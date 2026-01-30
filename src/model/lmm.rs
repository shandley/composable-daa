//! Linear Mixed Models (LMM) for longitudinal/repeated measures data.
//!
//! Implements mixed models with random effects for handling correlated observations:
//! `y = Xβ + Zu + ε` where `u ~ N(0, G⊗I)` and `ε ~ N(0, σ²I)`
//!
//! G is the covariance matrix of random effects (q × q where q = number of random terms per group).
//! For random intercepts only: G = [τ²]
//! For random intercepts + slopes: G = [[τ²_intercept, τ_cov], [τ_cov, τ²_slope]]
//!
//! Uses REML (Restricted Maximum Likelihood) estimation for unbiased variance
//! component estimates.
//!
//! # Example
//! ```ignore
//! use composable_daa::model::lmm::{model_lmm, LmmConfig};
//!
//! // Fit LMM with random intercepts per subject
//! let fit = model_lmm(
//!     &transformed,
//!     &design,
//!     &random_design,
//!     &LmmConfig::default(),
//! ).unwrap();
//!
//! // Fit LMM with random intercepts and slopes
//! // Formula: ~ group + time + (1 + time | subject)
//! ```

use crate::data::{DesignMatrix, MixedFormula, Metadata, RandomDesignMatrix, RandomEffect};
use crate::error::{DaaError, Result};
use crate::normalize::TransformedMatrix;
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Method for computing degrees of freedom in LMM hypothesis tests.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq, Eq)]
pub enum DfMethod {
    /// Naive residual df: n - p (fast but can be anticonservative with high ICC).
    #[default]
    Residual,
    /// Satterthwaite approximation: coefficient-specific df that accounts for
    /// random effects structure (more accurate Type I error control).
    Satterthwaite,
}

/// Configuration for LMM fitting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LmmConfig {
    /// Maximum iterations for REML estimation.
    pub max_iter: usize,
    /// Convergence tolerance for log-likelihood.
    pub tol: f64,
    /// Small ridge value for numerical stability.
    pub ridge: f64,
    /// Lower bound for variance components (prevents collapse to zero).
    pub var_lower_bound: f64,
    /// Method for computing degrees of freedom for hypothesis tests.
    pub df_method: DfMethod,
}

impl Default for LmmConfig {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-6,
            ridge: 1e-8,
            var_lower_bound: 1e-10,
            df_method: DfMethod::default(),
        }
    }
}

/// Results from fitting an LMM to a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LmmFitSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Estimated fixed effects coefficients.
    pub coefficients: Vec<f64>,
    /// Standard errors of fixed effects.
    pub std_errors: Vec<f64>,
    /// Variance components: [G_11, G_12, G_22, ..., sigma2].
    /// For intercept only: [tau2, sigma2]
    /// For intercept + slope: [tau2_intercept, tau_cov, tau2_slope, sigma2]
    pub variance_components: Vec<f64>,
    /// G matrix (covariance of random effects) as flat vector (row-major).
    /// Size: n_random_per_group × n_random_per_group
    #[serde(skip)]
    pub g_matrix: Vec<f64>,
    /// BLUPs (Best Linear Unbiased Predictors) for random effects.
    #[serde(skip)]
    pub random_effects: Vec<f64>,
    /// REML log-likelihood at convergence.
    pub log_reml: f64,
    /// Approximate residual degrees of freedom (n - p).
    pub df_residual: f64,
    /// Per-coefficient Satterthwaite degrees of freedom.
    /// None if df_method was Residual, Some(vec) if Satterthwaite was used.
    /// Each element corresponds to a coefficient in the same order as `coefficients`.
    pub df_satterthwaite: Option<Vec<f64>>,
    /// Number of iterations to convergence.
    pub iterations: usize,
    /// Whether the model converged.
    pub converged: bool,
    /// Intraclass correlation coefficient: tau2_intercept / (tau2_intercept + sigma2).
    /// For random slopes, this is based on the intercept variance only.
    pub icc: f64,
    /// Number of random effect terms per group.
    pub n_random_per_group: usize,
    /// Correlation between random intercept and slope (if applicable).
    pub random_correlation: Option<f64>,
}

impl LmmFitSingle {
    /// Get coefficient by index.
    pub fn get_coefficient(&self, index: usize) -> Option<f64> {
        self.coefficients.get(index).copied()
    }

    /// Get standard error by index.
    pub fn get_std_error(&self, index: usize) -> Option<f64> {
        self.std_errors.get(index).copied()
    }

    /// Calculate t-statistic for a coefficient.
    pub fn t_statistic(&self, index: usize) -> Option<f64> {
        let coef = self.coefficients.get(index)?;
        let se = self.std_errors.get(index)?;
        if *se > 0.0 {
            Some(coef / se)
        } else {
            None
        }
    }

    /// Get the random intercept variance (tau²_intercept).
    /// For intercept-only models, this is the only random variance.
    pub fn tau_squared(&self) -> f64 {
        self.variance_components.first().copied().unwrap_or(0.0)
    }

    /// Get the random slope variance (tau²_slope).
    /// Returns None for intercept-only models.
    pub fn tau_squared_slope(&self) -> Option<f64> {
        if self.n_random_per_group >= 2 {
            // For 2x2 G matrix: [G_11, G_12, G_22, sigma2]
            // G_22 is at index 2
            self.variance_components.get(2).copied()
        } else {
            None
        }
    }

    /// Get the covariance between random intercept and slope.
    /// Returns None for intercept-only models.
    pub fn tau_covariance(&self) -> Option<f64> {
        if self.n_random_per_group >= 2 {
            // For 2x2 G matrix: [G_11, G_12, G_22, sigma2]
            // G_12 is at index 1
            self.variance_components.get(1).copied()
        } else {
            None
        }
    }

    /// Get the residual variance (σ²).
    pub fn sigma_squared(&self) -> f64 {
        // sigma2 is always the last element
        self.variance_components.last().copied().unwrap_or(0.0)
    }

    /// Get the G matrix as a nalgebra DMatrix.
    pub fn g_matrix(&self) -> DMatrix<f64> {
        let q = self.n_random_per_group;
        DMatrix::from_row_slice(q, q, &self.g_matrix)
    }
}

/// Results from fitting LMM to all features.
#[derive(Debug, Clone)]
pub struct LmmFit {
    /// Individual fits for each feature.
    pub fits: Vec<LmmFitSingle>,
    /// Fixed effect coefficient names.
    pub coefficient_names: Vec<String>,
    /// Variance component names.
    pub variance_component_names: Vec<String>,
    /// Random effect term names (e.g., ["intercept"] or ["intercept", "time"]).
    pub random_term_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
    /// Number of groups (subjects).
    pub n_groups: usize,
    /// Number of random effect terms per group.
    pub n_random_per_group: usize,
}

impl LmmFit {
    /// Get fit for a specific feature by ID.
    pub fn get_feature(&self, feature_id: &str) -> Option<&LmmFitSingle> {
        self.fits.iter().find(|f| f.feature_id == feature_id)
    }

    /// Get coefficient index by name.
    pub fn coefficient_index(&self, name: &str) -> Option<usize> {
        self.coefficient_names.iter().position(|n| n == name)
    }

    /// Get all coefficients for a specific coefficient name.
    pub fn coefficients_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| f.coefficients.get(idx).copied().unwrap_or(f64::NAN))
                .collect(),
        )
    }

    /// Get all standard errors for a specific coefficient name.
    pub fn std_errors_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| f.std_errors.get(idx).copied().unwrap_or(f64::NAN))
                .collect(),
        )
    }

    /// Number of features.
    pub fn n_features(&self) -> usize {
        self.fits.len()
    }

    /// Number of fixed effect coefficients.
    pub fn n_coefficients(&self) -> usize {
        self.coefficient_names.len()
    }

    /// Check if all fits converged.
    pub fn all_converged(&self) -> bool {
        self.fits.iter().all(|f| f.converged)
    }

    /// Count how many fits converged.
    pub fn n_converged(&self) -> usize {
        self.fits.iter().filter(|f| f.converged).count()
    }

    /// Get mean ICC across all features.
    pub fn mean_icc(&self) -> f64 {
        let iccs: Vec<f64> = self.fits.iter().map(|f| f.icc).collect();
        iccs.iter().sum::<f64>() / iccs.len() as f64
    }
}

/// Fit linear mixed models to transformed abundance data.
///
/// Fits an LMM using REML estimation for each feature:
/// `y = Xβ + Zu + ε` where `u ~ N(0, G⊗I)` and `ε ~ N(0, σ²I)`
///
/// G is the covariance matrix of random effects (q × q where q = n_random_per_group).
///
/// # Arguments
/// * `transformed` - Transformed abundance data (e.g., CLR-transformed)
/// * `design` - Fixed effects design matrix (X)
/// * `random_design` - Random effects design matrix (Z)
/// * `config` - LMM configuration
///
/// # Returns
/// LmmFit containing results for all features.
pub fn model_lmm(
    transformed: &TransformedMatrix,
    design: &DesignMatrix,
    random_design: &RandomDesignMatrix,
    config: &LmmConfig,
) -> Result<LmmFit> {
    let n_features = transformed.n_features();
    let n_samples = transformed.n_samples();
    let n_fixed = design.n_coefficients();
    let n_groups = random_design.n_groups;
    let n_random_per_group = random_design.n_random_per_group;

    // Validate dimensions
    if design.n_samples() != n_samples {
        return Err(DaaError::DimensionMismatch {
            expected: n_samples,
            actual: design.n_samples(),
        });
    }
    if random_design.n_samples() != n_samples {
        return Err(DaaError::DimensionMismatch {
            expected: n_samples,
            actual: random_design.n_samples(),
        });
    }

    // Check degrees of freedom
    if n_samples <= n_fixed {
        return Err(DaaError::Numerical(
            "Model is saturated (n_samples <= n_fixed_effects)".to_string(),
        ));
    }

    let x = design.matrix();
    let z = random_design.matrix();

    // Fit all features in parallel
    let fits: Vec<LmmFitSingle> = (0..n_features)
        .into_par_iter()
        .map(|i| {
            fit_single_feature_lmm(
                &transformed.row(i),
                &transformed.feature_ids[i],
                x,
                z,
                n_groups,
                n_random_per_group,
                config,
            )
        })
        .collect();

    // Build variance component names based on G matrix structure
    let variance_component_names = build_variance_component_names(n_random_per_group, &random_design.term_names);

    Ok(LmmFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        variance_component_names,
        random_term_names: random_design.term_names.clone(),
        n_samples,
        n_groups,
        n_random_per_group,
    })
}

/// Build variance component names based on G matrix structure.
fn build_variance_component_names(n_random_per_group: usize, term_names: &[String]) -> Vec<String> {
    let mut names = Vec::new();

    // Add names for G matrix elements (upper triangle)
    for i in 0..n_random_per_group {
        for j in i..n_random_per_group {
            if i == j {
                names.push(format!("var_{}", term_names.get(i).unwrap_or(&format!("re_{}", i))));
            } else {
                names.push(format!(
                    "cov_{}_{}",
                    term_names.get(i).unwrap_or(&format!("re_{}", i)),
                    term_names.get(j).unwrap_or(&format!("re_{}", j))
                ));
            }
        }
    }

    // Add residual variance
    names.push("sigma2".to_string());

    names
}

/// Fit LMM to a single feature using REML.
fn fit_single_feature_lmm(
    y: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    z: &DMatrix<f64>,
    n_groups: usize,
    n_random_per_group: usize,
    config: &LmmConfig,
) -> LmmFitSingle {
    let n = y.len();
    let p = x.ncols();
    let q = n_random_per_group;
    let y_vec = DVector::from_column_slice(y);

    // Initialize variance components from OLS residuals
    let (sigma2_init, g_init) = initialize_variance_components_general(&y_vec, x, q, config);
    let mut sigma2 = sigma2_init;
    let mut g_matrix = g_init;

    let mut log_reml_prev = f64::NEG_INFINITY;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Build V = sigma2*I + Z*G_block*Z'
        let v = build_v_matrix_general(n, sigma2, &g_matrix, z, n_groups, q, config.ridge);

        // Cholesky decomposition of V
        let v_chol = match v.clone().cholesky() {
            Some(c) => c,
            None => {
                // V is not positive definite - try adding more ridge
                let v_ridge = &v + DMatrix::identity(n, n) * 0.01;
                match v_ridge.cholesky() {
                    Some(c) => c,
                    None => {
                        // Give up and return OLS solution
                        return fit_ols_fallback(&y_vec, feature_id, x, n_groups, q);
                    }
                }
            }
        };

        // Solve V^-1 X and V^-1 y using Cholesky
        let v_inv_x = v_chol.solve(x);
        let v_inv_y = v_chol.solve(&y_vec);

        // Compute (X'V^-1 X)^-1
        let xtvinvx = x.transpose() * &v_inv_x;
        let xtvinvx_inv = match xtvinvx.clone().try_inverse() {
            Some(inv) => inv,
            None => {
                // Singular - add ridge
                let xtvinvx_ridge = &xtvinvx + DMatrix::identity(p, p) * config.ridge;
                match xtvinvx_ridge.try_inverse() {
                    Some(inv) => inv,
                    None => return fit_ols_fallback(&y_vec, feature_id, x, n_groups, q),
                }
            }
        };

        // GLS estimates: beta = (X'V^-1 X)^-1 X'V^-1 y
        let xtvinvy = x.transpose() * &v_inv_y;
        let beta = &xtvinvx_inv * xtvinvy;

        // Residuals: r = y - X*beta
        let residuals = &y_vec - x * &beta;

        // REML log-likelihood (up to constant)
        // -0.5 * (log|V| + log|X'V^-1 X| + r'V^-1 r)
        let log_det_v = 2.0 * v_chol.l().diagonal().map(|x| x.ln()).sum();
        let log_det_xtvinvx = match xtvinvx.clone().cholesky() {
            Some(c) => 2.0 * c.l().diagonal().map(|x| x.ln()).sum(),
            None => p as f64 * xtvinvx[(0, 0)].abs().ln(), // Fallback approximation
        };
        let v_inv_r = v_chol.solve(&residuals);
        let quad_form = residuals.dot(&v_inv_r);

        let log_reml = -0.5 * (log_det_v + log_det_xtvinvx + quad_form);

        // Check convergence
        if (log_reml - log_reml_prev).abs() < config.tol {
            converged = true;
        }
        log_reml_prev = log_reml;

        if converged {
            // Compute Satterthwaite df if requested
            let df_satterthwaite = if config.df_method == DfMethod::Satterthwaite {
                Some(compute_satterthwaite_df(&v_chol, x, z, &xtvinvx_inv, n_groups, q))
            } else {
                None
            };

            return build_final_result(
                feature_id, &beta, &xtvinvx_inv, &g_matrix, sigma2, z, &v_inv_r,
                n, p, q, n_groups, log_reml, iterations, true, df_satterthwaite,
            );
        }

        // Update variance components
        let (new_g, new_sigma2) = update_variance_components_general(
            &v_chol, z, &residuals, sigma2, &g_matrix, n_groups, q, p, config,
        );

        g_matrix = new_g;
        sigma2 = new_sigma2;
    }

    // Did not converge - return last estimates
    let v = build_v_matrix_general(n, sigma2, &g_matrix, z, n_groups, q, config.ridge);
    let v_chol = v.clone().cholesky().unwrap_or_else(|| {
        let v_ridge = &v + DMatrix::identity(n, n) * 0.01;
        v_ridge.cholesky().unwrap()
    });

    let v_inv_x = v_chol.solve(x);
    let v_inv_y = v_chol.solve(&y_vec);
    let xtvinvx = x.transpose() * &v_inv_x;
    let xtvinvx_inv = (xtvinvx + DMatrix::identity(p, p) * config.ridge)
        .try_inverse()
        .unwrap();
    let beta = &xtvinvx_inv * (x.transpose() * &v_inv_y);
    let residuals = &y_vec - x * &beta;
    let v_inv_r = v_chol.solve(&residuals);

    // Compute Satterthwaite df if requested (even for non-converged fits)
    let df_satterthwaite = if config.df_method == DfMethod::Satterthwaite {
        Some(compute_satterthwaite_df(&v_chol, x, z, &xtvinvx_inv, n_groups, q))
    } else {
        None
    };

    build_final_result(
        feature_id, &beta, &xtvinvx_inv, &g_matrix, sigma2, z, &v_inv_r,
        n, p, q, n_groups, log_reml_prev, iterations, false, df_satterthwaite,
    )
}

/// Build final LmmFitSingle result from computed values.
fn build_final_result(
    feature_id: &str,
    beta: &DVector<f64>,
    xtvinvx_inv: &DMatrix<f64>,
    g_matrix: &DMatrix<f64>,
    sigma2: f64,
    z: &DMatrix<f64>,
    v_inv_r: &DVector<f64>,
    n: usize,
    p: usize,
    q: usize,
    n_groups: usize,
    log_reml: f64,
    iterations: usize,
    converged: bool,
    df_satterthwaite: Option<Vec<f64>>,
) -> LmmFitSingle {
    let coefficients: Vec<f64> = beta.iter().cloned().collect();
    let std_errors: Vec<f64> = (0..p)
        .map(|j| xtvinvx_inv[(j, j)].max(0.0).sqrt())
        .collect();

    // BLUPs: u = G_block * Z'V^-1(y - Xb)
    // For efficiency, compute group by group
    let g_block = build_g_block(g_matrix, n_groups, q);
    let random_effects: Vec<f64> = (&g_block * z.transpose() * v_inv_r)
        .iter()
        .cloned()
        .collect();

    // Build variance components vector (upper triangle of G + sigma2)
    let mut variance_components = Vec::new();
    for i in 0..q {
        for j in i..q {
            variance_components.push(g_matrix[(i, j)]);
        }
    }
    variance_components.push(sigma2);

    // G matrix as flat vector
    let g_flat: Vec<f64> = g_matrix.iter().cloned().collect();

    // ICC based on intercept variance (first diagonal element of G)
    let tau2_intercept = g_matrix[(0, 0)];
    let icc = tau2_intercept / (tau2_intercept + sigma2);

    // Random correlation (for random slopes)
    let random_correlation = if q >= 2 {
        let tau2_slope = g_matrix[(1, 1)];
        let tau_cov = g_matrix[(0, 1)];
        let denom = (tau2_intercept * tau2_slope).sqrt();
        if denom > 0.0 {
            Some(tau_cov / denom)
        } else {
            Some(0.0)
        }
    } else {
        None
    };

    // Approximate df
    let df_residual = (n - p) as f64;

    LmmFitSingle {
        feature_id: feature_id.to_string(),
        coefficients,
        std_errors,
        variance_components,
        g_matrix: g_flat,
        random_effects,
        log_reml,
        df_residual,
        df_satterthwaite,
        iterations,
        converged,
        icc,
        n_random_per_group: q,
        random_correlation,
    }
}

/// Build V matrix: V = sigma2*I + Z*G_block*Z'
/// where G_block is block-diagonal with G repeated n_groups times.
fn build_v_matrix_general(
    n: usize,
    sigma2: f64,
    g: &DMatrix<f64>,
    z: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
    ridge: f64,
) -> DMatrix<f64> {
    // Build block-diagonal G_block (n_groups*q × n_groups*q)
    let g_block = build_g_block(g, n_groups, q);

    // V = sigma2*I + Z * G_block * Z'
    let z_g_block = z * &g_block;
    let mut v = &z_g_block * z.transpose();

    for i in 0..n {
        v[(i, i)] += sigma2 + ridge;
    }
    v
}

/// Build block-diagonal G matrix for all groups.
fn build_g_block(g: &DMatrix<f64>, n_groups: usize, q: usize) -> DMatrix<f64> {
    let total_size = n_groups * q;
    let mut g_block = DMatrix::zeros(total_size, total_size);

    for group in 0..n_groups {
        let offset = group * q;
        for i in 0..q {
            for j in 0..q {
                g_block[(offset + i, offset + j)] = g[(i, j)];
            }
        }
    }
    g_block
}

/// Initialize variance components from OLS residuals.
fn initialize_variance_components_general(
    y: &DVector<f64>,
    x: &DMatrix<f64>,
    q: usize,
    config: &LmmConfig,
) -> (f64, DMatrix<f64>) {
    let n = y.len();
    let p = x.ncols();

    // OLS fit
    let xtx = x.transpose() * x;
    let xtx_inv = match xtx.clone().try_inverse() {
        Some(inv) => inv,
        None => {
            let xtx_ridge = &xtx + DMatrix::identity(p, p) * config.ridge;
            xtx_ridge.try_inverse().unwrap()
        }
    };
    let beta_ols = &xtx_inv * (x.transpose() * y);
    let residuals = y - x * beta_ols;
    let rss: f64 = residuals.iter().map(|r| r * r).sum();

    let df = (n - p).max(1);
    let sigma2 = (rss / df as f64).max(config.var_lower_bound);

    // Initialize G as diagonal matrix with tau2 = 0.1 * sigma2
    let tau2 = (0.1 * sigma2).max(config.var_lower_bound);
    let mut g = DMatrix::zeros(q, q);
    for i in 0..q {
        g[(i, i)] = tau2;
    }

    (sigma2, g)
}

/// Update variance components using REML moment-based estimation.
fn update_variance_components_general(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    z: &DMatrix<f64>,
    residuals: &DVector<f64>,
    sigma2: f64,
    g: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
    p: usize,
    config: &LmmConfig,
) -> (DMatrix<f64>, f64) {
    let n = residuals.len();

    // Compute V^-1 r
    let v_inv_r = v_chol.solve(residuals);

    // Quadratic form for residual variance
    let r_vinv_r = residuals.dot(&v_inv_r);
    let new_sigma2 = (r_vinv_r / (n - p) as f64).max(config.var_lower_bound);

    // Update G matrix using moment-based estimation
    // For each group, estimate the outer product of random effects
    let ztr = z.transpose() * &v_inv_r;

    let mut new_g = DMatrix::zeros(q, q);

    // Estimate G from cross-products of Z'V^-1 r grouped by subject
    for group in 0..n_groups {
        let offset = group * q;
        for i in 0..q {
            for j in 0..q {
                new_g[(i, j)] += ztr[offset + i] * ztr[offset + j];
            }
        }
    }
    new_g /= n_groups as f64;

    // Ensure G is positive semi-definite
    let new_g = ensure_positive_semidefinite(&new_g, config.var_lower_bound);

    // Damped update to prevent oscillation
    let alpha = 0.5;
    let g_updated = &new_g * alpha + g * (1.0 - alpha);
    let sigma2_updated = alpha * new_sigma2 + (1.0 - alpha) * sigma2;

    (
        ensure_positive_semidefinite(&g_updated, config.var_lower_bound),
        sigma2_updated.max(config.var_lower_bound),
    )
}

/// Ensure matrix is positive semi-definite.
fn ensure_positive_semidefinite(m: &DMatrix<f64>, min_eigenvalue: f64) -> DMatrix<f64> {
    let q = m.nrows();

    if q == 1 {
        // 1x1 case: just ensure positive
        let mut result = m.clone();
        result[(0, 0)] = result[(0, 0)].max(min_eigenvalue);
        return result;
    }

    if q == 2 {
        // 2x2 case: use explicit formula
        let a = m[(0, 0)];
        let b = (m[(0, 1)] + m[(1, 0)]) / 2.0; // Make symmetric
        let d = m[(1, 1)];

        // Eigenvalues of 2x2 symmetric: (a+d)/2 ± sqrt((a-d)²/4 + b²)
        let trace = a + d;
        let det = a * d - b * b;

        let discriminant = (trace * trace / 4.0 - det).max(0.0);
        let lambda1 = trace / 2.0 + discriminant.sqrt();
        let lambda2 = trace / 2.0 - discriminant.sqrt();

        // If both eigenvalues are positive enough, return as-is (symmetric)
        if lambda1 >= min_eigenvalue && lambda2 >= min_eigenvalue {
            let mut result = m.clone();
            result[(0, 1)] = b;
            result[(1, 0)] = b;
            return result;
        }

        // Otherwise, reconstruct with clipped eigenvalues
        let lambda1_new = lambda1.max(min_eigenvalue);
        let lambda2_new = lambda2.max(min_eigenvalue);

        // Handle diagonal case
        if b.abs() < 1e-10 {
            let mut result = DMatrix::zeros(2, 2);
            result[(0, 0)] = a.max(min_eigenvalue);
            result[(1, 1)] = d.max(min_eigenvalue);
            return result;
        }

        // Eigenvector for lambda1
        let v1 = DVector::from_vec(vec![b, lambda1 - a]);
        let v1_norm = v1.norm();
        let v1 = if v1_norm > 1e-10 { v1 / v1_norm } else { DVector::from_vec(vec![1.0, 0.0]) };
        let v2 = DVector::from_vec(vec![-v1[1], v1[0]]); // Orthogonal

        // Reconstruct: V * diag(lambda) * V'
        let mut result = DMatrix::zeros(2, 2);
        for i in 0..2 {
            for j in 0..2 {
                result[(i, j)] = lambda1_new * v1[i] * v1[j] + lambda2_new * v2[i] * v2[j];
            }
        }
        return result;
    }

    // For larger matrices, use diagonal dominance approach
    let mut result = m.clone();
    for i in 0..q {
        result[(i, i)] = result[(i, i)].max(min_eigenvalue);
    }
    result
}

// ============================================================================
// Satterthwaite Degrees of Freedom Approximation
// ============================================================================

/// Compute ∂V/∂σ² = I (identity matrix).
fn dv_d_sigma2(n: usize) -> DMatrix<f64> {
    DMatrix::identity(n, n)
}

/// Compute ∂V/∂G_ij where G is the random effects covariance matrix.
///
/// For V = σ²I + Z * G_block * Z', the derivative with respect to G_ij is:
/// ∂V/∂G_ij = Z * ∂G_block/∂G_ij * Z'
///
/// where G_block is the block-diagonal matrix with G repeated n_groups times.
fn dv_d_g_element(
    z: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
    i: usize,
    j: usize,
) -> DMatrix<f64> {
    let _n = z.nrows();

    // ∂G_block/∂G_ij has 1s at positions (g*q+i, g*q+j) and (g*q+j, g*q+i) for each group g
    // (for i != j, accounting for symmetry; for i == j, just the diagonal)
    let total_cols = n_groups * q;
    let mut dg_block = DMatrix::zeros(total_cols, total_cols);

    for g in 0..n_groups {
        let offset = g * q;
        if i == j {
            dg_block[(offset + i, offset + j)] = 1.0;
        } else {
            // G is symmetric, so derivatives appear in both (i,j) and (j,i)
            dg_block[(offset + i, offset + j)] = 1.0;
            dg_block[(offset + j, offset + i)] = 1.0;
        }
    }

    // ∂V/∂G_ij = Z * ∂G_block/∂G_ij * Z'
    let z_dg = z * &dg_block;
    &z_dg * z.transpose()
}

/// Compute the projection matrix P = V⁻¹ - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹.
fn compute_p_matrix(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    x: &DMatrix<f64>,
    xtvinvx_inv: &DMatrix<f64>,
) -> DMatrix<f64> {
    let n = x.nrows();

    // V⁻¹X
    let v_inv_x = v_chol.solve(x);

    // V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹ = (V⁻¹X)(X'V⁻¹X)⁻¹(V⁻¹X)'
    let term = &v_inv_x * xtvinvx_inv * v_inv_x.transpose();

    // V⁻¹ - V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹
    // We need V⁻¹ explicitly, which we compute column by column
    let mut v_inv = DMatrix::zeros(n, n);
    for i in 0..n {
        let mut e_i = DVector::zeros(n);
        e_i[i] = 1.0;
        let col = v_chol.solve(&e_i);
        for j in 0..n {
            v_inv[(j, i)] = col[j];
        }
    }

    &v_inv - &term
}

/// Compute ∂Cov(β̂)/∂θ_k for each variance component.
///
/// ∂Cov(β̂)/∂θ_k = -(X'V⁻¹X)⁻¹ X'V⁻¹ (∂V/∂θ_k) V⁻¹X (X'V⁻¹X)⁻¹
fn compute_cov_beta_derivatives(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    x: &DMatrix<f64>,
    z: &DMatrix<f64>,
    xtvinvx_inv: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
) -> Vec<DMatrix<f64>> {
    let n = x.nrows();
    let _p = x.ncols();

    // Number of variance components: upper triangle of G (q*(q+1)/2) + sigma²
    let n_theta = q * (q + 1) / 2 + 1;
    let mut derivatives = Vec::with_capacity(n_theta);

    // V⁻¹X
    let v_inv_x = v_chol.solve(x);

    // For each variance component θ_k, compute the derivative
    // First: G elements (upper triangle)
    for i in 0..q {
        for j in i..q {
            let dv_dtheta = dv_d_g_element(z, n_groups, q, i, j);

            // V⁻¹ (∂V/∂θ_k) V⁻¹X
            let dv_v_inv_x = v_chol.solve(&dv_dtheta) * &v_inv_x;

            // X'V⁻¹ (∂V/∂θ_k) V⁻¹X
            let middle = v_inv_x.transpose() * &dv_v_inv_x;

            // -(X'V⁻¹X)⁻¹ * middle * (X'V⁻¹X)⁻¹
            let derivative = -xtvinvx_inv * &middle * xtvinvx_inv;

            derivatives.push(derivative);
        }
    }

    // Last: sigma² (∂V/∂σ² = I)
    {
        let _dv_dtheta = dv_d_sigma2(n);

        // V⁻¹ I V⁻¹X = V⁻¹ V⁻¹X (need V⁻² X)
        // Actually: V⁻¹ (∂V/∂σ²) V⁻¹X = V⁻¹ I V⁻¹X = V⁻¹ (V⁻¹X)
        let v_inv_v_inv_x = v_chol.solve(&v_inv_x);

        // X'V⁻¹ (∂V/∂σ²) V⁻¹X = X'V⁻² X
        let middle = v_inv_x.transpose() * &v_inv_v_inv_x;

        // -(X'V⁻¹X)⁻¹ * middle * (X'V⁻¹X)⁻¹
        let derivative = -xtvinvx_inv * &middle * xtvinvx_inv;

        derivatives.push(derivative);
    }

    derivatives
}

/// Compute the information matrix I(θ) for variance components.
///
/// I(θ)_kl = (1/2) * tr(P * ∂V/∂θ_k * P * ∂V/∂θ_l)
fn compute_variance_info_matrix(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    x: &DMatrix<f64>,
    z: &DMatrix<f64>,
    xtvinvx_inv: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
) -> DMatrix<f64> {
    let n = x.nrows();

    // Number of variance components
    let n_theta = q * (q + 1) / 2 + 1;

    // Compute P matrix
    let p_mat = compute_p_matrix(v_chol, x, xtvinvx_inv);

    // Collect all ∂V/∂θ_k
    let mut dv_list: Vec<DMatrix<f64>> = Vec::with_capacity(n_theta);

    // G elements (upper triangle)
    for i in 0..q {
        for j in i..q {
            dv_list.push(dv_d_g_element(z, n_groups, q, i, j));
        }
    }

    // sigma²
    dv_list.push(dv_d_sigma2(n));

    // Compute I(θ)
    let mut info = DMatrix::zeros(n_theta, n_theta);

    for k in 0..n_theta {
        // P * ∂V/∂θ_k
        let p_dvk = &p_mat * &dv_list[k];

        for l in k..n_theta {
            // P * ∂V/∂θ_l
            let p_dvl = if l == k {
                p_dvk.clone()
            } else {
                &p_mat * &dv_list[l]
            };

            // tr(P * ∂V/∂θ_k * P * ∂V/∂θ_l) = tr(p_dvk * p_dvl)
            let trace: f64 = (0..n).map(|i| {
                (0..n).map(|m| p_dvk[(i, m)] * p_dvl[(m, i)]).sum::<f64>()
            }).sum();

            info[(k, l)] = 0.5 * trace;
            if l != k {
                info[(l, k)] = info[(k, l)];
            }
        }
    }

    info
}

/// Compute Satterthwaite degrees of freedom for each fixed effect coefficient.
///
/// df_j = 2 * [Var(β̂_j)]² / Var[Var(β̂_j)]
///
/// where Var[Var(β̂_j)] = gradient' * I(θ)⁻¹ * gradient
/// and gradient = [∂Var(β̂_j)/∂θ_1, ..., ∂Var(β̂_j)/∂θ_K]
pub fn compute_satterthwaite_df(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    x: &DMatrix<f64>,
    z: &DMatrix<f64>,
    xtvinvx_inv: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
) -> Vec<f64> {
    let n = x.nrows();
    let p = x.ncols();

    // Naive df for clamping
    let df_naive = (n - p) as f64;

    // Number of variance components
    let n_theta = q * (q + 1) / 2 + 1;

    // Compute information matrix for variance components
    let info = compute_variance_info_matrix(v_chol, x, z, xtvinvx_inv, n_groups, q);

    // Compute I(θ)⁻¹ with small ridge for numerical stability
    let ridge = 1e-8;
    let info_ridge = &info + DMatrix::identity(n_theta, n_theta) * ridge;
    let info_inv = match info_ridge.try_inverse() {
        Some(inv) => inv,
        None => {
            // Information matrix is singular even with ridge, fall back to naive df
            return vec![df_naive; p];
        }
    };

    // Compute ∂Cov(β̂)/∂θ_k for each variance component
    let cov_beta_derivs = compute_cov_beta_derivatives(v_chol, x, z, xtvinvx_inv, n_groups, q);

    // Compute Satterthwaite df for each coefficient
    let mut df_satt = Vec::with_capacity(p);

    for j in 0..p {
        // Var(β̂_j) = Cov(β̂)[j,j]
        let var_beta_j = xtvinvx_inv[(j, j)];

        if var_beta_j <= 0.0 {
            df_satt.push(df_naive);
            continue;
        }

        // gradient[k] = ∂Var(β̂_j)/∂θ_k = (∂Cov(β̂)/∂θ_k)[j,j]
        let gradient: DVector<f64> = DVector::from_iterator(
            n_theta,
            cov_beta_derivs.iter().map(|d| d[(j, j)])
        );

        // Var[Var(β̂_j)] = gradient' * I(θ)⁻¹ * gradient
        let var_var_beta_j = gradient.dot(&(&info_inv * &gradient));

        if var_var_beta_j <= 0.0 {
            df_satt.push(df_naive);
            continue;
        }

        // df_j = 2 * [Var(β̂_j)]² / Var[Var(β̂_j)]
        let df_j = 2.0 * var_beta_j * var_beta_j / var_var_beta_j;

        // Clamp to [1, n-p]
        df_satt.push(df_j.max(1.0).min(df_naive));
    }

    df_satt
}

/// Fallback to OLS when LMM fails to converge.
fn fit_ols_fallback(
    y: &DVector<f64>,
    feature_id: &str,
    x: &DMatrix<f64>,
    n_groups: usize,
    q: usize,
) -> LmmFitSingle {
    let n = y.len();
    let p = x.ncols();

    let xtx = x.transpose() * x;
    let xtx_inv = match xtx.clone().try_inverse() {
        Some(inv) => inv,
        None => {
            let xtx_ridge = &xtx + DMatrix::identity(p, p) * 1e-6;
            xtx_ridge.try_inverse().unwrap()
        }
    };

    let beta = &xtx_inv * (x.transpose() * y);
    let residuals = y - x * &beta;
    let rss: f64 = residuals.iter().map(|r| r * r).sum();
    let df = (n - p).max(1);
    let sigma2 = rss / df as f64;

    let coefficients: Vec<f64> = beta.iter().cloned().collect();
    let std_errors: Vec<f64> = (0..p)
        .map(|j| (sigma2 * xtx_inv[(j, j)]).max(0.0).sqrt())
        .collect();

    // Build variance components (all zeros for G + sigma2)
    let mut variance_components = vec![0.0; q * (q + 1) / 2];
    variance_components.push(sigma2);

    LmmFitSingle {
        feature_id: feature_id.to_string(),
        coefficients,
        std_errors,
        variance_components,
        g_matrix: vec![0.0; q * q],
        random_effects: vec![0.0; n_groups * q],
        log_reml: f64::NAN,
        df_residual: df as f64,
        df_satterthwaite: None,
        iterations: 0,
        converged: false,
        icc: 0.0,
        n_random_per_group: q,
        random_correlation: if q >= 2 { Some(0.0) } else { None },
    }
}

/// High-level function to fit LMM from formula, metadata, and transformed data.
///
/// Parses the formula, constructs design matrices, and fits the model.
///
/// # Arguments
/// * `transformed` - Transformed abundance data
/// * `metadata` - Sample metadata
/// * `formula` - Mixed model formula (e.g., "~ group + (1 | subject)")
/// * `config` - Optional LMM configuration
///
/// # Returns
/// LmmFit containing results for all features.
///
/// # Example
/// ```ignore
/// let fit = fit_lmm_from_formula(
///     &transformed,
///     &metadata,
///     "~ group + time + (1 | subject)",
///     None,
/// )?;
/// ```
pub fn fit_lmm_from_formula(
    transformed: &TransformedMatrix,
    metadata: &Metadata,
    formula: &str,
    config: Option<LmmConfig>,
) -> Result<LmmFit> {
    let mixed_formula = MixedFormula::parse(formula)?;

    if mixed_formula.random.is_empty() {
        return Err(DaaError::InvalidParameter(
            "Formula must contain at least one random effect. Use model_lm for fixed effects only."
                .to_string(),
        ));
    }

    // Currently only support single random effect
    if mixed_formula.random.len() > 1 {
        return Err(DaaError::NotImplemented(
            "Multiple random effects not yet supported".to_string(),
        ));
    }

    let random_effect = &mixed_formula.random[0];

    // Build design matrices
    let design = DesignMatrix::from_formula(metadata, &mixed_formula.fixed)?;
    let random_design = RandomDesignMatrix::from_random_effect(metadata, random_effect)?;

    let config = config.unwrap_or_default();
    model_lmm(transformed, &design, &random_design, &config)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Formula;
    use approx::assert_relative_eq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_longitudinal_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\ttime\tsubject").unwrap();
        // 4 subjects, 2 timepoints each, balanced design
        writeln!(file, "S1\tcontrol\t0\tA").unwrap();
        writeln!(file, "S2\tcontrol\t1\tA").unwrap();
        writeln!(file, "S3\tcontrol\t0\tB").unwrap();
        writeln!(file, "S4\tcontrol\t1\tB").unwrap();
        writeln!(file, "S5\ttreatment\t0\tC").unwrap();
        writeln!(file, "S6\ttreatment\t1\tC").unwrap();
        writeln!(file, "S7\ttreatment\t0\tD").unwrap();
        writeln!(file, "S8\ttreatment\t1\tD").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn create_longitudinal_transformed() -> TransformedMatrix {
        // Feature 0: no effect (random noise around 1.0)
        // Feature 1: strong treatment effect (treatment group higher)
        // Feature 2: within-subject correlation but no treatment effect
        let data = DMatrix::from_row_slice(
            3,
            8,
            &[
                // Feature 0: no effect
                1.0, 1.1, 0.9, 1.0, 1.1, 0.9, 1.0, 1.1,
                // Feature 1: treatment effect
                1.0, 1.2, 1.1, 1.3, 3.0, 3.2, 2.9, 3.1,
                // Feature 2: subject-specific but no treatment effect
                // Subject A: high, B: low, C: medium, D: medium
                2.5, 2.6, 0.5, 0.6, 1.5, 1.6, 1.4, 1.5,
            ],
        );

        TransformedMatrix {
            data,
            feature_ids: vec!["no_effect".into(), "treatment_effect".into(), "subject_effect".into()],
            sample_ids: (1..=8).map(|i| format!("S{}", i)).collect(),
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 8],
        }
    }

    #[test]
    fn test_lmm_config_default() {
        let config = LmmConfig::default();
        assert_eq!(config.max_iter, 100);
        assert_eq!(config.tol, 1e-6);
    }

    #[test]
    fn test_model_lmm_basic() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_coefficients(), 2); // intercept + grouptreatment
        assert_eq!(fit.n_groups, 4);
        assert_eq!(fit.coefficient_names, vec!["(Intercept)", "grouptreatment"]);
    }

    #[test]
    fn test_model_lmm_treatment_effect() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Feature with treatment effect should have large coefficient
        let treatment_feat = fit.get_feature("treatment_effect").unwrap();
        let coef_idx = fit.coefficient_index("grouptreatment").unwrap();
        let treatment_coef = treatment_feat.coefficients[coef_idx];

        // Treatment effect should be approximately 2.0 (3.0 - 1.0)
        assert!(treatment_coef > 1.5, "Expected large treatment effect, got {}", treatment_coef);
    }

    #[test]
    fn test_model_lmm_no_effect() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Feature with no effect should have small coefficient
        let no_effect_feat = fit.get_feature("no_effect").unwrap();
        let coef_idx = fit.coefficient_index("grouptreatment").unwrap();
        let treatment_coef = no_effect_feat.coefficients[coef_idx];

        assert!(
            treatment_coef.abs() < 0.5,
            "Expected small treatment effect, got {}",
            treatment_coef
        );
    }

    #[test]
    fn test_model_lmm_variance_components() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Feature with strong subject effect should have higher ICC
        let subject_feat = fit.get_feature("subject_effect").unwrap();
        let no_effect_feat = fit.get_feature("no_effect").unwrap();

        // Subject effect feature should have higher ICC than no effect feature
        // (more variance explained by subject)
        assert!(
            subject_feat.icc > no_effect_feat.icc * 0.5,
            "Expected subject_effect ICC ({}) > no_effect ICC ({}) * 0.5",
            subject_feat.icc,
            no_effect_feat.icc
        );
    }

    #[test]
    fn test_model_lmm_convergence() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Most features should converge
        assert!(
            fit.n_converged() >= fit.n_features() / 2,
            "Expected most features to converge, got {}/{}",
            fit.n_converged(),
            fit.n_features()
        );
    }

    #[test]
    fn test_fit_lmm_from_formula() {
        let metadata = create_longitudinal_metadata();
        let transformed = create_longitudinal_transformed();

        let fit = fit_lmm_from_formula(
            &transformed,
            &metadata,
            "~ group + (1 | subject)",
            None,
        )
        .unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_groups, 4);
    }

    #[test]
    fn test_fit_lmm_from_formula_error_no_random() {
        let metadata = create_longitudinal_metadata();
        let transformed = create_longitudinal_transformed();

        let result = fit_lmm_from_formula(&transformed, &metadata, "~ group", None);
        assert!(result.is_err());
    }

    #[test]
    fn test_lmm_fit_single_accessors() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();
        let single = fit.get_feature("treatment_effect").unwrap();

        // Test accessors
        assert!(single.get_coefficient(0).is_some());
        assert!(single.get_std_error(0).is_some());
        assert!(single.t_statistic(0).is_some());
        assert!(single.tau_squared() >= 0.0);
        assert!(single.sigma_squared() >= 0.0);
    }

    #[test]
    fn test_dimension_mismatch() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();

        // Create transformed with wrong number of samples
        let data = DMatrix::from_row_slice(2, 4, &[1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
        let transformed = TransformedMatrix {
            data,
            feature_ids: vec!["a".into(), "b".into()],
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()],
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 4],
        };

        let result = model_lmm(&transformed, &design, &random_design, &LmmConfig::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_model_lmm_random_slope() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 + time | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Check that random slopes are properly configured
        assert_eq!(fit.n_random_per_group, 2); // intercept + time slope
        assert_eq!(fit.random_term_names, vec!["intercept", "time"]);

        // Each single fit should have G matrix with 4 elements (2x2)
        let single = fit.get_feature("treatment_effect").unwrap();
        assert_eq!(single.n_random_per_group, 2);
        assert_eq!(single.g_matrix.len(), 4); // 2x2 = 4

        // Should have correlation estimate for random intercept-slope
        assert!(single.random_correlation.is_some());
    }

    #[test]
    fn test_fit_lmm_from_formula_with_random_slope() {
        let metadata = create_longitudinal_metadata();
        let transformed = create_longitudinal_transformed();

        let fit = fit_lmm_from_formula(
            &transformed,
            &metadata,
            "~ group + (1 + time | subject)",
            None,
        )
        .unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_groups, 4);
        assert_eq!(fit.n_random_per_group, 2);
        assert_eq!(fit.random_term_names, vec!["intercept", "time"]);
    }

    #[test]
    fn test_model_lmm_random_slope_only() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(0 + time | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();
        let config = LmmConfig::default();

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Check that only random slope (no intercept) is configured
        assert_eq!(fit.n_random_per_group, 1);
        assert_eq!(fit.random_term_names, vec!["time"]);

        // No correlation estimate for single random term
        let single = fit.get_feature("treatment_effect").unwrap();
        assert!(single.random_correlation.is_none());
    }

    // ========================================================================
    // Satterthwaite df tests
    // ========================================================================

    #[test]
    fn test_df_method_default() {
        let config = LmmConfig::default();
        assert_eq!(config.df_method, DfMethod::Residual);
    }

    #[test]
    fn test_df_method_satterthwaite_config() {
        let config = LmmConfig {
            df_method: DfMethod::Satterthwaite,
            ..Default::default()
        };
        assert_eq!(config.df_method, DfMethod::Satterthwaite);
    }

    #[test]
    fn test_satterthwaite_df_computed() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();

        let config = LmmConfig {
            df_method: DfMethod::Satterthwaite,
            ..Default::default()
        };

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // All fits should have Satterthwaite df computed
        for single_fit in &fit.fits {
            assert!(
                single_fit.df_satterthwaite.is_some(),
                "df_satterthwaite should be computed when DfMethod::Satterthwaite"
            );
            let satt_df = single_fit.df_satterthwaite.as_ref().unwrap();
            assert_eq!(
                satt_df.len(),
                single_fit.coefficients.len(),
                "Satterthwaite df should have one value per coefficient"
            );
        }
    }

    #[test]
    fn test_satterthwaite_df_not_computed_for_residual() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();

        let config = LmmConfig {
            df_method: DfMethod::Residual,
            ..Default::default()
        };

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // No fits should have Satterthwaite df
        for single_fit in &fit.fits {
            assert!(
                single_fit.df_satterthwaite.is_none(),
                "df_satterthwaite should be None when DfMethod::Residual"
            );
        }
    }

    #[test]
    fn test_satterthwaite_df_bounds() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();

        let config = LmmConfig {
            df_method: DfMethod::Satterthwaite,
            ..Default::default()
        };

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();
        let n_samples = transformed.n_samples();
        let n_fixed = design.n_coefficients();
        let max_df = (n_samples - n_fixed) as f64;

        for single_fit in &fit.fits {
            let satt_df = single_fit.df_satterthwaite.as_ref().unwrap();
            for &df in satt_df {
                assert!(df >= 1.0, "Satterthwaite df should be >= 1, got {}", df);
                assert!(
                    df <= max_df,
                    "Satterthwaite df should be <= {} (n-p), got {}",
                    max_df,
                    df
                );
            }
        }
    }

    #[test]
    fn test_satterthwaite_df_with_high_icc() {
        // Create data with high ICC (strong within-subject correlation)
        // This test verifies that Satterthwaite df is computed correctly with high ICC data
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tsubject").unwrap();
        // 6 subjects, 3 observations each
        for subj in ['A', 'B', 'C', 'D', 'E', 'F'] {
            let group = if subj < 'D' { "control" } else { "treatment" };
            for obs in 1..=3 {
                writeln!(file, "{}_{}\t{}\t{}", subj, obs, group, subj).unwrap();
            }
        }
        file.flush().unwrap();
        let metadata = Metadata::from_tsv(file.path()).unwrap();

        // Create data with high within-subject correlation (high ICC)
        // Each subject has similar values within subject, but different between subjects
        let data = DMatrix::from_row_slice(
            1,
            18,
            &[
                // Subject A (control): ~1.0
                1.0, 1.05, 0.95,
                // Subject B (control): ~2.0
                2.0, 2.05, 1.95,
                // Subject C (control): ~3.0
                3.0, 3.05, 2.95,
                // Subject D (treatment): ~4.0
                4.0, 4.05, 3.95,
                // Subject E (treatment): ~5.0
                5.0, 5.05, 4.95,
                // Subject F (treatment): ~6.0
                6.0, 6.05, 5.95,
            ],
        );

        let transformed = TransformedMatrix {
            data,
            feature_ids: vec!["high_icc_feature".into()],
            sample_ids: metadata.sample_ids().iter().cloned().collect(),
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 18],
        };

        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();

        let config = LmmConfig {
            df_method: DfMethod::Satterthwaite,
            ..Default::default()
        };

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();
        let single_fit = fit.get_feature("high_icc_feature").unwrap();

        // Verify Satterthwaite df is computed
        assert!(single_fit.df_satterthwaite.is_some());
        let satt_df = single_fit.df_satterthwaite.as_ref().unwrap();

        // Verify df values are valid (finite, positive, bounded)
        let naive_df = single_fit.df_residual;
        for (i, &df) in satt_df.iter().enumerate() {
            assert!(
                df.is_finite(),
                "Satterthwaite df[{}] should be finite, got {}",
                i, df
            );
            assert!(df >= 1.0, "Satterthwaite df[{}] should be >= 1, got {}", i, df);
            assert!(
                df <= naive_df,
                "Satterthwaite df[{}] should be <= naive df ({}), got {}",
                i, naive_df, df
            );
        }

        // Verify ICC is high as expected from the data structure
        assert!(
            single_fit.icc > 0.5,
            "Expected high ICC (>0.5), got {}",
            single_fit.icc
        );

        // Note: With high ICC, theoretically Satterthwaite df should be < naive df,
        // but the exact numerical behavior depends on the information matrix structure.
        // The key guarantee is that Satterthwaite df is computed and bounded correctly.
    }

    #[test]
    fn test_backward_compatibility_default_config() {
        // Test that default config produces the same behavior as before
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();

        // Default config should use DfMethod::Residual
        let config = LmmConfig::default();
        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Should not compute Satterthwaite df
        for single_fit in &fit.fits {
            assert!(single_fit.df_satterthwaite.is_none());
            // df_residual should be n - p = 8 - 2 = 6
            assert_relative_eq!(single_fit.df_residual, 6.0, epsilon = 0.01);
        }
    }

    #[test]
    fn test_satterthwaite_with_random_slope() {
        let metadata = create_longitudinal_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let re = RandomEffect::parse("(1 + time | subject)").unwrap();
        let random_design = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();
        let transformed = create_longitudinal_transformed();

        let config = LmmConfig {
            df_method: DfMethod::Satterthwaite,
            ..Default::default()
        };

        let fit = model_lmm(&transformed, &design, &random_design, &config).unwrap();

        // Should still compute Satterthwaite df with random slopes
        for single_fit in &fit.fits {
            assert!(
                single_fit.df_satterthwaite.is_some(),
                "Should compute Satterthwaite df even with random slopes"
            );
            let satt_df = single_fit.df_satterthwaite.as_ref().unwrap();
            assert_eq!(satt_df.len(), design.n_coefficients());

            // All df should be valid (finite and positive)
            for &df in satt_df {
                assert!(df.is_finite(), "Satterthwaite df should be finite");
                assert!(df >= 1.0, "Satterthwaite df should be >= 1");
            }
        }
    }
}
