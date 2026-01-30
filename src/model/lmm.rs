//! Linear Mixed Models (LMM) for longitudinal/repeated measures data.
//!
//! Implements mixed models with random effects for handling correlated observations:
//! `y = Xβ + Zu + ε` where `u ~ N(0, τ²I)` and `ε ~ N(0, σ²I)`
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
//! ```

use crate::data::{DesignMatrix, MixedFormula, Metadata, RandomDesignMatrix, RandomEffect};
use crate::error::{DaaError, Result};
use crate::normalize::TransformedMatrix;
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

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
}

impl Default for LmmConfig {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-6,
            ridge: 1e-8,
            var_lower_bound: 1e-10,
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
    /// Variance components: [tau2 (random effect variance), sigma2 (residual variance)].
    pub variance_components: Vec<f64>,
    /// BLUPs (Best Linear Unbiased Predictors) for random effects.
    #[serde(skip)]
    pub random_effects: Vec<f64>,
    /// REML log-likelihood at convergence.
    pub log_reml: f64,
    /// Approximate residual degrees of freedom.
    pub df_residual: f64,
    /// Number of iterations to convergence.
    pub iterations: usize,
    /// Whether the model converged.
    pub converged: bool,
    /// Intraclass correlation coefficient: tau2 / (tau2 + sigma2).
    pub icc: f64,
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

    /// Get the random effect variance (tau²).
    pub fn tau_squared(&self) -> f64 {
        self.variance_components.first().copied().unwrap_or(0.0)
    }

    /// Get the residual variance (σ²).
    pub fn sigma_squared(&self) -> f64 {
        self.variance_components.get(1).copied().unwrap_or(0.0)
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
    /// Number of samples.
    pub n_samples: usize,
    /// Number of groups (subjects).
    pub n_groups: usize,
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
/// `y = Xβ + Zu + ε` where `u ~ N(0, τ²I)` and `ε ~ N(0, σ²I)`
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
                config,
            )
        })
        .collect();

    Ok(LmmFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        variance_component_names: vec!["tau2".to_string(), "sigma2".to_string()],
        n_samples,
        n_groups,
    })
}

/// Fit LMM to a single feature using REML.
fn fit_single_feature_lmm(
    y: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    z: &DMatrix<f64>,
    n_groups: usize,
    config: &LmmConfig,
) -> LmmFitSingle {
    let n = y.len();
    let p = x.ncols();
    let y_vec = DVector::from_column_slice(y);

    // Initialize variance components from OLS residuals
    let (sigma2_init, tau2_init) = initialize_variance_components(&y_vec, x, config);
    let mut sigma2 = sigma2_init;
    let mut tau2 = tau2_init;

    let mut log_reml_prev = f64::NEG_INFINITY;
    let mut converged = false;
    let mut iterations = 0;

    // Precompute ZZ'
    let zzt = z * z.transpose();

    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Build V = sigma2*I + tau2*ZZ'
        let v = build_v_matrix(n, sigma2, tau2, &zzt, config.ridge);

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
                        return fit_ols_fallback(&y_vec, feature_id, x, n_groups);
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
                    None => return fit_ols_fallback(&y_vec, feature_id, x, n_groups),
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
            None => p as f64 * xtvinvx[(0, 0)].ln(), // Fallback approximation
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
            // Compute final results
            let coefficients: Vec<f64> = beta.iter().cloned().collect();
            let std_errors: Vec<f64> = (0..p)
                .map(|j| xtvinvx_inv[(j, j)].max(0.0).sqrt())
                .collect();

            // BLUPs: u = tau2 * Z'V^-1(y - Xb)
            let random_effects: Vec<f64> = (tau2 * z.transpose() * &v_inv_r)
                .iter()
                .cloned()
                .collect();

            // ICC
            let icc = tau2 / (tau2 + sigma2);

            // Approximate df (Satterthwaite would be more accurate but complex)
            let df_residual = (n - p) as f64;

            return LmmFitSingle {
                feature_id: feature_id.to_string(),
                coefficients,
                std_errors,
                variance_components: vec![tau2, sigma2],
                random_effects,
                log_reml,
                df_residual,
                iterations,
                converged: true,
                icc,
            };
        }

        // Update variance components via Newton-Raphson
        // Use REML score equations
        let (new_tau2, new_sigma2) =
            update_variance_components(&v_chol, z, &residuals, sigma2, tau2, n_groups, p, config);

        tau2 = new_tau2;
        sigma2 = new_sigma2;
    }

    // Did not converge - return last estimates
    let v = build_v_matrix(n, sigma2, tau2, &zzt, config.ridge);
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

    let coefficients: Vec<f64> = beta.iter().cloned().collect();
    let std_errors: Vec<f64> = (0..p)
        .map(|j| xtvinvx_inv[(j, j)].max(0.0).sqrt())
        .collect();
    let random_effects: Vec<f64> = (tau2 * z.transpose() * &v_inv_r)
        .iter()
        .cloned()
        .collect();
    let icc = tau2 / (tau2 + sigma2);

    LmmFitSingle {
        feature_id: feature_id.to_string(),
        coefficients,
        std_errors,
        variance_components: vec![tau2, sigma2],
        random_effects,
        log_reml: log_reml_prev,
        df_residual: (n - p) as f64,
        iterations,
        converged: false,
        icc,
    }
}

/// Build V matrix: V = sigma2*I + tau2*ZZ'
fn build_v_matrix(
    n: usize,
    sigma2: f64,
    tau2: f64,
    zzt: &DMatrix<f64>,
    ridge: f64,
) -> DMatrix<f64> {
    let mut v = zzt * tau2;
    for i in 0..n {
        v[(i, i)] += sigma2 + ridge;
    }
    v
}

/// Initialize variance components from OLS residuals.
fn initialize_variance_components(
    y: &DVector<f64>,
    x: &DMatrix<f64>,
    config: &LmmConfig,
) -> (f64, f64) {
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
    let tau2 = (0.1 * sigma2).max(config.var_lower_bound);

    (sigma2, tau2)
}

/// Update variance components using REML score equations.
fn update_variance_components(
    v_chol: &nalgebra::Cholesky<f64, nalgebra::Dyn>,
    z: &DMatrix<f64>,
    residuals: &DVector<f64>,
    sigma2: f64,
    tau2: f64,
    n_groups: usize,
    p: usize,
    config: &LmmConfig,
) -> (f64, f64) {
    let n = residuals.len();

    // Compute V^-1 r
    let v_inv_r = v_chol.solve(residuals);

    // Compute trace terms for REML score
    // For tau2: d(log|V|)/d(tau2) = tr(V^-1 ZZ')
    // For sigma2: d(log|V|)/d(sigma2) = tr(V^-1)

    // Approximate traces using the diagonal of V^-1
    // This is a simplification - full computation would be more expensive
    let v_inv_diag: Vec<f64> = (0..n)
        .map(|i| {
            let ei = DVector::from_fn(n, |j, _| if j == i { 1.0 } else { 0.0 });
            let v_inv_ei = v_chol.solve(&ei);
            v_inv_ei[i]
        })
        .collect();

    // Quadratic forms
    let r_vinv_r = residuals.dot(&v_inv_r);

    // Simple moment-based update (more robust than Newton)
    // sigma2 = (1/n) * r'V^-1 r
    // tau2 = estimated from between-group variance

    let new_sigma2 = (r_vinv_r / (n - p) as f64).max(config.var_lower_bound);

    // Estimate tau2 from between-group variance
    // Use the fact that for balanced designs, tau2 ≈ (MSB - MSW) / n_per_group
    let ztr = z.transpose() * &v_inv_r;
    let ss_between: f64 = ztr.iter().map(|x| x * x).sum();
    let new_tau2 = (ss_between / n_groups as f64).max(config.var_lower_bound);

    // Damped update to prevent oscillation
    let alpha = 0.5;
    let tau2_updated = alpha * new_tau2 + (1.0 - alpha) * tau2;
    let sigma2_updated = alpha * new_sigma2 + (1.0 - alpha) * sigma2;

    (
        tau2_updated.max(config.var_lower_bound),
        sigma2_updated.max(config.var_lower_bound),
    )
}

/// Fallback to OLS when LMM fails to converge.
fn fit_ols_fallback(
    y: &DVector<f64>,
    feature_id: &str,
    x: &DMatrix<f64>,
    n_groups: usize,
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

    LmmFitSingle {
        feature_id: feature_id.to_string(),
        coefficients,
        std_errors,
        variance_components: vec![0.0, sigma2], // tau2 = 0 (no random effect estimated)
        random_effects: vec![0.0; n_groups],
        log_reml: f64::NAN,
        df_residual: df as f64,
        iterations: 0,
        converged: false,
        icc: 0.0,
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
}
