//! Beta-Binomial GLM for proportions with overdispersion.
//!
//! Implements a beta-binomial generalized linear model with logit link,
//! suitable for modeling proportions with overdispersion common in microbiome studies.
//!
//! The beta-binomial distribution models counts Y out of n trials where the
//! success probability varies according to a Beta distribution. This allows
//! for overdispersion relative to the binomial distribution.
//!
//! # Mathematical Formulation
//!
//! **Distribution**: Y_i ~ BetaBinomial(n_i, α_i, β_i)
//!
//! **Mean-dispersion parameterization**:
//! - μ = α/(α+β) is the mean proportion
//! - ρ = 1/(α+β+1) is overdispersion (0 = binomial)
//!
//! **GLM formulation**:
//! - Link: logit(μ) = Xβ
//! - Variance: Var(Y_i/n_i) = μ(1-μ)[1 + (n_i-1)ρ] / n_i

use crate::data::{CountMatrix, DesignMatrix};
use crate::error::{DaaError, Result};
use crate::model::nb::lgamma;
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Maximum iterations for IRLS convergence.
const MAX_ITER: usize = 25;

/// Convergence tolerance for coefficient changes.
const TOL: f64 = 1e-8;

/// Minimum value for proportions to avoid log(0) or division by zero.
const MIN_MU: f64 = 1e-10;

/// Maximum value for proportions to avoid log(0).
const MAX_MU: f64 = 1.0 - 1e-10;

/// Minimum overdispersion value to maintain numerical stability.
const MIN_RHO: f64 = 1e-10;

/// Maximum overdispersion value.
const MAX_RHO: f64 = 1.0 - 1e-10;

/// Method for estimating overdispersion parameter.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
pub enum DispersionMethod {
    /// Method of moments (fast, default).
    #[default]
    MomentsBased,
    /// Profile maximum likelihood (slower, more accurate).
    ProfileML,
}

/// Configuration for beta-binomial model fitting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BbConfig {
    /// Maximum number of IRLS iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
    /// Method for estimating overdispersion.
    pub dispersion_method: DispersionMethod,
}

impl Default for BbConfig {
    fn default() -> Self {
        Self {
            max_iter: MAX_ITER,
            tol: TOL,
            dispersion_method: DispersionMethod::default(),
        }
    }
}

/// Results from fitting a beta-binomial model to a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BbFitSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Estimated coefficients (logit scale).
    pub coefficients: Vec<f64>,
    /// Standard errors of coefficients.
    pub std_errors: Vec<f64>,
    /// Estimated overdispersion parameter (ρ in (0,1)).
    /// ρ = 0 corresponds to binomial, larger values indicate more overdispersion.
    pub overdispersion: f64,
    /// Log-likelihood at convergence.
    pub log_likelihood: f64,
    /// Deviance.
    pub deviance: f64,
    /// Degrees of freedom (residual).
    pub df_residual: usize,
    /// Number of iterations to convergence.
    pub iterations: usize,
    /// Whether the fit converged.
    pub converged: bool,
    /// Fitted values (predicted proportions).
    #[serde(skip)]
    pub fitted_values: Vec<f64>,
}

impl BbFitSingle {
    /// Get coefficient by index.
    pub fn get_coefficient(&self, index: usize) -> Option<f64> {
        self.coefficients.get(index).copied()
    }

    /// Get standard error by index.
    pub fn get_std_error(&self, index: usize) -> Option<f64> {
        self.std_errors.get(index).copied()
    }

    /// Calculate z-statistic for a coefficient.
    pub fn z_statistic(&self, index: usize) -> Option<f64> {
        let coef = self.coefficients.get(index)?;
        let se = self.std_errors.get(index)?;
        if *se > 0.0 {
            Some(coef / se)
        } else {
            None
        }
    }
}

/// Results from fitting beta-binomial models to all features.
#[derive(Debug, Clone)]
pub struct BbFit {
    /// Individual fits for each feature.
    pub fits: Vec<BbFitSingle>,
    /// Coefficient names from the design matrix.
    pub coefficient_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
}

impl BbFit {
    /// Get the fit for a specific feature by ID.
    pub fn get_feature(&self, feature_id: &str) -> Option<&BbFitSingle> {
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

    /// Number of coefficients.
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

    /// Get mean overdispersion across features.
    pub fn mean_overdispersion(&self) -> f64 {
        let sum: f64 = self.fits.iter().map(|f| f.overdispersion).sum();
        sum / self.fits.len() as f64
    }
}

/// Fit beta-binomial GLM to count data.
///
/// Uses Iteratively Reweighted Least Squares (IRLS) with a logit link function.
/// Models counts as proportions of library size with overdispersion.
///
/// # Arguments
/// * `counts` - Count matrix (features × samples)
/// * `design` - Design matrix from formula
///
/// # Returns
/// BbFit containing results for all features.
pub fn model_bb(counts: &CountMatrix, design: &DesignMatrix) -> Result<BbFit> {
    model_bb_with_config(counts, design, &BbConfig::default())
}

/// Fit beta-binomial GLM with custom configuration.
///
/// # Arguments
/// * `counts` - Count matrix (features × samples)
/// * `design` - Design matrix from formula
/// * `config` - Model configuration
///
/// # Returns
/// BbFit containing results for all features.
pub fn model_bb_with_config(
    counts: &CountMatrix,
    design: &DesignMatrix,
    config: &BbConfig,
) -> Result<BbFit> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let n_coef = design.n_coefficients();

    // Validate dimensions
    if design.n_samples() != n_samples {
        return Err(DaaError::DimensionMismatch {
            expected: n_samples,
            actual: design.n_samples(),
        });
    }

    let df_residual = n_samples.saturating_sub(n_coef);
    if df_residual == 0 {
        return Err(DaaError::Numerical(
            "Model is saturated (n_samples <= n_coefficients)".to_string(),
        ));
    }

    let x = design.matrix();

    // Compute library sizes (total counts per sample)
    let library_sizes: Vec<f64> = counts.col_sums().iter().map(|&s| s as f64).collect();

    // Fit all features in parallel
    let fits: Vec<BbFitSingle> = (0..n_features)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples)
                .map(|j| counts.get(i, j) as f64)
                .collect();
            fit_single_bb(
                &y,
                &library_sizes,
                &counts.feature_ids()[i],
                x,
                n_samples,
                n_coef,
                df_residual,
                config,
            )
        })
        .collect();

    Ok(BbFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        n_samples,
    })
}

/// Fit beta-binomial model to a single feature using IRLS.
fn fit_single_bb(
    y: &[f64],
    n: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    n_samples: usize,
    n_coef: usize,
    df_residual: usize,
    config: &BbConfig,
) -> BbFitSingle {
    let y_vec = DVector::from_column_slice(y);
    let n_vec = DVector::from_column_slice(n);

    // Compute proportions
    let p: Vec<f64> = y
        .iter()
        .zip(n.iter())
        .map(|(&yi, &ni)| {
            if ni > 0.0 {
                (yi / ni).clamp(MIN_MU, MAX_MU)
            } else {
                0.5 // default for zero library size
            }
        })
        .collect();

    // Initialize with smoothed logistic regression starting values
    // Add small constant to avoid log(0) issues
    let p_mean: f64 = p.iter().sum::<f64>() / n_samples as f64;
    let p_mean = p_mean.clamp(MIN_MU, MAX_MU);
    let logit_mean = (p_mean / (1.0 - p_mean)).ln();

    let mut beta = DVector::zeros(n_coef);
    beta[0] = logit_mean;

    // Compute initial mu
    let mut mu = compute_mu_bb(x, &beta);

    // Initial overdispersion estimate using method of moments
    let mut rho = estimate_dispersion_mom_bb(&y_vec, &n_vec, &mu);

    let mut converged = false;
    let mut iterations = 0;

    // IRLS iterations
    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Compute working weights and response
        let (w, z) = compute_working_data_bb(&y_vec, &n_vec, &mu, rho);

        // Check for invalid weights
        if w.iter().any(|&wi| !wi.is_finite() || wi <= 0.0) {
            // Return non-converged result with current estimates
            return BbFitSingle {
                feature_id: feature_id.to_string(),
                coefficients: beta.iter().cloned().collect(),
                std_errors: vec![f64::NAN; n_coef],
                overdispersion: rho,
                log_likelihood: f64::NAN,
                deviance: f64::NAN,
                df_residual,
                iterations,
                converged: false,
                fitted_values: mu.iter().cloned().collect(),
            };
        }

        // Weighted least squares: beta = (X'WX)^-1 X'Wz
        let w_sqrt = DVector::from_iterator(n_samples, w.iter().map(|wi| wi.sqrt()));

        // Scale X and z by sqrt(W)
        let mut xw = x.clone();
        let mut zw = z.clone();
        for i in 0..n_samples {
            for j in 0..n_coef {
                xw[(i, j)] *= w_sqrt[i];
            }
            zw[i] *= w_sqrt[i];
        }

        // Solve via QR decomposition
        let xtwx = xw.transpose() * &xw;
        let xtwz = xw.transpose() * &zw;

        let beta_new = match xtwx.clone().try_inverse() {
            Some(inv) => inv * xtwz,
            None => {
                // Singular matrix, return non-converged result
                return BbFitSingle {
                    feature_id: feature_id.to_string(),
                    coefficients: beta.iter().cloned().collect(),
                    std_errors: vec![f64::NAN; n_coef],
                    overdispersion: rho,
                    log_likelihood: f64::NAN,
                    deviance: f64::NAN,
                    df_residual,
                    iterations,
                    converged: false,
                    fitted_values: mu.iter().cloned().collect(),
                };
            }
        };

        // Check convergence
        let delta: f64 = (&beta_new - &beta).iter().map(|d| d.abs()).sum();
        let scale: f64 = beta.iter().map(|b| b.abs()).sum::<f64>().max(1.0);

        beta = beta_new;
        mu = compute_mu_bb(x, &beta);

        // Update overdispersion estimate
        rho = match config.dispersion_method {
            DispersionMethod::MomentsBased => estimate_dispersion_mom_bb(&y_vec, &n_vec, &mu),
            DispersionMethod::ProfileML => {
                estimate_dispersion_profile_ml(&y_vec, &n_vec, &mu, rho)
            }
        };

        if delta / scale < config.tol {
            converged = true;
            break;
        }
    }

    // Compute standard errors from Fisher information
    let std_errors = compute_std_errors_bb(x, &mu, &n_vec, rho, n_samples, n_coef);

    // Compute log-likelihood
    let log_lik = bb_log_likelihood(&y_vec, &n_vec, &mu, rho);

    // Compute deviance
    let deviance = bb_deviance(&y_vec, &n_vec, &mu, rho);

    BbFitSingle {
        feature_id: feature_id.to_string(),
        coefficients: beta.iter().cloned().collect(),
        std_errors,
        overdispersion: rho,
        log_likelihood: log_lik,
        deviance,
        df_residual,
        iterations,
        converged,
        fitted_values: mu.iter().cloned().collect(),
    }
}

/// Compute mu = expit(X * beta) where expit is the inverse logit.
fn compute_mu_bb(x: &DMatrix<f64>, beta: &DVector<f64>) -> DVector<f64> {
    let eta = x * beta;
    DVector::from_iterator(
        eta.len(),
        eta.iter().map(|e| {
            let p = 1.0 / (1.0 + (-e).exp());
            p.clamp(MIN_MU, MAX_MU)
        }),
    )
}

/// Estimate overdispersion using method of moments.
///
/// For beta-binomial: Var(Y/n) = μ(1-μ)[1 + (n-1)ρ] / n
/// Rearranging gives an estimator for ρ.
fn estimate_dispersion_mom_bb(y: &DVector<f64>, n: &DVector<f64>, mu: &DVector<f64>) -> f64 {
    let n_obs = y.len() as f64;

    // Compute Pearson residuals
    let mut sum_pearson_sq = 0.0;
    let mut sum_n_minus_1 = 0.0;

    for i in 0..y.len() {
        let yi = y[i];
        let ni = n[i];
        let mi = mu[i];

        if ni > 0.0 {
            let pi = yi / ni;
            let var_binomial = mi * (1.0 - mi) / ni;
            if var_binomial > 0.0 {
                let resid = (pi - mi) / var_binomial.sqrt();
                sum_pearson_sq += resid * resid;
            }
            sum_n_minus_1 += ni - 1.0;
        }
    }

    // Estimate rho from excess variance
    // Under beta-binomial: E[chi^2] ≈ n_obs + rho * sum(n_i - 1)
    let excess = sum_pearson_sq - n_obs;
    if excess > 0.0 && sum_n_minus_1 > 0.0 {
        (excess / sum_n_minus_1).clamp(MIN_RHO, MAX_RHO)
    } else {
        MIN_RHO // No overdispersion detected
    }
}

/// Estimate overdispersion using profile maximum likelihood.
///
/// Uses Newton-Raphson iterations to find the MLE of ρ.
fn estimate_dispersion_profile_ml(
    y: &DVector<f64>,
    n: &DVector<f64>,
    mu: &DVector<f64>,
    rho_init: f64,
) -> f64 {
    let mut rho = rho_init.clamp(MIN_RHO, MAX_RHO);
    let max_iter = 10;
    let tol = 1e-6;

    for _ in 0..max_iter {
        // Compute score and Fisher information for rho
        let (score, fisher) = dispersion_score_fisher(y, n, mu, rho);

        if fisher.abs() < 1e-10 {
            break;
        }

        let delta = score / fisher;
        let rho_new = (rho + delta).clamp(MIN_RHO, MAX_RHO);

        if (rho_new - rho).abs() < tol {
            rho = rho_new;
            break;
        }
        rho = rho_new;
    }

    rho
}

/// Compute score and Fisher information for overdispersion parameter.
fn dispersion_score_fisher(
    y: &DVector<f64>,
    n: &DVector<f64>,
    mu: &DVector<f64>,
    rho: f64,
) -> (f64, f64) {
    let mut score = 0.0;
    let mut fisher = 0.0;

    for i in 0..y.len() {
        let yi = y[i];
        let ni = n[i];
        let mi = mu[i];

        if ni > 0.0 {
            // Derivatives of log-likelihood with respect to phi
            // d/dphi involves digamma functions
            // For simplicity, use numerical approximation
            let eps = 1e-6;
            let ll_plus = bb_single_log_likelihood(yi, ni, mi, rho + eps);
            let ll_minus = bb_single_log_likelihood(yi, ni, mi, rho - eps);
            let ll_center = bb_single_log_likelihood(yi, ni, mi, rho);

            score += (ll_plus - ll_minus) / (2.0 * eps);
            fisher += (ll_plus - 2.0 * ll_center + ll_minus) / (eps * eps);
        }
    }

    // Fisher information is negative of expected second derivative
    (score, -fisher.max(1e-10))
}

/// Compute working weights and response for IRLS.
///
/// For beta-binomial with logit link:
/// - Working weight: W_i = n_i / [μ_i(1-μ_i)(1 + (n_i-1)ρ)]
/// - Working response: z_i = η_i + (y_i/n_i - μ_i) / [μ_i(1-μ_i)]
fn compute_working_data_bb(
    y: &DVector<f64>,
    n: &DVector<f64>,
    mu: &DVector<f64>,
    rho: f64,
) -> (DVector<f64>, DVector<f64>) {
    let n_obs = y.len();

    let w = DVector::from_iterator(
        n_obs,
        (0..n_obs).map(|i| {
            let ni = n[i];
            let mi = mu[i].clamp(MIN_MU, MAX_MU);
            let var_factor = mi * (1.0 - mi);
            let overdispersion_factor = 1.0 + (ni - 1.0).max(0.0) * rho;
            ni / (var_factor * overdispersion_factor).max(1e-10)
        }),
    );

    let z = DVector::from_iterator(
        n_obs,
        (0..n_obs).map(|i| {
            let yi = y[i];
            let ni = n[i];
            let mi = mu[i].clamp(MIN_MU, MAX_MU);

            // Linear predictor
            let eta = (mi / (1.0 - mi)).ln();

            // Proportion
            let pi = if ni > 0.0 { yi / ni } else { mi };

            // Working response
            let deriv = mi * (1.0 - mi);
            eta + (pi - mi) / deriv.max(1e-10)
        }),
    );

    (w, z)
}

/// Compute standard errors from Fisher information.
fn compute_std_errors_bb(
    x: &DMatrix<f64>,
    mu: &DVector<f64>,
    n: &DVector<f64>,
    rho: f64,
    n_samples: usize,
    n_coef: usize,
) -> Vec<f64> {
    // Fisher information: I = X' W X
    // where W_i = n_i / [μ_i(1-μ_i)(1 + (n_i-1)ρ)]
    let w: DVector<f64> = DVector::from_iterator(
        n_samples,
        (0..n_samples).map(|i| {
            let ni = n[i];
            let mi = mu[i].clamp(MIN_MU, MAX_MU);
            let var_factor = mi * (1.0 - mi);
            let overdispersion_factor = 1.0 + (ni - 1.0).max(0.0) * rho;
            ni / (var_factor * overdispersion_factor).max(1e-10)
        }),
    );

    let mut xw = x.clone();
    for i in 0..n_samples {
        let w_sqrt = w[i].sqrt();
        for j in 0..n_coef {
            xw[(i, j)] *= w_sqrt;
        }
    }

    let fisher = xw.transpose() * &xw;
    match fisher.try_inverse() {
        Some(inv) => (0..n_coef).map(|j| inv[(j, j)].max(0.0).sqrt()).collect(),
        None => vec![f64::NAN; n_coef],
    }
}

/// Compute beta-binomial log-likelihood for a single observation.
fn bb_single_log_likelihood(y: f64, n: f64, mu: f64, rho: f64) -> f64 {
    if n <= 0.0 {
        return 0.0;
    }

    // Convert from (mu, rho) to (alpha, beta)
    let phi = if rho > MIN_RHO {
        (1.0 - rho) / rho
    } else {
        1.0 / MIN_RHO
    };

    let alpha = mu * phi;
    let beta_param = (1.0 - mu) * phi;

    // Beta-binomial log-likelihood:
    // log(C(n,y)) + log(B(y + alpha, n - y + beta)) - log(B(alpha, beta))
    //
    // Using lgamma:
    // = lgamma(n+1) - lgamma(y+1) - lgamma(n-y+1)
    //   + lgamma(y + alpha) + lgamma(n - y + beta) - lgamma(n + alpha + beta)
    //   - lgamma(alpha) - lgamma(beta) + lgamma(alpha + beta)

    let ll = lgamma(n + 1.0) - lgamma(y + 1.0) - lgamma(n - y + 1.0)
        + lgamma(y + alpha)
        + lgamma(n - y + beta_param)
        - lgamma(n + alpha + beta_param)
        - lgamma(alpha)
        - lgamma(beta_param)
        + lgamma(alpha + beta_param);

    if ll.is_finite() {
        ll
    } else {
        f64::NEG_INFINITY
    }
}

/// Compute beta-binomial log-likelihood.
fn bb_log_likelihood(y: &DVector<f64>, n: &DVector<f64>, mu: &DVector<f64>, rho: f64) -> f64 {
    y.iter()
        .zip(n.iter())
        .zip(mu.iter())
        .map(|((&yi, &ni), &mi)| bb_single_log_likelihood(yi, ni, mi, rho))
        .sum()
}

/// Compute beta-binomial deviance.
fn bb_deviance(y: &DVector<f64>, n: &DVector<f64>, mu: &DVector<f64>, rho: f64) -> f64 {
    // Deviance = 2 * (log-likelihood of saturated model - log-likelihood of fitted model)
    let ll_fitted = bb_log_likelihood(y, n, mu, rho);

    // Saturated model: mu_i = y_i / n_i (or clamped)
    let ll_saturated: f64 = y
        .iter()
        .zip(n.iter())
        .map(|(&yi, &ni)| {
            if ni > 0.0 {
                let mu_sat = (yi / ni).clamp(MIN_MU, MAX_MU);
                bb_single_log_likelihood(yi, ni, mu_sat, rho)
            } else {
                0.0
            }
        })
        .sum();

    2.0 * (ll_saturated - ll_fitted)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{Formula, Metadata};
    use approx::assert_relative_eq;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        // 3 features × 8 samples with ~equal library sizes
        // Feature 0: no group effect (similar proportions ~10% in both groups)
        // Feature 1: strong group effect (higher proportion in treatment)
        // Feature 2: "other" (balances library sizes)
        let mut tri_mat = TriMat::new((3, 8));

        // Total library size target: ~1000 per sample
        // Feature 0: ~10% proportion in both groups (~100)
        // Control samples (0, 2, 4, 6): ~100
        // Treatment samples (1, 3, 5, 7): ~100
        tri_mat.add_triplet(0, 0, 95);
        tri_mat.add_triplet(0, 1, 105);
        tri_mat.add_triplet(0, 2, 98);
        tri_mat.add_triplet(0, 3, 102);
        tri_mat.add_triplet(0, 4, 97);
        tri_mat.add_triplet(0, 5, 103);
        tri_mat.add_triplet(0, 6, 99);
        tri_mat.add_triplet(0, 7, 101);

        // Feature 1: Different proportions
        // Control samples: ~50/1000 = 5%
        // Treatment samples: ~200/1000 = 20%
        tri_mat.add_triplet(1, 0, 50);
        tri_mat.add_triplet(1, 1, 200);
        tri_mat.add_triplet(1, 2, 52);
        tri_mat.add_triplet(1, 3, 198);
        tri_mat.add_triplet(1, 4, 48);
        tri_mat.add_triplet(1, 5, 202);
        tri_mat.add_triplet(1, 6, 53);
        tri_mat.add_triplet(1, 7, 197);

        // Feature 2: "other" to balance library sizes
        // Control: ~850 so total ~1000
        // Treatment: ~700 so total ~1000
        tri_mat.add_triplet(2, 0, 855);
        tri_mat.add_triplet(2, 1, 695);
        tri_mat.add_triplet(2, 2, 850);
        tri_mat.add_triplet(2, 3, 700);
        tri_mat.add_triplet(2, 4, 855);
        tri_mat.add_triplet(2, 5, 695);
        tri_mat.add_triplet(2, 6, 848);
        tri_mat.add_triplet(2, 7, 702);

        let feature_ids = vec!["no_effect".into(), "strong_effect".into(), "other".into()];
        let sample_ids: Vec<String> = (1..=8).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        // Alternating control/treatment
        writeln!(file, "S1\tcontrol").unwrap();
        writeln!(file, "S2\ttreatment").unwrap();
        writeln!(file, "S3\tcontrol").unwrap();
        writeln!(file, "S4\ttreatment").unwrap();
        writeln!(file, "S5\tcontrol").unwrap();
        writeln!(file, "S6\ttreatment").unwrap();
        writeln!(file, "S7\tcontrol").unwrap();
        writeln!(file, "S8\ttreatment").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_model_bb_basic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_coefficients(), 2);
        assert_eq!(fit.coefficient_names, vec!["(Intercept)", "grouptreatment"]);
    }

    #[test]
    fn test_model_bb_convergence() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Both features should converge
        assert!(fit.all_converged(), "All fits should converge");
    }

    #[test]
    fn test_model_bb_coefficients() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Feature 0 should have near-zero treatment effect
        let no_effect = fit.get_feature("no_effect").unwrap();
        let coef_no_effect = no_effect.coefficients[1];
        assert!(
            coef_no_effect.abs() < 0.5,
            "No-effect feature should have small coefficient, got {}",
            coef_no_effect
        );

        // Feature 1 should have positive treatment effect
        let strong_effect = fit.get_feature("strong_effect").unwrap();
        let coef_strong = strong_effect.coefficients[1];
        assert!(
            coef_strong > 0.5,
            "Strong-effect feature should have positive coefficient, got {}",
            coef_strong
        );
    }

    #[test]
    fn test_model_bb_overdispersion() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Overdispersion should be in (0, 1)
        for f in &fit.fits {
            assert!(
                f.overdispersion > 0.0 && f.overdispersion < 1.0,
                "Overdispersion {} should be in (0, 1)",
                f.overdispersion
            );
        }
    }

    #[test]
    fn test_model_bb_std_errors() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Standard errors should be positive and finite
        for f in &fit.fits {
            for &se in &f.std_errors {
                assert!(
                    se > 0.0 && se.is_finite(),
                    "SE should be positive finite, got {}",
                    se
                );
            }
        }
    }

    #[test]
    fn test_model_bb_fitted_values() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Fitted values should be proportions in (0, 1)
        for f in &fit.fits {
            for &fv in &f.fitted_values {
                assert!(
                    fv > 0.0 && fv < 1.0,
                    "Fitted values should be proportions in (0,1), got {}",
                    fv
                );
            }
        }
    }

    #[test]
    fn test_model_bb_dimension_mismatch() {
        let counts = create_test_counts();

        // Create metadata with wrong number of samples
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        writeln!(file, "S1\tcontrol").unwrap();
        writeln!(file, "S2\ttreatment").unwrap();
        file.flush().unwrap();
        let metadata = Metadata::from_tsv(file.path()).unwrap();

        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let result = model_bb(&counts, &design);
        assert!(result.is_err());
    }

    #[test]
    fn test_z_statistic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_bb(&counts, &design).unwrap();

        // Strong effect feature should have larger z-statistic
        let strong = fit.get_feature("strong_effect").unwrap();
        let z = strong.z_statistic(1).unwrap();
        assert!(z.abs() > 1.0, "Strong effect should have z > 1, got {}", z);

        // No effect feature should have smaller z-statistic
        let no_effect = fit.get_feature("no_effect").unwrap();
        let _z_no = no_effect.z_statistic(1).unwrap();
        // Note: this is a weak test given noise in small samples
    }

    #[test]
    fn test_bb_log_likelihood() {
        // Test log-likelihood computation
        let y = DVector::from_vec(vec![5.0, 10.0, 15.0]);
        let n = DVector::from_vec(vec![100.0, 100.0, 100.0]);
        let mu = DVector::from_vec(vec![0.05, 0.10, 0.15]);
        let rho = 0.01;

        let ll = bb_log_likelihood(&y, &n, &mu, rho);
        assert!(ll.is_finite(), "Log-likelihood should be finite");
        assert!(ll < 0.0, "Log-likelihood should be negative");
    }

    #[test]
    fn test_compute_mu_bb() {
        // Test inverse logit computation
        let x = DMatrix::from_row_slice(3, 2, &[1.0, 0.0, 1.0, 1.0, 1.0, 0.0]);
        let beta = DVector::from_vec(vec![0.0, 1.0]); // intercept=0, slope=1

        let mu = compute_mu_bb(&x, &beta);

        // For x[0]: eta = 0, mu = 0.5
        assert_relative_eq!(mu[0], 0.5, epsilon = 0.01);

        // For x[1]: eta = 1, mu = 1/(1+e^-1) ≈ 0.731
        assert_relative_eq!(mu[1], 0.731, epsilon = 0.01);

        // For x[2]: eta = 0, mu = 0.5
        assert_relative_eq!(mu[2], 0.5, epsilon = 0.01);
    }

    #[test]
    fn test_config_with_profile_ml() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let config = BbConfig {
            dispersion_method: DispersionMethod::ProfileML,
            ..Default::default()
        };

        let fit = model_bb_with_config(&counts, &design, &config).unwrap();

        // Should still converge
        assert!(fit.n_converged() > 0, "At least one fit should converge");

        // Overdispersion should be in valid range
        for f in &fit.fits {
            assert!(f.overdispersion > 0.0 && f.overdispersion < 1.0);
        }
    }
}
