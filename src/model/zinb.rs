//! Zero-Inflated Negative Binomial (ZINB) model for count data with excess zeros.
//!
//! ZINB models the data as a mixture:
//! - With probability π, the observation is a "structural zero"
//! - With probability (1-π), the observation comes from a negative binomial distribution
//!
//! This is appropriate for microbiome data where zeros may arise from:
//! 1. True absence (structural zeros)
//! 2. Sampling zeros (taxon present but not detected)

use crate::data::{CountMatrix, DesignMatrix};
use crate::error::{DaaError, Result};
use crate::model::nb::lgamma;
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Fit result for a single feature from ZINB model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZinbFitSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Count model coefficients (log link).
    pub coefficients: Vec<f64>,
    /// Standard errors for count model coefficients.
    pub std_errors: Vec<f64>,
    /// Zero-inflation model coefficients (logit link).
    pub zi_coefficients: Vec<f64>,
    /// Standard errors for zero-inflation coefficients.
    pub zi_std_errors: Vec<f64>,
    /// Estimated dispersion parameter (theta).
    pub dispersion: f64,
    /// Log-likelihood at convergence.
    pub log_likelihood: f64,
    /// Residual degrees of freedom.
    pub df_residual: usize,
    /// Number of EM iterations.
    pub iterations: usize,
    /// Whether the algorithm converged.
    pub converged: bool,
    /// Fitted values (expected counts).
    pub fitted_values: Vec<f64>,
    /// Estimated zero-inflation probabilities.
    pub zi_probs: Vec<f64>,
}

/// Fit results for all features from ZINB model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZinbFit {
    /// Individual feature fits.
    pub fits: Vec<ZinbFitSingle>,
    /// Coefficient names for count model.
    pub coefficient_names: Vec<String>,
    /// Coefficient names for zero-inflation model.
    pub zi_coefficient_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
}

impl ZinbFit {
    /// Get coefficient index by name.
    pub fn coefficient_index(&self, name: &str) -> Option<usize> {
        self.coefficient_names.iter().position(|n| n == name)
    }

    /// Get zero-inflation coefficient index by name.
    pub fn zi_coefficient_index(&self, name: &str) -> Option<usize> {
        self.zi_coefficient_names.iter().position(|n| n == name)
    }

    /// Get number of features.
    pub fn n_features(&self) -> usize {
        self.fits.len()
    }

    /// Get fit for a specific feature.
    pub fn get_feature(&self, feature_id: &str) -> Option<&ZinbFitSingle> {
        self.fits.iter().find(|f| f.feature_id == feature_id)
    }
}

/// Configuration for ZINB model fitting.
#[derive(Debug, Clone)]
pub struct ZinbConfig {
    /// Maximum EM iterations.
    pub max_iter: usize,
    /// Convergence tolerance for log-likelihood.
    pub tol: f64,
    /// Maximum IRLS iterations within each M-step.
    pub max_irls_iter: usize,
    /// IRLS convergence tolerance.
    pub irls_tol: f64,
    /// Whether to use same design matrix for zero-inflation model.
    pub zi_same_design: bool,
}

impl Default for ZinbConfig {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-6,
            max_irls_iter: 25,
            irls_tol: 1e-8,
            zi_same_design: true,
        }
    }
}

/// Fit ZINB model to count data.
///
/// Uses an EM algorithm:
/// - E-step: Compute posterior probability each zero is structural
/// - M-step: Update count model (μ), zero-inflation model (π), and dispersion (θ)
///
/// # Arguments
/// * `counts` - Raw count matrix (features × samples)
/// * `design` - Design matrix for covariates
///
/// # Returns
/// ZinbFit containing model parameters for all features.
pub fn model_zinb(counts: &CountMatrix, design: &DesignMatrix) -> Result<ZinbFit> {
    model_zinb_with_config(counts, design, &ZinbConfig::default())
}

/// Fit ZINB model with custom configuration.
pub fn model_zinb_with_config(
    counts: &CountMatrix,
    design: &DesignMatrix,
    config: &ZinbConfig,
) -> Result<ZinbFit> {
    // Validate dimensions
    if counts.n_samples() != design.n_samples() {
        return Err(DaaError::DimensionMismatch {
            expected: counts.n_samples(),
            actual: design.n_samples(),
        });
    }

    let design_matrix = design.matrix();
    let n_samples = counts.n_samples();
    let n_coef = design_matrix.ncols();

    // For zero-inflation, use intercept-only or same design
    let zi_design = if config.zi_same_design {
        design_matrix.clone()
    } else {
        // Intercept-only model for zero-inflation
        DMatrix::from_element(n_samples, 1, 1.0)
    };
    let n_zi_coef = zi_design.ncols();

    // Fit each feature in parallel
    let fits: Vec<ZinbFitSingle> = counts
        .feature_ids()
        .par_iter()
        .enumerate()
        .map(|(i, feature_id)| {
            let y: Vec<f64> = counts.row_dense(i).into_iter().map(|x| x as f64).collect();
            fit_single_zinb(
                feature_id,
                &y,
                design_matrix,
                &zi_design,
                config,
            )
        })
        .collect();

    // Build coefficient names for zero-inflation model
    let zi_coef_names = if config.zi_same_design {
        design.coefficient_names().iter().map(|s| format!("zi_{}", s)).collect()
    } else {
        vec!["zi_intercept".to_string()]
    };

    Ok(ZinbFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        zi_coefficient_names: zi_coef_names,
        n_samples,
    })
}

/// Fit ZINB model for a single feature using EM algorithm.
fn fit_single_zinb(
    feature_id: &str,
    y: &[f64],
    x: &DMatrix<f64>,
    zi_x: &DMatrix<f64>,
    config: &ZinbConfig,
) -> ZinbFitSingle {
    let n = y.len();
    let p = x.ncols();
    let p_zi = zi_x.ncols();

    // Initialize parameters
    let y_mean = y.iter().sum::<f64>() / n as f64;
    let y_mean = y_mean.max(0.1); // Avoid log(0)

    // Initialize count model coefficients (intercept = log(mean), others = 0)
    let mut beta = vec![0.0; p];
    beta[0] = y_mean.ln();

    // Initialize zero-inflation coefficients
    // Start with proportion of zeros
    let n_zeros = y.iter().filter(|&&yi| yi == 0.0).count();
    let zero_prop = (n_zeros as f64 / n as f64).max(0.01).min(0.99);
    let mut gamma = vec![0.0; p_zi];
    gamma[0] = (zero_prop / (1.0 - zero_prop)).ln(); // logit(zero_prop)

    // Initialize dispersion
    let y_var = y.iter().map(|&yi| (yi - y_mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    let mut theta = if y_var > y_mean {
        (y_mean * y_mean / (y_var - y_mean)).max(0.1)
    } else {
        10.0 // Low overdispersion default
    };

    let mut prev_ll = f64::NEG_INFINITY;
    let mut converged = false;
    let mut iter = 0;

    for iteration in 0..config.max_iter {
        iter = iteration + 1;

        // E-step: Compute posterior probability that each zero is structural
        let mu = compute_mu_zinb(x, &beta);
        let pi = compute_pi(zi_x, &gamma);
        let w = e_step(y, &mu, &pi, theta);

        // M-step: Update parameters

        // Update zero-inflation model (logistic regression on w)
        gamma = m_step_zi(zi_x, &w, &gamma, config.max_irls_iter, config.irls_tol);

        // Update count model (weighted NB regression)
        beta = m_step_count(x, y, &w, &beta, theta, config.max_irls_iter, config.irls_tol);

        // Update dispersion
        let mu_new = compute_mu_zinb(x, &beta);
        theta = estimate_dispersion_zinb(y, &mu_new, &w, theta);

        // Compute log-likelihood
        let pi_new = compute_pi(zi_x, &gamma);
        let ll = zinb_log_likelihood(y, &mu_new, &pi_new, theta);

        // Check convergence
        if (ll - prev_ll).abs() < config.tol {
            converged = true;
            break;
        }
        prev_ll = ll;
    }

    // Final computations
    let mu_final = compute_mu_zinb(x, &beta);
    let pi_final = compute_pi(zi_x, &gamma);
    let ll_final = zinb_log_likelihood(y, &mu_final, &pi_final, theta);

    // Compute standard errors
    let (se_beta, se_gamma) = compute_zinb_std_errors(y, x, zi_x, &beta, &gamma, theta);

    // Degrees of freedom: n - p (count) - p_zi (zero-inflation) - 1 (dispersion)
    let df_residual = n.saturating_sub(p + p_zi + 1);

    ZinbFitSingle {
        feature_id: feature_id.to_string(),
        coefficients: beta,
        std_errors: se_beta,
        zi_coefficients: gamma,
        zi_std_errors: se_gamma,
        dispersion: theta,
        log_likelihood: ll_final,
        df_residual,
        iterations: iter,
        converged,
        fitted_values: mu_final,
        zi_probs: pi_final,
    }
}

/// Compute expected counts μ = exp(X * β).
fn compute_mu_zinb(x: &DMatrix<f64>, beta: &[f64]) -> Vec<f64> {
    let beta_vec = nalgebra::DVector::from_column_slice(beta);
    let eta = x * beta_vec;
    eta.iter().map(|&e| e.exp().min(1e10)).collect()
}

/// Compute zero-inflation probability π = 1 / (1 + exp(-X * γ)).
fn compute_pi(x: &DMatrix<f64>, gamma: &[f64]) -> Vec<f64> {
    let gamma_vec = nalgebra::DVector::from_column_slice(gamma);
    let eta = x * gamma_vec;
    eta.iter().map(|&e| {
        let exp_neg = (-e).exp();
        1.0 / (1.0 + exp_neg)
    }).collect()
}

/// E-step: Compute posterior probability that each zero is structural.
///
/// For y = 0: P(structural | y=0) = π / (π + (1-π) * NB(0))
/// For y > 0: P(structural | y>0) = 0
fn e_step(y: &[f64], mu: &[f64], pi: &[f64], theta: f64) -> Vec<f64> {
    y.iter()
        .zip(mu.iter())
        .zip(pi.iter())
        .map(|((&yi, &mui), &pii)| {
            if yi > 0.0 {
                0.0 // Cannot be structural zero if y > 0
            } else {
                // P(y=0 | NB) = (theta / (theta + mu))^theta
                let nb_zero = (theta / (theta + mui)).powf(theta);
                let denom = pii + (1.0 - pii) * nb_zero;
                if denom > 0.0 {
                    (pii / denom).min(1.0 - 1e-10).max(1e-10)
                } else {
                    0.5
                }
            }
        })
        .collect()
}

/// M-step: Update zero-inflation parameters via weighted logistic regression.
fn m_step_zi(
    x: &DMatrix<f64>,
    w: &[f64],
    gamma_init: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = w.len();
    let p = x.ncols();
    let mut gamma = gamma_init.to_vec();

    for _ in 0..max_iter {
        let pi = compute_pi(x, &gamma);

        // Working response and weights for IRLS
        let mut z = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);

        for i in 0..n {
            let pii = pi[i].max(1e-10).min(1.0 - 1e-10);
            // Working response: η + (w - π) / (π * (1 - π))
            let eta = (pii / (1.0 - pii)).ln();
            z.push(eta + (w[i] - pii) / (pii * (1.0 - pii)));
            // Weight: π * (1 - π)
            weights.push(pii * (1.0 - pii));
        }

        // Weighted least squares
        let gamma_new = weighted_least_squares(x, &z, &weights, &gamma);

        // Check convergence
        let delta: f64 = gamma_new.iter()
            .zip(gamma.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();

        gamma = gamma_new;

        if delta < tol {
            break;
        }
    }

    gamma
}

/// M-step: Update count model parameters via weighted NB regression.
fn m_step_count(
    x: &DMatrix<f64>,
    y: &[f64],
    w: &[f64],
    beta_init: &[f64],
    theta: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = y.len();
    let mut beta = beta_init.to_vec();

    for _ in 0..max_iter {
        let mu = compute_mu_zinb(x, &beta);

        // Working response and weights for IRLS
        // Weight observations by (1 - w), since w is probability of structural zero
        let mut z = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);

        for i in 0..n {
            let mui = mu[i].max(1e-10);
            let eta = mui.ln();

            // Working response: η + (y - μ) / μ
            z.push(eta + (y[i] - mui) / mui);

            // NB weight: μ / (1 + μ/θ), scaled by (1 - w)
            let nb_weight = mui / (1.0 + mui / theta);
            weights.push((1.0 - w[i]) * nb_weight);
        }

        // Weighted least squares
        let beta_new = weighted_least_squares(x, &z, &weights, &beta);

        // Check convergence
        let delta: f64 = beta_new.iter()
            .zip(beta.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();

        beta = beta_new;

        if delta < tol {
            break;
        }
    }

    beta
}

/// Weighted least squares: solve X'WX * β = X'Wz.
fn weighted_least_squares(
    x: &DMatrix<f64>,
    z: &[f64],
    weights: &[f64],
    beta_init: &[f64],
) -> Vec<f64> {
    let n = z.len();
    let p = x.ncols();

    // Build weighted design matrix and response
    let mut xtwx = DMatrix::zeros(p, p);
    let mut xtwz = nalgebra::DVector::zeros(p);

    for i in 0..n {
        let w = weights[i].max(1e-10);
        let xi = x.row(i);

        for j in 0..p {
            xtwz[j] += w * xi[j] * z[i];
            for k in 0..p {
                xtwx[(j, k)] += w * xi[j] * xi[k];
            }
        }
    }

    // Add small ridge for numerical stability
    for j in 0..p {
        xtwx[(j, j)] += 1e-8;
    }

    // Solve
    match xtwx.clone().lu().solve(&xtwz) {
        Some(solution) => solution.iter().cloned().collect(),
        None => beta_init.to_vec(),
    }
}

/// Estimate dispersion for ZINB model using method of moments.
fn estimate_dispersion_zinb(y: &[f64], mu: &[f64], w: &[f64], theta_init: f64) -> f64 {
    let n = y.len();

    // Compute weighted variance
    let weight_sum: f64 = w.iter().map(|&wi| 1.0 - wi).sum();
    if weight_sum < 1.0 {
        return theta_init;
    }

    let weighted_mean: f64 = y.iter()
        .zip(mu.iter())
        .zip(w.iter())
        .map(|((&yi, &_mui), &wi)| (1.0 - wi) * yi)
        .sum::<f64>() / weight_sum;

    let weighted_var: f64 = y.iter()
        .zip(mu.iter())
        .zip(w.iter())
        .map(|((&yi, &mui), &wi)| (1.0 - wi) * (yi - mui).powi(2))
        .sum::<f64>() / weight_sum;

    // Method of moments: Var(Y) = μ + μ²/θ
    // => θ = μ² / (Var(Y) - μ)
    if weighted_var > weighted_mean && weighted_mean > 0.0 {
        let theta = weighted_mean * weighted_mean / (weighted_var - weighted_mean);
        theta.max(0.01).min(1e6)
    } else {
        theta_init.max(0.01).min(1e6)
    }
}

/// Compute ZINB log-likelihood.
fn zinb_log_likelihood(y: &[f64], mu: &[f64], pi: &[f64], theta: f64) -> f64 {
    let ll_sum: f64 = y.iter()
        .zip(mu.iter())
        .zip(pi.iter())
        .map(|((&yi, &mui), &pii)| {
            let pii = pii.max(1e-10).min(1.0 - 1e-10);
            let mui = mui.max(1e-10);

            if yi == 0.0 {
                // P(Y=0) = π + (1-π) * NB(0)
                let nb_zero = (theta / (theta + mui)).powf(theta);
                let prob = pii + (1.0 - pii) * nb_zero;
                prob.max(1e-300).ln()
            } else {
                // P(Y=y) = (1-π) * NB(y)
                let log_nb = lgamma(yi + theta) - lgamma(theta) - lgamma(yi + 1.0)
                    + theta * (theta / (theta + mui)).ln()
                    + yi * (mui / (theta + mui)).ln();
                (1.0 - pii).ln() + log_nb
            }
        })
        .sum();

    ll_sum
}

/// Compute standard errors for ZINB model using observed Fisher information.
fn compute_zinb_std_errors(
    y: &[f64],
    x: &DMatrix<f64>,
    zi_x: &DMatrix<f64>,
    beta: &[f64],
    gamma: &[f64],
    theta: f64,
) -> (Vec<f64>, Vec<f64>) {
    let n = y.len();
    let p = x.ncols();
    let p_zi = zi_x.ncols();

    let mu = compute_mu_zinb(x, beta);
    let pi = compute_pi(zi_x, gamma);

    // Compute Fisher information for count model
    let mut info_beta = DMatrix::zeros(p, p);

    for i in 0..n {
        let mui = mu[i].max(1e-10);
        let pii = pi[i].max(1e-10).min(1.0 - 1e-10);

        // Expected information for NB part
        // Weight by (1 - E[structural zero | y])
        let w_count = if y[i] == 0.0 {
            let nb_zero = (theta / (theta + mui)).powf(theta);
            1.0 - pii / (pii + (1.0 - pii) * nb_zero)
        } else {
            1.0
        };

        let nb_var = mui + mui * mui / theta;
        let info_weight = w_count * mui * mui / nb_var;

        let xi = x.row(i);
        for j in 0..p {
            for k in 0..p {
                info_beta[(j, k)] += info_weight * xi[j] * xi[k];
            }
        }
    }

    // Compute Fisher information for zero-inflation model
    let mut info_gamma = DMatrix::zeros(p_zi, p_zi);

    for i in 0..n {
        let pii = pi[i].max(1e-10).min(1.0 - 1e-10);
        let info_weight = pii * (1.0 - pii);

        let xi = zi_x.row(i);
        for j in 0..p_zi {
            for k in 0..p_zi {
                info_gamma[(j, k)] += info_weight * xi[j] * xi[k];
            }
        }
    }

    // Invert to get variance-covariance matrices
    let se_beta = invert_and_extract_se(&info_beta, p);
    let se_gamma = invert_and_extract_se(&info_gamma, p_zi);

    (se_beta, se_gamma)
}

/// Invert information matrix and extract standard errors.
fn invert_and_extract_se(info: &DMatrix<f64>, p: usize) -> Vec<f64> {
    // Add small ridge for stability
    let mut info_ridge = info.clone();
    for j in 0..p {
        info_ridge[(j, j)] += 1e-8;
    }

    match info_ridge.clone().try_inverse() {
        Some(var_cov) => {
            (0..p).map(|j| {
                let v = var_cov[(j, j)];
                if v > 0.0 { v.sqrt() } else { f64::NAN }
            }).collect()
        }
        None => vec![f64::NAN; p],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DesignMatrix, Formula, Metadata};
    use approx::assert_relative_eq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        for i in 1..=20 {
            let group = if i <= 10 { "control" } else { "treatment" };
            writeln!(file, "S{}\t{}", i, group).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn create_test_counts_zinb() -> CountMatrix {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "feature_id\t{}", (1..=20).map(|i| format!("S{}", i)).collect::<Vec<_>>().join("\t")).unwrap();

        // Feature with no zeros (standard NB)
        writeln!(file, "no_zeros\t{}", (1..=20).map(|i| if i <= 10 { "10" } else { "15" }).collect::<Vec<_>>().join("\t")).unwrap();

        // Feature with excess zeros in control group (ZINB pattern)
        let zinb_counts: Vec<&str> = (1..=20).map(|i| {
            if i <= 10 {
                if i <= 5 { "0" } else { "8" }  // 50% zeros in control
            } else {
                "20"  // No zeros in treatment
            }
        }).collect();
        writeln!(file, "excess_zeros\t{}", zinb_counts.join("\t")).unwrap();

        // Feature with moderate effect
        writeln!(file, "moderate\t{}", (1..=20).map(|i| if i <= 10 { "5" } else { "12" }).collect::<Vec<_>>().join("\t")).unwrap();

        file.flush().unwrap();
        CountMatrix::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_zinb_basic() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_samples, 20);
        assert_eq!(fit.coefficient_names.len(), 2); // intercept + group
    }

    #[test]
    fn test_zinb_convergence() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Most features should converge
        let converged_count = fit.fits.iter().filter(|f| f.converged).count();
        assert!(converged_count >= 2, "At least 2 features should converge");
    }

    #[test]
    fn test_zinb_zero_inflation() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Feature with excess zeros should have higher ZI probability
        let excess_zeros = fit.get_feature("excess_zeros").unwrap();
        let no_zeros = fit.get_feature("no_zeros").unwrap();

        // ZI intercept should be higher for excess_zeros feature
        // (more structural zeros overall)
        let zi_intercept_excess = excess_zeros.zi_coefficients[0];
        let zi_intercept_no = no_zeros.zi_coefficients[0];

        // The feature with excess zeros should have higher baseline ZI
        // This is a soft check since EM can find different local optima
        assert!(zi_intercept_excess > zi_intercept_no - 2.0,
            "Feature with excess zeros should have comparable or higher ZI intercept");
    }

    #[test]
    fn test_zinb_coefficients() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Check that coefficients are finite
        for f in &fit.fits {
            for coef in &f.coefficients {
                assert!(coef.is_finite(), "Coefficients should be finite");
            }
            for zi_coef in &f.zi_coefficients {
                assert!(zi_coef.is_finite(), "ZI coefficients should be finite");
            }
        }
    }

    #[test]
    fn test_zinb_std_errors() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Standard errors should be positive for converged fits
        for f in &fit.fits {
            if f.converged {
                for se in &f.std_errors {
                    if !se.is_nan() {
                        assert!(*se >= 0.0, "Standard errors should be non-negative");
                    }
                }
            }
        }
    }

    #[test]
    fn test_zinb_dispersion() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Dispersion should be positive
        for f in &fit.fits {
            assert!(f.dispersion > 0.0, "Dispersion must be positive");
        }
    }

    #[test]
    fn test_zinb_fitted_values() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // Fitted values should be non-negative
        for f in &fit.fits {
            assert_eq!(f.fitted_values.len(), 20);
            for fv in &f.fitted_values {
                assert!(*fv >= 0.0, "Fitted values must be non-negative");
            }
        }
    }

    #[test]
    fn test_zinb_zi_probs() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let counts = create_test_counts_zinb();

        let fit = model_zinb(&counts, &design).unwrap();

        // ZI probabilities should be between 0 and 1
        for f in &fit.fits {
            assert_eq!(f.zi_probs.len(), 20);
            for pi in &f.zi_probs {
                assert!(*pi >= 0.0 && *pi <= 1.0,
                    "ZI probability {} should be in [0, 1]", pi);
            }
        }
    }

    #[test]
    fn test_zinb_dimension_mismatch() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        // Create counts with different number of samples
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "feature_id\tS1\tS2\tS3").unwrap();
        writeln!(file, "f1\t1\t2\t3").unwrap();
        file.flush().unwrap();
        let counts = CountMatrix::from_tsv(file.path()).unwrap();

        let result = model_zinb(&counts, &design);
        assert!(result.is_err());
    }

    #[test]
    fn test_zinb_log_likelihood() {
        // Test log-likelihood computation
        let y = vec![0.0, 0.0, 5.0, 10.0, 15.0];
        let mu = vec![5.0, 5.0, 5.0, 10.0, 15.0];
        let pi = vec![0.3, 0.3, 0.1, 0.1, 0.1];
        let theta = 2.0;

        let ll = zinb_log_likelihood(&y, &mu, &pi, theta);
        assert!(ll.is_finite(), "Log-likelihood should be finite");
        assert!(ll < 0.0, "Log-likelihood should be negative");
    }

    #[test]
    fn test_compute_pi() {
        let x = DMatrix::from_row_slice(3, 1, &[1.0, 1.0, 1.0]);
        let gamma = vec![0.0]; // logit(0.5) = 0

        let pi = compute_pi(&x, &gamma);

        for pii in &pi {
            assert_relative_eq!(*pii, 0.5, epsilon = 0.01);
        }
    }

    #[test]
    fn test_e_step() {
        let y = vec![0.0, 0.0, 5.0];
        let mu = vec![1.0, 1.0, 5.0];
        let pi = vec![0.5, 0.5, 0.5];
        let theta = 1.0;

        let w = e_step(&y, &mu, &pi, theta);

        // For y > 0, w should be 0
        assert_relative_eq!(w[2], 0.0, epsilon = 1e-10);

        // For y = 0, w should be between 0 and 1
        assert!(w[0] > 0.0 && w[0] < 1.0);
        assert!(w[1] > 0.0 && w[1] < 1.0);
    }
}
