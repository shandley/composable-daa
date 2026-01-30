//! Hurdle model for zero-inflated count data.
//!
//! A hurdle model is a two-part model:
//! 1. **Binary component**: Logistic regression for P(Y > 0)
//! 2. **Count component**: Truncated negative binomial for Y | Y > 0
//!
//! Unlike ZINB which models zeros as a mixture (structural + sampling),
//! hurdle models treat ALL zeros as coming from the binary process.
//! This is appropriate when zeros represent true absence rather than
//! a mixture of absence and under-detection.
//!
//! # Mathematical Formulation
//!
//! P(Y = 0) = 1 - p
//! P(Y = k) = p * f_TNB(k; μ, θ) for k > 0
//!
//! where:
//! - p = probability of non-zero (from logistic regression)
//! - f_TNB = truncated negative binomial density
//!
//! # When to Use Hurdle vs ZINB
//!
//! - **Hurdle**: When zeros have clear biological meaning (true absence)
//! - **ZINB**: When zeros may be a mixture of absence and under-sampling

use crate::data::{CountMatrix, DesignMatrix};
use crate::error::{DaaError, Result};
use crate::model::nb::lgamma;
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Maximum iterations for IRLS convergence.
const MAX_ITER: usize = 25;

/// Convergence tolerance.
const TOL: f64 = 1e-8;

/// Minimum probability to avoid log(0).
const MIN_P: f64 = 1e-10;

/// Maximum probability.
const MAX_P: f64 = 1.0 - 1e-10;

/// Minimum mean for count model.
const MIN_MU: f64 = 1e-10;

/// Type of count model for the positive counts.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default, PartialEq)]
pub enum CountModel {
    /// Truncated Poisson (no overdispersion).
    TruncatedPoisson,
    /// Truncated Negative Binomial (with overdispersion, default).
    #[default]
    TruncatedNB,
}

/// Configuration for hurdle model fitting.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HurdleConfig {
    /// Maximum IRLS iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
    /// Type of count model for positive counts.
    pub count_model: CountModel,
}

impl Default for HurdleConfig {
    fn default() -> Self {
        Self {
            max_iter: MAX_ITER,
            tol: TOL,
            count_model: CountModel::default(),
        }
    }
}

/// Results from fitting a hurdle model to a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HurdleFitSingle {
    /// Feature identifier.
    pub feature_id: String,

    // Binary component (logistic regression for P(Y > 0))
    /// Binary model coefficients (logit scale).
    pub binary_coefficients: Vec<f64>,
    /// Standard errors for binary model coefficients.
    pub binary_std_errors: Vec<f64>,

    // Count component (truncated NB for Y | Y > 0)
    /// Count model coefficients (log scale).
    pub count_coefficients: Vec<f64>,
    /// Standard errors for count model coefficients.
    pub count_std_errors: Vec<f64>,
    /// Dispersion parameter for count model (theta).
    /// Variance = mu + mu^2/theta. Large theta -> Poisson.
    pub dispersion: f64,

    // Model fit statistics
    /// Combined log-likelihood.
    pub log_likelihood: f64,
    /// Deviance.
    pub deviance: f64,
    /// Degrees of freedom (residual).
    pub df_residual: usize,
    /// Number of iterations (max of binary and count).
    pub iterations: usize,
    /// Whether both components converged.
    pub converged: bool,

    // Fitted values
    /// Fitted values (expected counts = P(Y>0) * E[Y|Y>0]).
    #[serde(skip)]
    pub fitted_values: Vec<f64>,
    /// Fitted probability of non-zero.
    #[serde(skip)]
    pub fitted_probs: Vec<f64>,
    /// Fitted mean of positive counts E[Y|Y>0].
    #[serde(skip)]
    pub fitted_positive_mean: Vec<f64>,
}

impl HurdleFitSingle {
    /// Get binary coefficient by index.
    pub fn get_binary_coefficient(&self, index: usize) -> Option<f64> {
        self.binary_coefficients.get(index).copied()
    }

    /// Get binary standard error by index.
    pub fn get_binary_std_error(&self, index: usize) -> Option<f64> {
        self.binary_std_errors.get(index).copied()
    }

    /// Get count coefficient by index.
    pub fn get_count_coefficient(&self, index: usize) -> Option<f64> {
        self.count_coefficients.get(index).copied()
    }

    /// Get count standard error by index.
    pub fn get_count_std_error(&self, index: usize) -> Option<f64> {
        self.count_std_errors.get(index).copied()
    }

    /// Calculate z-statistic for binary coefficient.
    pub fn binary_z_statistic(&self, index: usize) -> Option<f64> {
        let coef = self.binary_coefficients.get(index)?;
        let se = self.binary_std_errors.get(index)?;
        if *se > 0.0 && se.is_finite() {
            Some(coef / se)
        } else {
            None
        }
    }

    /// Calculate z-statistic for count coefficient.
    pub fn count_z_statistic(&self, index: usize) -> Option<f64> {
        let coef = self.count_coefficients.get(index)?;
        let se = self.count_std_errors.get(index)?;
        if *se > 0.0 && se.is_finite() {
            Some(coef / se)
        } else {
            None
        }
    }
}

/// Results from fitting hurdle models to all features.
#[derive(Debug, Clone)]
pub struct HurdleFit {
    /// Individual fits for each feature.
    pub fits: Vec<HurdleFitSingle>,
    /// Coefficient names for binary model.
    pub binary_coefficient_names: Vec<String>,
    /// Coefficient names for count model.
    pub count_coefficient_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
    /// Count model type used.
    pub count_model: CountModel,
}

impl HurdleFit {
    /// Get the fit for a specific feature by ID.
    pub fn get_feature(&self, feature_id: &str) -> Option<&HurdleFitSingle> {
        self.fits.iter().find(|f| f.feature_id == feature_id)
    }

    /// Get binary coefficient index by name.
    pub fn binary_coefficient_index(&self, name: &str) -> Option<usize> {
        self.binary_coefficient_names.iter().position(|n| n == name)
    }

    /// Get count coefficient index by name.
    pub fn count_coefficient_index(&self, name: &str) -> Option<usize> {
        self.count_coefficient_names.iter().position(|n| n == name)
    }

    /// Get all binary coefficients for a specific coefficient name.
    pub fn binary_coefficients_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.binary_coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| {
                    f.binary_coefficients
                        .get(idx)
                        .copied()
                        .unwrap_or(f64::NAN)
                })
                .collect(),
        )
    }

    /// Get all count coefficients for a specific coefficient name.
    pub fn count_coefficients_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.count_coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| f.count_coefficients.get(idx).copied().unwrap_or(f64::NAN))
                .collect(),
        )
    }

    /// Get all binary standard errors for a specific coefficient name.
    pub fn binary_std_errors_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.binary_coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| f.binary_std_errors.get(idx).copied().unwrap_or(f64::NAN))
                .collect(),
        )
    }

    /// Get all count standard errors for a specific coefficient name.
    pub fn count_std_errors_for(&self, coefficient_name: &str) -> Option<Vec<f64>> {
        let idx = self.count_coefficient_index(coefficient_name)?;
        Some(
            self.fits
                .iter()
                .map(|f| f.count_std_errors.get(idx).copied().unwrap_or(f64::NAN))
                .collect(),
        )
    }

    /// Number of features.
    pub fn n_features(&self) -> usize {
        self.fits.len()
    }

    /// Number of binary coefficients.
    pub fn n_binary_coefficients(&self) -> usize {
        self.binary_coefficient_names.len()
    }

    /// Number of count coefficients.
    pub fn n_count_coefficients(&self) -> usize {
        self.count_coefficient_names.len()
    }

    /// Check if all fits converged.
    pub fn all_converged(&self) -> bool {
        self.fits.iter().all(|f| f.converged)
    }

    /// Count how many fits converged.
    pub fn n_converged(&self) -> usize {
        self.fits.iter().filter(|f| f.converged).count()
    }

    /// Get mean dispersion across features.
    pub fn mean_dispersion(&self) -> f64 {
        let sum: f64 = self.fits.iter().map(|f| f.dispersion).sum();
        sum / self.fits.len() as f64
    }
}

/// Fit hurdle model to count data.
///
/// # Arguments
/// * `counts` - Count matrix (features × samples)
/// * `design` - Design matrix from formula
///
/// # Returns
/// HurdleFit containing results for all features.
pub fn model_hurdle(counts: &CountMatrix, design: &DesignMatrix) -> Result<HurdleFit> {
    model_hurdle_with_config(counts, design, &HurdleConfig::default())
}

/// Fit hurdle model with custom configuration.
///
/// # Arguments
/// * `counts` - Count matrix (features × samples)
/// * `design` - Design matrix from formula
/// * `config` - Model configuration
///
/// # Returns
/// HurdleFit containing results for all features.
pub fn model_hurdle_with_config(
    counts: &CountMatrix,
    design: &DesignMatrix,
    config: &HurdleConfig,
) -> Result<HurdleFit> {
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

    // Need at least some samples for both parts
    let df_residual = n_samples.saturating_sub(2 * n_coef);
    if n_samples <= 2 * n_coef {
        return Err(DaaError::Numerical(
            "Model is saturated (n_samples <= 2 * n_coefficients)".to_string(),
        ));
    }

    let x = design.matrix();

    // Fit all features in parallel
    let fits: Vec<HurdleFitSingle> = (0..n_features)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples)
                .map(|j| counts.get(i, j) as f64)
                .collect();
            fit_single_hurdle(
                &y,
                &counts.feature_ids()[i],
                x,
                n_samples,
                n_coef,
                df_residual,
                config,
            )
        })
        .collect();

    Ok(HurdleFit {
        fits,
        binary_coefficient_names: design.coefficient_names().to_vec(),
        count_coefficient_names: design.coefficient_names().to_vec(),
        n_samples,
        count_model: config.count_model,
    })
}

/// Fit hurdle model to a single feature.
fn fit_single_hurdle(
    y: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    n_samples: usize,
    n_coef: usize,
    df_residual: usize,
    config: &HurdleConfig,
) -> HurdleFitSingle {
    // Create binary response: z_i = I(y_i > 0)
    let z: Vec<f64> = y.iter().map(|&yi| if yi > 0.0 { 1.0 } else { 0.0 }).collect();

    // Identify positive observations
    let positive_idx: Vec<usize> = y
        .iter()
        .enumerate()
        .filter(|(_, &yi)| yi > 0.0)
        .map(|(i, _)| i)
        .collect();

    let n_positive = positive_idx.len();

    // Fit binary component (logistic regression on all samples)
    let (binary_coef, binary_se, binary_iter, binary_converged, fitted_probs) =
        fit_binary_component(&z, x, n_samples, n_coef, config);

    // Fit count component (truncated NB on positive samples only)
    let (count_coef, count_se, dispersion, count_iter, count_converged, fitted_pos_mean) =
        if n_positive >= n_coef {
            // Extract positive observations
            let y_pos: Vec<f64> = positive_idx.iter().map(|&i| y[i]).collect();
            let x_pos = DMatrix::from_fn(n_positive, n_coef, |i, j| x[(positive_idx[i], j)]);

            fit_count_component(&y_pos, &x_pos, n_positive, n_coef, config)
        } else {
            // Not enough positive observations
            (
                vec![f64::NAN; n_coef],
                vec![f64::NAN; n_coef],
                f64::NAN,
                0,
                false,
                vec![f64::NAN; n_samples],
            )
        };

    // Compute fitted values for all samples: E[Y] = P(Y>0) * E[Y|Y>0]
    let fitted_values: Vec<f64> = if n_positive >= n_coef {
        // Compute E[Y|Y>0] for all samples using count model
        let eta_count: DVector<f64> = x * &DVector::from_vec(count_coef.clone());
        let mu_all: Vec<f64> = eta_count.iter().map(|&e| e.exp().max(MIN_MU)).collect();

        // E[Y] = p * E[Y|Y>0] where E[Y|Y>0] needs adjustment for truncation
        // For truncated NB: E[Y|Y>0] = mu / (1 - P(Y=0|mu,theta))
        fitted_probs
            .iter()
            .zip(mu_all.iter())
            .map(|(&p, &mu)| {
                let p_zero = if config.count_model == CountModel::TruncatedNB && dispersion > 0.0 {
                    (dispersion / (dispersion + mu)).powf(dispersion)
                } else {
                    (-mu).exp() // Poisson P(Y=0)
                };
                let e_y_given_pos = mu / (1.0 - p_zero).max(MIN_P);
                p * e_y_given_pos
            })
            .collect()
    } else {
        vec![f64::NAN; n_samples]
    };

    // Compute log-likelihood
    let ll_binary = binary_log_likelihood(&z, &fitted_probs);
    let ll_count = if n_positive >= n_coef {
        let y_pos: Vec<f64> = positive_idx.iter().map(|&i| y[i]).collect();
        let mu_pos: Vec<f64> = positive_idx
            .iter()
            .map(|&i| {
                let eta: f64 = (0..n_coef).map(|j| x[(i, j)] * count_coef[j]).sum();
                eta.exp().max(MIN_MU)
            })
            .collect();
        truncated_nb_log_likelihood(&y_pos, &mu_pos, dispersion)
    } else {
        0.0
    };

    let log_likelihood = ll_binary + ll_count;

    // Compute deviance (simplified)
    let deviance = -2.0 * log_likelihood;

    HurdleFitSingle {
        feature_id: feature_id.to_string(),
        binary_coefficients: binary_coef,
        binary_std_errors: binary_se,
        count_coefficients: count_coef,
        count_std_errors: count_se,
        dispersion,
        log_likelihood,
        deviance,
        df_residual,
        iterations: binary_iter.max(count_iter),
        converged: binary_converged && count_converged,
        fitted_values,
        fitted_probs,
        fitted_positive_mean: fitted_pos_mean,
    }
}

/// Fit binary component using logistic regression (IRLS).
fn fit_binary_component(
    z: &[f64],
    x: &DMatrix<f64>,
    n_samples: usize,
    n_coef: usize,
    config: &HurdleConfig,
) -> (Vec<f64>, Vec<f64>, usize, bool, Vec<f64>) {
    let z_vec = DVector::from_column_slice(z);

    // Initialize with proportion
    let p_mean = z.iter().sum::<f64>() / n_samples as f64;
    let p_mean = p_mean.clamp(MIN_P, MAX_P);
    let logit_mean = (p_mean / (1.0 - p_mean)).ln();

    let mut beta = DVector::zeros(n_coef);
    beta[0] = logit_mean;

    let mut converged = false;
    let mut iterations = 0;

    // IRLS iterations
    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Compute fitted probabilities
        let eta = x * &beta;
        let p: DVector<f64> = DVector::from_iterator(
            n_samples,
            eta.iter().map(|&e| (1.0 / (1.0 + (-e).exp())).clamp(MIN_P, MAX_P)),
        );

        // Compute weights: W = p * (1 - p)
        let w: DVector<f64> =
            DVector::from_iterator(n_samples, p.iter().map(|&pi| (pi * (1.0 - pi)).max(MIN_P)));

        // Compute working response: z_work = eta + (y - p) / (p * (1-p))
        let z_work: DVector<f64> = DVector::from_iterator(
            n_samples,
            (0..n_samples).map(|i| {
                let pi = p[i];
                let wi = w[i];
                let eta_i = eta[i];
                eta_i + (z_vec[i] - pi) / wi.max(MIN_P)
            }),
        );

        // Weighted least squares
        let w_sqrt = DVector::from_iterator(n_samples, w.iter().map(|wi| wi.sqrt()));
        let mut xw = x.clone();
        let mut zw = z_work.clone();
        for i in 0..n_samples {
            for j in 0..n_coef {
                xw[(i, j)] *= w_sqrt[i];
            }
            zw[i] *= w_sqrt[i];
        }

        let xtwx = xw.transpose() * &xw;
        let xtwz = xw.transpose() * &zw;

        let beta_new = match xtwx.clone().try_inverse() {
            Some(inv) => inv * xtwz,
            None => {
                // Return with current estimates
                let p_final: Vec<f64> = {
                    let eta = x * &beta;
                    eta.iter()
                        .map(|&e| (1.0 / (1.0 + (-e).exp())).clamp(MIN_P, MAX_P))
                        .collect()
                };
                return (
                    beta.iter().cloned().collect(),
                    vec![f64::NAN; n_coef],
                    iterations,
                    false,
                    p_final,
                );
            }
        };

        // Check convergence
        let delta: f64 = (&beta_new - &beta).iter().map(|d| d.abs()).sum();
        let scale: f64 = beta.iter().map(|b| b.abs()).sum::<f64>().max(1.0);

        beta = beta_new;

        if delta / scale < config.tol {
            converged = true;
            break;
        }
    }

    // Compute final fitted probabilities
    let eta = x * &beta;
    let p_final: Vec<f64> = eta
        .iter()
        .map(|&e| (1.0 / (1.0 + (-e).exp())).clamp(MIN_P, MAX_P))
        .collect();

    // Compute standard errors from Fisher information
    let w_final: DVector<f64> = DVector::from_iterator(
        n_samples,
        p_final.iter().map(|&pi| (pi * (1.0 - pi)).max(MIN_P)),
    );

    let mut xw = x.clone();
    for i in 0..n_samples {
        let w_sqrt = w_final[i].sqrt();
        for j in 0..n_coef {
            xw[(i, j)] *= w_sqrt;
        }
    }

    let fisher = xw.transpose() * &xw;
    let std_errors = match fisher.try_inverse() {
        Some(inv) => (0..n_coef).map(|j| inv[(j, j)].max(0.0).sqrt()).collect(),
        None => vec![f64::NAN; n_coef],
    };

    (
        beta.iter().cloned().collect(),
        std_errors,
        iterations,
        converged,
        p_final,
    )
}

/// Fit count component using truncated NB (IRLS).
fn fit_count_component(
    y: &[f64],
    x: &DMatrix<f64>,
    n_obs: usize,
    n_coef: usize,
    config: &HurdleConfig,
) -> (Vec<f64>, Vec<f64>, f64, usize, bool, Vec<f64>) {
    let y_vec = DVector::from_column_slice(y);

    // Initialize with log of mean
    let y_mean = y.iter().sum::<f64>() / n_obs as f64;
    let log_y_mean = y_mean.max(1.0).ln();

    let mut beta = DVector::zeros(n_coef);
    beta[0] = log_y_mean;

    // Initialize dispersion (theta)
    let mut theta = 1.0;

    let mut converged = false;
    let mut iterations = 0;

    // IRLS iterations
    for iter in 0..config.max_iter {
        iterations = iter + 1;

        // Compute mu = exp(X * beta)
        let eta = x * &beta;
        let mu: DVector<f64> =
            DVector::from_iterator(n_obs, eta.iter().map(|&e| e.exp().max(MIN_MU)));

        // For truncated NB, we need adjusted mean
        // E[Y|Y>0] = mu / (1 - P(Y=0|mu,theta))
        let adj_factor: DVector<f64> = DVector::from_iterator(
            n_obs,
            mu.iter().map(|&m| {
                let p_zero = if config.count_model == CountModel::TruncatedNB {
                    (theta / (theta + m)).powf(theta)
                } else {
                    (-m).exp()
                };
                (1.0 - p_zero).max(MIN_P)
            }),
        );

        // Compute weights for truncated NB
        // Weight = mu / [var factor * truncation adjustment]
        let w: DVector<f64> = if config.count_model == CountModel::TruncatedNB {
            DVector::from_iterator(
                n_obs,
                (0..n_obs).map(|i| {
                    let m = mu[i];
                    let var = m + m * m / theta; // NB variance
                    (m / var.max(MIN_MU)).max(MIN_P)
                }),
            )
        } else {
            // Poisson: variance = mean
            DVector::from_iterator(n_obs, mu.iter().map(|_| 1.0))
        };

        // Working response: z = eta + (y - mu/adj) * adj / mu
        let z: DVector<f64> = DVector::from_iterator(
            n_obs,
            (0..n_obs).map(|i| {
                let m = mu[i];
                let adj = adj_factor[i];
                let y_i = y_vec[i];
                let expected = m / adj; // E[Y|Y>0]
                eta[i] + (y_i - expected) / m.max(MIN_MU)
            }),
        );

        // Weighted least squares
        let w_sqrt = DVector::from_iterator(n_obs, w.iter().map(|wi| wi.sqrt()));
        let mut xw = x.clone();
        let mut zw = z.clone();
        for i in 0..n_obs {
            for j in 0..n_coef {
                xw[(i, j)] *= w_sqrt[i];
            }
            zw[i] *= w_sqrt[i];
        }

        let xtwx = xw.transpose() * &xw;
        let xtwz = xw.transpose() * &zw;

        let beta_new = match xtwx.clone().try_inverse() {
            Some(inv) => inv * xtwz,
            None => {
                let mu_final: Vec<f64> = {
                    let eta = x * &beta;
                    eta.iter().map(|&e| e.exp().max(MIN_MU)).collect()
                };
                return (
                    beta.iter().cloned().collect(),
                    vec![f64::NAN; n_coef],
                    theta,
                    iterations,
                    false,
                    mu_final,
                );
            }
        };

        // Check convergence
        let delta: f64 = (&beta_new - &beta).iter().map(|d| d.abs()).sum();
        let scale: f64 = beta.iter().map(|b| b.abs()).sum::<f64>().max(1.0);

        beta = beta_new;

        // Update dispersion if using NB
        if config.count_model == CountModel::TruncatedNB {
            let eta_new = x * &beta;
            let mu_new: DVector<f64> =
                DVector::from_iterator(n_obs, eta_new.iter().map(|&e| e.exp().max(MIN_MU)));
            theta = estimate_truncated_dispersion(&y_vec, &mu_new);
        }

        if delta / scale < config.tol {
            converged = true;
            break;
        }
    }

    // Compute final fitted values
    let eta = x * &beta;
    let mu_final: Vec<f64> = eta.iter().map(|&e| e.exp().max(MIN_MU)).collect();

    // Compute standard errors
    let w_final: DVector<f64> = if config.count_model == CountModel::TruncatedNB {
        DVector::from_iterator(
            n_obs,
            mu_final.iter().map(|&m| {
                let var = m + m * m / theta;
                (m / var.max(MIN_MU)).max(MIN_P)
            }),
        )
    } else {
        DVector::from_element(n_obs, 1.0)
    };

    let mut xw = x.clone();
    for i in 0..n_obs {
        let w_sqrt = w_final[i].sqrt();
        for j in 0..n_coef {
            xw[(i, j)] *= w_sqrt;
        }
    }

    let fisher = xw.transpose() * &xw;
    let std_errors = match fisher.try_inverse() {
        Some(inv) => (0..n_coef).map(|j| inv[(j, j)].max(0.0).sqrt()).collect(),
        None => vec![f64::NAN; n_coef],
    };

    (
        beta.iter().cloned().collect(),
        std_errors,
        theta,
        iterations,
        converged,
        mu_final,
    )
}

/// Estimate dispersion for truncated NB using method of moments.
fn estimate_truncated_dispersion(y: &DVector<f64>, mu: &DVector<f64>) -> f64 {
    let n = y.len() as f64;

    // Compute Pearson residuals
    let pearson_chi_sq: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mi)| {
            let m = mi.max(MIN_MU);
            let resid = yi - m;
            resid * resid / m
        })
        .sum();

    // Estimate theta from excess variance
    let excess = pearson_chi_sq - n;
    if excess > 0.0 {
        let mean_mu_sq: f64 = mu.iter().map(|&m| m * m).sum::<f64>() / n;
        (n * mean_mu_sq / excess).max(0.1)
    } else {
        1e6 // No overdispersion
    }
}

/// Compute binary log-likelihood.
fn binary_log_likelihood(z: &[f64], p: &[f64]) -> f64 {
    z.iter()
        .zip(p.iter())
        .map(|(&zi, &pi)| {
            let pi = pi.clamp(MIN_P, MAX_P);
            if zi > 0.5 {
                pi.ln()
            } else {
                (1.0 - pi).ln()
            }
        })
        .sum()
}

/// Compute truncated NB log-likelihood.
fn truncated_nb_log_likelihood(y: &[f64], mu: &[f64], theta: f64) -> f64 {
    y.iter()
        .zip(mu.iter())
        .map(|(&yi, &mi)| {
            let m = mi.max(MIN_MU);
            // NB log-likelihood
            let ll_nb = lgamma(yi + theta) - lgamma(theta) - lgamma(yi + 1.0)
                + theta * (theta / (theta + m)).ln()
                + yi * (m / (theta + m)).ln();

            // Subtract log(1 - P(Y=0)) for truncation
            let p_zero = (theta / (theta + m)).powf(theta);
            let log_truncation = (1.0 - p_zero).max(MIN_P).ln();

            ll_nb - log_truncation
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{Formula, Metadata};
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        // 3 features × 8 samples
        // Feature 0: ~50% zeros, similar counts in both groups
        // Feature 1: ~50% zeros, different counts between groups
        // Feature 2: Few zeros, for comparison
        let mut tri_mat = TriMat::new((3, 8));

        // Feature 0: No effect, ~50% zeros
        // Control: 0, 10, 0, 12
        // Treatment: 0, 11, 0, 9
        tri_mat.add_triplet(0, 0, 0);
        tri_mat.add_triplet(0, 1, 0);
        tri_mat.add_triplet(0, 2, 10);
        tri_mat.add_triplet(0, 3, 11);
        tri_mat.add_triplet(0, 4, 0);
        tri_mat.add_triplet(0, 5, 0);
        tri_mat.add_triplet(0, 6, 12);
        tri_mat.add_triplet(0, 7, 9);

        // Feature 1: Strong effect, ~50% zeros
        // Control: 0, 5, 0, 6
        // Treatment: 0, 25, 0, 28
        tri_mat.add_triplet(1, 0, 0);
        tri_mat.add_triplet(1, 1, 0);
        tri_mat.add_triplet(1, 2, 5);
        tri_mat.add_triplet(1, 3, 25);
        tri_mat.add_triplet(1, 4, 0);
        tri_mat.add_triplet(1, 5, 0);
        tri_mat.add_triplet(1, 6, 6);
        tri_mat.add_triplet(1, 7, 28);

        // Feature 2: Few zeros, moderate effect
        tri_mat.add_triplet(2, 0, 50);
        tri_mat.add_triplet(2, 1, 80);
        tri_mat.add_triplet(2, 2, 55);
        tri_mat.add_triplet(2, 3, 75);
        tri_mat.add_triplet(2, 4, 48);
        tri_mat.add_triplet(2, 5, 82);
        tri_mat.add_triplet(2, 6, 52);
        tri_mat.add_triplet(2, 7, 78);

        let feature_ids = vec!["sparse_no_effect".into(), "sparse_effect".into(), "dense".into()];
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
    fn test_model_hurdle_basic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        assert_eq!(fit.n_features(), 3);
        assert_eq!(fit.n_binary_coefficients(), 2);
        assert_eq!(fit.n_count_coefficients(), 2);
        assert_eq!(
            fit.binary_coefficient_names,
            vec!["(Intercept)", "grouptreatment"]
        );
    }

    #[test]
    fn test_model_hurdle_convergence() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        // At least some features should converge
        assert!(fit.n_converged() > 0, "At least one feature should converge");
    }

    #[test]
    fn test_model_hurdle_coefficients() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        // Dense feature with positive effect should have positive count coefficient
        let dense = fit.get_feature("dense").unwrap();
        let count_effect = dense.count_coefficients[1];
        assert!(
            count_effect > 0.0,
            "Dense feature should have positive count effect, got {}",
            count_effect
        );
    }

    #[test]
    fn test_model_hurdle_dispersion() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        // Dispersion should be positive
        for f in &fit.fits {
            if f.dispersion.is_finite() {
                assert!(f.dispersion > 0.0, "Dispersion should be positive");
            }
        }
    }

    #[test]
    fn test_model_hurdle_fitted_probs() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        // Fitted probabilities should be in (0, 1)
        for f in &fit.fits {
            for &p in &f.fitted_probs {
                assert!(
                    p > 0.0 && p < 1.0,
                    "Fitted prob should be in (0,1), got {}",
                    p
                );
            }
        }
    }

    #[test]
    fn test_model_hurdle_dimension_mismatch() {
        let counts = create_test_counts();

        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        writeln!(file, "S1\tcontrol").unwrap();
        writeln!(file, "S2\ttreatment").unwrap();
        file.flush().unwrap();
        let metadata = Metadata::from_tsv(file.path()).unwrap();

        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let result = model_hurdle(&counts, &design);
        assert!(result.is_err());
    }

    #[test]
    fn test_binary_z_statistic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        let dense = fit.get_feature("dense").unwrap();
        let z = dense.binary_z_statistic(1);
        // Should return a valid z-statistic
        assert!(z.is_some() || dense.binary_std_errors[1].is_nan());
    }

    #[test]
    fn test_count_z_statistic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        let dense = fit.get_feature("dense").unwrap();
        let z = dense.count_z_statistic(1);
        // Should return a valid z-statistic
        assert!(z.is_some());
        if let Some(z_val) = z {
            assert!(z_val.is_finite());
        }
    }

    #[test]
    fn test_config_with_poisson() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let config = HurdleConfig {
            count_model: CountModel::TruncatedPoisson,
            ..Default::default()
        };

        let fit = model_hurdle_with_config(&counts, &design, &config).unwrap();

        assert_eq!(fit.count_model, CountModel::TruncatedPoisson);
        assert!(fit.n_converged() > 0);
    }

    #[test]
    fn test_binary_log_likelihood() {
        let z = vec![1.0, 0.0, 1.0, 0.0];
        let p = vec![0.8, 0.2, 0.9, 0.1];

        let ll = binary_log_likelihood(&z, &p);
        assert!(ll.is_finite());
        assert!(ll < 0.0); // Log-likelihood should be negative

        // Perfect prediction should give higher likelihood
        let p_perfect = vec![0.99, 0.01, 0.99, 0.01];
        let ll_perfect = binary_log_likelihood(&z, &p_perfect);
        assert!(ll_perfect > ll);
    }
}
