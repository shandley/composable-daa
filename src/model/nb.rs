//! Negative Binomial GLM for count data.
//!
//! Implements a negative binomial generalized linear model with log link,
//! suitable for overdispersed count data common in microbiome studies.

use crate::data::{CountMatrix, DesignMatrix};
use crate::error::{DaaError, Result};
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Maximum iterations for IRLS convergence.
const MAX_ITER: usize = 25;

/// Convergence tolerance for coefficient changes.
const TOL: f64 = 1e-8;

/// Minimum value for mean to avoid log(0).
const MIN_MU: f64 = 1e-10;

/// Results from fitting a negative binomial model to a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NbFitSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Estimated coefficients (log scale).
    pub coefficients: Vec<f64>,
    /// Standard errors of coefficients.
    pub std_errors: Vec<f64>,
    /// Estimated dispersion parameter (theta).
    /// Variance = mu + mu^2/theta.
    pub dispersion: f64,
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
    /// Fitted values (predicted means).
    #[serde(skip)]
    pub fitted_values: Vec<f64>,
}

impl NbFitSingle {
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

/// Results from fitting negative binomial models to all features.
#[derive(Debug, Clone)]
pub struct NbFit {
    /// Individual fits for each feature.
    pub fits: Vec<NbFitSingle>,
    /// Coefficient names from the design matrix.
    pub coefficient_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
}

impl NbFit {
    /// Get the fit for a specific feature by ID.
    pub fn get_feature(&self, feature_id: &str) -> Option<&NbFitSingle> {
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

    /// Get mean dispersion across features.
    pub fn mean_dispersion(&self) -> f64 {
        let sum: f64 = self.fits.iter().map(|f| f.dispersion).sum();
        sum / self.fits.len() as f64
    }
}

/// Fit negative binomial GLM to count data.
///
/// Uses Iteratively Reweighted Least Squares (IRLS) with a log link function.
/// Dispersion is estimated per-feature using method of moments with
/// optional MLE refinement.
///
/// # Arguments
/// * `counts` - Count matrix (features × samples)
/// * `design` - Design matrix from formula
///
/// # Returns
/// NbFit containing results for all features.
pub fn model_nb(counts: &CountMatrix, design: &DesignMatrix) -> Result<NbFit> {
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

    // Fit all features in parallel
    let fits: Vec<NbFitSingle> = (0..n_features)
        .into_par_iter()
        .map(|i| {
            let y: Vec<f64> = (0..n_samples)
                .map(|j| counts.get(i, j) as f64)
                .collect();
            fit_single_nb(
                &y,
                &counts.feature_ids()[i],
                x,
                n_samples,
                n_coef,
                df_residual,
            )
        })
        .collect();

    Ok(NbFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        n_samples,
    })
}

/// Fit negative binomial model to a single feature using IRLS.
fn fit_single_nb(
    y: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    n_samples: usize,
    n_coef: usize,
    df_residual: usize,
) -> NbFitSingle {
    let y_vec = DVector::from_column_slice(y);

    // Initialize with Poisson-like starting values
    let y_mean = y.iter().sum::<f64>() / n_samples as f64;
    let log_y_mean = (y_mean.max(MIN_MU)).ln();

    // Start with intercept-only model, then refine
    let mut beta = DVector::zeros(n_coef);
    beta[0] = log_y_mean;

    // Compute initial mu
    let mut mu = compute_mu(x, &beta);

    // Initial dispersion estimate using method of moments
    let mut theta = estimate_dispersion_mom(&y_vec, &mu);

    let mut converged = false;
    let mut iterations = 0;

    // IRLS iterations
    for iter in 0..MAX_ITER {
        iterations = iter + 1;

        // Compute working weights: W = mu / (1 + mu/theta)
        let w: DVector<f64> = DVector::from_iterator(
            n_samples,
            mu.iter().map(|&m| m / (1.0 + m / theta)),
        );

        // Compute working response: z = eta + (y - mu) / mu
        // where eta = X*beta = log(mu)
        let eta: DVector<f64> = DVector::from_iterator(n_samples, mu.iter().map(|m| m.ln()));
        let z: DVector<f64> = DVector::from_iterator(
            n_samples,
            (0..n_samples).map(|i| {
                let m = mu[i].max(MIN_MU);
                eta[i] + (y_vec[i] - m) / m
            }),
        );

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
                return NbFitSingle {
                    feature_id: feature_id.to_string(),
                    coefficients: beta.iter().cloned().collect(),
                    std_errors: vec![f64::NAN; n_coef],
                    dispersion: theta,
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
        mu = compute_mu(x, &beta);

        // Update dispersion estimate
        theta = estimate_dispersion_mom(&y_vec, &mu);

        if delta / scale < TOL {
            converged = true;
            break;
        }
    }

    // Compute standard errors from Fisher information
    // Fisher information: I = X' W X where W = mu / (1 + mu/theta)
    let w_final: DVector<f64> = DVector::from_iterator(
        n_samples,
        mu.iter().map(|&m| m / (1.0 + m / theta)),
    );

    let mut xw_final = x.clone();
    for i in 0..n_samples {
        let w_sqrt = w_final[i].sqrt();
        for j in 0..n_coef {
            xw_final[(i, j)] *= w_sqrt;
        }
    }

    let fisher = xw_final.transpose() * &xw_final;
    let std_errors = match fisher.try_inverse() {
        Some(inv) => (0..n_coef).map(|j| inv[(j, j)].max(0.0).sqrt()).collect(),
        None => vec![f64::NAN; n_coef],
    };

    // Compute log-likelihood
    let log_lik = nb_log_likelihood(&y_vec, &mu, theta);

    // Compute deviance
    let deviance = nb_deviance(&y_vec, &mu, theta);

    NbFitSingle {
        feature_id: feature_id.to_string(),
        coefficients: beta.iter().cloned().collect(),
        std_errors,
        dispersion: theta,
        log_likelihood: log_lik,
        deviance,
        df_residual,
        iterations,
        converged,
        fitted_values: mu.iter().cloned().collect(),
    }
}

/// Compute mu = exp(X * beta).
fn compute_mu(x: &DMatrix<f64>, beta: &DVector<f64>) -> DVector<f64> {
    let eta = x * beta;
    DVector::from_iterator(eta.len(), eta.iter().map(|e| e.exp().max(MIN_MU)))
}

/// Estimate dispersion using method of moments.
///
/// For NB: Var(Y) = mu + mu^2/theta
/// Rearranging: theta = mu^2 / (Var(Y) - mu)
fn estimate_dispersion_mom(y: &DVector<f64>, mu: &DVector<f64>) -> f64 {
    let n = y.len() as f64;

    // Compute Pearson residuals squared
    let pearson_chi_sq: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mi)| {
            let m = mi.max(MIN_MU);
            let resid = yi - m;
            resid * resid / m
        })
        .sum();

    // Method of moments: estimate theta from Pearson chi-squared
    // Under NB: E[chi^2] = n - p + (n-p)/theta approximately
    // Rearranging: theta ≈ (n-p) / (chi^2 - (n-p))
    let excess = pearson_chi_sq - n;
    if excess > 0.0 {
        n / excess
    } else {
        // No overdispersion detected, use large theta (approaches Poisson)
        1e6
    }
}

/// Compute negative binomial log-likelihood.
fn nb_log_likelihood(y: &DVector<f64>, mu: &DVector<f64>, theta: f64) -> f64 {
    y.iter()
        .zip(mu.iter())
        .map(|(&yi, &mi)| {
            let m = mi.max(MIN_MU);
            // NB log-likelihood: lgamma(y+theta) - lgamma(theta) - lgamma(y+1)
            //                    + theta*log(theta/(theta+mu)) + y*log(mu/(theta+mu))
            let t = theta;
            lgamma(yi + t) - lgamma(t) - lgamma(yi + 1.0)
                + t * (t / (t + m)).ln()
                + yi * (m / (t + m)).ln()
        })
        .sum()
}

/// Compute negative binomial deviance.
fn nb_deviance(y: &DVector<f64>, mu: &DVector<f64>, theta: f64) -> f64 {
    let dev_sum: f64 = y
        .iter()
        .zip(mu.iter())
        .map(|(&yi, &mi)| {
            let m = mi.max(MIN_MU);
            let y_safe = yi.max(MIN_MU);

            // Deviance contribution: y*log(y/mu) - (y+theta)*log((y+theta)/(mu+theta))
            let term1 = if yi > 0.0 {
                yi * (y_safe / m).ln()
            } else {
                0.0
            };
            let term2 = (yi + theta) * ((yi + theta) / (m + theta)).ln();
            term1 - term2
        })
        .sum();
    2.0 * dev_sum
}

/// Log gamma function (approximation using Stirling's formula for large values).
fn lgamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    if x < 12.0 {
        // Use recurrence relation for small values
        let mut result = 0.0;
        let mut z = x;
        while z < 12.0 {
            result -= z.ln();
            z += 1.0;
        }
        result + lgamma(z)
    } else {
        // Stirling's approximation for large values
        let z = x;
        let c = [
            1.0 / 12.0,
            -1.0 / 360.0,
            1.0 / 1260.0,
            -1.0 / 1680.0,
        ];
        let mut sum = 0.0;
        let mut zp = z;
        for &ci in &c {
            sum += ci / zp;
            zp *= z * z;
        }
        (z - 0.5) * z.ln() - z + 0.5 * (2.0 * std::f64::consts::PI).ln() + sum
    }
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
        // 2 features × 8 samples
        // Feature 0: no group effect (similar counts in both groups)
        // Feature 1: strong group effect (higher in treatment)
        let mut tri_mat = TriMat::new((2, 8));

        // Feature 0: ~100 in both groups
        tri_mat.add_triplet(0, 0, 95);
        tri_mat.add_triplet(0, 1, 105);
        tri_mat.add_triplet(0, 2, 98);
        tri_mat.add_triplet(0, 3, 102);
        tri_mat.add_triplet(0, 4, 97);
        tri_mat.add_triplet(0, 5, 103);
        tri_mat.add_triplet(0, 6, 99);
        tri_mat.add_triplet(0, 7, 101);

        // Feature 1: ~50 in control, ~200 in treatment
        tri_mat.add_triplet(1, 0, 48);
        tri_mat.add_triplet(1, 1, 195);
        tri_mat.add_triplet(1, 2, 52);
        tri_mat.add_triplet(1, 3, 205);
        tri_mat.add_triplet(1, 4, 47);
        tri_mat.add_triplet(1, 5, 198);
        tri_mat.add_triplet(1, 6, 53);
        tri_mat.add_triplet(1, 7, 202);

        let feature_ids = vec!["no_effect".into(), "strong_effect".into()];
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
    fn test_model_nb_basic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        assert_eq!(fit.n_features(), 2);
        assert_eq!(fit.n_coefficients(), 2);
        assert_eq!(fit.coefficient_names, vec!["(Intercept)", "grouptreatment"]);
    }

    #[test]
    fn test_model_nb_convergence() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Both features should converge
        assert!(fit.all_converged(), "All fits should converge");
    }

    #[test]
    fn test_model_nb_coefficients() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Feature 0 should have near-zero treatment effect
        let no_effect = fit.get_feature("no_effect").unwrap();
        let coef_no_effect = no_effect.coefficients[1];
        assert!(
            coef_no_effect.abs() < 0.2,
            "No-effect feature should have small coefficient, got {}",
            coef_no_effect
        );

        // Feature 1 should have positive treatment effect
        // log(200/50) ≈ 1.39
        let strong_effect = fit.get_feature("strong_effect").unwrap();
        let coef_strong = strong_effect.coefficients[1];
        assert!(
            coef_strong > 1.0,
            "Strong-effect feature should have large positive coefficient, got {}",
            coef_strong
        );
    }

    #[test]
    fn test_model_nb_dispersion() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Dispersion should be positive
        for f in &fit.fits {
            assert!(f.dispersion > 0.0, "Dispersion should be positive");
        }
    }

    #[test]
    fn test_model_nb_std_errors() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Standard errors should be positive and finite
        for f in &fit.fits {
            for &se in &f.std_errors {
                assert!(se > 0.0 && se.is_finite(), "SE should be positive finite");
            }
        }
    }

    #[test]
    fn test_model_nb_fitted_values() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Fitted values should be positive
        for f in &fit.fits {
            for &fv in &f.fitted_values {
                assert!(fv > 0.0, "Fitted values should be positive");
            }
        }
    }

    #[test]
    fn test_model_nb_dimension_mismatch() {
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

        let result = model_nb(&counts, &design);
        assert!(result.is_err());
    }

    #[test]
    fn test_lgamma() {
        // Test against known values
        assert_relative_eq!(lgamma(1.0), 0.0, epsilon = 1e-6);
        assert_relative_eq!(lgamma(2.0), 0.0, epsilon = 1e-6); // log(1!) = 0
        assert_relative_eq!(lgamma(3.0), (2.0_f64).ln(), epsilon = 1e-6); // log(2!) = log(2)
        assert_relative_eq!(lgamma(4.0), (6.0_f64).ln(), epsilon = 1e-6); // log(3!) = log(6)
    }

    #[test]
    fn test_z_statistic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        // Strong effect feature should have large z-statistic
        let strong = fit.get_feature("strong_effect").unwrap();
        let z = strong.z_statistic(1).unwrap();
        assert!(z.abs() > 2.0, "Strong effect should have large z-statistic");

        // No effect feature should have small z-statistic
        let no_effect = fit.get_feature("no_effect").unwrap();
        let z_no = no_effect.z_statistic(1).unwrap();
        assert!(
            z_no.abs() < 2.0,
            "No effect should have small z-statistic, got {}",
            z_no
        );
    }
}
