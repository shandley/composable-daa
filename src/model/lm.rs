//! Linear model fitting via OLS.

use crate::data::DesignMatrix;
use crate::error::{DaaError, Result};
use crate::normalize::TransformedMatrix;
use nalgebra::{DMatrix, DVector};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Results from fitting a linear model to a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LmFitSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Estimated coefficients.
    pub coefficients: Vec<f64>,
    /// Standard errors of coefficients.
    pub std_errors: Vec<f64>,
    /// Residuals.
    #[serde(skip)]
    pub residuals: Vec<f64>,
    /// Residual standard error (sigma).
    pub sigma: f64,
    /// R-squared.
    pub r_squared: f64,
    /// Degrees of freedom (residual).
    pub df_residual: usize,
    /// Whether the fit was successful.
    pub converged: bool,
}

impl LmFitSingle {
    /// Get coefficient by name (requires coefficient names from design matrix).
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
}

/// Results from fitting linear models to all features.
#[derive(Debug, Clone)]
pub struct LmFit {
    /// Individual fits for each feature.
    pub fits: Vec<LmFitSingle>,
    /// Coefficient names from the design matrix.
    pub coefficient_names: Vec<String>,
    /// Number of samples.
    pub n_samples: usize,
}

impl LmFit {
    /// Get the fit for a specific feature by ID.
    pub fn get_feature(&self, feature_id: &str) -> Option<&LmFitSingle> {
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

    /// Number of coefficients (including intercept if present).
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
}

/// Fit linear models to transformed abundance data.
///
/// Fits OLS regression for each feature against the design matrix.
/// Uses QR decomposition for numerical stability.
///
/// # Arguments
/// * `transformed` - CLR or other transformed abundance data
/// * `design` - Design matrix from formula
///
/// # Returns
/// LmFit containing results for all features.
pub fn model_lm(transformed: &TransformedMatrix, design: &DesignMatrix) -> Result<LmFit> {
    let n_features = transformed.n_features();
    let n_samples = transformed.n_samples();
    let n_coef = design.n_coefficients();

    // Validate dimensions
    if design.n_samples() != n_samples {
        return Err(DaaError::DimensionMismatch {
            expected: n_samples,
            actual: design.n_samples(),
        });
    }

    // Check rank
    let df_residual = n_samples.saturating_sub(n_coef);
    if df_residual == 0 {
        return Err(DaaError::Numerical(
            "Model is saturated (n_samples <= n_coefficients)".to_string(),
        ));
    }

    let x = design.matrix();
    let xtx = x.transpose() * x;

    // Try to compute (X'X)^-1
    let xtx_inv = xtx.clone().try_inverse().ok_or_else(|| {
        DaaError::Numerical("Design matrix is singular (X'X not invertible)".to_string())
    })?;

    // Fit all features in parallel
    let fits: Vec<LmFitSingle> = (0..n_features)
        .into_par_iter()
        .map(|i| {
            fit_single_feature(
                &transformed.row(i),
                &transformed.feature_ids[i],
                x,
                &xtx_inv,
                n_samples,
                n_coef,
                df_residual,
            )
        })
        .collect();

    Ok(LmFit {
        fits,
        coefficient_names: design.coefficient_names().to_vec(),
        n_samples,
    })
}

/// Fit a single feature using pre-computed (X'X)^-1.
fn fit_single_feature(
    y: &[f64],
    feature_id: &str,
    x: &DMatrix<f64>,
    xtx_inv: &DMatrix<f64>,
    n_samples: usize,
    n_coef: usize,
    df_residual: usize,
) -> LmFitSingle {
    let y_vec = DVector::from_column_slice(y);

    // Coefficients: beta = (X'X)^-1 X'y
    let xty = x.transpose() * &y_vec;
    let beta = xtx_inv * xty;
    let coefficients: Vec<f64> = beta.iter().cloned().collect();

    // Residuals: e = y - X*beta
    let y_hat = x * &beta;
    let residuals_vec = &y_vec - &y_hat;
    let residuals: Vec<f64> = residuals_vec.iter().cloned().collect();

    // Residual sum of squares
    let rss: f64 = residuals.iter().map(|e| e * e).sum();

    // Sigma (residual standard error)
    let sigma = if df_residual > 0 {
        (rss / df_residual as f64).sqrt()
    } else {
        f64::NAN
    };

    // Standard errors: SE = sigma * sqrt(diag((X'X)^-1))
    let std_errors: Vec<f64> = (0..n_coef)
        .map(|j| sigma * xtx_inv[(j, j)].sqrt())
        .collect();

    // R-squared
    let y_mean = y.iter().sum::<f64>() / n_samples as f64;
    let tss: f64 = y.iter().map(|yi| (yi - y_mean).powi(2)).sum();
    let r_squared = if tss > 0.0 { 1.0 - rss / tss } else { 0.0 };

    LmFitSingle {
        feature_id: feature_id.to_string(),
        coefficients,
        std_errors,
        residuals,
        sigma,
        r_squared,
        df_residual,
        converged: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{Formula, Metadata};
    use approx::assert_relative_eq;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        writeln!(file, "S1\tcontrol\t25").unwrap();
        writeln!(file, "S2\ttreatment\t30").unwrap();
        writeln!(file, "S3\tcontrol\t35").unwrap();
        writeln!(file, "S4\ttreatment\t28").unwrap();
        writeln!(file, "S5\tcontrol\t32").unwrap();
        writeln!(file, "S6\ttreatment\t27").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn create_test_transformed() -> TransformedMatrix {
        // Create transformed data with known effect
        // Feature 0: no group effect
        // Feature 1: clear group effect (treatment > control)
        let data = DMatrix::from_row_slice(2, 6, &[
            1.0, 1.2, 0.9, 1.1, 1.0, 1.1,  // feature 0: no effect
            1.0, 3.0, 1.2, 2.8, 0.9, 3.1,  // feature 1: treatment ~3x control
        ]);

        TransformedMatrix {
            data,
            feature_ids: vec!["feat_0".into(), "feat_1".into()],
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into(), "S4".into(), "S5".into(), "S6".into()],
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 6],
        }
    }

    #[test]
    fn test_model_lm_basic() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();

        assert_eq!(fit.n_features(), 2);
        assert_eq!(fit.n_coefficients(), 2);
        assert_eq!(fit.coefficient_names, vec!["(Intercept)", "grouptreatment"]);
        assert!(fit.all_converged());
    }

    #[test]
    fn test_model_lm_coefficients() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();

        // Feature 0 should have small treatment effect
        let feat0 = fit.get_feature("feat_0").unwrap();
        let treatment_effect_0 = feat0.coefficients[1];
        assert!(treatment_effect_0.abs() < 0.5);

        // Feature 1 should have large treatment effect (~2)
        let feat1 = fit.get_feature("feat_1").unwrap();
        let treatment_effect_1 = feat1.coefficients[1];
        assert!(treatment_effect_1 > 1.5);
    }

    #[test]
    fn test_model_lm_r_squared() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();

        // Feature 1 should have high R-squared (group explains variance)
        let feat1 = fit.get_feature("feat_1").unwrap();
        assert!(feat1.r_squared > 0.8);
    }

    #[test]
    fn test_model_lm_degrees_of_freedom() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group + age").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();

        // 6 samples - 3 coefficients = 3 df
        assert_eq!(fit.fits[0].df_residual, 3);
    }

    #[test]
    fn test_dimension_mismatch() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        // Create transformed with wrong number of samples
        let data = DMatrix::from_row_slice(2, 4, &[1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
        let transformed = TransformedMatrix {
            data,
            feature_ids: vec!["a".into(), "b".into()],
            sample_ids: vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()],
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 4],
        };

        let result = model_lm(&transformed, &design);
        assert!(result.is_err());
    }
}
