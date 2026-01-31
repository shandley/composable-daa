//! Wald test for coefficient significance.

use crate::error::{DaaError, Result};
use crate::model::{HurdleFit, LmFit, LmmFit, NbFit, ZinbFit};
use serde::{Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};

/// Result of Wald test for a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaldResultSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Coefficient name being tested.
    pub coefficient: String,
    /// Estimated coefficient value.
    pub estimate: f64,
    /// Standard error.
    pub std_error: f64,
    /// Wald statistic (t-statistic).
    pub statistic: f64,
    /// P-value (two-sided).
    pub p_value: f64,
    /// Degrees of freedom (f64 for Satterthwaite precision).
    pub df: f64,
}

/// Results of Wald tests for all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaldResult {
    /// Individual test results.
    pub results: Vec<WaldResultSingle>,
    /// Coefficient name being tested.
    pub coefficient: String,
}

impl WaldResult {
    /// Number of tests.
    pub fn len(&self) -> usize {
        self.results.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }

    /// Get p-values for all features.
    pub fn p_values(&self) -> Vec<f64> {
        self.results.iter().map(|r| r.p_value).collect()
    }

    /// Get estimates for all features.
    pub fn estimates(&self) -> Vec<f64> {
        self.results.iter().map(|r| r.estimate).collect()
    }

    /// Get standard errors for all features.
    pub fn std_errors(&self) -> Vec<f64> {
        self.results.iter().map(|r| r.std_error).collect()
    }

    /// Get feature IDs.
    pub fn feature_ids(&self) -> Vec<&str> {
        self.results.iter().map(|r| r.feature_id.as_str()).collect()
    }

    /// Get result for a specific feature.
    pub fn get_feature(&self, feature_id: &str) -> Option<&WaldResultSingle> {
        self.results.iter().find(|r| r.feature_id == feature_id)
    }
}

/// Perform Wald test on linear model coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the t-distribution.
/// The Wald statistic is t = β / SE(β), compared to a t-distribution
/// with df_residual degrees of freedom.
///
/// # Arguments
/// * `fit` - Linear model fit results
/// * `coefficient` - Name of coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald(fit: &LmFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Coefficient '{}' not found. Available: {:?}",
            coefficient, fit.coefficient_names
        ))
    })?;

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate t-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value
            let p_value = if !statistic.is_nan() && df > 0.0 {
                let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
                2.0 * (1.0 - t_dist.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Test multiple coefficients at once.
pub fn test_wald_multiple(fit: &LmFit, coefficients: &[&str]) -> Result<Vec<WaldResult>> {
    coefficients
        .iter()
        .map(|c| test_wald(fit, c))
        .collect()
}

/// Perform Wald test on negative binomial GLM coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the normal distribution (z-test).
/// For GLMs, the Wald statistic z = β / SE(β) is compared to a standard
/// normal distribution for large samples.
///
/// # Arguments
/// * `fit` - Negative binomial model fit results
/// * `coefficient` - Name of coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_nb(fit: &NbFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Coefficient '{}' not found. Available: {:?}",
            coefficient, fit.coefficient_names
        ))
    })?;

    let normal = Normal::new(0.0, 1.0).unwrap();

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate z-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using normal distribution
            let p_value = if !statistic.is_nan() {
                2.0 * (1.0 - normal.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Perform Wald test on ZINB count model coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the normal distribution (z-test).
/// This tests coefficients from the count component (negative binomial) of the ZINB model.
///
/// # Arguments
/// * `fit` - ZINB model fit results
/// * `coefficient` - Name of coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_zinb(fit: &ZinbFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Coefficient '{}' not found. Available: {:?}",
            coefficient, fit.coefficient_names
        ))
    })?;

    let normal = Normal::new(0.0, 1.0).unwrap();

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate z-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using normal distribution
            let p_value = if !statistic.is_nan() {
                2.0 * (1.0 - normal.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Perform Wald test on LMM fixed effect coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the t-distribution.
/// For LMMs, the Wald statistic t = β / SE(β) is compared to a t-distribution
/// with degrees of freedom determined by the df_method used during fitting.
///
/// Degrees of freedom priority:
/// 1. Kenward-Roger df (if available) - most conservative, uses adjusted SE
/// 2. Satterthwaite df (if available) - accounts for variance component uncertainty
/// 3. Residual df (n - p) - naive approach
///
/// When Kenward-Roger was used, this also uses the bias-corrected standard errors
/// (std_errors_kr) which are typically larger than the naive SE.
///
/// # Arguments
/// * `fit` - LMM fit results
/// * `coefficient` - Name of fixed effect coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_lmm(fit: &LmmFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Coefficient '{}' not found. Available: {:?}",
            coefficient, fit.coefficient_names
        ))
    })?;

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);

            // Get standard error: prefer KR adjusted SE > original SE
            let std_error = if let Some(ref kr_se) = f.std_errors_kr {
                kr_se.get(coef_idx).copied().unwrap_or_else(|| {
                    f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN)
                })
            } else {
                f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN)
            };

            // Get df: prefer KR df > Satterthwaite df > residual df
            let df = if let Some(ref kr_df) = f.df_kenward_roger {
                kr_df.get(coef_idx).copied().unwrap_or(f.df_residual)
            } else if let Some(ref satt_df) = f.df_satterthwaite {
                satt_df.get(coef_idx).copied().unwrap_or(f.df_residual)
            } else {
                f.df_residual
            };

            // Calculate t-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using t-distribution
            let p_value = if !statistic.is_nan() && df > 0.0 {
                let t_dist = StudentsT::new(0.0, 1.0, df).unwrap();
                2.0 * (1.0 - t_dist.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Perform Wald test on hurdle model count coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the normal distribution (z-test).
/// This tests coefficients from the count component (truncated NB) of the hurdle model.
///
/// # Arguments
/// * `fit` - Hurdle model fit results
/// * `coefficient` - Name of count coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_hurdle_count(fit: &HurdleFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.count_coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Count coefficient '{}' not found. Available: {:?}",
            coefficient, fit.count_coefficient_names
        ))
    })?;

    let normal = Normal::new(0.0, 1.0).unwrap();

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.count_coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.count_std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate z-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using normal distribution
            let p_value = if !statistic.is_nan() {
                2.0 * (1.0 - normal.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Perform Wald test on hurdle model binary coefficients.
///
/// Tests H0: β = 0 vs H1: β ≠ 0 using the normal distribution (z-test).
/// This tests coefficients from the binary component (logistic) of the hurdle model.
///
/// # Arguments
/// * `fit` - Hurdle model fit results
/// * `coefficient` - Name of binary coefficient to test
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_hurdle_binary(fit: &HurdleFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.binary_coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Binary coefficient '{}' not found. Available: {:?}",
            coefficient, fit.binary_coefficient_names
        ))
    })?;

    let normal = Normal::new(0.0, 1.0).unwrap();

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.binary_coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.binary_std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate z-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using normal distribution
            let p_value = if !statistic.is_nan() {
                2.0 * (1.0 - normal.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

/// Perform Wald test on ZINB zero-inflation model coefficients.
///
/// Tests H0: γ = 0 vs H1: γ ≠ 0 using the normal distribution (z-test).
/// This tests coefficients from the zero-inflation component (logistic) of the ZINB model.
///
/// # Arguments
/// * `fit` - ZINB model fit results
/// * `coefficient` - Name of zero-inflation coefficient to test (e.g., "zi_intercept", "zi_grouptreatment")
///
/// # Returns
/// WaldResult containing test statistics and p-values for all features.
pub fn test_wald_zinb_zi(fit: &ZinbFit, coefficient: &str) -> Result<WaldResult> {
    let coef_idx = fit.zi_coefficient_index(coefficient).ok_or_else(|| {
        DaaError::InvalidParameter(format!(
            "Zero-inflation coefficient '{}' not found. Available: {:?}",
            coefficient, fit.zi_coefficient_names
        ))
    })?;

    let normal = Normal::new(0.0, 1.0).unwrap();

    let results: Vec<WaldResultSingle> = fit
        .fits
        .iter()
        .map(|f| {
            let estimate = f.zi_coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let std_error = f.zi_std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            let df = f.df_residual as f64;

            // Calculate z-statistic
            let statistic = if std_error > 0.0 && !std_error.is_nan() {
                estimate / std_error
            } else {
                f64::NAN
            };

            // Calculate two-sided p-value using normal distribution
            let p_value = if !statistic.is_nan() {
                2.0 * (1.0 - normal.cdf(statistic.abs()))
            } else {
                f64::NAN
            };

            WaldResultSingle {
                feature_id: f.feature_id.clone(),
                coefficient: coefficient.to_string(),
                estimate,
                std_error,
                statistic,
                p_value,
                df,
            }
        })
        .collect();

    Ok(WaldResult {
        results,
        coefficient: coefficient.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DesignMatrix, Formula, Metadata};
    use crate::model::model_lm;
    use crate::normalize::TransformedMatrix;
    use approx::assert_relative_eq;
    use nalgebra::DMatrix;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
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

    fn create_test_transformed() -> TransformedMatrix {
        // Feature 0: no effect (random noise around 1.0)
        // Feature 1: strong effect (control ~1.0, treatment ~3.0)
        let data = DMatrix::from_row_slice(2, 8, &[
            1.1, 0.9, 1.0, 1.1, 0.9, 1.0, 1.1, 0.9,  // no effect
            1.0, 3.0, 1.1, 2.9, 0.9, 3.1, 1.0, 3.0,  // strong effect
        ]);

        TransformedMatrix {
            data,
            feature_ids: vec!["no_effect".into(), "strong_effect".into()],
            sample_ids: (1..=8).map(|i| format!("S{}", i)).collect(),
            transformation: "CLR".to_string(),
            geometric_means: vec![1.0; 8],
        }
    }

    #[test]
    fn test_wald_basic() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();
        let wald = test_wald(&fit, "grouptreatment").unwrap();

        assert_eq!(wald.len(), 2);
        assert_eq!(wald.coefficient, "grouptreatment");
    }

    #[test]
    fn test_wald_significance() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();
        let wald = test_wald(&fit, "grouptreatment").unwrap();

        // Feature with no effect should have large p-value
        let no_effect = wald.get_feature("no_effect").unwrap();
        assert!(no_effect.p_value > 0.1, "No-effect feature should not be significant");

        // Feature with strong effect should have small p-value
        let strong_effect = wald.get_feature("strong_effect").unwrap();
        assert!(strong_effect.p_value < 0.001, "Strong-effect feature should be significant");
    }

    #[test]
    fn test_wald_estimates() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();
        let wald = test_wald(&fit, "grouptreatment").unwrap();

        let strong_effect = wald.get_feature("strong_effect").unwrap();
        // Treatment effect should be approximately 2 (3.0 - 1.0)
        assert_relative_eq!(strong_effect.estimate, 2.0, epsilon = 0.2);
    }

    #[test]
    fn test_wald_invalid_coefficient() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();
        let result = test_wald(&fit, "nonexistent");

        assert!(result.is_err());
    }

    #[test]
    fn test_wald_p_value_bounds() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();
        let transformed = create_test_transformed();

        let fit = model_lm(&transformed, &design).unwrap();
        let wald = test_wald(&fit, "grouptreatment").unwrap();

        for result in &wald.results {
            assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
        }
    }
}
