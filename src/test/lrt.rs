//! Likelihood Ratio Test (LRT) for nested model comparison.
//!
//! LRT compares a full model to a reduced model (null model) by comparing
//! their log-likelihoods. The test statistic follows a chi-squared distribution.
//!
//! LRT is often preferred over Wald test because:
//! - More accurate for small samples
//! - Better behaved near parameter boundaries
//! - Can test multiple coefficients simultaneously

use crate::data::{CountMatrix, DesignMatrix};
use crate::error::{DaaError, Result};
use crate::model::nb::{model_nb, NbFit};
use crate::model::zinb::{model_zinb, ZinbFit};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use statrs::distribution::{ChiSquared, ContinuousCDF};

/// Result of LRT for a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LrtResultSingle {
    /// Feature identifier.
    pub feature_id: String,
    /// Coefficient(s) being tested.
    pub coefficients_tested: Vec<String>,
    /// Log-likelihood of full model.
    pub ll_full: f64,
    /// Log-likelihood of reduced model.
    pub ll_reduced: f64,
    /// LRT statistic: 2 * (ll_full - ll_reduced).
    pub statistic: f64,
    /// Degrees of freedom (number of coefficients tested).
    pub df: usize,
    /// P-value from chi-squared distribution.
    pub p_value: f64,
}

/// Results of LRT for all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LrtResult {
    /// Individual test results.
    pub results: Vec<LrtResultSingle>,
    /// Coefficient(s) being tested.
    pub coefficients_tested: Vec<String>,
    /// Degrees of freedom.
    pub df: usize,
}

impl LrtResult {
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

    /// Get LRT statistics for all features.
    pub fn statistics(&self) -> Vec<f64> {
        self.results.iter().map(|r| r.statistic).collect()
    }

    /// Get feature IDs.
    pub fn feature_ids(&self) -> Vec<&str> {
        self.results.iter().map(|r| r.feature_id.as_str()).collect()
    }

    /// Get result for a specific feature.
    pub fn get_feature(&self, feature_id: &str) -> Option<&LrtResultSingle> {
        self.results.iter().find(|r| r.feature_id == feature_id)
    }
}

/// Perform LRT on negative binomial model.
///
/// Tests H0: β = 0 for the specified coefficient(s) by comparing:
/// - Full model: includes all coefficients
/// - Reduced model: excludes the tested coefficient(s)
///
/// The LRT statistic is 2 * (LL_full - LL_reduced), compared to a
/// chi-squared distribution with df = number of coefficients tested.
///
/// # Arguments
/// * `full_fit` - Full NB model fit
/// * `counts` - Original count matrix
/// * `design` - Original design matrix
/// * `coefficients` - Names of coefficient(s) to test
///
/// # Returns
/// LrtResult containing test statistics and p-values for all features.
pub fn test_lrt_nb(
    full_fit: &NbFit,
    counts: &CountMatrix,
    design: &DesignMatrix,
    coefficients: &[&str],
) -> Result<LrtResult> {
    // Validate coefficients exist
    let coef_indices: Vec<usize> = coefficients
        .iter()
        .map(|c| {
            full_fit.coefficient_index(c).ok_or_else(|| {
                DaaError::InvalidParameter(format!(
                    "Coefficient '{}' not found. Available: {:?}",
                    c, full_fit.coefficient_names
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    // Cannot remove intercept
    for &idx in &coef_indices {
        if full_fit.coefficient_names.get(idx).map(|s| s.as_str()) == Some("(Intercept)") {
            return Err(DaaError::InvalidParameter(
                "Cannot test intercept with LRT; use a different coefficient".to_string(),
            ));
        }
    }

    let df = coefficients.len();

    // Create reduced design matrix by removing tested columns
    let reduced_design = create_reduced_design(design, &coef_indices)?;

    // Fit reduced model
    let reduced_fit = model_nb(counts, &reduced_design)?;

    // Compute LRT for each feature
    let chi_sq = ChiSquared::new(df as f64).unwrap();

    let results: Vec<LrtResultSingle> = full_fit
        .fits
        .iter()
        .zip(reduced_fit.fits.iter())
        .map(|(full, reduced)| {
            let ll_full = full.log_likelihood;
            let ll_reduced = reduced.log_likelihood;

            // LRT statistic
            let statistic = 2.0 * (ll_full - ll_reduced);
            let statistic = statistic.max(0.0); // Should be non-negative

            // P-value from chi-squared
            let p_value = if statistic.is_finite() {
                1.0 - chi_sq.cdf(statistic)
            } else {
                f64::NAN
            };

            LrtResultSingle {
                feature_id: full.feature_id.clone(),
                coefficients_tested: coefficients.iter().map(|s| s.to_string()).collect(),
                ll_full,
                ll_reduced,
                statistic,
                df,
                p_value,
            }
        })
        .collect();

    Ok(LrtResult {
        results,
        coefficients_tested: coefficients.iter().map(|s| s.to_string()).collect(),
        df,
    })
}

/// Perform LRT on ZINB count model coefficients.
///
/// Tests H0: β = 0 for the specified count model coefficient(s).
///
/// # Arguments
/// * `full_fit` - Full ZINB model fit
/// * `counts` - Original count matrix
/// * `design` - Original design matrix
/// * `coefficients` - Names of coefficient(s) to test
pub fn test_lrt_zinb(
    full_fit: &ZinbFit,
    counts: &CountMatrix,
    design: &DesignMatrix,
    coefficients: &[&str],
) -> Result<LrtResult> {
    // Validate coefficients exist
    let coef_indices: Vec<usize> = coefficients
        .iter()
        .map(|c| {
            full_fit.coefficient_index(c).ok_or_else(|| {
                DaaError::InvalidParameter(format!(
                    "Coefficient '{}' not found. Available: {:?}",
                    c, full_fit.coefficient_names
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    // Cannot remove intercept
    for &idx in &coef_indices {
        if full_fit.coefficient_names.get(idx).map(|s| s.as_str()) == Some("(Intercept)") {
            return Err(DaaError::InvalidParameter(
                "Cannot test intercept with LRT; use a different coefficient".to_string(),
            ));
        }
    }

    let df = coefficients.len();

    // Create reduced design matrix
    let reduced_design = create_reduced_design(design, &coef_indices)?;

    // Fit reduced ZINB model
    let reduced_fit = model_zinb(counts, &reduced_design)?;

    // Compute LRT for each feature
    let chi_sq = ChiSquared::new(df as f64).unwrap();

    let results: Vec<LrtResultSingle> = full_fit
        .fits
        .iter()
        .zip(reduced_fit.fits.iter())
        .map(|(full, reduced)| {
            let ll_full = full.log_likelihood;
            let ll_reduced = reduced.log_likelihood;

            let statistic = 2.0 * (ll_full - ll_reduced);
            let statistic = statistic.max(0.0);

            let p_value = if statistic.is_finite() {
                1.0 - chi_sq.cdf(statistic)
            } else {
                f64::NAN
            };

            LrtResultSingle {
                feature_id: full.feature_id.clone(),
                coefficients_tested: coefficients.iter().map(|s| s.to_string()).collect(),
                ll_full,
                ll_reduced,
                statistic,
                df,
                p_value,
            }
        })
        .collect();

    Ok(LrtResult {
        results,
        coefficients_tested: coefficients.iter().map(|s| s.to_string()).collect(),
        df,
    })
}

/// Create a reduced design matrix by removing specified columns.
fn create_reduced_design(
    design: &DesignMatrix,
    remove_indices: &[usize],
) -> Result<DesignMatrix> {
    let original = design.matrix();
    let n_rows = original.nrows();
    let n_cols = original.ncols();

    // Columns to keep
    let keep_cols: Vec<usize> = (0..n_cols)
        .filter(|i| !remove_indices.contains(i))
        .collect();

    if keep_cols.is_empty() {
        return Err(DaaError::InvalidParameter(
            "Cannot remove all coefficients from design matrix".to_string(),
        ));
    }

    // Build reduced matrix
    let mut reduced_data = Vec::with_capacity(n_rows * keep_cols.len());
    for row in 0..n_rows {
        for &col in &keep_cols {
            reduced_data.push(original[(row, col)]);
        }
    }

    let reduced_matrix = DMatrix::from_row_slice(n_rows, keep_cols.len(), &reduced_data);

    // Build reduced coefficient names
    let original_names = design.coefficient_names();
    let reduced_names: Vec<String> = keep_cols
        .iter()
        .map(|&i| original_names[i].clone())
        .collect();

    Ok(DesignMatrix::from_matrix(
        reduced_matrix,
        reduced_names,
        design.sample_ids().to_vec(),
    ))
}

/// Perform LRT comparing two pre-fitted NB models.
///
/// Use this when you have already fit both models and want to compare them.
/// The full model should be nested within (contain all parameters of) the reduced model.
///
/// # Arguments
/// * `full_fit` - Full NB model fit (more parameters)
/// * `reduced_fit` - Reduced NB model fit (fewer parameters)
/// * `df` - Degrees of freedom (difference in number of parameters)
pub fn test_lrt_nb_fitted(
    full_fit: &NbFit,
    reduced_fit: &NbFit,
    df: usize,
) -> Result<LrtResult> {
    if full_fit.fits.len() != reduced_fit.fits.len() {
        return Err(DaaError::DimensionMismatch {
            expected: full_fit.fits.len(),
            actual: reduced_fit.fits.len(),
        });
    }

    let chi_sq = ChiSquared::new(df as f64).unwrap();

    let results: Vec<LrtResultSingle> = full_fit
        .fits
        .iter()
        .zip(reduced_fit.fits.iter())
        .map(|(full, reduced)| {
            let ll_full = full.log_likelihood;
            let ll_reduced = reduced.log_likelihood;

            let statistic = 2.0 * (ll_full - ll_reduced);
            let statistic = statistic.max(0.0);

            let p_value = if statistic.is_finite() {
                1.0 - chi_sq.cdf(statistic)
            } else {
                f64::NAN
            };

            LrtResultSingle {
                feature_id: full.feature_id.clone(),
                coefficients_tested: vec!["(custom)".to_string()],
                ll_full,
                ll_reduced,
                statistic,
                df,
                p_value,
            }
        })
        .collect();

    Ok(LrtResult {
        results,
        coefficients_tested: vec!["(custom)".to_string()],
        df,
    })
}

/// Perform LRT comparing two pre-fitted ZINB models.
pub fn test_lrt_zinb_fitted(
    full_fit: &ZinbFit,
    reduced_fit: &ZinbFit,
    df: usize,
) -> Result<LrtResult> {
    if full_fit.fits.len() != reduced_fit.fits.len() {
        return Err(DaaError::DimensionMismatch {
            expected: full_fit.fits.len(),
            actual: reduced_fit.fits.len(),
        });
    }

    let chi_sq = ChiSquared::new(df as f64).unwrap();

    let results: Vec<LrtResultSingle> = full_fit
        .fits
        .iter()
        .zip(reduced_fit.fits.iter())
        .map(|(full, reduced)| {
            let ll_full = full.log_likelihood;
            let ll_reduced = reduced.log_likelihood;

            let statistic = 2.0 * (ll_full - ll_reduced);
            let statistic = statistic.max(0.0);

            let p_value = if statistic.is_finite() {
                1.0 - chi_sq.cdf(statistic)
            } else {
                f64::NAN
            };

            LrtResultSingle {
                feature_id: full.feature_id.clone(),
                coefficients_tested: vec!["(custom)".to_string()],
                ll_full,
                ll_reduced,
                statistic,
                df,
                p_value,
            }
        })
        .collect();

    Ok(LrtResult {
        results,
        coefficients_tested: vec!["(custom)".to_string()],
        df,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{Formula, Metadata};
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        for i in 1..=20 {
            let group = if i <= 10 { "control" } else { "treatment" };
            let age = 25 + i;
            writeln!(file, "S{}\t{}\t{}", i, group, age).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn create_test_counts() -> CountMatrix {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "feature_id\t{}",
            (1..=20).map(|i| format!("S{}", i)).collect::<Vec<_>>().join("\t")
        )
        .unwrap();

        // Feature with no effect
        writeln!(
            file,
            "no_effect\t{}",
            (1..=20).map(|_| "50").collect::<Vec<_>>().join("\t")
        )
        .unwrap();

        // Feature with treatment effect
        let effect_counts: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "20".to_string() } else { "80".to_string() })
            .collect();
        writeln!(file, "with_effect\t{}", effect_counts.join("\t")).unwrap();

        // Feature with moderate effect
        let mod_counts: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "30".to_string() } else { "50".to_string() })
            .collect();
        writeln!(file, "moderate\t{}", mod_counts.join("\t")).unwrap();

        file.flush().unwrap();
        CountMatrix::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_lrt_nb_basic() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let lrt = test_lrt_nb(&full_fit, &counts, &design, &["grouptreatment"]).unwrap();

        assert_eq!(lrt.len(), 3);
        assert_eq!(lrt.df, 1);
        assert_eq!(lrt.coefficients_tested, vec!["grouptreatment"]);
    }

    #[test]
    fn test_lrt_nb_significance() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let lrt = test_lrt_nb(&full_fit, &counts, &design, &["grouptreatment"]).unwrap();

        // Feature with no effect should have large p-value
        let no_effect = lrt.get_feature("no_effect").unwrap();
        assert!(
            no_effect.p_value > 0.05,
            "No-effect feature should not be significant, got p={}",
            no_effect.p_value
        );

        // Feature with effect should have small p-value
        let with_effect = lrt.get_feature("with_effect").unwrap();
        assert!(
            with_effect.p_value < 0.05,
            "Effect feature should be significant, got p={}",
            with_effect.p_value
        );
    }

    #[test]
    fn test_lrt_nb_statistic_positive() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let lrt = test_lrt_nb(&full_fit, &counts, &design, &["grouptreatment"]).unwrap();

        // LRT statistic should be non-negative
        for result in &lrt.results {
            assert!(
                result.statistic >= 0.0,
                "LRT statistic should be non-negative"
            );
        }
    }

    #[test]
    fn test_lrt_nb_p_value_bounds() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let lrt = test_lrt_nb(&full_fit, &counts, &design, &["grouptreatment"]).unwrap();

        for result in &lrt.results {
            assert!(
                result.p_value >= 0.0 && result.p_value <= 1.0,
                "P-value should be in [0, 1]"
            );
        }
    }

    #[test]
    fn test_lrt_nb_invalid_coefficient() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let result = test_lrt_nb(&full_fit, &counts, &design, &["nonexistent"]);

        assert!(result.is_err());
    }

    #[test]
    fn test_lrt_nb_cannot_test_intercept() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();
        let result = test_lrt_nb(&full_fit, &counts, &design, &["(Intercept)"]);

        assert!(result.is_err());
    }

    #[test]
    fn test_lrt_zinb_basic() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_zinb(&counts, &design).unwrap();
        let lrt = test_lrt_zinb(&full_fit, &counts, &design, &["grouptreatment"]).unwrap();

        assert_eq!(lrt.len(), 3);
        assert_eq!(lrt.df, 1);
    }

    #[test]
    fn test_lrt_nb_prefitted_models() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();

        // Full model with group
        let formula_full = Formula::parse("~ group").unwrap();
        let design_full = DesignMatrix::from_formula(&metadata, &formula_full).unwrap();
        let full_fit = model_nb(&counts, &design_full).unwrap();

        // Reduced model (intercept only)
        let formula_reduced = Formula::parse("~ 1").unwrap();
        let design_reduced = DesignMatrix::from_formula(&metadata, &formula_reduced).unwrap();
        let reduced_fit = model_nb(&counts, &design_reduced).unwrap();

        let lrt = super::test_lrt_nb_fitted(&full_fit, &reduced_fit, 1).unwrap();

        assert_eq!(lrt.len(), 3);
        assert_eq!(lrt.df, 1);

        // Check that results make sense
        for result in &lrt.results {
            assert!(result.statistic >= 0.0);
            assert!(result.p_value >= 0.0 && result.p_value <= 1.0);
        }
    }

    #[test]
    fn test_create_reduced_design() {
        let metadata = create_test_metadata();
        let formula = Formula::parse("~ group + age").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        // Original should have 3 columns: intercept, group, age
        assert_eq!(design.matrix().ncols(), 3);

        // Remove group (index 1)
        let reduced = create_reduced_design(&design, &[1]).unwrap();
        assert_eq!(reduced.matrix().ncols(), 2);
        assert_eq!(reduced.coefficient_names(), &["(Intercept)", "age"]);
    }

    #[test]
    fn test_lrt_multiple_coefficients() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group + age").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let full_fit = model_nb(&counts, &design).unwrap();

        // Test both group and age together
        let lrt = test_lrt_nb(&full_fit, &counts, &design, &["grouptreatment", "age"]).unwrap();

        assert_eq!(lrt.df, 2); // Testing 2 coefficients
        assert_eq!(lrt.coefficients_tested.len(), 2);
    }
}
