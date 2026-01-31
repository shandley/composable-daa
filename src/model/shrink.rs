//! Effect size shrinkage for log fold change estimates.
//!
//! Shrinkage stabilizes effect size estimates by pulling noisy estimates
//! toward zero. Features with high uncertainty (large standard errors)
//! are shrunk more than features with precise estimates.
//!
//! This is similar to DESeq2's `lfcShrink` function and improves
//! interpretability for sparse/low-count data.

use crate::model::hurdle::HurdleFit;
use crate::model::nb::NbFit;
use crate::model::zinb::ZinbFit;
use crate::test::WaldResult;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Shrinkage method to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ShrinkageMethod {
    /// Normal prior shrinkage (simple, fast).
    /// Assumes LFC ~ N(0, σ²) prior.
    Normal,
    /// Adaptive shrinkage using mixture model.
    /// Estimates prior from data, allows asymmetric shrinkage.
    Adaptive,
    /// No shrinkage, just compute confidence intervals.
    None,
}

impl Default for ShrinkageMethod {
    fn default() -> Self {
        ShrinkageMethod::Normal
    }
}

/// Configuration for shrinkage estimation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShrinkageConfig {
    /// Shrinkage method.
    pub method: ShrinkageMethod,
    /// Confidence level for intervals (default: 0.95).
    pub confidence_level: f64,
    /// Prior mean for LFC (default: 0.0, shrink toward zero).
    pub prior_mean: f64,
    /// Fixed prior variance (if None, estimated from data).
    pub prior_var: Option<f64>,
    /// Minimum prior variance to prevent over-shrinkage.
    pub min_prior_var: f64,
}

impl Default for ShrinkageConfig {
    fn default() -> Self {
        Self {
            method: ShrinkageMethod::Normal,
            confidence_level: 0.95,
            prior_mean: 0.0,
            prior_var: None,
            min_prior_var: 0.01,
        }
    }
}

/// Shrunk effect size estimate for a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShrunkEstimate {
    /// Feature identifier.
    pub feature_id: String,
    /// Original (unshrunk) log fold change estimate.
    pub lfc_unshrunk: f64,
    /// Shrunk log fold change estimate.
    pub lfc_shrunk: f64,
    /// Original standard error.
    pub se_unshrunk: f64,
    /// Posterior standard error (after shrinkage).
    pub se_posterior: f64,
    /// Lower bound of confidence interval.
    pub ci_lower: f64,
    /// Upper bound of confidence interval.
    pub ci_upper: f64,
    /// Shrinkage factor (0 = fully shrunk, 1 = no shrinkage).
    pub shrinkage_factor: f64,
    /// S-value: surprisal value (-log2(p-value)).
    pub s_value: f64,
}

impl ShrunkEstimate {
    /// Check if the confidence interval excludes zero.
    pub fn is_significant(&self) -> bool {
        self.ci_lower > 0.0 || self.ci_upper < 0.0
    }

    /// Get the fold change (exp of LFC).
    pub fn fold_change(&self) -> f64 {
        self.lfc_shrunk.exp()
    }
}

/// Results of shrinkage estimation for all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ShrinkageResult {
    /// Individual shrunk estimates.
    pub estimates: Vec<ShrunkEstimate>,
    /// Coefficient name that was shrunk.
    pub coefficient: String,
    /// Estimated prior variance (if empirical Bayes).
    pub prior_var: f64,
    /// Prior mean used.
    pub prior_mean: f64,
    /// Shrinkage method used.
    pub method: ShrinkageMethod,
    /// Confidence level used.
    pub confidence_level: f64,
}

impl ShrinkageResult {
    /// Get estimate for a specific feature.
    pub fn get_feature(&self, feature_id: &str) -> Option<&ShrunkEstimate> {
        self.estimates.iter().find(|e| e.feature_id == feature_id)
    }

    /// Get all shrunk LFC values.
    pub fn lfc_shrunk(&self) -> Vec<f64> {
        self.estimates.iter().map(|e| e.lfc_shrunk).collect()
    }

    /// Get all unshrunk LFC values.
    pub fn lfc_unshrunk(&self) -> Vec<f64> {
        self.estimates.iter().map(|e| e.lfc_unshrunk).collect()
    }

    /// Get features with significant effects (CI excludes zero).
    pub fn significant_features(&self) -> Vec<&ShrunkEstimate> {
        self.estimates.iter().filter(|e| e.is_significant()).collect()
    }

    /// Number of features.
    pub fn len(&self) -> usize {
        self.estimates.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.estimates.is_empty()
    }
}

/// Shrink log fold change estimates from Wald test results.
///
/// Uses empirical Bayes shrinkage to stabilize estimates:
/// 1. Estimates prior variance from the data
/// 2. Computes posterior mean (shrunk estimate)
/// 3. Computes posterior variance and confidence intervals
///
/// # Arguments
/// * `wald` - Wald test results containing estimates and standard errors
/// * `config` - Shrinkage configuration
///
/// # Returns
/// ShrinkageResult with shrunk estimates and confidence intervals.
pub fn shrink_lfc(wald: &WaldResult, config: &ShrinkageConfig) -> ShrinkageResult {
    let estimates: Vec<f64> = wald.estimates();
    let std_errors: Vec<f64> = wald.std_errors();
    let p_values: Vec<f64> = wald.p_values();

    // Filter out invalid estimates
    let valid: Vec<(usize, f64, f64, f64)> = estimates
        .iter()
        .zip(std_errors.iter())
        .zip(p_values.iter())
        .enumerate()
        .filter(|(_, ((e, se), _))| e.is_finite() && se.is_finite() && **se > 0.0)
        .map(|(i, ((e, se), p))| (i, *e, *se, *p))
        .collect();

    // Estimate prior variance if not provided
    let prior_var = match config.prior_var {
        Some(v) => v.max(config.min_prior_var),
        None => estimate_prior_variance(&valid, config),
    };

    // Compute critical value for confidence intervals
    let alpha = 1.0 - config.confidence_level;
    let z_crit = normal_quantile(1.0 - alpha / 2.0);

    // Compute shrunk estimates
    let shrunk_estimates: Vec<ShrunkEstimate> = wald
        .results
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let lfc = r.estimate;
            let se = r.std_error;
            let p_value = r.p_value;

            if !lfc.is_finite() || !se.is_finite() || se <= 0.0 {
                // Return unshrunk for invalid data
                return ShrunkEstimate {
                    feature_id: r.feature_id.clone(),
                    lfc_unshrunk: lfc,
                    lfc_shrunk: lfc,
                    se_unshrunk: se,
                    se_posterior: se,
                    ci_lower: f64::NEG_INFINITY,
                    ci_upper: f64::INFINITY,
                    shrinkage_factor: 1.0,
                    s_value: if p_value > 0.0 { -p_value.log2() } else { f64::INFINITY },
                };
            }

            let (lfc_shrunk, se_posterior, shrinkage_factor) = match config.method {
                ShrinkageMethod::Normal => {
                    shrink_normal(lfc, se, config.prior_mean, prior_var)
                }
                ShrinkageMethod::Adaptive => {
                    shrink_adaptive(lfc, se, config.prior_mean, prior_var)
                }
                ShrinkageMethod::None => {
                    (lfc, se, 1.0)
                }
            };

            let ci_lower = lfc_shrunk - z_crit * se_posterior;
            let ci_upper = lfc_shrunk + z_crit * se_posterior;

            let s_value = if p_value > 0.0 && p_value.is_finite() {
                -p_value.log2()
            } else if p_value == 0.0 {
                f64::INFINITY
            } else {
                f64::NAN
            };

            ShrunkEstimate {
                feature_id: r.feature_id.clone(),
                lfc_unshrunk: lfc,
                lfc_shrunk,
                se_unshrunk: se,
                se_posterior,
                ci_lower,
                ci_upper,
                shrinkage_factor,
                s_value,
            }
        })
        .collect();

    ShrinkageResult {
        estimates: shrunk_estimates,
        coefficient: wald.coefficient.clone(),
        prior_var,
        prior_mean: config.prior_mean,
        method: config.method,
        confidence_level: config.confidence_level,
    }
}

/// Shrink estimates from NB model fit.
///
/// # Arguments
/// * `fit` - NB model fit
/// * `coefficient` - Name of coefficient to shrink
/// * `config` - Shrinkage configuration
pub fn shrink_lfc_nb(
    fit: &NbFit,
    coefficient: &str,
    config: &ShrinkageConfig,
) -> Option<ShrinkageResult> {
    let coef_idx = fit.coefficient_index(coefficient)?;

    let estimates: Vec<(String, f64, f64)> = fit
        .fits
        .iter()
        .map(|f| {
            let est = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let se = f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            (f.feature_id.clone(), est, se)
        })
        .collect();

    Some(shrink_estimates(&estimates, coefficient, config))
}

/// Shrink estimates from ZINB model fit.
///
/// # Arguments
/// * `fit` - ZINB model fit
/// * `coefficient` - Name of coefficient to shrink
/// * `config` - Shrinkage configuration
pub fn shrink_lfc_zinb(
    fit: &ZinbFit,
    coefficient: &str,
    config: &ShrinkageConfig,
) -> Option<ShrinkageResult> {
    let coef_idx = fit.coefficient_index(coefficient)?;

    let estimates: Vec<(String, f64, f64)> = fit
        .fits
        .iter()
        .map(|f| {
            let est = f.coefficients.get(coef_idx).copied().unwrap_or(f64::NAN);
            let se = f.std_errors.get(coef_idx).copied().unwrap_or(f64::NAN);
            (f.feature_id.clone(), est, se)
        })
        .collect();

    Some(shrink_estimates(&estimates, coefficient, config))
}

/// Which component of the hurdle model to shrink.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HurdleComponent {
    /// Binary component (logistic regression for P(Y > 0)).
    Binary,
    /// Count component (truncated NB for E[Y | Y > 0]).
    Count,
}

impl Default for HurdleComponent {
    fn default() -> Self {
        HurdleComponent::Count
    }
}

/// Shrink estimates from hurdle model fit.
///
/// Hurdle models have two components:
/// - Binary: logistic regression for P(Y > 0), coefficients on logit scale
/// - Count: truncated NB for E[Y | Y > 0], coefficients on log scale
///
/// Use `component` to specify which component to shrink.
///
/// # Arguments
/// * `fit` - Hurdle model fit
/// * `coefficient` - Name of coefficient to shrink
/// * `component` - Which component to shrink (Binary or Count)
/// * `config` - Shrinkage configuration
pub fn shrink_lfc_hurdle(
    fit: &HurdleFit,
    coefficient: &str,
    component: HurdleComponent,
    config: &ShrinkageConfig,
) -> Option<ShrinkageResult> {
    let estimates: Vec<(String, f64, f64)> = match component {
        HurdleComponent::Binary => {
            let coef_idx = fit.binary_coefficient_index(coefficient)?;
            fit.fits
                .iter()
                .map(|f| {
                    let est = f
                        .binary_coefficients
                        .get(coef_idx)
                        .copied()
                        .unwrap_or(f64::NAN);
                    let se = f
                        .binary_std_errors
                        .get(coef_idx)
                        .copied()
                        .unwrap_or(f64::NAN);
                    (f.feature_id.clone(), est, se)
                })
                .collect()
        }
        HurdleComponent::Count => {
            let coef_idx = fit.count_coefficient_index(coefficient)?;
            fit.fits
                .iter()
                .map(|f| {
                    let est = f
                        .count_coefficients
                        .get(coef_idx)
                        .copied()
                        .unwrap_or(f64::NAN);
                    let se = f
                        .count_std_errors
                        .get(coef_idx)
                        .copied()
                        .unwrap_or(f64::NAN);
                    (f.feature_id.clone(), est, se)
                })
                .collect()
        }
    };

    // Add component suffix to coefficient name for clarity
    let coef_name = match component {
        HurdleComponent::Binary => format!("{} (binary)", coefficient),
        HurdleComponent::Count => format!("{} (count)", coefficient),
    };

    Some(shrink_estimates(&estimates, &coef_name, config))
}

/// Internal function to shrink a set of estimates.
fn shrink_estimates(
    estimates: &[(String, f64, f64)],
    coefficient: &str,
    config: &ShrinkageConfig,
) -> ShrinkageResult {
    // Filter valid estimates for prior estimation
    let valid: Vec<(usize, f64, f64, f64)> = estimates
        .iter()
        .enumerate()
        .filter(|(_, (_, e, se))| e.is_finite() && se.is_finite() && *se > 0.0)
        .map(|(i, (_, e, se))| (i, *e, *se, 1.0)) // dummy p-value
        .collect();

    let prior_var = match config.prior_var {
        Some(v) => v.max(config.min_prior_var),
        None => estimate_prior_variance(&valid, config),
    };

    let alpha = 1.0 - config.confidence_level;
    let z_crit = normal_quantile(1.0 - alpha / 2.0);

    let shrunk_estimates: Vec<ShrunkEstimate> = estimates
        .iter()
        .map(|(feature_id, lfc, se)| {
            if !lfc.is_finite() || !se.is_finite() || *se <= 0.0 {
                return ShrunkEstimate {
                    feature_id: feature_id.clone(),
                    lfc_unshrunk: *lfc,
                    lfc_shrunk: *lfc,
                    se_unshrunk: *se,
                    se_posterior: *se,
                    ci_lower: f64::NEG_INFINITY,
                    ci_upper: f64::INFINITY,
                    shrinkage_factor: 1.0,
                    s_value: f64::NAN,
                };
            }

            let (lfc_shrunk, se_posterior, shrinkage_factor) = match config.method {
                ShrinkageMethod::Normal => {
                    shrink_normal(*lfc, *se, config.prior_mean, prior_var)
                }
                ShrinkageMethod::Adaptive => {
                    shrink_adaptive(*lfc, *se, config.prior_mean, prior_var)
                }
                ShrinkageMethod::None => {
                    (*lfc, *se, 1.0)
                }
            };

            let ci_lower = lfc_shrunk - z_crit * se_posterior;
            let ci_upper = lfc_shrunk + z_crit * se_posterior;

            // Compute approximate p-value and s-value
            let z = lfc_shrunk / se_posterior;
            let p_value = 2.0 * (1.0 - normal_cdf(z.abs()));
            let s_value = if p_value > 0.0 { -p_value.log2() } else { f64::INFINITY };

            ShrunkEstimate {
                feature_id: feature_id.clone(),
                lfc_unshrunk: *lfc,
                lfc_shrunk,
                se_unshrunk: *se,
                se_posterior,
                ci_lower,
                ci_upper,
                shrinkage_factor,
                s_value,
            }
        })
        .collect();

    ShrinkageResult {
        estimates: shrunk_estimates,
        coefficient: coefficient.to_string(),
        prior_var,
        prior_mean: config.prior_mean,
        method: config.method,
        confidence_level: config.confidence_level,
    }
}

/// Estimate prior variance using empirical Bayes.
///
/// Uses method of moments: Var(β̂) = Var(β) + E[SE²]
/// So: Var(β) = Var(β̂) - mean(SE²)
fn estimate_prior_variance(
    valid: &[(usize, f64, f64, f64)],
    config: &ShrinkageConfig,
) -> f64 {
    if valid.is_empty() {
        return config.min_prior_var;
    }

    let estimates: Vec<f64> = valid.iter().map(|(_, e, _, _)| *e).collect();
    let variances: Vec<f64> = valid.iter().map(|(_, _, se, _)| se * se).collect();

    let n = estimates.len() as f64;

    // Variance of estimates
    let mean_est = estimates.iter().sum::<f64>() / n;
    let var_est = estimates.iter().map(|e| (e - mean_est).powi(2)).sum::<f64>() / (n - 1.0).max(1.0);

    // Mean of squared standard errors
    let mean_var = variances.iter().sum::<f64>() / n;

    // Prior variance = total variance - noise variance
    let prior_var = (var_est - mean_var).max(config.min_prior_var);

    prior_var
}

/// Normal shrinkage (shrink toward prior mean).
///
/// Posterior mean: μ_post = (μ_prior/σ²_prior + β̂/SE²) / (1/σ²_prior + 1/SE²)
/// Posterior var:  σ²_post = 1 / (1/σ²_prior + 1/SE²)
fn shrink_normal(
    estimate: f64,
    se: f64,
    prior_mean: f64,
    prior_var: f64,
) -> (f64, f64, f64) {
    let obs_var = se * se;

    // Posterior precision = prior precision + observation precision
    let post_precision = 1.0 / prior_var + 1.0 / obs_var;
    let post_var = 1.0 / post_precision;
    let post_se = post_var.sqrt();

    // Posterior mean is weighted average
    let post_mean = (prior_mean / prior_var + estimate / obs_var) / post_precision;

    // Shrinkage factor: how much of the original estimate is retained
    // 1 = no shrinkage, 0 = fully shrunk to prior
    let shrinkage_factor = (1.0 / obs_var) / post_precision;

    (post_mean, post_se, shrinkage_factor)
}

/// Adaptive shrinkage using scaled mixture.
///
/// Similar to ashr: allows for more shrinkage of estimates near zero
/// and less shrinkage of clearly non-null estimates.
fn shrink_adaptive(
    estimate: f64,
    se: f64,
    prior_mean: f64,
    prior_var: f64,
) -> (f64, f64, f64) {
    // Use a two-component mixture: null (point mass at 0) + non-null (normal)
    let obs_var = se * se;

    // Likelihood under null (estimate is noise around 0)
    let lik_null = normal_pdf(estimate, prior_mean, obs_var.sqrt());

    // Likelihood under alternative (estimate is from prior + noise)
    let alt_var = prior_var + obs_var;
    let lik_alt = normal_pdf(estimate, prior_mean, alt_var.sqrt());

    // Posterior probability of being non-null (assuming 50/50 prior)
    let prob_alt = lik_alt / (lik_null + lik_alt + 1e-300);

    // If likely null, shrink toward prior_mean
    // If likely non-null, shrink less
    let post_precision = 1.0 / prior_var + 1.0 / obs_var;
    let post_var = 1.0 / post_precision;

    // Weighted combination based on probability
    let shrunk_mean = prior_mean * (1.0 - prob_alt)
        + ((prior_mean / prior_var + estimate / obs_var) / post_precision) * prob_alt;

    // Posterior variance accounts for mixture uncertainty
    let post_se = (post_var * prob_alt + obs_var * (1.0 - prob_alt)).sqrt();

    let shrinkage_factor = prob_alt * (1.0 / obs_var) / post_precision;

    (shrunk_mean, post_se, shrinkage_factor)
}

/// Standard normal PDF.
fn normal_pdf(x: f64, mean: f64, sd: f64) -> f64 {
    let z = (x - mean) / sd;
    (-0.5 * z * z).exp() / (sd * (2.0 * PI).sqrt())
}

/// Standard normal CDF (approximation).
fn normal_cdf(x: f64) -> f64 {
    // Abramowitz and Stegun approximation
    let t = 1.0 / (1.0 + 0.2316419 * x.abs());
    let d = 0.3989423 * (-x * x / 2.0).exp();
    let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));

    if x > 0.0 {
        1.0 - p
    } else {
        p
    }
}

/// Standard normal quantile (approximation).
fn normal_quantile(p: f64) -> f64 {
    // Rational approximation
    if p <= 0.0 {
        return f64::NEG_INFINITY;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }

    let p = if p > 0.5 { 1.0 - p } else { p };
    let t = (-2.0 * p.ln()).sqrt();

    let c0 = 2.515517;
    let c1 = 0.802853;
    let c2 = 0.010328;
    let d1 = 1.432788;
    let d2 = 0.189269;
    let d3 = 0.001308;

    let z = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);

    if p > 0.5 {
        -z
    } else {
        z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{CountMatrix, DesignMatrix, Formula, Metadata};
    use crate::model::model_nb;
    use crate::test::test_wald_nb;
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

        // Feature with strong effect
        let strong: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "10".to_string() } else { "100".to_string() })
            .collect();
        writeln!(file, "strong_effect\t{}", strong.join("\t")).unwrap();

        // Feature with weak effect (noisy)
        let weak: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "40".to_string() } else { "60".to_string() })
            .collect();
        writeln!(file, "weak_effect\t{}", weak.join("\t")).unwrap();

        // Low count feature (high variance)
        let low: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "2".to_string() } else { "8".to_string() })
            .collect();
        writeln!(file, "low_count\t{}", low.join("\t")).unwrap();

        file.flush().unwrap();
        CountMatrix::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_shrink_lfc_basic() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());

        assert_eq!(shrunk.len(), 4);
        assert_eq!(shrunk.coefficient, "grouptreatment");
    }

    #[test]
    fn test_shrinkage_reduces_extreme_estimates() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());

        // Low count feature should be shrunk more (shrinkage_factor closer to 0)
        let low_count = shrunk.get_feature("low_count").unwrap();
        let strong = shrunk.get_feature("strong_effect").unwrap();

        // Low count has higher variance, so should be shrunk more
        // (shrinkage_factor = how much of original is retained)
        assert!(
            low_count.shrinkage_factor <= strong.shrinkage_factor,
            "Low count feature should be shrunk at least as much as strong effect"
        );
    }

    #[test]
    fn test_shrinkage_toward_zero() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());

        // Shrunk estimates should be closer to zero than unshrunk
        for est in &shrunk.estimates {
            if est.lfc_unshrunk.is_finite() && est.se_unshrunk > 0.0 {
                assert!(
                    est.lfc_shrunk.abs() <= est.lfc_unshrunk.abs() + 0.001,
                    "Shrunk LFC {} should be closer to zero than unshrunk {}",
                    est.lfc_shrunk,
                    est.lfc_unshrunk
                );
            }
        }
    }

    #[test]
    fn test_confidence_intervals() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());

        for est in &shrunk.estimates {
            if est.lfc_shrunk.is_finite() {
                assert!(est.ci_lower < est.lfc_shrunk);
                assert!(est.ci_upper > est.lfc_shrunk);
                assert!(est.ci_lower < est.ci_upper);
            }
        }
    }

    #[test]
    fn test_no_shrinkage_method() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let config = ShrinkageConfig {
            method: ShrinkageMethod::None,
            ..Default::default()
        };
        let shrunk = shrink_lfc(&wald, &config);

        // With no shrinkage, estimates should be unchanged
        for est in &shrunk.estimates {
            if est.lfc_unshrunk.is_finite() {
                assert_relative_eq!(est.lfc_shrunk, est.lfc_unshrunk, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_adaptive_shrinkage() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let config = ShrinkageConfig {
            method: ShrinkageMethod::Adaptive,
            ..Default::default()
        };
        let shrunk = shrink_lfc(&wald, &config);

        assert_eq!(shrunk.method, ShrinkageMethod::Adaptive);
        assert_eq!(shrunk.len(), 4);
    }

    #[test]
    fn test_shrink_lfc_nb() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        let shrunk = shrink_lfc_nb(&fit, "grouptreatment", &ShrinkageConfig::default());
        assert!(shrunk.is_some());

        let shrunk = shrunk.unwrap();
        assert_eq!(shrunk.len(), 4);
    }

    #[test]
    fn test_shrink_lfc_invalid_coefficient() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();

        let shrunk = shrink_lfc_nb(&fit, "nonexistent", &ShrinkageConfig::default());
        assert!(shrunk.is_none());
    }

    #[test]
    fn test_significant_features() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());
        let sig = shrunk.significant_features();

        // Strong effect should be significant
        assert!(sig.iter().any(|e| e.feature_id == "strong_effect"));
    }

    #[test]
    fn test_fold_change() {
        let est = ShrunkEstimate {
            feature_id: "test".to_string(),
            lfc_unshrunk: 1.0,
            lfc_shrunk: 0.693, // ln(2)
            se_unshrunk: 0.1,
            se_posterior: 0.1,
            ci_lower: 0.5,
            ci_upper: 0.9,
            shrinkage_factor: 0.9,
            s_value: 5.0,
        };

        assert_relative_eq!(est.fold_change(), 2.0, epsilon = 0.01);
    }

    #[test]
    fn test_normal_cdf() {
        assert_relative_eq!(normal_cdf(0.0), 0.5, epsilon = 0.001);
        assert!(normal_cdf(2.0) > 0.97);
        assert!(normal_cdf(-2.0) < 0.03);
    }

    #[test]
    fn test_normal_quantile() {
        assert_relative_eq!(normal_quantile(0.5), 0.0, epsilon = 0.01);
        assert_relative_eq!(normal_quantile(0.975), 1.96, epsilon = 0.01);
    }

    #[test]
    fn test_s_value() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let wald = test_wald_nb(&fit, "grouptreatment").unwrap();

        let shrunk = shrink_lfc(&wald, &ShrinkageConfig::default());

        // S-value should be positive (or infinity for p=0)
        for est in &shrunk.estimates {
            if est.s_value.is_finite() {
                assert!(est.s_value >= 0.0, "S-value should be non-negative");
            }
        }
    }

    // Helper to create sparse test data for hurdle model
    fn create_hurdle_test_counts() -> CountMatrix {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "feature_id\t{}",
            (1..=20).map(|i| format!("S{}", i)).collect::<Vec<_>>().join("\t")
        )
        .unwrap();

        // Feature with no effect (sparse)
        let no_effect: Vec<String> = (1..=20)
            .map(|i| if i % 3 == 0 { "0" } else { "50" }.to_string())
            .collect();
        writeln!(file, "no_effect\t{}", no_effect.join("\t")).unwrap();

        // Feature with strong effect
        let strong: Vec<String> = (1..=20)
            .map(|i| {
                if i <= 10 {
                    if i % 2 == 0 { "0" } else { "10" }
                } else {
                    "80"
                }
                .to_string()
            })
            .collect();
        writeln!(file, "strong_effect\t{}", strong.join("\t")).unwrap();

        // Dense feature
        let dense: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "30".to_string() } else { "60".to_string() })
            .collect();
        writeln!(file, "dense\t{}", dense.join("\t")).unwrap();

        file.flush().unwrap();
        CountMatrix::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_shrink_lfc_hurdle_count() {
        use crate::model::model_hurdle;

        let metadata = create_test_metadata();
        let counts = create_hurdle_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();
        let shrunk = shrink_lfc_hurdle(
            &fit,
            "grouptreatment",
            HurdleComponent::Count,
            &ShrinkageConfig::default(),
        );

        assert!(shrunk.is_some());
        let shrunk = shrunk.unwrap();
        assert_eq!(shrunk.len(), 3);
        assert!(shrunk.coefficient.contains("count"));
    }

    #[test]
    fn test_shrink_lfc_hurdle_binary() {
        use crate::model::model_hurdle;

        let metadata = create_test_metadata();
        let counts = create_hurdle_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();
        let shrunk = shrink_lfc_hurdle(
            &fit,
            "grouptreatment",
            HurdleComponent::Binary,
            &ShrinkageConfig::default(),
        );

        assert!(shrunk.is_some());
        let shrunk = shrunk.unwrap();
        assert_eq!(shrunk.len(), 3);
        assert!(shrunk.coefficient.contains("binary"));
    }

    #[test]
    fn test_shrink_lfc_hurdle_invalid_coefficient() {
        use crate::model::model_hurdle;

        let metadata = create_test_metadata();
        let counts = create_hurdle_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();

        // Test invalid coefficient for count component
        let shrunk = shrink_lfc_hurdle(
            &fit,
            "nonexistent",
            HurdleComponent::Count,
            &ShrinkageConfig::default(),
        );
        assert!(shrunk.is_none());

        // Test invalid coefficient for binary component
        let shrunk = shrink_lfc_hurdle(
            &fit,
            "nonexistent",
            HurdleComponent::Binary,
            &ShrinkageConfig::default(),
        );
        assert!(shrunk.is_none());
    }

    #[test]
    fn test_shrink_lfc_hurdle_confidence_intervals() {
        use crate::model::model_hurdle;

        let metadata = create_test_metadata();
        let counts = create_hurdle_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_hurdle(&counts, &design).unwrap();
        let shrunk = shrink_lfc_hurdle(
            &fit,
            "grouptreatment",
            HurdleComponent::Count,
            &ShrinkageConfig::default(),
        )
        .unwrap();

        for est in &shrunk.estimates {
            if est.lfc_shrunk.is_finite() {
                assert!(est.ci_lower < est.lfc_shrunk);
                assert!(est.ci_upper > est.lfc_shrunk);
                assert!(est.ci_lower < est.ci_upper);
            }
        }
    }

    #[test]
    fn test_hurdle_component_default() {
        assert_eq!(HurdleComponent::default(), HurdleComponent::Count);
    }
}
