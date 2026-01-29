//! Model comparison using information criteria.
//!
//! Provides AIC (Akaike Information Criterion) and BIC (Bayesian Information Criterion)
//! for comparing model fits and selecting the most appropriate model.
//!
//! - AIC = -2 * log_likelihood + 2 * k
//! - BIC = -2 * log_likelihood + k * log(n)
//!
//! Where k = number of parameters, n = number of observations.
//! Lower values indicate better models.

use crate::model::nb::NbFit;
use crate::model::zinb::ZinbFit;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Information criteria for a single feature's model fit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModelCriteria {
    /// Feature identifier.
    pub feature_id: String,
    /// Model name/type.
    pub model: String,
    /// Log-likelihood.
    pub log_likelihood: f64,
    /// Number of parameters.
    pub n_params: usize,
    /// Number of observations (samples).
    pub n_obs: usize,
    /// Akaike Information Criterion.
    pub aic: f64,
    /// Bayesian Information Criterion.
    pub bic: f64,
    /// AIC corrected for small samples (AICc).
    pub aicc: f64,
}

impl ModelCriteria {
    /// Create new model criteria from components.
    pub fn new(
        feature_id: &str,
        model: &str,
        log_likelihood: f64,
        n_params: usize,
        n_obs: usize,
    ) -> Self {
        let k = n_params as f64;
        let n = n_obs as f64;

        // AIC = -2 * LL + 2k
        let aic = -2.0 * log_likelihood + 2.0 * k;

        // BIC = -2 * LL + k * log(n)
        let bic = -2.0 * log_likelihood + k * n.ln();

        // AICc = AIC + (2k² + 2k) / (n - k - 1)
        // Corrected AIC for small sample sizes
        let aicc = if n > k + 1.0 {
            aic + (2.0 * k * k + 2.0 * k) / (n - k - 1.0)
        } else {
            f64::INFINITY
        };

        Self {
            feature_id: feature_id.to_string(),
            model: model.to_string(),
            log_likelihood,
            n_params,
            n_obs,
            aic,
            bic,
            aicc,
        }
    }
}

/// Comparison result for a single feature between two models.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureComparison {
    /// Feature identifier.
    pub feature_id: String,
    /// First model criteria.
    pub model1: ModelCriteria,
    /// Second model criteria.
    pub model2: ModelCriteria,
    /// Difference in AIC (model1 - model2). Negative = model1 better.
    pub delta_aic: f64,
    /// Difference in BIC (model1 - model2). Negative = model1 better.
    pub delta_bic: f64,
    /// Likelihood ratio statistic (if models are nested).
    pub lr_statistic: f64,
    /// Preferred model based on AIC.
    pub preferred_aic: String,
    /// Preferred model based on BIC.
    pub preferred_bic: String,
}

impl FeatureComparison {
    /// Create comparison from two model criteria.
    pub fn new(model1: ModelCriteria, model2: ModelCriteria) -> Self {
        let delta_aic = model1.aic - model2.aic;
        let delta_bic = model1.bic - model2.bic;

        // LR statistic = 2 * (LL_full - LL_reduced)
        let lr_statistic = 2.0 * (model1.log_likelihood - model2.log_likelihood);

        let preferred_aic = if delta_aic < 0.0 {
            model1.model.clone()
        } else {
            model2.model.clone()
        };

        let preferred_bic = if delta_bic < 0.0 {
            model1.model.clone()
        } else {
            model2.model.clone()
        };

        Self {
            feature_id: model1.feature_id.clone(),
            model1,
            model2,
            delta_aic,
            delta_bic,
            lr_statistic,
            preferred_aic,
            preferred_bic,
        }
    }

    /// Get the strength of evidence for the preferred model (AIC-based).
    ///
    /// Based on Burnham & Anderson guidelines:
    /// - |ΔAIC| < 2: Substantial support for both models
    /// - 2 ≤ |ΔAIC| < 4: Considerably less support for worse model
    /// - 4 ≤ |ΔAIC| < 7: Much less support for worse model
    /// - |ΔAIC| ≥ 10: Essentially no support for worse model
    pub fn evidence_strength_aic(&self) -> EvidenceStrength {
        let delta = self.delta_aic.abs();
        if delta < 2.0 {
            EvidenceStrength::Weak
        } else if delta < 4.0 {
            EvidenceStrength::Moderate
        } else if delta < 7.0 {
            EvidenceStrength::Strong
        } else {
            EvidenceStrength::VeryStrong
        }
    }

    /// Get the strength of evidence for the preferred model (BIC-based).
    ///
    /// Based on Kass & Raftery guidelines for Bayes factors:
    /// - |ΔBIC| < 2: Not worth more than a bare mention
    /// - 2 ≤ |ΔBIC| < 6: Positive evidence
    /// - 6 ≤ |ΔBIC| < 10: Strong evidence
    /// - |ΔBIC| ≥ 10: Very strong evidence
    pub fn evidence_strength_bic(&self) -> EvidenceStrength {
        let delta = self.delta_bic.abs();
        if delta < 2.0 {
            EvidenceStrength::Weak
        } else if delta < 6.0 {
            EvidenceStrength::Moderate
        } else if delta < 10.0 {
            EvidenceStrength::Strong
        } else {
            EvidenceStrength::VeryStrong
        }
    }
}

/// Strength of evidence for model preference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EvidenceStrength {
    /// Little difference between models.
    Weak,
    /// Moderate preference for one model.
    Moderate,
    /// Strong preference for one model.
    Strong,
    /// Very strong preference for one model.
    VeryStrong,
}

impl fmt::Display for EvidenceStrength {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EvidenceStrength::Weak => write!(f, "weak"),
            EvidenceStrength::Moderate => write!(f, "moderate"),
            EvidenceStrength::Strong => write!(f, "strong"),
            EvidenceStrength::VeryStrong => write!(f, "very strong"),
        }
    }
}

/// Summary of model comparison across all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComparisonSummary {
    /// Model 1 name.
    pub model1_name: String,
    /// Model 2 name.
    pub model2_name: String,
    /// Number of features.
    pub n_features: usize,
    /// Number preferring model 1 by AIC.
    pub n_prefer_model1_aic: usize,
    /// Number preferring model 2 by AIC.
    pub n_prefer_model2_aic: usize,
    /// Number preferring model 1 by BIC.
    pub n_prefer_model1_bic: usize,
    /// Number preferring model 2 by BIC.
    pub n_prefer_model2_bic: usize,
    /// Mean ΔAIC (model1 - model2).
    pub mean_delta_aic: f64,
    /// Mean ΔBIC (model1 - model2).
    pub mean_delta_bic: f64,
    /// Median ΔAIC.
    pub median_delta_aic: f64,
    /// Median ΔBIC.
    pub median_delta_bic: f64,
}

impl fmt::Display for ComparisonSummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Model Comparison Summary")?;
        writeln!(f, "========================")?;
        writeln!(f, "Models: {} vs {}", self.model1_name, self.model2_name)?;
        writeln!(f, "Features: {}", self.n_features)?;
        writeln!(f)?;
        writeln!(f, "AIC Preference:")?;
        writeln!(f, "  {}: {} features ({:.1}%)",
            self.model1_name,
            self.n_prefer_model1_aic,
            100.0 * self.n_prefer_model1_aic as f64 / self.n_features as f64
        )?;
        writeln!(f, "  {}: {} features ({:.1}%)",
            self.model2_name,
            self.n_prefer_model2_aic,
            100.0 * self.n_prefer_model2_aic as f64 / self.n_features as f64
        )?;
        writeln!(f, "  Mean ΔAIC: {:.2}", self.mean_delta_aic)?;
        writeln!(f)?;
        writeln!(f, "BIC Preference:")?;
        writeln!(f, "  {}: {} features ({:.1}%)",
            self.model1_name,
            self.n_prefer_model1_bic,
            100.0 * self.n_prefer_model1_bic as f64 / self.n_features as f64
        )?;
        writeln!(f, "  {}: {} features ({:.1}%)",
            self.model2_name,
            self.n_prefer_model2_bic,
            100.0 * self.n_prefer_model2_bic as f64 / self.n_features as f64
        )?;
        writeln!(f, "  Mean ΔBIC: {:.2}", self.mean_delta_bic)?;
        Ok(())
    }
}

/// Results of comparing two models across all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModelComparisonResult {
    /// Individual feature comparisons.
    pub comparisons: Vec<FeatureComparison>,
    /// Summary statistics.
    pub summary: ComparisonSummary,
}

impl ModelComparisonResult {
    /// Get comparison for a specific feature.
    pub fn get_feature(&self, feature_id: &str) -> Option<&FeatureComparison> {
        self.comparisons.iter().find(|c| c.feature_id == feature_id)
    }

    /// Get features where model 1 is preferred by AIC.
    pub fn features_preferring_model1_aic(&self) -> Vec<&str> {
        self.comparisons
            .iter()
            .filter(|c| c.delta_aic < 0.0)
            .map(|c| c.feature_id.as_str())
            .collect()
    }

    /// Get features where model 2 is preferred by AIC.
    pub fn features_preferring_model2_aic(&self) -> Vec<&str> {
        self.comparisons
            .iter()
            .filter(|c| c.delta_aic >= 0.0)
            .map(|c| c.feature_id.as_str())
            .collect()
    }

    /// Get features with strong evidence for one model (AIC-based).
    pub fn features_with_strong_preference(&self) -> Vec<&FeatureComparison> {
        self.comparisons
            .iter()
            .filter(|c| matches!(
                c.evidence_strength_aic(),
                EvidenceStrength::Strong | EvidenceStrength::VeryStrong
            ))
            .collect()
    }
}

impl fmt::Display for ModelComparisonResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.summary)
    }
}

/// Compute model criteria for a negative binomial fit.
///
/// Parameters counted:
/// - p regression coefficients
/// - 1 dispersion parameter
pub fn criteria_nb(fit: &NbFit) -> Vec<ModelCriteria> {
    let n_obs = fit.n_samples;
    let n_coef = fit.coefficient_names.len();

    fit.fits
        .iter()
        .map(|f| {
            // NB has: p coefficients + 1 dispersion
            let n_params = n_coef + 1;
            ModelCriteria::new(
                &f.feature_id,
                "NB",
                f.log_likelihood,
                n_params,
                n_obs,
            )
        })
        .collect()
}

/// Compute model criteria for a zero-inflated negative binomial fit.
///
/// Parameters counted:
/// - p count model coefficients
/// - p_zi zero-inflation coefficients
/// - 1 dispersion parameter
pub fn criteria_zinb(fit: &ZinbFit) -> Vec<ModelCriteria> {
    let n_obs = fit.n_samples;
    let n_coef = fit.coefficient_names.len();
    let n_zi_coef = fit.zi_coefficient_names.len();

    fit.fits
        .iter()
        .map(|f| {
            // ZINB has: p + p_zi coefficients + 1 dispersion
            let n_params = n_coef + n_zi_coef + 1;
            ModelCriteria::new(
                &f.feature_id,
                "ZINB",
                f.log_likelihood,
                n_params,
                n_obs,
            )
        })
        .collect()
}

/// Compare two sets of model criteria.
///
/// Both sets must have the same features in the same order.
pub fn compare_criteria(
    criteria1: &[ModelCriteria],
    criteria2: &[ModelCriteria],
) -> ModelComparisonResult {
    let comparisons: Vec<FeatureComparison> = criteria1
        .iter()
        .zip(criteria2.iter())
        .map(|(c1, c2)| FeatureComparison::new(c1.clone(), c2.clone()))
        .collect();

    let summary = summarize_comparisons(&comparisons);

    ModelComparisonResult {
        comparisons,
        summary,
    }
}

/// Compare NB and ZINB fits.
///
/// Returns comparison with NB as model1, ZINB as model2.
/// Negative delta values indicate NB is preferred.
pub fn compare_nb_zinb(nb_fit: &NbFit, zinb_fit: &ZinbFit) -> ModelComparisonResult {
    let nb_criteria = criteria_nb(nb_fit);
    let zinb_criteria = criteria_zinb(zinb_fit);
    compare_criteria(&nb_criteria, &zinb_criteria)
}

/// Summarize feature comparisons.
fn summarize_comparisons(comparisons: &[FeatureComparison]) -> ComparisonSummary {
    let n = comparisons.len();
    if n == 0 {
        return ComparisonSummary {
            model1_name: String::new(),
            model2_name: String::new(),
            n_features: 0,
            n_prefer_model1_aic: 0,
            n_prefer_model2_aic: 0,
            n_prefer_model1_bic: 0,
            n_prefer_model2_bic: 0,
            mean_delta_aic: 0.0,
            mean_delta_bic: 0.0,
            median_delta_aic: 0.0,
            median_delta_bic: 0.0,
        };
    }

    let model1_name = comparisons[0].model1.model.clone();
    let model2_name = comparisons[0].model2.model.clone();

    let n_prefer_model1_aic = comparisons.iter().filter(|c| c.delta_aic < 0.0).count();
    let n_prefer_model2_aic = n - n_prefer_model1_aic;

    let n_prefer_model1_bic = comparisons.iter().filter(|c| c.delta_bic < 0.0).count();
    let n_prefer_model2_bic = n - n_prefer_model1_bic;

    let delta_aics: Vec<f64> = comparisons.iter().map(|c| c.delta_aic).collect();
    let delta_bics: Vec<f64> = comparisons.iter().map(|c| c.delta_bic).collect();

    let mean_delta_aic = delta_aics.iter().sum::<f64>() / n as f64;
    let mean_delta_bic = delta_bics.iter().sum::<f64>() / n as f64;

    let median_delta_aic = median(&delta_aics);
    let median_delta_bic = median(&delta_bics);

    ComparisonSummary {
        model1_name,
        model2_name,
        n_features: n,
        n_prefer_model1_aic,
        n_prefer_model2_aic,
        n_prefer_model1_bic,
        n_prefer_model2_bic,
        mean_delta_aic,
        mean_delta_bic,
        median_delta_aic,
        median_delta_bic,
    }
}

/// Compute median of a slice.
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = sorted.len() / 2;
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) / 2.0
    } else {
        sorted[mid]
    }
}

/// Select the best model per feature based on a criterion.
///
/// Returns a vector of (feature_id, preferred_model_name) tuples.
pub fn select_best_model(
    comparisons: &[FeatureComparison],
    criterion: SelectionCriterion,
) -> Vec<(String, String)> {
    comparisons
        .iter()
        .map(|c| {
            let preferred = match criterion {
                SelectionCriterion::Aic => c.preferred_aic.clone(),
                SelectionCriterion::Bic => c.preferred_bic.clone(),
                SelectionCriterion::Aicc => {
                    if c.model1.aicc < c.model2.aicc {
                        c.model1.model.clone()
                    } else {
                        c.model2.model.clone()
                    }
                }
            };
            (c.feature_id.clone(), preferred)
        })
        .collect()
}

/// Criterion for model selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionCriterion {
    /// Akaike Information Criterion.
    Aic,
    /// Bayesian Information Criterion.
    Bic,
    /// Corrected AIC for small samples.
    Aicc,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{CountMatrix, DesignMatrix, Formula, Metadata};
    use crate::model::{model_nb, model_zinb};
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

        // Feature with no zeros (NB should be fine)
        writeln!(
            file,
            "no_zeros\t{}",
            (1..=20).map(|i| format!("{}", 10 + i)).collect::<Vec<_>>().join("\t")
        )
        .unwrap();

        // Feature with some zeros (might benefit from ZINB)
        let some_zeros: Vec<String> = (1..=20)
            .map(|i| if i % 4 == 0 { "0".to_string() } else { format!("{}", 15 + i) })
            .collect();
        writeln!(file, "some_zeros\t{}", some_zeros.join("\t")).unwrap();

        // Feature with many zeros (ZINB likely better)
        let many_zeros: Vec<String> = (1..=20)
            .map(|i| if i <= 10 { "0".to_string() } else { format!("{}", 20 + i) })
            .collect();
        writeln!(file, "many_zeros\t{}", many_zeros.join("\t")).unwrap();

        file.flush().unwrap();
        CountMatrix::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_model_criteria_creation() {
        let criteria = ModelCriteria::new("feature1", "NB", -100.0, 3, 20);

        assert_eq!(criteria.feature_id, "feature1");
        assert_eq!(criteria.model, "NB");
        assert_eq!(criteria.log_likelihood, -100.0);
        assert_eq!(criteria.n_params, 3);
        assert_eq!(criteria.n_obs, 20);

        // AIC = -2 * (-100) + 2 * 3 = 200 + 6 = 206
        assert!((criteria.aic - 206.0).abs() < 0.001);

        // BIC = -2 * (-100) + 3 * ln(20) = 200 + 3 * 2.996 ≈ 208.99
        assert!((criteria.bic - 208.99).abs() < 0.1);
    }

    #[test]
    fn test_feature_comparison() {
        let model1 = ModelCriteria::new("f1", "NB", -100.0, 3, 20);
        let model2 = ModelCriteria::new("f1", "ZINB", -95.0, 5, 20);

        let comp = FeatureComparison::new(model1, model2);

        assert_eq!(comp.feature_id, "f1");
        // ZINB has higher LL, but more parameters
        // NB: AIC = 206, ZINB: AIC = -2*(-95) + 2*5 = 190 + 10 = 200
        // delta_aic = 206 - 200 = 6 (positive means ZINB better)
        assert!(comp.delta_aic > 0.0);
        assert_eq!(comp.preferred_aic, "ZINB");
    }

    #[test]
    fn test_evidence_strength() {
        let model1 = ModelCriteria::new("f1", "NB", -100.0, 3, 20);
        let model2 = ModelCriteria::new("f1", "ZINB", -100.0, 5, 20);

        let comp = FeatureComparison::new(model1, model2);

        // Small difference should be weak evidence
        // With same LL, NB has lower AIC (fewer params)
        // delta_aic = 206 - 210 = -4
        let strength = comp.evidence_strength_aic();
        assert!(matches!(strength, EvidenceStrength::Moderate | EvidenceStrength::Strong));
    }

    #[test]
    fn test_criteria_nb() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_nb(&counts, &design).unwrap();
        let criteria = criteria_nb(&fit);

        assert_eq!(criteria.len(), 3);
        for c in &criteria {
            assert_eq!(c.model, "NB");
            assert_eq!(c.n_params, 3); // intercept + group + dispersion
            assert_eq!(c.n_obs, 20);
            assert!(c.aic.is_finite());
            assert!(c.bic.is_finite());
        }
    }

    #[test]
    fn test_criteria_zinb() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let fit = model_zinb(&counts, &design).unwrap();
        let criteria = criteria_zinb(&fit);

        assert_eq!(criteria.len(), 3);
        for c in &criteria {
            assert_eq!(c.model, "ZINB");
            // ZINB: 2 count coefs + 2 ZI coefs + 1 dispersion = 5
            assert_eq!(c.n_params, 5);
            assert_eq!(c.n_obs, 20);
            assert!(c.aic.is_finite());
            assert!(c.bic.is_finite());
        }
    }

    #[test]
    fn test_compare_nb_zinb() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let nb_fit = model_nb(&counts, &design).unwrap();
        let zinb_fit = model_zinb(&counts, &design).unwrap();

        let comparison = compare_nb_zinb(&nb_fit, &zinb_fit);

        assert_eq!(comparison.comparisons.len(), 3);
        assert_eq!(comparison.summary.n_features, 3);
        assert_eq!(comparison.summary.model1_name, "NB");
        assert_eq!(comparison.summary.model2_name, "ZINB");

        // Total preferences should equal number of features
        assert_eq!(
            comparison.summary.n_prefer_model1_aic + comparison.summary.n_prefer_model2_aic,
            3
        );
    }

    #[test]
    fn test_comparison_summary_display() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let nb_fit = model_nb(&counts, &design).unwrap();
        let zinb_fit = model_zinb(&counts, &design).unwrap();

        let comparison = compare_nb_zinb(&nb_fit, &zinb_fit);

        let display = format!("{}", comparison);
        assert!(display.contains("Model Comparison Summary"));
        assert!(display.contains("NB"));
        assert!(display.contains("ZINB"));
    }

    #[test]
    fn test_select_best_model() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let nb_fit = model_nb(&counts, &design).unwrap();
        let zinb_fit = model_zinb(&counts, &design).unwrap();

        let comparison = compare_nb_zinb(&nb_fit, &zinb_fit);
        let selections = select_best_model(&comparison.comparisons, SelectionCriterion::Aic);

        assert_eq!(selections.len(), 3);
        for (feature_id, model) in &selections {
            assert!(!feature_id.is_empty());
            assert!(model == "NB" || model == "ZINB");
        }
    }

    #[test]
    fn test_get_feature_comparison() {
        let metadata = create_test_metadata();
        let counts = create_test_counts();
        let formula = Formula::parse("~ group").unwrap();
        let design = DesignMatrix::from_formula(&metadata, &formula).unwrap();

        let nb_fit = model_nb(&counts, &design).unwrap();
        let zinb_fit = model_zinb(&counts, &design).unwrap();

        let comparison = compare_nb_zinb(&nb_fit, &zinb_fit);

        let feature = comparison.get_feature("no_zeros");
        assert!(feature.is_some());
        assert_eq!(feature.unwrap().feature_id, "no_zeros");

        let missing = comparison.get_feature("nonexistent");
        assert!(missing.is_none());
    }

    #[test]
    fn test_median() {
        assert_eq!(median(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[5.0]), 5.0);
        assert_eq!(median(&[]), 0.0);
    }
}
