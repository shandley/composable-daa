//! Spike-in evaluation: assess pipeline performance on spiked data.

use crate::data::DaResultSet;
use crate::spike::types::SpikeSpec;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Metrics for a specific prevalence tier.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TierMetrics {
    /// Number of spiked features in this tier.
    pub n_spiked: usize,
    /// Number of true positives (correctly detected spikes).
    pub true_positives: usize,
    /// Number of false negatives (missed spikes).
    pub false_negatives: usize,
    /// Sensitivity (TP / (TP + FN)).
    pub sensitivity: f64,
}

impl TierMetrics {
    fn compute_sensitivity(&mut self) {
        let total = self.true_positives + self.false_negatives;
        self.sensitivity = if total > 0 {
            self.true_positives as f64 / total as f64
        } else {
            0.0
        };
    }
}

/// Comprehensive evaluation of spike-in detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikeEvaluation {
    // Detection metrics
    /// Number of true positives (correctly detected spikes).
    pub true_positives: usize,
    /// Number of false positives (non-spiked features called significant).
    pub false_positives: usize,
    /// Number of false negatives (spiked features not detected).
    pub false_negatives: usize,
    /// Number of true negatives (non-spiked features correctly not called).
    pub true_negatives: usize,

    // Derived rates
    /// Sensitivity = TP / (TP + FN).
    pub sensitivity: f64,
    /// Specificity = TN / (TN + FP).
    pub specificity: f64,
    /// Precision = TP / (TP + FP).
    pub precision: f64,
    /// False discovery rate = FP / (TP + FP).
    pub fdr: f64,
    /// F1 score = 2 * precision * recall / (precision + recall).
    pub f1_score: f64,

    // Effect size recovery
    /// Correlation between estimated and true effect sizes.
    pub effect_correlation: f64,
    /// Mean bias (estimated - true).
    pub effect_bias: f64,
    /// Mean absolute error of effect estimates.
    pub effect_mae: f64,

    // By prevalence tier
    /// Metrics broken down by prevalence tier.
    pub by_tier: HashMap<String, TierMetrics>,

    // Metadata
    /// Number of spiked features.
    pub n_spiked: usize,
    /// Total number of features tested.
    pub n_tested: usize,
    /// FDR threshold used for calling significance.
    pub fdr_threshold: f64,
}

impl SpikeEvaluation {
    /// Check if this represents good performance.
    pub fn is_good(&self) -> bool {
        self.sensitivity >= 0.7 && self.fdr <= 0.15
    }

    /// Check if FDR is controlled.
    pub fn fdr_controlled(&self) -> bool {
        self.fdr <= self.fdr_threshold * 1.5 // Allow 50% slack
    }
}

impl std::fmt::Display for SpikeEvaluation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Spike-in Evaluation (n={} spiked, {} tested)", self.n_spiked, self.n_tested)?;
        writeln!(f, "  Detection:")?;
        writeln!(f, "    TP: {}, FP: {}, FN: {}, TN: {}",
                 self.true_positives, self.false_positives,
                 self.false_negatives, self.true_negatives)?;
        writeln!(f, "  Rates:")?;
        writeln!(f, "    Sensitivity: {:.1}%", self.sensitivity * 100.0)?;
        writeln!(f, "    Specificity: {:.1}%", self.specificity * 100.0)?;
        writeln!(f, "    Precision:   {:.1}%", self.precision * 100.0)?;
        writeln!(f, "    FDR:         {:.1}%", self.fdr * 100.0)?;
        writeln!(f, "    F1:          {:.3}", self.f1_score)?;
        writeln!(f, "  Effect size recovery:")?;
        writeln!(f, "    Correlation: {:.3}", self.effect_correlation)?;
        writeln!(f, "    Bias:        {:.3}", self.effect_bias)?;
        writeln!(f, "    MAE:         {:.3}", self.effect_mae)?;
        if !self.by_tier.is_empty() {
            writeln!(f, "  By tier:")?;
            for (tier, metrics) in &self.by_tier {
                writeln!(f, "    {}: {:.1}% sensitivity ({}/{})",
                         tier, metrics.sensitivity * 100.0,
                         metrics.true_positives, metrics.n_spiked)?;
            }
        }
        Ok(())
    }
}

/// Evaluate spike-in detection performance.
///
/// Compares the results from a differential abundance analysis to the
/// known spike-in specification to compute detection metrics.
///
/// # Arguments
/// * `results` - DA results from running the pipeline on spiked data
/// * `spike_spec` - Specification of what was spiked
/// * `fdr_threshold` - Threshold for calling significance (e.g., 0.05)
///
/// # Returns
/// SpikeEvaluation with comprehensive metrics.
pub fn evaluate_spikes(
    results: &DaResultSet,
    spike_spec: &SpikeSpec,
    fdr_threshold: f64,
) -> SpikeEvaluation {
    let spiked_set = spike_spec.spiked_set();
    let n_spiked = spike_spec.n_spiked();
    let n_tested = results.len();

    // Classify each result
    let mut true_positives = 0usize;
    let mut false_positives = 0usize;
    let mut false_negatives = 0usize;
    let mut true_negatives = 0usize;

    // For effect size comparison
    let mut estimated_effects: Vec<f64> = Vec::new();
    let mut true_effects: Vec<f64> = Vec::new();

    // By-tier tracking
    let mut tier_metrics: HashMap<String, TierMetrics> = HashMap::new();

    // Track which spiked features were tested
    let mut spiked_found: HashSet<&str> = HashSet::new();

    for result in &results.results {
        let is_spiked = spiked_set.contains(result.feature_id.as_str());
        let is_significant = result.q_value < fdr_threshold;

        // Get prevalence tier for this feature
        let tier_name = result.prevalence_tier.name().to_string();
        let tier_entry = tier_metrics.entry(tier_name.clone()).or_default();

        if is_spiked {
            spiked_found.insert(result.feature_id.as_str());
            tier_entry.n_spiked += 1;

            if is_significant {
                true_positives += 1;
                tier_entry.true_positives += 1;
            } else {
                false_negatives += 1;
                tier_entry.false_negatives += 1;
            }

            // Collect effect sizes for correlation
            if let Some(true_effect) = spike_spec.effect_size(&result.feature_id) {
                estimated_effects.push(result.estimate);
                true_effects.push(true_effect.ln()); // Convert fold change to log scale
            }
        } else {
            if is_significant {
                false_positives += 1;
            } else {
                true_negatives += 1;
            }
        }
    }

    // Check for spiked features not in results (filtered out?)
    for spiked_id in &spike_spec.spiked_features {
        if !spiked_found.contains(spiked_id.as_str()) {
            false_negatives += 1;
        }
    }

    // Compute rates
    let sensitivity = if n_spiked > 0 {
        true_positives as f64 / n_spiked as f64
    } else {
        0.0
    };

    let n_non_spiked = n_tested.saturating_sub(spiked_found.len());
    let specificity = if n_non_spiked > 0 {
        true_negatives as f64 / n_non_spiked as f64
    } else {
        1.0
    };

    let n_called = true_positives + false_positives;
    let precision = if n_called > 0 {
        true_positives as f64 / n_called as f64
    } else {
        1.0
    };

    let fdr = if n_called > 0 {
        false_positives as f64 / n_called as f64
    } else {
        0.0
    };

    let f1_score = if precision + sensitivity > 0.0 {
        2.0 * precision * sensitivity / (precision + sensitivity)
    } else {
        0.0
    };

    // Effect size metrics
    let (effect_correlation, effect_bias, effect_mae) = if estimated_effects.len() >= 2 {
        compute_effect_metrics(&estimated_effects, &true_effects)
    } else {
        (f64::NAN, f64::NAN, f64::NAN)
    };

    // Compute per-tier sensitivity
    for metrics in tier_metrics.values_mut() {
        metrics.compute_sensitivity();
    }

    SpikeEvaluation {
        true_positives,
        false_positives,
        false_negatives,
        true_negatives,
        sensitivity,
        specificity,
        precision,
        fdr,
        f1_score,
        effect_correlation,
        effect_bias,
        effect_mae,
        by_tier: tier_metrics,
        n_spiked,
        n_tested,
        fdr_threshold,
    }
}

/// Compute effect size recovery metrics.
fn compute_effect_metrics(estimated: &[f64], true_vals: &[f64]) -> (f64, f64, f64) {
    let n = estimated.len() as f64;
    if n < 2.0 {
        return (f64::NAN, f64::NAN, f64::NAN);
    }

    // Mean bias
    let bias: f64 = estimated
        .iter()
        .zip(true_vals.iter())
        .map(|(e, t)| e - t)
        .sum::<f64>() / n;

    // MAE
    let mae: f64 = estimated
        .iter()
        .zip(true_vals.iter())
        .map(|(e, t)| (e - t).abs())
        .sum::<f64>() / n;

    // Pearson correlation
    let mean_est = estimated.iter().sum::<f64>() / n;
    let mean_true = true_vals.iter().sum::<f64>() / n;

    let mut cov = 0.0;
    let mut var_est = 0.0;
    let mut var_true = 0.0;

    for (e, t) in estimated.iter().zip(true_vals.iter()) {
        let de = e - mean_est;
        let dt = t - mean_true;
        cov += de * dt;
        var_est += de * de;
        var_true += dt * dt;
    }

    let correlation = if var_est > 0.0 && var_true > 0.0 {
        cov / (var_est.sqrt() * var_true.sqrt())
    } else {
        f64::NAN
    };

    (correlation, bias, mae)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DaResult, PrevalenceTier};
    use crate::spike::types::SpikeType;

    fn create_test_results() -> DaResultSet {
        // 10 features: 3 spiked (feat_0, feat_1, feat_2), 7 not spiked
        let mut results = Vec::new();

        // feat_0: spiked, detected (TP)
        results.push(DaResult::new(
            "feat_0".into(), "group".into(),
            1.0, 0.2, 5.0, 0.001, 0.01, 0.8, 100.0,
        ));

        // feat_1: spiked, detected (TP)
        results.push(DaResult::new(
            "feat_1".into(), "group".into(),
            0.8, 0.2, 4.0, 0.005, 0.02, 0.7, 80.0,
        ));

        // feat_2: spiked, NOT detected (FN)
        results.push(DaResult::new(
            "feat_2".into(), "group".into(),
            0.3, 0.3, 1.0, 0.2, 0.3, 0.5, 50.0,
        ));

        // feat_3-6: not spiked, not detected (TN)
        for i in 3..7 {
            results.push(DaResult::new(
                format!("feat_{}", i), "group".into(),
                0.1, 0.2, 0.5, 0.5, 0.6, 0.6, 60.0,
            ));
        }

        // feat_7: not spiked, detected (FP)
        results.push(DaResult::new(
            "feat_7".into(), "group".into(),
            0.7, 0.2, 3.5, 0.01, 0.04, 0.7, 70.0,
        ));

        // feat_8-9: not spiked, not detected (TN)
        for i in 8..10 {
            results.push(DaResult::new(
                format!("feat_{}", i), "group".into(),
                0.05, 0.2, 0.25, 0.6, 0.7, 0.4, 40.0,
            ));
        }

        DaResultSet::new("test".into(), results)
    }

    fn create_test_spec() -> SpikeSpec {
        SpikeSpec::new(
            SpikeType::Abundance,
            vec!["feat_0".into(), "feat_1".into(), "feat_2".into()],
            vec![2.7, 2.2, 1.3], // fold changes (e^1.0, e^0.8, e^0.3)
            "treatment".into(),
            vec![0.8, 0.7, 0.5],
            42,
        )
    }

    #[test]
    fn test_evaluate_basic() {
        let results = create_test_results();
        let spec = create_test_spec();

        let eval = evaluate_spikes(&results, &spec, 0.05);

        assert_eq!(eval.true_positives, 2);  // feat_0, feat_1
        assert_eq!(eval.false_positives, 1); // feat_7
        assert_eq!(eval.false_negatives, 1); // feat_2
        assert_eq!(eval.true_negatives, 6);  // feat_3-6, feat_8-9
    }

    #[test]
    fn test_evaluate_rates() {
        let results = create_test_results();
        let spec = create_test_spec();

        let eval = evaluate_spikes(&results, &spec, 0.05);

        // Sensitivity = 2/3 ≈ 0.667
        assert!((eval.sensitivity - 2.0/3.0).abs() < 0.01);

        // Precision = 2/3 ≈ 0.667
        assert!((eval.precision - 2.0/3.0).abs() < 0.01);

        // FDR = 1/3 ≈ 0.333
        assert!((eval.fdr - 1.0/3.0).abs() < 0.01);
    }

    #[test]
    fn test_perfect_detection() {
        // All spiked detected, no false positives
        let mut results = Vec::new();
        results.push(DaResult::new(
            "feat_0".into(), "group".into(),
            1.0, 0.1, 10.0, 0.0001, 0.001, 0.8, 100.0,
        ));
        results.push(DaResult::new(
            "feat_1".into(), "group".into(),
            0.1, 0.1, 1.0, 0.5, 0.6, 0.8, 100.0,
        ));

        let result_set = DaResultSet::new("test".into(), results);
        let spec = SpikeSpec::new(
            SpikeType::Abundance,
            vec!["feat_0".into()],
            vec![2.7],
            "treatment".into(),
            vec![0.8],
            42,
        );

        let eval = evaluate_spikes(&result_set, &spec, 0.05);

        assert_eq!(eval.true_positives, 1);
        assert_eq!(eval.false_positives, 0);
        assert_eq!(eval.false_negatives, 0);
        assert!((eval.sensitivity - 1.0).abs() < 0.001);
        assert!((eval.fdr - 0.0).abs() < 0.001);
    }
}
