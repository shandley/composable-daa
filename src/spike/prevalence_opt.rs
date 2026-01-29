//! Prevalence threshold optimization via spike-in validation.
//!
//! This module provides tools to empirically determine optimal prevalence
//! filtering thresholds by sweeping across threshold values and measuring
//! detection performance via spike-in analysis.

use crate::data::{CountMatrix, DaResultSet, Metadata};
use crate::error::{DaaError, Result};
use crate::filter::{filter_prevalence_groupwise, filter_prevalence_overall, GroupwiseLogic};
use crate::spike::evaluate::evaluate_spikes;
use crate::spike::types::SpikeSelection;
use crate::spike::{spike_abundance_with_mode, SpikeMode};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

/// Simple LCG random number generator for reproducibility.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1);
        self.state
    }
}

/// Configuration for prevalence threshold optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceOptConfig {
    /// Thresholds to evaluate (e.g., [0.01, 0.02, 0.05, 0.10, 0.20, 0.30]).
    pub thresholds: Vec<f64>,
    /// Fold change for abundance spikes.
    pub fold_change: f64,
    /// Number of features to spike at each threshold.
    pub n_spike: usize,
    /// Number of replicates per threshold for stability.
    pub n_replicates: usize,
    /// Group column in metadata.
    pub group_column: String,
    /// Target group for spiking.
    pub target_group: String,
    /// FDR threshold for significance calls.
    pub fdr_threshold: f64,
    /// Criterion for selecting optimal threshold.
    pub selection_criterion: OptimizationCriterion,
    /// Filtering logic (overall vs. groupwise).
    pub filter_logic: PrevalenceFilterLogic,
    /// Spike mode to use.
    pub spike_mode: SpikeMode,
    /// Random seed.
    pub seed: u64,
}

impl Default for PrevalenceOptConfig {
    fn default() -> Self {
        Self {
            thresholds: vec![0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30],
            fold_change: 2.0,
            n_spike: 20,
            n_replicates: 5,
            group_column: "group".into(),
            target_group: "treatment".into(),
            fdr_threshold: 0.05,
            selection_criterion: OptimizationCriterion::MaxF1,
            filter_logic: PrevalenceFilterLogic::Overall,
            spike_mode: SpikeMode::Compositional,
            seed: 42,
        }
    }
}

impl PrevalenceOptConfig {
    /// Create a quick configuration for faster testing.
    pub fn quick() -> Self {
        Self {
            thresholds: vec![0.05, 0.10, 0.20],
            n_replicates: 2,
            ..Default::default()
        }
    }

    /// Create a thorough configuration for detailed analysis.
    pub fn thorough() -> Self {
        Self {
            thresholds: vec![
                0.01, 0.02, 0.03, 0.05, 0.07, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50,
            ],
            n_replicates: 10,
            ..Default::default()
        }
    }
}

/// How to filter by prevalence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PrevalenceFilterLogic {
    /// Overall prevalence across all samples.
    Overall,
    /// Must pass threshold in any group.
    AnyGroup,
    /// Must pass threshold in all groups.
    AllGroups,
}

/// Criterion for selecting optimal threshold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum OptimizationCriterion {
    /// Maximize F1 score.
    MaxF1,
    /// Maximize sensitivity while keeping FDR below target.
    MaxSensitivityAtFdr(f64),
    /// Maximize sensitivity per feature retained (efficiency).
    MaxEfficiency,
    /// Minimize FDR while maintaining minimum sensitivity.
    MinFdrAtSensitivity(f64),
}

impl Default for OptimizationCriterion {
    fn default() -> Self {
        OptimizationCriterion::MaxF1
    }
}

/// Result from a single replicate at a threshold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThresholdReplicate {
    /// Replicate number.
    pub replicate: usize,
    /// Random seed used.
    pub seed: u64,
    /// Sensitivity (TP / (TP + FN)).
    pub sensitivity: f64,
    /// False discovery rate (FP / (TP + FP)).
    pub fdr: f64,
    /// F1 score (harmonic mean of precision and sensitivity).
    pub f1_score: f64,
    /// Number of true positives.
    pub true_positives: usize,
    /// Number of false positives.
    pub false_positives: usize,
    /// Number of false negatives.
    pub false_negatives: usize,
    /// Number of features that were spiked (after filtering).
    pub n_spiked: usize,
}

/// Aggregated result for a single threshold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThresholdResult {
    /// The prevalence threshold.
    pub threshold: f64,
    /// Number of features retained after filtering.
    pub n_features_retained: usize,
    /// Percentage of original features retained.
    pub pct_features_retained: f64,
    /// Mean sensitivity across replicates.
    pub mean_sensitivity: f64,
    /// Standard deviation of sensitivity.
    pub std_sensitivity: f64,
    /// Mean FDR across replicates.
    pub mean_fdr: f64,
    /// Standard deviation of FDR.
    pub std_fdr: f64,
    /// Mean F1 score across replicates.
    pub mean_f1: f64,
    /// Standard deviation of F1.
    pub std_f1: f64,
    /// Individual replicate results.
    pub replicates: Vec<ThresholdReplicate>,
}

impl ThresholdResult {
    /// Calculate efficiency: sensitivity per feature retained.
    pub fn efficiency(&self) -> f64 {
        if self.n_features_retained == 0 {
            0.0
        } else {
            self.mean_sensitivity / (self.n_features_retained as f64)
        }
    }
}

/// Complete optimization result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceOptResult {
    /// Configuration used.
    pub config: PrevalenceOptConfig,
    /// Results for each threshold.
    pub results: Vec<ThresholdResult>,
    /// Optimal threshold based on selection criterion.
    pub optimal_threshold: f64,
    /// Index of optimal result in results vector.
    pub optimal_index: usize,
    /// Text recommendation.
    pub recommendation: String,
    /// Total runtime in seconds.
    pub runtime_seconds: f64,
}

impl PrevalenceOptResult {
    /// Get the optimal result.
    pub fn optimal_result(&self) -> &ThresholdResult {
        &self.results[self.optimal_index]
    }

    /// Get (threshold, sensitivity) pairs for plotting.
    pub fn sensitivity_curve(&self) -> Vec<(f64, f64)> {
        self.results
            .iter()
            .map(|r| (r.threshold, r.mean_sensitivity))
            .collect()
    }

    /// Get (threshold, fdr) pairs for plotting.
    pub fn fdr_curve(&self) -> Vec<(f64, f64)> {
        self.results.iter().map(|r| (r.threshold, r.mean_fdr)).collect()
    }

    /// Get (threshold, f1) pairs for plotting.
    pub fn f1_curve(&self) -> Vec<(f64, f64)> {
        self.results.iter().map(|r| (r.threshold, r.mean_f1)).collect()
    }

    /// Get (threshold, n_features) pairs for plotting.
    pub fn features_curve(&self) -> Vec<(f64, usize)> {
        self.results
            .iter()
            .map(|r| (r.threshold, r.n_features_retained))
            .collect()
    }

    /// Export as CSV.
    pub fn to_csv(&self) -> String {
        let mut lines = vec![
            "threshold,n_features,pct_features,sensitivity,std_sensitivity,fdr,std_fdr,f1,std_f1"
                .to_string(),
        ];
        for r in &self.results {
            lines.push(format!(
                "{:.3},{},{:.1},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}",
                r.threshold,
                r.n_features_retained,
                r.pct_features_retained * 100.0,
                r.mean_sensitivity,
                r.std_sensitivity,
                r.mean_fdr,
                r.std_fdr,
                r.mean_f1,
                r.std_f1,
            ));
        }
        lines.join("\n")
    }
}

impl fmt::Display for PrevalenceOptResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Prevalence Threshold Optimization")?;
        writeln!(f, "==================================")?;
        writeln!(f)?;
        writeln!(
            f,
            "{:>10}  {:>8}  {:>12}  {:>8}  {:>8}",
            "Threshold", "Features", "Sensitivity", "FDR", "F1"
        )?;
        writeln!(
            f,
            "{:>10}  {:>8}  {:>12}  {:>8}  {:>8}",
            "---------", "--------", "-----------", "-------", "------"
        )?;

        for (i, r) in self.results.iter().enumerate() {
            let marker = if i == self.optimal_index {
                " ← Optimal"
            } else {
                ""
            };
            writeln!(
                f,
                "{:>10.2}  {:>5} ({:>2.0}%)  {:>5.1} ± {:>4.1}%  {:>6.1}%  {:>6.2}{}",
                r.threshold,
                r.n_features_retained,
                r.pct_features_retained * 100.0,
                r.mean_sensitivity * 100.0,
                r.std_sensitivity * 100.0,
                r.mean_fdr * 100.0,
                r.mean_f1,
                marker
            )?;
        }

        writeln!(f)?;
        writeln!(f, "Recommendation: {}", self.recommendation)?;

        Ok(())
    }
}

/// Find optimal prevalence threshold via spike-in validation.
///
/// Sweeps across threshold values, applies prevalence filtering at each,
/// spikes known effects, runs the analysis pipeline, and measures detection
/// performance to identify the optimal threshold.
///
/// # Arguments
/// * `counts` - The count matrix
/// * `metadata` - Sample metadata
/// * `config` - Optimization configuration
/// * `run_pipeline` - Function that runs the analysis pipeline
///
/// # Returns
/// Optimization result with performance at each threshold and recommendation.
pub fn optimize_prevalence_threshold<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    config: &PrevalenceOptConfig,
    run_pipeline: F,
) -> Result<PrevalenceOptResult>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet> + Sync,
{
    let start = std::time::Instant::now();
    let n_original_features = counts.n_features();

    // Validate thresholds
    for &t in &config.thresholds {
        if !(0.0..=1.0).contains(&t) {
            return Err(DaaError::InvalidParameter(format!(
                "Threshold {} must be between 0 and 1",
                t
            )));
        }
    }

    // Evaluate each threshold
    let results: Vec<ThresholdResult> = config
        .thresholds
        .iter()
        .map(|&threshold| {
            evaluate_threshold(
                counts,
                metadata,
                threshold,
                n_original_features,
                config,
                &run_pipeline,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    // Select optimal threshold
    let (optimal_index, recommendation) = select_optimal(&results, &config.selection_criterion);
    let optimal_threshold = results[optimal_index].threshold;

    Ok(PrevalenceOptResult {
        config: config.clone(),
        results,
        optimal_threshold,
        optimal_index,
        recommendation,
        runtime_seconds: start.elapsed().as_secs_f64(),
    })
}

/// Evaluate a single threshold.
fn evaluate_threshold<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    threshold: f64,
    n_original_features: usize,
    config: &PrevalenceOptConfig,
    run_pipeline: &F,
) -> Result<ThresholdResult>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet> + Sync,
{
    // Apply prevalence filter
    let filtered = match &config.filter_logic {
        PrevalenceFilterLogic::Overall => filter_prevalence_overall(counts, threshold),
        PrevalenceFilterLogic::AnyGroup => filter_prevalence_groupwise(
            counts,
            metadata,
            &config.group_column,
            threshold,
            GroupwiseLogic::Any,
        ),
        PrevalenceFilterLogic::AllGroups => filter_prevalence_groupwise(
            counts,
            metadata,
            &config.group_column,
            threshold,
            GroupwiseLogic::All,
        ),
    };

    // Handle case where no features pass
    let filtered = match filtered {
        Ok(f) => f,
        Err(DaaError::EmptyData(_)) => {
            return Ok(ThresholdResult {
                threshold,
                n_features_retained: 0,
                pct_features_retained: 0.0,
                mean_sensitivity: 0.0,
                std_sensitivity: 0.0,
                mean_fdr: 0.0,
                std_fdr: 0.0,
                mean_f1: 0.0,
                std_f1: 0.0,
                replicates: vec![],
            });
        }
        Err(e) => return Err(e),
    };

    let n_features_retained = filtered.n_features();
    let pct_features_retained = n_features_retained as f64 / n_original_features as f64;

    // Adjust n_spike if needed
    let n_spike = config.n_spike.min(n_features_retained / 2).max(1);

    // Run replicates in parallel
    let replicates: Vec<ThresholdReplicate> = (0..config.n_replicates)
        .into_par_iter()
        .map(|rep| {
            let mut rng = SimpleRng::new(config.seed.wrapping_add(rep as u64 * 1000));
            let seed = rng.next_u64();

            evaluate_replicate(
                &filtered,
                metadata,
                n_spike,
                config,
                rep,
                seed,
                run_pipeline,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    // Aggregate results
    let sensitivities: Vec<f64> = replicates.iter().map(|r| r.sensitivity).collect();
    let fdrs: Vec<f64> = replicates.iter().map(|r| r.fdr).collect();
    let f1s: Vec<f64> = replicates.iter().map(|r| r.f1_score).collect();

    Ok(ThresholdResult {
        threshold,
        n_features_retained,
        pct_features_retained,
        mean_sensitivity: mean(&sensitivities),
        std_sensitivity: std_dev(&sensitivities),
        mean_fdr: mean(&fdrs),
        std_fdr: std_dev(&fdrs),
        mean_f1: mean(&f1s),
        std_f1: std_dev(&f1s),
        replicates,
    })
}

/// Evaluate a single replicate at a threshold.
fn evaluate_replicate<F>(
    filtered_counts: &CountMatrix,
    metadata: &Metadata,
    n_spike: usize,
    config: &PrevalenceOptConfig,
    replicate: usize,
    seed: u64,
    run_pipeline: &F,
) -> Result<ThresholdReplicate>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet>,
{
    // Spike abundance
    let spiked = spike_abundance_with_mode(
        filtered_counts,
        metadata,
        &config.group_column,
        n_spike,
        config.fold_change,
        &config.target_group,
        SpikeSelection::Random,
        seed,
        config.spike_mode,
    )?;

    // Run pipeline on spiked data
    let results = run_pipeline(&spiked.counts, metadata)?;

    // Evaluate
    let eval = evaluate_spikes(&results, &spiked.spec, config.fdr_threshold);

    // Calculate F1
    let precision = if eval.true_positives + eval.false_positives > 0 {
        eval.true_positives as f64 / (eval.true_positives + eval.false_positives) as f64
    } else {
        0.0
    };
    let f1_score = if precision + eval.sensitivity > 0.0 {
        2.0 * precision * eval.sensitivity / (precision + eval.sensitivity)
    } else {
        0.0
    };

    Ok(ThresholdReplicate {
        replicate,
        seed,
        sensitivity: eval.sensitivity,
        fdr: eval.fdr,
        f1_score,
        true_positives: eval.true_positives,
        false_positives: eval.false_positives,
        false_negatives: eval.false_negatives,
        n_spiked: spiked.spec.spiked_features.len(),
    })
}

/// Select optimal threshold based on criterion.
fn select_optimal(results: &[ThresholdResult], criterion: &OptimizationCriterion) -> (usize, String) {
    if results.is_empty() {
        return (0, "No valid thresholds evaluated".to_string());
    }

    let (best_idx, recommendation) = match criterion {
        OptimizationCriterion::MaxF1 => {
            let best = results
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.mean_f1.partial_cmp(&b.mean_f1).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);
            let r = &results[best];
            let rec = format!(
                "Use threshold {:.2} - Maximizes F1 score ({:.2}), retains {} features ({:.0}%), \
                 sensitivity {:.1}% at FDR {:.1}%",
                r.threshold,
                r.mean_f1,
                r.n_features_retained,
                r.pct_features_retained * 100.0,
                r.mean_sensitivity * 100.0,
                r.mean_fdr * 100.0
            );
            (best, rec)
        }
        OptimizationCriterion::MaxSensitivityAtFdr(max_fdr) => {
            // Find threshold with highest sensitivity where FDR <= max_fdr
            let valid: Vec<(usize, &ThresholdResult)> = results
                .iter()
                .enumerate()
                .filter(|(_, r)| r.mean_fdr <= *max_fdr)
                .collect();

            if valid.is_empty() {
                let best = results
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.mean_fdr.partial_cmp(&b.mean_fdr).unwrap())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                let r = &results[best];
                let rec = format!(
                    "WARNING: No threshold achieves FDR ≤ {:.1}%. Using {:.2} with lowest FDR ({:.1}%)",
                    max_fdr * 100.0,
                    r.threshold,
                    r.mean_fdr * 100.0
                );
                (best, rec)
            } else {
                let best = valid
                    .iter()
                    .max_by(|(_, a), (_, b)| a.mean_sensitivity.partial_cmp(&b.mean_sensitivity).unwrap())
                    .map(|(i, _)| *i)
                    .unwrap_or(0);
                let r = &results[best];
                let rec = format!(
                    "Use threshold {:.2} - Maximizes sensitivity ({:.1}%) while keeping FDR ≤ {:.1}% (actual: {:.1}%)",
                    r.threshold,
                    r.mean_sensitivity * 100.0,
                    max_fdr * 100.0,
                    r.mean_fdr * 100.0
                );
                (best, rec)
            }
        }
        OptimizationCriterion::MaxEfficiency => {
            let best = results
                .iter()
                .enumerate()
                .filter(|(_, r)| r.n_features_retained > 0)
                .max_by(|(_, a), (_, b)| a.efficiency().partial_cmp(&b.efficiency()).unwrap())
                .map(|(i, _)| i)
                .unwrap_or(0);
            let r = &results[best];
            let rec = format!(
                "Use threshold {:.2} - Maximizes efficiency (sensitivity per feature), \
                 retains {} features with {:.1}% sensitivity",
                r.threshold, r.n_features_retained, r.mean_sensitivity * 100.0
            );
            (best, rec)
        }
        OptimizationCriterion::MinFdrAtSensitivity(min_sens) => {
            // Find threshold with lowest FDR where sensitivity >= min_sens
            let valid: Vec<(usize, &ThresholdResult)> = results
                .iter()
                .enumerate()
                .filter(|(_, r)| r.mean_sensitivity >= *min_sens)
                .collect();

            if valid.is_empty() {
                let best = results
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| {
                        a.mean_sensitivity.partial_cmp(&b.mean_sensitivity).unwrap()
                    })
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                let r = &results[best];
                let rec = format!(
                    "WARNING: No threshold achieves sensitivity ≥ {:.1}%. Using {:.2} with highest sensitivity ({:.1}%)",
                    min_sens * 100.0,
                    r.threshold,
                    r.mean_sensitivity * 100.0
                );
                (best, rec)
            } else {
                let best = valid
                    .iter()
                    .min_by(|(_, a), (_, b)| a.mean_fdr.partial_cmp(&b.mean_fdr).unwrap())
                    .map(|(i, _)| *i)
                    .unwrap_or(0);
                let r = &results[best];
                let rec = format!(
                    "Use threshold {:.2} - Minimizes FDR ({:.1}%) while maintaining sensitivity ≥ {:.1}% (actual: {:.1}%)",
                    r.threshold,
                    r.mean_fdr * 100.0,
                    min_sens * 100.0,
                    r.mean_sensitivity * 100.0
                );
                (best, rec)
            }
        }
    };

    (best_idx, recommendation)
}

/// Calculate mean of a slice.
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Calculate standard deviation of a slice.
fn std_dev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let m = mean(values);
    let variance = values.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    variance.sqrt()
}

// ============================================================================
// Group-specific optimization
// ============================================================================

/// Configuration for per-group threshold optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerGroupOptConfig {
    /// Base optimization config.
    pub base: PrevalenceOptConfig,
    /// Groups to optimize separately.
    pub groups: Vec<String>,
}

/// Result from per-group optimization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerGroupOptResult {
    /// Optimal threshold for each group.
    pub group_thresholds: HashMap<String, f64>,
    /// Full results for each group.
    pub group_results: HashMap<String, PrevalenceOptResult>,
    /// Overall recommendation.
    pub recommendation: String,
}

impl fmt::Display for PerGroupOptResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Per-Group Prevalence Optimization")?;
        writeln!(f, "==================================")?;
        writeln!(f)?;

        for (group, threshold) in &self.group_thresholds {
            if let Some(result) = self.group_results.get(group) {
                let opt = result.optimal_result();
                writeln!(
                    f,
                    "{}: threshold {:.2} (F1: {:.2}, sensitivity: {:.1}%, FDR: {:.1}%)",
                    group,
                    threshold,
                    opt.mean_f1,
                    opt.mean_sensitivity * 100.0,
                    opt.mean_fdr * 100.0
                )?;
            }
        }

        writeln!(f)?;
        writeln!(f, "{}", self.recommendation)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{CountMatrix, DaResult};
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        // 100 features × 20 samples with varying sparsity
        let mut tri_mat = TriMat::new((100, 20));

        // First 30 features: high prevalence (present in 80% of samples = 16 samples)
        for feat in 0..30 {
            for sample in 0..16 {
                tri_mat.add_triplet(feat, sample, 100);
            }
        }

        // Next 40 features: moderate prevalence (present in 40% of samples = 8 samples)
        for feat in 30..70 {
            for sample in 0..8 {
                tri_mat.add_triplet(feat, sample, 50);
            }
        }

        // Last 30 features: low prevalence (present in 10% of samples = 2 samples)
        for feat in 70..100 {
            for sample in 0..2 {
                tri_mat.add_triplet(feat, sample, 20);
            }
        }

        let feature_ids: Vec<String> = (0..100).map(|i| format!("feature_{}", i)).collect();
        let sample_ids: Vec<String> = (0..20).map(|i| format!("sample_{}", i)).collect();

        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        for i in 0..10 {
            writeln!(file, "sample_{}\tcontrol", i).unwrap();
        }
        for i in 10..20 {
            writeln!(file, "sample_{}\ttreatment", i).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_prevalence_opt_config_default() {
        let config = PrevalenceOptConfig::default();
        assert_eq!(config.thresholds.len(), 8);
        assert_eq!(config.n_replicates, 5);
        assert_eq!(config.fold_change, 2.0);
    }

    #[test]
    fn test_prevalence_opt_config_quick() {
        let config = PrevalenceOptConfig::quick();
        assert_eq!(config.thresholds.len(), 3);
        assert_eq!(config.n_replicates, 2);
    }

    #[test]
    fn test_mean_std() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert!((mean(&values) - 3.0).abs() < 1e-10);
        assert!((std_dev(&values) - 1.5811388300841898).abs() < 1e-10);
    }

    #[test]
    fn test_threshold_result_efficiency() {
        let result = ThresholdResult {
            threshold: 0.1,
            n_features_retained: 100,
            pct_features_retained: 0.5,
            mean_sensitivity: 0.8,
            std_sensitivity: 0.05,
            mean_fdr: 0.05,
            std_fdr: 0.01,
            mean_f1: 0.85,
            std_f1: 0.03,
            replicates: vec![],
        };
        assert!((result.efficiency() - 0.008).abs() < 1e-10);
    }

    #[test]
    fn test_select_optimal_max_f1() {
        let results = vec![
            ThresholdResult {
                threshold: 0.05,
                n_features_retained: 100,
                pct_features_retained: 1.0,
                mean_sensitivity: 0.7,
                std_sensitivity: 0.05,
                mean_fdr: 0.1,
                std_fdr: 0.02,
                mean_f1: 0.75,
                std_f1: 0.03,
                replicates: vec![],
            },
            ThresholdResult {
                threshold: 0.10,
                n_features_retained: 80,
                pct_features_retained: 0.8,
                mean_sensitivity: 0.8,
                std_sensitivity: 0.04,
                mean_fdr: 0.05,
                std_fdr: 0.01,
                mean_f1: 0.85,
                std_f1: 0.02,
                replicates: vec![],
            },
            ThresholdResult {
                threshold: 0.20,
                n_features_retained: 50,
                pct_features_retained: 0.5,
                mean_sensitivity: 0.6,
                std_sensitivity: 0.06,
                mean_fdr: 0.03,
                std_fdr: 0.01,
                mean_f1: 0.72,
                std_f1: 0.04,
                replicates: vec![],
            },
        ];

        let (idx, _) = select_optimal(&results, &OptimizationCriterion::MaxF1);
        assert_eq!(idx, 1); // 0.10 threshold has highest F1
    }

    #[test]
    fn test_select_optimal_max_sensitivity_at_fdr() {
        let results = vec![
            ThresholdResult {
                threshold: 0.05,
                n_features_retained: 100,
                pct_features_retained: 1.0,
                mean_sensitivity: 0.9,
                std_sensitivity: 0.03,
                mean_fdr: 0.15, // Too high
                std_fdr: 0.02,
                mean_f1: 0.80,
                std_f1: 0.03,
                replicates: vec![],
            },
            ThresholdResult {
                threshold: 0.10,
                n_features_retained: 80,
                pct_features_retained: 0.8,
                mean_sensitivity: 0.75,
                std_sensitivity: 0.04,
                mean_fdr: 0.08, // Still high
                std_fdr: 0.02,
                mean_f1: 0.78,
                std_f1: 0.03,
                replicates: vec![],
            },
            ThresholdResult {
                threshold: 0.20,
                n_features_retained: 50,
                pct_features_retained: 0.5,
                mean_sensitivity: 0.65,
                std_sensitivity: 0.05,
                mean_fdr: 0.04, // Good
                std_fdr: 0.01,
                mean_f1: 0.75,
                std_f1: 0.04,
                replicates: vec![],
            },
        ];

        let (idx, _) = select_optimal(&results, &OptimizationCriterion::MaxSensitivityAtFdr(0.05));
        assert_eq!(idx, 2); // 0.20 is only one with FDR <= 5%
    }

    #[test]
    fn test_prevalence_opt_result_curves() {
        let results = vec![
            ThresholdResult {
                threshold: 0.05,
                n_features_retained: 100,
                pct_features_retained: 1.0,
                mean_sensitivity: 0.7,
                std_sensitivity: 0.05,
                mean_fdr: 0.1,
                std_fdr: 0.02,
                mean_f1: 0.75,
                std_f1: 0.03,
                replicates: vec![],
            },
            ThresholdResult {
                threshold: 0.10,
                n_features_retained: 80,
                pct_features_retained: 0.8,
                mean_sensitivity: 0.8,
                std_sensitivity: 0.04,
                mean_fdr: 0.05,
                std_fdr: 0.01,
                mean_f1: 0.85,
                std_f1: 0.02,
                replicates: vec![],
            },
        ];

        let opt_result = PrevalenceOptResult {
            config: PrevalenceOptConfig::default(),
            results,
            optimal_threshold: 0.10,
            optimal_index: 1,
            recommendation: "Use 0.10".to_string(),
            runtime_seconds: 1.0,
        };

        let sens_curve = opt_result.sensitivity_curve();
        assert_eq!(sens_curve.len(), 2);
        assert!((sens_curve[0].1 - 0.7).abs() < 1e-10);
        assert!((sens_curve[1].1 - 0.8).abs() < 1e-10);

        let feat_curve = opt_result.features_curve();
        assert_eq!(feat_curve[0].1, 100);
        assert_eq!(feat_curve[1].1, 80);
    }

    #[test]
    fn test_prevalence_opt_result_csv() {
        let results = vec![ThresholdResult {
            threshold: 0.10,
            n_features_retained: 80,
            pct_features_retained: 0.8,
            mean_sensitivity: 0.8,
            std_sensitivity: 0.04,
            mean_fdr: 0.05,
            std_fdr: 0.01,
            mean_f1: 0.85,
            std_f1: 0.02,
            replicates: vec![],
        }];

        let opt_result = PrevalenceOptResult {
            config: PrevalenceOptConfig::default(),
            results,
            optimal_threshold: 0.10,
            optimal_index: 0,
            recommendation: "Use 0.10".to_string(),
            runtime_seconds: 1.0,
        };

        let csv = opt_result.to_csv();
        assert!(csv.contains("threshold,n_features"));
        assert!(csv.contains("0.100,80,80.0"));
    }

    #[test]
    fn test_optimize_prevalence_threshold() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let config = PrevalenceOptConfig {
            thresholds: vec![0.05, 0.20, 0.50],
            fold_change: 3.0,
            n_spike: 5,
            n_replicates: 2,
            group_column: "group".into(),
            target_group: "treatment".into(),
            fdr_threshold: 0.10,
            selection_criterion: OptimizationCriterion::MaxF1,
            filter_logic: PrevalenceFilterLogic::Overall,
            spike_mode: SpikeMode::Raw,
            seed: 42,
        };

        // Simple pipeline: return mock results
        let pipeline = |counts: &CountMatrix, _metadata: &Metadata| -> Result<DaResultSet> {
            // Return results for all features with random p-values
            let results: Vec<DaResult> = counts
                .feature_ids()
                .iter()
                .enumerate()
                .map(|(i, fid)| {
                    // Make first few features significant
                    let q_value = if i < 5 { 0.01 } else { 0.5 };
                    DaResult::new(
                        fid.clone(),
                        "group".into(),
                        1.0,      // estimate
                        0.5,      // std_error
                        2.0,      // statistic
                        q_value * 0.1, // p_value
                        q_value,  // q_value
                        0.8,      // prevalence
                        100.0,    // mean_abundance
                    )
                })
                .collect();

            Ok(DaResultSet::new("mock".into(), results))
        };

        let result = optimize_prevalence_threshold(&counts, &metadata, &config, pipeline).unwrap();

        // Should have results for each threshold
        assert_eq!(result.results.len(), 3);

        // Features should decrease as threshold increases
        assert!(result.results[0].n_features_retained >= result.results[2].n_features_retained);

        // Should have a recommendation
        assert!(!result.recommendation.is_empty());
    }
}
