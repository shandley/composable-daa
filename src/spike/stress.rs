//! Compositional stress testing framework for evaluating differential abundance pipelines.
//!
//! This module provides systematic benchmarking to characterize how compositional effects
//! impact differential abundance detection across varying conditions.

use crate::data::{CountMatrix, DaResultSet, Metadata};
use crate::error::Result;
use crate::spike::evaluate::evaluate_spikes;
use crate::spike::types::{SpikeMode, SpikeSelection};
use crate::spike::abundance::spike_abundance_with_mode;
use crate::spike::validate::shuffle_labels;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Instant;

/// Simple LCG random number generator for reproducibility.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }
}

// ============================================================================
// Core Types
// ============================================================================

/// Configuration for a compositional stress test.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressConfig {
    /// Name of this stress test.
    pub name: String,
    /// Fractions of features to spike (e.g., [0.01, 0.05, 0.10, 0.25]).
    pub spike_fractions: Vec<f64>,
    /// Fold changes to test (e.g., [1.5, 2.0, 3.0, 5.0]).
    pub fold_changes: Vec<f64>,
    /// Spike modes to compare.
    pub spike_modes: Vec<SpikeMode>,
    /// Metadata column containing group labels.
    pub group_column: String,
    /// Group to receive the spike effect.
    pub target_group: String,
    /// FDR threshold for calling significance.
    pub fdr_threshold: f64,
    /// How to select features for spiking.
    pub selection: SpikeSelection,
    /// Number of replicates per parameter combination.
    pub n_replicates: usize,
    /// Number of permutations for FDR calibration.
    pub n_permutations: usize,
    /// Base random seed.
    pub seed: u64,
    /// Target power for power curve analysis (e.g., 0.80).
    pub power_target: f64,
}

impl Default for StressConfig {
    fn default() -> Self {
        Self {
            name: "stress_test".into(),
            spike_fractions: vec![0.01, 0.05, 0.10, 0.25],
            fold_changes: vec![1.5, 2.0, 3.0, 5.0],
            spike_modes: vec![SpikeMode::Raw, SpikeMode::Compositional, SpikeMode::Absolute],
            group_column: "group".into(),
            target_group: "treatment".into(),
            fdr_threshold: 0.05,
            selection: SpikeSelection::Random,
            n_replicates: 5,
            n_permutations: 5,
            seed: 42,
            power_target: 0.80,
        }
    }
}

impl StressConfig {
    /// Create a new stress config with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Create a quick config for fast testing (fewer replicates and permutations).
    pub fn quick() -> Self {
        Self {
            name: "quick_stress".into(),
            spike_fractions: vec![0.05, 0.10],
            fold_changes: vec![2.0, 3.0],
            spike_modes: vec![SpikeMode::Compositional],
            n_replicates: 2,
            n_permutations: 2,
            ..Default::default()
        }
    }

    /// Set spike fractions.
    pub fn with_spike_fractions(mut self, fractions: Vec<f64>) -> Self {
        self.spike_fractions = fractions;
        self
    }

    /// Set fold changes.
    pub fn with_fold_changes(mut self, fold_changes: Vec<f64>) -> Self {
        self.fold_changes = fold_changes;
        self
    }

    /// Set spike modes.
    pub fn with_modes(mut self, modes: Vec<SpikeMode>) -> Self {
        self.spike_modes = modes;
        self
    }

    /// Set group column and target.
    pub fn with_groups(mut self, column: &str, target: &str) -> Self {
        self.group_column = column.into();
        self.target_group = target.into();
        self
    }

    /// Set replication counts.
    pub fn with_replicates(mut self, n_replicates: usize, n_permutations: usize) -> Self {
        self.n_replicates = n_replicates;
        self.n_permutations = n_permutations;
        self
    }
}

/// Parameters for a single stress test run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressParams {
    /// Fraction of features to spike.
    pub spike_fraction: f64,
    /// Fold change to apply.
    pub fold_change: f64,
    /// Spike mode.
    pub mode: SpikeMode,
    /// Number of features being spiked.
    pub n_spiked: usize,
}

/// Result from a single stress test run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressRunResult {
    /// Parameters used for this run.
    pub params: StressParams,
    /// Replicate number.
    pub replicate: usize,
    /// Random seed used.
    pub seed: u64,

    // Detection metrics
    /// Sensitivity (true positive rate).
    pub sensitivity: f64,
    /// Observed false discovery rate.
    pub fdr: f64,
    /// F1 score.
    pub f1_score: f64,

    // Compositional diagnostics
    /// Ratio of geometric means (spiked / original).
    pub geometric_mean_ratio: f64,
    /// Effective CLR effect size.
    pub effective_clr_effect: f64,
    /// Effect attenuation: effective / nominal (1.0 = no attenuation).
    pub effect_attenuation: f64,

    // Permutation FPR
    /// Mean false positive rate from permutation tests.
    pub permuted_fpr_mean: f64,

    /// Runtime in seconds.
    pub runtime_seconds: f64,
}

/// Aggregated metrics across replicates for a parameter combination.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AggregatedMetrics {
    /// Parameters for this aggregate.
    pub params: StressParams,
    /// Number of replicates.
    pub n_replicates: usize,

    // Sensitivity
    /// Mean sensitivity.
    pub sensitivity_mean: f64,
    /// Standard deviation of sensitivity.
    pub sensitivity_std: f64,

    // FDR
    /// Mean FDR.
    pub fdr_mean: f64,
    /// Standard deviation of FDR.
    pub fdr_std: f64,

    // F1
    /// Mean F1 score.
    pub f1_mean: f64,
    /// Standard deviation of F1.
    pub f1_std: f64,

    // Effect attenuation
    /// Mean effect attenuation.
    pub attenuation_mean: f64,
    /// Standard deviation of attenuation.
    pub attenuation_std: f64,

    // Permuted FPR
    /// Mean of permuted FPR means.
    pub permuted_fpr_mean: f64,
}

/// A point on a power curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowerPoint {
    /// Spike fraction.
    pub spike_fraction: f64,
    /// Number of features spiked.
    pub n_spiked: usize,
    /// Fold change needed to achieve target power.
    pub fold_change_for_power: f64,
    /// Achieved sensitivity at this point.
    pub achieved_sensitivity: f64,
}

/// Comparison summary across spike modes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModeComparisonSummary {
    /// Best mode for balanced sensitivity/FDR.
    pub recommended_mode: SpikeMode,
    /// Mode with highest sensitivity.
    pub highest_sensitivity_mode: SpikeMode,
    /// Mode with lowest FDR.
    pub lowest_fdr_mode: SpikeMode,
    /// Per-mode statistics.
    pub mode_stats: HashMap<String, ModeStats>,
}

/// Statistics for a single mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModeStats {
    /// Mean sensitivity across all conditions.
    pub mean_sensitivity: f64,
    /// Mean FDR across all conditions.
    pub mean_fdr: f64,
    /// Mean effect attenuation.
    pub mean_attenuation: f64,
}

/// Complete results from a stress test.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StressSummary {
    /// Configuration used.
    pub config: StressConfig,
    /// All individual run results.
    pub runs: Vec<StressRunResult>,
    /// Aggregated metrics by parameter combination.
    pub aggregated: Vec<AggregatedMetrics>,
    /// Power curves by spike mode.
    pub power_curves: HashMap<String, Vec<PowerPoint>>,
    /// Mode comparison summary.
    pub mode_comparison: ModeComparisonSummary,
    /// Total runtime in seconds.
    pub total_runtime_seconds: f64,
}

// ============================================================================
// Grid Generation
// ============================================================================

/// Generate the parameter grid for stress testing.
pub fn generate_parameter_grid(config: &StressConfig, n_features: usize) -> Vec<StressParams> {
    let mut grid = Vec::new();

    for &fraction in &config.spike_fractions {
        let n_spiked = ((n_features as f64) * fraction).max(1.0) as usize;
        for &fold_change in &config.fold_changes {
            for &mode in &config.spike_modes {
                grid.push(StressParams {
                    spike_fraction: fraction,
                    fold_change,
                    mode,
                    n_spiked,
                });
            }
        }
    }

    grid
}

// ============================================================================
// Single Run Execution
// ============================================================================

/// Execute a single stress test run.
pub fn execute_stress_run<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    params: &StressParams,
    config: &StressConfig,
    replicate: usize,
    run_pipeline: &F,
) -> Result<StressRunResult>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet> + Sync,
{
    let start = Instant::now();

    // Generate seed for this run
    let mut rng = SimpleRng::new(config.seed);
    for _ in 0..(replicate * 1000) {
        rng.next_u64();
    }
    let run_seed = rng.next_u64();

    // Spike the data
    let spiked = spike_abundance_with_mode(
        counts,
        metadata,
        &config.group_column,
        params.n_spiked,
        params.fold_change,
        &config.target_group,
        config.selection.clone(),
        run_seed,
        params.mode,
    )?;

    // Run pipeline on spiked data
    let results = run_pipeline(&spiked.counts, metadata)?;

    // Evaluate spike detection
    let eval = evaluate_spikes(&results, &spiked.spec, config.fdr_threshold);

    // Extract compositional diagnostics
    let (geometric_mean_ratio, effective_clr_effect) = if let Some(diag) = &spiked.spec.diagnostics {
        (diag.geometric_mean_ratio, diag.effective_clr_effect)
    } else {
        (1.0, params.fold_change.ln())
    };

    // Calculate effect attenuation
    let nominal_effect = params.fold_change.ln();
    let effect_attenuation = if nominal_effect.abs() > 1e-10 {
        effective_clr_effect / nominal_effect
    } else {
        1.0
    };

    // Run permutation tests for FPR calibration
    let mut permuted_fprs = Vec::with_capacity(config.n_permutations);
    for perm_idx in 0..config.n_permutations {
        let perm_seed = run_seed.wrapping_add(10000 + perm_idx as u64);
        let shuffled = shuffle_labels(metadata, &config.group_column, perm_seed)?;
        let perm_results = run_pipeline(counts, &shuffled)?;

        let n_sig = perm_results
            .results
            .iter()
            .filter(|r| r.q_value < config.fdr_threshold)
            .count();
        let fpr = if perm_results.len() > 0 {
            n_sig as f64 / perm_results.len() as f64
        } else {
            0.0
        };
        permuted_fprs.push(fpr);
    }

    let permuted_fpr_mean = if !permuted_fprs.is_empty() {
        permuted_fprs.iter().sum::<f64>() / permuted_fprs.len() as f64
    } else {
        0.0
    };

    let runtime_seconds = start.elapsed().as_secs_f64();

    Ok(StressRunResult {
        params: params.clone(),
        replicate,
        seed: run_seed,
        sensitivity: eval.sensitivity,
        fdr: eval.fdr,
        f1_score: eval.f1_score,
        geometric_mean_ratio,
        effective_clr_effect,
        effect_attenuation,
        permuted_fpr_mean,
        runtime_seconds,
    })
}

// ============================================================================
// Aggregation
// ============================================================================

/// Aggregate results across replicates.
pub fn aggregate_results(runs: &[StressRunResult]) -> Vec<AggregatedMetrics> {
    // Group by parameter combination
    let mut groups: HashMap<String, Vec<&StressRunResult>> = HashMap::new();

    for run in runs {
        let key = format!(
            "{:.3}_{:.1}_{:?}",
            run.params.spike_fraction, run.params.fold_change, run.params.mode
        );
        groups.entry(key).or_default().push(run);
    }

    let mut aggregated = Vec::new();

    for (_key, group_runs) in groups {
        if group_runs.is_empty() {
            continue;
        }

        let n = group_runs.len() as f64;
        let params = group_runs[0].params.clone();

        // Compute means
        let sensitivity_mean = group_runs.iter().map(|r| r.sensitivity).sum::<f64>() / n;
        let fdr_mean = group_runs.iter().map(|r| r.fdr).sum::<f64>() / n;
        let f1_mean = group_runs.iter().map(|r| r.f1_score).sum::<f64>() / n;
        let attenuation_mean = group_runs.iter().map(|r| r.effect_attenuation).sum::<f64>() / n;
        let permuted_fpr_mean = group_runs.iter().map(|r| r.permuted_fpr_mean).sum::<f64>() / n;

        // Compute standard deviations
        let sensitivity_std = std_dev(group_runs.iter().map(|r| r.sensitivity), sensitivity_mean);
        let fdr_std = std_dev(group_runs.iter().map(|r| r.fdr), fdr_mean);
        let f1_std = std_dev(group_runs.iter().map(|r| r.f1_score), f1_mean);
        let attenuation_std = std_dev(group_runs.iter().map(|r| r.effect_attenuation), attenuation_mean);

        aggregated.push(AggregatedMetrics {
            params,
            n_replicates: group_runs.len(),
            sensitivity_mean,
            sensitivity_std,
            fdr_mean,
            fdr_std,
            f1_mean,
            f1_std,
            attenuation_mean,
            attenuation_std,
            permuted_fpr_mean,
        });
    }

    // Sort by spike fraction, then fold change, then mode
    aggregated.sort_by(|a, b| {
        a.params
            .spike_fraction
            .partial_cmp(&b.params.spike_fraction)
            .unwrap()
            .then_with(|| {
                a.params
                    .fold_change
                    .partial_cmp(&b.params.fold_change)
                    .unwrap()
            })
            .then_with(|| format!("{:?}", a.params.mode).cmp(&format!("{:?}", b.params.mode)))
    });

    aggregated
}

fn std_dev(values: impl Iterator<Item = f64>, mean: f64) -> f64 {
    let values: Vec<f64> = values.collect();
    if values.len() < 2 {
        return 0.0;
    }
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64;
    variance.sqrt()
}

// ============================================================================
// Power Curve Calculation
// ============================================================================

/// Calculate power curves for each spike mode.
pub fn calculate_power_curves(
    aggregated: &[AggregatedMetrics],
    power_target: f64,
) -> HashMap<String, Vec<PowerPoint>> {
    let mut curves: HashMap<String, Vec<PowerPoint>> = HashMap::new();

    // Group by mode and spike fraction
    for mode in [SpikeMode::Raw, SpikeMode::Compositional, SpikeMode::Absolute] {
        let mode_key = format!("{:?}", mode);
        let mut points = Vec::new();

        // Get unique spike fractions
        let mut fractions: Vec<f64> = aggregated
            .iter()
            .filter(|a| a.params.mode == mode)
            .map(|a| a.params.spike_fraction)
            .collect();
        fractions.sort_by(|a, b| a.partial_cmp(b).unwrap());
        fractions.dedup();

        for fraction in fractions {
            // Get all fold changes for this fraction/mode, sorted by fold change
            let mut candidates: Vec<_> = aggregated
                .iter()
                .filter(|a| a.params.mode == mode && (a.params.spike_fraction - fraction).abs() < 1e-6)
                .collect();
            candidates.sort_by(|a, b| a.params.fold_change.partial_cmp(&b.params.fold_change).unwrap());

            // Find the minimum fold change that achieves target power
            let mut fc_for_power = f64::INFINITY;
            let mut achieved = 0.0;

            for agg in &candidates {
                if agg.sensitivity_mean >= power_target {
                    fc_for_power = agg.params.fold_change;
                    achieved = agg.sensitivity_mean;
                    break;
                }
            }

            // If we didn't find one, use the max fold change and its sensitivity
            if fc_for_power.is_infinite() && !candidates.is_empty() {
                let last = candidates.last().unwrap();
                fc_for_power = last.params.fold_change;
                achieved = last.sensitivity_mean;
            }

            if !candidates.is_empty() {
                points.push(PowerPoint {
                    spike_fraction: fraction,
                    n_spiked: candidates[0].params.n_spiked,
                    fold_change_for_power: fc_for_power,
                    achieved_sensitivity: achieved,
                });
            }
        }

        if !points.is_empty() {
            curves.insert(mode_key, points);
        }
    }

    curves
}

// ============================================================================
// Mode Comparison
// ============================================================================

/// Compare performance across spike modes.
pub fn compare_modes(aggregated: &[AggregatedMetrics]) -> ModeComparisonSummary {
    let mut mode_stats: HashMap<String, ModeStats> = HashMap::new();

    for mode in [SpikeMode::Raw, SpikeMode::Compositional, SpikeMode::Absolute] {
        let mode_key = format!("{:?}", mode);
        let mode_results: Vec<_> = aggregated.iter().filter(|a| a.params.mode == mode).collect();

        if mode_results.is_empty() {
            continue;
        }

        let n = mode_results.len() as f64;
        let mean_sensitivity = mode_results.iter().map(|a| a.sensitivity_mean).sum::<f64>() / n;
        let mean_fdr = mode_results.iter().map(|a| a.fdr_mean).sum::<f64>() / n;
        let mean_attenuation = mode_results.iter().map(|a| a.attenuation_mean).sum::<f64>() / n;

        mode_stats.insert(
            mode_key,
            ModeStats {
                mean_sensitivity,
                mean_fdr,
                mean_attenuation,
            },
        );
    }

    // Find best modes
    let highest_sensitivity_mode = mode_stats
        .iter()
        .max_by(|a, b| a.1.mean_sensitivity.partial_cmp(&b.1.mean_sensitivity).unwrap())
        .map(|(k, _)| parse_mode(k))
        .unwrap_or(SpikeMode::Compositional);

    let lowest_fdr_mode = mode_stats
        .iter()
        .min_by(|a, b| a.1.mean_fdr.partial_cmp(&b.1.mean_fdr).unwrap())
        .map(|(k, _)| parse_mode(k))
        .unwrap_or(SpikeMode::Compositional);

    // Recommended mode: balance sensitivity and FDR using F1-like score
    let recommended_mode = mode_stats
        .iter()
        .max_by(|a, b| {
            let score_a = 2.0 * a.1.mean_sensitivity * (1.0 - a.1.mean_fdr)
                / (a.1.mean_sensitivity + (1.0 - a.1.mean_fdr) + 1e-10);
            let score_b = 2.0 * b.1.mean_sensitivity * (1.0 - b.1.mean_fdr)
                / (b.1.mean_sensitivity + (1.0 - b.1.mean_fdr) + 1e-10);
            score_a.partial_cmp(&score_b).unwrap()
        })
        .map(|(k, _)| parse_mode(k))
        .unwrap_or(SpikeMode::Compositional);

    ModeComparisonSummary {
        recommended_mode,
        highest_sensitivity_mode,
        lowest_fdr_mode,
        mode_stats,
    }
}

fn parse_mode(s: &str) -> SpikeMode {
    match s {
        "Raw" => SpikeMode::Raw,
        "Absolute" => SpikeMode::Absolute,
        _ => SpikeMode::Compositional,
    }
}

// ============================================================================
// Main Runner
// ============================================================================

/// Run a complete stress test.
///
/// This function systematically tests the pipeline across different combinations of:
/// - Spike fractions (how many features are affected)
/// - Fold changes (effect size)
/// - Spike modes (how compositional constraints are handled)
///
/// # Arguments
/// * `counts` - Original count matrix
/// * `metadata` - Sample metadata
/// * `config` - Stress test configuration
/// * `run_pipeline` - Function that takes (counts, metadata) and returns DA results
///
/// # Returns
/// StressSummary with comprehensive results and analysis.
pub fn run_stress_test<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    config: &StressConfig,
    run_pipeline: F,
) -> Result<StressSummary>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet> + Sync,
{
    let total_start = Instant::now();

    // Generate parameter grid
    let grid = generate_parameter_grid(config, counts.n_features());

    // Generate all (params, replicate) combinations
    let work_items: Vec<(StressParams, usize)> = grid
        .iter()
        .flat_map(|params| (0..config.n_replicates).map(move |rep| (params.clone(), rep)))
        .collect();

    // Run all combinations in parallel
    let runs: Vec<StressRunResult> = work_items
        .par_iter()
        .filter_map(|(params, replicate)| {
            execute_stress_run(counts, metadata, params, config, *replicate, &run_pipeline).ok()
        })
        .collect();

    // Aggregate results
    let aggregated = aggregate_results(&runs);

    // Calculate power curves
    let power_curves = calculate_power_curves(&aggregated, config.power_target);

    // Compare modes
    let mode_comparison = compare_modes(&aggregated);

    let total_runtime_seconds = total_start.elapsed().as_secs_f64();

    Ok(StressSummary {
        config: config.clone(),
        runs,
        aggregated,
        power_curves,
        mode_comparison,
        total_runtime_seconds,
    })
}

// ============================================================================
// Display Implementations
// ============================================================================

impl std::fmt::Display for StressSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Compositional Stress Test Results")?;
        writeln!(f, "==================================")?;
        writeln!(f)?;

        // Summary by spike fraction
        writeln!(f, "Summary by Spike Fraction:")?;
        let mut fractions: Vec<f64> = self
            .aggregated
            .iter()
            .map(|a| a.params.spike_fraction)
            .collect();
        fractions.sort_by(|a, b| a.partial_cmp(b).unwrap());
        fractions.dedup();

        for fraction in &fractions {
            let fraction_results: Vec<_> = self
                .aggregated
                .iter()
                .filter(|a| (a.params.spike_fraction - fraction).abs() < 1e-6)
                .collect();

            if fraction_results.is_empty() {
                continue;
            }

            let n_spiked = fraction_results[0].params.n_spiked;
            let mean_sens: f64 = fraction_results.iter().map(|a| a.sensitivity_mean).sum::<f64>()
                / fraction_results.len() as f64;
            let mean_std: f64 = fraction_results.iter().map(|a| a.sensitivity_std).sum::<f64>()
                / fraction_results.len() as f64;
            let mean_atten: f64 = fraction_results.iter().map(|a| a.attenuation_mean).sum::<f64>()
                / fraction_results.len() as f64;

            writeln!(
                f,
                "  {:.0}% spiked ({} features):",
                fraction * 100.0,
                n_spiked
            )?;
            writeln!(
                f,
                "    Sensitivity: {:.1}% +/- {:.1}%",
                mean_sens * 100.0,
                mean_std * 100.0
            )?;

            let attenuation_pct = (1.0 - mean_atten) * 100.0;
            let severity = if attenuation_pct < 10.0 {
                "minimal"
            } else if attenuation_pct < 30.0 {
                "moderate"
            } else {
                "substantial"
            };
            writeln!(
                f,
                "    Effect attenuation: {:.0}% ({})",
                attenuation_pct, severity
            )?;

            if attenuation_pct > 40.0 {
                writeln!(f, "    WARNING: Compositional closure reduces power")?;
            }
            writeln!(f)?;
        }

        // Power analysis
        writeln!(
            f,
            "Power Analysis (FC for {:.0}% sensitivity):",
            self.config.power_target * 100.0
        )?;
        for (mode, points) in &self.power_curves {
            writeln!(f, "  {} mode:", mode)?;
            for point in points {
                if point.fold_change_for_power.is_finite() && point.achieved_sensitivity >= self.config.power_target {
                    writeln!(
                        f,
                        "    {:.0}% spiked: {:.1}x",
                        point.spike_fraction * 100.0,
                        point.fold_change_for_power
                    )?;
                } else {
                    writeln!(
                        f,
                        "    {:.0}% spiked: >{:.1}x (max tested, achieved {:.0}%)",
                        point.spike_fraction * 100.0,
                        point.fold_change_for_power,
                        point.achieved_sensitivity * 100.0
                    )?;
                }
            }
        }
        writeln!(f)?;

        // Mode comparison
        writeln!(f, "Mode Comparison:")?;
        for (mode, stats) in &self.mode_comparison.mode_stats {
            writeln!(
                f,
                "  {}: sens={:.1}%, FDR={:.1}%, atten={:.0}%",
                mode,
                stats.mean_sensitivity * 100.0,
                stats.mean_fdr * 100.0,
                (1.0 - stats.mean_attenuation) * 100.0
            )?;
        }
        writeln!(
            f,
            "  Recommendation: Use '{:?}' for balanced sensitivity/FDR",
            self.mode_comparison.recommended_mode
        )?;
        writeln!(f)?;

        writeln!(
            f,
            "Total runtime: {:.1}s ({} runs)",
            self.total_runtime_seconds,
            self.runs.len()
        )?;

        Ok(())
    }
}

impl StressSummary {
    /// Convert to CSV format for export.
    pub fn to_csv(&self) -> String {
        let mut csv = String::new();

        // Header
        csv.push_str("spike_fraction,n_spiked,fold_change,mode,");
        csv.push_str("sensitivity_mean,sensitivity_std,");
        csv.push_str("fdr_mean,fdr_std,");
        csv.push_str("f1_mean,f1_std,");
        csv.push_str("attenuation_mean,attenuation_std,");
        csv.push_str("permuted_fpr_mean,n_replicates\n");

        // Data rows
        for agg in &self.aggregated {
            csv.push_str(&format!(
                "{:.4},{},{:.2},{:?},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{}\n",
                agg.params.spike_fraction,
                agg.params.n_spiked,
                agg.params.fold_change,
                agg.params.mode,
                agg.sensitivity_mean,
                agg.sensitivity_std,
                agg.fdr_mean,
                agg.fdr_std,
                agg.f1_mean,
                agg.f1_std,
                agg.attenuation_mean,
                agg.attenuation_std,
                agg.permuted_fpr_mean,
                agg.n_replicates,
            ));
        }

        csv
    }

    /// Convert to JSON format for export.
    pub fn to_json(&self) -> Result<String> {
        Ok(serde_json::to_string_pretty(self)?)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::DaResult;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        let mut tri_mat = TriMat::new((100, 20));
        for feat in 0..100 {
            for sample in 0..20 {
                // Varying prevalence and abundance
                if sample < 20 - feat % 5 {
                    tri_mat.add_triplet(feat, sample, 50 + (feat as u64 % 50));
                }
            }
        }
        let feature_ids: Vec<String> = (0..100).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..20).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        for i in 0..10 {
            writeln!(file, "S{}\tcontrol", i).unwrap();
        }
        for i in 10..20 {
            writeln!(file, "S{}\ttreatment", i).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn mock_pipeline(_counts: &CountMatrix, _metadata: &Metadata) -> Result<DaResultSet> {
        // Return mock results
        let mut results = Vec::new();
        for i in 0..100 {
            // First 10 features are "significant"
            let q_value = if i < 10 { 0.01 } else { 0.5 };
            results.push(DaResult::new(
                format!("feat_{}", i),
                "group".into(),
                0.5,
                0.1,
                5.0,
                q_value * 0.1,
                q_value,
                0.8,
                100.0,
            ));
        }
        Ok(DaResultSet::new("mock".into(), results))
    }

    #[test]
    fn test_stress_config_default() {
        let config = StressConfig::default();
        assert_eq!(config.spike_fractions.len(), 4);
        assert_eq!(config.fold_changes.len(), 4);
        assert_eq!(config.spike_modes.len(), 3);
        assert_eq!(config.n_replicates, 5);
    }

    #[test]
    fn test_stress_config_quick() {
        let config = StressConfig::quick();
        assert_eq!(config.n_replicates, 2);
        assert_eq!(config.n_permutations, 2);
        assert!(config.spike_fractions.len() < StressConfig::default().spike_fractions.len());
    }

    #[test]
    fn test_generate_parameter_grid() {
        let config = StressConfig::default()
            .with_spike_fractions(vec![0.05, 0.10])
            .with_fold_changes(vec![2.0, 3.0])
            .with_modes(vec![SpikeMode::Compositional]);

        let grid = generate_parameter_grid(&config, 100);

        // 2 fractions * 2 fold changes * 1 mode = 4 combinations
        assert_eq!(grid.len(), 4);

        // Check n_spiked calculation
        let five_pct: Vec<_> = grid.iter().filter(|p| (p.spike_fraction - 0.05).abs() < 1e-6).collect();
        assert!(five_pct.iter().all(|p| p.n_spiked == 5));

        let ten_pct: Vec<_> = grid.iter().filter(|p| (p.spike_fraction - 0.10).abs() < 1e-6).collect();
        assert!(ten_pct.iter().all(|p| p.n_spiked == 10));
    }

    #[test]
    fn test_aggregate_results() {
        let params = StressParams {
            spike_fraction: 0.05,
            fold_change: 2.0,
            mode: SpikeMode::Compositional,
            n_spiked: 5,
        };

        let runs = vec![
            StressRunResult {
                params: params.clone(),
                replicate: 0,
                seed: 42,
                sensitivity: 0.80,
                fdr: 0.10,
                f1_score: 0.84,
                geometric_mean_ratio: 1.05,
                effective_clr_effect: 0.65,
                effect_attenuation: 0.94,
                permuted_fpr_mean: 0.05,
                runtime_seconds: 1.0,
            },
            StressRunResult {
                params: params.clone(),
                replicate: 1,
                seed: 43,
                sensitivity: 0.70,
                fdr: 0.15,
                f1_score: 0.76,
                geometric_mean_ratio: 1.06,
                effective_clr_effect: 0.63,
                effect_attenuation: 0.91,
                permuted_fpr_mean: 0.06,
                runtime_seconds: 1.1,
            },
        ];

        let aggregated = aggregate_results(&runs);
        assert_eq!(aggregated.len(), 1);

        let agg = &aggregated[0];
        assert_eq!(agg.n_replicates, 2);
        assert!((agg.sensitivity_mean - 0.75).abs() < 0.01);
        assert!((agg.fdr_mean - 0.125).abs() < 0.01);
    }

    #[test]
    fn test_calculate_power_curves() {
        let aggregated = vec![
            AggregatedMetrics {
                params: StressParams {
                    spike_fraction: 0.05,
                    fold_change: 2.0,
                    mode: SpikeMode::Compositional,
                    n_spiked: 5,
                },
                n_replicates: 5,
                sensitivity_mean: 0.60,
                sensitivity_std: 0.05,
                fdr_mean: 0.10,
                fdr_std: 0.02,
                f1_mean: 0.70,
                f1_std: 0.03,
                attenuation_mean: 0.95,
                attenuation_std: 0.02,
                permuted_fpr_mean: 0.05,
            },
            AggregatedMetrics {
                params: StressParams {
                    spike_fraction: 0.05,
                    fold_change: 3.0,
                    mode: SpikeMode::Compositional,
                    n_spiked: 5,
                },
                n_replicates: 5,
                sensitivity_mean: 0.85,
                sensitivity_std: 0.04,
                fdr_mean: 0.08,
                fdr_std: 0.02,
                f1_mean: 0.88,
                f1_std: 0.02,
                attenuation_mean: 0.92,
                attenuation_std: 0.02,
                permuted_fpr_mean: 0.05,
            },
        ];

        let curves = calculate_power_curves(&aggregated, 0.80);

        assert!(curves.contains_key("Compositional"));
        let comp_curve = &curves["Compositional"];
        assert_eq!(comp_curve.len(), 1);

        // FC=2.0 gives 60% sens, FC=3.0 gives 85% sens
        // So FC for 80% power should be 3.0
        assert!((comp_curve[0].fold_change_for_power - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_compare_modes() {
        let aggregated = vec![
            AggregatedMetrics {
                params: StressParams {
                    spike_fraction: 0.05,
                    fold_change: 2.0,
                    mode: SpikeMode::Raw,
                    n_spiked: 5,
                },
                n_replicates: 5,
                sensitivity_mean: 0.90,
                sensitivity_std: 0.03,
                fdr_mean: 0.20,
                fdr_std: 0.05,
                f1_mean: 0.85,
                f1_std: 0.03,
                attenuation_mean: 1.00,
                attenuation_std: 0.01,
                permuted_fpr_mean: 0.05,
            },
            AggregatedMetrics {
                params: StressParams {
                    spike_fraction: 0.05,
                    fold_change: 2.0,
                    mode: SpikeMode::Compositional,
                    n_spiked: 5,
                },
                n_replicates: 5,
                sensitivity_mean: 0.80,
                sensitivity_std: 0.04,
                fdr_mean: 0.08,
                fdr_std: 0.02,
                f1_mean: 0.85,
                f1_std: 0.03,
                attenuation_mean: 0.95,
                attenuation_std: 0.02,
                permuted_fpr_mean: 0.05,
            },
        ];

        let comparison = compare_modes(&aggregated);

        // Raw has higher sensitivity
        assert_eq!(comparison.highest_sensitivity_mode, SpikeMode::Raw);
        // Compositional has lower FDR
        assert_eq!(comparison.lowest_fdr_mode, SpikeMode::Compositional);
    }

    #[test]
    fn test_execute_stress_run() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let config = StressConfig::quick().with_groups("group", "treatment");
        let params = StressParams {
            spike_fraction: 0.05,
            fold_change: 2.0,
            mode: SpikeMode::Compositional,
            n_spiked: 5,
        };

        let result = execute_stress_run(&counts, &metadata, &params, &config, 0, &mock_pipeline);
        assert!(result.is_ok());

        let run = result.unwrap();
        assert_eq!(run.replicate, 0);
        assert!(run.sensitivity >= 0.0 && run.sensitivity <= 1.0);
        assert!(run.fdr >= 0.0 && run.fdr <= 1.0);
        assert!(run.runtime_seconds > 0.0);
    }

    #[test]
    fn test_stress_summary_display() {
        let config = StressConfig::quick();
        let summary = StressSummary {
            config: config.clone(),
            runs: vec![],
            aggregated: vec![
                AggregatedMetrics {
                    params: StressParams {
                        spike_fraction: 0.05,
                        fold_change: 2.0,
                        mode: SpikeMode::Compositional,
                        n_spiked: 5,
                    },
                    n_replicates: 2,
                    sensitivity_mean: 0.80,
                    sensitivity_std: 0.05,
                    fdr_mean: 0.10,
                    fdr_std: 0.02,
                    f1_mean: 0.84,
                    f1_std: 0.03,
                    attenuation_mean: 0.95,
                    attenuation_std: 0.02,
                    permuted_fpr_mean: 0.05,
                },
            ],
            power_curves: HashMap::new(),
            mode_comparison: ModeComparisonSummary {
                recommended_mode: SpikeMode::Compositional,
                highest_sensitivity_mode: SpikeMode::Compositional,
                lowest_fdr_mode: SpikeMode::Compositional,
                mode_stats: HashMap::new(),
            },
            total_runtime_seconds: 10.0,
        };

        let display = format!("{}", summary);
        assert!(display.contains("Compositional Stress Test Results"));
        assert!(display.contains("5% spiked"));
    }

    #[test]
    fn test_stress_summary_csv() {
        let config = StressConfig::quick();
        let summary = StressSummary {
            config: config.clone(),
            runs: vec![],
            aggregated: vec![
                AggregatedMetrics {
                    params: StressParams {
                        spike_fraction: 0.05,
                        fold_change: 2.0,
                        mode: SpikeMode::Compositional,
                        n_spiked: 5,
                    },
                    n_replicates: 2,
                    sensitivity_mean: 0.80,
                    sensitivity_std: 0.05,
                    fdr_mean: 0.10,
                    fdr_std: 0.02,
                    f1_mean: 0.84,
                    f1_std: 0.03,
                    attenuation_mean: 0.95,
                    attenuation_std: 0.02,
                    permuted_fpr_mean: 0.05,
                },
            ],
            power_curves: HashMap::new(),
            mode_comparison: ModeComparisonSummary {
                recommended_mode: SpikeMode::Compositional,
                highest_sensitivity_mode: SpikeMode::Compositional,
                lowest_fdr_mode: SpikeMode::Compositional,
                mode_stats: HashMap::new(),
            },
            total_runtime_seconds: 10.0,
        };

        let csv = summary.to_csv();
        assert!(csv.contains("spike_fraction,n_spiked,fold_change,mode"));
        assert!(csv.contains("0.0500,5,2.00,Compositional"));
    }
}
