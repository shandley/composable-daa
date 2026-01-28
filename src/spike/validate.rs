//! Spike-in validation: run full validation experiments with positive and negative controls.

use crate::data::{CountMatrix, DaResultSet, Metadata};
use crate::error::{DaaError, Result};
use crate::spike::evaluate::{evaluate_spikes, SpikeEvaluation};
use crate::spike::types::SpikeSelection;
use crate::spike::{spike_abundance, spike_presence};
use serde::{Deserialize, Serialize};

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

    /// Shuffle a vector in place.
    fn shuffle<T>(&mut self, vec: &mut [T]) {
        for i in (1..vec.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            vec.swap(i, j);
        }
    }
}

/// Configuration for spike-in validation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationConfig {
    /// Number of features to spike.
    pub n_spike: usize,
    /// Fold change for abundance spikes.
    pub fold_change: f64,
    /// Target prevalence increase for presence spikes (0-1).
    pub prevalence_increase: f64,
    /// Group column in metadata.
    pub group_column: String,
    /// Target group for spiking.
    pub target_group: String,
    /// FDR threshold for calling significance.
    pub fdr_threshold: f64,
    /// How to select features for spiking.
    pub selection: SpikeSelection,
    /// Number of permutation iterations for null distribution.
    pub n_permutations: usize,
    /// Base random seed.
    pub seed: u64,
}

impl Default for ValidationConfig {
    fn default() -> Self {
        Self {
            n_spike: 50,
            fold_change: 2.0,
            prevalence_increase: 0.3,
            group_column: "group".into(),
            target_group: "treatment".into(),
            fdr_threshold: 0.05,
            selection: SpikeSelection::Random,
            n_permutations: 10,
            seed: 42,
        }
    }
}

/// Result from a single validation run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationResult {
    /// Evaluation from spiked data (positive control).
    pub spiked_eval: SpikeEvaluation,
    /// Evaluations from permuted data (negative controls).
    pub permuted_evals: Vec<PermutedEvaluation>,
    /// Configuration used.
    pub config: ValidationConfig,
}

/// Evaluation from a permutation (negative control).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PermutedEvaluation {
    /// Permutation iteration number.
    pub iteration: usize,
    /// Number of features called significant (should be ~FDR * n_tested).
    pub n_significant: usize,
    /// Empirical false positive rate.
    pub false_positive_rate: f64,
}

impl ValidationResult {
    /// Check if the validation passed basic criteria.
    pub fn passed(&self) -> bool {
        // Spiked should have reasonable sensitivity
        let good_sensitivity = self.spiked_eval.sensitivity >= 0.5;

        // Permuted should have low FP rate (FDR controlled)
        let fdr_controlled = self.mean_permuted_fpr() <= self.config.fdr_threshold * 2.0;

        good_sensitivity && fdr_controlled
    }

    /// Mean false positive rate across permutations.
    pub fn mean_permuted_fpr(&self) -> f64 {
        if self.permuted_evals.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.permuted_evals.iter().map(|e| e.false_positive_rate).sum();
        sum / self.permuted_evals.len() as f64
    }
}

impl std::fmt::Display for ValidationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Spike-in Validation Results")?;
        writeln!(f, "===========================")?;
        writeln!(f, "\nPositive Control (spiked data):")?;
        writeln!(f, "  Sensitivity: {:.1}%", self.spiked_eval.sensitivity * 100.0)?;
        writeln!(f, "  FDR:         {:.1}%", self.spiked_eval.fdr * 100.0)?;
        writeln!(f, "  Precision:   {:.1}%", self.spiked_eval.precision * 100.0)?;

        if !self.permuted_evals.is_empty() {
            writeln!(f, "\nNegative Controls ({} permutations):", self.permuted_evals.len())?;
            writeln!(f, "  Mean FP rate: {:.1}%", self.mean_permuted_fpr() * 100.0)?;
            let max_fpr = self.permuted_evals.iter()
                .map(|e| e.false_positive_rate)
                .fold(0.0, f64::max);
            writeln!(f, "  Max FP rate:  {:.1}%", max_fpr * 100.0)?;
        }

        writeln!(f, "\nValidation: {}", if self.passed() { "PASSED" } else { "FAILED" })?;
        Ok(())
    }
}

/// Summary across multiple validation runs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationSummary {
    /// Individual validation results.
    pub results: Vec<ValidationResult>,
    /// Mean sensitivity across runs.
    pub mean_sensitivity: f64,
    /// Mean FDR across runs.
    pub mean_fdr: f64,
    /// Mean permutation FP rate.
    pub mean_permuted_fpr: f64,
    /// Fraction of runs that passed.
    pub pass_rate: f64,
}

impl ValidationSummary {
    /// Create summary from validation results.
    pub fn from_results(results: Vec<ValidationResult>) -> Self {
        if results.is_empty() {
            return Self {
                results: vec![],
                mean_sensitivity: 0.0,
                mean_fdr: 0.0,
                mean_permuted_fpr: 0.0,
                pass_rate: 0.0,
            };
        }

        let n = results.len() as f64;
        let mean_sensitivity = results.iter().map(|r| r.spiked_eval.sensitivity).sum::<f64>() / n;
        let mean_fdr = results.iter().map(|r| r.spiked_eval.fdr).sum::<f64>() / n;
        let mean_permuted_fpr = results.iter().map(|r| r.mean_permuted_fpr()).sum::<f64>() / n;
        let pass_rate = results.iter().filter(|r| r.passed()).count() as f64 / n;

        Self {
            results,
            mean_sensitivity,
            mean_fdr,
            mean_permuted_fpr,
            pass_rate,
        }
    }
}

impl std::fmt::Display for ValidationSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Spike-in Validation Summary ({} runs)", self.results.len())?;
        writeln!(f, "=======================================")?;
        writeln!(f, "  Mean Sensitivity:   {:.1}%", self.mean_sensitivity * 100.0)?;
        writeln!(f, "  Mean FDR:           {:.1}%", self.mean_fdr * 100.0)?;
        writeln!(f, "  Mean Permuted FPR:  {:.1}%", self.mean_permuted_fpr * 100.0)?;
        writeln!(f, "  Pass Rate:          {:.1}%", self.pass_rate * 100.0)?;
        Ok(())
    }
}

/// Shuffle group labels in metadata for negative control analysis.
///
/// This creates a permuted version of the metadata where group assignments
/// are randomized, breaking any true association. Running the pipeline on
/// permuted data should yield ~FDR false positives.
///
/// # Arguments
/// * `metadata` - Original sample metadata
/// * `group_column` - Column containing group assignments to shuffle
/// * `seed` - Random seed for reproducibility
///
/// # Returns
/// New Metadata with shuffled group labels.
pub fn shuffle_labels(metadata: &Metadata, group_column: &str, seed: u64) -> Result<Metadata> {
    let mut rng = SimpleRng::new(seed);

    // Get all sample IDs and their group values
    let sample_ids = metadata.sample_ids();
    let mut group_values: Vec<_> = sample_ids
        .iter()
        .map(|id| metadata.get(id, group_column))
        .collect();

    // Verify all samples have the group column
    if group_values.iter().any(|v| v.is_none()) {
        return Err(DaaError::InvalidParameter(format!(
            "Not all samples have column '{}'",
            group_column
        )));
    }

    // Shuffle the group values
    rng.shuffle(&mut group_values);

    // Create new metadata with shuffled labels
    let shuffled = metadata.with_shuffled_column(group_column, &group_values)?;

    Ok(shuffled)
}

/// Run a complete spike-in validation experiment.
///
/// This function:
/// 1. Spikes known effects into the data
/// 2. Runs the analysis pipeline on spiked data
/// 3. Evaluates detection performance (positive control)
/// 4. Runs permutations as negative controls
///
/// # Arguments
/// * `counts` - Original count matrix
/// * `metadata` - Sample metadata
/// * `config` - Validation configuration
/// * `run_pipeline` - Function that takes (counts, metadata) and returns DA results
///
/// # Returns
/// ValidationResult with evaluation metrics.
pub fn run_spike_validation<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    config: &ValidationConfig,
    run_pipeline: F,
) -> Result<ValidationResult>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet>,
{
    // 1. Spike the data (positive control)
    let spiked = spike_abundance(
        counts,
        metadata,
        &config.group_column,
        config.n_spike,
        config.fold_change,
        &config.target_group,
        config.selection.clone(),
        config.seed,
    )?;

    // 2. Run pipeline on spiked data
    let spiked_results = run_pipeline(&spiked.counts, metadata)?;

    // 3. Evaluate spike detection
    let spiked_eval = evaluate_spikes(&spiked_results, &spiked.spec, config.fdr_threshold);

    // 4. Run permutations (negative controls)
    let mut permuted_evals = Vec::new();
    for i in 0..config.n_permutations {
        let perm_seed = config.seed.wrapping_add(1000 + i as u64);

        // Shuffle labels
        let shuffled_meta = shuffle_labels(metadata, &config.group_column, perm_seed)?;

        // Run pipeline on original data with shuffled labels
        let perm_results = run_pipeline(counts, &shuffled_meta)?;

        // Count significant results
        let n_significant = perm_results
            .results
            .iter()
            .filter(|r| r.q_value < config.fdr_threshold)
            .count();

        let n_tested = perm_results.len();
        let false_positive_rate = if n_tested > 0 {
            n_significant as f64 / n_tested as f64
        } else {
            0.0
        };

        permuted_evals.push(PermutedEvaluation {
            iteration: i,
            n_significant,
            false_positive_rate,
        });
    }

    Ok(ValidationResult {
        spiked_eval,
        permuted_evals,
        config: config.clone(),
    })
}

/// Run validation with presence spikes instead of abundance.
pub fn run_spike_validation_presence<F>(
    counts: &CountMatrix,
    metadata: &Metadata,
    config: &ValidationConfig,
    run_pipeline: F,
) -> Result<ValidationResult>
where
    F: Fn(&CountMatrix, &Metadata) -> Result<DaResultSet>,
{
    use crate::spike::types::AbundanceLevel;

    // 1. Spike the data with presence changes
    let spiked = spike_presence(
        counts,
        metadata,
        &config.group_column,
        config.n_spike,
        config.prevalence_increase,
        &config.target_group,
        AbundanceLevel::Median,
        config.selection.clone(),
        config.seed,
    )?;

    // 2. Run pipeline on spiked data
    let spiked_results = run_pipeline(&spiked.counts, metadata)?;

    // 3. Evaluate spike detection
    let spiked_eval = evaluate_spikes(&spiked_results, &spiked.spec, config.fdr_threshold);

    // 4. Run permutations (same as abundance validation)
    let mut permuted_evals = Vec::new();
    for i in 0..config.n_permutations {
        let perm_seed = config.seed.wrapping_add(1000 + i as u64);
        let shuffled_meta = shuffle_labels(metadata, &config.group_column, perm_seed)?;
        let perm_results = run_pipeline(counts, &shuffled_meta)?;

        let n_significant = perm_results
            .results
            .iter()
            .filter(|r| r.q_value < config.fdr_threshold)
            .count();

        let n_tested = perm_results.len();
        let false_positive_rate = if n_tested > 0 {
            n_significant as f64 / n_tested as f64
        } else {
            0.0
        };

        permuted_evals.push(PermutedEvaluation {
            iteration: i,
            n_significant,
            false_positive_rate,
        });
    }

    Ok(ValidationResult {
        spiked_eval,
        permuted_evals,
        config: config.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::Variable;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        let mut tri_mat = TriMat::new((10, 20));
        // Create features with varying prevalence
        for feat in 0..10 {
            for sample in 0..20 {
                // Higher feature index = lower prevalence
                if sample < (20 - feat * 2) {
                    tri_mat.add_triplet(feat, sample, 100 - feat as u64 * 5);
                }
            }
        }
        let feature_ids: Vec<String> = (0..10).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..20).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        for i in 0..10 {
            writeln!(file, "S{}\tcontrol\t{}", i, 25 + i).unwrap();
        }
        for i in 10..20 {
            writeln!(file, "S{}\ttreatment\t{}", i, 30 + i).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_shuffle_labels() {
        let metadata = create_test_metadata();

        let shuffled = shuffle_labels(&metadata, "group", 42).unwrap();

        // Same samples should exist
        assert_eq!(shuffled.sample_ids().len(), metadata.sample_ids().len());

        // Group values should be permuted (at least some different positions)
        let mut different_count = 0;
        for id in metadata.sample_ids() {
            let orig = metadata.get(id, "group");
            let shuf = shuffled.get(id, "group");
            if orig != shuf {
                different_count += 1;
            }
        }
        // With 20 samples, very unlikely all stay the same
        assert!(different_count > 0, "Shuffling should change at least some labels");
    }

    #[test]
    fn test_shuffle_labels_preserves_counts() {
        let metadata = create_test_metadata();

        let shuffled = shuffle_labels(&metadata, "group", 42).unwrap();

        // Count controls and treatments in original
        let orig_controls = metadata.sample_ids().iter()
            .filter(|id| {
                matches!(metadata.get(id, "group"), Some(Variable::Categorical(g)) if g == "control")
            })
            .count();

        // Count in shuffled
        let shuf_controls = shuffled.sample_ids().iter()
            .filter(|id| {
                matches!(shuffled.get(id, "group"), Some(Variable::Categorical(g)) if g == "control")
            })
            .count();

        // Should have same number of each group
        assert_eq!(orig_controls, shuf_controls);
    }

    #[test]
    fn test_validation_config_default() {
        let config = ValidationConfig::default();

        assert_eq!(config.n_spike, 50);
        assert_eq!(config.fold_change, 2.0);
        assert_eq!(config.fdr_threshold, 0.05);
        assert_eq!(config.n_permutations, 10);
    }

    #[test]
    fn test_validation_summary() {
        // Create mock validation results
        let config = ValidationConfig::default();

        let result1 = ValidationResult {
            spiked_eval: SpikeEvaluation {
                true_positives: 40,
                false_positives: 5,
                false_negatives: 10,
                true_negatives: 945,
                sensitivity: 0.8,
                specificity: 0.995,
                precision: 0.889,
                fdr: 0.111,
                f1_score: 0.842,
                effect_correlation: 0.9,
                effect_bias: 0.1,
                effect_mae: 0.2,
                by_tier: Default::default(),
                n_spiked: 50,
                n_tested: 1000,
                fdr_threshold: 0.05,
            },
            permuted_evals: vec![
                PermutedEvaluation { iteration: 0, n_significant: 30, false_positive_rate: 0.03 },
                PermutedEvaluation { iteration: 1, n_significant: 40, false_positive_rate: 0.04 },
            ],
            config: config.clone(),
        };

        let summary = ValidationSummary::from_results(vec![result1]);

        assert!((summary.mean_sensitivity - 0.8).abs() < 0.001);
        assert!((summary.mean_fdr - 0.111).abs() < 0.001);
        assert!((summary.mean_permuted_fpr - 0.035).abs() < 0.001);
        assert!((summary.pass_rate - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_permuted_evaluation() {
        let eval = PermutedEvaluation {
            iteration: 0,
            n_significant: 50,
            false_positive_rate: 0.05,
        };

        assert_eq!(eval.iteration, 0);
        assert_eq!(eval.n_significant, 50);
        assert!((eval.false_positive_rate - 0.05).abs() < 0.001);
    }
}
