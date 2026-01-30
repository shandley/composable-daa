//! Permutation tests for non-parametric hypothesis testing.
//!
//! Permutation tests provide a distribution-free approach to hypothesis testing
//! by shuffling group labels and computing a null distribution of test statistics.
//! This is particularly useful for microbiome data where distributional assumptions
//! may not hold.
//!
//! # Algorithm
//!
//! 1. Compute the observed test statistic for each feature
//! 2. Shuffle group labels n_permutations times
//! 3. For each permutation, recompute test statistics
//! 4. P-value = proportion of permuted statistics >= observed (two-sided)
//!
//! # Example
//!
//! ```ignore
//! use composable_daa::test::permutation::{test_permutation, PermutationConfig};
//!
//! let config = PermutationConfig::default();
//! let results = test_permutation(&transformed, &metadata, "group", "grouptreatment", &config)?;
//! ```

use crate::data::{DesignMatrix, Formula, Metadata, Variable};
use crate::error::{DaaError, Result};
use crate::model::lm::model_lm;
use crate::normalize::TransformedMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Configuration for permutation testing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PermutationConfig {
    /// Number of permutations to run.
    pub n_permutations: usize,
    /// Random seed for reproducibility.
    pub seed: u64,
    /// Whether to use parallel computation.
    pub parallel: bool,
    /// Minimum count of non-zero values to test a feature.
    pub min_nonzero: usize,
}

impl Default for PermutationConfig {
    fn default() -> Self {
        Self {
            n_permutations: 1000,
            seed: 42,
            parallel: true,
            min_nonzero: 3,
        }
    }
}

impl PermutationConfig {
    /// Create a quick configuration for testing (fewer permutations).
    pub fn quick() -> Self {
        Self {
            n_permutations: 100,
            ..Default::default()
        }
    }

    /// Create a thorough configuration (more permutations).
    pub fn thorough() -> Self {
        Self {
            n_permutations: 10000,
            ..Default::default()
        }
    }
}

/// Result of a permutation test for a single feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PermutationResult {
    /// Feature identifier.
    pub feature_id: String,
    /// Observed test statistic.
    pub observed_stat: f64,
    /// Two-sided permutation p-value.
    pub p_value: f64,
    /// Number of permutations used.
    pub n_permutations: usize,
    /// Number of permutations with |stat| >= |observed|.
    pub n_extreme: usize,
    /// Estimated effect size (coefficient).
    pub estimate: f64,
    /// Standard error of the estimate.
    pub std_error: f64,
}

/// Results of permutation testing across all features.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PermutationResults {
    /// Results for each feature.
    pub results: Vec<PermutationResult>,
    /// Configuration used.
    pub config: PermutationConfig,
    /// Coefficient that was tested.
    pub coefficient: String,
    /// Formula used.
    pub formula: String,
}

impl PermutationResults {
    /// Get the number of features tested.
    pub fn len(&self) -> usize {
        self.results.len()
    }

    /// Check if results are empty.
    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }

    /// Get results sorted by p-value.
    pub fn sorted_by_pvalue(&self) -> Vec<&PermutationResult> {
        let mut sorted: Vec<_> = self.results.iter().collect();
        sorted.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());
        sorted
    }

    /// Get significant results at a given threshold.
    pub fn significant(&self, alpha: f64) -> Vec<&PermutationResult> {
        self.results.iter().filter(|r| r.p_value < alpha).collect()
    }

    /// Apply Benjamini-Hochberg correction and return q-values.
    pub fn with_bh_correction(&self) -> Vec<(String, f64, f64)> {
        let mut indexed: Vec<(usize, f64)> = self
            .results
            .iter()
            .enumerate()
            .map(|(i, r)| (i, r.p_value))
            .collect();

        // Sort by p-value
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let n = indexed.len();
        let mut q_values = vec![0.0; n];

        // Calculate BH q-values
        let mut min_q: f64 = 1.0;
        for (rank, &(orig_idx, p)) in indexed.iter().rev().enumerate() {
            let adj_rank = n - rank;
            let q = (p * n as f64) / adj_rank as f64;
            min_q = min_q.min(q);
            q_values[orig_idx] = min_q.min(1.0);
        }

        // Return (feature_id, p_value, q_value)
        self.results
            .iter()
            .zip(q_values)
            .map(|(r, q)| (r.feature_id.clone(), r.p_value, q))
            .collect()
    }
}

/// Simple deterministic random number generator for permutations.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        // xorshift64
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    /// Fisher-Yates shuffle
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        let n = slice.len();
        for i in (1..n).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

/// Run permutation test on normalized data.
///
/// # Arguments
/// * `transformed` - Normalized data matrix
/// * `metadata` - Sample metadata
/// * `formula` - Model formula (e.g., "~ group")
/// * `coefficient` - Coefficient to test
/// * `config` - Permutation test configuration
///
/// # Returns
/// PermutationResults containing p-values for each feature.
pub fn test_permutation(
    transformed: &TransformedMatrix,
    metadata: &Metadata,
    formula: &str,
    coefficient: &str,
    config: &PermutationConfig,
) -> Result<PermutationResults> {
    let n_samples = transformed.n_samples();
    let n_features = transformed.n_features();

    if n_samples < 3 {
        return Err(DaaError::InvalidParameter(
            "Permutation test requires at least 3 samples".to_string(),
        ));
    }

    // Parse formula and create design matrix
    let parsed_formula = Formula::parse(formula)?;
    let design = DesignMatrix::from_formula(metadata, &parsed_formula)?;

    // Find the coefficient index
    let coef_idx = design
        .coefficient_names()
        .iter()
        .position(|name| name == coefficient)
        .ok_or_else(|| {
            DaaError::InvalidParameter(format!(
                "Coefficient '{}' not found. Available: {:?}",
                coefficient,
                design.coefficient_names()
            ))
        })?;

    // Fit original model to get observed statistics
    let original_fit = model_lm(transformed, &design)?;

    // Get observed statistics for all features
    let observed_stats: Vec<(f64, f64, f64)> = (0..n_features)
        .map(|i| {
            let fit = &original_fit.fits[i];
            let coef = fit.coefficients[coef_idx];
            let se = fit.std_errors[coef_idx];
            let t_stat = if se > 0.0 { coef / se } else { 0.0 };
            (t_stat, coef, se)
        })
        .collect();

    // Generate permuted design matrices
    let sample_indices: Vec<usize> = (0..n_samples).collect();

    // Run permutations
    let permuted_stats: Vec<Vec<f64>> = if config.parallel {
        (0..config.n_permutations)
            .into_par_iter()
            .map(|perm_idx| {
                let mut rng = SimpleRng::new(config.seed.wrapping_add(perm_idx as u64));
                let mut indices = sample_indices.clone();
                rng.shuffle(&mut indices);

                // Create permuted metadata
                let permuted_metadata = permute_metadata(metadata, &indices).unwrap();
                let permuted_design =
                    DesignMatrix::from_formula(&permuted_metadata, &parsed_formula).unwrap();

                // Fit model on permuted data
                let fit = model_lm(transformed, &permuted_design).unwrap();

                // Extract test statistics
                (0..n_features)
                    .map(|i| {
                        let f = &fit.fits[i];
                        let coef = f.coefficients[coef_idx];
                        let se = f.std_errors[coef_idx];
                        if se > 0.0 {
                            coef / se
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect()
    } else {
        let mut rng = SimpleRng::new(config.seed);
        (0..config.n_permutations)
            .map(|_| {
                let mut indices = sample_indices.clone();
                rng.shuffle(&mut indices);

                let permuted_metadata = permute_metadata(metadata, &indices).unwrap();
                let permuted_design =
                    DesignMatrix::from_formula(&permuted_metadata, &parsed_formula).unwrap();

                let fit = model_lm(transformed, &permuted_design).unwrap();

                (0..n_features)
                    .map(|i| {
                        let f = &fit.fits[i];
                        let coef = f.coefficients[coef_idx];
                        let se = f.std_errors[coef_idx];
                        if se > 0.0 {
                            coef / se
                        } else {
                            0.0
                        }
                    })
                    .collect()
            })
            .collect()
    };

    // Calculate p-values for each feature
    let results: Vec<PermutationResult> = (0..n_features)
        .map(|i| {
            let (obs_stat, estimate, std_error) = observed_stats[i];
            let abs_obs = obs_stat.abs();

            // Count permutations with |stat| >= |observed|
            let n_extreme: usize = permuted_stats
                .iter()
                .filter(|perm| perm[i].abs() >= abs_obs)
                .count();

            // Two-sided p-value with +1 correction for observed
            let p_value = (n_extreme as f64 + 1.0) / (config.n_permutations as f64 + 1.0);

            PermutationResult {
                feature_id: transformed.feature_ids[i].clone(),
                observed_stat: obs_stat,
                p_value,
                n_permutations: config.n_permutations,
                n_extreme,
                estimate,
                std_error,
            }
        })
        .collect();

    Ok(PermutationResults {
        results,
        config: config.clone(),
        coefficient: coefficient.to_string(),
        formula: formula.to_string(),
    })
}

/// Create metadata with permuted variable values.
///
/// This shuffles the values of all variables according to the given indices,
/// effectively permuting which sample has which variable values.
fn permute_metadata(metadata: &Metadata, indices: &[usize]) -> Result<Metadata> {
    let sample_ids = metadata.sample_ids();

    // For each column, create permuted values
    let mut result = metadata.clone();

    for col in metadata.column_names() {
        // Get original values in order
        let original_values: Vec<Option<&Variable>> = sample_ids
            .iter()
            .map(|id| metadata.get(id, col).map(|v| v))
            .collect();

        // Permute the values according to indices
        let permuted_values: Vec<Option<&Variable>> = indices
            .iter()
            .map(|&i| original_values.get(i).copied().flatten())
            .collect();

        result = result.with_shuffled_column(col, &permuted_values)?;
    }

    Ok(result)
}

/// Convenience function for quick permutation test.
pub fn test_permutation_quick(
    transformed: &TransformedMatrix,
    metadata: &Metadata,
    formula: &str,
    coefficient: &str,
) -> Result<PermutationResults> {
    test_permutation(
        transformed,
        metadata,
        formula,
        coefficient,
        &PermutationConfig::quick(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::CountMatrix;
    use crate::normalize::norm_clr;
    use crate::zero::pseudocount::add_pseudocount;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        writeln!(file, "S0\tcontrol").unwrap();
        writeln!(file, "S1\tcontrol").unwrap();
        writeln!(file, "S2\tcontrol").unwrap();
        writeln!(file, "S3\tcontrol").unwrap();
        writeln!(file, "S4\ttreatment").unwrap();
        writeln!(file, "S5\ttreatment").unwrap();
        writeln!(file, "S6\ttreatment").unwrap();
        writeln!(file, "S7\ttreatment").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    fn create_test_data() -> (TransformedMatrix, Metadata) {
        // Create transformed data directly with known effects
        // This avoids issues with CLR normalization artifacts
        use nalgebra::DMatrix;

        // 4 features × 8 samples (4 control, 4 treatment)
        // Feature 0: No group effect (values around 0)
        // Feature 1: Strong group effect (control ~-1, treatment ~+1)
        // Feature 2: Moderate group effect
        // Feature 3: No group effect
        let data = DMatrix::from_row_slice(4, 8, &[
            // F0: no effect
            0.1, -0.1, 0.2, -0.2, 0.0, 0.1, -0.1, 0.0,
            // F1: strong effect (control negative, treatment positive)
            -1.0, -1.2, -0.8, -1.1, 1.2, 1.0, 0.9, 1.1,
            // F2: moderate effect
            -0.3, -0.5, -0.4, -0.2, 0.4, 0.5, 0.3, 0.6,
            // F3: no effect
            -0.1, 0.2, 0.0, 0.1, -0.2, 0.1, 0.0, -0.1,
        ]);

        let transformed = TransformedMatrix {
            data,
            feature_ids: (0..4).map(|i| format!("F{}", i)).collect(),
            sample_ids: (0..8).map(|i| format!("S{}", i)).collect(),
            transformation: "test".to_string(),
            geometric_means: vec![1.0; 8],
        };

        // Create metadata
        let metadata = create_test_metadata();

        (transformed, metadata)
    }

    #[test]
    fn test_permutation_basic() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig {
            n_permutations: 100,
            seed: 42,
            parallel: false,
            min_nonzero: 1,
        };

        let results = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        assert_eq!(results.len(), 4);
        assert_eq!(results.coefficient, "grouptreatment");
    }

    #[test]
    fn test_permutation_pvalues_bounded() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig::quick();

        let results = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        // All p-values should be between 0 and 1
        for r in &results.results {
            assert!(r.p_value > 0.0 && r.p_value <= 1.0);
        }
    }

    #[test]
    fn test_permutation_detects_signal() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig {
            n_permutations: 500,
            seed: 42,
            parallel: true,
            min_nonzero: 1,
        };

        let results = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        // Feature 1 has strong signal, should have lower p-value than Feature 0
        let f1 = results.results.iter().find(|r| r.feature_id == "F1").unwrap();
        let f0 = results.results.iter().find(|r| r.feature_id == "F0").unwrap();

        // The feature with signal should have a lower p-value than the null feature
        assert!(
            f1.p_value < f0.p_value,
            "Feature with signal (p={}) should have lower p-value than null feature (p={})",
            f1.p_value, f0.p_value
        );

        // The observed statistic should be larger for the signal feature
        assert!(
            f1.observed_stat.abs() > f0.observed_stat.abs(),
            "Feature with signal should have larger test statistic"
        );
    }

    #[test]
    fn test_permutation_reproducible() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig {
            n_permutations: 100,
            seed: 12345,
            parallel: false,
            min_nonzero: 1,
        };

        let results1 = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        let results2 = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        // Same seed should give same p-values
        for (r1, r2) in results1.results.iter().zip(results2.results.iter()) {
            assert_eq!(r1.p_value, r2.p_value);
        }
    }

    #[test]
    fn test_bh_correction() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig::quick();

        let results = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        let corrected = results.with_bh_correction();

        // All q-values should be >= corresponding p-values
        for (_, p, q) in &corrected {
            assert!(*q >= *p, "q-value should be >= p-value");
            assert!(*q <= 1.0, "q-value should be <= 1.0");
        }
    }

    #[test]
    fn test_sorted_results() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig::quick();

        let results = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "grouptreatment",
            &config,
        )
        .unwrap();

        let sorted = results.sorted_by_pvalue();

        // Should be sorted by ascending p-value
        for i in 1..sorted.len() {
            assert!(sorted[i].p_value >= sorted[i - 1].p_value);
        }
    }

    #[test]
    fn test_permute_metadata() {
        // Create metadata with 4 samples
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        writeln!(file, "S0\tA").unwrap();
        writeln!(file, "S1\tA").unwrap();
        writeln!(file, "S2\tB").unwrap();
        writeln!(file, "S3\tB").unwrap();
        file.flush().unwrap();
        let metadata = Metadata::from_tsv(file.path()).unwrap();

        let indices = vec![3, 2, 1, 0]; // Reverse order
        let permuted = permute_metadata(&metadata, &indices).unwrap();

        // After permutation with indices [3, 2, 1, 0]:
        // S0 should have S3's group value (B)
        // S1 should have S2's group value (B)
        // S2 should have S1's group value (A)
        // S3 should have S0's group value (A)
        let s0_group = permuted.get("S0", "group").unwrap();
        let s3_group = permuted.get("S3", "group").unwrap();
        assert_eq!(s0_group.as_categorical(), Some("B"));
        assert_eq!(s3_group.as_categorical(), Some("A"));
    }

    #[test]
    fn test_invalid_coefficient() {
        let (transformed, metadata) = create_test_data();
        let config = PermutationConfig::quick();

        let result = test_permutation(
            &transformed,
            &metadata,
            "~ group",
            "nonexistent",
            &config,
        );

        assert!(result.is_err());
    }
}
