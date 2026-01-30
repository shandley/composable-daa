//! Spike-in normalization for absolute abundance estimation.
//!
//! Uses a constant spike-in (internal standard) to normalize counts to
//! absolute-like abundances. This addresses a fundamental limitation of
//! compositional data: the inability to distinguish relative from absolute
//! changes.
//!
//! # Theory
//!
//! When a spike-in organism is added at constant absolute amount to all samples:
//! - Its relative proportion varies inversely with total bacterial load
//! - Low spike proportion → high total load
//! - High spike proportion → low total load
//!
//! By normalizing to the spike-in, we recover absolute-like abundances:
//! ```text
//! normalized_ij = count_ij / spike_proportion_j
//!               = count_ij * library_size_j / spike_count_j
//! ```
//!
//! # Example
//!
//! ```ignore
//! use composable_daa::normalize::spikein::norm_spikein;
//!
//! // Normalize using Salinibacter as the spike-in
//! let result = norm_spikein(&counts, "Salinibacter_ruber")?;
//!
//! // The spike-in feature is removed from output by default
//! assert!(!result.feature_ids.contains(&"Salinibacter_ruber".to_string()));
//! ```
//!
//! # References
//!
//! Stammler et al. (2016) "Adjusting microbiome profiles for differences in
//! microbial load by spike-in bacteria" Microbiome 4:28

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// Result of spike-in normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikeinMatrix {
    /// The normalized data (features × samples).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers (spike-in removed unless keep_spikein=true).
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// The spike-in feature used for normalization.
    pub spikein_id: String,
    /// Spike-in counts per sample (before normalization).
    pub spikein_counts: Vec<u64>,
    /// Spike-in proportions per sample.
    pub spikein_proportions: Vec<f64>,
    /// Estimated relative total load per sample (inverse of spike proportion).
    pub estimated_load: Vec<f64>,
    /// Original library sizes.
    pub library_sizes: Vec<u64>,
}

impl SpikeinMatrix {
    /// Get the normalized value for a feature and sample.
    pub fn get(&self, feature: usize, sample: usize) -> f64 {
        self.data[(feature, sample)]
    }

    /// Number of features.
    pub fn n_features(&self) -> usize {
        self.data.nrows()
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.data.ncols()
    }

    /// Get a row (feature) as a vector.
    pub fn row(&self, feature: usize) -> Vec<f64> {
        self.data.row(feature).iter().cloned().collect()
    }

    /// Get a column (sample) as a vector.
    pub fn col(&self, sample: usize) -> Vec<f64> {
        self.data.column(sample).iter().cloned().collect()
    }

    /// Get reference to the underlying matrix.
    pub fn matrix(&self) -> &DMatrix<f64> {
        &self.data
    }

    /// Convert to TransformedMatrix for compatibility with downstream analysis.
    pub fn to_transformed(&self) -> super::clr::TransformedMatrix {
        super::clr::TransformedMatrix {
            data: self.data.clone(),
            feature_ids: self.feature_ids.clone(),
            sample_ids: self.sample_ids.clone(),
            transformation: format!("Spike-in normalized (ref={})", self.spikein_id),
            geometric_means: vec![],
        }
    }
}

/// Configuration for spike-in normalization.
#[derive(Debug, Clone)]
pub struct SpikeinConfig {
    /// ID of the spike-in feature to use as reference.
    pub spikein_id: String,
    /// Whether to keep the spike-in in the output (default: false).
    pub keep_spikein: bool,
    /// Minimum spike-in proportion to accept (default: 0.001 = 0.1%).
    /// Samples with lower proportions may have spike-in detection issues.
    pub min_proportion: f64,
    /// Whether to log-transform the normalized values (default: false).
    pub log_transform: bool,
    /// Pseudocount for log transformation (default: 1.0).
    pub log_pseudocount: f64,
}

impl Default for SpikeinConfig {
    fn default() -> Self {
        Self {
            spikein_id: String::new(),
            keep_spikein: false,
            min_proportion: 0.001,
            log_transform: false,
            log_pseudocount: 1.0,
        }
    }
}

impl SpikeinConfig {
    /// Create a new config with the spike-in feature ID.
    pub fn new(spikein_id: &str) -> Self {
        Self {
            spikein_id: spikein_id.to_string(),
            ..Default::default()
        }
    }

    /// Keep the spike-in feature in the output.
    pub fn keep_spikein(mut self) -> Self {
        self.keep_spikein = true;
        self
    }

    /// Set minimum acceptable spike-in proportion.
    pub fn min_proportion(mut self, min: f64) -> Self {
        self.min_proportion = min;
        self
    }

    /// Apply log transformation after normalization.
    pub fn log_transform(mut self, pseudocount: f64) -> Self {
        self.log_transform = true;
        self.log_pseudocount = pseudocount;
        self
    }
}

/// Apply spike-in normalization to a count matrix.
///
/// This uses a constant spike-in as an internal standard to convert
/// relative abundances to absolute-like abundances.
///
/// # Arguments
/// * `counts` - Count matrix
/// * `spikein_id` - Feature ID of the spike-in to use as reference
///
/// # Returns
/// A SpikeinMatrix with normalized values. The spike-in feature is removed
/// from the output by default.
pub fn norm_spikein(counts: &CountMatrix, spikein_id: &str) -> Result<SpikeinMatrix> {
    norm_spikein_with_config(counts, &SpikeinConfig::new(spikein_id))
}

/// Apply spike-in normalization with custom configuration.
pub fn norm_spikein_with_config(
    counts: &CountMatrix,
    config: &SpikeinConfig,
) -> Result<SpikeinMatrix> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply spike-in normalization to empty matrix".to_string(),
        ));
    }

    // Find the spike-in feature
    let spikein_idx = counts
        .feature_ids()
        .iter()
        .position(|id| id == &config.spikein_id)
        .ok_or_else(|| {
            DaaError::InvalidParameter(format!(
                "Spike-in feature '{}' not found in count matrix",
                config.spikein_id
            ))
        })?;

    // Get library sizes and spike-in counts
    let library_sizes = counts.col_sums();
    let spikein_counts: Vec<u64> = (0..n_samples)
        .map(|j| counts.get(spikein_idx, j))
        .collect();

    // Calculate spike-in proportions
    let spikein_proportions: Vec<f64> = spikein_counts
        .iter()
        .zip(library_sizes.iter())
        .map(|(&spike, &lib)| spike as f64 / lib as f64)
        .collect();

    // Check for low spike-in proportions
    for (j, &prop) in spikein_proportions.iter().enumerate() {
        if prop < config.min_proportion {
            return Err(DaaError::Numerical(format!(
                "Sample '{}' has spike-in proportion {:.4}% which is below minimum {:.4}%. \
                 This may indicate spike-in detection issues.",
                counts.sample_ids()[j],
                prop * 100.0,
                config.min_proportion * 100.0
            )));
        }
    }

    // Check for zero spike-in counts
    for (j, &count) in spikein_counts.iter().enumerate() {
        if count == 0 {
            return Err(DaaError::Numerical(format!(
                "Sample '{}' has zero spike-in counts, cannot normalize",
                counts.sample_ids()[j]
            )));
        }
    }

    // Calculate estimated relative total load (inverse of spike proportion)
    let estimated_load: Vec<f64> = spikein_proportions.iter().map(|&p| 1.0 / p).collect();

    // Determine which features to include in output
    let output_features: Vec<usize> = if config.keep_spikein {
        (0..n_features).collect()
    } else {
        (0..n_features).filter(|&i| i != spikein_idx).collect()
    };

    let output_feature_ids: Vec<String> = output_features
        .iter()
        .map(|&i| counts.feature_ids()[i].clone())
        .collect();

    // Normalize: divide by spike-in proportion (multiply by estimated load)
    let n_output_features = output_features.len();
    let mut data = DMatrix::zeros(n_output_features, n_samples);

    for (out_i, &in_i) in output_features.iter().enumerate() {
        for j in 0..n_samples {
            let count = counts.get(in_i, j) as f64;
            let mut normalized = count * estimated_load[j];

            if config.log_transform {
                normalized = (normalized + config.log_pseudocount).ln();
            }

            data[(out_i, j)] = normalized;
        }
    }

    Ok(SpikeinMatrix {
        data,
        feature_ids: output_feature_ids,
        sample_ids: counts.sample_ids().to_vec(),
        spikein_id: config.spikein_id.clone(),
        spikein_counts,
        spikein_proportions,
        estimated_load,
        library_sizes,
    })
}

/// Detect potential spike-in features based on naming patterns.
///
/// Returns a list of feature IDs that match common spike-in naming patterns.
pub fn detect_spikein_candidates(counts: &CountMatrix) -> Vec<String> {
    let patterns = [
        "spike", "spikein", "spike-in", "spike_in",
        "internal_standard", "synthetic",
        "salinibacter", "rhizobium", "alicyclobacillus",  // Common spike-in organisms
        "thermus", "imtechella",  // Other spike-ins
    ];

    counts
        .feature_ids()
        .iter()
        .filter(|id| {
            let lower = id.to_lowercase();
            patterns.iter().any(|p| lower.contains(p))
        })
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::TriMat;

    fn create_test_counts_with_spikein() -> CountMatrix {
        // 4 features × 3 samples
        // Feature 0: SPIKEIN (constant absolute amount, varies in proportion)
        // Feature 1-3: regular features
        let mut tri_mat = TriMat::new((4, 3));

        // Sample 0: low total load (spike is 10% = 100/1000)
        tri_mat.add_triplet(0, 0, 100);  // SPIKEIN
        tri_mat.add_triplet(1, 0, 400);
        tri_mat.add_triplet(2, 0, 300);
        tri_mat.add_triplet(3, 0, 200);
        // Total = 1000, spike = 10%

        // Sample 1: medium total load (spike is 5% = 100/2000)
        tri_mat.add_triplet(0, 1, 100);  // SPIKEIN (same absolute amount)
        tri_mat.add_triplet(1, 1, 900);
        tri_mat.add_triplet(2, 1, 600);
        tri_mat.add_triplet(3, 1, 400);
        // Total = 2000, spike = 5%

        // Sample 2: high total load (spike is 2% = 100/5000)
        tri_mat.add_triplet(0, 2, 100);  // SPIKEIN (same absolute amount)
        tri_mat.add_triplet(1, 2, 2400);
        tri_mat.add_triplet(2, 2, 1500);
        tri_mat.add_triplet(3, 2, 1000);
        // Total = 5000, spike = 2%

        let feature_ids = vec![
            "SPIKEIN".into(),
            "TaxonA".into(),
            "TaxonB".into(),
            "TaxonC".into(),
        ];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_spikein_normalization() {
        let counts = create_test_counts_with_spikein();
        let result = norm_spikein(&counts, "SPIKEIN").unwrap();

        // Spike-in should be removed from output
        assert_eq!(result.n_features(), 3);
        assert!(!result.feature_ids.contains(&"SPIKEIN".to_string()));

        // Check spike-in proportions
        assert_relative_eq!(result.spikein_proportions[0], 0.10, epsilon = 1e-10);
        assert_relative_eq!(result.spikein_proportions[1], 0.05, epsilon = 1e-10);
        assert_relative_eq!(result.spikein_proportions[2], 0.02, epsilon = 1e-10);

        // Check estimated load (inverse of proportion)
        assert_relative_eq!(result.estimated_load[0], 10.0, epsilon = 1e-10);
        assert_relative_eq!(result.estimated_load[1], 20.0, epsilon = 1e-10);
        assert_relative_eq!(result.estimated_load[2], 50.0, epsilon = 1e-10);
    }

    #[test]
    fn test_spikein_absolute_recovery() {
        let counts = create_test_counts_with_spikein();
        let result = norm_spikein(&counts, "SPIKEIN").unwrap();

        // TaxonA: If it truly has constant absolute abundance across samples,
        // after spike-in normalization, values should be similar
        //
        // Raw counts: 400, 900, 2400 (relative, varies with total load)
        // After normalization: 400*10, 900*20, 2400*50 = 4000, 18000, 120000
        //
        // Wait - in this test data, TaxonA is NOT constant absolute abundance.
        // It's varying with total load (compositional effect).
        //
        // Let's verify the math is correct:
        let taxa_idx = result.feature_ids.iter().position(|x| x == "TaxonA").unwrap();
        assert_relative_eq!(result.get(taxa_idx, 0), 400.0 * 10.0, epsilon = 1e-6);
        assert_relative_eq!(result.get(taxa_idx, 1), 900.0 * 20.0, epsilon = 1e-6);
        assert_relative_eq!(result.get(taxa_idx, 2), 2400.0 * 50.0, epsilon = 1e-6);
    }

    #[test]
    fn test_keep_spikein() {
        let counts = create_test_counts_with_spikein();
        let config = SpikeinConfig::new("SPIKEIN").keep_spikein();
        let result = norm_spikein_with_config(&counts, &config).unwrap();

        // Spike-in should be kept
        assert_eq!(result.n_features(), 4);
        assert!(result.feature_ids.contains(&"SPIKEIN".to_string()));
    }

    #[test]
    fn test_spikein_not_found() {
        let counts = create_test_counts_with_spikein();
        let result = norm_spikein(&counts, "NONEXISTENT");
        assert!(result.is_err());
    }

    #[test]
    fn test_detect_candidates() {
        let mut tri_mat = TriMat::new((5, 2));
        for i in 0..5 {
            tri_mat.add_triplet(i, 0, 10);
            tri_mat.add_triplet(i, 1, 20);
        }

        let feature_ids = vec![
            "Bacteroides".into(),
            "Salinibacter_ruber".into(),
            "spike_in_control".into(),
            "E_coli".into(),
            "Alicyclobacillus".into(),
        ];
        let sample_ids = vec!["S1".into(), "S2".into()];
        let counts = CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap();

        let candidates = detect_spikein_candidates(&counts);
        assert!(candidates.contains(&"Salinibacter_ruber".to_string()));
        assert!(candidates.contains(&"spike_in_control".to_string()));
        assert!(candidates.contains(&"Alicyclobacillus".to_string()));
        assert!(!candidates.contains(&"Bacteroides".to_string()));
        assert!(!candidates.contains(&"E_coli".to_string()));
    }

    #[test]
    fn test_log_transform() {
        let counts = create_test_counts_with_spikein();
        let config = SpikeinConfig::new("SPIKEIN").log_transform(1.0);
        let result = norm_spikein_with_config(&counts, &config).unwrap();

        // Values should be log-transformed
        let taxa_idx = result.feature_ids.iter().position(|x| x == "TaxonA").unwrap();
        let expected = (400.0_f64 * 10.0 + 1.0).ln();
        assert_relative_eq!(result.get(taxa_idx, 0), expected, epsilon = 1e-6);
    }

    #[test]
    fn test_to_transformed() {
        let counts = create_test_counts_with_spikein();
        let result = norm_spikein(&counts, "SPIKEIN").unwrap();
        let transformed = result.to_transformed();

        assert_eq!(transformed.n_features(), 3);
        assert!(transformed.transformation.contains("Spike-in"));
    }
}
