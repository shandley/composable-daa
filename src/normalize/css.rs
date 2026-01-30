//! Cumulative Sum Scaling (CSS) normalization from metagenomeSeq.
//!
//! CSS normalization addresses the sparse, zero-inflated nature of microbiome data
//! by calculating scaling factors based on the cumulative sum of counts up to a
//! percentile determined by the data. This makes it robust to:
//! - High sparsity (many zeros)
//! - Differential sequencing depth
//! - Highly abundant features that dominate library sizes
//!
//! # Algorithm
//!
//! 1. For each sample, sort non-zero counts in ascending order
//! 2. Compute cumulative sums
//! 3. Find the quantile where the cumulative sum captures most variation
//! 4. Use cumulative sum up to that quantile as the scaling factor
//! 5. Divide counts by scaling factor (and optionally log-transform)
//!
//! # Reference
//!
//! Paulson JN, Stine OC, Bravo HC, Pop M (2013). "Differential abundance analysis
//! for microbial marker-gene surveys." Nature Methods 10:1200-1202.

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use crate::normalize::clr::TransformedMatrix;
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Configuration for CSS normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CssConfig {
    /// Quantile for determining the scaling reference (default: 0.5).
    /// Higher values use more of the count distribution.
    pub quantile: f64,
    /// Whether to log-transform after scaling (default: true).
    pub log_transform: bool,
    /// Pseudocount to add before log transformation (default: 1.0).
    /// Only used if log_transform is true.
    pub log_pseudocount: f64,
    /// Minimum number of non-zero features required per sample.
    pub min_nonzero: usize,
}

impl Default for CssConfig {
    fn default() -> Self {
        Self {
            quantile: 0.5,
            log_transform: true,
            log_pseudocount: 1.0,
            min_nonzero: 1,
        }
    }
}

impl CssConfig {
    /// Create config for raw CSS (no log transform).
    pub fn raw() -> Self {
        Self {
            log_transform: false,
            ..Default::default()
        }
    }

    /// Create config with custom quantile.
    pub fn with_quantile(quantile: f64) -> Self {
        Self {
            quantile,
            ..Default::default()
        }
    }
}

/// Result of CSS normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CssMatrix {
    /// The normalized data (features x samples).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers.
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// CSS scaling factors per sample.
    pub scaling_factors: Vec<f64>,
    /// Quantile used for normalization.
    pub quantile: f64,
    /// Whether data is log-transformed.
    pub log_transformed: bool,
    /// Percentile cutoffs per sample (for diagnostics).
    pub percentile_cutoffs: Vec<f64>,
}

impl CssMatrix {
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

    /// Convert to TransformedMatrix for pipeline compatibility.
    pub fn to_transformed_matrix(&self) -> TransformedMatrix {
        TransformedMatrix {
            data: self.data.clone(),
            feature_ids: self.feature_ids.clone(),
            sample_ids: self.sample_ids.clone(),
            transformation: if self.log_transformed {
                format!("CSS(q={:.2},log)", self.quantile)
            } else {
                format!("CSS(q={:.2})", self.quantile)
            },
            geometric_means: self.scaling_factors.clone(),
        }
    }
}

/// Apply CSS normalization with default settings.
///
/// Uses quantile=0.5 and applies log transformation.
///
/// # Arguments
/// * `counts` - The count matrix to normalize
///
/// # Returns
/// A CssMatrix with normalized values.
///
/// # Example
/// ```ignore
/// use composable_daa::normalize::css::norm_css;
///
/// let normalized = norm_css(&counts)?;
/// println!("Scaling factors: {:?}", normalized.scaling_factors);
/// ```
pub fn norm_css(counts: &CountMatrix) -> Result<CssMatrix> {
    norm_css_with_config(counts, &CssConfig::default())
}

/// Apply CSS normalization with custom configuration.
///
/// # Arguments
/// * `counts` - The count matrix to normalize
/// * `config` - Configuration options
///
/// # Returns
/// A CssMatrix with normalized values.
pub fn norm_css_with_config(counts: &CountMatrix, config: &CssConfig) -> Result<CssMatrix> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply CSS to empty matrix".to_string(),
        ));
    }

    if !(0.0..=1.0).contains(&config.quantile) {
        return Err(DaaError::InvalidParameter(
            "CSS quantile must be between 0 and 1".to_string(),
        ));
    }

    // Calculate CSS scaling factors for each sample
    let css_results: Vec<Result<(f64, f64)>> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            // Get non-zero counts for this sample
            let mut nonzero_counts: Vec<u64> = (0..n_features)
                .filter_map(|i| {
                    let count = counts.get(i, j);
                    if count > 0 { Some(count) } else { None }
                })
                .collect();

            if nonzero_counts.len() < config.min_nonzero {
                return Err(DaaError::EmptyData(format!(
                    "Sample {} has fewer than {} non-zero features",
                    j, config.min_nonzero
                )));
            }

            // Sort counts in ascending order
            nonzero_counts.sort_unstable();

            // Compute cumulative sums
            let cumsum: Vec<u64> = nonzero_counts
                .iter()
                .scan(0u64, |acc, &x| {
                    *acc += x;
                    Some(*acc)
                })
                .collect();

            // Find the index corresponding to the quantile
            let quantile_idx = ((nonzero_counts.len() as f64 * config.quantile).ceil() as usize)
                .saturating_sub(1)
                .min(nonzero_counts.len() - 1);

            // Scaling factor is the cumulative sum at the quantile
            let scaling_factor = cumsum[quantile_idx] as f64;

            // Percentile cutoff (the count value at the quantile)
            let percentile_cutoff = nonzero_counts[quantile_idx] as f64;

            Ok((scaling_factor, percentile_cutoff))
        })
        .collect();

    // Check for errors and extract results
    let mut scaling_factors = Vec::with_capacity(n_samples);
    let mut percentile_cutoffs = Vec::with_capacity(n_samples);

    for result in css_results {
        let (sf, pc) = result?;
        scaling_factors.push(sf);
        percentile_cutoffs.push(pc);
    }

    // Compute median scaling factor for normalization reference
    let mut sorted_factors = scaling_factors.clone();
    sorted_factors.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_factor = if sorted_factors.len() % 2 == 0 {
        (sorted_factors[sorted_factors.len() / 2 - 1] + sorted_factors[sorted_factors.len() / 2])
            / 2.0
    } else {
        sorted_factors[sorted_factors.len() / 2]
    };

    // Normalize scaling factors relative to median
    let norm_factors: Vec<f64> = scaling_factors
        .iter()
        .map(|&sf| if median_factor > 0.0 { sf / median_factor } else { 1.0 })
        .collect();

    // Apply normalization
    let mut normalized = DMatrix::zeros(n_features, n_samples);

    for j in 0..n_samples {
        let factor = norm_factors[j];
        for i in 0..n_features {
            let count = counts.get(i, j) as f64;
            let scaled = if factor > 0.0 { count / factor } else { count };

            normalized[(i, j)] = if config.log_transform {
                (scaled + config.log_pseudocount).ln()
            } else {
                scaled
            };
        }
    }

    Ok(CssMatrix {
        data: normalized,
        feature_ids: counts.feature_ids().to_vec(),
        sample_ids: counts.sample_ids().to_vec(),
        scaling_factors: norm_factors,
        quantile: config.quantile,
        log_transformed: config.log_transform,
        percentile_cutoffs,
    })
}

/// Calculate CSS scaling factors without applying normalization.
///
/// Useful for diagnostics or custom normalization workflows.
pub fn css_factors(counts: &CountMatrix, quantile: f64) -> Result<Vec<f64>> {
    let config = CssConfig {
        quantile,
        log_transform: false,
        ..Default::default()
    };
    let result = norm_css_with_config(counts, &config)?;
    Ok(result.scaling_factors)
}

/// Estimate optimal quantile for CSS normalization.
///
/// Finds the quantile that maximizes stability of scaling factors across samples.
/// Based on the metagenomeSeq approach of finding where the cumulative sum
/// distribution stabilizes.
///
/// # Arguments
/// * `counts` - The count matrix
/// * `quantiles` - Quantiles to evaluate (default: 0.1 to 0.9 in steps of 0.05)
///
/// # Returns
/// The optimal quantile value.
pub fn estimate_css_quantile(counts: &CountMatrix, quantiles: Option<&[f64]>) -> Result<f64> {
    let default_quantiles: Vec<f64> = (1..=18).map(|i| i as f64 * 0.05).collect();
    let quantiles = quantiles.unwrap_or(&default_quantiles);

    if quantiles.is_empty() {
        return Err(DaaError::InvalidParameter(
            "Must provide at least one quantile to evaluate".to_string(),
        ));
    }

    // For each quantile, calculate the coefficient of variation of scaling factors
    let cvs: Vec<(f64, f64)> = quantiles
        .iter()
        .filter_map(|&q| {
            let factors = css_factors(counts, q).ok()?;
            let mean = factors.iter().sum::<f64>() / factors.len() as f64;
            let variance = factors.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                / factors.len() as f64;
            let cv = if mean > 0.0 { variance.sqrt() / mean } else { f64::INFINITY };
            Some((q, cv))
        })
        .collect();

    if cvs.is_empty() {
        return Err(DaaError::Numerical(
            "Could not calculate CV for any quantile".to_string(),
        ));
    }

    // Find quantile with minimum CV (most stable scaling factors)
    cvs.into_iter()
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .map(|(q, _)| q)
        .ok_or_else(|| DaaError::Numerical("Failed to find optimal quantile".to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::TriMat;

    fn create_test_counts() -> CountMatrix {
        // 5 features x 4 samples with varying sparsity
        let mut tri_mat = TriMat::new((5, 4));

        // Sample 0: moderate counts
        tri_mat.add_triplet(0, 0, 100);
        tri_mat.add_triplet(1, 0, 200);
        tri_mat.add_triplet(2, 0, 50);
        tri_mat.add_triplet(3, 0, 150);
        // Feature 4 is zero

        // Sample 1: higher counts
        tri_mat.add_triplet(0, 1, 200);
        tri_mat.add_triplet(1, 1, 400);
        tri_mat.add_triplet(2, 1, 100);
        tri_mat.add_triplet(3, 1, 300);
        tri_mat.add_triplet(4, 1, 50);

        // Sample 2: lower counts
        tri_mat.add_triplet(0, 2, 50);
        tri_mat.add_triplet(1, 2, 100);
        tri_mat.add_triplet(2, 2, 25);
        // Features 3, 4 are zero

        // Sample 3: very sparse
        tri_mat.add_triplet(0, 3, 20);
        tri_mat.add_triplet(1, 3, 30);
        // Features 2, 3, 4 are zero

        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..4).map(|i| format!("S{}", i)).collect();

        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_css_basic() {
        let counts = create_test_counts();
        let result = norm_css(&counts).unwrap();

        assert_eq!(result.n_features(), 5);
        assert_eq!(result.n_samples(), 4);
        assert_eq!(result.scaling_factors.len(), 4);
        assert!(result.log_transformed);
    }

    #[test]
    fn test_css_scaling_factors() {
        let counts = create_test_counts();
        let result = norm_css_with_config(&counts, &CssConfig::raw()).unwrap();

        // Scaling factors should be positive
        for sf in &result.scaling_factors {
            assert!(*sf > 0.0, "Scaling factor should be positive");
        }

        // Median scaling factor should be ~1.0 after normalization
        let mut sorted = result.scaling_factors.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = if sorted.len() % 2 == 0 {
            (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
        } else {
            sorted[sorted.len() / 2]
        };
        assert_relative_eq!(median, 1.0, epsilon = 0.1);
    }

    #[test]
    fn test_css_preserves_zeros() {
        let counts = create_test_counts();
        let result = norm_css_with_config(&counts, &CssConfig::raw()).unwrap();

        // Zero counts should remain zero (before log transform)
        // Feature 4 is zero in sample 0
        assert_eq!(result.get(4, 0), 0.0);
    }

    #[test]
    fn test_css_log_transform() {
        let counts = create_test_counts();

        // Without log transform
        let raw = norm_css_with_config(&counts, &CssConfig::raw()).unwrap();
        assert!(!raw.log_transformed);

        // With log transform (default)
        let logged = norm_css(&counts).unwrap();
        assert!(logged.log_transformed);

        // Log-transformed values should be smaller for non-zero counts
        // (since log(x+1) < x for x > ~1.7)
        let raw_val = raw.get(1, 1); // High count
        let log_val = logged.get(1, 1);
        assert!(log_val < raw_val, "Log transform should reduce large values");
    }

    #[test]
    fn test_css_quantile_effect() {
        let counts = create_test_counts();

        let low_q = norm_css_with_config(
            &counts,
            &CssConfig {
                quantile: 0.25,
                log_transform: false,
                ..Default::default()
            },
        )
        .unwrap();

        let high_q = norm_css_with_config(
            &counts,
            &CssConfig {
                quantile: 0.75,
                log_transform: false,
                ..Default::default()
            },
        )
        .unwrap();

        // Different quantiles should produce different scaling factors
        assert_ne!(low_q.scaling_factors, high_q.scaling_factors);
    }

    #[test]
    fn test_css_factors() {
        let counts = create_test_counts();
        let factors = css_factors(&counts, 0.5).unwrap();

        assert_eq!(factors.len(), 4);
        for f in &factors {
            assert!(*f > 0.0);
        }
    }

    #[test]
    fn test_estimate_css_quantile() {
        let counts = create_test_counts();
        let optimal_q = estimate_css_quantile(&counts, None).unwrap();

        assert!(optimal_q >= 0.05 && optimal_q <= 0.95);
    }

    #[test]
    fn test_css_to_transformed_matrix() {
        let counts = create_test_counts();
        let css = norm_css(&counts).unwrap();
        let transformed = css.to_transformed_matrix();

        assert_eq!(transformed.n_features(), 5);
        assert_eq!(transformed.n_samples(), 4);
        assert!(transformed.transformation.contains("CSS"));
    }

    #[test]
    fn test_css_invalid_quantile() {
        let counts = create_test_counts();

        let result = norm_css_with_config(
            &counts,
            &CssConfig {
                quantile: 1.5, // Invalid
                ..Default::default()
            },
        );
        assert!(result.is_err());

        let result = norm_css_with_config(
            &counts,
            &CssConfig {
                quantile: -0.1, // Invalid
                ..Default::default()
            },
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_css_config_presets() {
        let raw = CssConfig::raw();
        assert!(!raw.log_transform);
        assert_eq!(raw.quantile, 0.5);

        let custom = CssConfig::with_quantile(0.75);
        assert!(custom.log_transform);
        assert_eq!(custom.quantile, 0.75);
    }

    #[test]
    fn test_css_percentile_cutoffs() {
        let counts = create_test_counts();
        let result = norm_css(&counts).unwrap();

        // Percentile cutoffs should be stored
        assert_eq!(result.percentile_cutoffs.len(), 4);

        // Cutoffs should be positive
        for pc in &result.percentile_cutoffs {
            assert!(*pc > 0.0);
        }
    }

    #[test]
    fn test_css_sample_with_few_nonzeros() {
        // Create a very sparse sample
        let mut tri_mat = TriMat::new((5, 2));
        tri_mat.add_triplet(0, 0, 100);
        tri_mat.add_triplet(1, 0, 200);
        tri_mat.add_triplet(0, 1, 50); // Only one non-zero

        let counts = CountMatrix::new(
            tri_mat.to_csr(),
            vec!["A".into(), "B".into(), "C".into(), "D".into(), "E".into()],
            vec!["S1".into(), "S2".into()],
        )
        .unwrap();

        // Should still work with min_nonzero = 1
        let result = norm_css(&counts);
        assert!(result.is_ok());

        // Should fail with higher min_nonzero requirement
        let result = norm_css_with_config(
            &counts,
            &CssConfig {
                min_nonzero: 3,
                ..Default::default()
            },
        );
        assert!(result.is_err());
    }
}
