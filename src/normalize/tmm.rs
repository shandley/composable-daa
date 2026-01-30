//! Trimmed Mean of M-values (TMM) normalization.
//!
//! TMM is a normalization method from edgeR that calculates scaling factors
//! to account for compositional differences between samples. Unlike simple
//! library size normalization, TMM is robust to asymmetric differential
//! expression where a subset of features dominate expression changes.
//!
//! # Algorithm
//!
//! 1. Select a reference sample (by default, the one with upper quartile
//!    closest to the mean upper quartile)
//! 2. For each sample, calculate M-values (log-ratios) and A-values
//!    (average expression) for each feature
//! 3. Trim extreme values (default: 30% for M, 5% for A)
//! 4. Calculate weighted mean of M-values as the normalization factor
//!
//! # Reference
//!
//! Robinson MD, Oshlack A. A scaling normalization method for differential
//! expression analysis of RNA-seq data. Genome Biology 11, R25 (2010).

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Configuration for TMM normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TmmConfig {
    /// Fraction of M-values to trim from each tail (default: 0.30).
    pub trim_m: f64,
    /// Fraction of A-values to trim from each tail (default: 0.05).
    pub trim_a: f64,
    /// Minimum count for a feature to be used in TMM calculation.
    pub min_count: u64,
    /// Whether to log-transform the output (default: false).
    pub log_transform: bool,
    /// Pseudocount for log transformation.
    pub log_pseudocount: f64,
    /// Reference sample index (None = auto-select).
    pub reference_sample: Option<usize>,
}

impl Default for TmmConfig {
    fn default() -> Self {
        Self {
            trim_m: 0.30,
            trim_a: 0.05,
            min_count: 1,
            log_transform: false,
            log_pseudocount: 1.0,
            reference_sample: None,
        }
    }
}

/// Result of TMM normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TmmMatrix {
    /// The normalized data (features × samples).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers.
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// TMM normalization factors for each sample.
    pub norm_factors: Vec<f64>,
    /// Effective library sizes (library_size * norm_factor).
    pub effective_lib_sizes: Vec<f64>,
    /// Original library sizes.
    pub library_sizes: Vec<u64>,
    /// Index of the reference sample used.
    pub reference_sample: usize,
    /// Whether data is log-transformed.
    pub log_transformed: bool,
}

impl TmmMatrix {
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
            transformation: if self.log_transformed {
                "TMM (log)".to_string()
            } else {
                "TMM".to_string()
            },
            geometric_means: vec![],
        }
    }
}

/// Apply TMM normalization with default parameters.
///
/// # Example
/// ```ignore
/// let tmm = norm_tmm(&counts)?;
/// let factors = &tmm.norm_factors;  // Use these to scale counts
/// ```
pub fn norm_tmm(counts: &CountMatrix) -> Result<TmmMatrix> {
    norm_tmm_with_config(counts, &TmmConfig::default())
}

/// Apply TMM normalization with custom configuration.
pub fn norm_tmm_with_config(counts: &CountMatrix, config: &TmmConfig) -> Result<TmmMatrix> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply TMM to empty matrix".to_string(),
        ));
    }

    if n_samples < 2 {
        return Err(DaaError::InvalidParameter(
            "TMM requires at least 2 samples".to_string(),
        ));
    }

    if config.trim_m < 0.0 || config.trim_m >= 0.5 {
        return Err(DaaError::InvalidParameter(
            "trim_m must be in [0, 0.5)".to_string(),
        ));
    }

    if config.trim_a < 0.0 || config.trim_a >= 0.5 {
        return Err(DaaError::InvalidParameter(
            "trim_a must be in [0, 0.5)".to_string(),
        ));
    }

    // Get library sizes
    let library_sizes = counts.col_sums();

    // Check for zero library sizes
    for (j, &lib_size) in library_sizes.iter().enumerate() {
        if lib_size == 0 {
            return Err(DaaError::Numerical(format!(
                "Sample {} has zero total counts",
                counts.sample_ids()[j]
            )));
        }
    }

    // Select reference sample
    let ref_idx = config.reference_sample.unwrap_or_else(|| {
        select_reference_sample(counts, &library_sizes)
    });

    if ref_idx >= n_samples {
        return Err(DaaError::InvalidParameter(format!(
            "Reference sample index {} out of bounds (n_samples = {})",
            ref_idx, n_samples
        )));
    }

    // Calculate TMM factors for each sample
    let norm_factors: Vec<f64> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            if j == ref_idx {
                1.0 // Reference sample has factor 1.0
            } else {
                calculate_tmm_factor(
                    counts,
                    j,
                    ref_idx,
                    library_sizes[j],
                    library_sizes[ref_idx],
                    config,
                )
            }
        })
        .collect();

    // Calculate effective library sizes
    let effective_lib_sizes: Vec<f64> = library_sizes
        .iter()
        .zip(&norm_factors)
        .map(|(&lib, &factor)| lib as f64 * factor)
        .collect();

    // Normalize counts: count / effective_lib_size * mean_lib_size
    let mean_lib_size: f64 = library_sizes.iter().map(|&x| x as f64).sum::<f64>() / n_samples as f64;

    let normalized_cols: Vec<Vec<f64>> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            let eff_lib = effective_lib_sizes[j];
            (0..n_features)
                .map(|i| {
                    let count = counts.get(i, j) as f64;
                    let normalized = count / eff_lib * mean_lib_size;
                    if config.log_transform {
                        (normalized + config.log_pseudocount).ln()
                    } else {
                        normalized
                    }
                })
                .collect()
        })
        .collect();

    // Build the normalized matrix
    let mut data = DMatrix::zeros(n_features, n_samples);
    for (j, col) in normalized_cols.iter().enumerate() {
        for (i, &val) in col.iter().enumerate() {
            data[(i, j)] = val;
        }
    }

    Ok(TmmMatrix {
        data,
        feature_ids: counts.feature_ids().to_vec(),
        sample_ids: counts.sample_ids().to_vec(),
        norm_factors,
        effective_lib_sizes,
        library_sizes,
        reference_sample: ref_idx,
        log_transformed: config.log_transform,
    })
}

/// Select reference sample as the one with upper quartile closest to mean.
fn select_reference_sample(counts: &CountMatrix, library_sizes: &[u64]) -> usize {
    let n_samples = counts.n_samples();
    let n_features = counts.n_features();

    // Calculate upper quartile (75th percentile) for each sample
    let upper_quartiles: Vec<f64> = (0..n_samples)
        .map(|j| {
            let lib_size = library_sizes[j] as f64;
            let mut proportions: Vec<f64> = (0..n_features)
                .map(|i| counts.get(i, j) as f64 / lib_size)
                .filter(|&p| p > 0.0)
                .collect();

            if proportions.is_empty() {
                return 0.0;
            }

            proportions.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let q75_idx = (proportions.len() as f64 * 0.75) as usize;
            proportions[q75_idx.min(proportions.len() - 1)]
        })
        .collect();

    // Find sample closest to mean upper quartile
    let mean_uq: f64 = upper_quartiles.iter().sum::<f64>() / n_samples as f64;

    upper_quartiles
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let diff_a = (*a - mean_uq).abs();
            let diff_b = (*b - mean_uq).abs();
            diff_a.partial_cmp(&diff_b).unwrap()
        })
        .map(|(idx, _)| idx)
        .unwrap_or(0)
}

/// Calculate TMM factor for a sample relative to reference.
fn calculate_tmm_factor(
    counts: &CountMatrix,
    sample_idx: usize,
    ref_idx: usize,
    sample_lib: u64,
    ref_lib: u64,
    config: &TmmConfig,
) -> f64 {
    let n_features = counts.n_features();
    let sample_lib = sample_lib as f64;
    let ref_lib = ref_lib as f64;

    // Collect M and A values for features present in both samples
    let mut ma_values: Vec<(f64, f64, f64)> = Vec::new(); // (M, A, weight)

    for i in 0..n_features {
        let count_s = counts.get(i, sample_idx);
        let count_r = counts.get(i, ref_idx);

        // Skip features below minimum count in either sample
        if count_s < config.min_count || count_r < config.min_count {
            continue;
        }

        let prop_s = count_s as f64 / sample_lib;
        let prop_r = count_r as f64 / ref_lib;

        // M-value: log2 ratio
        let m = (prop_s / prop_r).log2();

        // A-value: average log2 expression
        let a = 0.5 * (prop_s * prop_r).log2();

        // Weight: inverse asymptotic variance
        // Var(M) ≈ (N-Y)/NY + (N'-Y')/N'Y' for counts Y, Y' and totals N, N'
        let var_m = (sample_lib - count_s as f64) / (sample_lib * count_s as f64)
            + (ref_lib - count_r as f64) / (ref_lib * count_r as f64);
        let weight = if var_m > 0.0 { 1.0 / var_m } else { 0.0 };

        if m.is_finite() && a.is_finite() && weight > 0.0 {
            ma_values.push((m, a, weight));
        }
    }

    if ma_values.is_empty() {
        return 1.0; // No valid features, return neutral factor
    }

    // Sort by M for M-trimming
    ma_values.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let n = ma_values.len();
    let trim_m_count = (n as f64 * config.trim_m) as usize;

    // Trim extreme M values
    let m_trimmed: Vec<(f64, f64, f64)> = if 2 * trim_m_count < n {
        ma_values[trim_m_count..(n - trim_m_count)].to_vec()
    } else {
        ma_values.clone()
    };

    if m_trimmed.is_empty() {
        return 1.0;
    }

    // Sort by A for A-trimming
    let mut a_sorted = m_trimmed.clone();
    a_sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let n = a_sorted.len();
    let trim_a_count = (n as f64 * config.trim_a) as usize;

    // Trim extreme A values
    let final_values: Vec<(f64, f64, f64)> = if 2 * trim_a_count < n {
        a_sorted[trim_a_count..(n - trim_a_count)].to_vec()
    } else {
        a_sorted
    };

    if final_values.is_empty() {
        return 1.0;
    }

    // Calculate weighted mean of M values
    let sum_weighted_m: f64 = final_values.iter().map(|(m, _, w)| m * w).sum();
    let sum_weights: f64 = final_values.iter().map(|(_, _, w)| w).sum();

    if sum_weights <= 0.0 {
        return 1.0;
    }

    let weighted_mean_m = sum_weighted_m / sum_weights;

    // TMM factor is 2^weighted_mean_M
    2.0_f64.powf(weighted_mean_m)
}

/// Get just the TMM normalization factors without normalizing the matrix.
///
/// This is useful when you want to apply the factors yourself or use them
/// with other tools like edgeR or DESeq2.
pub fn tmm_factors(counts: &CountMatrix) -> Result<Vec<f64>> {
    let tmm = norm_tmm(counts)?;
    Ok(tmm.norm_factors)
}

/// Get TMM factors with custom configuration.
pub fn tmm_factors_with_config(counts: &CountMatrix, config: &TmmConfig) -> Result<Vec<f64>> {
    let tmm = norm_tmm_with_config(counts, config)?;
    Ok(tmm.norm_factors)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::TriMat;

    fn create_test_counts() -> CountMatrix {
        // 4 features × 3 samples with different library sizes
        let mut tri_mat = TriMat::new((4, 3));

        // Sample 0: lib size 1000
        tri_mat.add_triplet(0, 0, 500);
        tri_mat.add_triplet(1, 0, 300);
        tri_mat.add_triplet(2, 0, 150);
        tri_mat.add_triplet(3, 0, 50);

        // Sample 1: lib size 2000 (double, same proportions)
        tri_mat.add_triplet(0, 1, 1000);
        tri_mat.add_triplet(1, 1, 600);
        tri_mat.add_triplet(2, 1, 300);
        tri_mat.add_triplet(3, 1, 100);

        // Sample 2: lib size 500 (half, same proportions)
        tri_mat.add_triplet(0, 2, 250);
        tri_mat.add_triplet(1, 2, 150);
        tri_mat.add_triplet(2, 2, 75);
        tri_mat.add_triplet(3, 2, 25);

        let feature_ids = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_asymmetric_counts() -> CountMatrix {
        // Create data where one sample has inflated high-abundance features
        let mut tri_mat = TriMat::new((4, 2));

        // Sample 0: baseline
        tri_mat.add_triplet(0, 0, 500);
        tri_mat.add_triplet(1, 0, 300);
        tri_mat.add_triplet(2, 0, 150);
        tri_mat.add_triplet(3, 0, 50);

        // Sample 1: feature A is 3x higher, others same (asymmetric)
        tri_mat.add_triplet(0, 1, 1500); // 3x
        tri_mat.add_triplet(1, 1, 300);
        tri_mat.add_triplet(2, 1, 150);
        tri_mat.add_triplet(3, 1, 50);

        let feature_ids = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let sample_ids = vec!["S1".into(), "S2".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_tmm_basic() {
        let counts = create_test_counts();
        let tmm = norm_tmm(&counts).unwrap();

        assert_eq!(tmm.n_features(), 4);
        assert_eq!(tmm.n_samples(), 3);
        assert_eq!(tmm.library_sizes, vec![1000, 2000, 500]);
    }

    #[test]
    fn test_tmm_factors_symmetric() {
        // When proportions are identical, TMM factors should be ~1.0
        let counts = create_test_counts();
        let tmm = norm_tmm(&counts).unwrap();

        for factor in &tmm.norm_factors {
            assert_relative_eq!(*factor, 1.0, epsilon = 0.1);
        }
    }

    #[test]
    fn test_tmm_factors_asymmetric() {
        // When one feature is inflated, TMM should detect this
        let counts = create_asymmetric_counts();
        let tmm = norm_tmm(&counts).unwrap();

        // Sample 1 has inflated feature A (1500 vs 500), causing higher lib size
        // TMM factor adjusts for this - the factor != 1.0 indicates asymmetry detected
        // Factor < 1 means sample has inflated counts relative to most features
        let factor_diff = (tmm.norm_factors[1] - 1.0).abs();
        assert!(factor_diff > 0.1, "TMM should detect asymmetric inflation");
    }

    #[test]
    fn test_tmm_reference_sample() {
        let counts = create_test_counts();

        // Force reference to be sample 0
        let config = TmmConfig {
            reference_sample: Some(0),
            ..Default::default()
        };
        let tmm = norm_tmm_with_config(&counts, &config).unwrap();

        assert_eq!(tmm.reference_sample, 0);
        assert_relative_eq!(tmm.norm_factors[0], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_tmm_effective_lib_sizes() {
        let counts = create_test_counts();
        let tmm = norm_tmm(&counts).unwrap();

        for j in 0..3 {
            let expected = tmm.library_sizes[j] as f64 * tmm.norm_factors[j];
            assert_relative_eq!(tmm.effective_lib_sizes[j], expected, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_tmm_with_log_transform() {
        let counts = create_test_counts();
        let config = TmmConfig {
            log_transform: true,
            log_pseudocount: 1.0,
            ..Default::default()
        };
        let tmm = norm_tmm_with_config(&counts, &config).unwrap();

        assert!(tmm.log_transformed);

        // All values should be finite (log of positive values)
        for i in 0..tmm.n_features() {
            for j in 0..tmm.n_samples() {
                assert!(tmm.get(i, j).is_finite());
            }
        }
    }

    #[test]
    fn test_tmm_config_validation() {
        let counts = create_test_counts();

        // Invalid trim_m
        let config = TmmConfig {
            trim_m: 0.6,
            ..Default::default()
        };
        assert!(norm_tmm_with_config(&counts, &config).is_err());

        // Invalid trim_a
        let config = TmmConfig {
            trim_a: -0.1,
            ..Default::default()
        };
        assert!(norm_tmm_with_config(&counts, &config).is_err());
    }

    #[test]
    fn test_tmm_single_sample_error() {
        let mut tri_mat = TriMat::new((2, 1));
        tri_mat.add_triplet(0, 0, 100);
        tri_mat.add_triplet(1, 0, 100);

        let counts = CountMatrix::new(
            tri_mat.to_csr(),
            vec!["A".into(), "B".into()],
            vec!["S1".into()],
        )
        .unwrap();

        assert!(norm_tmm(&counts).is_err());
    }

    #[test]
    fn test_tmm_factors_function() {
        let counts = create_test_counts();
        let factors = tmm_factors(&counts).unwrap();

        assert_eq!(factors.len(), 3);
        for factor in factors {
            assert!(factor > 0.0);
        }
    }

    #[test]
    fn test_tmm_to_transformed() {
        let counts = create_test_counts();
        let tmm = norm_tmm(&counts).unwrap();
        let transformed = tmm.to_transformed();

        assert_eq!(transformed.n_features(), 4);
        assert_eq!(transformed.n_samples(), 3);
        assert!(transformed.transformation.contains("TMM"));
    }

    #[test]
    fn test_select_reference_sample() {
        let counts = create_test_counts();
        let lib_sizes = counts.col_sums();
        let ref_idx = select_reference_sample(&counts, &lib_sizes);

        // Should select a valid sample
        assert!(ref_idx < 3);
    }
}
