//! Total Sum Scaling (TSS) normalization for compositional data.
//!
//! TSS converts counts to relative abundances by dividing each count by the
//! total counts in that sample. This is the simplest normalization for
//! compositional data and is widely used in microbiome analysis.

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Result of TSS normalization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssMatrix {
    /// The normalized data (features × samples).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers.
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// Scale factor applied (1.0 for proportions, 1e6 for CPM, etc.).
    pub scale_factor: f64,
    /// Library sizes (total counts per sample before normalization).
    pub library_sizes: Vec<u64>,
}

impl TssMatrix {
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
            transformation: format!("TSS (scale={})", self.scale_factor),
            geometric_means: vec![], // Not applicable to TSS
        }
    }
}

/// Apply Total Sum Scaling normalization to a count matrix.
///
/// TSS divides each count by the total counts in that sample, converting
/// to relative abundances (proportions). Optionally, a scale factor can
/// be applied (e.g., 1e6 for "counts per million" or CPM).
///
/// # Formula
/// For sample j: TSS(x_ij) = x_ij / sum(x_j) * scale_factor
///
/// # Arguments
/// * `counts` - Count matrix
/// * `scale_factor` - Multiplier for normalized values (1.0 for proportions,
///   1e6 for CPM, 1e4 for counts per 10k, etc.)
///
/// # Returns
/// A TssMatrix containing normalized values.
///
/// # Example
/// ```ignore
/// // Normalize to proportions (sum to 1)
/// let tss = norm_tss(&counts, 1.0)?;
///
/// // Normalize to counts per million (CPM)
/// let cpm = norm_tss(&counts, 1e6)?;
/// ```
pub fn norm_tss(counts: &CountMatrix, scale_factor: f64) -> Result<TssMatrix> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply TSS to empty matrix".to_string(),
        ));
    }

    if scale_factor <= 0.0 {
        return Err(DaaError::InvalidParameter(
            "Scale factor must be positive".to_string(),
        ));
    }

    // Get library sizes (column sums)
    let library_sizes = counts.col_sums();

    // Check for zero library sizes
    for (j, &lib_size) in library_sizes.iter().enumerate() {
        if lib_size == 0 {
            return Err(DaaError::Numerical(format!(
                "Sample {} has zero total counts, cannot normalize",
                counts.sample_ids()[j]
            )));
        }
    }

    // Apply TSS normalization in parallel
    let normalized_cols: Vec<Vec<f64>> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            let lib_size = library_sizes[j] as f64;
            (0..n_features)
                .map(|i| (counts.get(i, j) as f64 / lib_size) * scale_factor)
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

    Ok(TssMatrix {
        data,
        feature_ids: counts.feature_ids().to_vec(),
        sample_ids: counts.sample_ids().to_vec(),
        scale_factor,
        library_sizes,
    })
}

/// Apply TSS normalization with pseudocount added first.
///
/// This is useful when downstream analysis requires non-zero values
/// (e.g., log transformation).
///
/// # Arguments
/// * `counts` - Count matrix
/// * `pseudocount` - Value to add before normalization
/// * `scale_factor` - Multiplier for normalized values
pub fn norm_tss_with_pseudocount(
    counts: &CountMatrix,
    pseudocount: f64,
    scale_factor: f64,
) -> Result<TssMatrix> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply TSS to empty matrix".to_string(),
        ));
    }

    if scale_factor <= 0.0 {
        return Err(DaaError::InvalidParameter(
            "Scale factor must be positive".to_string(),
        ));
    }

    if pseudocount < 0.0 {
        return Err(DaaError::InvalidParameter(
            "Pseudocount must be non-negative".to_string(),
        ));
    }

    // Get original library sizes
    let library_sizes = counts.col_sums();

    // Apply TSS with pseudocount
    let normalized_cols: Vec<Vec<f64>> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            // Add pseudocount to each value, then compute new total
            let values_with_pseudo: Vec<f64> = (0..n_features)
                .map(|i| counts.get(i, j) as f64 + pseudocount)
                .collect();
            let total: f64 = values_with_pseudo.iter().sum();

            values_with_pseudo
                .iter()
                .map(|&v| (v / total) * scale_factor)
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

    Ok(TssMatrix {
        data,
        feature_ids: counts.feature_ids().to_vec(),
        sample_ids: counts.sample_ids().to_vec(),
        scale_factor,
        library_sizes,
    })
}

/// Common scale factors for TSS normalization.
pub mod scale {
    /// Proportions (sum to 1.0 per sample).
    pub const PROPORTION: f64 = 1.0;
    /// Counts per million (CPM).
    pub const CPM: f64 = 1_000_000.0;
    /// Counts per 10,000.
    pub const CP10K: f64 = 10_000.0;
    /// Counts per 100 (percentages).
    pub const PERCENT: f64 = 100.0;
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use sprs::TriMat;

    fn create_test_counts() -> CountMatrix {
        // 3 features × 4 samples
        let mut tri_mat = TriMat::new((3, 4));

        // Sample 0: total = 100
        tri_mat.add_triplet(0, 0, 50);  // 50%
        tri_mat.add_triplet(1, 0, 30);  // 30%
        tri_mat.add_triplet(2, 0, 20);  // 20%

        // Sample 1: total = 200
        tri_mat.add_triplet(0, 1, 100); // 50%
        tri_mat.add_triplet(1, 1, 60);  // 30%
        tri_mat.add_triplet(2, 1, 40);  // 20%

        // Sample 2: total = 50
        tri_mat.add_triplet(0, 2, 25);  // 50%
        tri_mat.add_triplet(1, 2, 15);  // 30%
        tri_mat.add_triplet(2, 2, 10);  // 20%

        // Sample 3: total = 1000
        tri_mat.add_triplet(0, 3, 500); // 50%
        tri_mat.add_triplet(1, 3, 300); // 30%
        tri_mat.add_triplet(2, 3, 200); // 20%

        let feature_ids = vec!["A".into(), "B".into(), "C".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_tss_proportions() {
        let counts = create_test_counts();
        let tss = norm_tss(&counts, scale::PROPORTION).unwrap();

        assert_eq!(tss.n_features(), 3);
        assert_eq!(tss.n_samples(), 4);
        assert_eq!(tss.scale_factor, 1.0);

        // All samples should have same proportions (50%, 30%, 20%)
        for j in 0..4 {
            assert_relative_eq!(tss.get(0, j), 0.50, epsilon = 1e-10);
            assert_relative_eq!(tss.get(1, j), 0.30, epsilon = 1e-10);
            assert_relative_eq!(tss.get(2, j), 0.20, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_tss_column_sums() {
        let counts = create_test_counts();
        let tss = norm_tss(&counts, scale::PROPORTION).unwrap();

        // Each column should sum to 1.0
        for j in 0..tss.n_samples() {
            let col_sum: f64 = (0..tss.n_features()).map(|i| tss.get(i, j)).sum();
            assert_relative_eq!(col_sum, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_tss_cpm() {
        let counts = create_test_counts();
        let cpm = norm_tss(&counts, scale::CPM).unwrap();

        // Each column should sum to 1,000,000
        for j in 0..cpm.n_samples() {
            let col_sum: f64 = (0..cpm.n_features()).map(|i| cpm.get(i, j)).sum();
            assert_relative_eq!(col_sum, 1_000_000.0, epsilon = 1e-6);
        }

        // First feature should be 500,000 CPM (50%)
        assert_relative_eq!(cpm.get(0, 0), 500_000.0, epsilon = 1e-6);
    }

    #[test]
    fn test_tss_library_sizes() {
        let counts = create_test_counts();
        let tss = norm_tss(&counts, scale::PROPORTION).unwrap();

        assert_eq!(tss.library_sizes, vec![100, 200, 50, 1000]);
    }

    #[test]
    fn test_tss_with_pseudocount() {
        let counts = create_test_counts();
        let tss = norm_tss_with_pseudocount(&counts, 1.0, scale::PROPORTION).unwrap();

        // With pseudocount=1, sample 0 totals become 103 instead of 100
        // Feature A: (50+1)/103 ≈ 0.495
        assert_relative_eq!(tss.get(0, 0), 51.0 / 103.0, epsilon = 1e-10);

        // Columns should still sum to 1.0
        for j in 0..tss.n_samples() {
            let col_sum: f64 = (0..tss.n_features()).map(|i| tss.get(i, j)).sum();
            assert_relative_eq!(col_sum, 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_tss_to_transformed() {
        let counts = create_test_counts();
        let tss = norm_tss(&counts, scale::CPM).unwrap();
        let transformed = tss.to_transformed();

        assert_eq!(transformed.n_features(), 3);
        assert_eq!(transformed.n_samples(), 4);
        assert!(transformed.transformation.contains("TSS"));
    }

    #[test]
    fn test_tss_empty_matrix() {
        let tri_mat: TriMat<u64> = TriMat::new((0, 0));
        let counts = CountMatrix::new(tri_mat.to_csr(), vec![], vec![]).unwrap();
        let result = norm_tss(&counts, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_tss_invalid_scale_factor() {
        let counts = create_test_counts();
        let result = norm_tss(&counts, 0.0);
        assert!(result.is_err());

        let result = norm_tss(&counts, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_tss_zero_library_size() {
        // Create matrix with one sample having all zeros
        let mut tri_mat = TriMat::new((2, 2));
        tri_mat.add_triplet(0, 0, 10);
        tri_mat.add_triplet(1, 0, 10);
        // Sample 1 has no counts (all zeros)

        let feature_ids = vec!["A".into(), "B".into()];
        let sample_ids = vec!["S1".into(), "S2".into()];
        let counts = CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap();

        let result = norm_tss(&counts, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn test_tss_accessors() {
        let counts = create_test_counts();
        let tss = norm_tss(&counts, scale::PROPORTION).unwrap();

        // Test row accessor
        let row0 = tss.row(0);
        assert_eq!(row0.len(), 4);
        assert!(row0.iter().all(|&v| (v - 0.5).abs() < 1e-10));

        // Test col accessor
        let col0 = tss.col(0);
        assert_eq!(col0.len(), 3);
        assert_relative_eq!(col0[0], 0.50, epsilon = 1e-10);
        assert_relative_eq!(col0[1], 0.30, epsilon = 1e-10);
        assert_relative_eq!(col0[2], 0.20, epsilon = 1e-10);
    }
}
