//! Abundance-based filtering for count matrices.

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Filter features by relative abundance thresholds.
///
/// Keeps features whose mean relative abundance falls within the specified range.
/// Relative abundance is calculated as the proportion of total counts.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `min_abundance` - Minimum mean relative abundance (0.0 to 1.0)
/// * `max_abundance` - Maximum mean relative abundance (0.0 to 1.0), or None for no upper limit
///
/// # Returns
/// A new CountMatrix containing only features meeting the criteria.
pub fn filter_abundance(
    counts: &CountMatrix,
    min_abundance: f64,
    max_abundance: Option<f64>,
) -> Result<CountMatrix> {
    if !(0.0..=1.0).contains(&min_abundance) {
        return Err(DaaError::InvalidParameter(
            "min_abundance must be between 0 and 1".to_string(),
        ));
    }
    if let Some(max) = max_abundance {
        if !(0.0..=1.0).contains(&max) {
            return Err(DaaError::InvalidParameter(
                "max_abundance must be between 0 and 1".to_string(),
            ));
        }
        if max < min_abundance {
            return Err(DaaError::InvalidParameter(
                "max_abundance cannot be less than min_abundance".to_string(),
            ));
        }
    }

    // Calculate total counts per sample (library sizes)
    let col_sums = counts.col_sums();
    let total_reads: u64 = col_sums.iter().sum();

    if total_reads == 0 {
        return Err(DaaError::EmptyData("All counts are zero".to_string()));
    }

    // Calculate mean relative abundance for each feature
    let row_sums = counts.row_sums();
    let max_abund = max_abundance.unwrap_or(1.0);

    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter(|&row| {
            let rel_abund = row_sums[row] as f64 / total_reads as f64;
            rel_abund >= min_abundance && rel_abund <= max_abund
        })
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No features have relative abundance between {:.2e} and {:.2e}",
            min_abundance, max_abund
        )));
    }

    counts.subset_features(&keep_indices)
}

/// Filter features by mean abundance per sample.
///
/// Keeps features whose mean count per sample (where present) meets the threshold.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `min_mean` - Minimum mean count when present
///
/// # Returns
/// A new CountMatrix containing only features meeting the threshold.
pub fn filter_mean_abundance(counts: &CountMatrix, min_mean: f64) -> Result<CountMatrix> {
    if min_mean < 0.0 {
        return Err(DaaError::InvalidParameter(
            "min_mean must be non-negative".to_string(),
        ));
    }

    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter(|&row| {
            let row_data = counts.row_dense(row);
            let non_zero: Vec<u64> = row_data.iter().copied().filter(|&v| v > 0).collect();
            if non_zero.is_empty() {
                return false;
            }
            let mean = non_zero.iter().sum::<u64>() as f64 / non_zero.len() as f64;
            mean >= min_mean
        })
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No features have mean abundance >= {}",
            min_mean
        )));
    }

    counts.subset_features(&keep_indices)
}

/// Filter features by total count threshold.
///
/// Keeps features whose total count across all samples meets the threshold.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `min_total` - Minimum total count
///
/// # Returns
/// A new CountMatrix containing only features meeting the threshold.
pub fn filter_total_count(counts: &CountMatrix, min_total: u64) -> Result<CountMatrix> {
    let row_sums = counts.row_sums();

    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter(|&row| row_sums[row] >= min_total)
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No features have total count >= {}",
            min_total
        )));
    }

    counts.subset_features(&keep_indices)
}

/// Result of abundance filtering with statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AbundanceFilterResult {
    /// Number of features before filtering.
    pub n_before: usize,
    /// Number of features after filtering.
    pub n_after: usize,
    /// Number of features removed.
    pub n_removed: usize,
    /// Proportion of features retained.
    pub retention_rate: f64,
    /// Proportion of total reads retained.
    pub reads_retained: f64,
}

impl std::fmt::Display for AbundanceFilterResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Abundance Filter Result")?;
        writeln!(f, "  Features before:  {}", self.n_before)?;
        writeln!(f, "  Features after:   {}", self.n_after)?;
        writeln!(f, "  Features removed: {}", self.n_removed)?;
        writeln!(f, "  Feature retention: {:.1}%", self.retention_rate * 100.0)?;
        writeln!(f, "  Reads retained:    {:.1}%", self.reads_retained * 100.0)?;
        Ok(())
    }
}

/// Filter with statistics about what was filtered.
pub fn filter_abundance_with_stats(
    counts: &CountMatrix,
    min_abundance: f64,
    max_abundance: Option<f64>,
) -> Result<(CountMatrix, AbundanceFilterResult)> {
    let n_before = counts.n_features();
    let total_reads_before: u64 = counts.row_sums().iter().sum();

    let filtered = filter_abundance(counts, min_abundance, max_abundance)?;

    let n_after = filtered.n_features();
    let total_reads_after: u64 = filtered.row_sums().iter().sum();

    let result = AbundanceFilterResult {
        n_before,
        n_after,
        n_removed: n_before - n_after,
        retention_rate: n_after as f64 / n_before as f64,
        reads_retained: if total_reads_before > 0 {
            total_reads_after as f64 / total_reads_before as f64
        } else {
            0.0
        },
    };

    Ok((filtered, result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        // 5 features × 4 samples
        // Total counts: 1000
        let mut tri_mat = TriMat::new((5, 4));

        // Feature 0: high abundance (400 total = 40%)
        for col in 0..4 {
            tri_mat.add_triplet(0, col, 100);
        }
        // Feature 1: medium abundance (300 total = 30%)
        for col in 0..4 {
            tri_mat.add_triplet(1, col, 75);
        }
        // Feature 2: low abundance (200 total = 20%)
        for col in 0..4 {
            tri_mat.add_triplet(2, col, 50);
        }
        // Feature 3: very low abundance (90 total = 9%)
        for col in 0..3 {
            tri_mat.add_triplet(3, col, 30);
        }
        // Feature 4: rare (10 total = 1%)
        tri_mat.add_triplet(4, 0, 10);

        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..4).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_filter_abundance_min() {
        let counts = create_test_matrix();

        // Filter to keep only features with >= 10% relative abundance
        let filtered = filter_abundance(&counts, 0.10, None).unwrap();

        // Should keep features 0 (40%), 1 (30%), 2 (20%)
        assert_eq!(filtered.n_features(), 3);
        assert_eq!(filtered.feature_ids(), &["feat_0", "feat_1", "feat_2"]);
    }

    #[test]
    fn test_filter_abundance_range() {
        let counts = create_test_matrix();

        // Filter to keep only features between 15% and 35%
        let filtered = filter_abundance(&counts, 0.15, Some(0.35)).unwrap();

        // Should keep features 1 (30%), 2 (20%)
        assert_eq!(filtered.n_features(), 2);
        assert_eq!(filtered.feature_ids(), &["feat_1", "feat_2"]);
    }

    #[test]
    fn test_filter_mean_abundance() {
        let counts = create_test_matrix();

        // Filter to keep features with mean count >= 50
        let filtered = filter_mean_abundance(&counts, 50.0).unwrap();

        // Feature 0: mean 100, Feature 1: mean 75, Feature 2: mean 50
        assert_eq!(filtered.n_features(), 3);
    }

    #[test]
    fn test_filter_total_count() {
        let counts = create_test_matrix();

        // Filter to keep features with total count >= 200
        let filtered = filter_total_count(&counts, 200).unwrap();

        // Features 0 (400), 1 (300), 2 (200)
        assert_eq!(filtered.n_features(), 3);
    }

    #[test]
    fn test_filter_abundance_with_stats() {
        let counts = create_test_matrix();

        let (filtered, stats) = filter_abundance_with_stats(&counts, 0.10, None).unwrap();

        assert_eq!(stats.n_before, 5);
        assert_eq!(stats.n_after, 3);
        assert_eq!(stats.n_removed, 2);
        assert!((stats.retention_rate - 0.6).abs() < 1e-10);
        // Reads retained: (400+300+200)/1000 = 90%
        assert!((stats.reads_retained - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_parameters() {
        let counts = create_test_matrix();

        assert!(filter_abundance(&counts, -0.1, None).is_err());
        assert!(filter_abundance(&counts, 1.1, None).is_err());
        assert!(filter_abundance(&counts, 0.5, Some(0.3)).is_err()); // max < min
        assert!(filter_mean_abundance(&counts, -1.0).is_err());
    }
}
