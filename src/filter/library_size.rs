//! Library size-based filtering for samples.

use crate::data::{CountMatrix, Metadata};
use crate::error::{DaaError, Result};
use serde::{Deserialize, Serialize};

/// Filter samples by library size (total counts).
///
/// Removes samples with total counts below or above specified thresholds.
/// This is important for quality control as very low library sizes indicate
/// poor sequencing depth.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `min_reads` - Minimum total reads per sample (None for no minimum)
/// * `max_reads` - Maximum total reads per sample (None for no maximum)
///
/// # Returns
/// A new CountMatrix containing only samples meeting the criteria.
pub fn filter_library_size(
    counts: &CountMatrix,
    min_reads: Option<u64>,
    max_reads: Option<u64>,
) -> Result<CountMatrix> {
    if let (Some(min), Some(max)) = (min_reads, max_reads) {
        if max < min {
            return Err(DaaError::InvalidParameter(
                "max_reads cannot be less than min_reads".to_string(),
            ));
        }
    }

    let col_sums = counts.col_sums();
    let min = min_reads.unwrap_or(0);
    let max = max_reads.unwrap_or(u64::MAX);

    let keep_indices: Vec<usize> = col_sums
        .iter()
        .enumerate()
        .filter(|(_, &sum)| sum >= min && sum <= max)
        .map(|(i, _)| i)
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No samples have library size between {} and {}",
            min, max
        )));
    }

    counts.subset_samples(&keep_indices)
}

/// Filter samples by library size using quantile thresholds.
///
/// Removes samples in the lower or upper quantiles of library size distribution.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `lower_quantile` - Remove samples below this quantile (0.0 to 1.0)
/// * `upper_quantile` - Remove samples above this quantile (0.0 to 1.0)
///
/// # Returns
/// A new CountMatrix containing only samples within the quantile range.
pub fn filter_library_size_quantile(
    counts: &CountMatrix,
    lower_quantile: f64,
    upper_quantile: f64,
) -> Result<CountMatrix> {
    if !(0.0..=1.0).contains(&lower_quantile) || !(0.0..=1.0).contains(&upper_quantile) {
        return Err(DaaError::InvalidParameter(
            "Quantiles must be between 0 and 1".to_string(),
        ));
    }
    if upper_quantile < lower_quantile {
        return Err(DaaError::InvalidParameter(
            "upper_quantile cannot be less than lower_quantile".to_string(),
        ));
    }

    let col_sums = counts.col_sums();
    let mut sorted_sums: Vec<u64> = col_sums.clone();
    sorted_sums.sort_unstable();

    let n = sorted_sums.len();
    let lower_idx = ((lower_quantile * n as f64).floor() as usize).min(n - 1);
    let upper_idx = ((upper_quantile * n as f64).ceil() as usize).min(n - 1);

    let min_threshold = sorted_sums[lower_idx];
    let max_threshold = sorted_sums[upper_idx];

    filter_library_size(counts, Some(min_threshold), Some(max_threshold))
}

/// Filter samples by minimum library size relative to median.
///
/// A common approach is to remove samples with library size below
/// some fraction of the median library size.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `min_fraction_of_median` - Minimum library size as fraction of median
///
/// # Returns
/// A new CountMatrix containing only samples meeting the threshold.
pub fn filter_library_size_relative(
    counts: &CountMatrix,
    min_fraction_of_median: f64,
) -> Result<CountMatrix> {
    if min_fraction_of_median < 0.0 {
        return Err(DaaError::InvalidParameter(
            "min_fraction_of_median must be non-negative".to_string(),
        ));
    }

    let col_sums = counts.col_sums();
    let mut sorted_sums: Vec<u64> = col_sums.clone();
    sorted_sums.sort_unstable();

    let n = sorted_sums.len();
    let median = if n % 2 == 0 {
        (sorted_sums[n / 2 - 1] + sorted_sums[n / 2]) / 2
    } else {
        sorted_sums[n / 2]
    };

    let min_threshold = (median as f64 * min_fraction_of_median).round() as u64;

    filter_library_size(counts, Some(min_threshold), None)
}

/// Result of library size filtering with statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibrarySizeFilterResult {
    /// Number of samples before filtering.
    pub n_before: usize,
    /// Number of samples after filtering.
    pub n_after: usize,
    /// Number of samples removed.
    pub n_removed: usize,
    /// IDs of removed samples.
    pub removed_samples: Vec<String>,
    /// Minimum library size in retained samples.
    pub min_retained: u64,
    /// Maximum library size in retained samples.
    pub max_retained: u64,
    /// Median library size in retained samples.
    pub median_retained: u64,
}

impl std::fmt::Display for LibrarySizeFilterResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Library Size Filter Result")?;
        writeln!(f, "  Samples before:  {}", self.n_before)?;
        writeln!(f, "  Samples after:   {}", self.n_after)?;
        writeln!(f, "  Samples removed: {}", self.n_removed)?;
        if !self.removed_samples.is_empty() {
            writeln!(f, "  Removed: {:?}", self.removed_samples)?;
        }
        writeln!(f, "  Retained range:  {} - {}", self.min_retained, self.max_retained)?;
        writeln!(f, "  Retained median: {}", self.median_retained)?;
        Ok(())
    }
}

/// Filter with statistics about what was filtered.
pub fn filter_library_size_with_stats(
    counts: &CountMatrix,
    min_reads: Option<u64>,
    max_reads: Option<u64>,
) -> Result<(CountMatrix, LibrarySizeFilterResult)> {
    let n_before = counts.n_samples();
    let col_sums = counts.col_sums();
    let min = min_reads.unwrap_or(0);
    let max = max_reads.unwrap_or(u64::MAX);

    // Find which samples to remove
    let removed_samples: Vec<String> = counts
        .sample_ids()
        .iter()
        .enumerate()
        .filter(|(i, _)| {
            let sum = col_sums[*i];
            sum < min || sum > max
        })
        .map(|(_, id)| id.clone())
        .collect();

    let filtered = filter_library_size(counts, min_reads, max_reads)?;

    let n_after = filtered.n_samples();
    let filtered_sums = filtered.col_sums();
    let mut sorted_sums: Vec<u64> = filtered_sums.clone();
    sorted_sums.sort_unstable();

    let n = sorted_sums.len();
    let median = if n % 2 == 0 {
        (sorted_sums[n / 2 - 1] + sorted_sums[n / 2]) / 2
    } else {
        sorted_sums[n / 2]
    };

    let result = LibrarySizeFilterResult {
        n_before,
        n_after,
        n_removed: n_before - n_after,
        removed_samples,
        min_retained: *sorted_sums.first().unwrap_or(&0),
        max_retained: *sorted_sums.last().unwrap_or(&0),
        median_retained: median,
    };

    Ok((filtered, result))
}

/// Filter both counts and metadata to keep only retained samples.
///
/// This ensures counts and metadata stay synchronized after sample filtering.
pub fn filter_library_size_with_metadata(
    counts: &CountMatrix,
    metadata: &Metadata,
    min_reads: Option<u64>,
    max_reads: Option<u64>,
) -> Result<(CountMatrix, Metadata)> {
    let filtered_counts = filter_library_size(counts, min_reads, max_reads)?;
    let filtered_metadata = metadata.subset_samples(filtered_counts.sample_ids())?;
    Ok((filtered_counts, filtered_metadata))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        // 3 features × 6 samples with varying library sizes
        let mut tri_mat = TriMat::new((3, 6));

        // Sample library sizes: 100, 200, 500, 1000, 2000, 5000
        let lib_sizes = [100u64, 200, 500, 1000, 2000, 5000];

        for (col, &lib_size) in lib_sizes.iter().enumerate() {
            // Distribute counts across features
            tri_mat.add_triplet(0, col, lib_size / 2);
            tri_mat.add_triplet(1, col, lib_size / 3);
            tri_mat.add_triplet(2, col, lib_size - lib_size / 2 - lib_size / 3);
        }

        let feature_ids: Vec<String> = (0..3).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..6).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_filter_library_size_min() {
        let counts = create_test_matrix();

        // Filter to keep only samples with >= 500 reads
        let filtered = filter_library_size(&counts, Some(500), None).unwrap();

        // Should keep S2 (500), S3 (1000), S4 (2000), S5 (5000)
        assert_eq!(filtered.n_samples(), 4);
        assert_eq!(filtered.sample_ids(), &["S2", "S3", "S4", "S5"]);
    }

    #[test]
    fn test_filter_library_size_range() {
        let counts = create_test_matrix();

        // Filter to keep samples with 200-2000 reads
        let filtered = filter_library_size(&counts, Some(200), Some(2000)).unwrap();

        // Should keep S1 (200), S2 (500), S3 (1000), S4 (2000)
        assert_eq!(filtered.n_samples(), 4);
        assert_eq!(filtered.sample_ids(), &["S1", "S2", "S3", "S4"]);
    }

    #[test]
    fn test_filter_library_size_quantile() {
        let counts = create_test_matrix();

        // Filter to keep middle 50% (remove bottom 25% and top 25%)
        // With 6 samples, 25th percentile is index 1, 75th is index 4-5
        let filtered = filter_library_size_quantile(&counts, 0.25, 0.75).unwrap();

        // Should keep samples between the 25th and 75th percentile thresholds
        assert!(filtered.n_samples() >= 2);
        assert!(filtered.n_samples() <= 6);
    }

    #[test]
    fn test_filter_library_size_relative() {
        let counts = create_test_matrix();

        // Median library size is between 500 and 1000 (~750)
        // Filter to keep samples with >= 50% of median
        let filtered = filter_library_size_relative(&counts, 0.5).unwrap();

        // Should keep most samples except the very small ones
        assert!(filtered.n_samples() >= 3);
    }

    #[test]
    fn test_filter_library_size_with_stats() {
        let counts = create_test_matrix();

        let (filtered, stats) = filter_library_size_with_stats(&counts, Some(500), None).unwrap();

        assert_eq!(stats.n_before, 6);
        assert_eq!(stats.n_after, 4);
        assert_eq!(stats.n_removed, 2);
        assert_eq!(stats.removed_samples, vec!["S0", "S1"]);
    }

    #[test]
    fn test_invalid_parameters() {
        let counts = create_test_matrix();

        assert!(filter_library_size(&counts, Some(1000), Some(500)).is_err()); // max < min
        assert!(filter_library_size_quantile(&counts, 0.8, 0.2).is_err()); // upper < lower
        assert!(filter_library_size_quantile(&counts, -0.1, 0.9).is_err()); // invalid quantile
        assert!(filter_library_size_relative(&counts, -0.5).is_err()); // negative fraction
    }
}
