//! Prevalence-based filtering for count matrices.

use crate::data::{CountMatrix, Metadata, Variable};
use crate::error::{DaaError, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Logic for combining group-wise prevalence filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GroupwiseLogic {
    /// Feature must pass threshold in ANY group.
    Any,
    /// Feature must pass threshold in ALL groups.
    All,
    /// Feature must pass threshold in at least N groups.
    AtLeast(usize),
}

/// Filter features by overall prevalence threshold.
///
/// Keeps features that are present in at least `threshold` proportion of samples.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `threshold` - Minimum prevalence (0.0 to 1.0)
///
/// # Returns
/// A new CountMatrix containing only features meeting the threshold.
pub fn filter_prevalence_overall(counts: &CountMatrix, threshold: f64) -> Result<CountMatrix> {
    if !(0.0..=1.0).contains(&threshold) {
        return Err(DaaError::InvalidParameter(
            "Prevalence threshold must be between 0 and 1".to_string(),
        ));
    }

    let n_samples = counts.n_samples();
    let min_samples = (threshold * n_samples as f64).ceil() as usize;

    // Calculate prevalence for each feature in parallel
    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter_map(|row| {
            let nnz = counts.data().outer_view(row)
                .map(|v| v.nnz())
                .unwrap_or(0);
            if nnz >= min_samples {
                Some(row)
            } else {
                None
            }
        })
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No features pass prevalence threshold of {:.1}%",
            threshold * 100.0
        )));
    }

    counts.subset_features(&keep_indices)
}

/// Filter features by group-wise prevalence threshold.
///
/// Evaluates prevalence within each group defined by a categorical variable,
/// then combines results using the specified logic.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `metadata` - Sample metadata
/// * `group_column` - Name of the categorical column defining groups
/// * `threshold` - Minimum prevalence within group (0.0 to 1.0)
/// * `logic` - How to combine group-wise results
///
/// # Returns
/// A new CountMatrix containing only features meeting the criteria.
pub fn filter_prevalence_groupwise(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    threshold: f64,
    logic: GroupwiseLogic,
) -> Result<CountMatrix> {
    if !(0.0..=1.0).contains(&threshold) {
        return Err(DaaError::InvalidParameter(
            "Prevalence threshold must be between 0 and 1".to_string(),
        ));
    }

    // Build group membership
    let groups = build_group_indices(counts, metadata, group_column)?;
    let n_groups = groups.len();

    if n_groups == 0 {
        return Err(DaaError::EmptyData(
            "No groups found in metadata".to_string(),
        ));
    }

    // For AtLeast logic, validate n
    if let GroupwiseLogic::AtLeast(n) = logic {
        if n > n_groups {
            return Err(DaaError::InvalidParameter(format!(
                "AtLeast({}) requires at least {} groups, but only {} found",
                n, n, n_groups
            )));
        }
    }

    // Calculate group-wise prevalence for each feature
    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter(|&row| {
            let row_data = counts.row_dense(row);
            let groups_passing = groups
                .iter()
                .filter(|(_, indices)| {
                    let group_size = indices.len();
                    if group_size == 0 {
                        return false;
                    }
                    let min_samples = (threshold * group_size as f64).ceil() as usize;
                    let nnz = indices.iter().filter(|&&i| row_data[i] > 0).count();
                    nnz >= min_samples
                })
                .count();

            match logic {
                GroupwiseLogic::Any => groups_passing > 0,
                GroupwiseLogic::All => groups_passing == n_groups,
                GroupwiseLogic::AtLeast(n) => groups_passing >= n,
            }
        })
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(format!(
            "No features pass group-wise prevalence threshold of {:.1}%",
            threshold * 100.0
        )));
    }

    counts.subset_features(&keep_indices)
}

/// Build a mapping from group levels to sample indices.
fn build_group_indices(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
) -> Result<HashMap<String, Vec<usize>>> {
    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();

    for (idx, sample_id) in counts.sample_ids().iter().enumerate() {
        let var = metadata
            .get(sample_id, group_column)
            .ok_or_else(|| DaaError::SampleMismatch(format!(
                "Sample '{}' not found in metadata",
                sample_id
            )))?;

        if let Variable::Categorical(level) = var {
            groups.entry(level.clone()).or_default().push(idx);
        } else if !var.is_missing() {
            return Err(DaaError::InvalidVariableType {
                column: group_column.to_string(),
                reason: "Expected categorical variable for grouping".to_string(),
            });
        }
    }

    Ok(groups)
}

/// Result of prevalence filtering with statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilterResult {
    /// Number of features before filtering.
    pub n_before: usize,
    /// Number of features after filtering.
    pub n_after: usize,
    /// Number of features removed.
    pub n_removed: usize,
    /// Proportion of features retained.
    pub retention_rate: f64,
}

impl std::fmt::Display for FilterResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Filter Result")?;
        writeln!(f, "  Before:    {} features", self.n_before)?;
        writeln!(f, "  After:     {} features", self.n_after)?;
        writeln!(f, "  Removed:   {} features", self.n_removed)?;
        writeln!(f, "  Retained:  {:.1}%", self.retention_rate * 100.0)?;
        Ok(())
    }
}

/// Filter with statistics about what was filtered.
pub fn filter_prevalence_overall_with_stats(
    counts: &CountMatrix,
    threshold: f64,
) -> Result<(CountMatrix, FilterResult)> {
    let n_before = counts.n_features();
    let filtered = filter_prevalence_overall(counts, threshold)?;
    let n_after = filtered.n_features();

    let result = FilterResult {
        n_before,
        n_after,
        n_removed: n_before - n_after,
        retention_rate: n_after as f64 / n_before as f64,
    };

    Ok((filtered, result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_matrix() -> CountMatrix {
        // 5 features × 8 samples
        let mut tri_mat = TriMat::new((5, 8));
        // Feature 0: present in all 8 samples (100%)
        for col in 0..8 {
            tri_mat.add_triplet(0, col, 10);
        }
        // Feature 1: present in 6 samples (75%)
        for col in 0..6 {
            tri_mat.add_triplet(1, col, 20);
        }
        // Feature 2: present in 4 samples (50%)
        for col in 0..4 {
            tri_mat.add_triplet(2, col, 30);
        }
        // Feature 3: present in 2 samples (25%)
        for col in 0..2 {
            tri_mat.add_triplet(3, col, 40);
        }
        // Feature 4: present in 1 sample (12.5%)
        tri_mat.add_triplet(4, 0, 50);

        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..8).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        // Group A: S0-S3 (4 samples)
        for i in 0..4 {
            writeln!(file, "S{}\tA", i).unwrap();
        }
        // Group B: S4-S7 (4 samples)
        for i in 4..8 {
            writeln!(file, "S{}\tB", i).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_filter_prevalence_overall() {
        let counts = create_test_matrix();

        // 50% threshold should keep features 0, 1, 2
        let filtered = filter_prevalence_overall(&counts, 0.5).unwrap();
        assert_eq!(filtered.n_features(), 3);
        assert_eq!(filtered.feature_ids(), &["feat_0", "feat_1", "feat_2"]);

        // 75% threshold should keep features 0, 1
        let filtered = filter_prevalence_overall(&counts, 0.75).unwrap();
        assert_eq!(filtered.n_features(), 2);
        assert_eq!(filtered.feature_ids(), &["feat_0", "feat_1"]);
    }

    #[test]
    fn test_filter_prevalence_groupwise_any() {
        let counts = create_test_matrix();
        let metadata = create_test_metadata();

        // Feature 2 is present in 4 samples (all in group A, 100% in A, 0% in B)
        // Feature 3 is present in 2 samples (both in group A, 50% in A, 0% in B)
        // With 50% threshold and Any logic, features meeting threshold in either group pass
        let filtered = filter_prevalence_groupwise(
            &counts,
            &metadata,
            "group",
            0.5,
            GroupwiseLogic::Any,
        ).unwrap();

        // Features 0, 1 pass in both groups; 2, 3 pass in group A
        assert_eq!(filtered.n_features(), 4);
    }

    #[test]
    fn test_filter_prevalence_groupwise_all() {
        let counts = create_test_matrix();
        let metadata = create_test_metadata();

        // With 50% threshold and All logic, only features meeting threshold in both groups pass
        let filtered = filter_prevalence_groupwise(
            &counts,
            &metadata,
            "group",
            0.5,
            GroupwiseLogic::All,
        ).unwrap();

        // Only features 0 and 1 (75%+ overall) pass in both groups
        assert_eq!(filtered.n_features(), 2);
    }

    #[test]
    fn test_filter_with_stats() {
        let counts = create_test_matrix();
        let (filtered, stats) = filter_prevalence_overall_with_stats(&counts, 0.5).unwrap();

        assert_eq!(stats.n_before, 5);
        assert_eq!(stats.n_after, 3);
        assert_eq!(stats.n_removed, 2);
        assert!((stats.retention_rate - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_threshold() {
        let counts = create_test_matrix();
        assert!(filter_prevalence_overall(&counts, -0.1).is_err());
        assert!(filter_prevalence_overall(&counts, 1.1).is_err());
    }
}
