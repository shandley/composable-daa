//! Sparsity profiling for count matrices.

use crate::data::CountMatrix;
use serde::{Deserialize, Serialize};

/// Profile of sparsity characteristics in a count matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparsityProfile {
    /// Total number of entries (features × samples).
    pub total_entries: usize,
    /// Number of non-zero entries.
    pub nonzero_entries: usize,
    /// Number of zero entries.
    pub zero_entries: usize,
    /// Overall sparsity (proportion of zeros).
    pub sparsity: f64,
    /// Sparsity per feature (row).
    pub feature_sparsity: Vec<f64>,
    /// Sparsity per sample (column).
    pub sample_sparsity: Vec<f64>,
    /// Mean sparsity across features.
    pub mean_feature_sparsity: f64,
    /// Mean sparsity across samples.
    pub mean_sample_sparsity: f64,
    /// Median sparsity across features.
    pub median_feature_sparsity: f64,
    /// Median sparsity across samples.
    pub median_sample_sparsity: f64,
}

impl SparsityProfile {
    /// Check if the data is highly sparse (> 50% zeros).
    pub fn is_highly_sparse(&self) -> bool {
        self.sparsity > 0.5
    }

    /// Check if the data is ultra-sparse (> 90% zeros).
    pub fn is_ultra_sparse(&self) -> bool {
        self.sparsity > 0.9
    }
}

impl std::fmt::Display for SparsityProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Sparsity Profile")?;
        writeln!(f, "  Total entries:     {}", self.total_entries)?;
        writeln!(f, "  Non-zero entries:  {}", self.nonzero_entries)?;
        writeln!(f, "  Zero entries:      {}", self.zero_entries)?;
        writeln!(f, "  Overall sparsity:  {:.2}%", self.sparsity * 100.0)?;
        writeln!(f, "  Mean feature sparsity:   {:.2}%", self.mean_feature_sparsity * 100.0)?;
        writeln!(f, "  Median feature sparsity: {:.2}%", self.median_feature_sparsity * 100.0)?;
        writeln!(f, "  Mean sample sparsity:    {:.2}%", self.mean_sample_sparsity * 100.0)?;
        writeln!(f, "  Median sample sparsity:  {:.2}%", self.median_sample_sparsity * 100.0)?;
        Ok(())
    }
}

/// Profile sparsity characteristics of a count matrix.
pub fn profile_sparsity(counts: &CountMatrix) -> SparsityProfile {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let total_entries = n_features * n_samples;
    let nonzero_entries = counts.nnz();
    let zero_entries = total_entries - nonzero_entries;
    let sparsity = zero_entries as f64 / total_entries as f64;

    // Per-feature sparsity
    let feature_sparsity: Vec<f64> = (0..n_features)
        .map(|row| {
            let row_nnz = counts.data().outer_view(row)
                .map(|v| v.nnz())
                .unwrap_or(0);
            (n_samples - row_nnz) as f64 / n_samples as f64
        })
        .collect();

    // Per-sample sparsity
    let mut sample_nnz = vec![0usize; n_samples];
    for row_vec in counts.data().outer_iterator() {
        for (col, _) in row_vec.iter() {
            sample_nnz[col] += 1;
        }
    }
    let sample_sparsity: Vec<f64> = sample_nnz
        .iter()
        .map(|&nnz| (n_features - nnz) as f64 / n_features as f64)
        .collect();

    // Statistics
    let mean_feature_sparsity = feature_sparsity.iter().sum::<f64>() / n_features as f64;
    let mean_sample_sparsity = sample_sparsity.iter().sum::<f64>() / n_samples as f64;

    let median_feature_sparsity = median(&feature_sparsity);
    let median_sample_sparsity = median(&sample_sparsity);

    SparsityProfile {
        total_entries,
        nonzero_entries,
        zero_entries,
        sparsity,
        feature_sparsity,
        sample_sparsity,
        mean_feature_sparsity,
        mean_sample_sparsity,
        median_feature_sparsity,
        median_sample_sparsity,
    }
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
    } else {
        sorted[n / 2]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        // 3 features × 4 samples, some zeros
        let mut tri_mat = TriMat::new((3, 4));
        tri_mat.add_triplet(0, 0, 10);
        tri_mat.add_triplet(0, 1, 20);
        // (0,2) and (0,3) are zero
        tri_mat.add_triplet(0, 3, 5);
        tri_mat.add_triplet(1, 0, 100);
        tri_mat.add_triplet(1, 1, 200);
        tri_mat.add_triplet(1, 2, 150);
        tri_mat.add_triplet(1, 3, 175);
        tri_mat.add_triplet(2, 0, 1);
        // feature 2 is very sparse

        let feature_ids = vec!["A".into(), "B".into(), "C".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_profile_sparsity() {
        let counts = create_test_matrix();
        let profile = profile_sparsity(&counts);

        assert_eq!(profile.total_entries, 12);
        assert_eq!(profile.nonzero_entries, 8);
        assert_eq!(profile.zero_entries, 4);
        assert!((profile.sparsity - 4.0/12.0).abs() < 1e-10);
    }

    #[test]
    fn test_feature_sparsity() {
        let counts = create_test_matrix();
        let profile = profile_sparsity(&counts);

        // Feature 0: 3 non-zero, 1 zero → 25% sparse
        assert!((profile.feature_sparsity[0] - 0.25).abs() < 1e-10);
        // Feature 1: 4 non-zero, 0 zero → 0% sparse
        assert!((profile.feature_sparsity[1] - 0.0).abs() < 1e-10);
        // Feature 2: 1 non-zero, 3 zero → 75% sparse
        assert!((profile.feature_sparsity[2] - 0.75).abs() < 1e-10);
    }
}
