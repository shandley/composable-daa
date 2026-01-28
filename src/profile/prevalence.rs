//! Prevalence profiling for count matrices.

use crate::data::CountMatrix;
use serde::{Deserialize, Serialize};

/// Profile of prevalence characteristics in a count matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceProfile {
    /// Number of features.
    pub n_features: usize,
    /// Number of samples.
    pub n_samples: usize,
    /// Prevalence (proportion of non-zero samples) per feature.
    pub feature_prevalence: Vec<f64>,
    /// Mean prevalence across features.
    pub mean_prevalence: f64,
    /// Median prevalence across features.
    pub median_prevalence: f64,
    /// Minimum prevalence.
    pub min_prevalence: f64,
    /// Maximum prevalence.
    pub max_prevalence: f64,
    /// Number of features present in all samples.
    pub n_ubiquitous: usize,
    /// Number of features present in only one sample.
    pub n_singletons: usize,
    /// Number of features below 10% prevalence.
    pub n_rare: usize,
    /// Number of features below 5% prevalence.
    pub n_very_rare: usize,
}

impl PrevalenceProfile {
    /// Get features above a prevalence threshold.
    pub fn features_above(&self, threshold: f64) -> Vec<usize> {
        self.feature_prevalence
            .iter()
            .enumerate()
            .filter(|(_, &p)| p >= threshold)
            .map(|(i, _)| i)
            .collect()
    }

    /// Get features below a prevalence threshold.
    pub fn features_below(&self, threshold: f64) -> Vec<usize> {
        self.feature_prevalence
            .iter()
            .enumerate()
            .filter(|(_, &p)| p < threshold)
            .map(|(i, _)| i)
            .collect()
    }

    /// Get features within a prevalence range.
    pub fn features_in_range(&self, min: f64, max: f64) -> Vec<usize> {
        self.feature_prevalence
            .iter()
            .enumerate()
            .filter(|(_, &p)| p >= min && p <= max)
            .map(|(i, _)| i)
            .collect()
    }
}

impl std::fmt::Display for PrevalenceProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Prevalence Profile")?;
        writeln!(f, "  Features:      {}", self.n_features)?;
        writeln!(f, "  Samples:       {}", self.n_samples)?;
        writeln!(f, "  Mean prevalence:   {:.2}%", self.mean_prevalence * 100.0)?;
        writeln!(f, "  Median prevalence: {:.2}%", self.median_prevalence * 100.0)?;
        writeln!(f, "  Min prevalence:    {:.2}%", self.min_prevalence * 100.0)?;
        writeln!(f, "  Max prevalence:    {:.2}%", self.max_prevalence * 100.0)?;
        writeln!(f, "  Ubiquitous (100%): {}", self.n_ubiquitous)?;
        writeln!(f, "  Singletons (1 sample): {}", self.n_singletons)?;
        writeln!(f, "  Rare (<10%):  {}", self.n_rare)?;
        writeln!(f, "  Very rare (<5%): {}", self.n_very_rare)?;
        Ok(())
    }
}

/// Profile prevalence characteristics of a count matrix.
pub fn profile_prevalence(counts: &CountMatrix) -> PrevalenceProfile {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    // Calculate prevalence per feature
    let feature_prevalence: Vec<f64> = (0..n_features)
        .map(|row| {
            let nnz = counts.data().outer_view(row)
                .map(|v| v.nnz())
                .unwrap_or(0);
            nnz as f64 / n_samples as f64
        })
        .collect();

    // Statistics
    let mean_prevalence = if n_features > 0 {
        feature_prevalence.iter().sum::<f64>() / n_features as f64
    } else {
        0.0
    };

    let median_prevalence = median(&feature_prevalence);
    let min_prevalence = feature_prevalence.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_prevalence = feature_prevalence.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    let n_ubiquitous = feature_prevalence.iter().filter(|&&p| p >= 1.0).count();
    let singleton_threshold = 1.0 / n_samples as f64 + 1e-10;
    let n_singletons = feature_prevalence.iter().filter(|&&p| p <= singleton_threshold && p > 0.0).count();
    let n_rare = feature_prevalence.iter().filter(|&&p| p < 0.10).count();
    let n_very_rare = feature_prevalence.iter().filter(|&&p| p < 0.05).count();

    PrevalenceProfile {
        n_features,
        n_samples,
        feature_prevalence,
        mean_prevalence,
        median_prevalence,
        min_prevalence: if min_prevalence.is_infinite() { 0.0 } else { min_prevalence },
        max_prevalence: if max_prevalence.is_infinite() { 0.0 } else { max_prevalence },
        n_ubiquitous,
        n_singletons,
        n_rare,
        n_very_rare,
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
        // 4 features × 4 samples
        let mut tri_mat = TriMat::new((4, 4));
        // Feature 0: present in all 4 samples (ubiquitous)
        tri_mat.add_triplet(0, 0, 10);
        tri_mat.add_triplet(0, 1, 20);
        tri_mat.add_triplet(0, 2, 15);
        tri_mat.add_triplet(0, 3, 12);
        // Feature 1: present in 2 samples (50%)
        tri_mat.add_triplet(1, 0, 5);
        tri_mat.add_triplet(1, 2, 8);
        // Feature 2: present in 1 sample (singleton, 25%)
        tri_mat.add_triplet(2, 1, 3);
        // Feature 3: present in 3 samples (75%)
        tri_mat.add_triplet(3, 0, 100);
        tri_mat.add_triplet(3, 1, 200);
        tri_mat.add_triplet(3, 3, 150);

        let feature_ids = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_profile_prevalence() {
        let counts = create_test_matrix();
        let profile = profile_prevalence(&counts);

        assert_eq!(profile.n_features, 4);
        assert_eq!(profile.n_samples, 4);
        assert!((profile.feature_prevalence[0] - 1.0).abs() < 1e-10);
        assert!((profile.feature_prevalence[1] - 0.5).abs() < 1e-10);
        assert!((profile.feature_prevalence[2] - 0.25).abs() < 1e-10);
        assert!((profile.feature_prevalence[3] - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_prevalence_counts() {
        let counts = create_test_matrix();
        let profile = profile_prevalence(&counts);

        assert_eq!(profile.n_ubiquitous, 1); // Feature 0
        assert_eq!(profile.n_singletons, 1); // Feature 2
    }

    #[test]
    fn test_features_above_threshold() {
        let counts = create_test_matrix();
        let profile = profile_prevalence(&counts);

        let above_50 = profile.features_above(0.5);
        assert_eq!(above_50, vec![0, 1, 3]);

        let above_75 = profile.features_above(0.75);
        assert_eq!(above_75, vec![0, 3]);
    }
}
