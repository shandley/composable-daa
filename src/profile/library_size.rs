//! Library size profiling for count matrices.

use crate::data::CountMatrix;
use serde::{Deserialize, Serialize};

/// Profile of library size characteristics in a count matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibrarySizeProfile {
    /// Number of samples.
    pub n_samples: usize,
    /// Library size (total counts) per sample.
    pub library_sizes: Vec<u64>,
    /// Mean library size.
    pub mean: f64,
    /// Median library size.
    pub median: f64,
    /// Standard deviation of library sizes.
    pub std_dev: f64,
    /// Minimum library size.
    pub min: u64,
    /// Maximum library size.
    pub max: u64,
    /// Coefficient of variation (std_dev / mean).
    pub cv: f64,
    /// Log2 fold change between max and min.
    pub log2_fold_range: f64,
    /// Number of samples below mean - 2*std_dev.
    pub n_low_depth: usize,
    /// Number of samples above mean + 2*std_dev.
    pub n_high_depth: usize,
}

impl LibrarySizeProfile {
    /// Check if library sizes are highly variable (CV > 0.5).
    pub fn is_highly_variable(&self) -> bool {
        self.cv > 0.5
    }

    /// Check if there's a large range in library sizes (> 4-fold).
    pub fn has_large_range(&self) -> bool {
        self.log2_fold_range > 2.0
    }

    /// Get indices of samples with library size below a threshold.
    pub fn samples_below(&self, threshold: u64) -> Vec<usize> {
        self.library_sizes
            .iter()
            .enumerate()
            .filter(|(_, &s)| s < threshold)
            .map(|(i, _)| i)
            .collect()
    }

    /// Get indices of samples with library size above a threshold.
    pub fn samples_above(&self, threshold: u64) -> Vec<usize> {
        self.library_sizes
            .iter()
            .enumerate()
            .filter(|(_, &s)| s > threshold)
            .map(|(i, _)| i)
            .collect()
    }
}

impl std::fmt::Display for LibrarySizeProfile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Library Size Profile")?;
        writeln!(f, "  Samples: {}", self.n_samples)?;
        writeln!(f, "  Mean:    {:.0}", self.mean)?;
        writeln!(f, "  Median:  {:.0}", self.median)?;
        writeln!(f, "  Std Dev: {:.0}", self.std_dev)?;
        writeln!(f, "  Min:     {}", self.min)?;
        writeln!(f, "  Max:     {}", self.max)?;
        writeln!(f, "  CV:      {:.2}", self.cv)?;
        writeln!(f, "  Log2 fold range: {:.2}", self.log2_fold_range)?;
        writeln!(f, "  Low depth samples:  {}", self.n_low_depth)?;
        writeln!(f, "  High depth samples: {}", self.n_high_depth)?;
        Ok(())
    }
}

/// Profile library size characteristics of a count matrix.
pub fn profile_library_size(counts: &CountMatrix) -> LibrarySizeProfile {
    let library_sizes = counts.col_sums();
    let n_samples = library_sizes.len();

    if n_samples == 0 {
        return LibrarySizeProfile {
            n_samples: 0,
            library_sizes: vec![],
            mean: 0.0,
            median: 0.0,
            std_dev: 0.0,
            min: 0,
            max: 0,
            cv: 0.0,
            log2_fold_range: 0.0,
            n_low_depth: 0,
            n_high_depth: 0,
        };
    }

    let mean = library_sizes.iter().sum::<u64>() as f64 / n_samples as f64;

    let variance = library_sizes
        .iter()
        .map(|&x| {
            let diff = x as f64 - mean;
            diff * diff
        })
        .sum::<f64>() / n_samples as f64;
    let std_dev = variance.sqrt();

    let median = median_u64(&library_sizes);
    let min = *library_sizes.iter().min().unwrap_or(&0);
    let max = *library_sizes.iter().max().unwrap_or(&0);

    let cv = if mean > 0.0 { std_dev / mean } else { 0.0 };
    let log2_fold_range = if min > 0 {
        (max as f64 / min as f64).log2()
    } else {
        f64::INFINITY
    };

    let low_threshold = mean - 2.0 * std_dev;
    let high_threshold = mean + 2.0 * std_dev;
    let n_low_depth = library_sizes
        .iter()
        .filter(|&&x| (x as f64) < low_threshold)
        .count();
    let n_high_depth = library_sizes
        .iter()
        .filter(|&&x| (x as f64) > high_threshold)
        .count();

    LibrarySizeProfile {
        n_samples,
        library_sizes,
        mean,
        median,
        std_dev,
        min,
        max,
        cv,
        log2_fold_range: if log2_fold_range.is_infinite() { f64::MAX } else { log2_fold_range },
        n_low_depth,
        n_high_depth,
    }
}

fn median_u64(values: &[u64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort();
    let n = sorted.len();
    if n % 2 == 0 {
        (sorted[n / 2 - 1] + sorted[n / 2]) as f64 / 2.0
    } else {
        sorted[n / 2] as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        // 3 features × 4 samples with varying library sizes
        let mut tri_mat = TriMat::new((3, 4));
        // Sample 0: library size = 100
        tri_mat.add_triplet(0, 0, 40);
        tri_mat.add_triplet(1, 0, 50);
        tri_mat.add_triplet(2, 0, 10);
        // Sample 1: library size = 200
        tri_mat.add_triplet(0, 1, 80);
        tri_mat.add_triplet(1, 1, 100);
        tri_mat.add_triplet(2, 1, 20);
        // Sample 2: library size = 150
        tri_mat.add_triplet(0, 2, 60);
        tri_mat.add_triplet(1, 2, 75);
        tri_mat.add_triplet(2, 2, 15);
        // Sample 3: library size = 50
        tri_mat.add_triplet(0, 3, 20);
        tri_mat.add_triplet(1, 3, 25);
        tri_mat.add_triplet(2, 3, 5);

        let feature_ids = vec!["A".into(), "B".into(), "C".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_profile_library_size() {
        let counts = create_test_matrix();
        let profile = profile_library_size(&counts);

        assert_eq!(profile.n_samples, 4);
        assert_eq!(profile.library_sizes, vec![100, 200, 150, 50]);
        assert_eq!(profile.min, 50);
        assert_eq!(profile.max, 200);
        assert!((profile.mean - 125.0).abs() < 1e-10);
    }

    #[test]
    fn test_library_size_variation() {
        let counts = create_test_matrix();
        let profile = profile_library_size(&counts);

        // 4-fold range (200/50 = 4), so log2 = 2.0
        assert!((profile.log2_fold_range - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_samples_below_threshold() {
        let counts = create_test_matrix();
        let profile = profile_library_size(&counts);

        let below_100 = profile.samples_below(100);
        assert_eq!(below_100, vec![3]); // Only sample 3 has size 50
    }
}
