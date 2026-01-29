//! Compositionality profiling: dominance and evenness analysis.
//!
//! This module analyzes the compositional structure of count data,
//! measuring how dominated the community is by a few taxa and
//! the overall evenness of the distribution.

use crate::data::CountMatrix;
use serde::{Deserialize, Serialize};

/// A dominant feature with its mean proportion.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DominantFeature {
    /// Feature identifier.
    pub feature_id: String,
    /// Mean proportion across samples.
    pub mean_proportion: f64,
    /// Rank (1 = most dominant).
    pub rank: usize,
}

/// Compositionality profile for a count matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompositionalityProfile {
    /// Proportion of reads from top 1 feature.
    pub dominance_top_1: f64,
    /// Proportion of reads from top 3 features.
    pub dominance_top_3: f64,
    /// Proportion of reads from top 10 features.
    pub dominance_top_10: f64,
    /// Shannon evenness (0-1 scale, 1 = perfectly even).
    pub evenness: f64,
    /// Gini coefficient (0-1, 0 = perfectly even, 1 = maximally uneven).
    pub gini: f64,
    /// Top dominant features.
    pub dominant_features: Vec<DominantFeature>,
    /// Number of features analyzed.
    pub n_features: usize,
}

impl CompositionalityProfile {
    /// Categorize dominance level.
    pub fn dominance_category(&self) -> DominanceCategory {
        if self.dominance_top_3 < 0.20 {
            DominanceCategory::Low
        } else if self.dominance_top_3 < 0.40 {
            DominanceCategory::Moderate
        } else if self.dominance_top_3 < 0.60 {
            DominanceCategory::High
        } else {
            DominanceCategory::Extreme
        }
    }

    /// Get interpretation text for dominance level.
    pub fn dominance_interpretation(&self) -> &'static str {
        match self.dominance_category() {
            DominanceCategory::Low => {
                "Low dominance - compositional effects are minimal"
            }
            DominanceCategory::Moderate => {
                "Moderate dominance - some compositional artifacts possible"
            }
            DominanceCategory::High => {
                "High dominance - significant compositional artifacts likely"
            }
            DominanceCategory::Extreme => {
                "Extreme dominance - results need very careful interpretation"
            }
        }
    }
}

/// Dominance category based on top-3 proportion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DominanceCategory {
    /// <20% from top 3
    Low,
    /// 20-40% from top 3
    Moderate,
    /// 40-60% from top 3
    High,
    /// >60% from top 3
    Extreme,
}

impl DominanceCategory {
    /// Get string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            DominanceCategory::Low => "low",
            DominanceCategory::Moderate => "moderate",
            DominanceCategory::High => "high",
            DominanceCategory::Extreme => "extreme",
        }
    }
}

/// Profile compositionality of a count matrix.
///
/// Calculates dominance (proportion from top features) and evenness
/// metrics to characterize the compositional structure.
pub fn profile_compositionality(counts: &CountMatrix) -> CompositionalityProfile {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();

    if n_features == 0 || n_samples == 0 {
        return CompositionalityProfile {
            dominance_top_1: 0.0,
            dominance_top_3: 0.0,
            dominance_top_10: 0.0,
            evenness: 0.0,
            gini: 0.0,
            dominant_features: vec![],
            n_features: 0,
        };
    }

    // Calculate mean proportion for each feature across samples
    let mut feature_proportions: Vec<(usize, f64)> = Vec::with_capacity(n_features);

    for feat_idx in 0..n_features {
        let mut total_proportion = 0.0;
        let mut n_valid = 0;

        for sample_idx in 0..n_samples {
            // Get library size for this sample
            let lib_size: u64 = (0..n_features)
                .map(|f| counts.get(f, sample_idx))
                .sum();

            if lib_size > 0 {
                let count = counts.get(feat_idx, sample_idx);
                total_proportion += count as f64 / lib_size as f64;
                n_valid += 1;
            }
        }

        let mean_prop = if n_valid > 0 {
            total_proportion / n_valid as f64
        } else {
            0.0
        };

        feature_proportions.push((feat_idx, mean_prop));
    }

    // Sort by proportion descending
    feature_proportions.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    // Calculate dominance metrics
    let dominance_top_1 = feature_proportions.first().map(|(_, p)| *p).unwrap_or(0.0);

    let dominance_top_3: f64 = feature_proportions.iter().take(3).map(|(_, p)| p).sum();

    let dominance_top_10: f64 = feature_proportions.iter().take(10).map(|(_, p)| p).sum();

    // Get dominant features
    let dominant_features: Vec<DominantFeature> = feature_proportions
        .iter()
        .take(10)
        .enumerate()
        .filter(|(_, (_, p))| *p > 0.01) // Only include if >1%
        .map(|(rank, (feat_idx, prop))| DominantFeature {
            feature_id: counts.feature_ids()[*feat_idx].clone(),
            mean_proportion: *prop,
            rank: rank + 1,
        })
        .collect();

    // Calculate evenness metrics
    let proportions: Vec<f64> = feature_proportions.iter().map(|(_, p)| *p).collect();
    let evenness = shannon_evenness(&proportions);
    let gini = gini_coefficient(&proportions);

    CompositionalityProfile {
        dominance_top_1,
        dominance_top_3,
        dominance_top_10,
        evenness,
        gini,
        dominant_features,
        n_features,
    }
}

/// Calculate Shannon evenness (Pielou's J).
///
/// E = H / H_max = H / ln(S)
/// where H is Shannon diversity and S is species richness.
///
/// Returns value between 0 (completely uneven) and 1 (perfectly even).
fn shannon_evenness(proportions: &[f64]) -> f64 {
    if proportions.is_empty() {
        return 0.0;
    }

    // Count non-zero species
    let s = proportions.iter().filter(|&&p| p > 0.0).count();
    if s <= 1 {
        return 1.0; // Single species is "perfectly even" by definition
    }

    // Calculate Shannon diversity H = -sum(p_i * ln(p_i))
    let h: f64 = proportions
        .iter()
        .filter(|&&p| p > 0.0)
        .map(|&p| -p * p.ln())
        .sum();

    // Maximum diversity H_max = ln(S)
    let h_max = (s as f64).ln();

    if h_max > 0.0 {
        h / h_max
    } else {
        1.0
    }
}

/// Calculate Gini coefficient.
///
/// Measures inequality in distribution. 0 = perfectly equal, 1 = maximally unequal.
fn gini_coefficient(proportions: &[f64]) -> f64 {
    let n = proportions.len();
    if n <= 1 {
        return 0.0;
    }

    // Sort proportions
    let mut sorted: Vec<f64> = proportions.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Calculate Gini
    let mut sum_weighted = 0.0;
    let mut sum_props = 0.0;

    for (i, &p) in sorted.iter().enumerate() {
        sum_weighted += (i + 1) as f64 * p;
        sum_props += p;
    }

    if sum_props == 0.0 {
        return 0.0;
    }

    let gini = (2.0 * sum_weighted) / (n as f64 * sum_props) - (n as f64 + 1.0) / n as f64;
    gini.max(0.0).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_even_counts() -> CountMatrix {
        // 10 features, all with equal counts
        let mut tri_mat = TriMat::new((10, 5));
        for feat in 0..10 {
            for sample in 0..5 {
                tri_mat.add_triplet(feat, sample, 100);
            }
        }
        let feature_ids: Vec<String> = (0..10).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..5).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_dominated_counts() -> CountMatrix {
        // 10 features, one dominates
        let mut tri_mat = TriMat::new((10, 5));
        for sample in 0..5 {
            // First feature dominates (80% of reads)
            tri_mat.add_triplet(0, sample, 800);
            // Rest share 20%
            for feat in 1..10 {
                tri_mat.add_triplet(feat, sample, 22);
            }
        }
        let feature_ids: Vec<String> = (0..10).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..5).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_even_distribution() {
        let counts = create_even_counts();
        let profile = profile_compositionality(&counts);

        // Each of 10 features should have ~10%
        assert!((profile.dominance_top_1 - 0.1).abs() < 0.01);
        assert!((profile.dominance_top_3 - 0.3).abs() < 0.01);

        // Evenness should be high
        assert!(profile.evenness > 0.95);

        // Gini should be low
        assert!(profile.gini < 0.1);

        // With 10 even features, top 3 = 30%, which is Moderate (20-40%)
        assert_eq!(profile.dominance_category(), DominanceCategory::Moderate);
    }

    #[test]
    fn test_dominated_distribution() {
        let counts = create_dominated_counts();
        let profile = profile_compositionality(&counts);

        // First feature should dominate
        assert!(profile.dominance_top_1 > 0.75);
        assert!(profile.dominance_top_3 > 0.80);

        // Evenness should be low
        assert!(profile.evenness < 0.5);

        // Gini should be high
        assert!(profile.gini > 0.5);

        assert_eq!(profile.dominance_category(), DominanceCategory::Extreme);
    }

    #[test]
    fn test_dominant_features_list() {
        let counts = create_dominated_counts();
        let profile = profile_compositionality(&counts);

        // Should identify the dominant feature
        assert!(!profile.dominant_features.is_empty());
        assert_eq!(profile.dominant_features[0].feature_id, "feat_0");
        assert_eq!(profile.dominant_features[0].rank, 1);
        assert!(profile.dominant_features[0].mean_proportion > 0.75);
    }

    #[test]
    fn test_shannon_evenness_perfect() {
        // Perfect evenness
        let props = vec![0.25, 0.25, 0.25, 0.25];
        let e = shannon_evenness(&props);
        assert!((e - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_shannon_evenness_uneven() {
        // Very uneven
        let props = vec![0.97, 0.01, 0.01, 0.01];
        let e = shannon_evenness(&props);
        assert!(e < 0.3);
    }

    #[test]
    fn test_gini_perfect_equality() {
        let props = vec![0.2, 0.2, 0.2, 0.2, 0.2];
        let g = gini_coefficient(&props);
        assert!(g < 0.01);
    }

    #[test]
    fn test_gini_inequality() {
        let props = vec![0.0, 0.0, 0.0, 0.0, 1.0];
        let g = gini_coefficient(&props);
        assert!(g > 0.7);
    }

    #[test]
    fn test_dominance_categories() {
        assert_eq!(DominanceCategory::Low.as_str(), "low");
        assert_eq!(DominanceCategory::Extreme.as_str(), "extreme");
    }
}
