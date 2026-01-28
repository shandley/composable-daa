//! Stratified filtering by prevalence tier.
//!
//! Allows different filtering thresholds for features in different prevalence tiers,
//! useful for applying more lenient criteria to rare features.

use crate::data::{CountMatrix, PrevalenceTier};
use crate::error::{DaaError, Result};
use crate::profile::profile_prevalence;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for per-tier filtering thresholds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TierThresholds {
    /// Thresholds by prevalence tier (prevalence, min_abundance, min_samples).
    thresholds: HashMap<PrevalenceTier, TierConfig>,
    /// Default config for tiers not explicitly specified.
    default: TierConfig,
}

/// Configuration for a single tier.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TierConfig {
    /// Minimum prevalence to keep feature (0.0-1.0).
    pub min_prevalence: f64,
    /// Minimum mean abundance when present.
    pub min_abundance: f64,
    /// Minimum number of samples with non-zero counts.
    pub min_samples: usize,
}

impl Default for TierConfig {
    fn default() -> Self {
        Self {
            min_prevalence: 0.0,
            min_abundance: 0.0,
            min_samples: 1,
        }
    }
}

impl TierThresholds {
    /// Create new tier thresholds with default config.
    pub fn new() -> Self {
        Self {
            thresholds: HashMap::new(),
            default: TierConfig::default(),
        }
    }

    /// Set threshold for a specific tier.
    pub fn with_tier(mut self, tier: PrevalenceTier, config: TierConfig) -> Self {
        self.thresholds.insert(tier, config);
        self
    }

    /// Set the default threshold for unspecified tiers.
    pub fn with_default(mut self, config: TierConfig) -> Self {
        self.default = config;
        self
    }

    /// Get config for a tier.
    pub fn get(&self, tier: &PrevalenceTier) -> &TierConfig {
        self.thresholds.get(tier).unwrap_or(&self.default)
    }

    /// Create a commonly-used configuration that's lenient on rare features.
    ///
    /// This is useful when you want to retain rare features that might be
    /// biologically interesting while still filtering out noise.
    pub fn lenient_rare() -> Self {
        Self::new()
            .with_tier(
                PrevalenceTier::VeryHigh,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 1.0,
                    min_samples: 2,
                },
            )
            .with_tier(
                PrevalenceTier::High,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 1.0,
                    min_samples: 2,
                },
            )
            .with_tier(
                PrevalenceTier::Medium,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 2.0,
                    min_samples: 2,
                },
            )
            .with_tier(
                PrevalenceTier::Low,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 5.0,
                    min_samples: 2,
                },
            )
            .with_tier(
                PrevalenceTier::Rare,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 10.0,
                    min_samples: 2,
                },
            )
    }

    /// Create a strict configuration that removes most rare features.
    pub fn strict() -> Self {
        Self::new()
            .with_tier(
                PrevalenceTier::VeryHigh,
                TierConfig {
                    min_prevalence: 0.5,
                    min_abundance: 1.0,
                    min_samples: 5,
                },
            )
            .with_tier(
                PrevalenceTier::High,
                TierConfig {
                    min_prevalence: 0.3,
                    min_abundance: 2.0,
                    min_samples: 5,
                },
            )
            .with_tier(
                PrevalenceTier::Medium,
                TierConfig {
                    min_prevalence: 0.2,
                    min_abundance: 5.0,
                    min_samples: 5,
                },
            )
            .with_tier(
                PrevalenceTier::Low,
                TierConfig {
                    min_prevalence: 0.1,
                    min_abundance: 10.0,
                    min_samples: 3,
                },
            )
            .with_tier(
                PrevalenceTier::Rare,
                TierConfig {
                    min_prevalence: 0.05,
                    min_abundance: 20.0,
                    min_samples: 2,
                },
            )
    }
}

impl Default for TierThresholds {
    fn default() -> Self {
        Self::new()
    }
}

/// Filter features using stratified thresholds by prevalence tier.
///
/// Each feature is assigned to a prevalence tier, then filtered according
/// to the tier-specific thresholds.
///
/// # Arguments
/// * `counts` - The count matrix to filter
/// * `thresholds` - Per-tier filtering configuration
///
/// # Returns
/// A new CountMatrix containing only features meeting their tier's criteria.
pub fn filter_stratified(
    counts: &CountMatrix,
    thresholds: &TierThresholds,
) -> Result<CountMatrix> {
    let prevalence_profile = profile_prevalence(counts);

    let keep_indices: Vec<usize> = (0..counts.n_features())
        .into_par_iter()
        .filter(|&row| {
            let prevalence = prevalence_profile.feature_prevalence[row];
            let tier = PrevalenceTier::from_prevalence(prevalence);
            let config = thresholds.get(&tier);

            // Check prevalence threshold
            if prevalence < config.min_prevalence {
                return false;
            }

            // Check minimum samples
            let row_data = counts.row_dense(row);
            let nnz = row_data.iter().filter(|&&v| v > 0).count();
            if nnz < config.min_samples {
                return false;
            }

            // Check mean abundance
            if nnz > 0 {
                let sum: u64 = row_data.iter().filter(|&&v| v > 0).sum();
                let mean = sum as f64 / nnz as f64;
                if mean < config.min_abundance {
                    return false;
                }
            }

            true
        })
        .collect();

    if keep_indices.is_empty() {
        return Err(DaaError::EmptyData(
            "No features pass stratified filtering criteria".to_string(),
        ));
    }

    counts.subset_features(&keep_indices)
}

/// Result of stratified filtering with per-tier statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StratifiedFilterResult {
    /// Number of features before filtering.
    pub n_before: usize,
    /// Number of features after filtering.
    pub n_after: usize,
    /// Breakdown by tier.
    pub by_tier: HashMap<String, TierFilterStats>,
}

/// Statistics for a single tier.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TierFilterStats {
    /// Features in this tier before filtering.
    pub n_before: usize,
    /// Features in this tier after filtering.
    pub n_after: usize,
    /// Retention rate for this tier.
    pub retention_rate: f64,
}

impl std::fmt::Display for StratifiedFilterResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Stratified Filter Result")?;
        writeln!(f, "  Total before: {}", self.n_before)?;
        writeln!(f, "  Total after:  {}", self.n_after)?;
        writeln!(f, "  By tier:")?;
        for (tier, stats) in &self.by_tier {
            writeln!(
                f,
                "    {}: {} → {} ({:.1}%)",
                tier,
                stats.n_before,
                stats.n_after,
                stats.retention_rate * 100.0
            )?;
        }
        Ok(())
    }
}

/// Filter with per-tier statistics.
pub fn filter_stratified_with_stats(
    counts: &CountMatrix,
    thresholds: &TierThresholds,
) -> Result<(CountMatrix, StratifiedFilterResult)> {
    let n_before = counts.n_features();
    let prevalence_profile = profile_prevalence(counts);

    // Count features per tier before filtering
    let mut before_by_tier: HashMap<String, usize> = HashMap::new();
    for &prev in &prevalence_profile.feature_prevalence {
        let tier = PrevalenceTier::from_prevalence(prev);
        *before_by_tier.entry(tier.name().to_string()).or_default() += 1;
    }

    let filtered = filter_stratified(counts, thresholds)?;

    // Count features per tier after filtering
    let filtered_prevalence = profile_prevalence(&filtered);
    let mut after_by_tier: HashMap<String, usize> = HashMap::new();
    for &prev in &filtered_prevalence.feature_prevalence {
        let tier = PrevalenceTier::from_prevalence(prev);
        *after_by_tier.entry(tier.name().to_string()).or_default() += 1;
    }

    // Build per-tier stats
    let mut by_tier = HashMap::new();
    for (tier_name, &n_before_tier) in &before_by_tier {
        let n_after_tier = after_by_tier.get(tier_name).copied().unwrap_or(0);
        by_tier.insert(
            tier_name.clone(),
            TierFilterStats {
                n_before: n_before_tier,
                n_after: n_after_tier,
                retention_rate: if n_before_tier > 0 {
                    n_after_tier as f64 / n_before_tier as f64
                } else {
                    1.0
                },
            },
        );
    }

    let result = StratifiedFilterResult {
        n_before,
        n_after: filtered.n_features(),
        by_tier,
    };

    Ok((filtered, result))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        // 10 features × 20 samples with varying prevalence
        let mut tri_mat = TriMat::new((10, 20));

        // Features 0-1: very high prevalence (100%)
        for feat in 0..2 {
            for sample in 0..20 {
                tri_mat.add_triplet(feat, sample, 100);
            }
        }
        // Features 2-3: high prevalence (80%)
        for feat in 2..4 {
            for sample in 0..16 {
                tri_mat.add_triplet(feat, sample, 50);
            }
        }
        // Features 4-5: medium prevalence (50%)
        for feat in 4..6 {
            for sample in 0..10 {
                tri_mat.add_triplet(feat, sample, 30);
            }
        }
        // Features 6-7: low prevalence (20%)
        for feat in 6..8 {
            for sample in 0..4 {
                tri_mat.add_triplet(feat, sample, 20);
            }
        }
        // Features 8-9: rare prevalence (5%)
        tri_mat.add_triplet(8, 0, 10);
        tri_mat.add_triplet(9, 0, 5);

        let feature_ids: Vec<String> = (0..10).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..20).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_tier_thresholds_default() {
        let thresholds = TierThresholds::new();
        let config = thresholds.get(&PrevalenceTier::High);

        // Default should be very permissive
        assert_eq!(config.min_prevalence, 0.0);
        assert_eq!(config.min_abundance, 0.0);
        assert_eq!(config.min_samples, 1);
    }

    #[test]
    fn test_filter_stratified_basic() {
        let counts = create_test_matrix();

        // Use default thresholds (very permissive)
        let thresholds = TierThresholds::new();
        let filtered = filter_stratified(&counts, &thresholds).unwrap();

        // All features should pass with default thresholds
        assert_eq!(filtered.n_features(), 10);
    }

    #[test]
    fn test_filter_stratified_strict() {
        let counts = create_test_matrix();

        // Use strict thresholds
        let thresholds = TierThresholds::strict();
        let filtered = filter_stratified(&counts, &thresholds).unwrap();

        // Should filter out rare features
        assert!(filtered.n_features() < 10);
        // High prevalence features should remain
        assert!(filtered.n_features() >= 2);
    }

    #[test]
    fn test_filter_stratified_custom() {
        let counts = create_test_matrix();

        // Custom thresholds: require min 5 samples for rare features
        let thresholds = TierThresholds::new()
            .with_tier(
                PrevalenceTier::Rare,
                TierConfig {
                    min_prevalence: 0.0,
                    min_abundance: 0.0,
                    min_samples: 5, // Rare features need 5 samples
                },
            );

        let filtered = filter_stratified(&counts, &thresholds).unwrap();

        // Rare features (8, 9) have only 1 sample, should be filtered
        assert_eq!(filtered.n_features(), 8);
    }

    #[test]
    fn test_filter_stratified_with_stats() {
        let counts = create_test_matrix();

        let thresholds = TierThresholds::lenient_rare();
        let (filtered, stats) = filter_stratified_with_stats(&counts, &thresholds).unwrap();

        assert_eq!(stats.n_before, 10);
        assert!(stats.n_after <= 10);
        assert!(!stats.by_tier.is_empty());

        // Check that very_high tier has high retention
        if let Some(vh_stats) = stats.by_tier.get("very_high") {
            assert!(vh_stats.retention_rate >= 0.5);
        }
    }

    #[test]
    fn test_lenient_rare_preset() {
        let thresholds = TierThresholds::lenient_rare();

        // Rare tier should have higher abundance requirement
        let rare_config = thresholds.get(&PrevalenceTier::Rare);
        let high_config = thresholds.get(&PrevalenceTier::High);

        assert!(rare_config.min_abundance > high_config.min_abundance);
    }
}
