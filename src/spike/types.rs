//! Core types for the spike-in framework.

use crate::data::{CountMatrix, PrevalenceTier};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

/// Type of spike-in effect applied.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SpikeType {
    /// Fold-change in abundance (multiplicative effect on non-zeros).
    Abundance,
    /// Change in presence/absence (zeros become non-zeros).
    Presence,
    /// Combined presence and abundance effect.
    Hurdle,
}

/// How to select features for spiking.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SpikeSelection {
    /// Select randomly from all eligible features.
    Random,
    /// Select from a specific prevalence tier.
    ByPrevalenceTier(PrevalenceTier),
    /// Select features within an abundance range (proportion of total).
    ByAbundance { min: f64, max: f64 },
    /// Select specific features by ID.
    Specific(Vec<String>),
}

/// Abundance level for presence spikes.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum AbundanceLevel {
    /// 10th percentile of non-zero values.
    Low,
    /// Median of non-zero values.
    Median,
    /// 90th percentile of non-zero values.
    High,
    /// Fixed count value.
    Fixed(u64),
}

/// Specification of a spike-in experiment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikeSpec {
    /// Type of spike applied.
    pub spike_type: SpikeType,
    /// IDs of features that were spiked.
    pub spiked_features: Vec<String>,
    /// True fold changes (for abundance spikes) or prevalence increases (for presence spikes).
    pub effect_sizes: Vec<f64>,
    /// Group that received the increased signal.
    pub affected_group: String,
    /// Original prevalence of each spiked feature in the affected group.
    pub original_prevalence: Vec<f64>,
    /// Random seed used.
    pub seed: u64,
}

impl SpikeSpec {
    /// Create a new spike specification.
    pub fn new(
        spike_type: SpikeType,
        spiked_features: Vec<String>,
        effect_sizes: Vec<f64>,
        affected_group: String,
        original_prevalence: Vec<f64>,
        seed: u64,
    ) -> Self {
        Self {
            spike_type,
            spiked_features,
            effect_sizes,
            affected_group,
            original_prevalence,
            seed,
        }
    }

    /// Get the set of spiked feature IDs for quick lookup.
    pub fn spiked_set(&self) -> HashSet<&str> {
        self.spiked_features.iter().map(|s| s.as_str()).collect()
    }

    /// Check if a feature was spiked.
    pub fn is_spiked(&self, feature_id: &str) -> bool {
        self.spiked_features.iter().any(|f| f == feature_id)
    }

    /// Get the effect size for a specific feature.
    pub fn effect_size(&self, feature_id: &str) -> Option<f64> {
        self.spiked_features
            .iter()
            .position(|f| f == feature_id)
            .map(|i| self.effect_sizes[i])
    }

    /// Number of spiked features.
    pub fn n_spiked(&self) -> usize {
        self.spiked_features.len()
    }
}

/// Result of a spike-in operation: modified counts plus specification.
#[derive(Debug, Clone)]
pub struct SpikedData {
    /// Modified count matrix with spikes injected.
    pub counts: CountMatrix,
    /// Specification of what was spiked.
    pub spec: SpikeSpec,
    /// Original count matrix (before spiking).
    pub original_counts: CountMatrix,
}

impl SpikedData {
    /// Create new spiked data.
    pub fn new(counts: CountMatrix, spec: SpikeSpec, original_counts: CountMatrix) -> Self {
        Self {
            counts,
            spec,
            original_counts,
        }
    }

    /// Get the spiked feature IDs.
    pub fn spiked_features(&self) -> &[String] {
        &self.spec.spiked_features
    }

    /// Check if a feature was spiked.
    pub fn is_spiked(&self, feature_id: &str) -> bool {
        self.spec.is_spiked(feature_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spike_spec() {
        let spec = SpikeSpec::new(
            SpikeType::Abundance,
            vec!["feat_1".into(), "feat_2".into()],
            vec![2.0, 3.0],
            "treatment".into(),
            vec![0.8, 0.6],
            42,
        );

        assert_eq!(spec.n_spiked(), 2);
        assert!(spec.is_spiked("feat_1"));
        assert!(!spec.is_spiked("feat_3"));
        assert_eq!(spec.effect_size("feat_1"), Some(2.0));
        assert_eq!(spec.effect_size("feat_2"), Some(3.0));
        assert_eq!(spec.effect_size("feat_3"), None);
    }
}
