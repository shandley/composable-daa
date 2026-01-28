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

/// How the spike-in handles compositional constraints.
///
/// Microbiome sequencing data is inherently compositional - we only observe
/// relative abundances, not absolute. Different spike modes model different
/// assumptions about what a "true" biological effect looks like in the data.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum SpikeMode {
    /// Raw multiplication of counts (default).
    ///
    /// Simply multiplies non-zero counts by the fold change. This increases
    /// library size for affected samples, which may not reflect biological
    /// reality but tests whether methods can detect inflated counts.
    #[default]
    Raw,

    /// Compositional mode: spike then renormalize to original library size.
    ///
    /// After applying the fold change, scales all counts so that the total
    /// library size remains unchanged. This models the scenario where one
    /// taxon increases in relative abundance while others decrease proportionally,
    /// as would be observed in real sequencing data with fixed depth.
    Compositional,

    /// Absolute mode: model what true absolute changes look like post-sequencing.
    ///
    /// Simulates a scenario where the spiked taxa truly increased in absolute
    /// abundance while others stayed constant, then the sample was sequenced
    /// at a fixed depth. This is the most biologically realistic model but
    /// results in smaller observed fold changes due to the compositional closure.
    Absolute,
}

/// Diagnostic information about the compositional effects of spiking.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpikeDiagnostics {
    /// Original geometric mean of counts (per affected sample, averaged).
    pub original_geometric_mean: f64,
    /// Geometric mean after spiking.
    pub spiked_geometric_mean: f64,
    /// Ratio of geometric means (indicates compositional shift).
    pub geometric_mean_ratio: f64,
    /// Nominal fold change requested.
    pub nominal_fold_change: f64,
    /// Effective CLR effect size (accounting for geometric mean shift).
    pub effective_clr_effect: f64,
    /// Library size change factor (1.0 for Compositional mode).
    pub library_size_factor: f64,
    /// Number of features spiked.
    pub n_spiked: usize,
    /// Total features in matrix.
    pub n_total_features: usize,
    /// Warning if compositional effects may dominate.
    pub compositional_warning: Option<String>,
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
    /// Spike mode used (Raw, Compositional, or Absolute).
    pub mode: SpikeMode,
    /// Diagnostic information about compositional effects.
    pub diagnostics: Option<SpikeDiagnostics>,
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
            mode: SpikeMode::Raw,
            diagnostics: None,
        }
    }

    /// Create a new spike specification with mode and diagnostics.
    pub fn with_diagnostics(
        spike_type: SpikeType,
        spiked_features: Vec<String>,
        effect_sizes: Vec<f64>,
        affected_group: String,
        original_prevalence: Vec<f64>,
        seed: u64,
        mode: SpikeMode,
        diagnostics: SpikeDiagnostics,
    ) -> Self {
        Self {
            spike_type,
            spiked_features,
            effect_sizes,
            affected_group,
            original_prevalence,
            seed,
            mode,
            diagnostics: Some(diagnostics),
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
