//! Spike-in validation framework for empirical pipeline evaluation.
//!
//! This module provides tools to inject known signals into real data
//! and evaluate how well analysis pipelines recover those signals.

mod types;
mod abundance;
mod presence;
mod evaluate;
mod validate;

pub use types::{
    AbundanceLevel, SpikeDiagnostics, SpikeMode, SpikeSelection, SpikeSpec, SpikedData, SpikeType,
};
pub use abundance::{spike_abundance, spike_abundance_with_mode};
pub use presence::{spike_presence, spike_hurdle};
pub use evaluate::{evaluate_spikes, SpikeEvaluation, TierMetrics};
pub use validate::{
    run_spike_validation, run_spike_validation_presence, shuffle_labels,
    ValidationConfig, ValidationResult, ValidationSummary, PermutedEvaluation,
};
