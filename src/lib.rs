//! Composable Differential Abundance Analysis (DAA) Library
//!
//! This library provides modular primitives for differential abundance analysis
//! of microbiome and other compositional data.
//!
//! # Overview
//!
//! The library is organized into composable modules:
//!
//! - **data**: Core data structures (CountMatrix, Metadata, Results)
//! - **profile**: Data profiling (sparsity, prevalence, library size)
//! - **filter**: Data filtering (prevalence-based)
//! - **zero**: Zero handling strategies (pseudocount)
//! - **normalize**: Normalization methods (CLR)
//! - **model**: Statistical models (linear model)
//! - **test**: Hypothesis testing (Wald test)
//! - **correct**: Multiple testing correction (Benjamini-Hochberg)
//! - **pipeline**: Pipeline composition and execution
//!
//! # Example
//!
//! ```no_run
//! use composable_daa::prelude::*;
//!
//! // Load data
//! let counts = CountMatrix::from_tsv("counts.tsv").unwrap();
//! let metadata = Metadata::from_tsv("metadata.tsv").unwrap();
//!
//! // Run analysis pipeline
//! let results = Pipeline::new()
//!     .filter_prevalence(0.1)
//!     .add_pseudocount(0.5)
//!     .normalize_clr()
//!     .model_lm("~ group")
//!     .test_wald("grouptreatment")
//!     .correct_bh()
//!     .run(&counts, &metadata)
//!     .unwrap();
//! ```

pub mod correct;
pub mod data;
pub mod error;
pub mod filter;
pub mod model;
pub mod normalize;
pub mod pipeline;
pub mod profile;
pub mod spike;
pub mod test;
pub mod zero;

/// Convenient re-exports for common usage.
pub mod prelude {
    pub use crate::correct::bh::correct_bh;
    pub use crate::data::{
        Confidence, CountMatrix, DaResult, DaResultSet, DesignMatrix, Formula, Metadata,
        PrevalenceTier, Term, Variable,
    };
    pub use crate::error::{DaaError, Result};
    pub use crate::filter::{
        // Prevalence filtering
        filter_prevalence_groupwise, filter_prevalence_overall, GroupwiseLogic, FilterResult,
        // Abundance filtering
        filter_abundance, filter_mean_abundance, filter_total_count, AbundanceFilterResult,
        // Library size filtering
        filter_library_size, filter_library_size_quantile, filter_library_size_relative,
        LibrarySizeFilterResult,
        // Stratified filtering
        filter_stratified, TierConfig, TierThresholds, StratifiedFilterResult,
    };
    pub use crate::model::lm::{model_lm, LmFit};
    pub use crate::model::nb::{model_nb, NbFit};
    pub use crate::normalize::clr::{norm_clr, TransformedMatrix};
    pub use crate::normalize::tss::{norm_tss, norm_tss_with_pseudocount, TssMatrix, scale as tss_scale};
    pub use crate::pipeline::{Pipeline, PipelineStep, StratifiedPreset};
    pub use crate::profile::{
        LibrarySizeProfile, PrevalenceProfile, SparsityProfile,
        profile_library_size, profile_prevalence, profile_sparsity,
    };
    pub use crate::test::wald::{test_wald, test_wald_nb, WaldResult};
    pub use crate::zero::pseudocount::add_pseudocount;
    pub use crate::spike::{
        AbundanceLevel, SpikeDiagnostics, SpikeEvaluation, SpikeMode, SpikeSelection,
        SpikeSpec, SpikeType, SpikedData, TierMetrics, ValidationConfig, ValidationResult,
        ValidationSummary, evaluate_spikes, run_spike_validation, shuffle_labels,
        spike_abundance, spike_abundance_with_mode, spike_presence,
        // Stress testing
        StressConfig, StressParams, StressRunResult, StressSummary,
        AggregatedMetrics, PowerPoint, ModeComparisonSummary, ModeStats,
        run_stress_test, generate_parameter_grid, aggregate_results,
        calculate_power_curves, compare_modes,
    };
}
