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

pub mod benchmark;
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
    pub use crate::benchmark::{
        generate_synthetic, Direction, EffectType, GroundTruth, SyntheticConfig, SyntheticData,
        // Dataset fetching
        fetch_dataset, list_datasets, clear_cache, BenchmarkDataset, DatasetInfo, FetchedDataset,
    };
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
    pub use crate::model::compare::{
        compare_nb_zinb, criteria_nb, criteria_zinb, select_best_model,
        ComparisonSummary, EvidenceStrength, FeatureComparison, ModelComparisonResult,
        ModelCriteria, SelectionCriterion,
    };
    pub use crate::model::lm::{model_lm, LmFit};
    pub use crate::model::nb::{model_nb, NbFit};
    pub use crate::model::shrink::{
        shrink_lfc, shrink_lfc_nb, shrink_lfc_zinb, ShrinkageConfig, ShrinkageMethod,
        ShrinkageResult, ShrunkEstimate,
    };
    pub use crate::model::zinb::{model_zinb, model_zinb_with_config, ZinbConfig, ZinbFit};
    pub use crate::normalize::alr::{
        norm_alr, norm_alr_default, norm_alr_with_pseudocount, AlrMatrix, ReferenceSelection,
    };
    pub use crate::normalize::clr::{norm_clr, TransformedMatrix};
    pub use crate::normalize::css::{
        css_factors, estimate_css_quantile, norm_css, norm_css_with_config, CssConfig, CssMatrix,
    };
    pub use crate::normalize::tmm::{norm_tmm, norm_tmm_with_config, tmm_factors, TmmConfig, TmmMatrix};
    pub use crate::normalize::tss::{norm_tss, norm_tss_with_pseudocount, TssMatrix, scale as tss_scale};
    pub use crate::normalize::spikein::{
        norm_spikein, norm_spikein_with_config, detect_spikein_candidates,
        SpikeinConfig, SpikeinMatrix,
    };
    pub use crate::pipeline::{Pipeline, PipelineStep, StratifiedPreset};
    pub use crate::profile::{
        LibrarySizeProfile, PrevalenceProfile, SparsityProfile,
        profile_library_size, profile_prevalence, profile_sparsity,
        // Compositionality
        CompositionalityProfile, DominanceCategory, DominantFeature,
        profile_compositionality,
        // LLM profiling
        LlmProfile, profile_for_llm,
    };
    pub use crate::test::lrt::{test_lrt_nb, test_lrt_zinb, test_lrt_nb_fitted, test_lrt_zinb_fitted, LrtResult};
    pub use crate::test::permutation::{
        test_permutation, test_permutation_quick, PermutationConfig, PermutationResult,
        PermutationResults,
    };
    pub use crate::test::wald::{test_wald, test_wald_nb, test_wald_zinb, test_wald_zinb_zi, WaldResult};
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
        // Prevalence optimization
        optimize_prevalence_threshold, PrevalenceOptConfig, PrevalenceOptResult,
        PrevalenceFilterLogic, OptimizationCriterion, ThresholdResult,
    };
}
