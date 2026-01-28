//! Filtering primitives for count matrices.
//!
//! This module provides various filtering strategies:
//!
//! - **Prevalence filtering**: Remove features present in too few samples
//! - **Abundance filtering**: Remove features by relative abundance
//! - **Library size filtering**: Remove samples with too few/many reads
//! - **Stratified filtering**: Apply different thresholds per prevalence tier

pub mod abundance;
pub mod library_size;
pub mod prevalence;
pub mod stratified;

pub use abundance::{
    filter_abundance, filter_abundance_with_stats, filter_mean_abundance, filter_total_count,
    AbundanceFilterResult,
};
pub use library_size::{
    filter_library_size, filter_library_size_quantile, filter_library_size_relative,
    filter_library_size_with_metadata, filter_library_size_with_stats, LibrarySizeFilterResult,
};
pub use prevalence::{
    filter_prevalence_groupwise, filter_prevalence_overall, filter_prevalence_overall_with_stats,
    FilterResult, GroupwiseLogic,
};
pub use stratified::{
    filter_stratified, filter_stratified_with_stats, StratifiedFilterResult, TierConfig,
    TierFilterStats, TierThresholds,
};
