//! Filtering primitives for count matrices.

pub mod prevalence;

pub use prevalence::{filter_prevalence_groupwise, filter_prevalence_overall, GroupwiseLogic};
