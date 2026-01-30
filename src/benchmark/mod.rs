//! Benchmarking utilities for evaluating differential abundance methods.
//!
//! This module provides tools for:
//! - Generating synthetic datasets with known ground truth
//! - Fetching classic benchmark datasets from Zenodo
//! - Evaluating method performance

mod datasets;
mod generate;

pub use datasets::{
    fetch_dataset, list_datasets, clear_cache,
    BenchmarkDataset, DatasetInfo, FetchedDataset,
};
pub use generate::{
    generate_synthetic, Direction, EffectType, GroundTruth, SyntheticConfig, SyntheticData,
};
