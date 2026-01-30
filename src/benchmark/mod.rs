//! Benchmarking utilities for evaluating differential abundance methods.
//!
//! This module provides tools for generating synthetic datasets with known
//! ground truth and evaluating method performance.

mod generate;

pub use generate::{
    generate_synthetic, Direction, EffectType, GroundTruth, SyntheticConfig, SyntheticData,
};
