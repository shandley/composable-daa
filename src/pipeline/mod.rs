//! Pipeline composition and execution for differential abundance analysis.

mod runner;

pub use runner::{Pipeline, PipelineConfig, PipelineStep, run_linda};
