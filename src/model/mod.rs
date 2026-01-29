//! Statistical models for differential abundance analysis.

pub mod lm;
pub mod nb;

pub use lm::{model_lm, LmFit, LmFitSingle};
pub use nb::{model_nb, NbFit, NbFitSingle};
