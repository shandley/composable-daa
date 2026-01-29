//! Statistical models for differential abundance analysis.

pub mod lm;
pub mod nb;
pub mod zinb;

pub use lm::{model_lm, LmFit, LmFitSingle};
pub use nb::{model_nb, NbFit, NbFitSingle};
pub use zinb::{model_zinb, model_zinb_with_config, ZinbConfig, ZinbFit, ZinbFitSingle};
