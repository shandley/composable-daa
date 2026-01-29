//! Statistical models for differential abundance analysis.

pub mod compare;
pub mod lm;
pub mod nb;
pub mod zinb;

pub use compare::{
    compare_criteria, compare_nb_zinb, criteria_nb, criteria_zinb, select_best_model,
    ComparisonSummary, EvidenceStrength, FeatureComparison, ModelComparisonResult,
    ModelCriteria, SelectionCriterion,
};
pub use lm::{model_lm, LmFit, LmFitSingle};
pub use nb::{model_nb, NbFit, NbFitSingle};
pub use zinb::{model_zinb, model_zinb_with_config, ZinbConfig, ZinbFit, ZinbFitSingle};
