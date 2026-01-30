//! Statistical models for differential abundance analysis.

pub mod bb;
pub mod compare;
pub mod lm;
pub mod lmm;
pub mod nb;
pub mod shrink;
pub mod zinb;

pub use compare::{
    compare_criteria, compare_nb_zinb, criteria_nb, criteria_zinb, select_best_model,
    ComparisonSummary, EvidenceStrength, FeatureComparison, ModelComparisonResult,
    ModelCriteria, SelectionCriterion,
};
pub use lm::{model_lm, LmFit, LmFitSingle};
pub use lmm::{fit_lmm_from_formula, model_lmm, DfMethod, LmmConfig, LmmFit, LmmFitSingle};
pub use nb::{model_nb, NbFit, NbFitSingle};
pub use shrink::{
    shrink_lfc, shrink_lfc_nb, shrink_lfc_zinb, ShrinkageConfig, ShrinkageMethod,
    ShrinkageResult, ShrunkEstimate,
};
pub use zinb::{model_zinb, model_zinb_with_config, ZinbConfig, ZinbFit, ZinbFitSingle};
pub use bb::{model_bb, model_bb_with_config, BbConfig, BbFit, BbFitSingle, DispersionMethod};
