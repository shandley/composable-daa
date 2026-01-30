//! Data structures for differential abundance analysis.

mod count_matrix;
mod design_matrix;
mod formula;
mod metadata;
pub mod random_effects;
mod result;

pub use count_matrix::CountMatrix;
pub use design_matrix::DesignMatrix;
pub use formula::{Formula, Term};
pub use metadata::{Metadata, Variable, VariableType};
pub use random_effects::{MixedFormula, RandomDesignMatrix, RandomEffect};
pub use result::{Confidence, DaResult, DaResultSet, PrevalenceTier};
