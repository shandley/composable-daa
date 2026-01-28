//! Data structures for differential abundance analysis.

mod count_matrix;
mod design_matrix;
mod formula;
mod metadata;
mod result;

pub use count_matrix::CountMatrix;
pub use design_matrix::DesignMatrix;
pub use formula::{Formula, Term};
pub use metadata::{Metadata, Variable, VariableType};
pub use result::{Confidence, DaResult, DaResultSet, PrevalenceTier};
