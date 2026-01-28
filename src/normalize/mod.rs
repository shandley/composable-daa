//! Normalization methods for compositional data.

pub mod clr;
pub mod tss;

pub use clr::{norm_clr, TransformedMatrix};
pub use tss::{norm_tss, norm_tss_with_pseudocount, TssMatrix, scale};
