//! Statistical hypothesis testing for differential abundance.

pub mod wald;

pub use wald::{test_wald, test_wald_nb, WaldResult, WaldResultSingle};
