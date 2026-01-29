//! Statistical hypothesis testing for differential abundance.

pub mod wald;

pub use wald::{test_wald, test_wald_nb, test_wald_zinb, test_wald_zinb_zi, WaldResult, WaldResultSingle};
