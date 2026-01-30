//! Statistical hypothesis testing for differential abundance.

pub mod lrt;
pub mod permutation;
pub mod wald;

pub use lrt::{test_lrt_nb, test_lrt_nb_fitted, test_lrt_zinb, test_lrt_zinb_fitted, LrtResult, LrtResultSingle};
pub use permutation::{
    test_permutation, test_permutation_quick, PermutationConfig, PermutationResult,
    PermutationResults,
};
pub use wald::{test_wald, test_wald_nb, test_wald_zinb, test_wald_zinb_zi, WaldResult, WaldResultSingle};
