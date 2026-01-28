//! Data profiling primitives for understanding count matrix characteristics.

mod library_size;
mod prevalence;
mod sparsity;

pub use library_size::{profile_library_size, LibrarySizeProfile};
pub use prevalence::{profile_prevalence, PrevalenceProfile};
pub use sparsity::{profile_sparsity, SparsityProfile};
