//! Data profiling primitives for understanding count matrix characteristics.

pub mod compositionality;
mod library_size;
pub mod llm;
mod prevalence;
mod sparsity;

pub use compositionality::{
    profile_compositionality, CompositionalityProfile, DominanceCategory, DominantFeature,
};
pub use library_size::{profile_library_size, LibrarySizeProfile};
pub use llm::{profile_for_llm, LlmProfile};
pub use prevalence::{profile_prevalence, PrevalenceProfile};
pub use sparsity::{profile_sparsity, SparsityProfile};
