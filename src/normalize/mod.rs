//! Normalization methods for compositional data.
//!
//! This module provides several normalization approaches:
//!
//! - **CLR**: Centered log-ratio transformation (compositional)
//! - **TSS**: Total sum scaling / relative abundance
//! - **Spike-in**: Absolute abundance estimation using internal standards

pub mod clr;
pub mod spikein;
pub mod tss;

pub use clr::{norm_clr, TransformedMatrix};
pub use spikein::{
    detect_spikein_candidates, norm_spikein, norm_spikein_with_config,
    SpikeinConfig, SpikeinMatrix,
};
pub use tss::{norm_tss, norm_tss_with_pseudocount, scale, TssMatrix};
