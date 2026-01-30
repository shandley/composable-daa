//! Normalization methods for compositional data.
//!
//! This module provides several normalization approaches:
//!
//! - **CLR**: Centered log-ratio transformation (compositional)
//! - **TSS**: Total sum scaling / relative abundance
//! - **TMM**: Trimmed mean of M-values (robust to asymmetric changes)
//! - **Spike-in**: Absolute abundance estimation using internal standards

pub mod clr;
pub mod spikein;
pub mod tmm;
pub mod tss;

pub use clr::{norm_clr, TransformedMatrix};
pub use spikein::{
    detect_spikein_candidates, norm_spikein, norm_spikein_with_config, SpikeinConfig, SpikeinMatrix,
};
pub use tmm::{norm_tmm, norm_tmm_with_config, tmm_factors, tmm_factors_with_config, TmmConfig, TmmMatrix};
pub use tss::{norm_tss, norm_tss_with_pseudocount, scale, TssMatrix};
