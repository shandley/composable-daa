//! Normalization methods for compositional data.
//!
//! This module provides several normalization approaches:
//!
//! - **CLR**: Centered log-ratio transformation (compositional)
//! - **ALR**: Additive log-ratio transformation (relative to reference taxon)
//! - **TSS**: Total sum scaling / relative abundance
//! - **TMM**: Trimmed mean of M-values (robust to asymmetric changes)
//! - **CSS**: Cumulative sum scaling (metagenomeSeq-style, robust to sparsity)
//! - **Spike-in**: Absolute abundance estimation using internal standards

pub mod alr;
pub mod clr;
pub mod css;
pub mod spikein;
pub mod tmm;
pub mod tss;

pub use alr::{
    norm_alr, norm_alr_default, norm_alr_with_pseudocount, AlrMatrix, ReferenceSelection,
};
pub use clr::{norm_clr, TransformedMatrix};
pub use css::{
    css_factors, estimate_css_quantile, norm_css, norm_css_with_config, CssConfig, CssMatrix,
};
pub use spikein::{
    detect_spikein_candidates, norm_spikein, norm_spikein_with_config, SpikeinConfig, SpikeinMatrix,
};
pub use tmm::{norm_tmm, norm_tmm_with_config, tmm_factors, tmm_factors_with_config, TmmConfig, TmmMatrix};
pub use tss::{norm_tss, norm_tss_with_pseudocount, scale, TssMatrix};
