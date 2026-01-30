# Composable DAA - Project Context

## Purpose

A Rust library for **differential abundance analysis (DAA)** of sparse count data (microbiome/virome). The core innovation is exposing statistical primitives that can be composed into custom pipelines, validated empirically via spike-in analysis, and optimized for specific data characteristics.

## The Problem We're Solving

1. **Method selection is arbitrary** - Researchers pick DESeq2, LinDA, ANCOM-BC, etc. without knowing which works best for their data structure
2. **Prevalence filtering is a shot in the dark** - Thresholds like 0.1 or 0.05 are chosen arbitrarily with no empirical basis
3. **Sparse data breaks assumptions** - 70-90% zeros in virome data violates many statistical assumptions
4. **No empirical validation** - Methods are trusted without testing on the user's actual data structure

## Our Solution

### 1. Composable Primitives
All DAA methods are combinations of: filter -> zero-handle -> normalize -> model -> test -> correct

We expose each step as independent primitives that can be mixed and matched.

### 2. Spike-in Validation
Inject known effects into real data, run the pipeline, measure sensitivity/FDR. This empirically validates method performance on the user's actual data structure.

### 3. AI-Assisted Pipeline Design
Profile the data, generate structured summaries, and use LLM assistance to recommend appropriate pipeline configurations.

### 4. Prevalence Threshold Optimization
Use spike-in validation to find optimal prevalence thresholds rather than guessing. Support group-specific thresholds (healthy vs. disease may need different filters).

## Architecture

```
src/
├── data/           # Core data structures
│   ├── count_matrix.rs   # Sparse count matrix (features x samples)
│   ├── metadata.rs       # Sample metadata with type inference
│   ├── formula.rs        # R-style formula parsing
│   ├── design_matrix.rs  # Design matrix construction
│   ├── random_effects.rs # Random effects for mixed models (MixedFormula, RandomDesignMatrix)
│   └── result.rs         # DAResult, DaResultSet, Confidence, PrevalenceTier
├── profile/        # Data characterization
│   ├── sparsity.rs       # Zero patterns
│   ├── prevalence.rs     # Feature prevalence by group
│   ├── library_size.rs   # Sequencing depth analysis
│   ├── compositionality.rs # Dominance and evenness analysis
│   └── llm.rs            # LLM-friendly profiling for AI-assisted design
├── filter/         # Data filtering
│   ├── prevalence.rs     # Overall and groupwise prevalence filtering
│   ├── abundance.rs      # Count-based filtering
│   ├── library_size.rs   # Sample filtering by depth
│   └── stratified.rs     # Tiered filtering by prevalence
├── zero/           # Zero handling
│   └── pseudocount.rs    # Add pseudocount (fixed or adaptive)
├── normalize/      # Normalization methods
│   ├── clr.rs            # Centered log-ratio
│   ├── alr.rs            # Additive log-ratio (relative to reference)
│   ├── tss.rs            # Total sum scaling
│   ├── tmm.rs            # Trimmed mean of M-values (edgeR-style)
│   ├── css.rs            # Cumulative sum scaling (metagenomeSeq-style)
│   └── spikein.rs        # Spike-in normalization for absolute abundance
├── model/          # Statistical models
│   ├── lm.rs             # Linear model (QR decomposition)
│   ├── lmm.rs            # Linear mixed model (REML estimation)
│   ├── nb.rs             # Negative binomial GLM (IRLS)
│   ├── zinb.rs           # Zero-inflated NB (EM algorithm)
│   ├── compare.rs        # AIC/BIC model comparison
│   └── shrink.rs         # Empirical Bayes effect size shrinkage
├── test/           # Hypothesis testing
│   ├── wald.rs           # Wald test for coefficients
│   ├── lrt.rs            # Likelihood ratio test
│   └── permutation.rs    # Permutation tests (distribution-free)
├── correct/        # Multiple testing correction
│   └── bh.rs             # Benjamini-Hochberg FDR
├── spike/          # Spike-in validation framework
│   ├── types.rs          # SpikeSpec, SpikeMode, etc.
│   ├── abundance.rs      # Abundance spike-ins
│   ├── presence.rs       # Presence/absence spike-ins
│   ├── evaluate.rs       # Sensitivity/FDR evaluation
│   ├── validate.rs       # Full validation runs
│   ├── stress.rs         # Compositional stress testing
│   └── prevalence_opt.rs # Prevalence threshold optimization
├── benchmark/      # Benchmarking utilities
│   ├── generate.rs       # Synthetic data generation with ground truth
│   └── datasets.rs       # Classic benchmark dataset fetcher (Zenodo)
├── pipeline/       # Pipeline composition
│   └── runner.rs         # Pipeline builder and executor
├── error.rs        # Error types
├── lib.rs          # Public API and prelude
└── bin/
    └── daa.rs      # CLI tool
```

## Key Patterns

### Consistent Return Types
- Filtering returns `FilterResult { kept, removed, kept_ids, removed_ids }`
- Models return `*Fit` structs with coefficients, std errors, fit statistics
- Tests return `WaldResult`, `LrtResult`, or `PermutationResults` with p-values and statistics

### Prevalence-Aware Design
- `PrevalenceTier` enum: Ubiquitous, Moderate, GroupSpecific, Rare
- `Confidence` enum: High, Medium, Low, Untestable
- Results annotated with prevalence context

### Spike Modes
- `SpikeMode::Raw` - Multiply counts directly (unrealistic)
- `SpikeMode::Compositional` - Realistic: spiked features take from others
- `SpikeMode::Absolute` - Adds reads (increases library size)

### Pipeline Builder
```rust
// Standard linear model pipeline
Pipeline::new()
    .filter_prevalence(0.1)
    .add_pseudocount(0.5)
    .normalize_clr()
    .model_lm("~ group")
    .test_wald("grouptreatment")
    .correct_bh()
    .run(&counts, &metadata)

// Linear mixed model for longitudinal data
Pipeline::new()
    .filter_prevalence(0.1)
    .add_pseudocount(0.5)
    .normalize_clr()
    .model_lmm("~ group + time + (1 | subject)")  // Random intercept per subject
    .test_wald("grouptreatment")
    .correct_bh()
    .run(&counts, &metadata)
```

## Current Status (306 tests passing)

### Implemented
- Core data structures (CountMatrix, Metadata, Formula, DesignMatrix, MixedFormula, RandomDesignMatrix)
- Profiling (sparsity, prevalence, library size, LLM-friendly output)
- Filtering (prevalence, abundance, library size, stratified)
- Zero handling (pseudocount)
- Normalization (CLR, ALR, TSS, TMM, CSS, spike-in)
- Models (LM, LMM, NB, ZINB with model comparison and shrinkage)
- Linear mixed models (LMM) with REML estimation for longitudinal/repeated measures data
  - Random intercepts: `(1 | subject)`
  - Supports lme4-style formula syntax
  - Variance component estimation (G matrix, sigma²)
  - Intraclass correlation coefficient (ICC)
  - BLUPs for random effects
  - Random slopes: `(1 + time | subject)`, `(0 + time | subject)`
  - Correlation estimation between random intercept and slope
- Testing (Wald, LRT, permutation) - including test_wald_lmm for mixed models
- Correction (Benjamini-Hochberg)
- Spike-in validation (abundance, presence, stress testing, threshold optimization)
- Pipeline composition with YAML serialization
- Benchmarking (synthetic data generation, classic dataset fetcher)
- CLI tool with full feature coverage

### Future Work
- Satterthwaite/Kenward-Roger df approximation for LMM
- Beta-binomial models
- Hurdle models

## Build and Test

```bash
cargo build
cargo test --lib
cargo run --bin daa -- --help
```

## CLI Usage

```bash
# Profile data
daa profile -c counts.tsv

# Generate LLM-friendly profile for AI-assisted pipeline design
daa profile-llm -c counts.tsv -m metadata.tsv -g group

# Run LinDA-style analysis
daa linda -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Run with permutation tests (non-parametric)
daa permutation -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Run pipeline from YAML config
daa run -c counts.tsv -m metadata.tsv --config pipeline.yaml -o results.tsv

# Validate with spike-ins
daa validate -c counts.tsv -m metadata.tsv -g group -t treatment -f "~ group" --test-coef grouptreatment

# Stress test compositional effects
daa stress -c counts.tsv -m metadata.tsv -g group -t treatment -f "~ group" --test-coef grouptreatment

# Optimize prevalence threshold via spike-in
daa optimize-prevalence -c counts.tsv -m metadata.tsv -g group -t treatment -f "~ group" --test-coef grouptreatment

# Generate synthetic benchmark data
daa generate -p typical_16s -o ./synthetic_data/

# Fetch classic benchmark datasets
daa fetch --list
daa fetch -d ravel -o ./ravel_data/
```

## Design Philosophy

1. **Primitives over methods** - Expose building blocks, not opinionated workflows
2. **Empirical over theoretical** - Validate on user's data, don't trust assumptions
3. **Sparse-first** - Design for 70-90% zeros from the start
4. **Composable** - Consistent interfaces, arbitrary chaining
5. **Prevalence-aware** - Groupwise prevalence as first-class concept

## Experiments Directory

The `experiments/` directory documents research opportunities identified during development:

1. **BV Compositional Analysis** - Ravel 2011 reanalysis showing compositional constraints
2. **Spike-in Load Estimation** - Stammler 2016 demonstrating 9.8x load variation
3. **Compositional Artifact Audit** - Proposed meta-analysis of published findings

These serve as validation case studies and potential standalone publications. See `experiments/README.md` for details.
