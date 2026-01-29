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
All DAA methods are combinations of: filter → zero-handle → normalize → model → test → correct

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
│   ├── count_matrix.rs   # Sparse count matrix (features × samples)
│   ├── metadata.rs       # Sample metadata with type inference
│   ├── formula.rs        # R-style formula parsing
│   ├── design_matrix.rs  # Design matrix construction
│   └── result.rs         # DAResult, DaResultSet, Confidence, PrevalenceTier
├── profile/        # Data characterization
│   ├── sparsity.rs       # Zero patterns
│   ├── prevalence.rs     # Feature prevalence by group
│   └── library_size.rs   # Sequencing depth analysis
├── filter/         # Data filtering
│   ├── prevalence.rs     # Overall and groupwise prevalence filtering
│   ├── abundance.rs      # Count-based filtering
│   ├── library_size.rs   # Sample filtering by depth
│   └── stratified.rs     # Tiered filtering by prevalence
├── zero/           # Zero handling
│   └── pseudocount.rs    # Add pseudocount (fixed or adaptive)
├── normalize/      # Normalization methods
│   ├── clr.rs            # Centered log-ratio
│   └── tss.rs            # Total sum scaling
├── model/          # Statistical models
│   ├── lm.rs             # Linear model (QR decomposition)
│   ├── nb.rs             # Negative binomial GLM (IRLS)
│   ├── zinb.rs           # Zero-inflated NB (EM algorithm)
│   ├── compare.rs        # AIC/BIC model comparison
│   └── shrink.rs         # Empirical Bayes effect size shrinkage
├── test/           # Hypothesis testing
│   ├── wald.rs           # Wald test for coefficients
│   └── lrt.rs            # Likelihood ratio test
├── correct/        # Multiple testing correction
│   └── bh.rs             # Benjamini-Hochberg FDR
├── spike/          # Spike-in validation framework
│   ├── types.rs          # SpikeSpec, SpikeMode, etc.
│   ├── abundance.rs      # Abundance spike-ins
│   ├── presence.rs       # Presence/absence spike-ins
│   ├── evaluate.rs       # Sensitivity/FDR evaluation
│   ├── validate.rs       # Full validation runs
│   └── stress.rs         # Compositional stress testing
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
- Tests return `WaldResult` or `LrtResult` with p-values and statistics

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
Pipeline::new()
    .filter_prevalence(0.1)
    .add_pseudocount(0.5)
    .normalize_clr()
    .model_lm("~ group")
    .test_wald("grouptreatment")
    .correct_bh()
    .run(&counts, &metadata)
```

## Current Status (194 tests passing)

### Implemented
- Core data structures (CountMatrix, Metadata, Formula, DesignMatrix)
- Profiling (sparsity, prevalence, library size)
- Filtering (prevalence, abundance, library size, stratified)
- Zero handling (pseudocount)
- Normalization (CLR, TSS)
- Models (LM, NB, ZINB with model comparison and shrinkage)
- Testing (Wald, LRT)
- Correction (Benjamini-Hochberg)
- Spike-in validation (abundance, presence, stress testing)
- Pipeline composition with YAML serialization
- CLI tool

### Next Priorities
1. **Prevalence threshold optimization** - Find optimal thresholds via spike-in
2. **LLM-friendly profiling** - Structured output for AI-assisted design
3. **Group-specific prevalence handling** - Different thresholds per group

### Future
- TMM normalization
- Mixed models (LMM)
- Permutation tests
- Beta-binomial models
- Hurdle models

## Build & Test

```bash
cargo build
cargo test --lib
cargo run --bin daa -- --help
```

## CLI Usage

```bash
# Profile data
daa profile -c counts.tsv -m metadata.tsv -g group

# Run pipeline
daa run -c counts.tsv -m metadata.tsv -f "~ group" --test-coef grouptreatment

# Validate with spike-ins
daa validate -c counts.tsv -m metadata.tsv -g group -t treatment

# Stress test compositional effects
daa stress -c counts.tsv -m metadata.tsv -g group -t treatment
```

## Design Philosophy

1. **Primitives over methods** - Expose building blocks, not opinionated workflows
2. **Empirical over theoretical** - Validate on user's data, don't trust assumptions
3. **Sparse-first** - Design for 70-90% zeros from the start
4. **Composable** - Consistent interfaces, arbitrary chaining
5. **Prevalence-aware** - Groupwise prevalence as first-class concept
