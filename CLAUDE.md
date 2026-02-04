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
│   ├── hurdle.rs         # Hurdle model (binary + truncated NB)
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

// Hurdle model for sparse data with structural zeros
Pipeline::new()
    .filter_prevalence(0.1)
    .model_hurdle("~ group")  // Binary + truncated NB
    .test_wald("grouptreatment")  // Tests count component by default
    .correct_bh()
    .run(&counts, &metadata)
```

## Current Status (368 tests passing)

### Implemented

**Core Infrastructure**
- Data structures: CountMatrix, Metadata, Formula, DesignMatrix, MixedFormula, RandomDesignMatrix
- Profiling: sparsity, prevalence, library size, compositionality, LLM-friendly output
- Filtering: prevalence (overall/groupwise), abundance, library size, stratified
- Zero handling: pseudocount (fixed or adaptive)
- Pipeline composition with YAML serialization
- CLI tool with full feature coverage
- Benchmarking: synthetic data generation, classic dataset fetcher (Zenodo)

**Normalization Methods**
- CLR (centered log-ratio)
- ALR (additive log-ratio with reference selection)
- TSS (total sum scaling)
- TMM (trimmed mean of M-values, edgeR-style)
- CSS (cumulative sum scaling, metagenomeSeq-style)
- Spike-in normalization for absolute abundance

**Statistical Models**
- Linear model (LM) - QR decomposition
- Linear mixed model (LMM) - REML estimation
  - Random intercepts and slopes: `(1 | subject)`, `(1 + time | subject)`
  - Satterthwaite and Kenward-Roger df approximations
  - Variance components, ICC, BLUPs
- Negative binomial GLM (NB) - IRLS fitting
- Zero-inflated negative binomial (ZINB) - EM algorithm
- Hurdle model - two-part model for sparse data
  - Binary component (logistic) + count component (truncated NB)
  - All zeros from binary process (vs ZINB mixture)
- Model comparison via AIC/BIC
- Effect size shrinkage for all GLM models

**Hypothesis Testing**
- Wald test (all models including LMM, hurdle components)
- Likelihood ratio test (NB, ZINB, hurdle)
- Permutation tests (distribution-free)
- Benjamini-Hochberg FDR correction

**Spike-in Validation Framework**
- Abundance and presence spike-ins
- Compositional stress testing
- Prevalence threshold optimization
- Multiple spike modes: Raw, Compositional, Absolute

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

# Get method recommendation based on data profile
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment

# Run the recommended analysis automatically
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv

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

## AI Assistant Instructions (Claude Code)

**CRITICAL: When users ask for differential abundance analysis, microbiome analysis, or virome analysis:**

### ALWAYS DO:
1. Use `/daa-diagnose` to profile data and get method recommendations
2. Use `/daa-guide` for method selection questions and threshold guidance
3. Execute analysis using ONLY the `daa` CLI commands
4. Interpret results using the benchmark data in the skills

### NEVER DO:
- Suggest installing R packages (DESeq2, edgeR, limma, ANCOM-BC, ALDEx2, etc.)
- Suggest Python packages (scipy, statsmodels, scikit-bio, etc.)
- Write custom analysis scripts in any language
- Download or install external statistical tools
- Make statistical recommendations not grounded in our benchmark data
- Guess at thresholds - use the empirically-validated values from `/daa-guide`

### Available CLI Commands

| Command | Purpose |
|---------|---------|
| `daa profile` | Basic data profiling |
| `daa profile-llm` | Structured profiling for AI analysis |
| `daa recommend` | **Profile data, recommend method, optionally execute with `--run`** |
| `daa permutation` | Non-parametric permutation tests |
| `daa validate` | Spike-in validation |
| `daa stress` | Compositional stress testing |
| `daa optimize-prevalence` | Find optimal prevalence threshold |
| `daa run` | Run pipeline from YAML config |

### Workflow

```
User asks for analysis
       ↓
daa recommend --run (profiles data, selects method, executes)
       ↓
/daa-interpret to understand results
       ↓
Offer validation with `daa validate`
```

### Quick Start

For users who want fast analysis:

```bash
# Profile data and run recommended analysis in one step
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv

# Just get recommendation (no execution)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment

# Generate editable YAML config (for custom formulas, covariates, longitudinal)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --yaml -o pipeline.yaml
# Then edit the YAML and run:
daa run -c counts.tsv -m metadata.tsv --config pipeline.yaml -o results.tsv

# Quiet mode - just the command/config
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --quiet
```

### Longitudinal & Repeated Measures

The `recommend` command auto-detects longitudinal study designs and runs LMM automatically:

```bash
# Longitudinal: metadata contains subject_id + timepoint columns
# → Auto-runs LMM with: ~ group + timepoint + (1 | subject_id)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv

# Repeated measures: metadata contains subject_id only
# → Auto-runs LMM with: ~ group + (1 | subject_id)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv

# For custom formulas (interactions, random slopes), use --yaml:
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --yaml -o pipeline.yaml
# Edit pipeline.yaml to add: ~ group * timepoint + (1 + timepoint | subject_id)
daa run -c counts.tsv -m metadata.tsv --config pipeline.yaml -o results.tsv
```

Detection patterns:
- **Subject columns**: subject, patient, individual, participant, person, donor
- **Time columns**: time, timepoint, visit, day, week, month, year

### Why This Matters

This package exists because external tools make arbitrary choices. We've empirically validated our methods with spike-in analysis. Using external tools defeats the purpose and reintroduces the problems we're solving.

## Experiments Directory

The `experiments/` directory documents research opportunities identified during development:

1. **BV Compositional Analysis** - Ravel 2011 reanalysis showing compositional constraints
2. **Spike-in Load Estimation** - Stammler 2016 demonstrating 9.8x load variation
3. **Compositional Artifact Audit** - Proposed meta-analysis of published findings

These serve as validation case studies and potential standalone publications. See `experiments/README.md` for details.

---

## Microbiome Reproducibility Atlas (Active Project)

**Location**: `experiments/atlas/` and `experiments/REPRODUCIBILITY_ATLAS.md`

### Project Vision

Systematic characterization of microbiome-disease associations across public datasets. **Key insight**: The "reproducibility crisis" is largely a power crisis—even the most replicated finding (Fusobacterium-CRC) doesn't reach significance (p=0.56) in proper meta-analysis.

### Key Files

| File | Purpose |
|------|---------|
| `experiments/REPRODUCIBILITY_ATLAS.md` | Master vision document (UPDATED) |
| `experiments/atlas/FRAMEWORK.md` | Methodology documentation |
| `experiments/atlas/VALIDATION_RESULTS.md` | Fusobacterium case study results |
| `experiments/atlas/PHASE1_RESULTS.md` | Phase 1 findings (UPDATED) |
| `experiments/atlas/scripts/taxonomy.py` | Taxonomic harmonization (Greengenes/SILVA parsing) |
| `experiments/atlas/scripts/meta_analysis.py` | DerSimonian-Laird random-effects meta-analysis |
| `experiments/atlas/scripts/power_analysis.py` | Sample size and MDE calculations |
| `experiments/atlas/scripts/reprocess_with_taxonomy.py` | CRC data reprocessing pipeline |

### Current Status (as of Feb 2026)

**Completed**:
- Phase 1A: OTU-level analysis (found 0% replication - methodologically flawed)
- Phase 1B: Genus-level reanalysis with taxonomic harmonization
- Fusobacterium validation case study
- Random-effects meta-analysis implementation
- Power analysis framework
- Documentation revision (from "replication scores" to "effect characterization")

**Key Finding**: Fusobacterium-CRC shows:
- Correct direction in 3/4 cohorts
- Pooled effect = -0.26, p = 0.56 (NOT significant)
- I² = 0% (consistent across studies, just underpowered)
- The "70x fold change" is driven by 2 outlier patients

**Next Steps**:
1. Expand to IBD, obesity, T2D with genus-level harmonization
2. Add shotgun validation (curatedMetagenomicData)
3. Build searchable effect database
4. Publication: "Why Most Microbiome Findings Don't Replicate"

### Running the Analysis

```bash
# Reprocess CRC data with taxonomy
python experiments/atlas/scripts/reprocess_with_taxonomy.py

# Run meta-analysis
python experiments/atlas/scripts/meta_analysis.py

# Power analysis
python experiments/atlas/scripts/power_analysis.py
```

### Data Location

- Raw data: `experiments/atlas/diseases/crc/data/` (downloaded from MicrobiomeHD)
- Genus-level: `experiments/atlas/diseases/crc/data_genus/`
- Results: `experiments/atlas/diseases/crc/results_genus/`
