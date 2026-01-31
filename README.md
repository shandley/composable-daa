# composable-daa

A Rust library for **differential abundance analysis (DAA)** of sparse count data from microbiome and virome studies.

## Overview

Existing DAA methods (DESeq2, LinDA, ANCOM-BC, ALDEx2, etc.) are monolithic pipelines that make implicit choices about filtering, normalization, and modeling. This library takes a different approach: **expose the statistical primitives** and let users compose them into pipelines that can be **empirically validated** on their own data.

### Key Features

- **Automatic method selection** - Profiles your data and selects the optimal method
- **Composable primitives** - Mix and match filtering, normalization, modeling, and testing steps
- **Spike-in validation** - Inject known effects into your data to measure pipeline sensitivity and FDR
- **Longitudinal support** - Auto-detects repeated measures and applies appropriate mixed models
- **Sparse-first design** - Built for 70-90% zero counts typical in virome data

## Installation

```bash
cargo install composable-daa
```

Or build from source:

```bash
git clone https://github.com/shandley/composable-daa
cd composable-daa
cargo build --release
```

## Quick Start

### The Recommended Workflow

The simplest way to analyze your data is with the unified `recommend` command:

```bash
# Profile data, select method, and run analysis in one step
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv
```

This will:
1. Profile your data (sparsity, library size, group balance)
2. Detect study design (cross-sectional, longitudinal, repeated measures)
3. Select the appropriate method (LinDA, Hurdle, ZINB, or LMM)
4. Run the analysis and write results

### Example Output

```
=== Data Summary ===
Features:        500
Samples:         60
Sparsity:        72.4%
Min group size:  30

=== Recommendation ===
Method:    HURDLE
Threshold: q < 0.05

Rationale: High sparsity (>70%) - Hurdle model handles structural zeros well

=== Running HURDLE Analysis ===
  Formula: ~ group
  Testing: grouptreatment
Writing results to "results.tsv"...

=== Results Summary ===
Features tested:     342
Significant q<0.05:  15
  Up in target:      8
  Down in target:    7
```

### Longitudinal & Repeated Measures

The `recommend` command auto-detects longitudinal study designs:

```bash
# If metadata contains subject_id + timepoint columns:
# → Auto-runs LMM with: ~ group + timepoint + (1 | subject_id)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv

# If metadata contains subject_id only (repeated measures):
# → Auto-runs LMM with: ~ group + (1 | subject_id)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --run -o results.tsv
```

Detection patterns:
- **Subject columns**: subject, patient, individual, participant, person, donor
- **Time columns**: time, timepoint, visit, day, week, month, year

### Custom Pipelines

For more control, generate an editable YAML config:

```bash
# Generate YAML with detected covariates and comments
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --yaml -o pipeline.yaml

# Edit the YAML to customize (add interactions, change thresholds, etc.)
# Then run:
daa run -c counts.tsv -m metadata.tsv --config pipeline.yaml -o results.tsv
```

Example generated YAML:

```yaml
name: recommended_pipeline
steps:
  - filter_prevalence:
      threshold: 0.1
  - add_pseudocount:
      value: 0.5
  - normalize_clr: {}
  - model_lm:
      formula: "~ group"
  - test_wald:
      coefficient: "grouptreatment"
  - correct_bh:
      alpha: 0.05
```

## Method Selection

The `recommend` command uses these evidence-based rules:

| Sparsity | Method | Threshold | Rationale |
|----------|--------|-----------|-----------|
| >70% | Hurdle | q < 0.05 | Best for structural zeros |
| 50-70% | ZINB | q < 0.05 | Handles excess zeros well |
| <30% | LinDA | q < 0.10 | CLR + linear model works well |
| Longitudinal | LMM | q < 0.05 | Random effects for repeated measures |

**Important**: LinDA requires q < 0.10 (not 0.05) due to CLR effect attenuation.

## Validation

### Spike-in Validation

Test your pipeline's performance by injecting known effects:

```bash
daa validate -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment --n-iterations 100
```

This measures:
- **Sensitivity** - Can the pipeline detect true effects?
- **FDR** - How many false positives does it produce?
- **Effect size recovery** - Are estimates accurate?

### Compositional Stress Test

Quantify how dominant taxa affect results:

```bash
daa stress -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment
```

### Prevalence Threshold Optimization

Find the optimal filtering threshold empirically:

```bash
daa optimize-prevalence -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment
```

## CLI Commands

| Command | Description |
|---------|-------------|
| `recommend` | **Primary workflow** - Profile data, select method, optionally run analysis |
| `profile` | Analyze data characteristics (sparsity, prevalence, library size) |
| `profile-llm` | Generate structured data profile for AI-assisted analysis |
| `run` | Execute pipeline from YAML configuration |
| `permutation` | Run analysis with permutation tests (non-parametric) |
| `validate` | Spike-in validation to measure pipeline performance |
| `stress` | Compositional stress testing |
| `optimize-prevalence` | Find optimal prevalence threshold via spike-in |
| `generate` | Create synthetic benchmark data with known ground truth |
| `fetch` | Download classic benchmark datasets (Ravel BV, Stammler spike-in) |

### Command Examples

```bash
# Just get recommendation (no execution)
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment

# Quiet mode - just the command
daa recommend -c counts.tsv -m metadata.tsv -g group -t treatment --quiet

# Non-parametric permutation test
daa permutation -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Generate synthetic test data
daa generate -p typical_16s -o ./synthetic_data/

# Fetch classic benchmark datasets
daa fetch --list
daa fetch -d ravel -o ./ravel_data/
```

## Rust API

```rust
use composable_daa::prelude::*;

// Load data
let counts = CountMatrix::from_tsv("counts.tsv")?;
let metadata = Metadata::from_tsv("metadata.tsv")?;

// Build and run a pipeline
let results = Pipeline::new()
    .name("my-analysis")
    .filter_prevalence(0.1)
    .add_pseudocount(0.5)
    .normalize_clr()
    .model_lm("~ group")
    .test_wald("grouptreatment")
    .correct_bh()
    .run(&counts, &metadata)?;

// Get significant results
for result in results.significant() {
    println!("{}: estimate={:.2}, q={:.4}",
        result.feature_id, result.estimate, result.q_value);
}

// For longitudinal data with random effects
let results = Pipeline::new()
    .name("longitudinal-analysis")
    .filter_prevalence(0.1)
    .add_pseudocount(0.5)
    .normalize_clr()
    .model_lmm("~ group + time + (1 | subject)")
    .test_wald("grouptreatment")
    .correct_bh()
    .run(&counts, &metadata)?;
```

## Pipeline Components

### Filtering
- **Prevalence filtering** - Overall, group-wise (any/all/differential), or stratified
- **Abundance filtering** - Minimum count, mean abundance thresholds
- **Library size filtering** - Remove low-depth samples

### Normalization
- **CLR** - Centered log-ratio (compositional)
- **ALR** - Additive log-ratio (relative to reference)
- **TSS** - Total sum scaling (relative abundance)
- **TMM** - Trimmed mean of M-values (edgeR-style)
- **CSS** - Cumulative sum scaling (metagenomeSeq-style)

### Models
- **Linear model (LM)** - For CLR-transformed data
- **Linear mixed model (LMM)** - For longitudinal/repeated measures
- **Negative binomial GLM (NB)** - For raw counts with overdispersion
- **Zero-inflated NB (ZINB)** - For excess zeros beyond NB expectation
- **Hurdle model** - Two-part model for sparse counts

### Testing
- **Wald test** - Coefficient significance (all models)
- **Likelihood ratio test (LRT)** - Nested model comparison
- **Permutation test** - Distribution-free, non-parametric

### Multiple Testing
- **Benjamini-Hochberg** - FDR control

## Data Formats

### Count Matrix (TSV)
```
feature_id	sample1	sample2	sample3
taxa_1	100	0	50
taxa_2	0	200	30
```

### Metadata (TSV)
```
sample_id	group	subject	timepoint
sample1	control	subj1	0
sample2	treatment	subj1	1
sample3	control	subj2	0
```

## Interpreting Results

### Effect Sizes

**LinDA/LMM (CLR-transformed)**:
- Effect sizes are attenuated by ~75% due to CLR
- Multiply CLR estimate by ~4 to approximate true log2 fold change
- Use q < 0.10 threshold

**Hurdle/ZINB (count models)**:
- Effect sizes are on natural log scale
- `fold_change = exp(estimate)`
- Use q < 0.05 threshold

### Result Columns

| Column | Description |
|--------|-------------|
| `feature_id` | Feature identifier |
| `estimate` | Effect size (interpretation depends on method) |
| `std_error` | Standard error |
| `statistic` | Test statistic |
| `p_value` | Raw p-value |
| `q_value` | FDR-corrected q-value |
| `prevalence` | Proportion of samples with feature |
| `confidence` | high, moderate, suggestive, not_significant |

## Testing

```bash
cargo test           # Run all tests
cargo test --lib     # Library tests only
```

368+ tests covering all components.

## License

MIT

## Citation

If you use this software, please cite:

```
@software{composable_daa,
  title = {composable-daa: Composable Differential Abundance Analysis},
  author = {Handley, Scott},
  year = {2025},
  url = {https://github.com/shandley/composable-daa}
}
```

## Related Work

This library implements concepts from:
- LinDA (Zhou et al., 2022)
- DESeq2 (Love et al., 2014)
- ANCOM-BC (Lin & Peddada, 2020)
- ALDEx2 (Fernandes et al., 2014)
- edgeR TMM (Robinson & Oshlack, 2010)
- metagenomeSeq CSS (Paulson et al., 2013)

The spike-in validation framework is inspired by experimental spike-in studies like Stammler et al. (2016).
