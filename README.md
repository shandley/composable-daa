# composable-daa

A Rust library for **differential abundance analysis (DAA)** of sparse count data from microbiome and virome studies.

## Overview

Existing DAA methods (DESeq2, LinDA, ANCOM-BC, ALDEx2, etc.) are monolithic pipelines that make implicit choices about filtering, normalization, and modeling. This library takes a different approach: **expose the statistical primitives** and let users compose them into pipelines that can be **empirically validated** on their own data.

### Key Features

- **Composable primitives** - Mix and match filtering, normalization, modeling, and testing steps
- **Spike-in validation** - Inject known effects into your data to measure pipeline sensitivity and FDR
- **Prevalence threshold optimization** - Empirically find optimal filtering thresholds via spike-in analysis
- **AI-assisted pipeline design** - Generate structured data profiles for LLM-assisted method selection
- **Sparse-first design** - Built for 70-90% zero counts typical in virome data

## Installation

```bash
cargo install composable-daa
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
composable-daa = "0.1"
```

## Quick Start

### Command Line

```bash
# Run a LinDA-style analysis
daa linda -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Validate your pipeline with spike-ins
daa validate -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment

# Find optimal prevalence threshold
daa optimize-prevalence -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment
```

### Rust API

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
    println!("{}: log2FC={:.2}, q={:.4}",
        result.feature_id, result.estimate, result.q_value);
}
```

## Pipeline Components

### Filtering
- **Prevalence filtering** - Overall, group-wise (any/all/differential), or stratified by tier
- **Abundance filtering** - Minimum count, mean abundance, total count thresholds
- **Library size filtering** - Remove low-depth samples

### Normalization
- **CLR** - Centered log-ratio (compositional)
- **TSS** - Total sum scaling (relative abundance)
- **TMM** - Trimmed mean of M-values (edgeR-style)
- **Spike-in** - Absolute abundance via experimental spike-ins

### Models
- **Linear model** - For CLR-transformed data
- **Negative binomial GLM** - For raw counts with overdispersion
- **Zero-inflated NB** - For excess zeros beyond NB expectation

### Testing
- **Wald test** - Coefficient significance
- **Likelihood ratio test** - Nested model comparison
- **Permutation test** - Distribution-free, non-parametric

### Multiple Testing
- **Benjamini-Hochberg** - FDR control

## Spike-in Validation

The key innovation is empirical validation. Instead of trusting a method works, test it:

```bash
# Inject known effects, measure detection performance
daa validate -c counts.tsv -m metadata.tsv -g group -t treatment \
    -f "~ group" --test-coef grouptreatment --n-iterations 100

# Output includes:
# - Sensitivity (true positive rate)
# - FDR (false discovery rate)
# - Effect size recovery accuracy
# - Performance by prevalence tier
```

### Spike Modes

- **Compositional** (default) - Realistic: spiked features take reads from others
- **Absolute** - Adds reads, increases library size
- **Raw** - Direct multiplication (unrealistic but simple)

## CLI Commands

| Command | Description |
|---------|-------------|
| `profile` | Analyze data characteristics (sparsity, prevalence, library size) |
| `profile-llm` | Generate AI-friendly data profile for pipeline recommendations |
| `linda` | Run LinDA-style analysis (CLR + LM + Wald) |
| `permutation` | Run analysis with permutation tests |
| `run` | Execute pipeline from YAML configuration |
| `validate` | Spike-in validation to measure pipeline performance |
| `stress` | Compositional stress testing |
| `optimize-prevalence` | Find optimal prevalence threshold via spike-in |
| `generate` | Create synthetic benchmark data with known ground truth |
| `fetch` | Download classic benchmark datasets (Ravel BV, Stammler spike-in) |

## Data Formats

### Count Matrix (TSV)
```
feature_id	sample1	sample2	sample3
taxa_1	100	0	50
taxa_2	0	200	30
```

### Metadata (TSV)
```
sample_id	group	age
sample1	control	25
sample2	treatment	30
sample3	control	28
```

## Performance

- Sparse matrix representation (CSR format)
- Parallel processing via rayon
- Memory-efficient for large datasets

## Testing

```bash
cargo test           # Run all tests
cargo test --lib     # Library tests only
```

254 tests covering all components.

## License

MIT

## Citation

If you use this software, please cite:

```
@software{composable_daa,
  title = {composable-daa: Composable Differential Abundance Analysis},
  author = {Handley, Scott},
  year = {2025},
  url = {https://github.com/[username]/composable-daa}
}
```

## Related Work

This library implements concepts from:
- LinDA (Zhou et al., 2022)
- DESeq2 (Love et al., 2014)
- ANCOM-BC (Lin & Peddada, 2020)
- ALDEx2 (Fernandes et al., 2014)
- edgeR TMM (Robinson & Oshlack, 2010)

The spike-in validation framework is inspired by experimental spike-in studies like Stammler et al. (2016).
