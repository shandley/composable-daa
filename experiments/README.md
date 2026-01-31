# Experimental Opportunities

This directory documents potential research studies that emerged during development of the composable-daa toolkit. These serve as case studies demonstrating the toolkit's capabilities and could form the basis of methods paper validation or standalone publications.

## Overview

| # | Experiment | Status | Priority | Data Available |
|---|------------|--------|----------|----------------|
| 01 | [BV Compositional Analysis](01-bv-compositional-analysis.md) | Exploratory complete | High | Yes (Ravel 2011) |
| 02 | [Spike-in Load Estimation](02-spikein-load-estimation.md) | Proof-of-concept | High | Yes (Stammler 2016) |
| 03 | [Compositional Artifact Audit](03-compositional-artifact-audit.md) | Conceptual | Medium | Requires meta-analysis |
| 04 | [Variance Estimation](04-variance-estimation.md) | **New** - Motivated by BB fix | High | Synthetic |
| 05 | [LinDA Sensitivity Analysis](05-linda-sensitivity.md) | **New** - Benchmarking finding | High | Synthetic |
| 06 | [Effect Size Recovery](06-effect-size-recovery.md) | **New** | Medium | Synthetic |
| 07 | [Library Size Confounding](07-library-size-confounding.md) | **New** | Medium | Synthetic |
| 08 | [Mixed Model Validation](08-mixed-model-validation.md) | **New** | Medium | Synthetic |

## Key Findings So Far

### 1. Total Load Variation is Massive
Using Stammler spike-in data, we demonstrated that samples vary **9.8x in total bacterial load**. This creates potential log2FC artifacts of **±3.3** - comparable to or larger than most reported effect sizes in the literature.

### 2. BV Signature is Compositionally Constrained
Analysis of Ravel 2011 BV data shows:
- 47/55 taxa significantly associated with BV (q<0.05)
- Sum of log2FC = 0.17 (essentially zero due to CLR closure)
- Classic findings (Lactobacillus decrease, Gardnerella increase) cannot be distinguished from total load effects

### 3. Statistical Robustness ≠ Biological Interpretability
Permutation tests show 0% FPR (well-calibrated statistics), but compositional constraints mean we cannot determine:
- Whether Lactobacillus truly decreases in absolute terms
- Whether BV-associated bacteria truly increase
- Or whether total load changes create the appearance of both

### 4. **NEW: Theoretical Variance Can Fail Catastrophically**
Benchmarking revealed that beta-binomial model-based standard errors underestimate variance by 100-1000x for microbiome data with large library sizes. This caused 98.5% FPR on null data (fixed by switching to empirical variance).

### 5. **NEW: LinDA May Be Underpowered**
On synthetic data with 1.0 log2FC effect (2-fold), LinDA showed 0% sensitivity while maintaining excellent FDR control. This suggests power analysis is critical for study design.

## New Experiments (from Benchmarking)

### Variance Estimation (Exp 04)
Why did theoretical beta-binomial variance fail so badly? When can we trust model-based SEs vs empirical variance? This is fundamental to all GLM-based methods.

**Key question**: Can we develop diagnostic tests to detect variance misspecification?

### LinDA Sensitivity (Exp 05)
Why does CLR + linear model fail to detect moderate effects? What sample sizes are needed? This is critical for study design guidance.

**Key question**: What is the power curve for LinDA across sparsity levels?

### Effect Size Recovery (Exp 06)
How accurately do methods recover true log2FC? Are estimates biased? Are CIs calibrated? This matters for meta-analysis and biological interpretation.

**Key question**: Does shrinkage improve effect size estimation?

### Library Size Confounding (Exp 07)
How do methods handle systematic library size differences between groups? This extends the spike-in load work to ask about DAA robustness.

**Key question**: Which normalization methods protect against confounding?

### Mixed Model Validation (Exp 08)
Does our LMM implementation properly control FPR for longitudinal data? How important is correct random effect specification?

**Key question**: What random effect structures are appropriate for microbiome studies?

## Directory Structure

```
experiments/
├── README.md                          # This file
├── 01-bv-compositional-analysis.md    # BV case study
├── 02-spikein-load-estimation.md      # Stammler spike-in
├── 03-compositional-artifact-audit.md # Meta-analysis proposal
├── 04-variance-estimation.md          # NEW: SE estimation
├── 05-linda-sensitivity.md            # NEW: Power analysis
├── 06-effect-size-recovery.md         # NEW: Effect calibration
├── 07-library-size-confounding.md     # NEW: Confounding
├── 08-mixed-model-validation.md       # NEW: LMM validation
├── preliminary-results.md             # Raw findings
└── scripts/
    ├── 01-bv-analysis/
    ├── 02-spikein-analysis/
    ├── 03-artifact-audit/
    ├── 04-variance-estimation/
    ├── 05-linda-sensitivity/
    ├── 06-effect-size-recovery/
    ├── 07-library-size-confounding/
    └── 08-mixed-model-validation/
```

## Relationship to Methods Paper

These experiments could serve as:
1. **Validation case studies** demonstrating toolkit capabilities
2. **Motivating examples** for why compositional-aware analysis matters
3. **Standalone findings** worthy of separate publication
4. **Benchmarking evidence** for method recommendations

The variance estimation and LinDA sensitivity experiments are particularly important for establishing when to use which method.

## Next Steps

1. ~~Complete main package architecture~~ DONE (368 tests passing)
2. ~~Fix critical BB model FPR issue~~ DONE (commit d60f594)
3. Formalize these analyses with reproducible scripts
4. Generate publication-quality figures
5. Consider integration with absolute quantification data (qPCR, flow cytometry)
6. Write methods paper demonstrating toolkit on these case studies

## Running Experiments

Most experiments use synthetic data generated by the toolkit:

```bash
# Generate typical 16S data
daa generate -p typical_16s -o data/typical_16s

# Generate sparse virome data
daa generate -p sparse_virome -o data/sparse_virome

# Generate confounded data (library size imbalance)
daa generate -p confounded -o data/confounded

# Fetch real benchmark data
daa fetch -d ravel -o data/ravel
daa fetch -d stammler -o data/stammler
```

See individual experiment files for specific analysis plans.
