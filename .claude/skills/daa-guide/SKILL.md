---
name: daa-guide
description: Guide for differential abundance analysis method selection, threshold choices, and result interpretation. Use when user asks about LinDA, ZINB, Hurdle, NB methods, q-value thresholds, effect sizes, sensitivity, FDR, or how to interpret DAA results. Also use when user runs an analysis and gets unexpected results (like 0 significant features).
---

# Differential Abundance Analysis Guide

This skill provides evidence-based guidance for method selection, threshold choices, and result interpretation based on comprehensive benchmarking.

## Quick Method Selection

| Goal | Method | Threshold | Sensitivity | FDR | Use Case |
|------|--------|-----------|-------------|-----|----------|
| **Discovery** | ZINB | q < 0.05 | 83% | 29% | Maximize true positives |
| **Discovery** | Hurdle | q < 0.05 | 83% | 25% | Sparse data with structural zeros |
| **Confirmation** | LinDA | q < 0.10 | 39% | 12.5% | High-confidence findings |
| **Non-parametric** | Permutation | p < 0.05 | - | - | Unknown distributions |
| **Longitudinal** | LMM | p < 0.05 | - | - | Repeated measures |

## CRITICAL: LinDA Threshold

**LinDA requires q < 0.10, NOT q < 0.05**

LinDA uses CLR (Centered Log-Ratio) transformation which attenuates effect sizes by ~75%:
- True 4.0 log2FC (16x fold change) becomes ~1.0 observed
- At q < 0.05: **0% sensitivity** (nothing detected)
- At q < 0.10: **39% sensitivity** with excellent **12.5% FDR**

This is by design, not a bug. CLR centers by geometric mean to handle compositional data.

### When User Gets 0 Significant Results with LinDA

1. First, check if they used q < 0.05 → suggest q < 0.10
2. If still nothing at q < 0.10, the effects may be too small
3. LinDA needs >8x fold changes to detect anything reliably
4. Suggest ZINB or Hurdle for discovery if FDR control is less critical

## Effect Size Requirements

| Method | Minimum Detectable Effect | Notes |
|--------|---------------------------|-------|
| LinDA | >8x fold change (3 log2FC) | Due to CLR attenuation |
| ZINB | >2x fold change (1 log2FC) | Good sensitivity |
| Hurdle | >2x fold change (1 log2FC) | Best for sparse data |
| NB | >16x fold change (4 log2FC) | Conservative |

## Method Details

### LinDA (CLR + Linear Model)
- **Best for**: High-confidence findings, FDR-controlled discovery
- **Threshold**: q < 0.10 (NOT 0.05)
- **Pros**: Excellent FDR control (12.5%), handles compositionality
- **Cons**: Low sensitivity, only detects very large effects
- **Effect sizes**: NOT directly interpretable as fold changes (attenuated by ~75%)

### ZINB (Zero-Inflated Negative Binomial)
- **Best for**: Discovery, count data with excess zeros
- **Threshold**: q < 0.05
- **Pros**: High sensitivity (83%), models zero-inflation
- **Cons**: Moderate FDR (29%), assumes NB distribution

### Hurdle Model
- **Best for**: Sparse data with structural zeros, two-part analysis
- **Threshold**: q < 0.05
- **Pros**: High sensitivity (83%), good FDR (25%), separates presence/abundance
- **Cons**: More complex interpretation (binary + count components)

### NB (Negative Binomial)
- **Best for**: Low-sparsity data, simple overdispersion
- **Threshold**: q < 0.05
- **Pros**: Simple, well-understood
- **Cons**: Low sensitivity (6%), doesn't handle excess zeros

### Permutation Test
- **Best for**: Unknown distributions, non-parametric inference
- **Threshold**: p < 0.05
- **Pros**: Distribution-free, robust
- **Cons**: Computationally intensive, may be conservative

### LMM (Linear Mixed Model)
- **Best for**: Longitudinal data, repeated measures
- **Threshold**: p < 0.05
- **Pros**: Handles within-subject correlation
- **Cons**: Requires CLR transformation, same attenuation as LinDA

## Interpreting Results

### LinDA Results
- Effect sizes are CLR-transformed, NOT fold changes
- Observed estimate of 1.0 may represent true 4x fold change
- Focus on significance (q-value), not effect magnitude
- Use q < 0.10 threshold

### ZINB/Hurdle Results
- Effect sizes are on log scale (interpretable as log fold change)
- Higher sensitivity means more discoveries but also more false positives
- Consider biological plausibility of findings

### Sample Size Considerations
- n=20 per group: Only large effects (>4x) reliably detectable
- n=50 per group: Moderate effects (>2x) become detectable
- Power analysis recommended before study

## Decision Tree

```
Is your data longitudinal/repeated measures?
├── YES → Use LMM with q < 0.05
└── NO → Continue...

What is your primary goal?
├── DISCOVERY (find as many true effects as possible)
│   ├── High sparsity (>50% zeros)? → Hurdle (q < 0.05)
│   └── Moderate sparsity? → ZINB (q < 0.05)
├── CONFIRMATION (high-confidence findings only)
│   └── LinDA (q < 0.10)
└── UNSURE about distributional assumptions
    └── Permutation test (p < 0.05)
```

## For More Details

- Method comparison benchmarks: see [method-comparison.md](method-comparison.md)
- Q-value threshold analysis: see [thresholds.md](thresholds.md)
- Result interpretation guide: see [interpretation.md](interpretation.md)

## CLI Quick Reference

```bash
# Discovery (high sensitivity)
daa zinb -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv
daa hurdle -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Confirmation (low FDR) - NOTE: use q < 0.10 when filtering results
daa linda -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

# Non-parametric
daa permutation -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv
```
