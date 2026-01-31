# Experiment 07: Library Size Confounding Analysis

## Summary

Systematic evaluation of how DAA methods handle confounding by library size (sequencing depth) differences between groups, and assessment of normalization effectiveness.

## Background

Library size (total reads per sample) often differs systematically between groups due to:
- Biological differences (e.g., bacterial load in disease)
- Technical variation (sample processing, extraction efficiency)
- Sequencing batch effects

### The Problem

If treatment group has 3x higher library sizes on average:
- All taxa appear more abundant in treatment (false positives)
- Normalization should correct this, but may fail
- Compositional methods may be more robust, or may amplify artifacts

### Our Synthetic Data

We generated "confounded" dataset with:
- 3x library size ratio between groups
- 15 true differential features (1.0 log2FC)
- 50% sparsity

## Research Questions

1. **How robust are methods to library size confounding?**
   - Which methods fail when library sizes differ?
   - What magnitude of imbalance causes problems?
   - Does sample size help or hurt?

2. **How effective are normalization methods?**
   - Does CLR correct for library size bias?
   - Do TMM/CSS/RLE help?
   - Is spike-in normalization superior?

3. **Can we detect confounding?**
   - What diagnostics reveal library size artifacts?
   - Can we correct post-hoc?
   - When should we trust results?

## Proposed Analysis

### Phase 1: Confounding Severity

1. Generate data with varying library size ratios:
   - Ratios: 1.0, 1.5, 2.0, 3.0, 5.0, 10.0
   - Balanced sample sizes (n=20 per group)
   - True effect: 1.0 log2FC

2. Run all methods without special handling:
   - Measure FPR on null features
   - Measure direction of bias
   - Assess sensitivity on true effects

### Phase 2: Normalization Comparison

For each confounding level:

1. Apply different normalizations:
   - None (raw proportions)
   - TSS (relative abundance)
   - CLR (centered log-ratio)
   - TMM (trimmed mean of M-values)
   - CSS (cumulative sum scaling)
   - Spike-in (if available)

2. Compare:
   - FPR control
   - Sensitivity preservation
   - Bias in effect estimates

### Phase 3: Diagnostic Development

1. Develop confounding diagnostics:
   - Correlation of effect with library size
   - Global shift detection
   - Asymmetry in significant taxa

2. Create correction procedures:
   - Library size as covariate
   - Post-hoc adjustment
   - Sensitivity analysis

## Data Requirements

- Synthetic data with controlled library size confounding
- Range of confounding magnitudes
- Real data examples with known confounding

## Expected Outcomes

1. Guidelines for handling library size confounding
2. Normalization method recommendations
3. Diagnostic tools for detecting confounding
4. Methods paper section on technical variation

## Relationship to Experiment 02

This extends the spike-in load estimation work by asking:
- Beyond measuring load variation, how does it affect DAA?
- Can normalization correct for load-driven artifacts?
- When is spike-in normalization necessary vs optional?

## Code Location

Scripts: `experiments/scripts/07-library-size-confounding/`

## References

- McMurdie & Holmes (2014) - Waste not, want not: normalization
- Weiss et al. (2017) - Normalization and microbial differential abundance
- Vandeputte et al. (2017) - Quantitative microbiome profiling
