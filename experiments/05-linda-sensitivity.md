# Experiment 05: LinDA Sensitivity Analysis

## Summary

Investigation of why LinDA (CLR + linear model) shows 0% sensitivity on typical microbiome data with moderate effect sizes, and development of power guidelines.

## Background

During benchmarking on synthetic data with known ground truth:

| Dataset | Effect Size | LinDA Sensitivity | BB Sensitivity |
|---------|-------------|-------------------|----------------|
| typical_16s (64% sparse) | 1.0 log2FC | 0% | 60% |
| sparse_virome (89% sparse) | 1.5 log2FC | 0% | 40% |
| group_specific | 1.0 log2FC | 45% | 70% |

LinDA detected nothing on two datasets despite true effects being present.

### Possible Explanations

1. **CLR transformation dampens signal**
   - CLR centers on geometric mean
   - Effects on low-abundance taxa are compressed
   - High-prevalence features dominate

2. **Effect size too small for sample size**
   - n=20 per group may be underpowered
   - 1.0 log2FC (2-fold) requires larger n

3. **Sparsity interaction**
   - Many zeros create CLR artifacts
   - Pseudocount choice affects detection

4. **Compositional closure**
   - CLR preserves closure constraint
   - True effects may cancel in log-ratio space

## Research Questions

1. **What is the power curve for LinDA?**
   - At what effect sizes is detection reliable?
   - How does sample size affect power?
   - What is the minimum detectable effect size?

2. **How does sparsity affect power?**
   - Does high sparsity reduce sensitivity?
   - Is there a prevalence threshold for detection?
   - How does zero handling affect results?

3. **Is LinDA's conservatism appropriate?**
   - Does low sensitivity come with low FDR?
   - Is the sensitivity/FDR tradeoff optimal?
   - When should we prefer LinDA vs. other methods?

## Proposed Analysis

### Phase 1: Power Curves

1. Generate synthetic data:
   - Effect sizes: 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0 log2FC
   - Sample sizes: n = 10, 20, 30, 50, 100 per group
   - Sparsity levels: 30%, 50%, 70%, 90%

2. Run LinDA on each configuration (100 replicates)

3. Compute:
   - Power at each effect size
   - Minimum detectable effect (80% power)
   - False positive rate under null

### Phase 2: Component Analysis

1. Isolate each step's contribution:
   - Raw proportions vs CLR-transformed
   - With/without pseudocount
   - Different prevalence filters

2. Compare alternatives:
   - ALR vs CLR
   - TSS vs CLR
   - No transformation

### Phase 3: Comparison

1. Compare to other methods:
   - BB model (proportion-based)
   - Hurdle model (zero-explicit)
   - Permutation test (non-parametric)
   - DESeq2/edgeR-style (if implemented)

2. ROC/PR curves across configurations

## Data Requirements

- Synthetic data with known ground truth
- Range of configurations (sparsity, effect size, n)
- Null data for FPR validation

## Expected Outcomes

1. Power tables and guidelines for sample size planning
2. Recommendations for when to use LinDA vs. alternatives
3. Understanding of CLR transformation effects on power
4. Methods paper section on power analysis

## Code Location

Scripts: `experiments/scripts/05-linda-sensitivity/`

## References

- Zhou et al. (2022) - LinDA paper
- Mandal et al. (2015) - ANCOM
- Gloor et al. (2017) - ALDEx2 and compositional data analysis
