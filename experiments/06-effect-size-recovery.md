# Experiment 06: Effect Size Recovery Accuracy

## Summary

Assessment of how accurately different DAA methods recover true effect sizes (log2 fold changes), including bias, variance, and calibration across prevalence tiers.

## Background

In our benchmarking, we observed that methods detect effects with varying accuracy:

```
Ground truth: 1.0 log2FC (2-fold change)
Observed estimates ranged from -3.1 to +2.6 log2FC
```

This variation arises from:
1. Sampling noise (expected)
2. Compositional effects (CLR closure)
3. Model misspecification
4. Shrinkage (if applied)

### Why Effect Size Accuracy Matters

- **Meta-analysis**: Combining studies requires comparable effect sizes
- **Biological interpretation**: "2-fold increase" has meaning; "significant" does not
- **Power analysis**: Sample size planning requires expected effect sizes
- **Replication**: Effect size estimates should be reproducible

## Research Questions

1. **How biased are effect size estimates?**
   - Is there systematic over/underestimation?
   - Does bias vary by prevalence tier?
   - How does transformation affect estimates?

2. **How variable are estimates?**
   - What is the distribution of estimates around truth?
   - Are confidence intervals well-calibrated?
   - Does shrinkage help?

3. **Are estimates comparable across methods?**
   - Do CLR-based and proportion-based estimates correlate?
   - Can we convert between scales?
   - Which scale is most interpretable?

## Proposed Analysis

### Phase 1: Estimation Accuracy

1. Generate synthetic data with known effects:
   - True log2FC: -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3
   - Multiple prevalence tiers
   - 100 replicates per configuration

2. For each method, compute:
   - Mean estimated effect
   - Bias = E[estimate] - truth
   - RMSE = sqrt(E[(estimate - truth)²])
   - Coverage of 95% CIs

3. Stratify by:
   - Prevalence tier (rare, low, moderate, high)
   - Effect direction (up vs down)
   - Sparsity level

### Phase 2: Shrinkage Effects

1. Compare shrunk vs unshrunk estimates:
   - DESeq2-style adaptive shrinkage
   - Normal prior shrinkage
   - No shrinkage

2. Evaluate:
   - Bias-variance tradeoff
   - Improvement in RMSE
   - Effect on ranking

### Phase 3: Cross-Method Comparison

1. Convert between scales:
   - CLR log2FC vs proportion log2FC
   - Absolute vs relative interpretation

2. Assess concordance:
   - Correlation of estimates
   - Agreement on direction
   - Agreement on magnitude

## Data Requirements

- Synthetic data with exact known effects
- Multiple prevalence tiers
- Range of effect sizes (including null)

## Expected Outcomes

1. Effect size calibration guidelines
2. Recommendations for shrinkage
3. Conversion factors between methods
4. Methods paper section on effect size interpretation

## Code Location

Scripts: `experiments/scripts/06-effect-size-recovery/`

## References

- Stephens (2017) - False discovery rates and effect size shrinkage
- Zhu et al. (2019) - Effect sizes in microbiome studies
- Weiss et al. (2017) - Normalization effects on microbiome analysis
