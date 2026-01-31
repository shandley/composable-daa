# Power Analysis Results

**Date**: 2025-01-31
**Data**: Typical 16S profile (65% sparsity, n=20 per group, 200 features, 20 differential)

## Summary

This analysis evaluates statistical power (sensitivity) and false discovery rate (FDR) control across five models at four effect sizes (0.5, 1.0, 2.0, 4.0 log2FC corresponding to 1.4x, 2x, 4x, 16x fold changes).

### Key Findings

1. **LinDA is extremely conservative** but has excellent FDR control when detections occur
2. **Beta-binomial (BB) has unacceptably high FDR** (~85%) at all effect sizes
3. **ZINB offers best balance** of sensitivity and FDR for large effects
4. **Hurdle model excels at large effects** with good FDR control
5. **NB is very conservative** - useful for high-confidence discoveries

## Power Curves (Sensitivity at q < 0.05)

| Model  | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA  | 0%         | 0%         | 0%         | 0%         |
| NB     | 0%         | 0%         | 0%         | 6%         |
| ZINB   | 11%        | 11%        | 58%        | 83%        |
| BB     | 61%        | 58%        | 95%        | 100%       |
| Hurdle | 0%         | 0%         | 26%        | 83%        |

## FDR at q < 0.05

| Model  | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA  | n/a        | n/a        | n/a        | n/a        |
| NB     | n/a        | n/a        | n/a        | 0%         |
| ZINB   | 50%        | 78%        | 45%        | 29%        |
| BB     | 87%        | 86%        | 81%        | 81%        |
| Hurdle | 100%       | 100%       | 17%        | 25%        |

## Detailed Analysis by Model

### LinDA (CLR + Linear Model)
- **Strength**: Excellent FDR control when using relaxed threshold (q < 0.10: 39% sensitivity, 12.5% FDR at 16-fold change)
- **Weakness**: Zero detections at q < 0.05 even for large effects
- **Why**: CLR transformation attenuates effect sizes; 4.0 log2FC true effects appear as 0.5-1.5 log2FC after CLR
- **Recommendation**: Use q < 0.10 threshold OR increase sample size for moderate effects

### Negative Binomial (NB)
- **Strength**: Very conservative, minimal false positives
- **Weakness**: Extremely low power - misses most true effects
- **Why**: Standard errors may be inflated for sparse count data
- **Recommendation**: Use only when false positive cost is very high

### Zero-Inflated Negative Binomial (ZINB)
- **Strength**: Good sensitivity at large effects (83% at 4.0 log2FC)
- **Weakness**: Elevated FDR at small effects (50-78%)
- **Why**: Models both zero-inflation and overdispersion, but may overfit
- **Recommendation**: Best choice for detecting strong effects (>4x fold change)

### Beta-Binomial (BB)
- **Strength**: Highest sensitivity across all effect sizes
- **Weakness**: **CRITICAL: 81-87% FDR at all effect sizes**
- **Why**: Despite FPR fix on null data, the model produces many false positives when real effects exist (compositional artifacts, residual variance issues)
- **Recommendation**: **DO NOT USE** for inference until FDR issue is resolved

### Hurdle Model
- **Strength**: Good FDR control (17-25%) at large effects with high sensitivity (83%)
- **Weakness**: Poor power at small effects (<2x fold change)
- **Why**: Binary component filters low-signal features; count component tests survivors
- **Recommendation**: Excellent choice for detecting presence/absence shifts and large abundance changes

## Recommendations by Use Case

### Exploratory Analysis (Hypothesis Generation)
- **Use**: ZINB or Hurdle at q < 0.10
- **Expect**: ~60-95% sensitivity, ~25-40% FDR
- **Validate**: Follow up with targeted experiments

### Confirmatory Analysis (Publication)
- **Use**: LinDA at q < 0.10 or Hurdle at q < 0.05
- **Expect**: ~40-80% sensitivity, <25% FDR
- **Note**: Accept lower power for higher confidence

### High-Stakes Decisions (Clinical)
- **Use**: NB or Permutation test
- **Expect**: Low sensitivity but minimal false positives
- **Note**: May miss real effects but discoveries are reliable

## Methodological Notes

### Effect Size Attenuation
The CLR transformation (used in LinDA) preserves ratios but not absolute differences. A true 16-fold change (4.0 log2FC) may appear as only 1-2 log2FC after transformation, requiring larger sample sizes or more permissive thresholds for detection.

### Sample Size Considerations
With n=20 per group and 65% sparsity:
- 2-fold changes are essentially undetectable by conservative methods
- 4-fold changes require relaxed thresholds (q < 0.10)
- 16-fold changes are detectable by most methods

### Multiple Testing
The Benjamini-Hochberg correction becomes increasingly stringent with many borderline-significant features. This particularly affects LinDA where many features have moderate p-values.

## Files

- `effect_*/` - Synthetic data at each effect size
- `*_pipeline.yaml` - Pipeline configurations for each model
- `analyze_power.py` - Analysis script (q < 0.05)
- `analyze_power_multi.py` - Multi-threshold analysis script
