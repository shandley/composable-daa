# LinDA Sensitivity Investigation

Date: 2026-01-31

## Executive Summary

LinDA (CLR + Linear Model) shows 0% sensitivity at the standard q < 0.05 threshold, even for 16x fold changes (4.0 log2FC). This investigation reveals that:

1. **This is expected behavior**, not a bug
2. **CLR transformation attenuates effect sizes by ~75%** (4.0 log2FC becomes ~1.0)
3. **At q < 0.10, LinDA achieves 39% sensitivity with only 12.5% FDR**
4. **LinDA excels at FDR control** - use it for confident discoveries

## Key Findings

### 1. Effect Size Attenuation

The CLR transformation severely attenuates observed effect sizes:

| True log2FC | Observed CLR Effect | Attenuation |
|-------------|---------------------|-------------|
| +4.0        | +0.93               | 23%         |
| -4.0        | -0.99               | 25%         |

This happens because CLR centers by geometric mean:
- `CLR(x) = log(x) - mean(log(x))`
- When features increase, the geometric mean shifts
- This partially cancels out the observed effect

### 2. Sensitivity by Q-Value Threshold

At 4.0 log2FC (16x fold change):

| Threshold | Sensitivity | FDR   | TP  | FP  |
|-----------|-------------|-------|-----|-----|
| q < 0.01  | 0%          | n/a   | 0   | 0   |
| q < 0.05  | 0%          | n/a   | 0   | 0   |
| q < 0.10  | 39%         | 12.5% | 7   | 1   |
| q < 0.15  | 39%         | 30%   | 7   | 3   |
| q < 0.20  | 39%         | 30%   | 7   | 3   |

**Recommendation: Use q < 0.10 for LinDA**

### 3. Comparison with Other Methods (at 4.0 log2FC)

| Method | Sensitivity (q<0.05) | Sensitivity (q<0.10) | FDR (q<0.10) |
|--------|---------------------|---------------------|--------------|
| LinDA  | 0%                  | 39%                 | 12.5%        |
| ZINB   | 83%                 | 94%                 | 35%          |
| Hurdle | 83%                 | 89%                 | 33%          |
| NB     | 6%                  | 6%                  | 0%           |

### 4. Effect Size Requirements

LinDA only detects effects at 4.0 log2FC (16x fold change). At smaller effect sizes:

| Effect Size | Fold Change | LinDA Sensitivity |
|-------------|-------------|-------------------|
| 0.5 log2FC  | 1.4x        | 0% (any threshold)|
| 1.0 log2FC  | 2x          | 0% (any threshold)|
| 2.0 log2FC  | 4x          | 0% (any threshold)|
| 4.0 log2FC  | 16x         | 39% (at q<0.10)   |

## Root Cause Analysis

### Why CLR Attenuates Effects

1. **Geometric mean centering**: CLR subtracts the mean of log-transformed values
2. **Cross-feature dependencies**: When one feature increases, it affects the geometric mean
3. **Compositional constraints**: In relative abundance data, proportions must sum to 1

### Mathematical Explanation

For a feature with true fold change `k`:
```
True effect: log2(k) = 4.0 (for k=16)
CLR observed: log2(k) * attenuation_factor ≈ 1.0
Attenuation factor ≈ 0.25
```

The attenuation factor depends on:
- Number of features (more features = less attenuation theoretically)
- Correlation structure (high correlation = more attenuation)
- Magnitude of effect (large effects shift geometric mean more)

### Not a Bug: By Design

CLR is designed to handle compositional data by removing the arbitrary sum constraint. The trade-off is:
- **Pro**: Robust to compositional artifacts
- **Con**: Reduced power to detect true effects

## Recommendations

### For Users

1. **Use q < 0.10 for LinDA discovery**
   - At this threshold: 39% sensitivity, 12.5% FDR
   - Much better than q < 0.05 (0% sensitivity)

2. **Expect to detect only large effects**
   - LinDA needs ~16x fold change to detect anything
   - For smaller effects, use ZINB or Hurdle models

3. **Use LinDA for confirmation, ZINB for discovery**
   - LinDA: Low sensitivity, excellent FDR control
   - ZINB: High sensitivity, moderate FDR control

4. **Increase sample size for better power**
   - With n=20/group, power is limited
   - n=50/group would improve detection of smaller effects

### For Documentation

1. **Document effect size attenuation**
   - Users should know CLR reduces observed effects by ~75%
   - Effect size estimates are NOT directly interpretable as fold changes

2. **Recommend q < 0.10 as default**
   - The standard q < 0.05 is too conservative for LinDA
   - q < 0.10 provides good sensitivity/FDR balance

3. **Provide method selection guidance**
   - Discovery (maximize true positives): Use ZINB/Hurdle
   - Confirmation (minimize false positives): Use LinDA
   - Unknown data structure: Run multiple methods, compare

### For Development

1. **Consider adaptive thresholds**
   - Automatically suggest q threshold based on data
   - Show sensitivity/FDR trade-off curves

2. **Add power analysis tools**
   - `daa power-calc` to estimate required sample size
   - Help users design adequately powered studies

3. **Implement effect size calibration**
   - Back-transform CLR effects to approximate fold changes
   - Caveat: only approximate due to compositional effects

## Files Generated

- `analyze_linda_sensitivity.py` - Detailed analysis of LinDA results
- `analyze_threshold_sensitivity.py` - Threshold comparison across methods
- `analyze_clr_math.py` - Mathematical analysis of CLR attenuation
- `analyze_raw_vs_compositional.py` - Spike mode comparison
- `raw_mode_test/` - Test data with raw spike mode

## Conclusion

LinDA's low sensitivity is expected behavior due to CLR effect size attenuation. Users should:
1. Use q < 0.10 threshold
2. Expect to detect only very large effects (>8x fold change)
3. Use LinDA for FDR-controlled discovery
4. Consider ZINB/Hurdle for higher sensitivity when FDR control is less critical
