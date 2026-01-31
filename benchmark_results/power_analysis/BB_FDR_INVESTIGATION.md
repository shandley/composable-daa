# Beta-Binomial FDR Investigation

**Date**: 2025-01-31

## Executive Summary

The beta-binomial (BB) model shows ~85% false discovery rate across all effect sizes. Investigation reveals this is **NOT a bug** but rather a **fundamental limitation of proportion-based inference** for compositional data.

**Key Finding**: BB's "false positives" have **real proportion changes** in the data (due to compositional closure). BB correctly detects these changes, but they are artifacts of the closed-sum constraint, not independent biological differences.

## Root Cause Analysis

### 1. Compositional Closure Creates Real Proportion Changes

When we spike UP some features by 16x fold change:

| True Effect | Count Change | Proportion Change |
|-------------|--------------|-------------------|
| Spiked UP features | +16x | Increase |
| Null features | No change | **DECREASE** (compositional) |

Because proportions must sum to 1, increasing some features' proportions forces others to decrease. These are **real** proportion changes, just not independent biological changes.

### 2. BB Correctly Detects These Changes

Example false positive Feature_0061:
- Control proportion: 0.00165 (0.165%)
- Treatment proportion: 0.00017 (0.017%)
- **Actual log2 fold change: -3.24** (10x decrease in proportion!)

BB correctly identifies this 10x proportion decrease with q < 3e-9.

### 3. LinDA/CLR Avoids This Problem

The same feature with CLR transformation:
- CLR difference: -1.52 (smaller than -3.24 raw)
- LinDA q-value: 0.09 (not significant)

CLR transformation partially corrects for compositional artifacts by using log-ratios relative to the geometric mean.

## Quantitative Evidence

### Direction of False Positives

At 4.0 log2FC effect:
- True positives: 10 spiked UP, 8 spiked DOWN (net UP)
- False positives: 21 positive, **57 negative**

The false positives are predominantly negative (DOWN), exactly as expected from compositional closure when true effects are net positive.

### Magnitude of Proportion Changes

Of 78 BB "false positives":
- 77 have |actual log2FC| > 0.5 (substantial real proportion changes)
- Mean |estimate| = 1.08 (large detected effects)

These are **not statistical artifacts** - they are real proportion changes being correctly detected.

### Model Comparison

| Model | FP Direction | Mechanism |
|-------|--------------|-----------|
| BB | 21+:57- (opposite to TP) | Compositional artifacts |
| ZINB | 6+:0- (same as TP) | Different mechanism |
| Hurdle | 4+:1- (same as TP) | Different mechanism |
| LinDA | 0 FPs | CLR handles compositionality |

## Conclusions

### BB is NOT Suitable for Standard DAA

The beta-binomial model:
1. Models Y/n (counts as proportion of library size)
2. Detects **proportion changes**, not absolute abundance changes
3. Will always detect compositional artifacts as significant
4. FDR will remain high (~85%) regardless of fixes

### This is By Design, Not a Bug

BB does exactly what it's supposed to do: model and detect proportion differences. The problem is that **proportions are not the appropriate target** for differential abundance analysis due to the closed-sum constraint.

## Recommendations

### Immediate Actions

1. **Add prominent warning** to BB documentation:
   ```
   WARNING: The beta-binomial model detects proportion changes, which
   include compositional artifacts. For compositionally-aware analysis,
   use LinDA (CLR + LM) or ZINB with log-link.
   ```

2. **Update CLI help text** with warning about BB limitations

3. **Update BENCHMARK_FINDINGS.md** with explanation

### Long-term Options

1. **Deprecate BB model** for standard DAA use
2. **Keep BB for specialized use cases** where proportion changes are the target
3. **Implement compositional correction** (CLR-transform proportions before BB)
4. **Add BB-CLR hybrid** model that transforms before modeling

### User Guidance

| Use Case | Recommended Model |
|----------|-------------------|
| Standard DAA | LinDA, ZINB, or Hurdle |
| Detecting proportion changes | BB (but understand limitations) |
| High FDR tolerance | ZINB at q<0.05 |
| Conservative inference | LinDA at q<0.10 |
| Large effects only | Hurdle at q<0.05 |

## Files

- `analyze_bb_fdr.py` - FDR pattern analysis
- `compare_models_fdr.py` - Cross-model comparison
- `check_proportions.py` - Verify real proportion changes
- `check_linda_proportions.py` - Compare BB vs LinDA handling
