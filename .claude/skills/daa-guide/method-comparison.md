# Method Comparison Benchmarks

Comprehensive benchmarking results from synthetic data with known ground truth.

## Test Conditions

- **Data**: Typical 16S profile (65% sparsity)
- **Samples**: n=20 per group (40 total)
- **Effect sizes tested**: 0.5, 1.0, 2.0, 4.0 log2FC (1.4x to 16x fold changes)
- **Differential features**: 18 (10 up, 8 down)

## Sensitivity at q < 0.05

| Method | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA | 0% | 0% | 0% | 0% |
| NB | 0% | 0% | 0% | 6% |
| ZINB | 11% | 11% | 58% | 83% |
| Hurdle | 0% | 0% | 26% | 83% |

**Key finding**: LinDA detects nothing at q < 0.05, even at 16x fold changes.

## LinDA at q < 0.10 (Recommended)

| Effect Size | Sensitivity | FDR | True Positives | False Positives |
|-------------|-------------|-----|----------------|-----------------|
| 0.5 log2FC | 0% | n/a | 0 | 0 |
| 1.0 log2FC | 0% | n/a | 0 | 0 |
| 2.0 log2FC | 0% | n/a | 0 | 0 |
| 4.0 log2FC | 39% | 12.5% | 7 | 1 |

**Key finding**: At q < 0.10, LinDA achieves 39% sensitivity with excellent 12.5% FDR.

## FDR at q < 0.05

| Method | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA | n/a | n/a | n/a | n/a |
| NB | n/a | n/a | n/a | 0% |
| ZINB | 50% | 78% | 45% | 29% |
| Hurdle | 100% | 100% | 17% | 25% |

**Key finding**: ZINB and Hurdle have moderate FDR at large effects but poor FDR at small effects.

## False Positive Rate (Null Data)

Testing on data with NO true effects:

| Method | False Positives | FPR | Notes |
|--------|-----------------|-----|-------|
| LinDA | 0/200 | 0% | Conservative |
| Hurdle | 4/200 | 2.0% | Well-calibrated |
| Permutation | 0/200 | 0% | Conservative |
| NB | 0/200 | 0% | Conservative |
| ZINB | 0/200 | 0% | Conservative |
| LMM | 1/100 | 1% | Well-calibrated |

**Key finding**: All methods show proper type I error control.

## Why LinDA Has Low Sensitivity

### The CLR Attenuation Problem

CLR transformation attenuates effect sizes by ~75%:

| True log2FC | Observed CLR Effect | Attenuation |
|-------------|---------------------|-------------|
| +4.0 | +0.93 | 23% |
| -4.0 | -0.99 | 25% |

### Mathematical Explanation

CLR formula: `CLR(x) = log(x) - mean(log(x))`

When a feature increases:
1. Its log value increases
2. The geometric mean (mean of logs) also increases
3. The CLR value increases less than the raw log change
4. Net effect: ~25% of true effect observed

### This Is By Design

CLR is designed to handle compositional data constraints:
- **Pro**: Robust to compositional artifacts
- **Con**: Reduced power to detect true effects

## Method Selection Summary

### For Discovery (Maximize True Positives)
- **Use**: ZINB or Hurdle at q < 0.05
- **Expected**: 83% sensitivity, 25-29% FDR
- **Trade-off**: More false positives, but find most true effects

### For Confirmation (Minimize False Positives)
- **Use**: LinDA at q < 0.10
- **Expected**: 39% sensitivity, 12.5% FDR
- **Trade-off**: Miss many true effects, but high confidence in findings

### When Unsure
- Run multiple methods and compare overlap
- Features significant in both LinDA and ZINB/Hurdle are high confidence
- Features only significant in ZINB/Hurdle need biological validation
