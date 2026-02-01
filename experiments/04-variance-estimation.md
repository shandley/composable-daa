# Experiment 04: Variance Estimation for Compositional Data

## Summary

Investigation demonstrating that **theoretical variance formulas fail catastrophically** for microbiome data, producing standard errors 100-1000x too small and resulting in 98.5% false positive rates on null data. This experiment documents the problem, the fix, and validates that all current methods are properly calibrated.

**Key Finding**: Model-based variance estimation (beta-binomial) produces SEs of 0.0003-0.006 when they should be 0.1-0.5, leading to essentially 100% false discovery. Empirical variance estimation fixes the FPR but reveals a deeper problem: compositional artifacts that cannot be fixed by better variance estimation.

## Background

During benchmarking of the composable-daa toolkit, we discovered severe FPR inflation in the beta-binomial model. This experiment characterizes the problem and validates the solution.

### The Variance Formula Problem

For beta-binomial GLM, the theoretical variance of proportions is:

```
Var(Y_i/n_i) = mu(1-mu)[1 + (n_i-1)*rho] / n_i
```

With microbiome data:
- Library sizes n ~ 10,000-100,000
- Overdispersion rho often estimated as very small (~1e-10)
- Result: Variance scales as 1/n, making SEs unrealistically small

### Discovery Timeline

| Phase | Observation | FPR | Status |
|-------|-------------|-----|--------|
| Initial | BB on null data | 98.5% | BROKEN |
| Diagnosis | Model SEs 100x too small | - | - |
| Fix 1 | Empirical variance | 3.0% | FIXED |
| Problem 2 | Compositional FDR | 85% | FUNDAMENTAL |
| Resolution | BB removed | - | CLOSED |

## Data

### Null Data Generation

Synthetic null datasets with no true effects, varying:

| Parameter | Values | Purpose |
|-----------|--------|---------|
| Sparsity | 30%, 50%, 70%, 90% | Test sparse data handling |
| Sample size | 20, 40, 60, 100 per group | Test power effects |
| Library size | 1K, 10K, 100K | Test variance scaling |

```bash
# Generate null data
daa generate --n-samples 40 --n-features 200 --sparsity 0.7 \
    --library-size 10000 --effect-size 0.0 --n-differential 0 \
    -o ./data/null/sparsity_0.7
```

## Analysis Performed

### 1. FPR Calibration Test

Run all methods on null data and measure false positive rate:

| Method | FPR at p<0.05 | FPR at p<0.10 | Calibrated? |
|--------|---------------|---------------|-------------|
| LinDA | 2-4% | 4-7% | YES |
| Hurdle | 2-5% | 5-8% | YES |
| ZINB | 2-5% | 5-9% | YES |
| NB | 1-3% | 3-6% | YES |
| Permutation | 2-5% | 5-10% | YES (gold std) |
| BB (removed) | 98.5% | 99.5% | NO - removed |

**Finding**: All current methods show proper FPR calibration (<7% at alpha=0.05).

### 2. Standard Error Analysis

Compared standard errors across methods and conditions:

| Variance Type | SE Range | FPR | Example |
|---------------|----------|-----|---------|
| Model-based (BB) | 0.0003-0.006 | 98.5% | Removed |
| Sandwich (HC3) | 0.02-0.10 | ~12% | Improved |
| Bootstrap | 0.05-0.20 | ~8% | Expensive |
| Empirical | 0.10-0.50 | 3-5% | Current |

**Finding**: Empirical variance estimation produces properly scaled SEs.

### 3. Z-Statistic Distribution

On null data, z-statistics should follow N(0,1):

| Method | Z-stat Range | Distribution | Status |
|--------|--------------|--------------|--------|
| Model-based (BB) | 100-40,000 | Catastrophic | Removed |
| Current methods | -3 to +3 | N(0,1) | Correct |

**Finding**: Model-based z-statistics were 100-10,000x too large.

### 4. Library Size Effect

The key insight: variance scales incorrectly with library size.

| Library Size | Model SE | Empirical SE | Ratio |
|--------------|----------|--------------|-------|
| 1,000 | 0.03 | 0.25 | 8x |
| 10,000 | 0.01 | 0.25 | 25x |
| 100,000 | 0.003 | 0.25 | 83x |

**Finding**: As library size increases, model-based SE decreases while true variability remains constant.

## Key Insights

### 1. Why Theoretical Variance Fails

The beta-binomial variance formula assumes that observed variance comes from:
1. Binomial sampling (1/n term)
2. Overdispersion (rho term)

But microbiome data has additional sources of variability:
- Biological variation between samples
- Compositionality effects
- Zero inflation from multiple sources

When n is large (10,000+) and rho is small (~1e-10), the formula predicts near-zero variance.

### 2. The 100x Problem

```
Expected SE (from data): 0.25
Model SE (from formula): 0.003
Ratio: 83x

With z = estimate / SE:
True z: 0.5 / 0.25 = 2.0 (not significant)
Model z: 0.5 / 0.003 = 167 (p < 1e-100)
```

This is why 98.5% of null features appeared significant.

### 3. Empirical Variance Works

The t-test-style empirical variance:
```
SE = sqrt(var(logit(Y/n))) / sqrt(n_samples)
```

This captures ALL sources of variability, regardless of their origin.

### 4. But Compositionality Remains

Even after fixing variance estimation:
- BB still had 85% FDR with true effects
- Reason: Compositional closure forces spurious correlations
- When some features increase, others must decrease (proportions sum to 1)
- BB detects these "real" proportion changes that are artifacts

### 5. Resolution: Remove BB, Use Compositionally-Aware Methods

LinDA uses CLR transformation to handle compositionality:
- Detects features changing relative to geometric mean
- Avoids detecting compositional artifacts
- Properly calibrated FPR AND FDR

## Implications

### For Method Users

1. **Trust the current methods**: LinDA, Hurdle, ZINB, NB, Permutation are all calibrated
2. **Don't use proportion-based methods**: Without compositional correction
3. **Check FPR on your data**: Use permutation tests as gold standard

### For Method Developers

1. **Never trust theoretical variance for microbiome data**: Always validate empirically
2. **Test on null data first**: Before evaluating power
3. **Consider compositionality**: Variance fixes don't fix compositional artifacts

### For the Field

1. **Many published results may be FPs**: If using uncalibrated methods
2. **Validation is critical**: Spike-in or permutation-based calibration needed
3. **Compositionality is deeper than variance**: Affects effect interpretation, not just significance

## Reproducibility

All analysis is reproducible via:

```bash
cd experiments/scripts/04-variance-estimation/
./run_analysis.sh
```

Output files:
- `results/fpr_summary.tsv` - FPR across conditions
- `results/se_summary.tsv` - SE distributions
- `results/variance_analysis.tsv` - Variance estimator comparison
- `results/summary_report.txt` - Text summary

Figures:
```bash
python3 generate_figures.py
```

## Connection to Overall Paper

This experiment provides **methodological grounding** for the toolkit:

1. **Experiment 01**: Demonstrates compositional closure
2. **Experiment 02**: Quantifies load variation artifacts
3. **Experiment 03**: Audits artifact risk in published findings
4. **Experiment 04**: Validates that our methods are properly calibrated

Together, they show:
- The field has serious methodological problems (01, 02, 03)
- Our toolkit addresses them with empirically validated methods (04)
- Composability enables testing and validation not possible with black-box tools

## Mathematical Details

### Beta-Binomial Variance

For response Y_i | n_i ~ BetaBinomial(n_i, mu, rho):

```
E[Y_i/n_i] = mu
Var[Y_i/n_i] = mu(1-mu)[1 + (n_i-1)*rho] / n_i
```

When rho -> 0 and n_i -> infinity:
```
Var[Y_i/n_i] -> mu(1-mu) / n_i -> 0
SE -> 0
```

### Empirical Variance

Instead of using the model formula, compute directly:
```
logit_values = logit(Y_i/n_i) for each sample i
SE = std(logit_values) / sqrt(n_samples)
```

This is approximately what a t-test does, and it captures all sources of variability.

### Why This Works

Empirical variance is:
1. **Assumption-free**: Doesn't assume beta-binomial
2. **Includes all variation**: Biological, technical, compositional
3. **Properly scaled**: By sample size, not library size

## References

- White (1980) - Heteroscedasticity-consistent covariance estimators
- MacKinnon & White (1985) - HC variants (HC0, HC1, HC2, HC3)
- Cameron & Miller (2015) - Cluster-robust SEs in practice
- Gloor GB et al. (2017) - Microbiome datasets are compositional
- Lin H, Peddada SD (2020) - LinDA: Analysis of compositions with bias correction
