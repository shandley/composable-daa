# Experiment 04: Variance Estimation for Compositional Data

## Summary

Investigation of why theoretical variance formulas fail for microbiome data, and development of robust alternatives.

## Background

During benchmarking, we discovered that the beta-binomial model's theoretical variance formula produces severely underestimated standard errors for microbiome data, leading to 98.5% false positive rates on null data.

### The Problem

For beta-binomial GLM, the theoretical variance of proportions is:

```
Var(Y_i/n_i) = μ(1-μ)[1 + (n_i-1)ρ] / n_i
```

With microbiome data:
- Library sizes n ~ 10,000-100,000
- Overdispersion ρ often estimated as very small (~1e-10)
- Result: Variance scales as 1/n, making SEs unrealistically small

### Key Finding

| Approach | SE Range | FPR on Null |
|----------|----------|-------------|
| Model-based | 0.0003-0.006 | 98.5% |
| Empirical (t-test style) | 0.1-0.5 | 3.0% |

The model-based approach fails because:
1. Overdispersion estimation is unreliable for sparse features
2. The beta-binomial model may not capture all sources of biological variation
3. Large library sizes amplify any misspecification

## Research Questions

1. **When does theoretical variance fail?**
   - At what sample sizes/sparsity levels?
   - For which data distributions?
   - How does it vary by prevalence tier?

2. **What are the alternatives?**
   - Sandwich (robust) estimators
   - Bootstrap-based variance
   - Empirical Bayes shrinkage of variance
   - Permutation-based calibration

3. **Can we detect when model assumptions fail?**
   - Diagnostic tests for variance misspecification
   - Automatic switching between approaches

## Proposed Analysis

### Phase 1: Characterization

1. Generate synthetic data across parameter space:
   - Sparsity: 20%, 50%, 70%, 90%, 95%
   - Sample sizes: n = 10, 20, 50, 100 per group
   - Effect sizes: 0.5, 1.0, 2.0, 3.0 log2FC
   - Library sizes: 1K, 10K, 100K

2. Compare variance estimators:
   - Model-based (Fisher information)
   - Sandwich (HC0, HC1, HC2, HC3)
   - Bootstrap (parametric, nonparametric)
   - Empirical (sample variance)

3. Evaluate:
   - Coverage of 95% CIs
   - FPR at α = 0.05
   - Power for detecting true effects

### Phase 2: Development

1. Implement alternative variance estimators in the toolkit
2. Create diagnostic tools for detecting variance misspecification
3. Develop adaptive approach that selects best estimator

### Phase 3: Validation

1. Apply to real datasets (Ravel, HMP, etc.)
2. Compare to permutation test p-values (gold standard)
3. Assess impact on biological conclusions

## Data Requirements

- Synthetic data (generated via toolkit)
- Null datasets for FPR calibration
- Real datasets for validation

## Expected Outcomes

1. Guidelines for when to use different variance estimators
2. New variance estimation functions in the toolkit
3. Diagnostic tools for model misspecification
4. Methods paper section on robust inference for compositional data

## Code Location

Scripts: `experiments/scripts/04-variance-estimation/`

## References

- White (1980) - Heteroscedasticity-consistent covariance estimators
- MacKinnon & White (1985) - HC variants
- Cameron & Miller (2015) - Cluster-robust SEs in practice
