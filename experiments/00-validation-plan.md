# Validation Plan for Publication

## Overview

This document outlines the testing and validation required to support a methods paper submission and build user confidence in the composable-daa toolkit.

## Part 1: Comparison to Established Methods

### 1.1 Methods to Compare

| Method | R Package | Model Type | Our Equivalent |
|--------|-----------|------------|----------------|
| DESeq2 | DESeq2 | NB GLM + shrinkage | `model_nb` + `shrink_lfc` |
| edgeR | edgeR | NB GLM (TMM) | `norm_tmm` + `model_nb` |
| LinDA | MicrobiomeStat | CLR + LM | `norm_clr` + `model_lm` |
| ANCOM-BC | ANCOMBC | Log-linear + bias correction | Not directly implemented |
| ALDEx2 | ALDEx2 | CLR + Welch's t | `norm_clr` + Wald test |
| Maaslin2 | Maaslin2 | LM with transformations | Various pipelines |

### 1.2 Comparison Approach

**Synthetic Data Comparison:**
1. Generate data with known ground truth
2. Run all methods
3. Compare: sensitivity, FDR, effect size accuracy, ranking

**Real Data Comparison:**
1. Run same real datasets through all methods
2. Compare: overlap in significant features, effect size correlation, ranking concordance

### 1.3 Deliverables

- [ ] R script to run DESeq2, edgeR, LinDA, ALDEx2 on test data
- [ ] Comparison tables and figures
- [ ] Statistical tests for method agreement

---

## Part 2: Statistical Validation

### 2.1 FPR Calibration (Type I Error)

**Test**: Under null hypothesis (no true effects), proportion of p < 0.05 should be ~5%

| Model | Current Status | Needed |
|-------|----------------|--------|
| LM (LinDA) | ✅ 0% on null | Verify across more scenarios |
| LMM | ❓ Not tested | Test with correlated null data |
| NB | ❓ Not tested | Test on count-scale null data |
| ZINB | ❓ Not tested | Test on zero-inflated null |
| BB | ✅ 3% after fix | Verify stable across scenarios |
| Hurdle | ✅ 2% | Verify across scenarios |
| Permutation | ✅ 0% | Gold standard reference |

**Scenarios to test:**
- [ ] Varying sample sizes: n = 10, 20, 50, 100 per group
- [ ] Varying sparsity: 30%, 50%, 70%, 90%
- [ ] Varying library sizes: balanced vs imbalanced
- [ ] With/without batch effects

### 2.2 Power Analysis

**Test**: Proportion of true effects detected at various effect sizes

**Deliverables:**
- [ ] Power curves: Power vs effect size for each method
- [ ] Power tables: Sample size needed for 80% power
- [ ] Comparison across sparsity levels

### 2.3 Confidence Interval Coverage

**Test**: Do 95% CIs contain true value 95% of the time?

**Approach:**
1. Generate 1000 datasets with known effect = 1.0 log2FC
2. Compute 95% CI for each
3. Check coverage: should be ~95%

**Deliverables:**
- [ ] Coverage rates for each method
- [ ] Identification of under/over-coverage issues

### 2.4 Effect Size Calibration

**Test**: Are effect estimates unbiased?

**Metrics:**
- Bias: E[estimate] - truth
- RMSE: sqrt(E[(estimate - truth)²])
- Correlation: cor(estimate, truth)

**Deliverables:**
- [ ] Calibration plots (estimated vs true)
- [ ] Bias tables by prevalence tier
- [ ] RMSE comparison across methods

---

## Part 3: Real Data Validation

### 3.1 Datasets with Ground Truth

| Dataset | Ground Truth | Analysis |
|---------|--------------|----------|
| Stammler spike-in | Known spike taxa | Should detect spike-ins as differential if comparing spike levels |
| Synthetic spike-in | Injected effects | Validate spike-in framework |

### 3.2 Datasets with Expected Results

| Dataset | Expected Findings | Analysis |
|---------|-------------------|----------|
| Ravel BV | Lactobacillus ↓, Gardnerella ↑ | Confirm known associations |
| HMP body sites | Site-specific taxa | Confirm obvious differences |
| IBD studies | Known IBD taxa | Compare to published results |

### 3.3 Reproducibility Checks

**Test**: Does our implementation match R package results?

**Approach:**
1. Run LinDA R package on Ravel data
2. Run our LinDA pipeline on same data
3. Compare: coefficients, SEs, p-values, significant features

**Acceptance criteria:**
- Coefficients within 1% (allowing for numerical precision)
- Same features significant at q < 0.05
- Ranking correlation > 0.99

---

## Part 4: Robustness Testing

### 4.1 Edge Cases

| Scenario | Expected Behavior | Test |
|----------|-------------------|------|
| Feature absent in one group | Flag or handle gracefully | ✅ |
| Single sample per group | Error with message | ❓ |
| All zeros for a feature | Skip or flag | ❓ |
| Extreme outlier sample | Robust to influence | ❓ |
| Collinear covariates | Error with message | ❓ |

### 4.2 Numerical Stability

| Scenario | Expected Behavior | Test |
|----------|-------------------|------|
| Very small counts | No underflow | ❓ |
| Very large counts | No overflow | ❓ |
| Near-singular design matrix | Warning + fallback | ❓ |
| Non-convergence | Clear warning | ✅ |

### 4.3 Scalability

| Metric | Target | Test |
|--------|--------|------|
| Features | 10,000+ | ❓ |
| Samples | 1,000+ | ❓ |
| Memory usage | < 8GB for typical data | ❓ |
| Runtime | < 1 min for typical data | ❓ |

---

## Part 5: Documentation Validation

### 5.1 User Guidance

- [ ] Decision tree: which method to use when
- [ ] Parameter selection guide
- [ ] Interpretation guide for results
- [ ] Common pitfalls and solutions

### 5.2 Worked Examples

- [ ] Basic two-group comparison
- [ ] Longitudinal analysis with LMM
- [ ] Multi-factor design
- [ ] Spike-in validation workflow

---

## Priority Order

### Phase 1: Critical for Submission
1. FPR calibration for all models (Part 2.1)
2. Comparison to DESeq2/edgeR/LinDA (Part 1)
3. Power analysis (Part 2.2)
4. LinDA reproducibility check (Part 3.3)

### Phase 2: Strengthens Paper
5. CI coverage (Part 2.3)
6. Effect size calibration (Part 2.4)
7. Real data validation (Part 3.1, 3.2)

### Phase 3: User Confidence
8. Edge case testing (Part 4.1)
9. Scalability benchmarks (Part 4.3)
10. Documentation (Part 5)

---

## Timeline Estimate

| Phase | Tasks | Estimate |
|-------|-------|----------|
| Phase 1 | FPR, comparison, power | 2-3 weeks |
| Phase 2 | CI, calibration, real data | 2 weeks |
| Phase 3 | Edge cases, docs | 1-2 weeks |

Total: ~6 weeks for comprehensive validation

---

## Immediate Next Steps

1. **Set up R comparison environment**
   - Install DESeq2, edgeR, MicrobiomeStat, ALDEx2
   - Create wrapper scripts for standardized comparison

2. **Generate comprehensive test datasets**
   - Null datasets for FPR (multiple configurations)
   - Effect datasets for power (range of effect sizes)
   - Edge case datasets

3. **Implement automated testing framework**
   - CI/CD for regression testing
   - Automated benchmark reports
