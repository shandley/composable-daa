# Benchmark Findings Summary

Date: 2025-01-31

## Critical Bug Fix: Beta-Binomial FPR Inflation

### Problem Discovered
- Beta-binomial model had **98.5% false positive rate** on null data
- T-test on same data: 4.5% FPR (expected)
- LinDA: 0% FPR (well-calibrated)
- Permutation test: 0% FPR (well-calibrated)

### Root Cause
The theoretical beta-binomial variance formula underestimates variance for microbiome data:
- Standard errors were 0.0003-0.006 (should be 0.1-0.5)
- With large library sizes (n~10,000), model-based SEs scale as 1/sqrt(n)
- When overdispersion ρ is small, variance is severely underestimated
- Resulted in massive z-statistics (100-40,000) even for null effects

### Solution Implemented
Replaced model-based SEs with empirical variance-based SEs (similar to t-test):
- Uses observed logit-scale variance across samples
- Sets minimum variance floor (0.5) based on typical microbiome variability
- Location: `src/model/bb.rs:compute_std_errors_bb()`

### Result
- **Before fix**: 98.5% FPR on null data (197/200 false positives)
- **After fix**: 3.0% FPR on null data (6/200 false positives)

This is within acceptable range for α=0.05.

---

## Critical Finding: Beta-Binomial FDR Is Fundamental, Not a Bug

### Problem
Despite fixing the FPR issue on null data (98.5% → 3%), BB still shows ~85% FDR when true effects exist.

### Investigation Findings

**The "false positives" have REAL proportion changes in the data.**

When we spike UP some features by 16x, the proportions of OTHER features must decrease (compositional closure). BB correctly detects these decreases:

| Feature | Control Prop | Treatment Prop | Actual log2FC | BB Detects? |
|---------|--------------|----------------|---------------|-------------|
| Feature_0061 | 0.00165 | 0.00017 | -3.24 | YES (q<1e-8) |
| Feature_0107 | 0.00148 | 0.00019 | -2.94 | YES (q<1e-8) |
| Feature_0048 | 0.00484 | 0.00085 | -2.50 | YES (q<1e-8) |

These are **real** proportion changes (up to 10x decrease), not statistical artifacts.

### Root Cause: Proportion-Based Inference on Compositional Data

BB models proportions directly: Y_i / library_size. For compositional data:
- Proportions must sum to 1 (closed-sum constraint)
- Increasing some features' proportions FORCES others to decrease
- BB correctly detects these forced decreases
- This is expected behavior, not a bug

### Evidence

**Direction of false positives matches compositional prediction:**
- True positives: 10 spiked UP, 8 spiked DOWN (net UP)
- False positives: 21 positive, 57 negative (opposite to TP direction)

**LinDA (CLR-based) handles this correctly:**
- Same features that BB detects with q<1e-8
- LinDA reports q > 0.09 (not significant)
- CLR transformation corrects for compositional artifacts

### Conclusion

BB is fundamentally unsuitable for standard differential abundance analysis due to compositional artifacts. **This is not fixable** without changing the model's target (from proportions to something else).

### Updated Recommendation

**DO NOT USE BB for standard DAA.** Use instead:
- LinDA (q<0.10) for compositionally-aware analysis
- ZINB for count-based analysis with zero-inflation
- Hurdle for presence/absence + abundance modeling

See `benchmark_results/power_analysis/BB_FDR_INVESTIGATION.md` for full details.

---

## Model Performance Comparison

### On Synthetic Data (effect size = 1.0 log2FC, n=20 per group)

| Dataset | Model | True Diff | Significant | TP | FP | Sensitivity | FDR |
|---------|-------|-----------|-------------|----|----|-------------|-----|
| typical_16s (64% sparse) | LinDA | 15 | 0 | 0 | 0 | 0% | - |
| typical_16s | BB (fixed) | 15 | 80 | 9 | 71 | 60% | 89% |
| typical_16s | Hurdle | 15 | 3 | 1 | 2 | 7% | 67% |
| sparse_virome (89% sparse) | LinDA | 10 | 0 | 0 | 0 | 0% | - |
| sparse_virome | BB (fixed) | 10 | 81 | 4 | 77 | 40% | 95% |
| group_specific | LinDA | 20 | 11 | 9 | 2 | 45% | 18% |
| group_specific | BB (fixed) | 20 | 64 | 14 | 50 | 70% | 78% |

### Key Observations

1. **LinDA is conservative** - 0% sensitivity on typical_16s and sparse_virome, but good FDR control when it does detect (18% on group_specific)

2. **BB (fixed) is sensitive but has high FDR** - Detects more true effects but also many false positives. The high FDR suggests the synthetic data generation may create spurious correlations.

3. **Hurdle underperforms** - Very low sensitivity, may need investigation

4. **Effect size of 1.0 log2FC is challenging** - With n=20 per group, power is limited for 2-fold changes

### Null Data FPR Summary

**Cross-sectional null data (n=40, 50% control/treatment, no effect)**

| Model | False Positives | FPR | Notes |
|-------|-----------------|-----|-------|
| LinDA | 0/200 | 0% | Conservative |
| BB (before fix) | 197/200 | 98.5% | Critical bug |
| BB (after fix) | 6/200 | 3.0% | ✓ Fixed |
| Hurdle | 4/200 | 2.0% | Conservative |
| Permutation | 0/200 | 0% | Conservative |
| NB | 0/200 | 0% | Conservative |
| ZINB | 0/200 | 0% | Conservative |

**Longitudinal null data (20 subjects x 3 timepoints, ICC~0.3, no group effect)**

| Model | Features | p < 0.01 | p < 0.05 | p < 0.10 | Notes |
|-------|----------|----------|----------|----------|-------|
| LMM | 100 | 0% | 0% | 1% | Slightly conservative, proper FPR control |

**Interpretation**: All models show proper type I error control after the BB fix. The models tend to be conservative (fewer false positives than nominal α), which is preferable to anti-conservative behavior.

---

## New Experiment Opportunities Identified

### 1. Variance Estimation Deep-Dive
**Question**: Why does the theoretical BB variance fail so dramatically?

**Proposed analysis**:
- Compare model-based vs empirical variance across different data structures
- Investigate when overdispersion estimation fails
- Test sandwich estimators vs empirical approach

**Potential publication**: "Variance estimation for compositional data: When theory fails"

### 2. LinDA Sensitivity Analysis
**Question**: Why does LinDA have 0% sensitivity on typical microbiome data?

**Proposed analysis**:
- Profile the CLR transformation effect on signal detection
- Compare to untransformed proportion tests
- Test with larger effect sizes (2, 3, 4 log2FC)

**Potential publication**: "Power analysis of compositional data methods"

### 3. Effect Size Recovery Accuracy
**Question**: How well do different methods recover true effect sizes?

**Proposed analysis**:
- Compare estimated log2FC to true log2FC
- Assess shrinkage effects
- Test calibration across prevalence tiers

### 4. Confounding by Library Size
**Question**: How do methods handle systematic library size differences?

**Proposed analysis**:
- Use confounded synthetic data (3x library size difference)
- Compare methods with/without normalization
- Test spike-in normalization effectiveness

### 5. Mixed Model Validation
**Question**: Does LMM properly handle longitudinal/repeated measures?

**Proposed analysis**:
- Generate synthetic repeated measures data
- Compare LMM to naive LM
- Test random effect specification impact

---

## Power Analysis Results

Comprehensive power analysis across effect sizes 0.5, 1.0, 2.0, 4.0 log2FC (1.4x to 16x fold changes). Data: typical 16S profile (65% sparsity, n=20 per group).

### Sensitivity at q < 0.05

| Model  | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA  | 0%         | 0%         | 0%         | 0%         |
| NB     | 0%         | 0%         | 0%         | 6%         |
| ZINB   | 11%        | 11%        | 58%        | 83%        |
| BB     | 61%        | 58%        | 95%        | 100%       |
| Hurdle | 0%         | 0%         | 26%        | 83%        |

### FDR at q < 0.05

| Model  | 0.5 log2FC | 1.0 log2FC | 2.0 log2FC | 4.0 log2FC |
|--------|------------|------------|------------|------------|
| LinDA  | n/a        | n/a        | n/a        | n/a        |
| NB     | n/a        | n/a        | n/a        | 0%         |
| ZINB   | 50%        | 78%        | 45%        | 29%        |
| BB     | **87%**    | **86%**    | **81%**    | **81%**    |
| Hurdle | 100%       | 100%       | 17%        | 25%        |

### Key Insights

1. **LinDA detects nothing at q<0.05** due to CLR effect size attenuation; at q<0.10 achieves 39% sensitivity with only 12.5% FDR
2. **BB has ~85% FDR at ALL effect sizes** - the FPR fix didn't solve broader FDR issues
3. **ZINB best for large effects** - 83% sensitivity, 29% FDR at 4.0 log2FC
4. **Hurdle excels at large effects** - 83% sensitivity, 25% FDR with good FDR control

See `benchmark_results/power_analysis/POWER_ANALYSIS_RESULTS.md` for full details.

---

## Recommendations

### For Users
1. **Use ZINB or Hurdle** for detecting large effects (>4x fold change)
2. **Use LinDA at q < 0.10** for conservative inference with excellent FDR control
3. **Use permutation tests** when distributional assumptions are uncertain
4. **AVOID BB model** - 85% FDR makes discoveries unreliable

### For Development
1. **Investigate BB FDR issue** - root cause differs from FPR issue (may be compositional artifacts)
2. **Increase sample size recommendations** - n=20 insufficient for detecting 2-fold changes
3. **Implement power calculation tools** to help users design adequately powered studies
4. **Consider adaptive thresholds** - q<0.10 may be appropriate for exploratory analysis

---

## Files Generated

### Synthetic Data
- `benchmark_results/typical_16s/` - Typical 16S synthetic data (64% sparse)
- `benchmark_results/sparse_virome/` - Sparse virome synthetic data (89% sparse)
- `benchmark_results/extreme_sparse/` - Extreme sparsity edge case
- `benchmark_results/confounded/` - Library size confounded data
- `benchmark_results/group_specific/` - Presence/absence effects
- `benchmark_results/stammler/` - Stammler spike-in dataset (Zenodo)

### FPR Validation Data
- `benchmark_results/lmm_null_data/` - Longitudinal null data (20 subjects x 3 timepoints, ICC~0.3)
  - `counts.tsv` - 100 features x 60 samples with within-subject correlation
  - `metadata.tsv` - Subject, group, time columns

### FPR Validation Results
- `benchmark_results/null_nb_results.tsv` - NB model on null data (0% FPR)
- `benchmark_results/null_zinb_results.tsv` - ZINB model on null data (0% FPR)
- `benchmark_results/null_lmm_results.tsv` - LMM model on longitudinal null data (0% FPR at α=0.05)

### Pipeline Configurations
- `benchmark_results/nb_pipeline.yaml` - Negative binomial FPR test pipeline
- `benchmark_results/zinb_pipeline.yaml` - Zero-inflated NB FPR test pipeline
- `benchmark_results/lmm_pipeline.yaml` - Linear mixed model FPR test pipeline

### Power Analysis
- `benchmark_results/power_analysis/` - Power analysis directory
  - `effect_0.5/`, `effect_1.0/`, `effect_2.0/`, `effect_4.0/` - Data at each effect size
  - `*_pipeline.yaml` - Pipeline configurations for each model
  - `POWER_ANALYSIS_RESULTS.md` - Comprehensive power analysis documentation
  - `analyze_power.py`, `analyze_power_multi.py` - Analysis scripts
