# Figure Legends - Variance Estimation for Compositional Data

## Figure 1: FPR Calibration Heatmap

**File**: `fig1_fpr_heatmap.png/pdf`

Three-panel heatmap showing false positive rate (FPR) across all implemented methods and data conditions on null data (no true effects).

**Panel A (Sparsity)**: FPR at different sparsity levels (30%, 50%, 70%, 90%). All methods maintain FPR near the nominal 0.05 level across sparsity conditions, indicating proper calibration even for highly sparse data.

**Panel B (Sample Size)**: FPR at different sample sizes (n = 20, 40, 60, 100 per group). Methods show consistent FPR regardless of sample size.

**Panel C (Library Size)**: FPR at different library sizes (1K, 10K, 100K reads). This is the critical test - model-based variance would fail catastrophically at high library sizes, but all current methods use empirical variance and remain calibrated.

Color scale: Green indicates calibrated FPR (< 0.07), red indicates inflated FPR (> 0.10).

---

## Figure 2: Standard Error Comparison

**File**: `fig2_se_comparison.png/pdf`

**Panel A**: Standard error estimates by library size for each method. Bars show median SE with error bars indicating min-max range. Empirical methods (LinDA, Hurdle) show consistent SEs across library sizes, while model-based methods show similar stability after fixes. Gray dashed line indicates minimum reasonable SE (~0.1) based on typical microbiome variability.

**Panel B**: The variance estimation problem illustrated. Left bar shows model-based SE from the removed beta-binomial model (0.003) - approximately 100x too small. Right bar shows empirical SE (0.25) - properly calibrated. The green shaded region indicates the expected SE range (0.1-0.5) for microbiome data.

---

## Figure 3: The Variance Problem

**File**: `fig3_variance_problem.png/pdf`

Detailed visualization of why theoretical variance formulas fail for microbiome data.

**Panel A (SE Scaling)**: Log-log plot showing how standard errors scale with library size. Model-based SE (red) decreases as 1/sqrt(n), becoming unrealistically small at high library sizes (0.003 at n=10,000). Empirical SE (green) remains relatively constant (~0.1-0.3), reflecting true biological variability.

**Panel B (FPR by Estimator)**: Bar chart comparing FPR for different variance estimators on null data:
- Model-based (BB): 98.5% FPR (catastrophic failure, model removed)
- Sandwich (HC3): ~12% FPR (improved but still inflated)
- Bootstrap: ~8% FPR (better but computationally expensive)
- Empirical: ~3.5% FPR (current approach, properly calibrated)

Red dashed line = nominal 0.05 level; orange dashed line = acceptable 0.10 threshold.

**Panel C (Z-statistic Distribution)**: Histogram of z-statistics on null data. Empirical method produces z~N(0,1) as expected (green). Model-based method produced z-statistics of 100-40,000 (shown clipped), leading to essentially all features appearing significant.

**Panel D (Fix Timeline)**: Timeline showing the discovery and resolution process:
1. Discovery: BB shows 98.5% FPR on null data
2. Diagnosis: SEs 100-1000x too small (0.003 vs 0.3)
3. Fix Attempt: Empirical variance reduces FPR to 3%
4. New Problem: FDR remains 85% (compositional artifacts)
5. Resolution: BB model removed, recommend LinDA/ZINB/Hurdle

---

## Figure 4: Method FPR Summary

**File**: `fig4_method_fpr_summary.png/pdf`

Bar chart summarizing FPR by method across all tested conditions. Bars show mean FPR with error bars indicating standard deviation across conditions. All implemented methods (LinDA, Hurdle, ZINB, NB, Permutation) show FPR within the acceptable range (3-7%). Red dashed line indicates the nominal 0.05 level. Green shaded region indicates the acceptable range (0.03-0.07).

This figure demonstrates that all current methods in the toolkit are properly calibrated for Type I error control after the variance estimation fixes.

---

## Figure 5: Summary Figure

**File**: `fig5_summary.png/pdf`

Three-panel summary suitable for graphical abstract.

**Panel A (SE Estimation)**: Bar chart showing model-based SE (0.003) vs empirical SE (0.25) on log scale. Demonstrates the 100x underestimation problem.

**Panel B (FPR Calibration)**: Before/after comparison showing FPR reduction from 98.5% to 3% after switching to empirical variance estimation.

**Panel C (Current Methods)**: Status chart showing all current methods (LinDA, Hurdle, ZINB, NB, Permutation) are properly calibrated.

**Key Message Box**: Summary stating that theoretical variance formulas severely underestimate uncertainty for microbiome data, and empirical variance estimation is the solution.

---

## Technical Notes

- All figures generated at 300 DPI for publication quality
- PDF versions available for vector graphics
- Color scheme:
  - LinDA: #9b59b6 (purple)
  - Hurdle: #3498db (blue)
  - ZINB: #27ae60 (green)
  - NB: #e67e22 (orange)
  - Permutation: #95a5a6 (gray)
  - Calibrated: #27ae60 (green)
  - Uncalibrated: #e74c3c (red)
  - Warning: #f39c12 (yellow)
- Nominal FPR level: 0.05 (alpha)
- Acceptable FPR range: 0.03-0.07
- Generated by `generate_figures.py` using matplotlib

## Key Findings

1. **Model-based variance fails**: Theoretical beta-binomial variance produces SEs 100-1000x too small for microbiome data, leading to 98.5% FPR

2. **Root cause**: Large library sizes (n~10,000-100,000) combined with small overdispersion estimates cause variance to scale as 1/n

3. **Solution**: Empirical variance estimation (similar to t-test) produces properly calibrated SEs (0.1-0.5)

4. **BB removed**: Even after FPR fix, beta-binomial had 85% FDR due to compositional artifacts (not variance-related)

5. **Current methods calibrated**: LinDA, Hurdle, ZINB, NB, and Permutation all show FPR < 5% on null data
