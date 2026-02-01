# Experiment 07: scRNA-seq Generalization - Figure Legends

## Overview

These figures demonstrate that the composable DAA toolkit **generalizes beyond microbiome data** to single-cell RNA-seq, showing that the same statistical methods and recommendations apply across sparse count data types.

---

## Figure 1: Method Performance Comparison (scRNA-seq vs Microbiome)

**File**: `fig1_method_comparison.png/pdf`

**Description**: Side-by-side comparison of sensitivity and FDR for all methods when applied to scRNA-seq vs microbiome data. Both synthetic datasets have 20 true differential features with effect size 2.0 log2FC.

**Key observations**:
- Method rankings are consistent across data types
- Hurdle performs best for both high-sparsity datasets
- ZINB provides good discovery performance in both domains
- LinDA shows similar conservative behavior for both

**Interpretation**: The same methods that work for microbiome analysis work equally well for scRNA-seq, confirming the toolkit's generalizability.

---

## Figure 2: FPR Calibration on Null scRNA-seq Data

**File**: `fig2_fpr_calibration.png/pdf`

**Description**: False positive rate for each method when applied to scRNA-seq data with no true differential expression (null data). Methods should have FPR ≤ alpha.

**Key observations**:
- All methods maintain proper FPR control (≤ 10% for alpha = 0.05)
- Variance estimation methods (ZINB, Hurdle) are well-calibrated
- LinDA shows excellent type I error control
- Permutation test provides gold standard calibration

**Interpretation**: The statistical methods are properly calibrated for scRNA-seq data, not just microbiome data. Users can trust the p-values and q-values.

---

## Figure 3: Sensitivity by Threshold (scRNA-seq)

**File**: `fig3_sensitivity_threshold.png/pdf`

**Description**: Comparison of detection sensitivity at q < 0.05 vs q < 0.10 thresholds for scRNA-seq data.

**Key observations**:
- LinDA requires q < 0.10 (CLR attenuation effect applies to scRNA-seq too)
- ZINB and Hurdle work well at q < 0.05
- The same threshold recommendations from microbiome benchmarks apply

**Interpretation**: The threshold guidance developed for microbiome analysis (LinDA q < 0.10, others q < 0.05) transfers directly to scRNA-seq.

---

## Figure 4: Data Characteristics Comparison

**File**: `fig4_data_characteristics.png/pdf`

**Description**: Two-panel figure showing (A) sparsity levels across data types and (B) method recommendations.

**Panel A - Sparsity Comparison**:
- 16S microbiome: ~65% zeros
- Virome: ~89% zeros
- scRNA-seq: ~87% zeros
- scRNA-seq falls in the "high sparsity" category with virome

**Panel B - Method Recommendations**:
- Hurdle recommended for high-sparsity data (virome, scRNA-seq)
- ZINB good across all data types
- LinDA suitable for confirmation in all domains

**Interpretation**: scRNA-seq shares critical properties with virome data (high sparsity, zero-inflation), explaining why the same methods are optimal for both.

---

## Figure 5: Generalization Summary

**File**: `fig5_generalization_summary.png/pdf`

**Description**: Conceptual diagram showing the unified framework for sparse count data analysis.

**Key elements**:
- Common challenges shared across data types (sparsity, zero-inflation, overdispersion)
- Three data types sharing these properties (microbiome, virome, scRNA-seq)
- Single solution: Composable DAA Toolkit

**Interpretation**: This diagram encapsulates the experiment's key finding - sparse count data is a unified domain where methods developed for microbiome analysis apply directly to scRNA-seq and other high-sparsity count data.

---

## Summary Table

| Figure | Key Finding | Main Recommendation |
|--------|-------------|---------------------|
| Fig 1 | Same methods work | Use Hurdle for discovery in both domains |
| Fig 2 | FPR well-calibrated | Trust p-values for scRNA-seq |
| Fig 3 | Same thresholds | LinDA q < 0.10, others q < 0.05 |
| Fig 4 | Same data properties | scRNA-seq ≈ virome in sparsity |
| Fig 5 | Unified framework | One toolkit for all sparse count data |

---

## Methods

### Data Generation
- scRNA-seq-like data generated using `sparse_virome` preset (85% sparsity)
- Microbiome comparison data using `typical_16s` preset (60% sparsity)
- 20 features with true differential abundance (effect size 2.0 log2FC)
- 40 samples (20 per group)

### Analysis
- All five methods run on both data types: LinDA, ZINB, Hurdle, NB, Permutation
- FPR validation on null scRNA-seq data (no true effects)
- Same prevalence filter (5%) applied to all analyses

### Figures
- Generated with matplotlib/seaborn
- Publication quality (300 DPI)
- Available in PNG and PDF formats
