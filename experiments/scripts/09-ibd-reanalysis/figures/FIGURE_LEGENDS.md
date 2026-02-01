# Experiment 09: IBD Cohort Reanalysis - Figure Legends

## Overview

These figures demonstrate the application of artifact risk assessment to IBD microbiome findings, showing that many associations fall within the artifact-risk zone.

---

## Figure 1: Artifact Risk by Effect Size and Method

**File**: `fig1_artifact_risk.png/pdf`

**Description**: Comparison of robust vs at-risk findings for moderate (2.0 log2FC) and large (4.0 log2FC) effect datasets.

**Panel A - Moderate Effect**:
- Most findings fall below artifact threshold
- Higher proportion at-risk across all methods

**Panel B - Large Effect**:
- More findings exceed artifact threshold
- Better separation of robust findings

**Key observations**:
- Effect size determines artifact risk classification
- Moderate effects are largely indistinguishable from artifacts
- Only large effects (>3.3 log2FC) are robust

---

## Figure 2: Effect Size Distribution with Artifact Threshold

**File**: `fig2_effect_distribution.png/pdf`

**Description**: Histogram of effect sizes with the 3.3 log2FC artifact threshold marked.

**Key observations**:
- Most taxa have effect sizes below the threshold
- Threshold derived from 9.8x load variation (Stammler data)
- Effects within threshold could be artifacts

---

## Figure 3: FPR Calibration on Null IBD Data

**File**: `fig3_fpr_validation.png/pdf`

**Description**: False positive rate for each method on null data (no true effects).

**Key observations**:
- All methods maintain proper FPR control
- LinDA, Hurdle, and Permutation well-calibrated
- Confirms methods work correctly on IBD-like data

---

## Figure 4: Clinical Implications Diagram

**File**: `fig4_clinical_implications.png/pdf`

**Description**: Conceptual diagram showing IBD-specific confounders and recommended approach.

**Key elements**:
- IBD confounders: diarrhea, inflammation, disease activity
- Consequence: 10x load variation → ±3.3 log2FC artifacts
- Solution: Report effect sizes, flag small effects, use spike-ins

---

## Summary Table

| Figure | Key Finding |
|--------|-------------|
| Fig 1 | Moderate effects mostly at-risk |
| Fig 2 | Distribution shows most effects below threshold |
| Fig 3 | Methods properly calibrated |
| Fig 4 | Clinical recommendations for IBD studies |

---

## Implications for IBD Research

1. **Many published associations may be artifacts**
   - Effect sizes < 3.3 log2FC cannot be distinguished from load variation
   - IBD patients have increased load variation (diarrhea, inflammation)

2. **Robust findings exist**
   - Some taxa show very large, consistent effects
   - These are likely biologically real

3. **Recommendations**
   - Report effect sizes, not just p-values
   - Use spike-in controls when possible
   - Flag findings < 3.3 log2FC as needing validation
