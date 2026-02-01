# Experiment 08: HMP Gingival Analysis - Figure Legends

## Overview

These figures demonstrate that the composable DAA toolkit generalizes across body sites, from vaginal (Ravel BV study) to oral (HMP gingival) microbiome data.

---

## Figure 1: Method Performance on HMP Gingival Data

**File**: `fig1_method_summary.png/pdf`

**Description**: Number of significant taxa detected by each method at q < 0.05 and q < 0.10 thresholds on HMP gingival microbiome data.

**Key observations**:
- Method rankings are consistent with vaginal microbiome results
- LinDA shows fewer detections due to CLR attenuation
- Hurdle and ZINB show good detection rates

---

## Figure 2: Cross-Body-Site Comparison

**File**: `fig2_body_site_comparison.png/pdf`

**Description**: Comparison of relative method performance across three body sites: vaginal (Ravel), oral (HMP gingival), and gut (synthetic/IBD).

**Key observations**:
- Method rankings are consistent across all body sites
- ZINB and Hurdle perform well universally
- LinDA remains conservative across all sites
- Toolkit is body-site-agnostic

---

## Figure 3: Compositional Closure Verification

**File**: `fig3_compositional_closure.png/pdf`

**Description**: Verification that CLR compositional closure (sum of estimates = 0) applies across all body sites.

**Key observations**:
- Sum of CLR estimates is exactly 0 for all body sites
- This is a mathematical constraint, not a biological finding
- Demonstrates toolkit correctly implements CLR transformation

---

## Figure 4: Cross-Body-Site Generalization Diagram

**File**: `fig4_generalization_diagram.png/pdf`

**Description**: Conceptual diagram showing how the toolkit's validation extends across body sites.

**Key elements**:
- Three body sites: Vaginal, Oral, Gut
- Same methods, thresholds, and interpretation
- Validates body-site-agnostic design

---

## Summary Table

| Figure | Key Finding |
|--------|-------------|
| Fig 1 | Methods work on oral microbiome |
| Fig 2 | Rankings consistent across body sites |
| Fig 3 | CLR closure applies universally |
| Fig 4 | Toolkit is body-site-agnostic |

---

## Methods

### Data
- HMP V1-V3 gingival microbiome (~300 samples)
- Fetched via `daa fetch -d hmp_v13`

### Analysis
- All methods run with prevalence filter 10%
- Comparison based on sex (if available in metadata)
- Same pipeline configurations as other experiments
