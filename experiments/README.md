# Experimental Opportunities

This directory documents potential research studies that emerged during development of the composable-daa toolkit. These serve as case studies demonstrating the toolkit's capabilities and could form the basis of methods paper validation or standalone publications.

## Overview

| Experiment | Status | Priority | Data Available |
|------------|--------|----------|----------------|
| [BV Compositional Analysis](01-bv-compositional-analysis.md) | Exploratory complete | High | Yes (Ravel 2011) |
| [Spike-in Load Estimation](02-spikein-load-estimation.md) | Proof-of-concept | High | Yes (Stammler 2016) |
| [Compositional Artifact Audit](03-compositional-artifact-audit.md) | Conceptual | Medium | Requires meta-analysis |

## Key Findings So Far

### 1. Total Load Variation is Massive
Using Stammler spike-in data, we demonstrated that samples vary **9.8x in total bacterial load**. This creates potential log2FC artifacts of **±3.3** - comparable to or larger than most reported effect sizes in the literature.

### 2. BV Signature is Compositionally Constrained
Analysis of Ravel 2011 BV data shows:
- 47/55 taxa significantly associated with BV (q<0.05)
- Sum of log2FC = 0.17 (essentially zero due to CLR closure)
- Classic findings (Lactobacillus decrease, Gardnerella increase) cannot be distinguished from total load effects

### 3. Statistical Robustness ≠ Biological Interpretability
Permutation tests show 0% FPR (well-calibrated statistics), but compositional constraints mean we cannot determine:
- Whether Lactobacillus truly decreases in absolute terms
- Whether BV-associated bacteria truly increase
- Or whether total load changes create the appearance of both

## Relationship to Methods Paper

These experiments could serve as:
1. **Validation case studies** demonstrating toolkit capabilities
2. **Motivating examples** for why compositional-aware analysis matters
3. **Standalone findings** worthy of separate publication

The BV analysis is particularly compelling given the direct collaboration with the Ravel lab.

## Next Steps

1. ~~Complete main package architecture~~ DONE (254 tests passing)
2. Formalize these analyses with reproducible scripts
3. Generate publication-quality figures
4. Consider integration with absolute quantification data (qPCR, flow cytometry)
5. Write methods paper demonstrating toolkit on these case studies
