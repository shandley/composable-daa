# Experiment 03: Compositional Artifact Audit

## Summary

A systematic meta-analysis to quantify how many published microbiome differential abundance findings may be artifacts of compositional data analysis, and to develop robustness scores for existing datasets.

## Motivation

From our preliminary analyses:
- Stammler spike-in data shows **9.8x variation in total load** across samples
- This creates **±3.3 log2FC artifact potential**
- Many published effect sizes are within this range
- The BV "Lactobacillus decrease" signature sums to zero (compositional closure)

**Central question**: What fraction of published microbiome findings would survive correction for total bacterial load?

## Proposed Approach

### Phase 1: Literature Survey

1. **Identify target studies**: Select 20-50 high-impact microbiome papers across disease areas:
   - Inflammatory bowel disease (IBD)
   - Bacterial vaginosis (BV)
   - Obesity/metabolic syndrome
   - Colorectal cancer
   - Autism spectrum disorder

2. **Extract reported findings**:
   - Taxa reported as differentially abundant
   - Effect sizes (log2FC or equivalent)
   - Sample sizes and statistical methods
   - Whether absolute quantification was used

3. **Categorize by load potential**:
   - Studies with load measurements (rare)
   - Studies with spike-ins (rare)
   - Studies with relative abundance only (most)

### Phase 2: Stress Testing

For studies with available raw data (via SRA/ENA):

1. **Reanalyze with composable-daa toolkit**
2. **Run stress testing** to assess compositional sensitivity
3. **Calculate robustness scores**

```bash
daa stress \
    --counts study_counts.tsv \
    --metadata study_metadata.tsv \
    --formula "~ condition" \
    --test-coef conditiondisease \
    --spike-fractions "0.01,0.05,0.10,0.25" \
    --fold-changes "1.5,2.0,3.0,5.0"
```

### Phase 3: Meta-analysis

1. **Compare effect sizes** across studies with/without load correction
2. **Identify consistent findings** that survive stress testing
3. **Flag potentially artifactual findings** for follow-up

## Robustness Score Definition

For each taxon-disease association, calculate:

```
Robustness = Sensitivity_at_25%_spike / Sensitivity_at_1%_spike
```

- **Robustness > 0.8**: Finding likely real (survives compositional stress)
- **Robustness 0.5-0.8**: Moderate concern
- **Robustness < 0.5**: High artifact risk

## Expected Findings

Based on preliminary analysis, we predict:

1. **Many findings will show low robustness** (<0.5), especially:
   - Small effect sizes (<2 log2FC)
   - Disease states known to alter total load (IBD, BV)
   - Studies using older normalization methods (TSS, rarefaction)

2. **Some findings will be robust** (>0.8):
   - Large effect sizes (>3 log2FC)
   - Organisms with known biology (pathogens, probiotics)
   - Studies with independent validation

3. **Disease areas will differ**:
   - BV: High artifact risk (biofilm = increased load)
   - IBD: High artifact risk (inflammation = dysbiosis)
   - Obesity: Moderate risk
   - Localized infections: May be more robust

## Diagnostic Metrics

For each study, calculate:

| Metric | What It Measures |
|--------|------------------|
| Sum of log2FC | Compositional closure (should be ~0 for CLR) |
| Fraction unidirectional | Whether changes are balanced |
| Effect size distribution | Whether effects are plausible |
| Permutation FPR | Statistical calibration |
| Stress test sensitivity | Robustness to compositional artifacts |

## Therapeutic Implications

If many findings are load artifacts:

1. **Drug development** based on microbiome signatures may be misdirected
2. **Probiotic formulations** may work for wrong reasons
3. **Diagnostic biomarkers** may not replicate
4. **Personalized medicine** based on microbiome may need recalibration

## Case Study: Vandeputte et al. (2017)

This Nature paper demonstrated that the Bacteroides/Prevotella trade-off in IBD is partially a load artifact. Their QMP method (flow cytometry + 16S) showed:

- Bacteroides doesn't actually "increase" in healthy
- Total bacterial load is lower in IBD
- The apparent trade-off is compositional closure

**This is our model for what a broader audit could reveal.**

## Data Requirements

1. **Raw counts** (not rarefied or normalized)
2. **Sample metadata** (disease status, covariates)
3. **Ideally**: qPCR total load, spike-in data, or flow cytometry

## Collaboration Opportunities

- Ravel lab: Extensive BV data with potential absolute quantification
- IBD consortia: Large cohorts with detailed metadata
- Human Microbiome Project: Reference data for healthy variation

## Timeline

| Phase | Scope | Deliverable |
|-------|-------|-------------|
| 1 | Literature survey | Curated list of 50 studies with extracted findings |
| 2 | Stress testing | Robustness scores for 20 datasets |
| 3 | Meta-analysis | Publication-ready figures and tables |

## Code Location

Analysis scripts will be added to `experiments/scripts/03-artifact-audit/` when formalized.

## References

- Vandeputte et al. (2017) Nature - QMP method paper
- Gloor et al. (2017) - Compositional data analysis review
- Morton et al. (2019) - Differential ranking methods
- Lin & Peddada (2020) - LinDA method paper
