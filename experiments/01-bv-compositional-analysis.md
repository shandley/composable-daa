# Experiment 01: BV Compositional Analysis

## Summary

Reanalysis of the Ravel 2011 vaginal microbiome dataset to demonstrate that the classic BV signature (Lactobacillus decrease, Gardnerella/Prevotella increase) is compositionally constrained and cannot be distinguished from total bacterial load effects.

## Background

Bacterial vaginosis (BV) is characterized by:
- Decreased Lactobacillus dominance
- Increased diversity (Gardnerella, Prevotella, Atopobium, etc.)
- Biofilm formation

The standard interpretation is that Lactobacillus "decreases" and pathogens "increase." However, BV involves biofilm formation which substantially increases total bacterial load. Without absolute quantification, we cannot distinguish:

1. **Scenario A**: Lactobacillus dies, pathogens fill the niche (constant load)
2. **Scenario B**: Lactobacillus unchanged, pathogens bloom on top (increased load)
3. **Scenario C**: Some combination

All three scenarios produce identical 16S relative abundance profiles.

## Data

**Source**: MicrobiomeBenchmarkData (Zenodo 6911027)
**Dataset**: Ravel 2011 vaginal microbiome
**Samples**: 394 total (248 healthy, 97 BV, 49 intermediate)
**Features**: 247 taxa

```bash
# Fetch data
daa fetch ravel-bv
```

## Analysis Performed

### 1. Standard LinDA/CLR Analysis

```bash
daa linda \
    -c counts_hbv.tsv \
    -m metadata_hbv.tsv \
    --formula "~ study_condition" \
    --test-coef study_conditionhealthy \
    --output linda_results.tsv
```

**Results**:
- 55 taxa tested (prevalence >10%)
- 47 significant at q<0.05 (85%)
- 28 higher in healthy, 19 lower in healthy

### 2. Compositional Diagnostic

**Sum of log2FC = 0.17** (essentially zero)

This is the signature of compositional closure under CLR normalization. Increases must be balanced by decreases - the data cannot tell us about absolute changes.

### 3. Top Differential Taxa

| Taxon | log2FC (healthy) | Interpretation |
|-------|------------------|----------------|
| L. crispatus | +4.44 | Higher in healthy |
| L. jensenii | +3.14 | Higher in healthy |
| Prevotella | -3.86 | Higher in BV |
| Megasphaera | -3.63 | Higher in BV |
| Gardnerella | -1.29 | Higher in BV |

### 4. Permutation Testing

5 permutation runs: **0 false positives** in all runs

Confirms statistical calibration is correct - the findings are not type I errors.

### 5. Sensitivity Analysis

If total bacterial load increases 10x in BV:

| Taxon | Observed | Load-corrected |
|-------|----------|----------------|
| L. crispatus | +4.44 | +7.76 (increases more!) |
| Prevotella | -3.86 | -0.54 (modest increase) |
| Gardnerella | -1.29 | +2.03 (actually increases) |

This completely changes the biological interpretation.

## Key Insight

The observation that "Lactobacillus decreases in BV" may be partially or entirely an artifact of total load increase. If true:

- Lactobacillus probiotics might work by displacing pathogens (not restoring Lactobacillus per se)
- Antimicrobial approaches targeting specific pathogens may be more effective
- Total bacterial load itself may be a better biomarker than composition

## Therapeutic Implications

| If True Mechanism Is... | Therapeutic Approach |
|------------------------|---------------------|
| Lactobacillus death | Lactobacillus probiotics |
| Pathogen bloom (load increase) | Targeted antimicrobials |
| Biofilm disruption needed | Biofilm-disrupting agents |

Current treatments may work for the wrong reasons, or the right treatment may depend on which mechanism dominates in a given patient.

## Next Steps

1. **Integrate absolute quantification**: Find BV datasets with qPCR total load measurements
2. **Apply spike-in normalization**: Use internal standards where available
3. **Meta-analysis**: Compare findings across multiple BV studies
4. **Collaboration opportunity**: Ravel lab has extensive BV data with potential absolute quantification

## Code Location

Analysis scripts will be added to `experiments/scripts/01-bv-analysis/` when formalized.

## References

- Ravel et al. (2011) PNAS - Original vaginal microbiome CST paper
- Vandeputte et al. (2017) Nature - QMP demonstration of compositional artifacts
- Stammler et al. (2016) - Spike-in methodology for absolute quantification
