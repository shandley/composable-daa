# Experiment 01: BV Compositional Analysis

## Summary

Reanalysis of the Ravel 2011 vaginal microbiome dataset demonstrating that the classic BV signature (Lactobacillus decrease, Gardnerella/Prevotella increase) is compositionally constrained and cannot be distinguished from total bacterial load effects without absolute quantification.

**Key Finding**: The sum of CLR log fold changes equals exactly zero (0.0000), proving compositional closure. This means that observed "increases" and "decreases" are mathematically coupled - they are not independent biological observations.

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
**Analysis Samples**: 345 (248 healthy, 97 BV) - intermediate excluded
**Features Tested**: 55 taxa (after 10% prevalence filtering)
**Sparsity**: 92.1%

```bash
# Fetch and filter data
daa fetch -d ravel -o ./data/
# Filter to healthy vs BV in Python (see run_analysis.sh)
```

## Analysis Performed

### 1. Data Profiling

```bash
daa profile-llm -c counts_hbv.tsv -m metadata_hbv.tsv -g study_condition
```

Key metrics:
- Overall sparsity: 92.1% (very sparse → Hurdle model recommended)
- Library size CV: 0.47 (moderate variation)
- Samples: 248 healthy, 97 BV

### 2. Hurdle Model Analysis (Recommended for Sparse Data)

```bash
daa recommend -c counts_hbv.tsv -m metadata_hbv.tsv \
    -g study_condition -t healthy --run -o hurdle_results.tsv
```

**Results**:
- 55 taxa tested (10% prevalence threshold)
- 28 significant at q<0.05
- Sum of estimates: **-30.45** (not zero - count-based model)

### 3. LinDA (CLR) Analysis for Compositional Demonstration

```yaml
# linda_clr_pipeline.yaml
name: LinDA_CLR
steps:
- !FilterPrevalence
  threshold: 0.1
- !AddPseudocount
  value: 0.5
- !NormalizeCLR
- !ModelLM
  formula: "~ study_condition"
- !TestWald
  coefficient: study_conditionhealthy
- CorrectBH
```

```bash
daa run -c counts_hbv.tsv -m metadata_hbv.tsv \
    --config linda_clr_pipeline.yaml -o linda_clr_results.tsv
```

**Results**:
- 55 taxa tested
- 47 significant at q<0.05 (85%)
- **Sum of estimates: 0.0000** ← KEY FINDING

### 4. Compositional Closure Proof

| Method | Sum of Estimates | Interpretation |
|--------|------------------|----------------|
| Hurdle (count-based) | -30.45 | No compositional constraint |
| LinDA (CLR-based) | **0.0000** | Perfect compositional closure |

The CLR transformation forces the sum of log fold changes to zero. This is a mathematical property, not a biological finding.

**Implication**: If Lactobacillus shows log2FC = +4.44 in healthy, other taxa must show a combined log2FC = -4.44 to maintain the sum of zero. The "increase" in pathogens is mathematically required, not independently observed.

### 5. Top Differential Taxa

| Taxon | log2FC (healthy) | Std Error | q-value | Interpretation |
|-------|------------------|-----------|---------|----------------|
| L. crispatus | +4.44 | 0.37 | <0.001 | Higher in healthy |
| L. jensenii | +3.14 | 0.29 | <0.001 | Higher in healthy |
| L. iners | +1.87 | 0.39 | <0.001 | Higher in healthy |
| Prevotella | -3.86 | 0.18 | <0.001 | Higher in BV |
| Megasphaera | -3.63 | 0.18 | <0.001 | Higher in BV |
| Atopobium | -3.14 | 0.16 | <0.001 | Higher in BV |
| Sneathia | -3.05 | 0.18 | <0.001 | Higher in BV |
| Gardnerella | -1.29 | 0.11 | <0.001 | Higher in BV |

### 6. Sensitivity Analysis: Effect of Total Load

Formula: `True_log2FC = Observed_CLR + log2(Load_BV / Load_Healthy)`

Reference: Stammler et al. (2016) showed 9.8x variation in bacterial load via spike-in normalization.

| Taxon | Observed CLR | 1x Load | 2x Load | 5x Load | 10x Load | 20x Load |
|-------|--------------|---------|---------|---------|----------|----------|
| L. crispatus | +4.44 | +4.44 | +5.44 | +6.76 | +7.76 | +8.76 |
| L. iners | +1.87 | +1.87 | +2.87 | +4.20 | +5.20 | +6.20 |
| L. jensenii | +3.14 | +3.14 | +4.14 | +5.46 | +6.46 | +7.46 |
| Prevotella | -3.86 | -3.86 | -2.86 | -1.54 | **-0.54** | +0.46 |
| Gardnerella | -1.29 | -1.29 | -0.29 | +1.03 | **+2.03** | +3.03 |
| Megasphaera | -3.63 | -3.63 | -2.63 | -1.31 | **-0.31** | +0.69 |
| Sneathia | -3.05 | -3.05 | -2.05 | -0.73 | +0.27 | +1.27 |
| Atopobium | -3.14 | -3.14 | -2.14 | -0.82 | +0.18 | +1.18 |

**At 10x load (plausible for BV biofilm)**:
- **Prevotella**: Nearly UNCHANGED in absolute terms (-0.54 ≈ 0)
- **Gardnerella**: Actually HIGHER in BV (+2.03)
- **Lactobacillus**: Even MORE protective (+7.76 vs +4.44)

## Key Insights

### 1. Compositional Closure is Real

The sum of CLR estimates = 0.0000 exactly. This proves:
- Taxa cannot change independently under CLR normalization
- Every "increase" requires a compensating "decrease"
- Classic BV findings are mathematically coupled, not independent

### 2. The Lactobacillus Paradox

The observed "decrease" in Lactobacillus during BV may be:
1. **Real decrease**: Lactobacillus cells die or are displaced
2. **Dilution artifact**: More total bacteria → lower Lactobacillus fraction
3. **Compositional artifact**: Sum-to-zero constraint forces apparent decrease

Without absolute quantification, these are indistinguishable.

### 3. BV "Pathogens" May Be Neutral

If BV involves a 10x increase in total load (plausible for biofilm):
- Prevotella, Megasphaera, Atopobium may be UNCHANGED in absolute abundance
- Their apparent "increase" is entirely due to Lactobacillus fraction dropping
- They may be opportunistic colonizers, not drivers of dysbiosis

### 4. Gardnerella May Be Higher in BV (Absolutely)

Even with load correction, Gardnerella shows a positive effect at 10x load. This suggests Gardnerella genuinely increases in absolute terms, making it a stronger candidate for a causative role than Prevotella.

## Therapeutic Implications

| If True Mechanism Is... | Therapeutic Approach |
|------------------------|---------------------|
| Lactobacillus death | Lactobacillus probiotics |
| Pathogen bloom (load increase) | Targeted antimicrobials |
| Biofilm disruption needed | Biofilm-disrupting agents |
| Gardnerella is driver | Gardnerella-specific therapy |

Current treatments may work for the wrong reasons, or the right treatment may depend on which mechanism dominates.

## Composability Enables Discovery

This analysis would be difficult with monolithic tools like DESeq2 or ANCOM-BC because:

1. **Hurdle vs LinDA comparison**: Requires running two different model types
2. **Sum of estimates**: Not a standard output from any tool
3. **Sensitivity analysis**: Requires manual calculation on raw results
4. **Pipeline customization**: CLR + LM combination not pre-packaged

The composable-daa approach allowed us to:
- Mix normalization (CLR) with simple linear model
- Extract raw estimates for post-hoc analysis
- Compare count-based (Hurdle) vs ratio-based (CLR) methods
- Quantify compositional closure directly

## Reproducibility

All analysis is reproducible via:

```bash
cd experiments/scripts/01-bv-analysis/
./run_analysis.sh
```

Output files:
- `results/hurdle_results.tsv` - Hurdle model results
- `results/linda_clr_results.tsv` - LinDA CLR results
- `results/sensitivity_analysis.tsv` - Load-corrected estimates
- `results/compositional_summary.txt` - Closure verification
- `results/analysis_report.md` - Summary report

## Next Steps

1. **Integrate absolute quantification**: Find BV datasets with qPCR total load
2. **Validate with spike-in data**: Apply to Stammler et al. dataset
3. **Cross-validate methods**: Compare Hurdle vs ZINB vs LinDA systematically
4. **Extend to other conditions**: Apply framework to IBD, obesity, etc.
5. **Publication**: Write up as methods paper with BV as case study

## References

- Ravel J, et al. (2011). Vaginal microbiome of reproductive-age women. PNAS.
- Stammler F, et al. (2016). Adjustment of 16S rRNA gene sequence data by spike-in. Sci Rep.
- Vandeputte D, et al. (2017). Quantitative microbiome profiling links gut community variation to microbial load. Nature.
- Aitchison J (1986). The Statistical Analysis of Compositional Data. Chapman & Hall.
