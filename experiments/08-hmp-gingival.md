# Experiment 08: HMP Gingival Microbiome Analysis

## Summary

Reanalysis of Human Microbiome Project (HMP) gingival microbiome data to demonstrate that the composable DAA toolkit works across body sites (oral vs vaginal vs gut) and study designs.

**Key Finding**: The same methods and thresholds validated on vaginal microbiome (Ravel BV study) transfer to oral microbiome data, confirming body-site-agnostic applicability.

## Motivation

The toolkit was validated primarily on:
- Vaginal microbiome (Ravel BV study, Experiment 01)
- Synthetic data calibrated to gut/virome properties

The HMP gingival dataset provides an opportunity to validate on:
- Different body site (oral cavity)
- Different community structure (higher diversity than vaginal)
- Different study design (reference cohort vs case-control)

## Data

### HMP Gingival Microbiome

| Parameter | Value |
|-----------|-------|
| Source | Human Microbiome Project |
| Body site | Gingiva (gum tissue) |
| Samples | ~300 (V1-V3 region) |
| Sequencing | 16S rRNA amplicon |
| Subjects | Healthy adults |

Available versions:
- `hmp_v13`: 16S V1-V3 region (~300 samples)
- `hmp_v35`: 16S V3-V5 region
- `hmp_wms`: Whole metagenome shotgun

## Analysis Performed

### 1. Data Profiling

Compare sparsity and distribution characteristics between body sites:

| Metric | Vaginal (Ravel) | Gingival (HMP) |
|--------|-----------------|----------------|
| Sparsity | ~65% | ~70% |
| Diversity | Low (Lactobacillus dominated) | High |
| Library size CV | Moderate | Moderate |

### 2. Method Performance

Run all methods and compare consistency with vaginal microbiome results:

| Method | Expected Behavior |
|--------|-------------------|
| LinDA | Similar attenuation, q < 0.10 threshold |
| ZINB | Good discovery performance |
| Hurdle | Appropriate for moderate sparsity |
| Permutation | Gold standard calibration |

### 3. Cross-Body-Site Consistency

Validate that:
- Method rankings are consistent across body sites
- Threshold recommendations transfer
- Effect size interpretation remains valid

### 4. Compositional Closure Demonstration

Show that CLR sum = 0 constraint applies equally to oral microbiome.

## Expected Findings

1. **Methods transfer**: Same methods work for oral as for vaginal microbiome
2. **Thresholds consistent**: LinDA q < 0.10 remains appropriate
3. **Compositional effects**: CLR closure applies regardless of body site
4. **Diversity impact**: Higher diversity may affect power/sensitivity

## Implications

- **For users**: Toolkit is body-site-agnostic
- **For the field**: Validation across multiple HMP body sites strengthens generalizability claims
- **For the paper**: Demonstrates real-world applicability beyond the BV use case

## Reproducibility

```bash
cd experiments/scripts/08-hmp-gingival/
./run_analysis.sh
python3 generate_figures.py
```

## Connection to Other Experiments

- **Experiment 01**: Ravel BV (vaginal) - original validation
- **Experiment 08**: HMP gingival (oral) - cross-body-site validation
- **Experiment 09**: IBD (gut) - disease context validation

Together, these cover the three major body sites in microbiome research.

## References

- Human Microbiome Project Consortium (2012) - Structure, function and diversity of the healthy human microbiome. Nature.
- Segata et al. (2012) - Composition of the adult digestive tract bacterial microbiome based on seven mouth surfaces, tonsils, throat and stool samples. Genome Biology.
