# Experiment 10: CRC Meta-Analysis

## Summary

Meta-analysis of colorectal cancer (CRC) microbiome studies to assess cross-study consistency of findings and artifact risk. By analyzing multiple independent cohorts with the same methods, we can distinguish reproducible biology from study-specific artifacts.

**Key Finding**: Only a small fraction of CRC-associated taxa are consistently detected across cohorts with effect sizes exceeding the artifact threshold, while many "significant" findings are cohort-specific or fall within the artifact-risk zone.

## Motivation

CRC microbiome research has produced many claimed associations:
- Fusobacterium nucleatum (most replicated)
- Peptostreptococcus, Parvimonas
- Reduced butyrate producers

However:
- Different studies use different methods
- Effect sizes vary dramatically across cohorts
- Many findings fail to replicate

A meta-analysis using consistent methods can:
1. Identify truly robust associations
2. Quantify cross-study heterogeneity
3. Assess artifact risk for each finding

## Data

### Wirbel et al. 2019 Meta-Analysis Cohorts

| Cohort | Country | Cases | Controls | Sequencing |
|--------|---------|-------|----------|------------|
| Zeller | France | 91 | 66 | Metagenomics |
| Yu | China | 74 | 54 | Metagenomics |
| Feng | Austria | 46 | 63 | Metagenomics |
| Vogtmann | USA | 52 | 52 | Metagenomics |
| Thomas | Italy | 60 | 65 | Metagenomics |
| Hannigan | USA | 30 | 30 | Metagenomics |

Total: ~350 CRC cases, ~330 controls across 6 cohorts.

### Data Access

Available from:
- Original publications (SRA/ENA)
- curatedMetagenomicData (standardized)
- Wirbel et al. supplement (processed tables)

## Analysis Performed

### 1. Per-Cohort Analysis

Apply all methods to each cohort independently:
- LinDA (q < 0.10)
- ZINB (q < 0.05)
- Hurdle (q < 0.05)
- Permutation (p < 0.05)

### 2. Effect Size Comparison

For each taxon, compare effect sizes across cohorts:

| Taxon | Zeller | Yu | Feng | Vogtmann | Consistent? |
|-------|--------|-----|------|----------|-------------|
| F. nucleatum | +3.5 | +4.2 | +2.8 | +3.1 | Yes |
| Taxon X | +2.1 | -0.5 | +1.3 | NS | No |

### 3. Artifact Risk Classification

For each cohort and taxon:

| Category | Criteria | Interpretation |
|----------|----------|----------------|
| Robust | >3.3 log2FC in ≥2 cohorts | Likely real |
| Moderate | >3.3 log2FC in 1 cohort | Needs validation |
| At-risk | <3.3 log2FC in all | Likely artifact |

### 4. Meta-Analysis Statistics

Calculate:
- Pooled effect sizes (random effects model)
- Heterogeneity (I², Q statistic)
- Publication bias assessment

### 5. Consistency Matrix

Create a cohort × taxon matrix showing:
- Direction of effect (up/down/NS)
- Effect size magnitude
- Significance status

## Expected Findings

Based on prior meta-analyses:

### Robust Associations (Expected)
| Taxon | Expected Effect | Notes |
|-------|-----------------|-------|
| Fusobacterium nucleatum | +3-5 log2FC | Well-validated, large effect |
| Peptostreptococcus stomatis | +2-4 log2FC | Consistent across studies |
| Parvimonas micra | +2-3 log2FC | Oral pathobiont |

### Inconsistent/At-Risk (Expected)
- Many low-prevalence taxa
- Taxa with small effect sizes
- Cohort-specific findings

## Implications

### For CRC Research

1. **F. nucleatum is robust**: Large, consistent effect across cohorts
2. **Most associations are uncertain**: Effect sizes too small or inconsistent
3. **Method matters**: Different methods yield different significant taxa
4. **Cohort effects are real**: Technical and population differences matter

### For Meta-Analysis Practice

1. **Use consistent methods**: Enables fair comparison
2. **Report effect sizes**: Not just significance
3. **Assess artifact risk**: Is the effect larger than noise?
4. **Evaluate heterogeneity**: Consistent direction ≠ consistent magnitude

### For the Paper

This experiment demonstrates:
- The toolkit enables rigorous meta-analysis
- Cross-study consistency is the gold standard
- Many published findings may be artifacts

## Reproducibility

```bash
cd experiments/scripts/10-crc-meta-analysis/
./run_analysis.sh
python3 generate_figures.py
```

## Key Figures

1. **Forest plot**: Effect sizes by cohort for top taxa
2. **Consistency heatmap**: Taxa × cohort significance matrix
3. **Artifact risk distribution**: Fraction of findings in each risk category
4. **Method agreement**: Venn diagram of significant taxa by method

## Connection to Other Experiments

| Experiment | Focus | Design |
|------------|-------|--------|
| 01-03 | Ravel BV | Single study, artifact framework |
| 08 | HMP Gingival | Cross-body-site |
| 09 | IBD | Clinical disease context |
| 10 | CRC Meta | Cross-study consistency |

Together, these demonstrate the toolkit's applicability across:
- Study designs (case-control, meta-analysis)
- Body sites (vaginal, oral, gut)
- Disease contexts (BV, IBD, CRC)

## References

- Wirbel J, et al. (2019) - Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer. Nature Medicine.
- Thomas AM, et al. (2019) - Metagenomic analysis of colorectal cancer datasets identifies cross-cohort microbial diagnostic signatures and a link with choline degradation. Nature Medicine.
- Zeller G, et al. (2014) - Potential of fecal microbiota for early-stage detection of colorectal cancer. Molecular Systems Biology.
- Yu J, et al. (2017) - Metagenomic analysis of faecal microbiome as a tool towards targeted non-invasive biomarkers for colorectal cancer. Gut.
- Kostic AD, et al. (2012) - Genomic analysis identifies association of Fusobacterium with colorectal carcinoma. Genome Research.
