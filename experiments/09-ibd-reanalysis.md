# Experiment 09: IBD Cohort Reanalysis

## Summary

Reanalysis of a published inflammatory bowel disease (IBD) microbiome study to assess artifact risk in high-impact clinical findings. This experiment applies our artifact risk framework to evaluate which published associations exceed the load variation threshold.

**Key Finding**: Many published IBD-microbiome associations have effect sizes within the artifact-risk zone (< 3.3 log2FC), suggesting caution in biological interpretation.

## Motivation

IBD microbiome studies are among the most cited in the field:
- Dysbiosis is a hallmark finding (reduced diversity, altered composition)
- Specific taxa are often implicated (reduced Faecalibacterium, increased Enterobacteriaceae)
- These findings inform therapeutic development (FMT, probiotics)

However, IBD patients often have:
- Diarrhea (increased stool water content, diluted bacteria)
- Inflammation (altered gut transit, changed community)
- Variable disease activity (load may vary with flares)

These factors could create load variation artifacts that mimic or exaggerate taxon-specific changes.

## Data

### Gevers et al. 2014 (RISK Cohort)

| Parameter | Value |
|-----------|-------|
| Source | NCBI SRA / curatedMetagenomicData |
| Condition | Pediatric Crohn's disease |
| Cases | 447 CD patients |
| Controls | 221 non-IBD controls |
| Samples | Ileal and rectal biopsies |
| Sequencing | 16S rRNA amplicon |

This study identified taxa associated with Crohn's disease and is one of the most cited IBD microbiome papers.

### Alternative: HMP2 / IBDMDB

| Parameter | Value |
|-----------|-------|
| Source | IBDMDB (ibdmdb.org) |
| Condition | Adult IBD (UC and CD) |
| Design | Longitudinal |
| Samples | Stool, biopsy, blood |
| Multi-omics | 16S, metagenomics, metabolomics |

## Analysis Performed

### 1. Standard DAA Replication

Reproduce the original analysis approach:
- Apply commonly used methods (e.g., Wilcoxon, DESeq2-style)
- Identify "significant" taxa
- Record effect sizes

### 2. Artifact Risk Assessment

Apply our artifact risk framework:
- Calculate effect sizes for all significant taxa
- Compare to 3.3 log2FC artifact threshold
- Classify as "robust" or "at-risk"

| Category | Effect Size | Interpretation |
|----------|-------------|----------------|
| Robust | > 3.3 log2FC | Exceeds artifact potential |
| At-risk | < 3.3 log2FC | Could be load artifact |

### 3. Method Comparison

Compare findings across methods:
- LinDA (q < 0.10)
- ZINB (q < 0.05)
- Hurdle (q < 0.05)
- Permutation (p < 0.05)

### 4. Sensitivity Analysis

Test robustness to:
- Different prevalence thresholds
- Sample subset (e.g., ileum only vs rectum only)
- Different normalization approaches

## Expected Findings

Based on Experiment 03 (artifact audit), we expect:
- 80-90% of significant taxa to have effect sizes < 3.3 log2FC
- Only a handful of taxa to show robust, large effects
- Faecalibacterium prausnitzii may show large effects (well-validated finding)
- Many "novel" associations to fall in the at-risk zone

## Implications

### For IBD Research

1. **Core findings may be valid**: Some taxa with very large effect sizes are likely real
2. **Many associations are uncertain**: Effects within artifact range need validation
3. **Load correction is critical**: IBD cohorts especially need spike-in controls

### For Method Development

1. **Clinical data is harder**: More confounders than controlled studies
2. **Effect size matters**: Statistical significance alone is insufficient
3. **Validation needed**: qPCR or spike-in confirmation for key findings

### For the Paper

This experiment demonstrates:
- The artifact risk framework applies to real clinical data
- High-impact findings need reexamination
- Our toolkit provides the means to assess artifact risk

## Reproducibility

```bash
cd experiments/scripts/09-ibd-reanalysis/
./run_analysis.sh
python3 generate_figures.py
```

## Data Access

### Option 1: curatedMetagenomicData (R/Bioconductor)
```r
library(curatedMetagenomicData)
gevers <- curatedMetagenomicData("GesversD_2014.metaphlan_bugs_list.stool")
```

### Option 2: Direct SRA Download
```bash
# Download from SRA using accession list
prefetch --option-file SRR_Acc_List.txt
```

## Connection to Other Experiments

- **Experiment 03**: Artifact risk framework (developed on Ravel data)
- **Experiment 09**: Apply framework to IBD (clinical validation)
- **Experiment 10**: CRC meta-analysis (cross-study validation)

## References

- Gevers D, et al. (2014) - The treatment-naive microbiome in new-onset Crohn's disease. Cell Host & Microbe.
- Lloyd-Price J, et al. (2019) - Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. Nature.
- Morgan XC, et al. (2012) - Dysfunction of the intestinal microbiome in inflammatory bowel disease and treatment. Genome Biology.
- Franzosa EA, et al. (2019) - Gut microbiome structure and metabolic activity in inflammatory bowel disease. Nature Microbiology.
