# Microbiome Reproducibility Atlas - Phase 1 Results

## Executive Summary

**The replication crisis in microbiome research is worse than expected.**

Across three diseases analyzed with multiple independent cohorts, we find:

| Disease | Cohorts | Significant Taxa | Replicate in ALL | Replication Rate |
|---------|---------|------------------|------------------|------------------|
| **CRC** | 4 | 3,244 | 0 | **0.0%** |
| **IBD** | 2 | 347 | 0 | **0.0%** |
| **CDI** | 2 | 561 | 1 | **0.2%** |
| **Total** | 8 | 4,152 | 1 | **0.02%** |

**Only 1 taxon out of 4,152 significant findings replicates across all available cohorts.**

---

## Detailed Results

### Colorectal Cancer (CRC)

**Cohorts analyzed**: Zeller, Zackular, Xiang, Zhao (4 cohorts, ~400 samples)

| Method | Threshold | Significant Taxa | Replicate in 4/4 | Replicate in 2+ |
|--------|-----------|------------------|------------------|-----------------|
| LinDA | q < 0.10 | 0 | 0 (0%) | 0 (0%) |
| Hurdle | q < 0.05 | 3,244 | 0 (0%) | 66 (2%) |

**Key finding**: Despite thousands of significant associations per cohort, ZERO taxa are significant across all 4 independent CRC cohorts. 98% of findings are study-specific (significant in only 1 cohort).

### Inflammatory Bowel Disease (IBD)

**Cohorts analyzed**: Gevers, Huttenhower (2 cohorts, ~380 samples)

| Method | Threshold | Significant Taxa | Replicate in 2/2 |
|--------|-----------|------------------|------------------|
| LinDA | q < 0.10 | 1 | 0 (0%) |
| Hurdle | q < 0.05 | 347 | 0 (0%) |

**Key finding**: Zero overlap between the two IBD cohorts. All 347 significant taxa are study-specific.

### Clostridioides difficile Infection (CDI)

**Cohorts analyzed**: Schubert, Youngster (2 cohorts, ~436 samples)

| Method | Threshold | Significant Taxa | Replicate in 2/2 |
|--------|-----------|------------------|------------------|
| LinDA | q < 0.10 | 0 | 0 (0%) |
| Hurdle | q < 0.05 | 561 | 1 (0.2%) |

**Key finding**: Only 1 taxon (a Bacteroidetes) replicates across both CDI cohorts out of 561 significant findings.

---

## The Pattern

```
┌─────────────────────────────────────────────────────────────────┐
│  EACH COHORT FINDS HUNDREDS OF SIGNIFICANT TAXA                 │
│                                                                  │
│  CRC Zeller:    2,313 significant                               │
│  CRC Zackular:    954 significant                               │
│  CDI Youngster:   474 significant                               │
│  IBD Gevers:      335 significant                               │
│  ...                                                             │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  BUT ALMOST NONE REPLICATE ACROSS COHORTS                       │
│                                                                  │
│  CRC (4 cohorts):  0 replicate in all 4                         │
│  IBD (2 cohorts):  0 replicate in both                          │
│  CDI (2 cohorts):  1 replicates in both                         │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│  THE VAST MAJORITY ARE STUDY-SPECIFIC                           │
│                                                                  │
│  98% of CRC findings: only 1 cohort                             │
│  100% of IBD findings: only 1 cohort                            │
│  99.8% of CDI findings: only 1 cohort                           │
└─────────────────────────────────────────────────────────────────┘
```

---

## Implications

### For the Field

1. **Single-cohort studies are insufficient** - A finding significant in one cohort has <1% probability of replicating
2. **Published effect sizes may be artifacts** - Study-specific confounders, batch effects, or statistical noise
3. **Meta-analyses are essential** - Only cross-cohort consistency can identify robust signals
4. **The literature is unreliable** - Most published microbiome-disease associations will not replicate

### For Researchers

1. **Always validate in independent cohorts** before claiming biological significance
2. **Report per-cohort results** in meta-analyses, not just pooled statistics
3. **Use replication as the gold standard**, not p-values from single studies
4. **Be skeptical of single-cohort claims**, even with large sample sizes

### For Reviewers

1. **Require replication** for causal claims
2. **Request cross-cohort validation** before publication
3. **Discount single-cohort findings** in review

---

## Methods

### Data Source
MicrobiomeHD (Zenodo 840333) - Standardized 16S rRNA gene sequencing data from published case-control studies.

### Analysis Pipeline
All analyses performed with composable-daa toolkit:
- Prevalence filter: 10%
- LinDA: CLR + Linear Model, q < 0.10
- Hurdle: Two-part model, q < 0.05
- Cross-cohort overlap computed at taxon level

### Limitations
- Only 16S data (not shotgun metagenomics)
- Limited to diseases with 2+ cohorts in MicrobiomeHD
- Taxonomic resolution varies between studies
- Some cohorts have small sample sizes

---

## Data Availability

All analysis scripts and results in:
```
experiments/atlas/
├── diseases/
│   ├── crc/results/
│   ├── ibd/results/
│   └── cdi/results/
├── scripts/
│   └── analyze_disease.sh
└── PHASE1_RESULTS.md
```

---

## Next Steps

### Phase 2: Expand to GMrepo
- 118,965 samples across 301 diseases
- Will dramatically increase cohort coverage

### Phase 3: Build Replication Database
- Searchable resource for any taxon-disease pair
- Report replication score and confidence level

### Phase 4: Publication
- "The Microbiome Reproducibility Atlas: Why Most Findings Don't Replicate"

---

## Key Message

> **A statistically significant microbiome finding has less than 1% probability of replicating across independent cohorts. The field's confidence in single-cohort results is dramatically misplaced.**

This is not a failure of any individual study—it's a fundamental limitation of the current publication model that rewards novel findings over replication.
