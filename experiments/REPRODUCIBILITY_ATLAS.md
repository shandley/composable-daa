# The Microbiome Reproducibility Atlas

## Executive Summary

A systematic initiative to calculate **replication scores** for every published microbiome-disease association by reanalyzing all publicly available standardized microbiome datasets using the composable-daa toolkit.

**Core Question**: What fraction of published microbiome findings actually replicate across independent cohorts?

**Preliminary Answer**: Based on our CRC meta-analysis (Experiment 11), only ~15% of significant findings replicate across all available cohorts. This initiative will determine if this pattern holds across all diseases.

---

## Motivation

### The Problem

1. **Single-cohort studies dominate** - Most microbiome papers report findings from one cohort
2. **Replication is rare** - Few studies validate findings in independent populations
3. **The field lacks a reproducibility metric** - No systematic assessment exists
4. **Our Experiment 11 finding** - Only 15% of CRC-microbiome associations replicate across 3 cohorts

### The Opportunity

Massive standardized databases now exist:

| Database | Samples | Diseases | Status |
|----------|---------|----------|--------|
| GMrepo v3 | 118,965 | 301 | Standardized |
| QIITA | 168,000+ | Diverse | Standardized |
| mBodyMap | 63,148 | 56 | Standardized |
| PRIME | 53,449 | 101 | Standardized |
| curatedMetagenomicData | ~10,000+ | Multiple | Standardized |
| MicrobiomeHD | ~2,500 | 10 | Standardized |
| **Total** | **~400,000+** | **300+** | |

### The Solution

Use composable-daa to systematically reanalyze ALL available cohorts and calculate cross-study replication rates for every taxon-disease pair.

---

## Methodology

### Analysis Pipeline

For each disease with 2+ available cohorts:

```yaml
# Standard analysis pipeline
steps:
- !FilterPrevalence
  threshold: 0.10
- !AddPseudocount
  value: 0.5
- NormalizeCLR
- !ModelLM
  formula: "~ group"
- !TestWald
  coefficient: group{disease}
- CorrectBH
```

Plus Hurdle model for discovery and Permutation for validation.

### Replication Score Definition

For each taxon-disease pair:

```
Replication Score = (# cohorts where q < threshold) / (# cohorts tested)

Classification:
- Score = 1.0 (100%): ROBUST - significant in all cohorts
- Score ≥ 0.67 (67%+): HIGH - significant in most cohorts
- Score ≥ 0.33 (33%+): MODERATE - significant in some cohorts
- Score < 0.33: LOW - mostly study-specific
- Score = 0: NOT REPLICATED - never significant
```

### Additional Metrics

1. **Direction Consistency**: Do effect sizes agree in direction across cohorts?
2. **Effect Size Range**: What is the range of observed effect sizes?
3. **Method Agreement**: Do LinDA and Hurdle agree?
4. **Artifact Risk**: Is the effect size above the 3.3 log2FC artifact threshold?

---

## Phase 1: Top 10 Diseases

### Disease Selection Criteria

1. High public health impact
2. Multiple independent cohorts available
3. Prior literature on microbiome associations
4. Mix of gut and other body sites

### Selected Diseases

| Priority | Disease | Expected Cohorts | Primary Database |
|----------|---------|------------------|------------------|
| 1 | **Colorectal Cancer (CRC)** | 3+ | MicrobiomeHD, GMrepo |
| 2 | **Inflammatory Bowel Disease (IBD)** | 5+ | MicrobiomeHD, GMrepo |
| 3 | **Type 2 Diabetes (T2D)** | 3+ | GMrepo, curatedMetagenomicData |
| 4 | **Obesity** | 5+ | MicrobiomeHD, GMrepo |
| 5 | **Crohn's Disease** | 3+ | MicrobiomeHD, GMrepo |
| 6 | **Ulcerative Colitis** | 3+ | GMrepo |
| 7 | **Liver Cirrhosis** | 2+ | GMrepo |
| 8 | **Parkinson's Disease** | 2+ | GMrepo |
| 9 | **Depression/Anxiety** | 2+ | GMrepo |
| 10 | **Rheumatoid Arthritis** | 2+ | GMrepo |

### Data Sources

#### Primary: MicrobiomeHD (Zenodo 840333)
- Standardized 16S OTU counts
- 28 case-control studies
- 10 diseases
- Consistent processing pipeline

**Available datasets:**
- CRC: Baxter, Zeller, Zackular (DONE in Exp 11)
- IBD: Multiple cohorts
- Obesity: Multiple cohorts
- CDI (C. difficile): Multiple cohorts
- And more...

#### Secondary: GMrepo v3
- 118,965 samples
- 301 diseases
- Both 16S and WGS
- API access available

#### Tertiary: curatedMetagenomicData
- Shotgun metagenomics
- Functional data
- Bioconductor access

---

## Implementation Plan

### Phase 1A: MicrobiomeHD Complete Analysis

MicrobiomeHD contains these disease datasets:

| Disease | Studies | Samples | Status |
|---------|---------|---------|--------|
| CRC | 3 | 666 | **DONE (Exp 11)** |
| CDI | 4 | ~500 | Pending |
| IBD | 4 | ~800 | Pending |
| Obesity | 3 | ~400 | Pending |
| T1D | 2 | ~100 | Pending |
| HIV | 2 | ~200 | Pending |
| Cirrhosis | 1 | ~200 | Pending |
| NASH | 1 | ~100 | Pending |
| Autism | 1 | ~50 | Pending |
| RA | 1 | ~100 | Pending |

### Phase 1B: Expand with GMrepo

For diseases with only 1-2 cohorts in MicrobiomeHD, supplement with GMrepo data.

### Directory Structure

```
experiments/
├── REPRODUCIBILITY_ATLAS.md          # This document
├── atlas/
│   ├── README.md                     # Atlas overview
│   ├── diseases/
│   │   ├── crc/                      # Colorectal cancer
│   │   │   ├── cohorts.yaml          # List of cohorts
│   │   │   ├── run_analysis.sh       # Analysis script
│   │   │   ├── results/              # Per-cohort results
│   │   │   └── replication_scores.tsv
│   │   ├── ibd/                      # Inflammatory bowel disease
│   │   ├── obesity/
│   │   ├── t2d/
│   │   └── ...
│   ├── scripts/
│   │   ├── download_microbiomehd.sh
│   │   ├── download_gmrepo.py
│   │   ├── run_all_analyses.sh
│   │   └── calculate_replication.py
│   └── results/
│       ├── disease_summary.tsv       # Replication rates by disease
│       ├── taxon_scores.tsv          # All taxon-disease scores
│       └── figures/
```

---

## Expected Outcomes

### Primary Deliverable

**The Replication Score Database**: A searchable resource where researchers can look up any taxon-disease pair and see:
- How many cohorts tested it
- How many found it significant
- Effect size range and direction consistency
- Confidence classification

### Key Findings (Hypotheses)

1. **Overall replication rate**: We predict 10-20% of findings replicate across all cohorts
2. **Disease variation**: Some diseases (IBD?) may have higher replication than others
3. **Taxon variation**: Some taxa (Fusobacterium in CRC?) may be consistently robust
4. **Method matters**: LinDA vs Hurdle may show different replication patterns

### Publication

**Proposed Title**: "The Microbiome Reproducibility Atlas: Systematic Cross-Cohort Validation of Disease Associations"

**Key Message**: The majority of published microbiome-disease associations do not replicate across independent cohorts. We provide replication scores for X,XXX taxon-disease pairs across Y diseases.

---

## Technical Notes

### Composable-DAA Advantages

1. **Speed**: Rust implementation handles large datasets efficiently
2. **Consistency**: Same pipeline on all cohorts ensures fair comparison
3. **Reproducibility**: YAML configs document exact methods
4. **Validation**: Built-in spike-in validation available

### Challenges

1. **Taxonomic harmonization**: Different studies use different databases (Greengenes, SILVA)
2. **Metadata standardization**: Disease labels vary across studies
3. **Batch effects**: Technical differences between studies
4. **Sample size variation**: Some cohorts are small (n<30)

### Quality Filters

- Minimum 20 samples per group
- Minimum 10% prevalence filter
- Require raw counts (not pre-normalized)
- Document any exclusions

---

## Progress Log

### 2024-02-01: Phase 1 Complete

**STRIKING FINDING: <1% of microbiome findings replicate across cohorts**

| Disease | Cohorts | Significant Taxa | Replicate in ALL | Rate |
|---------|---------|------------------|------------------|------|
| CRC | 4 | 3,244 | 0 | 0.0% |
| IBD | 2 | 347 | 0 | 0.0% |
| CDI | 2 | 561 | 1 | 0.2% |
| **Total** | **8** | **4,152** | **1** | **0.02%** |

This is far worse than our initial Experiment 11 finding (15% with 3 CRC cohorts). As we add more cohorts and analyze more diseases, the replication rate drops to near zero.

See `atlas/PHASE1_RESULTS.md` for full details.

### Phase 1 Completed

- Analyzed 3 diseases (CRC, IBD, CDI) with 8 total cohorts
- All analyses run with composable-daa toolkit
- LinDA and Hurdle methods on each cohort
- Cross-cohort replication calculated

### Next Steps

1. Expand to GMrepo (118,965 samples, 301 diseases)
2. Build searchable replication score database
3. Write publication: "The Microbiome Reproducibility Atlas"

---

## References

1. Duvallet C, et al. (2017). Meta-analysis of gut microbiome studies identifies disease-specific and shared responses. Nat Commun.
2. MicrobiomeHD: https://zenodo.org/record/840333
3. GMrepo: https://gmrepo.humangut.info/
4. PRIME: https://prime.microbiome.cloud/
5. curatedMetagenomicData: https://waldronlab.io/curatedMetagenomicData/
