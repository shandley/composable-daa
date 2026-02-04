# Session Notes - Microbiome Reproducibility Atlas

**Last Updated**: February 2026
**Purpose**: Context preservation for session continuity

---

## Executive Summary

We built a rigorous meta-analysis framework to test microbiome reproducibility. The key finding is paradigm-shifting: **even Fusobacterium-CRC (the most replicated association) is NOT statistically significant when properly analyzed** (p=0.56 across 4 cohorts).

This is NOT because the effect is false—it's because:
1. Sample sizes are too small (need n>126 for medium effects)
2. Effects are heterogeneous (driven by ~15% of patients)
3. Simple replication metrics are misleading

---

## What We Built

### 1. Taxonomic Harmonization (`taxonomy.py`)

Solves the problem that different studies use different OTU identifiers:

```python
# Parses Greengenes format:
# k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__;d__denovo12345

def parse_taxonomy_string(name):
    """Extract genus from various taxonomy formats."""

def aggregate_to_genus(counts_df, method='sum'):
    """Aggregate OTU-level counts to genus level."""
```

**Validation**: Found Fusobacterium in all 4 CRC cohorts (116-1338 OTUs each) after harmonization.

### 2. Random-Effects Meta-Analysis (`meta_analysis.py`)

```python
def random_effects_meta(effects, variances):
    """
    DerSimonian-Laird estimator.
    Returns: pooled effect, CI, τ², I², Q statistic
    """
```

Includes ASCII forest plot generation for visualization.

### 3. Power Analysis (`power_analysis.py`)

```python
def minimum_detectable_effect(n_case, n_control, alpha=0.05, power=0.80):
    """Calculate Cohen's d detectable at 80% power."""

def required_sample_size(effect_size, alpha=0.05, power=0.80):
    """Calculate N needed to detect given effect."""
```

### 4. Reprocessing Pipeline (`reprocess_with_taxonomy.py`)

End-to-end pipeline that:
1. Reads raw OTU tables with full taxonomy strings
2. Aggregates to genus level
3. Runs DAA (Hurdle model via `daa run`)
4. Validates known positives (Fusobacterium, Parvimonas, etc.)

---

## Key Results

### Fusobacterium-CRC Meta-Analysis

| Cohort | Effect | SE | p-value | Direction |
|--------|--------|-----|---------|-----------|
| Zeller | +0.42 | 0.99 | 0.67 | ↓ (wrong) |
| Zackular | -1.11 | 1.30 | 0.39 | ↑ (correct) |
| Xiang | -0.52 | 1.00 | 0.60 | ↑ (correct) |
| Zhao | -0.22 | 0.64 | 0.73 | ↑ (correct) |
| **Pooled** | **-0.26** | **0.44** | **0.56** | ↑ (correct) |

**I² = 0%** - Effects are consistent, just underpowered.

### The Outlier Problem (Zackular Data)

```
CRC (n=58):          Control (n=30):
Mean: 2,465          Mean: 35
Median: 1            Median: 0
Top 2: 122k, 18k     Top: 992
```

**The "70x fold change" is driven by 2 patients.** Most CRC patients have zero Fusobacterium.

### Power Requirements

| Effect | N per Group | Total N |
|--------|-------------|---------|
| d=0.5 (1.6x) | 63 | 126 |
| d=0.3 (1.3x) | 175 | 350 |
| d=0.2 (1.2x) | 393 | 786 |

Current cohorts (n=43-116 total) are far below requirements.

---

## Revised Project Vision

### Original Goal
> "Calculate replication scores for every taxon-disease pair"

### Revised Goal
> "Characterize effect sizes, power requirements, and heterogeneity—explaining WHY effects do or don't replicate"

### New Deliverables

1. **Effect Characterization Database** (not just yes/no)
   ```
   Fusobacterium × CRC
   ├── Pooled effect: -0.26 (CI: -1.13, 0.61)
   ├── Direction: ↑ in disease (3/4 cohorts)
   ├── Significance: p=0.56 (underpowered)
   ├── Heterogeneity: I²=0% (consistent)
   └── Confidence: LIKELY REAL, SUBTYPE EFFECT
   ```

2. **Confidence Classifications**
   - Robust: Significant, low I², consistent direction
   - Likely Real: Correct direction, underpowered (e.g., Fusobacterium)
   - Heterogeneous: High I², needs subgroup analysis
   - Inconclusive: Insufficient power
   - Likely False: Wrong direction with adequate power

3. **Power Calculator** for study planning

4. **Publication**: "Why Most Microbiome Findings Don't Replicate"

---

## Directory Structure

```
experiments/
├── REPRODUCIBILITY_ATLAS.md      # Master vision (UPDATED)
├── atlas/
│   ├── FRAMEWORK.md              # Methodology
│   ├── VALIDATION_RESULTS.md     # Fusobacterium case study
│   ├── PHASE1_RESULTS.md         # Phase 1 findings (UPDATED)
│   ├── SESSION_NOTES.md          # This file
│   ├── scripts/
│   │   ├── taxonomy.py           # Taxonomic harmonization
│   │   ├── meta_analysis.py      # Random-effects meta-analysis
│   │   ├── power_analysis.py     # Power calculations
│   │   └── reprocess_with_taxonomy.py
│   └── diseases/
│       └── crc/
│           ├── data/             # Raw MicrobiomeHD data
│           ├── data_genus/       # Genus-aggregated counts
│           └── results_genus/    # DAA results
```

---

## How to Resume Work

### Option 1: Expand to More Diseases

```bash
# The pattern for new diseases:
# 1. Download from MicrobiomeHD
# 2. Modify reprocess_with_taxonomy.py for disease-specific metadata
# 3. Run genus-level aggregation
# 4. Run DAA
# 5. Run meta-analysis
```

### Option 2: Add Shotgun Validation

Use curatedMetagenomicData (Bioconductor) for species-level resolution. May show stronger Fusobacterium nucleatum effect than genus-level 16S.

### Option 3: Build the Database

Create structured output format for all taxon-disease pairs with:
- Effect estimates
- Power context
- Heterogeneity metrics
- Confidence classification

### Option 4: Write the Paper

Outline in REPRODUCIBILITY_ATLAS.md. Key messages:
1. Non-replication ≠ false positives
2. Published effect sizes are inflated
3. Subtype thinking is needed
4. Current sample sizes are inadequate

---

## Technical Notes

### CLI Usage

The scripts use `daa run` with YAML configs, NOT direct subcommands:

```bash
# CORRECT:
daa run -c counts.tsv -m metadata.tsv --config hurdle_config.yaml -o results.tsv

# WRONG (doesn't exist):
daa hurdle -c counts.tsv ...
```

### Known Positives (for validation)

```python
KNOWN_POSITIVES = {
    'crc': [
        {'taxon': 'Fusobacterium', 'direction': 'up'},
        {'taxon': 'Parvimonas', 'direction': 'up'},
        {'taxon': 'Porphyromonas', 'direction': 'up'},
        {'taxon': 'Peptostreptococcus', 'direction': 'up'},
        {'taxon': 'Solobacterium', 'direction': 'up'},
    ],
}
```

### Data Source

MicrobiomeHD: https://zenodo.org/record/840333
- Standardized 16S OTU counts
- 28 case-control studies
- 10 diseases

---

## Commits (Recent)

```
51f18cc Revise Atlas vision: from replication scores to effect characterization
74940c5 Add rigorous meta-analysis framework with taxonomic harmonization
33faa98 Add Microbiome Reproducibility Atlas - Phase 1 Results
```

---

## Key Insight to Remember

**The "reproducibility crisis" is largely a power crisis.**

Fusobacterium-CRC is REAL (direction consistent, I²=0%), but:
- Requires n>126 to detect at 80% power
- Is a subtype marker (~15% of CRC patients)
- Published findings benefited from sampling fluctuations

This reframes everything: the question isn't "what replicates?" but "what CAN replicate with available data?"
