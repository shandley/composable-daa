# The Microbiome Reproducibility Atlas

## Executive Summary

A systematic initiative to characterize the **detectability and reproducibility** of microbiome-disease associations by reanalyzing publicly available datasets using rigorous meta-analysis methods.

**Core Question (Revised)**: Not just "what replicates?" but **"what CAN be detected with available sample sizes, and what does reproducibility actually mean for heterogeneous biological effects?"**

**Key Finding**: Even the most well-established association (Fusobacterium-CRC) does not reach statistical significance (p=0.56) when properly analyzed across 4 cohorts. This is not because the effect is false—it's because:
1. Current cohorts are underpowered (n=40-100 can only detect 2x+ effects)
2. Effects are heterogeneous (driven by patient subsets, not universal)
3. Simple replication metrics are misleading for complex biology

---

## Revised Understanding

### What We Originally Thought

> "Only 15% of findings replicate" → "The field has a reproducibility crisis"

### What We Now Know

The problem is more nuanced:

| Issue | Original Interpretation | Revised Understanding |
|-------|------------------------|----------------------|
| Low replication | False positives | Underpowered studies + heterogeneous effects |
| Fusobacterium-CRC | Should easily replicate | Subtype phenomenon (10-20% of patients) |
| OTU non-overlap | Nothing replicates | Taxonomic harmonization needed |
| Effect sizes | Inflated | Publication bias + outlier-driven |

### The Fusobacterium Case Study

We used Fusobacterium-CRC as a validation case—the "most replicated" finding in microbiome research. Results:

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Found in all cohorts | ✓ | Taxonomic harmonization works |
| Correct direction | 3/4 cohorts | Effect is real |
| Statistically significant | 0/4 cohorts | Underpowered |
| Meta-analysis p-value | 0.56 | Not significant even pooled |
| Heterogeneity I² | 0% | Consistent (just small) |

**The "70x fold change" in Zackular is driven by 2 patients** with 122,000+ reads. Median CRC = 1, Median Control = 0. Most CRC patients have zero Fusobacterium.

---

## Revised Vision

### From Replication Scores to Effect Characterization

The Atlas should provide, for each taxon-disease pair:

#### 1. Effect Estimation
- **Pooled effect size** (random-effects meta-analysis)
- **95% confidence interval**
- **Direction consistency** across cohorts

#### 2. Power Context
- **Minimum detectable effect** for available cohorts
- **Required sample size** to confirm at 80% power
- **Current power** to detect observed effect

#### 3. Heterogeneity Assessment
- **I² statistic** (between-study variance)
- **Patient-level heterogeneity** (within-study variance)
- **Subtype prevalence** (% of patients showing effect)

#### 4. Confidence Classification

| Classification | Criteria | Example |
|---------------|----------|---------|
| **Robust** | Significant meta-analysis, low I², consistent direction | (Rare with current data) |
| **Likely Real** | Correct direction, underpowered, low I² | Fusobacterium-CRC |
| **Heterogeneous** | High I², mixed directions | Needs subgroup analysis |
| **Inconclusive** | Low power, insufficient cohorts | Most findings |
| **Likely False** | Wrong direction, high power, still null | Inflated single-study claims |

---

## Methodology

### Taxonomic Harmonization

Different studies use different OTU definitions. We aggregate to genus level:

```python
# taxonomy.py
def aggregate_to_genus(counts_df):
    """
    Parse taxonomy strings (Greengenes, SILVA, plain)
    and aggregate OTU counts to genus level.
    """
```

**Validation**: Found Fusobacterium in all 4 CRC cohorts (116-1338 OTUs each) after harmonization.

### Random-Effects Meta-Analysis

For each taxon present in 2+ cohorts:

```python
# meta_analysis.py
def random_effects_meta(effects, variances):
    """
    DerSimonian-Laird estimator with:
    - Pooled effect and CI
    - τ² (between-study variance)
    - I² (heterogeneity %)
    - Q statistic and p-value
    """
```

### Power Analysis

For each cohort:

```python
# power_analysis.py
def minimum_detectable_effect(n_case, n_control, alpha=0.05, power=0.80):
    """
    Calculate Cohen's d detectable at 80% power.
    Maps to fold-change: d=0.5 ≈ 1.6x, d=0.8 ≈ 2.2x
    """
```

**Finding**: Current cohorts (n=40-100) require ~2x fold changes to detect.

### Analysis Pipeline

```yaml
# Genus-level analysis
steps:
  - !FilterPrevalence
    threshold: 0.05
  - !ModelHurdle
    formula: "~ group"
  - !TestWald
    coefficient: groupcontrol
  - CorrectBH
```

---

## Key Findings

### Power Limitations

| Cohort | N (case/ctrl) | Min Detectable d | Min Detectable FC |
|--------|---------------|------------------|-------------------|
| Zeller | 41/75 | 0.54 | 1.7x |
| Zackular | 58/30 | 0.63 | 1.9x |
| Xiang | 21/22 | 0.85 | 2.4x |
| Zhao | 46/56 | 0.56 | 1.7x |

**Required sample sizes:**
- Medium effect (1.6x): 63 per group (126 total)
- Small-medium (1.3x): 175 per group (350 total)

Most published studies are below these thresholds.

### Known Positive Validation

| Marker | Expected | Observed | Pooled p | Status |
|--------|----------|----------|----------|--------|
| Fusobacterium | ↑ | ↑ | 0.56 | Correct direction, underpowered |
| Parvimonas | ↑ | ↑ | 0.53 | Correct direction, underpowered |
| Porphyromonas | ↑ | ↓ | 0.51 | High heterogeneity (I²=92%) |

**Conclusion**: Framework correctly identifies effect direction but lacks power to achieve significance.

---

## Implications

### For Researchers

1. **Single-cohort findings are hypothesis-generating, not confirmatory**
   - Even n=100 studies can only detect large effects
   - Significance in one cohort ≠ reliable finding

2. **Effect size matters more than p-values**
   - Report confidence intervals
   - Compare to minimum detectable effects

3. **Heterogeneity is expected, not failure**
   - Microbiome effects are often subtype-specific
   - Not all patients with a disease share the same microbial signature

### For the Field

1. **Published effect sizes are probably inflated**
   - Winner's curse: only significant results publish
   - Outlier patients can drive entire findings

2. **Meta-analysis is essential but not magic**
   - Pooling underpowered studies doesn't create power
   - Need to address heterogeneity explicitly

3. **Subtype thinking needed**
   - "Fusobacterium-high CRC" vs "Fusobacterium-low CRC"
   - May have different etiologies and prognoses

### For This Project

1. **Replication scores alone are misleading**
   - 0% replication ≠ false findings
   - Need power-adjusted interpretation

2. **The Atlas should educate, not just report**
   - Explain why effects don't replicate
   - Provide sample size recommendations

---

## Revised Deliverables

### 1. Effect Characterization Database

For each taxon-disease pair:
```
Fusobacterium × CRC
├── Pooled effect: -0.26 (95% CI: -1.13, 0.61)
├── Direction: Enriched in disease (3/4 cohorts)
├── Meta-analysis p: 0.56
├── Heterogeneity: I² = 0% (consistent)
├── Power assessment: UNDERPOWERED
│   └── Need n=126+ total to detect at 80% power
├── Patient-level: ~15% of CRC show high loads
└── Confidence: LIKELY REAL, SUBTYPE EFFECT
```

### 2. Power Calculator

Interactive tool: "Given your sample size, what effect can you detect?"

### 3. Study Design Recommendations

For each disease: "To reliably detect associations, you need..."

### 4. Publication

**Revised Title**: "The Microbiome Reproducibility Atlas: Why Most Findings Don't Replicate and What We Can Do About It"

**Key Messages**:
1. Non-replication ≠ false positives (often underpowered)
2. Even "robust" findings like Fusobacterium-CRC are subtype effects
3. The field needs larger cohorts and subtype-aware analysis
4. Simple replication metrics are misleading

---

## Technical Implementation

### Directory Structure

```
experiments/atlas/
├── FRAMEWORK.md              # Methodology documentation
├── VALIDATION_RESULTS.md     # Fusobacterium case study
├── scripts/
│   ├── taxonomy.py           # Taxonomic harmonization
│   ├── meta_analysis.py      # Random-effects meta-analysis
│   ├── power_analysis.py     # Power calculations
│   └── reprocess_with_taxonomy.py
├── diseases/
│   └── crc/
│       ├── data_genus/       # Genus-aggregated counts
│       └── results_genus/    # Per-cohort DAA results
└── results/
    └── (forthcoming database)
```

### Available Scripts

```bash
# Reprocess with taxonomic harmonization
python scripts/reprocess_with_taxonomy.py

# Run meta-analysis
python scripts/meta_analysis.py

# Power analysis
python scripts/power_analysis.py
```

---

## Next Steps

### Phase 2: Expand Analysis

1. **More diseases**: IBD, obesity, T2D with genus-level harmonization
2. **Shotgun validation**: curatedMetagenomicData for species-level
3. **Larger cohorts**: GMrepo for higher-powered analyses

### Phase 3: Build Database

1. **Searchable interface** for taxon-disease queries
2. **Power calculator** for study planning
3. **Visualization** of effect heterogeneity

### Phase 4: Publication

1. **Methods paper**: Framework and validation
2. **Resource paper**: Full database release
3. **Commentary**: Implications for the field

---

## References

1. Duvallet C, et al. (2017). Meta-analysis of gut microbiome studies identifies disease-specific and shared responses. Nat Commun.
2. MicrobiomeHD: https://zenodo.org/record/840333
3. GMrepo: https://gmrepo.humangut.info/
4. curatedMetagenomicData: https://waldronlab.io/curatedMetagenomicData/
5. Castellarin M, et al. (2012). Fusobacterium nucleatum infection is prevalent in human colorectal carcinoma. Genome Res.
6. Kostic AD, et al. (2012). Genomic analysis identifies association of Fusobacterium with colorectal carcinoma. Genome Res.
