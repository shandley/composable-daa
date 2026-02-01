# Microbiome Reproducibility Atlas - Rigorous Framework

## Overview

A methodologically rigorous framework for cross-cohort meta-analysis that addresses the key limitations of naive OTU-level comparison:

1. **Taxonomic Harmonization** - Map disparate OTU IDs to common taxonomy
2. **Power Analysis** - Exclude or weight underpowered cohorts appropriately
3. **Formal Meta-Analysis** - Random-effects models with heterogeneity assessment
4. **Shotgun Validation** - Extend framework to metagenomic data
5. **Known Positive Validation** - Verify detection of established associations

---

## 1. Taxonomic Harmonization

### The Problem

Different studies use different:
- Clustering thresholds (97% vs 99% OTUs vs ASVs)
- Reference databases (Greengenes, SILVA, RDP)
- 16S regions (V1-V3, V3-V4, V4, V4-V5)
- Taxonomic classifiers (RDP, BLAST, exact matching)

Result: The same organism gets different identifiers across studies.

### The Solution

#### Level 1: Genus-Level Aggregation
Simplest approach - aggregate OTUs to genus level before comparison.

```
OTU_123 (Fusobacterium nucleatum) → Fusobacterium
OTU_456 (Fusobacterium sp.) → Fusobacterium
denovo_789 (Fusobacterium) → Fusobacterium
```

Pros: Simple, robust to sequence variation
Cons: Loses species-level resolution

#### Level 2: Taxonomic String Matching
Parse RDP/SILVA classifications and match at each level.

```python
def harmonize_taxonomy(otu_name):
    """Extract standardized taxonomy from various formats."""
    # Handle: "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium"
    # Handle: "Fusobacterium nucleatum"
    # Handle: "denovo123_Fusobacterium"
    # Return: {"kingdom": "Bacteria", "phylum": "Fusobacteria", ..., "genus": "Fusobacterium"}
```

#### Level 3: Sequence-Based Matching (Gold Standard)
Use representative sequences to match OTUs across studies.

```bash
# Cluster all representative sequences at 97%
vsearch --cluster_fast all_rep_seqs.fasta --id 0.97 --centroids centroids.fasta --uc clusters.uc

# Map original OTUs to unified clusters
# OTU_123 (study1) → Cluster_001
# denovo_456 (study2) → Cluster_001
```

Pros: Most accurate
Cons: Requires sequence data (not always available in processed datasets)

### Implementation

```python
class TaxonomyHarmonizer:
    def __init__(self, level='genus'):
        self.level = level  # 'genus', 'family', 'species'
        self.mapping = {}

    def parse_taxonomy_string(self, name):
        """Parse various taxonomy formats."""
        # k__X;p__Y;c__Z format
        if ';' in name and '__' in name:
            return self._parse_greengenes(name)
        # Plain name
        return self._parse_plain(name)

    def harmonize(self, taxa_list):
        """Map list of taxa to standardized names at specified level."""
        return [self.parse_taxonomy_string(t).get(self.level, 'Unknown') for t in taxa_list]
```

---

## 2. Power Analysis

### The Problem

Small cohorts (n < 30 per group) may lack power to detect real effects, leading to false negatives that inflate "non-replication."

### The Solution

#### Pre-Analysis Power Calculation

For each cohort, estimate detectable effect size:

```python
def minimum_detectable_effect(n_case, n_control, alpha=0.05, power=0.80):
    """Calculate minimum detectable log2 fold change."""
    from scipy import stats

    # Two-sample t-test power calculation
    # For 80% power at alpha=0.05 with n1=n2=20: d ≈ 0.9 (large effect)
    # For 80% power at alpha=0.05 with n1=n2=50: d ≈ 0.57 (medium effect)

    n_total = n_case + n_control
    harmonic_n = 2 * n_case * n_control / n_total

    # Effect size (Cohen's d) detectable at given power
    ncp = stats.norm.ppf(1 - alpha/2) + stats.norm.ppf(power)
    d = ncp / np.sqrt(harmonic_n / 2)

    return d  # Approximate log2 fold change
```

#### Power-Based Cohort Weighting

```python
def cohort_weight(n_case, n_control, target_effect=1.0):
    """Weight cohorts by their power to detect target effect."""
    mde = minimum_detectable_effect(n_case, n_control)
    if mde > target_effect:
        # Cohort cannot detect the target effect
        return 0.0
    else:
        # Weight by effective sample size
        return 2 * n_case * n_control / (n_case + n_control)
```

#### Interpretation Guidelines

| n per group | Minimum Detectable Effect | Use For |
|-------------|---------------------------|---------|
| < 15 | > 1.5 log2FC (3x) | Exclude or discovery only |
| 15-30 | 1.0-1.5 log2FC | Large effects only |
| 30-50 | 0.7-1.0 log2FC | Medium-large effects |
| 50-100 | 0.5-0.7 log2FC | Medium effects |
| > 100 | < 0.5 log2FC | Small effects |

---

## 3. Formal Meta-Analysis

### The Problem

"Vote counting" (counting how many studies find significance) ignores:
- Effect size magnitude
- Confidence intervals
- Study heterogeneity
- Direction of effects

### The Solution

#### Random-Effects Meta-Analysis

For each taxon present in 2+ cohorts:

```python
def random_effects_meta(effects, variances):
    """
    Random-effects meta-analysis using DerSimonian-Laird estimator.

    Parameters:
    - effects: list of effect sizes (log2FC or similar)
    - variances: list of variance estimates (SE^2)

    Returns:
    - pooled_effect: weighted mean effect
    - pooled_se: standard error of pooled effect
    - tau2: between-study variance
    - I2: heterogeneity statistic
    - Q: Cochran's Q
    """
    import numpy as np

    effects = np.array(effects)
    variances = np.array(variances)
    weights = 1 / variances

    # Fixed-effect estimate (for Q calculation)
    fe_effect = np.sum(weights * effects) / np.sum(weights)

    # Cochran's Q
    Q = np.sum(weights * (effects - fe_effect)**2)
    df = len(effects) - 1

    # DerSimonian-Laird tau^2
    C = np.sum(weights) - np.sum(weights**2) / np.sum(weights)
    tau2 = max(0, (Q - df) / C)

    # Random-effects weights
    re_weights = 1 / (variances + tau2)

    # Pooled estimate
    pooled_effect = np.sum(re_weights * effects) / np.sum(re_weights)
    pooled_var = 1 / np.sum(re_weights)
    pooled_se = np.sqrt(pooled_var)

    # I^2 heterogeneity
    I2 = max(0, (Q - df) / Q) if Q > 0 else 0

    return {
        'effect': pooled_effect,
        'se': pooled_se,
        'ci_lower': pooled_effect - 1.96 * pooled_se,
        'ci_upper': pooled_effect + 1.96 * pooled_se,
        'p_value': 2 * (1 - stats.norm.cdf(abs(pooled_effect / pooled_se))),
        'tau2': tau2,
        'I2': I2,
        'Q': Q,
        'Q_pvalue': 1 - stats.chi2.cdf(Q, df),
        'n_studies': len(effects)
    }
```

#### Interpretation

| I² | Heterogeneity | Interpretation |
|----|---------------|----------------|
| 0-25% | Low | Consistent across studies |
| 25-50% | Moderate | Some variation |
| 50-75% | Substantial | Consider subgroup analysis |
| > 75% | High | Studies measuring different things |

#### Forest Plot Output

```
Taxon: Fusobacterium
─────────────────────────────────────────────────────────────────
                              Effect [95% CI]          Weight
Baxter 2016    ────●────      1.23 [0.45, 2.01]        24.3%
Zeller 2014    ─────●───      1.56 [0.89, 2.23]        28.1%
Zackular 2014  ────●────      0.98 [0.12, 1.84]        18.7%
Zhao 2012      ───●─────      1.42 [0.67, 2.17]        28.9%
─────────────────────────────────────────────────────────────────
Pooled (RE)    ────◆────      1.32 [0.91, 1.73]
─────────────────────────────────────────────────────────────────
Heterogeneity: I² = 12%, Q = 3.4, p = 0.34
Test for effect: z = 6.3, p < 0.0001
```

---

## 4. Shotgun Validation

### The Problem

16S has limited taxonomic resolution. Shotgun metagenomics provides:
- Species/strain-level identification
- Functional gene content
- Better quantification

### The Solution

#### Data Source: curatedMetagenomicData

```r
# Bioconductor package with standardized shotgun data
library(curatedMetagenomicData)

# Get CRC studies
crc_data <- curatedMetagenomicData("*CRC*", dryrun = FALSE)
```

#### Equivalent Analysis Pipeline

```yaml
# Shotgun version of our pipeline
name: shotgun_linda
steps:
- !FilterPrevalence
  threshold: 0.10
- !AddPseudocount
  value: 0.5
- NormalizeCLR
- !ModelLM
  formula: "~ disease + age + sex"  # Can include covariates!
- !TestWald
  coefficient: diseaseCRC
- CorrectBH
```

#### Advantages

1. Species-level taxonomy (no harmonization needed)
2. Consistent processing (MetaPhlAn)
3. Rich metadata for covariate adjustment
4. Functional validation possible

---

## 5. Known Positive Validation

### The Concept

Before claiming "nothing replicates," verify we can detect established associations.

### Case Study: Fusobacterium in CRC

**Known biology:**
- Fusobacterium nucleatum is enriched in CRC tumors
- Detected in tissue, stool, and blood
- Mechanistic studies show it promotes tumorigenesis
- One of the most replicated findings in CRC microbiome research

**Expected finding:**
- Fusobacterium should be significantly enriched in CRC across most/all cohorts
- If we don't find it, our methodology has a problem

### Validation Framework

```python
KNOWN_POSITIVES = {
    'crc': [
        {'taxon': 'Fusobacterium', 'direction': 'up', 'confidence': 'high'},
        {'taxon': 'Parvimonas', 'direction': 'up', 'confidence': 'high'},
        {'taxon': 'Porphyromonas', 'direction': 'up', 'confidence': 'medium'},
        {'taxon': 'Peptostreptococcus', 'direction': 'up', 'confidence': 'medium'},
    ],
    'ibd': [
        {'taxon': 'Faecalibacterium', 'direction': 'down', 'confidence': 'high'},
        {'taxon': 'Roseburia', 'direction': 'down', 'confidence': 'high'},
        {'taxon': 'Enterobacteriaceae', 'direction': 'up', 'confidence': 'high'},
    ],
    'cdi': [
        {'taxon': 'Clostridioides', 'direction': 'up', 'confidence': 'high'},
    ],
}

def validate_known_positives(results, disease):
    """Check if known associations are detected."""
    expected = KNOWN_POSITIVES.get(disease, [])

    for kp in expected:
        taxon = kp['taxon']
        direction = kp['direction']

        # Search for taxon in results (genus-level match)
        matches = results[results.index.str.contains(taxon, case=False)]

        if len(matches) == 0:
            print(f"WARNING: {taxon} not found in results")
            continue

        # Check significance and direction
        best_match = matches.loc[matches['q_value'].idxmin()]
        is_sig = best_match['q_value'] < 0.05
        is_correct_dir = (best_match['estimate'] > 0) == (direction == 'up')

        status = "✓" if (is_sig and is_correct_dir) else "✗"
        print(f"{status} {taxon}: q={best_match['q_value']:.3f}, effect={best_match['estimate']:.2f}")
```

### If Validation Fails

If we don't detect Fusobacterium in CRC:
1. **Check taxonomy** - Is it present under a different name?
2. **Check prevalence** - Is it filtered out by prevalence threshold?
3. **Check power** - Are cohorts too small?
4. **Check methodology** - Is our pipeline biased?

This is the "sanity check" that validates our entire framework.

---

## Implementation Plan

### Phase 1: Fusobacterium Validation (Today)
- Search for Fusobacterium in CRC results
- If not found, diagnose why
- Implement genus-level aggregation if needed

### Phase 2: Taxonomic Harmonization
- Build genus-level aggregator
- Re-run CRC meta-analysis at genus level
- Compare to OTU-level results

### Phase 3: Power Analysis
- Calculate power for each cohort
- Weight or exclude underpowered cohorts
- Re-assess replication with power adjustment

### Phase 4: Formal Meta-Analysis
- Implement random-effects meta-analysis
- Generate forest plots
- Calculate heterogeneity statistics

### Phase 5: Shotgun Validation
- Download curatedMetagenomicData CRC studies
- Run equivalent analysis
- Compare 16S vs shotgun results

---

## Success Criteria

The framework is validated when:

1. **Fusobacterium is detected** as enriched in CRC across cohorts
2. **Genus-level replication** is higher than OTU-level
3. **Meta-analysis** shows significant pooled effect with low-moderate heterogeneity
4. **Shotgun data** confirms 16S findings at species level
5. **Known positives** are detected for IBD and CDI as well
