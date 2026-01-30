# Preliminary Results

Raw findings from exploratory analysis performed during toolkit development. These will be formalized into proper analyses when we return to the experiments.

## Date: 2025-01-29

### Ravel BV Analysis

**Data**: 345 samples (248 healthy, 97 BV), 247 taxa

**LinDA Results**:
```
Total features tested: 55
Significant (q<0.05): 47

Direction of effects:
  Higher in healthy: 28 (59.6%)
  Lower in healthy: 19 (40.4%)

Sum of log fold changes: 0.17
Mean log fold change: 0.00
```

**Top Hits**:
```
Lactobacillus crispatus: estimate=4.442, q=0.0000
Lactobacillus jensenii: estimate=3.140, q=0.0000
Prevotella: estimate=-3.858, q=0.0000
Megasphaera: estimate=-3.629, q=0.0000
Sneathia: estimate=-3.051, q=0.0000
```

**Permutation Tests**: 0/5 false positives (FPR well-calibrated)

---

### Stammler Spike-in Analysis

**Data**: 17 samples, 4036 taxa, 3 spike-in organisms

**Spike-in Feature IDs**:
```
Salinibacter ruber: AF323500XXXX (constant spike)
Alicyclobacillus acidiphilus: AB076660XXXX
Rhizobium radiobacter: AB247615XXXX
```

**Spike-in Variation**:
```
Salinibacter (constant spike):
  Mean proportion: 7.112%
  Range: 1.622% - 15.938%
  CV: 64.7%

Alicyclobacillus:
  Mean proportion: 0.610%
  Range: 0.163% - 1.435%
  CV: 61.5%

Rhizobium:
  Mean proportion: 34.268%
  Range: 6.285% - 82.663%
  CV: 69.6%
```

**Estimated Total Load** (from Salinibacter):
```
Load range: 0.27x - 2.70x of average
Variation: 9.8x across samples
Log2FC artifact potential: ±3.3
```

---

### Key Quantitative Findings

1. **9.8x variation in total bacterial load** (Stammler)
2. **±3.3 log2FC artifact potential** from load alone
3. **Sum of log2FC = 0** in BV analysis (compositional closure)
4. **47/55 taxa significant** in BV (85%) - but interpretation uncertain

---

### Commands Used

```bash
# Fetch datasets
daa fetch stammler-spikein
daa fetch ravel-bv

# Run BV analysis
daa linda \
    -c /tmp/ravel/counts_hbv_correct.tsv \
    -m /tmp/ravel/metadata_hbv.tsv \
    --formula "~ study_condition" \
    --test-coef study_conditionhealthy \
    --output /tmp/ravel/linda_results.tsv

# Data location
# Ravel: /tmp/ravel/
# Stammler: /tmp/stammler/
```

---

### Files Generated

```
/tmp/ravel/
  counts_hbv_correct.tsv    # Filtered counts (H vs BV only)
  metadata_hbv.tsv          # Filtered metadata
  linda_results.tsv         # LinDA output
  metadata_perm_*.tsv       # Permutation test metadata

/tmp/stammler/
  counts.tsv
  metadata.tsv
  taxonomy.tsv              # Maps feature IDs to taxonomy
```
