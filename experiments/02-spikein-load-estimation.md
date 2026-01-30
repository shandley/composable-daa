# Experiment 02: Spike-in Load Estimation

## Summary

Using the Stammler 2016 spike-in dataset to demonstrate that total bacterial load varies dramatically across samples, and that this variation is sufficient to explain many reported effect sizes in the microbiome literature.

## Background

Stammler et al. (2016) added known quantities of three bacterial species to fecal samples before DNA extraction:

1. **Salinibacter ruber** - Constant amount across all samples
2. **Rhizobium radiobacter** - Variable amounts (simulating differential abundance)
3. **Alicyclobacillus acidiphilus** - Variable amounts

Because Salinibacter was added at constant absolute abundance, its relative abundance reflects the inverse of total bacterial load. Higher Salinibacter proportion = lower total bacteria.

## Data

**Source**: MicrobiomeBenchmarkData (Zenodo 6911027)
**Dataset**: Stammler 2016 spike-in
**Samples**: 17
**Features**: 4036

```bash
# Fetch data
daa fetch stammler-spikein
```

## Analysis Performed

### 1. Spike-in Identification

| Organism | Feature ID | Role |
|----------|-----------|------|
| Salinibacter ruber | AF323500XXXX | Constant spike (internal standard) |
| Alicyclobacillus acidiphilus | AB076660XXXX | Variable spike |
| Rhizobium radiobacter | AB247615XXXX | Variable spike |

### 2. Spike-in Proportion Variation

| Spike-in | Mean % | Range % | CV |
|----------|--------|---------|-----|
| Salinibacter | 7.1% | 1.6% - 15.9% | 64.7% |
| Alicyclobacillus | 0.6% | 0.2% - 1.4% | 61.5% |
| Rhizobium | 34.3% | 6.3% - 82.7% | 69.6% |

**Key observation**: All spike-ins show ~65% CV in relative abundance, despite being added at known amounts. This variation reflects sample-to-sample differences in total bacterial load.

### 3. Estimated Total Load

Using Salinibacter as internal standard:

```
Estimated_Load = 1 / Salinibacter_Proportion
```

**Results**:
- Load range: 0.27x to 2.70x of average
- **9.8x variation** in total bacterial load across samples
- Log2 artifact potential: **±3.3**

### 4. Spike-in Correlation

All three spike-ins are highly correlated (r > 0.9) despite being added at different amounts. This correlation arises because all are diluted by the same total load factor.

```
Salinibacter ~ Alicyclobacillus: r = 0.90
Salinibacter ~ Rhizobium: r = 0.85
Alicyclobacillus ~ Rhizobium: r = 0.88
```

## Key Insight

**A 9.8x variation in total load creates ±3.3 log2FC artifacts.**

This is comparable to or larger than most reported effect sizes in microbiome differential abundance studies. Many published findings could be entirely explained by total load differences between groups.

## Implications

### For Study Design
- Include internal standards (spike-ins) whenever possible
- Measure total load independently (qPCR, flow cytometry)
- Report both relative and absolute abundances

### For Interpretation
- Effect sizes <3 log2FC could be load artifacts
- Bidirectional changes (some taxa up, others down) are expected under compositional closure
- "Significant" findings may reflect load, not taxa-specific biology

### For Meta-analysis
- Studies without load correction cannot be directly compared
- Heterogeneity across studies may reflect load variation
- Need standardized spike-in or absolute quantification protocols

## Spike-in Normalization Method

We implemented `norm_spikein()` in the toolkit:

```rust
use composable_daa::norm_spikein;

// Normalize counts using spike-in as internal standard
let normalized = norm_spikein(&counts, "AF323500XXXX")?;

// Access estimated total load
let load = normalized.estimated_load;
```

## Next Steps

1. **Validate on additional spike-in datasets**
2. **Compare spike-in normalization to other methods** (RLE, TMM, CSS)
3. **Develop guidelines** for spike-in experimental design
4. **Integrate with differential abundance** to produce load-corrected results

## Code Location

Analysis scripts will be added to `experiments/scripts/02-spikein-analysis/` when formalized.

## References

- Stammler et al. (2016) - Original spike-in methodology paper
- Vandeputte et al. (2017) Nature - QMP method using flow cytometry
- Tkacz et al. (2018) - Spike-in protocols for plant microbiome
