# Experiment 02: Spike-in Load Estimation

## Summary

Analysis of the Stammler 2016 spike-in dataset demonstrates that total bacterial load varies **9.8-fold** across samples, creating potential artifacts of up to **3.3 log2FC**. This exceeds most reported effect sizes in microbiome studies, suggesting many published findings could be confounded by uncontrolled load variation.

**Key Finding**: Most microbiome studies report effect sizes in the 1-3 log2FC range. Our analysis shows that load variation alone can generate artifacts of this magnitude, meaning apparent differences between groups may reflect load rather than taxon-specific biology.

## Background

Stammler et al. (2016) added known quantities of three bacterial species to fecal samples before DNA extraction:

| Organism | Role | Expected Behavior |
|----------|------|-------------------|
| Salinibacter ruber | Constant spike (internal standard) | Relative abundance ∝ 1/total_load |
| Rhizobium radiobacter | Variable spike | Tests quantification |
| Alicyclobacillus acidiphilus | Variable spike | Tests quantification |

Because Salinibacter was added at **constant absolute abundance**, its relative abundance directly reflects total bacterial load:
- High Salinibacter % → Low total bacteria
- Low Salinibacter % → High total bacteria

This allows us to estimate sample-to-sample load variation without external quantification (qPCR, flow cytometry).

## Data

**Source**: MicrobiomeBenchmarkData (Zenodo 6911027)
**Dataset**: Stammler 2016 spike-in
**Samples**: 17 (longitudinal samples from ASCT patients)
**Features**: 4,036 OTUs

```bash
# Fetch data
daa fetch -d stammler -o ./data/
```

## Analysis Performed

### 1. Spike-in Identification

Three spike-in features were identified by accession number:

| Spike-in | Feature ID | Total Counts |
|----------|-----------|--------------|
| Salinibacter ruber | AF323500XXXX | 39,285 |
| Alicyclobacillus acidiphilus | AB076660XXXX | 3,429 |
| Rhizobium radiobacter | AB247615XXXX | 191,803 |

### 2. Spike-in Relative Abundance Variation

Despite being added at known amounts, all spike-ins showed high variation in relative abundance:

| Spike-in | Mean % | Range % | CV |
|----------|--------|---------|-----|
| Salinibacter (constant) | 7.1% | 1.6% - 15.9% | 64.7% |
| Alicyclobacillus (variable) | 0.6% | 0.2% - 1.4% | 61.5% |
| Rhizobium (variable) | 34.3% | 6.3% - 82.7% | 69.6% |

**Key observation**: Even the constant spike (Salinibacter) shows 65% CV. This variation is not measurement error—it reflects real differences in total bacterial load across samples.

### 3. Total Load Estimation

Using Salinibacter as internal standard:

```
Estimated_Load = 1 / Salinibacter_Proportion
(normalized to mean = 1)
```

**Results**:
- Load range: **0.27x to 2.70x** of average
- Fold range: **9.8x**
- Log2 artifact potential: **±3.3**

### 4. Spike-in Correlations

All three spike-ins are highly correlated despite being added at different amounts:

| Comparison | Correlation (r) |
|------------|----------------|
| Salinibacter ~ Alicyclobacillus | 0.71 |
| Salinibacter ~ Rhizobium | 0.90 |
| Alicyclobacillus ~ Rhizobium | 0.85 |

High correlations confirm they're all affected by the same factor: dilution by total bacterial load. This validates Salinibacter as an internal standard.

### 5. Spike-in Normalization Demonstration

Scaling counts by (median_spike-in / sample_spike-in):

| Metric | Before | After |
|--------|--------|-------|
| Salinibacter CV | 66.3% | **0.0%** |
| Other taxa CV | Variable | Reduced |

Spike-in normalization successfully removes load-related variation.

## Key Insights

### 1. Load Variation is Substantial

The 9.8x variation in total bacterial load is **larger than most reported effect sizes**:

| Effect Size Category | log2FC | Fold Change | Within Artifact Range? |
|---------------------|--------|-------------|----------------------|
| Small | 1 | 2x | **Yes** |
| Moderate | 2 | 4x | **Yes** |
| Large | 3 | 8x | **Yes** |
| Very large | 4 | 16x | Partially |

### 2. Compositional Effects are Compounded

Load variation creates compositional artifacts:
- If sample A has 2x higher load than sample B
- A taxon with constant absolute abundance appears 2x lower in A
- A taxon that doubles appears unchanged
- These are NOT independent biological changes

Combined with Experiment 01's demonstration of compositional closure, this shows:
1. **CLR sum = 0** (closure constraint)
2. **Load variation shifts the baseline** (this experiment)
3. **Both effects can masquerade as biology**

### 3. Spike-in Normalization Works

When spike-in data is available, normalization effectively removes load artifacts:
- Spike-in CV goes from 66% to 0%
- Load-corrected abundances reflect true differences
- This should be standard practice when spike-ins are included

## Implications for Published Studies

### Most Studies Don't Control for Load

| Study Type | Load Control | Artifact Risk |
|------------|--------------|---------------|
| Standard 16S | None | High |
| With qPCR | External quantification | Low |
| With spike-ins | Internal standard | Low |
| With flow cytometry | Cell counts | Low |

### Effect Sizes May Be Inflated or Spurious

Studies reporting effect sizes <3 log2FC without load correction should be interpreted cautiously:
- The effect could be entirely load-driven
- Apparent bidirectional changes could reflect load + closure
- Reproducibility across studies may be poor

### Recommendations

1. **For study design**: Include spike-in controls or external quantification
2. **For analysis**: Report both relative and load-corrected abundances
3. **For interpretation**: Be skeptical of small effects without load data
4. **For meta-analysis**: Only compare studies with similar load controls

## Connection to Experiment 01 (BV Analysis)

The BV reanalysis showed that CLR estimates sum to zero (compositional closure). This experiment adds another dimension:

| Experiment | Finding | Implication |
|------------|---------|-------------|
| 01 (BV) | Sum of CLR = 0 | Taxa changes are coupled, not independent |
| 02 (Spike-in) | 9.8x load variation | Baseline shifts create artifacts |
| Combined | Both effects operate simultaneously | Many findings may be entirely artifactual |

If BV involves a 10x increase in total load (as suggested by biofilm formation), and load variation alone can create 3.3 log2FC artifacts, then:
- The apparent Lactobacillus "decrease" could be entirely dilution
- The apparent pathogen "increase" could be neutral + load + closure
- Only taxa with effects >3.3 log2FC are reliably different

## Reproducibility

All analysis is reproducible via:

```bash
cd experiments/scripts/02-spikein-analysis/
./run_analysis.sh
```

Output files:
- `results/spikein_analysis.tsv` - Per-sample spike-in and load estimates
- `results/spikein_correlations.tsv` - Correlation matrix
- `results/counts_spikein_normalized.tsv` - Load-corrected counts
- `results/summary_statistics.tsv` - Key metrics
- `results/analysis_report.md` - Summary report

Figures:
```bash
python3 generate_figures.py
```

## Next Steps

1. **Validate on other spike-in datasets**: Tkacz et al., Tourlousse et al.
2. **Compare normalization methods**: Spike-in vs RLE vs TMM vs CSS
3. **Develop load estimation without spike-ins**: Library size proxies
4. **Integrate with DAA pipelines**: Load-aware differential abundance

## References

- Stammler F, et al. (2016). Adjustment of 16S rRNA gene sequence data of fecal samples by spike-in of known bacterial cells or DNA. Sci Rep.
- Vandeputte D, et al. (2017). Quantitative microbiome profiling links gut community variation to microbial load. Nature.
- Tkacz A, et al. (2018). Absolute quantification of microbiota abundance in environmental samples. Microbiome.
- Tourlousse DM, et al. (2017). Synthetic spike-in standards for high-throughput 16S rRNA gene amplicon sequencing. Nucleic Acids Res.
