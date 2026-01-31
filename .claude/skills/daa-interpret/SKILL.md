---
name: daa-interpret
description: Interpret differential abundance analysis results. Use when user has DAA output and wants to understand significant features, effect sizes, confidence levels, or next steps.
argument-hint: "[results.tsv]"
allowed-tools: Bash, Read, Glob
---

# DAA Results Interpretation Workflow

This skill helps users interpret differential abundance analysis results from the `daa` CLI.

## Step 1: Load Results

If the user provided a results file, read it. Otherwise, find recent results:

```bash
# Find recent results files
ls -lt *.tsv *results*.tsv 2>/dev/null | head -10
```

Read the results file:
```bash
cat {RESULTS_FILE}
```

## Step 2: Identify the Method Used

Look at the output structure to determine which method was used:

| Column Pattern | Method |
|---------------|--------|
| `estimate`, `std_error`, `statistic` | LinDA (linear model) |
| `count_estimate`, `zero_estimate` | Hurdle model |
| `mu_estimate`, `zi_estimate` | ZINB |
| `nb_estimate`, `dispersion` | Negative binomial |

## Step 3: Summarize Key Findings

### Count Significant Features

```
Total features tested: {n_total}
Significant at q < 0.05: {n_sig_05}
Significant at q < 0.10: {n_sig_10}
```

### Categorize by Direction

```
Up in treatment/disease: {n_up}
Down in treatment/disease: {n_down}
```

### Categorize by Confidence

Use the `confidence` column if present:
- **high**: Strong evidence, large effect
- **moderate**: Good evidence, medium effect
- **suggestive**: Weak evidence, worth following up
- **not_significant**: No evidence

## Step 4: Interpret Effect Sizes

### For LinDA (CLR-transformed)

CLR effects are attenuated by ~75%. To estimate true fold change:

```
Approximate true log2FC = CLR_estimate × 4
```

| CLR Estimate | Approx True FC |
|--------------|----------------|
| 0.25 | ~2x |
| 0.50 | ~4x |
| 1.00 | ~16x |
| 1.50 | ~64x |

### For ZINB/Hurdle/NB (Count Models)

Estimates are in log scale. Convert to fold change:

```
Fold change = exp(estimate)
```

| Estimate | Fold Change |
|----------|-------------|
| 0.69 | 2x |
| 1.10 | 3x |
| 1.39 | 4x |
| 2.30 | 10x |

## Step 5: Assess Result Quality

### Warning Signs

1. **All dominant taxa significant**
   - Top features by abundance are all significant
   - May indicate compositional artifacts
   - Recommend: Run `daa stress` to quantify

2. **Only rare taxa significant**
   - Only low-prevalence features significant
   - May be noise or biological signal
   - Recommend: Check prevalence tier

3. **Many false discoveries expected**
   - At q < 0.05 with 100 significant features
   - Expect ~5 false positives
   - Recommend: Use q < 0.01 for high-confidence subset

4. **Effect sizes seem too large**
   - >100x fold changes are suspicious
   - May indicate data issues
   - Recommend: Check raw data for outliers

### Quality Indicators

1. **Good result characteristics**
   - Mix of prevalence tiers
   - Effect sizes 2-20x
   - FDR < 20%
   - Some biological coherence

2. **Concerning result characteristics**
   - All from one prevalence tier
   - Extreme effect sizes (>50x)
   - FDR > 30%
   - No biological pattern

## Step 6: Generate Summary Report

### Template Output

```
## Results Summary

**Method**: {method}
**Features tested**: {n_total}
**Significant (q < 0.05)**: {n_sig}

### Top Significant Features

| Feature | Effect | q-value | Prevalence | Interpretation |
|---------|--------|---------|------------|----------------|
| {feature1} | {effect}x {direction} | {q} | {prev}% | {interp} |
| ... | | | | |

### Effect Size Distribution

- Large effects (>4x): {n_large}
- Moderate effects (2-4x): {n_moderate}
- Small effects (<2x): {n_small}

### Prevalence Distribution

- Very high (>75%): {n_vhigh}
- High (50-75%): {n_high}
- Medium (25-50%): {n_medium}
- Low (<25%): {n_low}

### Quality Assessment

{quality_statement}

### Recommended Next Steps

1. {next_step_1}
2. {next_step_2}
3. {next_step_3}
```

## Step 7: Offer Follow-up Actions

Based on results, offer:

1. **If many significant features**:
   - "Run cross-validation with a second method?"
   - "Generate high-confidence subset at q < 0.01?"

2. **If few/no significant features**:
   - "Try relaxing threshold to q < 0.10 (for LinDA)?"
   - "Run ZINB or Hurdle for better sensitivity?"

3. **If compositional concerns**:
   - "Run `daa stress` to quantify effect attenuation?"

4. **For publication**:
   - "Generate volcano plot data?"
   - "Export for visualization?"

## Effect Size Reference

See [effect-sizes.md](effect-sizes.md) for detailed effect size interpretation.

## Example Interpretation

**User**: "I ran hurdle and got these results"

**Claude**:
1. Reads results file
2. Counts: 8 significant at q < 0.05
3. Categorizes: 5 up, 3 down
4. Interprets effect sizes: range 2.1x to 8.3x
5. Checks quality: mix of prevalence tiers, reasonable FDR
6. Summarizes: "8 features show significant differential abundance..."
7. Offers: "Would you like to validate with spike-in testing?"
