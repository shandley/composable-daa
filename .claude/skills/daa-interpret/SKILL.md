---
name: daa-interpret
description: Interpret differential abundance analysis results. Use when user has DAA output and wants to understand significant features, effect sizes, confidence levels, or next steps.
argument-hint: "[results.tsv]"
allowed-tools: Read, Glob, Bash
---

# DAA Results Interpretation Workflow

This skill helps users interpret differential abundance analysis results from the `daa` CLI.

## Step 1: Load Results

If the user provided a results file path, use the Read tool to load it directly.

If no file was provided, search for recent results:

```bash
ls -lt *.tsv *results*.tsv 2>/dev/null | head -10
```

Then use the Read tool to load the identified results file.

## Step 2: Parse Results Structure

The results TSV has these standard columns:

| Column | Description |
|--------|-------------|
| `feature_id` | Feature identifier |
| `coefficient` | Which coefficient was tested (e.g., `grouptreatment`) |
| `estimate` | Effect size estimate |
| `std_error` | Standard error of estimate |
| `statistic` | Test statistic (t or z) |
| `p_value` | Raw p-value |
| `q_value` | FDR-corrected q-value |
| `prevalence` | Proportion of samples with feature |
| `mean_abundance` | Mean CLR or log abundance |
| `prevalence_tier` | very_high, high, moderate, low, rare |
| `confidence` | high, moderate, suggestive, not_significant |

### Identify the Method Used

Check the output structure:

| Column Pattern | Method |
|---------------|--------|
| `estimate`, `std_error`, `statistic` (simple) | LinDA (LM) or LMM |
| `count_estimate`, `zero_estimate` | Hurdle model |
| `mu_estimate`, `zi_estimate` | ZINB |
| `nb_estimate`, `dispersion` | Negative binomial |

## Step 3: Calculate Summary Statistics

From the loaded data, calculate:

```
Total features tested: count all rows
Significant at q < 0.05: count where q_value < 0.05
Significant at q < 0.10: count where q_value < 0.10

Direction (based on estimate sign):
- Positive estimate = Up in treatment/target group
- Negative estimate = Down in treatment/target group
```

### Apply Method-Specific Thresholds

**Critical**: Different methods require different q-value thresholds based on empirical benchmarks:

| Method | Recommended Threshold | Reason |
|--------|----------------------|--------|
| LinDA/LMM | q < 0.10 | CLR attenuation reduces power |
| Hurdle | q < 0.05 | Standard threshold works well |
| ZINB | q < 0.05 | Standard threshold works well |
| NB | q < 0.05 | Standard threshold works well |

## Step 4: Interpret Effect Sizes

### For LinDA/LMM (CLR-transformed)

CLR effects are attenuated by ~75%. To estimate true fold change:

```
Approximate true log2FC = CLR_estimate × 4
True fold change = 2^(CLR_estimate × 4)
```

| CLR Estimate | Approx True FC | Interpretation |
|--------------|----------------|----------------|
| 0.25 | ~2x | Small effect |
| 0.50 | ~4x | Moderate effect |
| 0.75 | ~8x | Large effect |
| 1.00 | ~16x | Very large effect |
| 1.50 | ~64x | Check for artifacts |

### For ZINB/Hurdle/NB (Count Models)

Estimates are in natural log scale:

```
Fold change = exp(estimate)
```

| Estimate | Fold Change | Interpretation |
|----------|-------------|----------------|
| 0.69 | 2x | Small effect |
| 1.10 | 3x | Moderate effect |
| 1.39 | 4x | Large effect |
| 2.30 | 10x | Very large effect |
| 4.61 | 100x | Check for artifacts |

## Step 5: Detect Compositional Artifacts

**This is critical for quality assessment.** Compositional data can produce spurious associations.

### Check 1: Dominant Taxa Bias

Sort features by `mean_abundance` (highest first). Check if top abundant features are disproportionately significant:

```
IF >50% of top 10 most abundant features are significant:
    → WARNING: Possible compositional artifact
    → Dominant taxa shifts may cascade through all features
    → Recommend: Run `daa stress` to quantify
```

### Check 2: Direction Imbalance

Count significant features by direction:

```
IF >80% of significant features go same direction:
    → WARNING: Possible compositional artifact
    → A dominant feature changing may push others opposite direction
    → Recommend: Check if a single dominant taxon is driving results
```

### Check 3: Effect Size Correlation with Abundance

```
IF negative correlation between significance and abundance:
    → May indicate real biological signal (rare taxa responding)
IF positive correlation between significance and abundance:
    → May indicate compositional artifact
```

### Check 4: Extreme Effect Sizes

```
IF any effect >100x fold change:
    → WARNING: Suspicious effect size
    → Check raw data for outliers or zero-inflation issues
```

## Step 6: Cross-Reference with Data Profile (if available)

If a data profile is available (from `daa profile-llm`), cross-reference:

| Profile Metric | Results Check |
|----------------|---------------|
| High sparsity (>70%) | Are significant features less sparse? |
| Unbalanced groups | Did smaller group drive significance? |
| High library size CV | Are effects confounded with depth? |
| Batch present | Were batch effects controlled? |

Run profile if not available:
```bash
daa profile-llm -c {COUNTS_FILE} -m {METADATA_FILE} -g {GROUP_COLUMN}
```

## Step 7: Generate Interpretation Report

### Template Output

```
## Results Interpretation

**Method**: {method}
**Threshold applied**: q < {threshold}
**Features tested**: {n_total}
**Significant features**: {n_sig}
**Up in {target_group}**: {n_up}
**Down in {target_group}**: {n_down}

### Top Significant Features

| Feature | Effect | Fold Change | q-value | Prevalence | Confidence |
|---------|--------|-------------|---------|------------|------------|
| {feature1} | {estimate} | {fc}x {dir} | {q} | {prev}% | {conf} |
| ... | ... | ... | ... | ... | ... |

### Effect Size Distribution

- Large effects (>4x): {n_large}
- Moderate effects (2-4x): {n_moderate}
- Small effects (<2x): {n_small}

### Prevalence Distribution of Significant Features

- Very high (>75%): {n_vhigh}
- High (50-75%): {n_high}
- Moderate (25-50%): {n_moderate}
- Low (<25%): {n_low}

### Quality Assessment

{quality_assessment based on artifact checks}

### Biological Interpretation

{Brief interpretation of what the results suggest biologically}
```

## Step 8: Recommend Validation Steps

Based on results, recommend specific validation commands:

### If Many Significant Features (>10)

```bash
# Validate with spike-in testing
daa validate -c {counts} -m {metadata} -g {group} -t {target} -f "{formula}" --test-coef {coefficient}
```

**Why**: Confirms the method has expected sensitivity/FDR on your data structure.

### If Compositional Concerns Detected

```bash
# Run compositional stress test
daa stress -c {counts} -m {metadata} -g {group} -t {target} -f "{formula}" --test-coef {coefficient}
```

**Why**: Quantifies how much dominant taxa shifts affect other features.

### If Few/No Significant Features

1. **For LinDA**: Was q < 0.10 used? (Required due to CLR attenuation)
2. **Consider alternative method**:
   ```bash
   daa recommend -c {counts} -m {metadata} -g {group} -t {target} --run
   ```
3. **Check power**:
   - Sample size <20/group has very limited power
   - Only >4x effects detectable with n<20

### If Effect Sizes Seem Extreme

```bash
# Profile data to check for outliers
daa profile-llm -c {counts} -m {metadata} -g {group}
```

### For Cross-Validation

```bash
# Run with alternative method for comparison
daa recommend -c {counts} -m {metadata} -g {group} -t {target} --yaml -o alt_pipeline.yaml
# Edit to use different method, then:
daa run -c {counts} -m {metadata} --config alt_pipeline.yaml -o alt_results.tsv
```

## Step 9: Final Summary

Provide a clear, actionable summary:

```
## Summary

Your analysis identified {n_sig} differentially abundant features at q < {threshold}.

**Key findings**:
- {top_feature_1}: {effect}x {direction} (q = {q})
- {top_feature_2}: {effect}x {direction} (q = {q})
- {top_feature_3}: {effect}x {direction} (q = {q})

**Quality assessment**: {Good/Moderate/Concerning}
{If concerning: specific issues identified}

**Recommended next steps**:
1. {specific_recommendation_1}
2. {specific_recommendation_2}
```

## Effect Size Reference

See [effect-sizes.md](effect-sizes.md) for detailed effect size interpretation tables.

## Example Interpretations

### Example 1: Clean Results

**Results**: 8 significant at q < 0.05, mix of up/down, effect sizes 2-8x, various prevalence tiers

**Interpretation**: "Results look clean. Mix of prevalence tiers and balanced directions suggest real biological signal rather than compositional artifacts. Recommend spike-in validation before publication."

### Example 2: Compositional Concern

**Results**: 15 significant, 14 are down, top 3 abundant taxa are all significant

**Interpretation**: "Warning: Strong directional bias and dominant taxa significance suggests possible compositional artifact. One abundant taxon may be driving these associations. Recommend running `daa stress` to quantify compositional effects before interpreting results."

### Example 3: No Significant Results

**Results**: 0 significant at q < 0.05, method was LinDA

**Interpretation**: "No significant features at q < 0.05. However, LinDA requires q < 0.10 due to CLR attenuation. Checking at q < 0.10... [3 features]. Also note sample size (n=15/group) limits power to detect effects <4x."
