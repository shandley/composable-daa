---
name: daa-diagnose
description: Diagnose microbiome/virome count data and recommend differential abundance analysis methods. Use when user wants to analyze their data, asks which method to use for their specific dataset, or provides count data files for analysis.
argument-hint: "[counts.tsv] [metadata.tsv] [group_column]"
allowed-tools: Bash, Read, Glob
---

# DAA Data Diagnosis Workflow

This skill diagnoses user data and recommends appropriate differential abundance analysis methods based on empirical benchmarks.

## Step 1: Locate Data Files

If the user provided file paths, use those. Otherwise, search for common patterns:

```bash
# Look for count/abundance files
ls -la *.tsv *.csv 2>/dev/null | head -20

# Common naming patterns
ls -la *count* *abundance* *otu* *asv* *feature* 2>/dev/null | head -10
```

Ask the user to confirm which file is the count matrix and which is metadata if unclear.

## Step 2: Profile the Data

Run the LLM-friendly profiler to get structured diagnostics:

```bash
daa profile-llm -c {COUNTS_FILE} -m {METADATA_FILE} -g {GROUP_COLUMN}
```

If the group column is unknown, first inspect the metadata:
```bash
head -5 {METADATA_FILE}
```

## Step 3: Parse Profile Output

The profile-llm command outputs structured sections. Extract key metrics:

### Critical Metrics for Method Selection

| Metric | Where to Find | Decision Impact |
|--------|---------------|-----------------|
| Overall sparsity | `sparsity.overall_zero_fraction` | >70% → Hurdle; 30-70% → ZINB |
| Samples per group | `samples` section | <20/group → warn about power |
| Library size CV | `library_size.cv` | >1.0 → check for confounding |
| Features | `dimensions` | >500 → more stringent correction |
| Group-specific prevalence | `prevalence` section | Varies by group → consider groupwise filtering |

## Step 4: Apply Decision Rules

Based on the profile, apply these evidence-based rules:

### Primary Method Selection

```
IF sparsity > 70%:
    → Recommend HURDLE (best for structural zeros)

ELIF sparsity 50-70%:
    → Recommend ZINB (handles excess zeros well)

ELIF sparsity 30-50%:
    → Recommend ZINB or NB

ELSE (sparsity < 30%):
    → Recommend LinDA or NB (standard methods work)
```

### Sample Size Considerations

```
IF n_per_group < 10:
    → WARN: Very low power, only huge effects detectable
    → Consider pooling or different study design

ELIF n_per_group 10-20:
    → WARN: Limited power, need >4x fold changes
    → Recommend ZINB/Hurdle for discovery

ELIF n_per_group 20-50:
    → Moderate power, >2x fold changes detectable
    → Full method suite available

ELSE (n_per_group > 50):
    → Good power, even moderate effects detectable
```

### Library Size Variation

```
IF library_size_cv > 1.5:
    → WARN: High library size variation
    → Check for batch effects or technical issues
    → Ensure proper normalization

IF library_size differs by group (>2x):
    → WARN: Potential confounding
    → Consider library size as covariate
```

## Step 5: Generate Recommendations

Output a structured recommendation:

### Template Output

```
## Data Diagnosis Summary

**Dataset**: {n_features} features × {n_samples} samples
**Sparsity**: {sparsity}%
**Samples per group**: {group_sizes}
**Library size CV**: {cv}

## Recommended Analysis

### Primary Method: {METHOD}
Rationale: {why this method based on data characteristics}

```bash
daa {method} -c {counts} -m {metadata} -f "~ {group}" -t {group}{level} -o results_{method}.tsv
```

### Threshold: {q_threshold}
- For {METHOD}: use q < {threshold}
{If LinDA: explain q < 0.10 requirement}

### Secondary Validation (Optional)
Run a second method to validate top hits:
```bash
daa {validation_method} -c {counts} -m {metadata} -f "~ {group}" -t {group}{level} -o results_{validation}.tsv
```

## Warnings
{Any warnings about power, confounding, etc.}

## Expected Performance
Based on benchmarks with similar data:
- Sensitivity: ~{X}% at {effect_size} fold change
- FDR: ~{Y}% at q < {threshold}
```

## Step 6: Offer to Run Analysis

After presenting recommendations, offer:

1. "Should I run the recommended analysis now?"
2. "Would you like me to run multiple methods for comparison?"
3. "Do you want to adjust the prevalence filter threshold?"

## Decision Rules Reference

See [decision-rules.md](decision-rules.md) for complete decision tree and benchmark data.

## Example Session

**User**: "Analyze my microbiome data in data/counts.tsv"

**Claude**:
1. Runs `daa profile-llm -c data/counts.tsv -m data/metadata.tsv -g treatment`
2. Sees: 68% sparsity, n=25/group, CV=0.8
3. Recommends: "Based on 68% sparsity and n=25/group, I recommend **Hurdle model** at q < 0.05"
4. Provides ready-to-run command
5. Offers to execute
