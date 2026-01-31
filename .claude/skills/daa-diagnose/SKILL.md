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

### Study Design Detection

Inspect metadata for study design indicators:

```bash
head -5 {METADATA_FILE}
```

Look for these columns:

| Column Pattern | Study Design | Recommended Model |
|----------------|--------------|-------------------|
| `subject`, `patient`, `individual`, `id` | Repeated measures | LMM with random intercept |
| `time`, `timepoint`, `visit`, `day` | Longitudinal | LMM with time + random intercept |
| `batch`, `run`, `plate`, `lane` | Batch effects | Include as covariate or random effect |
| `age`, `bmi`, `sex` | Covariates | Include in formula |

### Longitudinal Data Rules

```
IF has_subject_column AND has_time_column:
    → Study type: Longitudinal
    → Recommend: LMM with random intercept per subject
    → Formula: "~ group * time + (1 | subject)"

IF has_subject_column AND NOT has_time_column:
    → Study type: Repeated measures (paired)
    → Recommend: LMM with random intercept per subject
    → Formula: "~ group + (1 | subject)"

IF multiple_samples_per_subject:
    → MUST use LMM to account for non-independence
    → Standard models will have inflated Type I error
```

### Covariate Recommendations

```
IF has_batch_column:
    → Always include batch in model
    → Consider as random effect if many batches: (1 | batch)
    → Or fixed effect if few batches: ~ group + batch

IF has_continuous_covariates (age, bmi, etc.):
    → Include if potentially confounding
    → Formula: "~ group + age + bmi"

IF library_size_imbalance > 2x AND correlated_with_group:
    → Include log(library_size) as covariate
    → Formula: "~ group + log_library_size"
```

## Step 5: Generate Recommendations

Output a structured recommendation:

### Template Output

```
## Data Diagnosis Summary

**Dataset**: {n_features} features × {n_samples} samples
**Sparsity**: {sparsity}%
**Study design**: {cross-sectional | longitudinal | repeated-measures}
**Samples per group**: {group_sizes}
**Library size CV**: {cv}
**Covariates detected**: {list of potential covariates}

## Recommended Analysis

### Primary Method: {METHOD}
Rationale: {why this method based on data characteristics}

### Formula: {formula}
{If longitudinal: explain random effects}
{If covariates: explain why included}

```bash
daa recommend -c {counts} -m {metadata} -g {group_column} -t {target_level} --run -o results.tsv
```

### Threshold: {q_threshold}
- For {METHOD}: use q < {threshold}
{If LinDA: explain q < 0.10 requirement}

### Secondary Validation (Optional)
Run spike-in validation to empirically test the method on your data:
```bash
daa validate -c {counts} -m {metadata} -g {group_column} -t {target_level} -f "{formula}" --test-coef {coefficient}
```

## Warnings
{Any warnings about power, confounding, non-independence, etc.}

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

## Example Sessions

### Example 1: Cross-sectional Study

**User**: "Analyze my microbiome data in data/counts.tsv"

**Claude**:
1. Runs `daa profile-llm -c data/counts.tsv -m data/metadata.tsv -g treatment`
2. Sees: 68% sparsity, n=25/group, CV=0.8
3. Recommends: "Based on 68% sparsity and n=25/group, I recommend **Hurdle model** at q < 0.05"
4. Provides ready-to-run command
5. Offers to execute

### Example 2: Longitudinal Study

**User**: "Analyze my longitudinal microbiome study"

**Claude**:
1. Inspects metadata: sees `subject_id`, `timepoint`, `treatment` columns
2. Detects: Longitudinal design with repeated measures
3. Runs `daa profile-llm -c counts.tsv -m metadata.tsv -g treatment`
4. Sees: 55% sparsity, n=30 subjects, 3 timepoints each
5. Recommends:
   ```
   Study design: Longitudinal (3 timepoints per subject)

   Primary Method: LMM (Linear Mixed Model)
   Rationale: Repeated measures require random effects to account for
   within-subject correlation. Moderate sparsity allows CLR transformation.

   Formula: "~ treatment * timepoint + (1 | subject_id)"
   - treatment: main effect of interest
   - timepoint: time effect
   - treatment * timepoint: interaction (different trajectories?)
   - (1 | subject_id): random intercept per subject

   Note: For longitudinal designs with random effects, use daa run with
   a YAML pipeline config, as recommend doesn't yet auto-detect longitudinal designs.
   ```
6. Warns: "Use LMM to avoid inflated Type I error from non-independent samples"

### Example 3: Study with Batch Effects

**User**: "My data has batch effects from different sequencing runs"

**Claude**:
1. Inspects metadata: sees `batch`, `group` columns
2. Detects: 4 batches, potential confounding
3. Recommends:
   ```
   Covariates detected: batch (4 levels)

   Formula: "~ group + batch"
   - Include batch as fixed effect (few levels)

   Alternative if many batches:
   Formula: "~ group + (1 | batch)"
   - Include batch as random effect
   ```
