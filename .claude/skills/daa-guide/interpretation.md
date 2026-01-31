# Interpreting DAA Results

Guide for understanding and communicating differential abundance analysis results.

## Understanding Effect Sizes

### LinDA (CLR-Transformed)

**Effect sizes are NOT directly interpretable as fold changes.**

| Observed CLR Effect | Approximate True Fold Change |
|--------------------|------------------------------|
| 0.5 | ~2x (but uncertain) |
| 1.0 | ~4x (but uncertain) |
| 1.5 | ~6-8x (but uncertain) |
| 2.0 | ~8-16x (but uncertain) |

**Why?** CLR transformation attenuates effects by ~75%. An observed effect of 1.0
may represent a true 4x fold change, but the exact conversion depends on:
- Number of features
- Correlation structure
- Magnitude of changes in other features

**Recommendation**: Focus on significance (q-value), not effect magnitude.

### ZINB/Hurdle/NB (Log-Scale)

Effect sizes are on natural log scale:

| Estimate | Fold Change | Interpretation |
|----------|-------------|----------------|
| 0.69 | 2x | Doubled |
| 1.10 | 3x | Tripled |
| 1.39 | 4x | Quadrupled |
| 2.30 | 10x | 10-fold increase |
| -0.69 | 0.5x | Halved |
| -2.30 | 0.1x | 10-fold decrease |

**Conversion**: Fold change = exp(estimate)

### Hurdle Model Components

Hurdle models have two components:

1. **Binary component**: Log-odds of presence
   - Positive = more likely present in treatment
   - Negative = more likely absent in treatment

2. **Count component**: Log fold change in abundance (given presence)
   - Interpretation same as NB/ZINB

**Example interpretation**:
- Binary estimate = 1.5 → 4.5x more likely to be present
- Count estimate = 1.0 → 2.7x higher abundance when present

## Interpreting P-Values and Q-Values

### P-Value vs Q-Value

| Metric | What It Means | When to Use |
|--------|---------------|-------------|
| p-value | Probability of result under null hypothesis | Single comparisons |
| q-value | Expected FDR at this threshold | Multiple testing |

### Q-Value Interpretation

| Q-Value | Interpretation |
|---------|----------------|
| q < 0.01 | Very strong evidence, <1% expected false positives |
| q < 0.05 | Strong evidence, <5% expected false positives |
| q < 0.10 | Moderate evidence, <10% expected false positives |
| q < 0.20 | Suggestive, <20% expected false positives |
| q > 0.20 | Weak evidence, high false positive risk |

### Method-Specific Q-Value Guidance

| Method | Recommended | Why |
|--------|-------------|-----|
| LinDA | q < 0.10 | CLR attenuation makes 0.05 too strict |
| ZINB | q < 0.05 | Standard threshold works |
| Hurdle | q < 0.05 | Standard threshold works |
| NB | q < 0.05 | Standard threshold works |

## Confidence Levels

Results include a `confidence` column based on prevalence:

| Confidence | Meaning | Action |
|------------|---------|--------|
| high | High prevalence, reliable estimate | Trust the result |
| medium | Moderate prevalence, reasonable estimate | Interpret with care |
| low | Low prevalence, uncertain estimate | Validate independently |
| suggestive | Below threshold but interesting | Hypothesis only |
| not_significant | Above threshold | No evidence of effect |

## Common Interpretation Mistakes

### Mistake 1: LinDA Effect = Fold Change

**Wrong**: "LinDA estimate of 1.5 means 1.5-fold increase"

**Right**: "LinDA estimate of 1.5 is a CLR-transformed value. The true fold change
is likely much larger (~6x), but CLR effects should be interpreted for significance,
not magnitude."

### Mistake 2: High Q-Value = No Effect

**Wrong**: "q = 0.15 means there's no difference"

**Right**: "q = 0.15 means we don't have strong statistical evidence at conventional
thresholds. The feature may still be biologically different, but our study lacks
power to detect it confidently."

### Mistake 3: Ignoring Multiple Testing

**Wrong**: "p < 0.05 so it's significant"

**Right**: "With 200 features tested, we expect ~10 false positives at p < 0.05.
Use q-values (FDR-adjusted) for interpretation."

### Mistake 4: Comparing Effect Sizes Across Methods

**Wrong**: "LinDA effect of 1.0 < ZINB effect of 2.0, so ZINB found a bigger effect"

**Right**: "LinDA and ZINB use different scales. CLR effects are attenuated;
log-scale effects are not. They cannot be directly compared."

## Reporting Results

### For Publications

```
Differential abundance analysis was performed using [METHOD] with
Benjamini-Hochberg correction for multiple testing. Features with
q < [THRESHOLD] were considered significant.

[For LinDA]: Note that effect sizes are CLR-transformed and represent
relative changes after compositional normalization, not direct fold
changes.

[For ZINB/Hurdle/NB]: Effect sizes are reported on the natural log
scale; fold changes can be computed as exp(estimate).
```

### For Collaborators

```
Key findings:
- X features significantly different (q < [THRESHOLD])
- Y features increased in treatment
- Z features decreased in treatment

Top hits:
1. Feature_A: [fold change]x [higher/lower], q = [value]
2. Feature_B: [fold change]x [higher/lower], q = [value]

Caveats:
- [For LinDA]: Only very large effects (>8x) detectable; smaller
  differences may exist but weren't detected
- FDR of ~[X]% expected; some findings may be false positives
```

## When Results Are Unexpected

### No Significant Features

Possible reasons:
1. **Threshold too strict** → Try q < 0.10 (especially for LinDA)
2. **Effects too small** → Need more samples or larger biological effect
3. **High variability** → Consider covariates or batch effects
4. **No true effect** → Groups may not differ

### Too Many Significant Features

Possible reasons:
1. **Batch effects** → Check for confounding
2. **Compositional artifacts** → Consider absolute abundance methods
3. **Threshold too lenient** → Use stricter threshold
4. **True biological signal** → Large perturbation has broad effects

### Different Methods Disagree

| Scenario | Interpretation |
|----------|----------------|
| LinDA sig, ZINB not sig | Unusual; check data carefully |
| LinDA not sig, ZINB sig | Common; LinDA more conservative |
| Both sig, different direction | Error; check for bugs |
| Both sig, same direction | High confidence finding |

## Sample Output Interpretation

### Example LinDA Output

```
feature_id      estimate    std_error    statistic    p_value    q_value
Feature_A       1.94        0.56         3.47         0.001      0.091
Feature_B       -2.02       0.59         -3.43        0.001      0.091
Feature_C       0.69        0.82         0.84         0.408      0.894
```

**Interpretation**:
- Feature_A: Significantly higher in treatment (q=0.091 < 0.10)
  - CLR effect of 1.94 suggests substantial increase (likely >4x)
- Feature_B: Significantly lower in treatment (q=0.091 < 0.10)
  - CLR effect of -2.02 suggests substantial decrease
- Feature_C: Not significant (q=0.894 > 0.10)
  - No evidence of differential abundance

### Example ZINB Output

```
feature_id      estimate    std_error    statistic    p_value    q_value
Feature_A       2.30        0.45         5.11         0.000      0.001
Feature_B       -1.61       0.38         -4.24        0.000      0.003
Feature_C       0.35        0.52         0.67         0.502      0.721
```

**Interpretation**:
- Feature_A: 10-fold higher (exp(2.30) = 10x), highly significant
- Feature_B: 5-fold lower (exp(-1.61) = 0.2x), highly significant
- Feature_C: Not significant, ~1.4x change but uncertain
