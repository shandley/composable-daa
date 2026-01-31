# Effect Size Interpretation Guide

## Understanding DAA Effect Sizes

Effect sizes in differential abundance analysis depend on the method and transformation used.

## Method-Specific Interpretation

### LinDA (CLR-Transformed Linear Model)

LinDA uses Centered Log-Ratio transformation, which **attenuates** effect sizes.

#### The Attenuation Problem

| True log2FC | Observed CLR Effect | Attenuation |
|-------------|---------------------|-------------|
| 1.0 (2x) | ~0.25 | 75% |
| 2.0 (4x) | ~0.50 | 75% |
| 3.0 (8x) | ~0.75 | 75% |
| 4.0 (16x) | ~1.00 | 75% |

#### Why This Happens

CLR centers each sample around its geometric mean. When one feature increases, the geometric mean shifts, reducing the apparent effect on all features.

#### Correcting for Attenuation

To estimate true fold change from LinDA output:

```
Approximate true log2FC = CLR_estimate × 4
True fold change = 2^(CLR_estimate × 4)
```

| CLR Estimate | Approx log2FC | Approx Fold Change |
|--------------|---------------|-------------------|
| 0.25 | 1.0 | 2x |
| 0.50 | 2.0 | 4x |
| 0.75 | 3.0 | 8x |
| 1.00 | 4.0 | 16x |
| 1.25 | 5.0 | 32x |
| 1.50 | 6.0 | 64x |

**Important**: This is an approximation. Actual attenuation varies by:
- Number of features
- Dominance structure
- Proportion of truly differential features

### ZINB and Hurdle Models

These models work on raw counts with log link. Estimates are in natural log scale.

#### Direct Interpretation

```
Fold change = exp(estimate)
```

| Estimate | Fold Change | Direction |
|----------|-------------|-----------|
| -2.30 | 0.10x (10x down) | Decrease |
| -1.39 | 0.25x (4x down) | Decrease |
| -0.69 | 0.50x (2x down) | Decrease |
| 0.00 | 1x (no change) | None |
| 0.69 | 2x | Increase |
| 1.10 | 3x | Increase |
| 1.39 | 4x | Increase |
| 1.61 | 5x | Increase |
| 2.30 | 10x | Increase |
| 3.00 | 20x | Increase |
| 4.61 | 100x | Increase |

#### Hurdle Two-Part Interpretation

Hurdle models have two components:

1. **Binary component** (zero vs non-zero)
   - Estimate in log-odds scale
   - exp(estimate) = odds ratio
   - Positive = more likely to be present

2. **Count component** (abundance when present)
   - Estimate in log scale
   - exp(estimate) = fold change
   - Positive = higher abundance

### Negative Binomial GLM

Same interpretation as ZINB count component:

```
Fold change = exp(estimate)
```

## Biological Significance vs Statistical Significance

### What's a Meaningful Effect?

| Fold Change | Biological Relevance |
|-------------|---------------------|
| <1.5x | Usually noise, rarely meaningful |
| 1.5-2x | Small effect, may be meaningful in some contexts |
| 2-4x | Moderate effect, often biologically relevant |
| 4-10x | Large effect, strong biological signal |
| >10x | Very large effect, check for data issues |

### Context Matters

- **Dominant taxa**: Even 2x changes can have large ecosystem impact
- **Rare taxa**: May need >10x changes to matter functionally
- **Pathogens**: Any significant change may be clinically relevant

## Confidence Categories

The `confidence` column in results indicates:

| Confidence | Meaning |
|------------|---------|
| high | Large effect, low q-value, good prevalence |
| moderate | Medium effect or higher q-value |
| suggestive | Borderline significance, worth following up |
| not_significant | No evidence of differential abundance |

### How Confidence is Determined

```
IF q_value < 0.01 AND |effect| > 2x AND prevalence > 50%:
    confidence = "high"
ELIF q_value < 0.05 AND |effect| > 1.5x:
    confidence = "moderate"
ELIF q_value < 0.10 OR (q_value < 0.15 AND |effect| > 2x):
    confidence = "suggestive"
ELSE:
    confidence = "not_significant"
```

## Comparing Effect Sizes Across Methods

### Same Data, Different Methods

| Method | Reported Effect | True Effect Estimate |
|--------|-----------------|---------------------|
| LinDA | 0.5 CLR units | ~4x fold change |
| ZINB | 1.39 log units | 4x fold change |
| Hurdle | 1.39 log units | 4x fold change |
| NB | 1.39 log units | 4x fold change |

### Why Methods Report Different Values

1. **LinDA**: CLR attenuation reduces apparent effect
2. **ZINB**: Models zeros separately, count effect is cleaner
3. **Hurdle**: Separates presence from abundance
4. **NB**: Direct count model, no zero handling

## Red Flags in Effect Sizes

### Suspiciously Large Effects

| Warning | Likely Cause |
|---------|--------------|
| >100x fold change | Data error, outliers, or very rare taxa |
| All effects same direction | Compositional artifact |
| Only rare taxa significant | Noise or real biology (check carefully) |
| Only dominant taxa significant | Compositional artifact |

### What to Do

1. **Check raw data** for outliers or errors
2. **Run `daa stress`** to quantify compositional effects
3. **Validate with second method** for concordance
4. **Check biological plausibility** with domain knowledge

## Summary Table

| Method | Effect Scale | Conversion | Attenuation |
|--------|--------------|------------|-------------|
| LinDA | CLR | ×4 to get log2FC | ~75% |
| ZINB | log | exp() for FC | None |
| Hurdle | log | exp() for FC | None |
| NB | log | exp() for FC | None |
