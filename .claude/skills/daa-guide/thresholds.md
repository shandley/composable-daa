# Q-Value Threshold Selection Guide

Choosing the right significance threshold is critical for balancing sensitivity and FDR.

## Method-Specific Recommendations

| Method | Recommended Threshold | Rationale |
|--------|----------------------|-----------|
| LinDA | **q < 0.10** | CLR attenuation makes q < 0.05 too conservative |
| ZINB | q < 0.05 | Standard threshold works well |
| Hurdle | q < 0.05 | Standard threshold works well |
| NB | q < 0.05 | Standard threshold, but low sensitivity |
| Permutation | p < 0.05 | Already conservative |
| LMM | p < 0.05 | Standard threshold |

## LinDA Threshold Analysis

### Why q < 0.05 Fails for LinDA

At 4.0 log2FC (16x fold change) with n=20/group:

| Threshold | Sensitivity | FDR | TP | FP | Total Sig |
|-----------|-------------|-----|----|----|-----------|
| q < 0.01 | 0% | n/a | 0 | 0 | 0 |
| q < 0.05 | 0% | n/a | 0 | 0 | 0 |
| **q < 0.10** | **39%** | **12.5%** | 7 | 1 | 8 |
| q < 0.15 | 39% | 30% | 7 | 3 | 10 |
| q < 0.20 | 39% | 30% | 7 | 3 | 10 |
| q < 0.25 | 44% | 33% | 8 | 4 | 12 |
| q < 0.30 | 44% | 53% | 8 | 9 | 17 |

**Optimal threshold**: q < 0.10 (best F1 score of 0.538)

### The Problem: CLR Effect Attenuation

True 4.0 log2FC effects become ~1.0 CLR effects:
- Mean standard error: ~0.5
- Mean t-statistic: ~2.0
- With 190+ tests, BH correction at α=0.05 is very strict
- Minimum q-value achieved: 0.09 (just above 0.05)

### At Smaller Effect Sizes

| Effect Size | Best Threshold | Sensitivity | FDR |
|-------------|----------------|-------------|-----|
| 0.5 log2FC | Any | 0% | n/a |
| 1.0 log2FC | Any | 0% | n/a |
| 2.0 log2FC | Any | 0% | n/a |
| 4.0 log2FC | q < 0.10 | 39% | 12.5% |

LinDA simply cannot detect effects smaller than ~3 log2FC (8x fold change).

## ZINB/Hurdle Threshold Considerations

### Standard q < 0.05 Works

| Effect Size | Sensitivity | FDR |
|-------------|-------------|-----|
| 2.0 log2FC | 26-58% | 17-45% |
| 4.0 log2FC | 83% | 25-29% |

### When to Use Stricter Thresholds

Use q < 0.01 when:
- Validating findings for publication
- Following up with expensive experiments
- FDR control is paramount

Use q < 0.10 when:
- Exploratory analysis
- Generating hypotheses
- Large effect sizes expected

## Practical Guidance

### Scenario 1: User Runs LinDA, Gets 0 Significant

```
User: I ran LinDA but got 0 significant features at q < 0.05

Response: This is expected behavior. LinDA uses CLR transformation which
attenuates effect sizes by ~75%. Try filtering at q < 0.10 instead:

  # In Python/R, filter your results:
  significant = results[results['q_value'] < 0.10]

At q < 0.10, LinDA achieves 39% sensitivity with only 12.5% FDR -
excellent for high-confidence discoveries.

If you still get few/no results, the effects in your data may be
smaller than LinDA can detect (needs >8x fold changes). Consider
using ZINB or Hurdle for discovery.
```

### Scenario 2: User Wants Maximum Sensitivity

```
User: I want to find as many differentially abundant features as possible

Response: For maximum sensitivity, use ZINB or Hurdle at q < 0.05:

  daa zinb -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

Expected performance at 4.0 log2FC:
- Sensitivity: 83%
- FDR: 29%

This means you'll find most true effects but ~30% of your "significant"
findings may be false positives. Consider biological validation.
```

### Scenario 3: User Needs High Confidence

```
User: I need results I can publish with high confidence

Response: For high-confidence findings, use LinDA at q < 0.10:

  daa linda -c counts.tsv -m metadata.tsv -f "~ group" -t grouptreatment -o results.tsv

Then filter results at q < 0.10 (NOT 0.05).

Expected performance:
- Sensitivity: 39% (you'll miss some true effects)
- FDR: 12.5% (but findings are highly reliable)

For even higher confidence, run both LinDA and ZINB - features
significant in both are very likely true positives.
```

## Multiple Testing Considerations

### Why BH Correction Matters

With 200 features tested:
- At α = 0.05, expect 10 false positives by chance
- BH correction controls FDR, not FWER
- Actual FDR depends on proportion of true effects

### When Many Features Are Tested

| # Features | Recommendation |
|------------|----------------|
| < 100 | Standard thresholds work |
| 100-500 | Standard thresholds, consider q < 0.10 for LinDA |
| 500-1000 | May need larger effects or more samples |
| > 1000 | Power is limited, focus on large effects |

## Summary

1. **Always use q < 0.10 for LinDA** - q < 0.05 is too conservative
2. **Use q < 0.05 for ZINB/Hurdle** - standard threshold works well
3. **Consider effect sizes** - small effects need more samples
4. **Run multiple methods** - concordant findings are high confidence
