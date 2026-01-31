# DAA Method Selection Decision Rules

Evidence-based rules derived from comprehensive benchmarking.

## Master Decision Tree

```
START
  │
  ├─► Is data longitudinal/repeated measures?
  │     YES ──► LMM (Linear Mixed Model)
  │             • Formula: ~ group + time + (1 | subject)
  │             • Threshold: p < 0.05
  │             • Note: Same CLR attenuation as LinDA
  │     NO ───► Continue...
  │
  ├─► What is the sparsity level?
  │     │
  │     ├─► >70% (very sparse, typical virome)
  │     │     └─► HURDLE MODEL
  │     │         • Separates presence/absence from abundance
  │     │         • Best for structural zeros
  │     │         • Threshold: q < 0.05
  │     │         • Expected: 83% sensitivity, 25% FDR at 4x FC
  │     │
  │     ├─► 50-70% (moderately sparse, typical 16S)
  │     │     └─► ZINB
  │     │         • Models zero-inflation as mixture
  │     │         • Good for excess sampling zeros
  │     │         • Threshold: q < 0.05
  │     │         • Expected: 83% sensitivity, 29% FDR at 4x FC
  │     │
  │     ├─► 30-50% (moderate sparsity)
  │     │     └─► ZINB or NB
  │     │         • ZINB if zeros seem inflated
  │     │         • NB if zeros are sampling-expected
  │     │         • Threshold: q < 0.05
  │     │
  │     └─► <30% (low sparsity)
  │           └─► LinDA or NB
  │               • LinDA: q < 0.10 (handles compositionality)
  │               • NB: q < 0.05 (simpler model)
  │
  ├─► What is the primary goal?
  │     │
  │     ├─► DISCOVERY (maximize true positives)
  │     │     └─► Use ZINB or Hurdle at q < 0.05
  │     │         • Accept ~25-30% FDR
  │     │         • Validate top hits independently
  │     │
  │     ├─► CONFIRMATION (minimize false positives)
  │     │     └─► Use LinDA at q < 0.10
  │     │         • Accept lower sensitivity (~39%)
  │     │         • High confidence in findings (~12% FDR)
  │     │
  │     └─► UNSURE
  │           └─► Run multiple methods, report concordance
  │               • High confidence: significant in LinDA + ZINB
  │               • Medium confidence: significant in ZINB only
  │
  └─► Sample size check
        │
        ├─► n < 10 per group
        │     └─► WARNING: Very underpowered
        │         • Only massive effects (>16x) detectable
        │         • Consider: pooling, different design
        │         • If must proceed: Permutation test
        │
        ├─► n = 10-20 per group
        │     └─► WARNING: Limited power
        │         • Need >4x fold changes
        │         • Use ZINB/Hurdle for best sensitivity
        │         • Don't expect to find subtle effects
        │
        ├─► n = 20-50 per group
        │     └─► Moderate power
        │         • Can detect >2x fold changes
        │         • Full method suite appropriate
        │
        └─► n > 50 per group
              └─► Good power
                  • Can detect moderate effects
                  • LinDA becomes more viable
```

## Sparsity-Based Selection

| Sparsity | Recommended | Alternative | Avoid |
|----------|-------------|-------------|-------|
| >80% | Hurdle | Permutation | LinDA, NB |
| 70-80% | Hurdle | ZINB | NB |
| 50-70% | ZINB | Hurdle | NB |
| 30-50% | ZINB | NB, LinDA | - |
| <30% | LinDA, NB | ZINB | - |

## Sample Size Guidelines

| n per group | Power Level | Detectable Effect | Recommended Method |
|-------------|-------------|-------------------|-------------------|
| <10 | Very low | >16x only | Permutation (or reconsider study) |
| 10-20 | Low | >4x | ZINB, Hurdle |
| 20-30 | Moderate | >2-4x | Any method |
| 30-50 | Good | >2x | Any method |
| >50 | High | >1.5x | Any method, LinDA viable |

## Method Performance Summary

### At 4.0 log2FC (16x fold change), n=20/group

| Method | Threshold | Sensitivity | FDR | Best For |
|--------|-----------|-------------|-----|----------|
| Hurdle | q < 0.05 | 83% | 25% | Very sparse data |
| ZINB | q < 0.05 | 83% | 29% | Moderately sparse |
| LinDA | q < 0.10 | 39% | 12.5% | FDR control |
| NB | q < 0.05 | 6% | 0% | Low sparsity |
| Permutation | p < 0.05 | ~30% | ~5% | Unknown distributions |

### At 2.0 log2FC (4x fold change), n=20/group

| Method | Threshold | Sensitivity | FDR |
|--------|-----------|-------------|-----|
| Hurdle | q < 0.05 | 26% | 17% |
| ZINB | q < 0.05 | 58% | 45% |
| LinDA | q < 0.10 | 0% | n/a |
| NB | q < 0.05 | 0% | n/a |

## Warning Triggers

### High Priority Warnings

| Condition | Warning | Action |
|-----------|---------|--------|
| n < 10/group | "Very underpowered study" | Suggest pooling or redesign |
| Sparsity > 90% | "Extreme sparsity" | Use Hurdle, consider presence-only |
| Library size CV > 2.0 | "High technical variation" | Check for batch effects |
| Groups have very different n | "Unbalanced design" | Note reduced power |

### Medium Priority Warnings

| Condition | Warning | Action |
|-----------|---------|--------|
| Library size differs by group | "Potential confounding" | Add as covariate or check |
| Sparsity differs by group | "Group-specific sparsity" | Use groupwise prevalence filter |
| <50 features pass filter | "Few testable features" | Consider relaxing filter |

## Formula Construction

### Basic Two-Group Comparison
```
Formula: ~ group
Test coefficient: group{treatment_level}
Example: daa zinb -f "~ treatment" -t treatmentdisease
```

### With Covariates
```
Formula: ~ group + age + sex
Test coefficient: group{treatment_level}
Example: daa zinb -f "~ treatment + age + sex" -t treatmentdisease
```

### Longitudinal Data
```
Formula: ~ group + time + (1 | subject)
Test coefficient: group{treatment_level}
Example: daa lmm -f "~ treatment + timepoint + (1 | subject_id)" -t treatmentdisease
```

## Prevalence Filter Guidelines

| Data Type | Recommended Threshold | Rationale |
|-----------|----------------------|-----------|
| 16S (typical) | 10% | Balance detection vs noise |
| Virome (sparse) | 5% | Don't lose rare viruses |
| Shotgun metagenomics | 10-20% | Higher depth allows stricter |
| Low sample size (<30) | 20% | Need more observations |

## Output Interpretation Guide

### LinDA Results
- Effect sizes are CLR-transformed (NOT fold changes)
- Multiply observed effect by ~4 for approximate true fold change
- Use q < 0.10 threshold
- Focus on significance, not magnitude

### ZINB/Hurdle Results
- Effect sizes are log-scale (natural log)
- exp(estimate) = fold change
- estimate of 1.0 ≈ 2.7x fold change
- estimate of 2.3 ≈ 10x fold change

### Concordance Analysis
When running multiple methods:
- **High confidence**: Significant in LinDA (q<0.10) AND ZINB/Hurdle (q<0.05)
- **Medium confidence**: Significant in ZINB/Hurdle only
- **Low confidence**: Significant in only one sensitive method
- **Likely false positive**: Only significant in most sensitive method at lenient threshold
