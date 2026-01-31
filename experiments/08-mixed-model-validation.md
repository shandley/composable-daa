# Experiment 08: Mixed Model Validation

## Summary

Validation of the linear mixed model (LMM) implementation for longitudinal and repeated measures microbiome data, comparing to standard linear models and assessing random effect specification.

## Background

The toolkit implements LMM with:
- REML estimation for variance components
- Random intercepts: `(1 | subject)`
- Random slopes: `(1 + time | subject)`
- Satterthwaite and Kenward-Roger degrees of freedom

### Why Mixed Models Matter

Microbiome studies often have:
- **Longitudinal designs**: Samples at multiple time points
- **Repeated measures**: Multiple samples per subject
- **Clustering**: Subjects nested within sites/batches

Ignoring correlation structure leads to:
- Inflated Type I error (pseudoreplication)
- Loss of power (inefficient estimation)
- Biased effect estimates

## Research Questions

1. **Does LMM properly handle correlation?**
   - Is FPR controlled with correlated data?
   - Does ignoring correlation inflate FPR?
   - How much power is lost with naive LM?

2. **How do different random effect structures perform?**
   - Random intercepts only vs intercepts + slopes
   - Correct vs misspecified random effects
   - Convergence across data structures

3. **Are degrees of freedom methods appropriate?**
   - Satterthwaite vs Kenward-Roger vs residual
   - Which is most conservative?
   - Impact on inference

## Proposed Analysis

### Phase 1: Simulation Study

1. Generate longitudinal data:
   - n = 10, 20, 50 subjects
   - t = 3, 5, 10 time points per subject
   - ICC = 0.1, 0.3, 0.5, 0.7 (within-subject correlation)
   - True effect on time or treatment

2. Compare methods:
   - LMM with correct specification
   - LMM with wrong specification
   - Naive LM (ignoring correlation)
   - GEE (generalized estimating equations, if implemented)

3. Evaluate:
   - FPR under null (should be ~5%)
   - Power for detecting true effects
   - Effect estimate accuracy

### Phase 2: Random Effect Specification

1. Test sensitivity to misspecification:
   - True random intercepts, fit random slopes
   - True random slopes, fit random intercepts only
   - True crossed effects, fit nested

2. Assess:
   - Bias in fixed effect estimates
   - SE accuracy
   - Convergence rates

### Phase 3: Degrees of Freedom

1. Compare DF methods across scenarios:
   - Balanced vs unbalanced designs
   - Small vs large clusters
   - Few vs many clusters

2. Evaluate:
   - Coverage of confidence intervals
   - FPR calibration
   - Conservatism

## Data Requirements

- Synthetic longitudinal data with controlled ICC
- Varying cluster sizes and numbers
- Both balanced and unbalanced designs

## Expected Outcomes

1. Guidelines for random effect specification in microbiome LMM
2. DF method recommendations
3. Power analysis for longitudinal designs
4. Methods paper section on mixed models for microbiome data

## Relevance to Microbiome Studies

Common designs that require LMM:
- **Diet interventions**: Pre/post samples from same subjects
- **Disease progression**: Samples over time
- **Multi-site studies**: Subjects nested in sites
- **Technical replicates**: Multiple extractions per sample

## Code Location

Scripts: `experiments/scripts/08-mixed-model-validation/`

## References

- Bates et al. (2015) - lme4: Linear mixed models in R
- Barr et al. (2013) - Random effects structure in psycholinguistics
- Bolker et al. (2009) - GLMMs for ecology
- La Rosa et al. (2012) - Mixed models for microbiome studies
