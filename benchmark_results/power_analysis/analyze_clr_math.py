#!/usr/bin/env python3
"""Analyze the mathematical reason for CLR effect size attenuation."""

import pandas as pd
import numpy as np

def analyze_clr_math():
    """Explain and demonstrate CLR effect size attenuation."""

    print("=" * 80)
    print("CLR EFFECT SIZE ATTENUATION: MATHEMATICAL ANALYSIS")
    print("=" * 80)

    # Load the count data
    counts = pd.read_csv("effect_4.0/counts.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv("effect_4.0/metadata.tsv", sep="\t")
    gt = pd.read_csv("effect_4.0/ground_truth.tsv", sep="\t")

    # Identify control and treatment samples
    control_samples = metadata[metadata["group"] == "control"]["sample_id"].tolist()
    treatment_samples = metadata[metadata["group"] == "treatment"]["sample_id"].tolist()

    # Get spiked features
    up_spiked = gt[gt["direction"] == "up"]["feature_id"].tolist()
    down_spiked = gt[gt["direction"] == "down"]["feature_id"].tolist()

    print("\n1. RAW COUNT DIFFERENCES FOR SPIKED FEATURES")
    print("-" * 50)

    # Calculate raw fold changes
    for feature in up_spiked[:3]:  # Show first 3
        if feature in counts.index:
            control_mean = counts.loc[feature, control_samples].mean()
            treatment_mean = counts.loc[feature, treatment_samples].mean()
            if control_mean > 0:
                raw_fc = treatment_mean / control_mean
                raw_log2fc = np.log2(raw_fc) if raw_fc > 0 else np.nan
                print(f"  {feature}: control={control_mean:.1f}, treatment={treatment_mean:.1f}")
                print(f"    Raw FC={raw_fc:.2f}, log2FC={raw_log2fc:.2f}")

    print("\n2. WHY CLR ATTENUATES EFFECTS")
    print("-" * 50)
    print("""
CLR transformation: CLR(x_i) = log(x_i) - mean(log(x))

When feature i increases by factor k:
- x_i' = k * x_i
- log(x_i') = log(k) + log(x_i)
- Geometric mean also increases (because x_i contributes)
- CLR(x_i') = log(k) + log(x_i) - [mean(log(x)) + log(k)/D]
- Net effect: log(k) * (1 - 1/D) where D = number of features

For D=200 features:
- Theoretical attenuation factor = (1 - 1/200) = 0.995 (almost no attenuation)

BUT in practice, the attenuation is much worse because:
1. Many features are correlated (not independent)
2. Compositional closure: when some features increase, others must decrease
3. The spike-in data uses compositional mode (redistributes reads)
""")

    # Actually compute CLR and compare
    print("\n3. COMPUTING CLR EFFECT SIZES")
    print("-" * 50)

    # Add pseudocount and compute CLR
    pseudocount = 0.5
    counts_pc = counts + pseudocount

    # CLR: log(x) - mean(log(x)) for each sample
    log_counts = np.log(counts_pc)
    clr = log_counts.sub(log_counts.mean(axis=0), axis=1)

    # Compare CLR effect sizes
    print("\nUP-spiked features (true log2FC = +4.0):")
    up_effects = []
    for feature in up_spiked:
        if feature in clr.index:
            control_mean = clr.loc[feature, control_samples].mean()
            treatment_mean = clr.loc[feature, treatment_samples].mean()
            clr_diff = treatment_mean - control_mean
            up_effects.append(clr_diff)
            print(f"  {feature}: CLR diff = {clr_diff:.3f}")

    print(f"\n  Mean CLR effect: {np.mean(up_effects):.3f}")
    print(f"  Attenuation: {np.mean(up_effects)/4.0:.1%}")

    print("\nDOWN-spiked features (true log2FC = -4.0):")
    down_effects = []
    for feature in down_spiked:
        if feature in clr.index:
            control_mean = clr.loc[feature, control_samples].mean()
            treatment_mean = clr.loc[feature, treatment_samples].mean()
            clr_diff = treatment_mean - control_mean
            down_effects.append(clr_diff)
            print(f"  {feature}: CLR diff = {clr_diff:.3f}")

    print(f"\n  Mean CLR effect: {np.mean(down_effects):.3f}")
    print(f"  Attenuation: {np.mean(down_effects)/-4.0:.1%}")

    # Key insight
    print("\n4. THE KEY INSIGHT: COMPOSITIONAL NORMALIZATION")
    print("-" * 50)
    print("""
The spike-in uses COMPOSITIONAL mode, which:
1. Increases spiked features by 16x
2. Decreases all OTHER features proportionally (to maintain library size)

This creates a "double penalty" for CLR:
1. The spiked feature goes UP
2. But ALL other features go DOWN (slightly)
3. The geometric mean shifts toward the spiked features
4. This partially cancels out the spike effect

In absolute abundance studies (with spike-in controls), you'd see
the full 4.0 log2FC. But in compositional data, the relative change
is dampened by the constraints.
""")

    # Compare with other normalization
    print("\n5. COMPARISON: CLR vs TSS vs RAW")
    print("-" * 50)

    # TSS (proportions)
    tss = counts.div(counts.sum(axis=0), axis=1)

    print("\nFor first 3 UP-spiked features:")
    for feature in up_spiked[:3]:
        if feature in counts.index:
            # Raw counts
            raw_control = counts.loc[feature, control_samples].mean()
            raw_treatment = counts.loc[feature, treatment_samples].mean()
            raw_log2fc = np.log2(raw_treatment / raw_control) if raw_control > 0 else np.nan

            # TSS (proportions)
            tss_control = tss.loc[feature, control_samples].mean()
            tss_treatment = tss.loc[feature, treatment_samples].mean()
            tss_log2fc = np.log2(tss_treatment / tss_control) if tss_control > 0 else np.nan

            # CLR
            clr_control = clr.loc[feature, control_samples].mean()
            clr_treatment = clr.loc[feature, treatment_samples].mean()
            clr_diff = clr_treatment - clr_control

            print(f"\n  {feature}:")
            print(f"    Raw counts log2FC:  {raw_log2fc:.3f}")
            print(f"    TSS log2FC:         {tss_log2fc:.3f}")
            print(f"    CLR diff:           {clr_diff:.3f}")

    print("\n6. RECOMMENDATIONS")
    print("-" * 50)
    print("""
For LinDA users:
1. Use q < 0.10 instead of q < 0.05 for discovery
2. Expect to detect only very large effects (>8x fold change)
3. For absolute abundance studies, use spike-in normalization
4. LinDA excels at FDR control - use it for confident discoveries

For method comparison:
- LinDA is CONSERVATIVE (low sensitivity, excellent FDR control)
- ZINB is SENSITIVE (high sensitivity, moderate FDR)
- Choose based on research goals:
  * Discovery (minimize false negatives): use ZINB
  * Confirmation (minimize false positives): use LinDA
""")

if __name__ == "__main__":
    analyze_clr_math()
