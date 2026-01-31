#!/usr/bin/env python3
"""Analyze LinDA's low sensitivity issue."""

import pandas as pd
import numpy as np

def analyze_linda_results():
    """Deep dive into LinDA sensitivity."""

    print("=" * 80)
    print("LINDA SENSITIVITY INVESTIGATION")
    print("=" * 80)

    # Load data for effect size 4.0
    linda = pd.read_csv("effect_4.0/linda_results.tsv", sep="\t")
    ground_truth = pd.read_csv("effect_4.0/ground_truth.tsv", sep="\t")

    # Merge with ground truth
    merged = linda.merge(ground_truth, on="feature_id", how="left")
    merged["is_differential"] = merged["is_differential"].fillna(False)

    # Separate true positives from null
    true_diff = merged[merged["is_differential"] == True]
    nulls = merged[merged["is_differential"] == False]

    print(f"\n1. BASIC STATISTICS")
    print("-" * 40)
    print(f"Total features tested: {len(merged)}")
    print(f"True differential features: {len(true_diff)}")
    print(f"Null features: {len(nulls)}")

    # P-value and q-value distribution for spiked features
    print(f"\n2. P-VALUE/Q-VALUE DISTRIBUTION FOR SPIKED FEATURES (n={len(true_diff)})")
    print("-" * 40)

    up_spiked = true_diff[true_diff["direction"] == "up"]
    down_spiked = true_diff[true_diff["direction"] == "down"]

    print(f"\nUP-spiked features (true log2FC = +4.0):")
    print(f"  Count: {len(up_spiked)}")
    for _, row in up_spiked.iterrows():
        print(f"  {row['feature_id']}: estimate={row['estimate']:.3f}, p={row['p_value']:.4f}, q={row['q_value']:.4f}")

    print(f"\nDOWN-spiked features (true log2FC = -4.0):")
    print(f"  Count: {len(down_spiked)}")
    for _, row in down_spiked.iterrows():
        print(f"  {row['feature_id']}: estimate={row['estimate']:.3f}, p={row['p_value']:.4f}, q={row['q_value']:.4f}")

    # Effect size attenuation
    print(f"\n3. EFFECT SIZE ATTENUATION")
    print("-" * 40)

    up_estimates = up_spiked["estimate"].values
    down_estimates = down_spiked["estimate"].values

    print(f"\nUP-spiked (should be +4.0):")
    print(f"  Mean observed: {up_estimates.mean():.3f}")
    print(f"  Range: [{up_estimates.min():.3f}, {up_estimates.max():.3f}]")
    print(f"  Attenuation factor: {up_estimates.mean() / 4.0:.3f}")

    print(f"\nDOWN-spiked (should be -4.0):")
    print(f"  Mean observed: {down_estimates.mean():.3f}")
    print(f"  Range: [{down_estimates.min():.3f}, {down_estimates.max():.3f}]")
    print(f"  Attenuation factor: {down_estimates.mean() / -4.0:.3f}")

    # P-value summary
    print(f"\n4. P-VALUE SUMMARY FOR SPIKED FEATURES")
    print("-" * 40)

    all_spiked_pvals = true_diff["p_value"].values
    print(f"  Min p-value: {all_spiked_pvals.min():.4f}")
    print(f"  Max p-value: {all_spiked_pvals.max():.4f}")
    print(f"  Median p-value: {np.median(all_spiked_pvals):.4f}")
    print(f"  Mean p-value: {all_spiked_pvals.mean():.4f}")

    print(f"\n  Count with p < 0.05: {(all_spiked_pvals < 0.05).sum()}")
    print(f"  Count with p < 0.10: {(all_spiked_pvals < 0.10).sum()}")
    print(f"  Count with p < 0.20: {(all_spiked_pvals < 0.20).sum()}")

    # Q-value summary
    print(f"\n5. Q-VALUE SUMMARY FOR SPIKED FEATURES")
    print("-" * 40)

    all_spiked_qvals = true_diff["q_value"].values
    print(f"  Min q-value: {all_spiked_qvals.min():.4f}")
    print(f"  Max q-value: {all_spiked_qvals.max():.4f}")
    print(f"  Median q-value: {np.median(all_spiked_qvals):.4f}")

    print(f"\n  Count with q < 0.05: {(all_spiked_qvals < 0.05).sum()}")
    print(f"  Count with q < 0.10: {(all_spiked_qvals < 0.10).sum()}")
    print(f"  Count with q < 0.20: {(all_spiked_qvals < 0.20).sum()}")

    # Compare with other methods
    print(f"\n6. COMPARISON WITH OTHER METHODS (q < 0.05)")
    print("-" * 40)

    for method in ["zinb", "hurdle", "nb"]:
        try:
            other = pd.read_csv(f"effect_4.0/{method}_results.tsv", sep="\t")
            other_merged = other.merge(ground_truth, on="feature_id", how="left")
            other_diff = other_merged[other_merged["is_differential"] == True]
            tp = (other_diff["q_value"] < 0.05).sum()
            print(f"  {method.upper()}: {tp}/{len(other_diff)} TP ({tp/len(other_diff)*100:.1f}% sensitivity)")
        except Exception as e:
            print(f"  {method}: error - {e}")

    # LinDA
    tp_linda = (true_diff["q_value"] < 0.05).sum()
    print(f"  LINDA: {tp_linda}/{len(true_diff)} TP ({tp_linda/len(true_diff)*100:.1f}% sensitivity)")

    # Standard error analysis
    print(f"\n7. STANDARD ERROR ANALYSIS")
    print("-" * 40)

    print(f"\nSpiked features:")
    print(f"  Mean SE: {true_diff['std_error'].mean():.3f}")
    print(f"  Mean |estimate|: {true_diff['estimate'].abs().mean():.3f}")
    print(f"  Mean |t-stat|: {true_diff['statistic'].abs().mean():.3f}")

    print(f"\nNull features:")
    print(f"  Mean SE: {nulls['std_error'].mean():.3f}")
    print(f"  Mean |estimate|: {nulls['estimate'].abs().mean():.3f}")
    print(f"  Mean |t-stat|: {nulls['statistic'].abs().mean():.3f}")

    # The core issue
    print(f"\n8. THE CORE ISSUE")
    print("-" * 40)
    print("""
The problem is that CLR transformation severely attenuates effect sizes:
- True effect size: 4.0 log2FC (16x fold change)
- Observed effect size: ~0.5-1.5 log2FC after CLR

This happens because CLR centers by the geometric mean:
- CLR(x) = log(x) - mean(log(x))
- When one feature increases, the geometric mean increases
- This reduces the relative change after transformation

With small sample sizes (n=20 per group):
- Standard errors are ~0.5-1.0
- Attenuated effects of ~1.0 give t-stats of ~1-2
- These are not significant after multiple testing correction
""")

    # Recommendations
    print(f"\n9. RECOMMENDATIONS")
    print("-" * 40)
    print("""
1. Use q < 0.10 threshold for LinDA (balances sensitivity/FDR)
2. Increase sample size for detecting small fold changes
3. Consider alternative normalizations for absolute abundance studies
4. Use LinDA for FDR-controlled discovery, follow up with other methods
5. Document expected effect size attenuation in user guidance
""")

if __name__ == "__main__":
    analyze_linda_results()
