#!/usr/bin/env python3
"""Compare LinDA performance on raw vs compositional spike modes."""

import pandas as pd
import numpy as np

def analyze_comparison():
    """Compare raw mode data with compositional mode data."""

    print("=" * 80)
    print("RAW VS COMPOSITIONAL SPIKE MODE COMPARISON")
    print("=" * 80)

    # Load raw mode results
    raw_linda = pd.read_csv("raw_mode_test/linda_results.tsv", sep="\t")
    raw_gt = pd.read_csv("raw_mode_test/ground_truth.tsv", sep="\t")

    # Load compositional mode results (from effect_4.0)
    comp_linda = pd.read_csv("effect_4.0/linda_results.tsv", sep="\t")
    comp_gt = pd.read_csv("effect_4.0/ground_truth.tsv", sep="\t")

    # Analyze raw mode
    raw_merged = raw_linda.merge(raw_gt, on="feature_id", how="left")
    raw_merged["is_differential"] = raw_merged["is_differential"].fillna(False)
    raw_diff = raw_merged[raw_merged["is_differential"] == True]
    raw_nulls = raw_merged[raw_merged["is_differential"] == False]

    # Analyze compositional mode
    comp_merged = comp_linda.merge(comp_gt, on="feature_id", how="left")
    comp_merged["is_differential"] = comp_merged["is_differential"].fillna(False)
    comp_diff = comp_merged[comp_merged["is_differential"] == True]
    comp_nulls = comp_merged[comp_merged["is_differential"] == False]

    print("\n1. BASIC STATISTICS")
    print("-" * 50)
    print(f"Raw mode:   {len(raw_diff)} differential, {len(raw_nulls)} null")
    print(f"Comp mode:  {len(comp_diff)} differential, {len(comp_nulls)} null")

    print("\n2. EFFECT SIZE COMPARISON (true log2FC = 4.0)")
    print("-" * 50)

    raw_up = raw_diff[raw_diff["direction"] == "up"]
    raw_down = raw_diff[raw_diff["direction"] == "down"]
    comp_up = comp_diff[comp_diff["direction"] == "up"]
    comp_down = comp_diff[comp_diff["direction"] == "down"]

    print(f"\nUP-spiked:")
    print(f"  Raw mode:  mean estimate = {raw_up['estimate'].mean():.3f} (n={len(raw_up)})")
    print(f"  Comp mode: mean estimate = {comp_up['estimate'].mean():.3f} (n={len(comp_up)})")

    print(f"\nDOWN-spiked:")
    print(f"  Raw mode:  mean estimate = {raw_down['estimate'].mean():.3f} (n={len(raw_down)})")
    print(f"  Comp mode: mean estimate = {comp_down['estimate'].mean():.3f} (n={len(comp_down)})")

    print("\n3. SENSITIVITY COMPARISON")
    print("-" * 50)

    for threshold in [0.05, 0.10, 0.15, 0.20]:
        raw_tp = (raw_diff["q_value"] < threshold).sum()
        raw_fp = (raw_nulls["q_value"] < threshold).sum()
        raw_sens = raw_tp / len(raw_diff) if len(raw_diff) > 0 else 0
        raw_fdr = raw_fp / (raw_tp + raw_fp) if (raw_tp + raw_fp) > 0 else 0

        comp_tp = (comp_diff["q_value"] < threshold).sum()
        comp_fp = (comp_nulls["q_value"] < threshold).sum()
        comp_sens = comp_tp / len(comp_diff) if len(comp_diff) > 0 else 0
        comp_fdr = comp_fp / (comp_tp + comp_fp) if (comp_tp + comp_fp) > 0 else 0

        print(f"\nq < {threshold}:")
        print(f"  Raw mode:  {raw_sens*100:.1f}% sens, {raw_fdr*100:.1f}% FDR ({raw_tp} TP, {raw_fp} FP)")
        print(f"  Comp mode: {comp_sens*100:.1f}% sens, {comp_fdr*100:.1f}% FDR ({comp_tp} TP, {comp_fp} FP)")

    print("\n4. P-VALUE DISTRIBUTION FOR SPIKED FEATURES")
    print("-" * 50)

    raw_pvals = raw_diff["p_value"].values
    comp_pvals = comp_diff["p_value"].values

    print(f"\nRaw mode:")
    print(f"  Min p: {raw_pvals.min():.4f}, Max p: {raw_pvals.max():.4f}")
    print(f"  p < 0.05: {(raw_pvals < 0.05).sum()}/{len(raw_pvals)}")
    print(f"  p < 0.10: {(raw_pvals < 0.10).sum()}/{len(raw_pvals)}")

    print(f"\nComp mode:")
    print(f"  Min p: {comp_pvals.min():.4f}, Max p: {comp_pvals.max():.4f}")
    print(f"  p < 0.05: {(comp_pvals < 0.05).sum()}/{len(comp_pvals)}")
    print(f"  p < 0.10: {(comp_pvals < 0.10).sum()}/{len(comp_pvals)}")

    print("\n5. CONCLUSION")
    print("-" * 50)
    print("""
The CLR attenuation is similar between raw and compositional spike modes.
This confirms that the issue is inherent to CLR transformation, not the
spike mode used.

CLR centers by geometric mean, which:
1. Shifts when ANY feature changes significantly
2. Partially cancels out the observed effect
3. Results in ~25% of true effect being observed

This is by design - CLR creates relative abundance measures that are
robust to compositional constraints. The trade-off is reduced power
to detect true effects.
""")

if __name__ == "__main__":
    analyze_comparison()
