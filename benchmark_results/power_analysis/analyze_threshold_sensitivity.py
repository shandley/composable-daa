#!/usr/bin/env python3
"""Analyze LinDA sensitivity at different q-value thresholds."""

import pandas as pd
import numpy as np

def analyze_thresholds():
    """Compare sensitivity and FDR at different q-value thresholds."""

    print("=" * 80)
    print("LINDA SENSITIVITY BY Q-VALUE THRESHOLD")
    print("=" * 80)

    effect_sizes = ["0.5", "1.0", "2.0", "4.0"]
    thresholds = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]

    results = []

    for effect in effect_sizes:
        try:
            linda = pd.read_csv(f"effect_{effect}/linda_results.tsv", sep="\t")
            gt = pd.read_csv(f"effect_{effect}/ground_truth.tsv", sep="\t")

            merged = linda.merge(gt, on="feature_id", how="left")
            merged["is_differential"] = merged["is_differential"].fillna(False)

            n_true_diff = merged["is_differential"].sum()
            n_null = (~merged["is_differential"]).sum()

            for threshold in thresholds:
                sig = merged[merged["q_value"] < threshold]
                tp = sig["is_differential"].sum()
                fp = (~sig["is_differential"]).sum()

                sensitivity = tp / n_true_diff if n_true_diff > 0 else 0
                fdr = fp / len(sig) if len(sig) > 0 else 0

                results.append({
                    "effect_size": effect,
                    "threshold": threshold,
                    "n_significant": len(sig),
                    "tp": tp,
                    "fp": fp,
                    "sensitivity": sensitivity,
                    "fdr": fdr
                })

        except Exception as e:
            print(f"Error for effect {effect}: {e}")

    # Print results table
    print("\nSensitivity (True Positive Rate) by threshold:")
    print("-" * 70)
    print(f"{'Threshold':<12}", end="")
    for effect in effect_sizes:
        print(f"{effect} log2FC".rjust(12), end="")
    print()
    print("-" * 70)

    for threshold in thresholds:
        print(f"q < {threshold:<8}", end="")
        for effect in effect_sizes:
            row = [r for r in results if r["effect_size"] == effect and r["threshold"] == threshold]
            if row:
                print(f"{row[0]['sensitivity']*100:>10.1f}%", end="")
            else:
                print(f"{'n/a':>11}", end="")
        print()

    print("\nFalse Discovery Rate by threshold:")
    print("-" * 70)
    print(f"{'Threshold':<12}", end="")
    for effect in effect_sizes:
        print(f"{effect} log2FC".rjust(12), end="")
    print()
    print("-" * 70)

    for threshold in thresholds:
        print(f"q < {threshold:<8}", end="")
        for effect in effect_sizes:
            row = [r for r in results if r["effect_size"] == effect and r["threshold"] == threshold]
            if row:
                if row[0]['n_significant'] > 0:
                    print(f"{row[0]['fdr']*100:>10.1f}%", end="")
                else:
                    print(f"{'n/a':>11}", end="")
            else:
                print(f"{'n/a':>11}", end="")
        print()

    print("\nNumber of Discoveries by threshold:")
    print("-" * 70)
    print(f"{'Threshold':<12}", end="")
    for effect in effect_sizes:
        print(f"{effect} log2FC".rjust(12), end="")
    print()
    print("-" * 70)

    for threshold in thresholds:
        print(f"q < {threshold:<8}", end="")
        for effect in effect_sizes:
            row = [r for r in results if r["effect_size"] == effect and r["threshold"] == threshold]
            if row:
                print(f"{row[0]['tp']}/{row[0]['n_significant']:>9}", end="")
            else:
                print(f"{'n/a':>11}", end="")
        print()

    # Summary and recommendations
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    # Best threshold analysis
    print("\nBest threshold for each effect size (maximizing F1-score):")
    print("-" * 50)

    for effect in effect_sizes:
        best_f1 = 0
        best_thresh = None
        best_stats = None

        for threshold in thresholds:
            row = [r for r in results if r["effect_size"] == effect and r["threshold"] == threshold]
            if row and row[0]['n_significant'] > 0:
                sens = row[0]['sensitivity']
                prec = 1 - row[0]['fdr']
                if sens + prec > 0:
                    f1 = 2 * sens * prec / (sens + prec)
                    if f1 > best_f1:
                        best_f1 = f1
                        best_thresh = threshold
                        best_stats = row[0]

        if best_thresh:
            print(f"  {effect} log2FC: q < {best_thresh} → {best_stats['sensitivity']*100:.1f}% sens, {best_stats['fdr']*100:.1f}% FDR, F1={best_f1:.3f}")
        else:
            print(f"  {effect} log2FC: No discoveries at any threshold")

    # Compare with ZINB at q<0.10
    print("\n\nComparison with ZINB at q < 0.10:")
    print("-" * 50)

    for effect in effect_sizes:
        try:
            zinb = pd.read_csv(f"effect_{effect}/zinb_results.tsv", sep="\t")
            gt = pd.read_csv(f"effect_{effect}/ground_truth.tsv", sep="\t")
            merged = zinb.merge(gt, on="feature_id", how="left")
            merged["is_differential"] = merged["is_differential"].fillna(False)

            n_true_diff = merged["is_differential"].sum()
            sig = merged[merged["q_value"] < 0.10]
            tp = sig["is_differential"].sum()
            fp = (~sig["is_differential"]).sum()
            fdr = fp / len(sig) if len(sig) > 0 else 0

            linda_row = [r for r in results if r["effect_size"] == effect and r["threshold"] == 0.10]
            if linda_row:
                print(f"  {effect} log2FC: LinDA {linda_row[0]['sensitivity']*100:.1f}% sens, {linda_row[0]['fdr']*100:.1f}% FDR")
                print(f"            ZINB  {tp/n_true_diff*100:.1f}% sens, {fdr*100:.1f}% FDR")
        except:
            pass

if __name__ == "__main__":
    analyze_thresholds()
