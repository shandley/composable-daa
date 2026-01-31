#!/usr/bin/env python3
"""Analyze power curves from benchmark results."""

import os
import sys

def load_ground_truth(path):
    """Load ground truth file and return set of differential features."""
    differential = set()
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                feature_id = parts[0]
                is_diff = parts[4].lower() == 'true'
                if is_diff:
                    differential.add(feature_id)
    return differential

def load_results(path):
    """Load results file and return set of significant features (q < 0.05)."""
    significant = set()
    all_features = set()
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                feature_id = parts[0]
                all_features.add(feature_id)
                try:
                    q_value = float(parts[6])
                    if q_value < 0.05:
                        significant.add(feature_id)
                except ValueError:
                    pass
    return significant, all_features

def calculate_metrics(differential, significant, total_tested):
    """Calculate sensitivity, FDR, and precision."""
    tp = len(differential & significant)
    fp = len(significant - differential)
    fn = len(differential - significant)

    # Only count features that were tested
    total_diff_tested = len(differential & total_tested)

    sensitivity = tp / total_diff_tested if total_diff_tested > 0 else 0
    fdr = fp / len(significant) if len(significant) > 0 else 0
    precision = tp / len(significant) if len(significant) > 0 else 0

    return {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'total_diff': total_diff_tested,
        'total_sig': len(significant),
        'sensitivity': sensitivity,
        'fdr': fdr,
        'precision': precision
    }

def main():
    base_dir = 'benchmark_results/power_analysis'
    effect_sizes = ['0.5', '1.0', '2.0', '4.0']
    models = ['linda', 'nb', 'zinb', 'bb', 'hurdle']

    print("=" * 80)
    print("POWER ANALYSIS RESULTS")
    print("=" * 80)
    print()

    # Collect all results for summary table
    all_results = {}

    for effect in effect_sizes:
        effect_dir = f'{base_dir}/effect_{effect}'
        gt_path = f'{effect_dir}/ground_truth.tsv'

        if not os.path.exists(gt_path):
            print(f"Missing ground truth for effect {effect}")
            continue

        differential = load_ground_truth(gt_path)

        print(f"Effect Size: {effect} log2FC ({2**float(effect):.1f}-fold change)")
        print(f"  True differential features: {len(differential)}")
        print()
        print(f"  {'Model':<10} {'TP':>4} {'FP':>4} {'FN':>4} {'Sens':>8} {'FDR':>8} {'Prec':>8}")
        print(f"  {'-'*10} {'-'*4} {'-'*4} {'-'*4} {'-'*8} {'-'*8} {'-'*8}")

        all_results[effect] = {}

        for model in models:
            results_path = f'{effect_dir}/{model}_results.tsv'
            if not os.path.exists(results_path):
                print(f"  {model:<10} -- missing results --")
                continue

            significant, all_features = load_results(results_path)
            metrics = calculate_metrics(differential, significant, all_features)
            all_results[effect][model] = metrics

            print(f"  {model:<10} {metrics['tp']:>4} {metrics['fp']:>4} {metrics['fn']:>4} "
                  f"{metrics['sensitivity']*100:>7.1f}% {metrics['fdr']*100:>7.1f}% "
                  f"{metrics['precision']*100:>7.1f}%")

        print()

    # Summary table: Sensitivity by effect size
    print("=" * 80)
    print("POWER CURVE SUMMARY (Sensitivity %)")
    print("=" * 80)
    print()
    print(f"{'Model':<10}", end="")
    for effect in effect_sizes:
        print(f" {effect:>8}", end="")
    print()
    print("-" * 50)

    for model in models:
        print(f"{model:<10}", end="")
        for effect in effect_sizes:
            if effect in all_results and model in all_results[effect]:
                sens = all_results[effect][model]['sensitivity'] * 100
                print(f" {sens:>7.1f}%", end="")
            else:
                print(f" {'--':>8}", end="")
        print()

    print()
    print("=" * 80)
    print("FDR SUMMARY (%)")
    print("=" * 80)
    print()
    print(f"{'Model':<10}", end="")
    for effect in effect_sizes:
        print(f" {effect:>8}", end="")
    print()
    print("-" * 50)

    for model in models:
        print(f"{model:<10}", end="")
        for effect in effect_sizes:
            if effect in all_results and model in all_results[effect]:
                fdr = all_results[effect][model]['fdr'] * 100
                print(f" {fdr:>7.1f}%", end="")
            else:
                print(f" {'--':>8}", end="")
        print()

if __name__ == '__main__':
    main()
