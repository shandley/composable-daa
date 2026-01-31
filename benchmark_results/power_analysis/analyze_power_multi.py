#!/usr/bin/env python3
"""Analyze power curves from benchmark results with multiple thresholds."""

import os

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

def load_results(path, q_threshold=0.05):
    """Load results file and return set of significant features."""
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
                    if q_value < q_threshold:
                        significant.add(feature_id)
                except ValueError:
                    pass
    return significant, all_features

def calculate_metrics(differential, significant, total_tested):
    """Calculate sensitivity, FDR, and precision."""
    tp = len(differential & significant)
    fp = len(significant - differential)
    fn = len(differential - significant)
    total_diff_tested = len(differential & total_tested)

    sensitivity = tp / total_diff_tested if total_diff_tested > 0 else 0
    fdr = fp / len(significant) if len(significant) > 0 else 0
    precision = tp / len(significant) if len(significant) > 0 else 0

    return {
        'tp': tp, 'fp': fp, 'fn': fn,
        'total_diff': total_diff_tested,
        'total_sig': len(significant),
        'sensitivity': sensitivity,
        'fdr': fdr, 'precision': precision
    }

def main():
    base_dir = 'benchmark_results/power_analysis'
    effect_sizes = ['0.5', '1.0', '2.0', '4.0']
    models = ['linda', 'nb', 'zinb', 'bb', 'hurdle']
    thresholds = [0.05, 0.10, 0.20]

    for q_thresh in thresholds:
        print("=" * 80)
        print(f"POWER CURVE SUMMARY at q < {q_thresh} (Sensitivity %)")
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
                effect_dir = f'{base_dir}/effect_{effect}'
                gt_path = f'{effect_dir}/ground_truth.tsv'
                results_path = f'{effect_dir}/{model}_results.tsv'

                if not os.path.exists(results_path):
                    print(f" {'--':>8}", end="")
                    continue

                differential = load_ground_truth(gt_path)
                significant, all_features = load_results(results_path, q_thresh)
                metrics = calculate_metrics(differential, significant, all_features)
                sens = metrics['sensitivity'] * 100
                print(f" {sens:>7.1f}%", end="")
            print()

        print()
        print(f"FDR at q < {q_thresh} (%)")
        print("-" * 50)

        for model in models:
            print(f"{model:<10}", end="")
            for effect in effect_sizes:
                effect_dir = f'{base_dir}/effect_{effect}'
                gt_path = f'{effect_dir}/ground_truth.tsv'
                results_path = f'{effect_dir}/{model}_results.tsv'

                if not os.path.exists(results_path):
                    print(f" {'--':>8}", end="")
                    continue

                differential = load_ground_truth(gt_path)
                significant, all_features = load_results(results_path, q_thresh)
                metrics = calculate_metrics(differential, significant, all_features)
                fdr = metrics['fdr'] * 100
                n_sig = metrics['total_sig']
                if n_sig > 0:
                    print(f" {fdr:>7.1f}%", end="")
                else:
                    print(f" {'n/a':>8}", end="")
            print()
        print()

if __name__ == '__main__':
    main()
