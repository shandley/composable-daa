#!/usr/bin/env python3
"""Compare FDR patterns across models."""

import os

def load_ground_truth(path):
    """Load ground truth."""
    features = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                feat_id = parts[0]
                log2fc = float(parts[1])
                is_diff = parts[4].lower() == 'true'
                features[feat_id] = {'true_log2fc': log2fc, 'is_differential': is_diff}
    return features

def load_results(path):
    """Load results."""
    features = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                feat_id = parts[0]
                try:
                    estimate = float(parts[2]) if parts[2] not in ['NA', 'NaN'] else None
                    q_value = float(parts[6]) if parts[6] not in ['NA', 'NaN'] else None
                except:
                    estimate = None
                    q_value = None
                features[feat_id] = {
                    'estimate': estimate,
                    'q_value': q_value,
                    'significant': q_value is not None and q_value < 0.05
                }
    return features

def analyze_model(gt, results):
    """Analyze FDR pattern for a model."""
    fp_pos = 0
    fp_neg = 0
    tp = 0
    fp = 0

    for feat_id, truth in gt.items():
        if feat_id not in results:
            continue
        res = results[feat_id]

        is_diff = truth['is_differential']
        is_sig = res['significant']
        estimate = res['estimate']

        if is_diff and is_sig:
            tp += 1
        elif not is_diff and is_sig:
            fp += 1
            if estimate and estimate > 0:
                fp_pos += 1
            elif estimate and estimate < 0:
                fp_neg += 1

    fdr = fp / (tp + fp) * 100 if (tp + fp) > 0 else 0
    return {
        'tp': tp, 'fp': fp, 'fdr': fdr,
        'fp_pos': fp_pos, 'fp_neg': fp_neg
    }

def main():
    base_dir = 'benchmark_results/power_analysis'
    effect = '4.0'
    models = ['linda', 'nb', 'zinb', 'bb', 'hurdle']

    gt_path = f'{base_dir}/effect_{effect}/ground_truth.tsv'
    gt = load_ground_truth(gt_path)

    print(f"Model Comparison at {effect} log2FC")
    print("=" * 70)
    print()
    print(f"{'Model':<10} {'TP':>4} {'FP':>4} {'FDR':>8} {'FP+':>6} {'FP-':>6} {'Ratio':>10}")
    print("-" * 70)

    for model in models:
        results_path = f'{base_dir}/effect_{effect}/{model}_results.tsv'
        if not os.path.exists(results_path):
            continue

        results = load_results(results_path)
        stats = analyze_model(gt, results)

        ratio = f"{stats['fp_pos']}:{stats['fp_neg']}" if stats['fp'] > 0 else "n/a"
        fdr_str = f"{stats['fdr']:.1f}%" if stats['fdr'] > 0 else "n/a"

        print(f"{model:<10} {stats['tp']:>4} {stats['fp']:>4} {fdr_str:>8} "
              f"{stats['fp_pos']:>6} {stats['fp_neg']:>6} {ratio:>10}")

    print()
    print("FP+ = False positives with positive estimate (higher in treatment)")
    print("FP- = False positives with negative estimate (lower in treatment)")
    print()
    print("If compositional artifact is the cause, FP- >> FP+ (since net true effect is UP)")

if __name__ == '__main__':
    main()
