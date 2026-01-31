#!/usr/bin/env python3
"""Analyze the BB model FDR issue in detail."""

import os
import sys

def load_ground_truth(path):
    """Load ground truth with full details."""
    features = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                feat_id = parts[0]
                log2fc = float(parts[1])
                direction = parts[2]
                effect_type = parts[3]
                is_diff = parts[4].lower() == 'true'
                features[feat_id] = {
                    'true_log2fc': log2fc,
                    'direction': direction,
                    'effect_type': effect_type,
                    'is_differential': is_diff
                }
    return features

def load_results(path):
    """Load BB results with full details."""
    features = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                feat_id = parts[0]
                estimate = float(parts[2]) if parts[2] != 'NA' else None
                std_error = float(parts[3]) if parts[3] != 'NA' else None
                statistic = float(parts[4]) if parts[4] != 'NA' else None
                p_value = float(parts[5]) if parts[5] != 'NA' else None
                q_value = float(parts[6]) if parts[6] != 'NA' else None
                features[feat_id] = {
                    'estimate': estimate,
                    'std_error': std_error,
                    'statistic': statistic,
                    'p_value': p_value,
                    'q_value': q_value,
                    'significant': q_value is not None and q_value < 0.05
                }
    return features

def main():
    effect = sys.argv[1] if len(sys.argv) > 1 else '4.0'
    base_dir = 'benchmark_results/power_analysis'

    gt_path = f'{base_dir}/effect_{effect}/ground_truth.tsv'
    bb_path = f'{base_dir}/effect_{effect}/bb_results.tsv'

    ground_truth = load_ground_truth(gt_path)
    bb_results = load_results(bb_path)

    print(f"=" * 80)
    print(f"BB FDR ANALYSIS - Effect Size {effect} log2FC")
    print(f"=" * 80)
    print()

    # Categorize results
    true_pos = []  # Differential and significant
    false_pos = []  # Not differential but significant
    true_neg = []  # Not differential and not significant
    false_neg = []  # Differential but not significant

    for feat_id, gt in ground_truth.items():
        if feat_id not in bb_results:
            continue
        bb = bb_results[feat_id]

        is_diff = gt['is_differential']
        is_sig = bb['significant']

        entry = {
            'feature': feat_id,
            'true_log2fc': gt['true_log2fc'],
            'estimate': bb['estimate'],
            'std_error': bb['std_error'],
            'q_value': bb['q_value'],
            'direction': gt['direction']
        }

        if is_diff and is_sig:
            true_pos.append(entry)
        elif not is_diff and is_sig:
            false_pos.append(entry)
        elif not is_diff and not is_sig:
            true_neg.append(entry)
        else:
            false_neg.append(entry)

    print(f"Summary:")
    print(f"  True Positives:  {len(true_pos)}")
    print(f"  False Positives: {len(false_pos)}")
    print(f"  True Negatives:  {len(true_neg)}")
    print(f"  False Negatives: {len(false_neg)}")
    print(f"  FDR: {len(false_pos) / (len(true_pos) + len(false_pos)) * 100:.1f}%")
    print()

    # Analyze false positives
    print("=" * 80)
    print("FALSE POSITIVE ANALYSIS")
    print("=" * 80)
    print()

    # Direction of false positives
    fp_positive = [f for f in false_pos if f['estimate'] and f['estimate'] > 0]
    fp_negative = [f for f in false_pos if f['estimate'] and f['estimate'] < 0]

    print(f"Direction of false positives:")
    print(f"  Positive estimates (higher in treatment): {len(fp_positive)}")
    print(f"  Negative estimates (lower in treatment):  {len(fp_negative)}")
    print()

    # Compare to true positive directions
    tp_positive = [f for f in true_pos if f['estimate'] and f['estimate'] > 0]
    tp_negative = [f for f in true_pos if f['estimate'] and f['estimate'] < 0]

    print(f"Direction of true positives:")
    print(f"  Positive estimates: {len(tp_positive)}")
    print(f"  Negative estimates: {len(tp_negative)}")
    print()

    # Effect size comparison
    print("=" * 80)
    print("EFFECT SIZE COMPARISON")
    print("=" * 80)
    print()

    print("True Positives (top 10):")
    print(f"  {'Feature':<15} {'True log2FC':>12} {'Estimate':>12} {'SE':>10} {'q-value':>10}")
    print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
    for entry in sorted(true_pos, key=lambda x: abs(x['estimate'] or 0), reverse=True)[:10]:
        print(f"  {entry['feature']:<15} {entry['true_log2fc']:>12.2f} {entry['estimate']:>12.4f} "
              f"{entry['std_error']:>10.4f} {entry['q_value']:>10.2e}")

    print()
    print("False Positives (top 10 by |estimate|):")
    print(f"  {'Feature':<15} {'True log2FC':>12} {'Estimate':>12} {'SE':>10} {'q-value':>10}")
    print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
    for entry in sorted(false_pos, key=lambda x: abs(x['estimate'] or 0), reverse=True)[:10]:
        print(f"  {entry['feature']:<15} {entry['true_log2fc']:>12.2f} {entry['estimate']:>12.4f} "
              f"{entry['std_error']:>10.4f} {entry['q_value']:>10.2e}")

    print()
    print("=" * 80)
    print("STANDARD ERROR COMPARISON")
    print("=" * 80)
    print()

    tp_se = [f['std_error'] for f in true_pos if f['std_error']]
    fp_se = [f['std_error'] for f in false_pos if f['std_error']]
    tn_se = [f['std_error'] for f in true_neg if f['std_error']]

    print(f"Mean SE for True Positives:  {sum(tp_se)/len(tp_se):.4f}" if tp_se else "  N/A")
    print(f"Mean SE for False Positives: {sum(fp_se)/len(fp_se):.4f}" if fp_se else "  N/A")
    print(f"Mean SE for True Negatives:  {sum(tn_se)/len(tn_se):.4f}" if tn_se else "  N/A")

    tp_est = [abs(f['estimate']) for f in true_pos if f['estimate']]
    fp_est = [abs(f['estimate']) for f in false_pos if f['estimate']]
    tn_est = [abs(f['estimate']) for f in true_neg if f['estimate']]

    print()
    print(f"Mean |estimate| for True Positives:  {sum(tp_est)/len(tp_est):.4f}" if tp_est else "  N/A")
    print(f"Mean |estimate| for False Positives: {sum(fp_est)/len(fp_est):.4f}" if fp_est else "  N/A")
    print(f"Mean |estimate| for True Negatives:  {sum(tn_est)/len(tn_est):.4f}" if tn_est else "  N/A")

    print()
    print("=" * 80)
    print("COMPOSITIONAL ARTIFACT HYPOTHESIS")
    print("=" * 80)
    print()

    # Check if FP direction is opposite to overall TP direction
    # If we spike UP features, compositional closure means others go DOWN
    n_tp_up = len([f for f in true_pos if f['true_log2fc'] > 0])
    n_tp_down = len([f for f in true_pos if f['true_log2fc'] < 0])

    print(f"True differential features:")
    print(f"  Spiked UP:   {n_tp_up}")
    print(f"  Spiked DOWN: {n_tp_down}")
    print()

    # For each FP, check if its direction is consistent with compositional artifact
    # If a TP went UP, null features should go DOWN (compositional closure)
    print("Compositional artifact check:")
    print("  If this is compositional, FP direction should be OPPOSITE to net TP effect")

    net_tp_direction = 'UP' if n_tp_up > n_tp_down else 'DOWN' if n_tp_down > n_tp_up else 'BALANCED'
    print(f"  Net TP direction: {net_tp_direction}")

    expected_fp_direction = 'negative' if net_tp_direction == 'UP' else 'positive' if net_tp_direction == 'DOWN' else 'mixed'
    actual_fp_direction = 'positive' if len(fp_positive) > len(fp_negative) else 'negative' if len(fp_negative) > len(fp_positive) else 'mixed'

    print(f"  Expected FP direction (if compositional): {expected_fp_direction}")
    print(f"  Actual FP direction: {actual_fp_direction}")

    if expected_fp_direction == actual_fp_direction or (net_tp_direction == 'BALANCED'):
        print(f"  --> CONSISTENT with compositional artifact hypothesis")
    else:
        print(f"  --> NOT consistent with compositional artifact hypothesis")

if __name__ == '__main__':
    main()
