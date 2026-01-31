#!/usr/bin/env python3
"""Check if false positives have actually changed proportions in the data."""

import os
import sys

def load_counts(path):
    """Load count matrix as dict of feature -> sample -> count."""
    counts = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        sample_ids = header[1:]  # First column is feature_id
        for line in f:
            parts = line.strip().split('\t')
            feat_id = parts[0]
            counts[feat_id] = {}
            for i, count_str in enumerate(parts[1:]):
                counts[feat_id][sample_ids[i]] = int(count_str)
    return counts, sample_ids

def load_metadata(path):
    """Load metadata to get group assignments."""
    groups = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            sample_id = parts[0]
            group = parts[1]
            groups[sample_id] = group
    return groups

def load_ground_truth(path):
    """Load ground truth."""
    differential = set()
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5 and parts[4].lower() == 'true':
                differential.add(parts[0])
    return differential

def load_bb_results(path):
    """Load BB results."""
    results = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                feat_id = parts[0]
                try:
                    q_value = float(parts[6])
                    estimate = float(parts[2])
                except:
                    continue
                results[feat_id] = {
                    'q_value': q_value,
                    'estimate': estimate,
                    'significant': q_value < 0.05
                }
    return results

def main():
    effect = sys.argv[1] if len(sys.argv) > 1 else '4.0'
    base_dir = f'benchmark_results/power_analysis/effect_{effect}'

    counts, sample_ids = load_counts(f'{base_dir}/counts.tsv')
    groups = load_metadata(f'{base_dir}/metadata.tsv')
    differential = load_ground_truth(f'{base_dir}/ground_truth.tsv')
    bb_results = load_bb_results(f'{base_dir}/bb_results.tsv')

    # Calculate library sizes
    lib_sizes = {}
    for sample in sample_ids:
        lib_sizes[sample] = sum(counts[feat][sample] for feat in counts)

    # Calculate proportions for each feature
    def calc_group_mean_prop(feat_id, group_name):
        group_samples = [s for s in sample_ids if groups.get(s) == group_name]
        props = []
        for s in group_samples:
            count = counts[feat_id].get(s, 0)
            lib = lib_sizes[s]
            if lib > 0:
                props.append(count / lib)
        return sum(props) / len(props) if props else 0

    # Find false positives
    false_positives = []
    for feat_id, res in bb_results.items():
        if feat_id in differential:
            continue  # Skip true positives
        if res['significant']:
            control_prop = calc_group_mean_prop(feat_id, 'control')
            treatment_prop = calc_group_mean_prop(feat_id, 'treatment')
            log2fc = 0
            if control_prop > 0 and treatment_prop > 0:
                log2fc = (treatment_prop / control_prop)
                import math
                log2fc = math.log2(log2fc) if log2fc > 0 else 0

            false_positives.append({
                'feature': feat_id,
                'control_prop': control_prop,
                'treatment_prop': treatment_prop,
                'actual_log2fc': log2fc,
                'bb_estimate': res['estimate'],
                'q_value': res['q_value']
            })

    print(f"False Positive Proportion Analysis - Effect Size {effect}")
    print("=" * 90)
    print()
    print(f"Found {len(false_positives)} false positives")
    print()

    # Sort by absolute actual log2fc
    false_positives.sort(key=lambda x: abs(x['actual_log2fc']), reverse=True)

    print(f"{'Feature':<15} {'Ctrl Prop':>12} {'Treat Prop':>12} {'Actual log2FC':>14} {'BB Est':>10}")
    print("-" * 90)

    for fp in false_positives[:20]:
        print(f"{fp['feature']:<15} {fp['control_prop']:>12.6f} {fp['treatment_prop']:>12.6f} "
              f"{fp['actual_log2fc']:>14.4f} {fp['bb_estimate']:>10.4f}")

    print()
    print("Analysis:")

    # Check how many FPs have actual proportion changes
    significant_changes = [fp for fp in false_positives if abs(fp['actual_log2fc']) > 0.5]
    moderate_changes = [fp for fp in false_positives if 0.1 < abs(fp['actual_log2fc']) <= 0.5]
    small_changes = [fp for fp in false_positives if abs(fp['actual_log2fc']) <= 0.1]

    print(f"  FPs with |actual log2FC| > 0.5:  {len(significant_changes)}")
    print(f"  FPs with 0.1 < |actual log2FC| <= 0.5: {len(moderate_changes)}")
    print(f"  FPs with |actual log2FC| <= 0.1: {len(small_changes)}")
    print()

    # Direction check
    neg_actual = [fp for fp in false_positives if fp['actual_log2fc'] < -0.1]
    pos_actual = [fp for fp in false_positives if fp['actual_log2fc'] > 0.1]

    print(f"  FPs with negative actual change: {len(neg_actual)} (decreased in treatment)")
    print(f"  FPs with positive actual change: {len(pos_actual)} (increased in treatment)")
    print()
    print("Conclusion: If FPs have REAL proportion changes (due to compositional closure),")
    print("            they are detected correctly by BB but shouldn't be called 'false' positives!")

if __name__ == '__main__':
    main()
