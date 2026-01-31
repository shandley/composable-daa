#!/usr/bin/env python3
"""Check how LinDA handles the same features that BB calls significant."""

import os
import math

def load_counts(path):
    """Load count matrix."""
    counts = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        sample_ids = header[1:]
        for line in f:
            parts = line.strip().split('\t')
            feat_id = parts[0]
            counts[feat_id] = {sample_ids[i]: int(parts[i+1]) for i in range(len(sample_ids))}
    return counts, sample_ids

def load_metadata(path):
    """Load metadata."""
    groups = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            groups[parts[0]] = parts[1]
    return groups

def load_results(path):
    """Load results file."""
    results = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 7:
                try:
                    results[parts[0]] = {
                        'estimate': float(parts[2]),
                        'p_value': float(parts[5]),
                        'q_value': float(parts[6]),
                        'significant': float(parts[6]) < 0.05
                    }
                except:
                    pass
    return results

def calc_clr(counts, sample_ids, pseudocount=0.5):
    """Calculate CLR-transformed values."""
    clr = {}
    for sample in sample_ids:
        # Get all counts for this sample
        sample_counts = [counts[f].get(sample, 0) + pseudocount for f in counts]
        # Geometric mean
        log_sum = sum(math.log(c) for c in sample_counts)
        geo_mean = math.exp(log_sum / len(sample_counts))
        # CLR transform
        for feat in counts:
            if feat not in clr:
                clr[feat] = {}
            count = counts[feat].get(sample, 0) + pseudocount
            clr[feat][sample] = math.log(count / geo_mean)
    return clr

def main():
    base_dir = 'benchmark_results/power_analysis/effect_4.0'

    counts, sample_ids = load_counts(f'{base_dir}/counts.tsv')
    groups = load_metadata(f'{base_dir}/metadata.tsv')
    bb_results = load_results(f'{base_dir}/bb_results.tsv')
    linda_results = load_results(f'{base_dir}/linda_results.tsv')

    # Calculate CLR values
    clr = calc_clr(counts, sample_ids)

    # Library sizes
    lib_sizes = {s: sum(counts[f].get(s, 0) for f in counts) for s in sample_ids}

    # Get features that BB calls significant but aren't truly differential
    # (from earlier analysis we know Feature_0061, Feature_0107, etc.)
    bb_fps = ['Feature_0189', 'Feature_0061', 'Feature_0107', 'Feature_0048',
              'Feature_0166', 'Feature_0120', 'Feature_0073', 'Feature_0032']

    print("Comparison: BB vs LinDA on BB's 'false positives'")
    print("=" * 100)
    print()
    print("These features have REAL proportion changes (compositional artifacts)")
    print("BB detects them; does LinDA?")
    print()

    ctrl_samples = [s for s in sample_ids if groups[s] == 'control']
    treat_samples = [s for s in sample_ids if groups[s] == 'treatment']

    print(f"{'Feature':<15} {'Prop Change':>12} {'BB q':>10} {'BB Sig':>8} {'LinDA q':>10} {'LinDA Sig':>10} {'CLR Diff':>10}")
    print("-" * 100)

    for feat in bb_fps:
        # Proportion change
        ctrl_prop = sum(counts[feat].get(s, 0) / lib_sizes[s] for s in ctrl_samples) / len(ctrl_samples)
        treat_prop = sum(counts[feat].get(s, 0) / lib_sizes[s] for s in treat_samples) / len(treat_samples)
        prop_change = math.log2(treat_prop / ctrl_prop) if ctrl_prop > 0 and treat_prop > 0 else 0

        # CLR difference
        ctrl_clr = sum(clr[feat].get(s, 0) for s in ctrl_samples) / len(ctrl_samples)
        treat_clr = sum(clr[feat].get(s, 0) for s in treat_samples) / len(treat_samples)
        clr_diff = treat_clr - ctrl_clr

        bb = bb_results.get(feat, {})
        linda = linda_results.get(feat, {})

        bb_q = f"{bb.get('q_value', 'NA'):.2e}" if bb.get('q_value') else 'NA'
        bb_sig = 'YES' if bb.get('significant') else 'no'
        linda_q = f"{linda.get('q_value', 'NA'):.2e}" if linda.get('q_value') else 'NA'
        linda_sig = 'YES' if linda.get('significant') else 'no'

        print(f"{feat:<15} {prop_change:>12.2f} {bb_q:>10} {bb_sig:>8} {linda_q:>10} {linda_sig:>10} {clr_diff:>10.2f}")

    print()
    print("Interpretation:")
    print("  - 'Prop Change' = log2(treatment_prop / control_prop) - REAL proportion change")
    print("  - 'CLR Diff' = mean(treatment_CLR) - mean(control_CLR) - compositionally-corrected")
    print()
    print("  If CLR Diff ≈ 0, the feature is stable after compositional correction")
    print("  LinDA uses CLR, so it should NOT detect these compositional artifacts")

if __name__ == '__main__':
    main()
