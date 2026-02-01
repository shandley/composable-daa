#!/usr/bin/env python3
"""
Random-Effects Meta-Analysis for Microbiome Reproducibility Atlas

Implements DerSimonian-Laird random-effects meta-analysis to pool
effect sizes across cohorts.
"""

import pandas as pd
import numpy as np
from scipy import stats
import os
import sys


def random_effects_meta(effects, variances):
    """
    Random-effects meta-analysis using DerSimonian-Laird estimator.

    Parameters:
    - effects: list of effect sizes (log scale)
    - variances: list of variance estimates (SE^2)

    Returns:
    - dict with pooled effect, CI, heterogeneity stats
    """
    effects = np.array(effects)
    variances = np.array(variances)

    # Remove any NaN/inf values
    valid = np.isfinite(effects) & np.isfinite(variances) & (variances > 0)
    if valid.sum() < 2:
        return None

    effects = effects[valid]
    variances = variances[valid]
    weights = 1 / variances

    # Fixed-effect estimate (for Q calculation)
    fe_effect = np.sum(weights * effects) / np.sum(weights)

    # Cochran's Q
    Q = np.sum(weights * (effects - fe_effect)**2)
    df = len(effects) - 1

    # DerSimonian-Laird tau^2
    C = np.sum(weights) - np.sum(weights**2) / np.sum(weights)
    tau2 = max(0, (Q - df) / C) if C > 0 else 0

    # Random-effects weights
    re_weights = 1 / (variances + tau2)

    # Pooled estimate
    pooled_effect = np.sum(re_weights * effects) / np.sum(re_weights)
    pooled_var = 1 / np.sum(re_weights)
    pooled_se = np.sqrt(pooled_var)

    # I^2 heterogeneity
    I2 = max(0, (Q - df) / Q) if Q > df else 0

    # P-value
    z = pooled_effect / pooled_se if pooled_se > 0 else 0
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    return {
        'effect': pooled_effect,
        'se': pooled_se,
        'ci_lower': pooled_effect - 1.96 * pooled_se,
        'ci_upper': pooled_effect + 1.96 * pooled_se,
        'z': z,
        'p_value': p_value,
        'tau2': tau2,
        'I2': I2,
        'Q': Q,
        'Q_pvalue': 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0,
        'n_studies': len(effects)
    }


def format_forest_row(cohort, effect, se, weight_pct):
    """Format a single row for forest plot display."""
    ci_lower = effect - 1.96 * se
    ci_upper = effect + 1.96 * se

    # Create ASCII forest plot bar
    bar_center = 30  # center position
    bar_scale = 10   # pixels per unit
    pos = int(bar_center + effect * bar_scale)
    pos = max(5, min(55, pos))

    bar = ['-'] * 60
    bar[bar_center] = '|'  # zero line
    bar[pos] = '●'

    return f"{cohort:15} {''.join(bar[:35])} {effect:6.2f} [{ci_lower:5.2f}, {ci_upper:5.2f}] {weight_pct:5.1f}%"


def forest_plot(taxon, cohort_data, meta_result):
    """Generate ASCII forest plot."""
    lines = []
    lines.append(f"\nTaxon: {taxon}")
    lines.append("─" * 70)
    lines.append(f"{'Cohort':15} {'Effect Plot':35} Effect [95% CI]      Weight")
    lines.append("─" * 70)

    total_weight = sum(1/d['variance'] for d in cohort_data if d['variance'] > 0)

    for d in cohort_data:
        if d['variance'] > 0:
            weight_pct = 100 * (1/d['variance']) / total_weight
            lines.append(format_forest_row(d['cohort'], d['effect'], d['se'], weight_pct))

    lines.append("─" * 70)

    if meta_result:
        bar = ['-'] * 60
        bar[30] = '|'
        pos = int(30 + meta_result['effect'] * 10)
        pos = max(5, min(55, pos))
        bar[pos] = '◆'

        lines.append(f"{'Pooled (RE)':15} {''.join(bar[:35])} {meta_result['effect']:6.2f} [{meta_result['ci_lower']:5.2f}, {meta_result['ci_upper']:5.2f}]")
        lines.append("─" * 70)
        lines.append(f"Heterogeneity: I² = {100*meta_result['I2']:.0f}%, Q = {meta_result['Q']:.1f}, p = {meta_result['Q_pvalue']:.3f}")
        lines.append(f"Test for effect: z = {meta_result['z']:.2f}, p = {meta_result['p_value']:.4f}")

    return "\n".join(lines)


def run_fusobacterium_meta(results_dir):
    """
    Run meta-analysis for Fusobacterium across CRC cohorts.
    """
    print("=" * 70)
    print("META-ANALYSIS: Fusobacterium in CRC")
    print("=" * 70)

    cohorts = ['zeller', 'zackular', 'xiang', 'zhao']
    cohort_data = []

    for cohort in cohorts:
        result_file = os.path.join(results_dir, f'{cohort}_hurdle.tsv')
        if not os.path.exists(result_file):
            continue

        df = pd.read_csv(result_file, sep='\t')

        # Find Fusobacterium
        fuso = df[df['feature_id'].str.contains('Fusobacterium', case=False, na=False)]

        if len(fuso) > 0:
            row = fuso.iloc[0]
            effect = row['estimate']
            se = row['std_error']

            if pd.notna(effect) and pd.notna(se) and se > 0:
                cohort_data.append({
                    'cohort': cohort,
                    'effect': effect,
                    'se': se,
                    'variance': se**2,
                    'p_value': row['p_value'],
                    'q_value': row['q_value']
                })
                print(f"\n{cohort}:")
                print(f"  Effect: {effect:.3f} (SE: {se:.3f})")
                print(f"  p-value: {row['p_value']:.4f}, q-value: {row['q_value']:.4f}")
            else:
                print(f"\n{cohort}: Invalid effect or SE")
        else:
            print(f"\n{cohort}: Fusobacterium not found")

    if len(cohort_data) < 2:
        print("\nNot enough cohorts for meta-analysis")
        return

    # Run meta-analysis
    effects = [d['effect'] for d in cohort_data]
    variances = [d['variance'] for d in cohort_data]

    result = random_effects_meta(effects, variances)

    if result:
        print("\n" + "=" * 70)
        print("POOLED RESULT")
        print("=" * 70)
        print(f"Pooled effect: {result['effect']:.3f} (SE: {result['se']:.3f})")
        print(f"95% CI: [{result['ci_lower']:.3f}, {result['ci_upper']:.3f}]")
        print(f"z = {result['z']:.2f}, p = {result['p_value']:.4f}")
        print(f"\nHeterogeneity:")
        print(f"  τ² = {result['tau2']:.3f}")
        print(f"  I² = {100*result['I2']:.1f}%")
        print(f"  Q = {result['Q']:.2f}, p = {result['Q_pvalue']:.3f}")

        # Interpretation
        print("\n" + "=" * 70)
        print("INTERPRETATION")
        print("=" * 70)

        if result['effect'] < 0:
            direction = "enriched in CRC (higher than controls)"
        else:
            direction = "depleted in CRC (lower than controls)"

        print(f"Direction: Fusobacterium is {direction}")

        if result['p_value'] < 0.05:
            print(f"Significance: SIGNIFICANT at p < 0.05")
        elif result['p_value'] < 0.10:
            print(f"Significance: Suggestive (p < 0.10)")
        else:
            print(f"Significance: Not significant (p = {result['p_value']:.3f})")

        if result['I2'] < 0.25:
            het = "LOW - effects consistent across cohorts"
        elif result['I2'] < 0.50:
            het = "MODERATE - some variation"
        elif result['I2'] < 0.75:
            het = "SUBSTANTIAL - consider subgroup analysis"
        else:
            het = "HIGH - studies may be measuring different things"
        print(f"Heterogeneity: {het}")

        # Forest plot
        print(forest_plot("Fusobacterium", cohort_data, result))

    return result


def run_known_positives_meta(results_dir):
    """
    Run meta-analysis for all known CRC positives.
    """
    from taxonomy import KNOWN_POSITIVES

    print("\n\n" + "=" * 70)
    print("META-ANALYSIS: All Known CRC Markers")
    print("=" * 70)

    cohorts = ['zeller', 'zackular', 'xiang', 'zhao']

    # Load all result files
    results = {}
    for cohort in cohorts:
        result_file = os.path.join(results_dir, f'{cohort}_hurdle.tsv')
        if os.path.exists(result_file):
            results[cohort] = pd.read_csv(result_file, sep='\t')

    # For each known positive
    summary = []
    for kp in KNOWN_POSITIVES['crc']:
        taxon = kp['taxon']
        expected_dir = kp['direction']

        print(f"\n--- {taxon} (expected: {expected_dir}) ---")

        cohort_data = []
        for cohort, df in results.items():
            matches = df[df['feature_id'].str.contains(taxon, case=False, na=False)]
            if len(matches) > 0:
                row = matches.iloc[0]
                effect = row['estimate']
                se = row['std_error']

                if pd.notna(effect) and pd.notna(se) and se > 0 and se < 1000:
                    cohort_data.append({
                        'cohort': cohort,
                        'effect': effect,
                        'se': se,
                        'variance': se**2
                    })
                    print(f"  {cohort}: effect={effect:.2f}, SE={se:.2f}")

        if len(cohort_data) >= 2:
            effects = [d['effect'] for d in cohort_data]
            variances = [d['variance'] for d in cohort_data]

            result = random_effects_meta(effects, variances)
            if result:
                observed_dir = 'up' if result['effect'] < 0 else 'down'
                correct = observed_dir == expected_dir
                sig = result['p_value'] < 0.10

                status = "✓" if (sig and correct) else "~" if correct else "✗"

                print(f"  POOLED: effect={result['effect']:.2f}, p={result['p_value']:.4f}")
                print(f"  Direction: observed={observed_dir}, expected={expected_dir} {'✓' if correct else '✗'}")
                print(f"  Status: {status}")

                summary.append({
                    'taxon': taxon,
                    'expected': expected_dir,
                    'observed': observed_dir,
                    'pooled_effect': result['effect'],
                    'p_value': result['p_value'],
                    'I2': result['I2'],
                    'n_cohorts': result['n_studies'],
                    'correct_direction': correct,
                    'significant': sig
                })
        else:
            print(f"  Not enough cohorts for meta-analysis")

    # Summary table
    if summary:
        print("\n\n" + "=" * 70)
        print("SUMMARY: Known Positive Validation via Meta-Analysis")
        print("=" * 70)

        n_correct = sum(1 for s in summary if s['correct_direction'])
        n_sig = sum(1 for s in summary if s['significant'])
        n_both = sum(1 for s in summary if s['correct_direction'] and s['significant'])

        print(f"\nTested: {len(summary)} known markers")
        print(f"Correct direction: {n_correct}/{len(summary)} ({100*n_correct/len(summary):.0f}%)")
        print(f"Significant (p<0.10): {n_sig}/{len(summary)} ({100*n_sig/len(summary):.0f}%)")
        print(f"Both correct + significant: {n_both}/{len(summary)} ({100*n_both/len(summary):.0f}%)")

        print("\nPer-taxon results:")
        for s in summary:
            status = "✓" if (s['correct_direction'] and s['significant']) else "~" if s['correct_direction'] else "✗"
            print(f"  {status} {s['taxon']}: pooled effect={s['pooled_effect']:.2f}, p={s['p_value']:.4f}, I²={100*s['I2']:.0f}%")

    return summary


if __name__ == '__main__':
    # Add parent directory for imports
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    results_dir = os.path.join(os.path.dirname(__file__), '..', 'diseases', 'crc', 'results_genus')
    results_dir = os.path.abspath(results_dir)

    # Run Fusobacterium meta-analysis
    fuso_result = run_fusobacterium_meta(results_dir)

    # Run all known positives meta-analysis
    summary = run_known_positives_meta(results_dir)
