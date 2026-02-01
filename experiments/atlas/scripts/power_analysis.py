#!/usr/bin/env python3
"""
Power Analysis for Microbiome Reproducibility Atlas

Calculates power and minimum detectable effect sizes for each cohort.
"""

import numpy as np
from scipy import stats
import pandas as pd
import os


def minimum_detectable_effect(n_case, n_control, alpha=0.05, power=0.80):
    """
    Calculate minimum detectable effect size (Cohen's d).

    For sparse count data, effects are on log scale, so:
    - d = 0.5 (medium) ≈ 1.6x fold change
    - d = 0.8 (large) ≈ 2.2x fold change
    - d = 1.0 ≈ 2.7x fold change
    """
    n_total = n_case + n_control
    harmonic_n = 2 * n_case * n_control / n_total

    # Effect size (Cohen's d) detectable at given power
    # Using approximation: d ≈ (z_alpha + z_beta) / sqrt(n/2)
    z_alpha = stats.norm.ppf(1 - alpha/2)
    z_beta = stats.norm.ppf(power)

    d = (z_alpha + z_beta) / np.sqrt(harmonic_n / 2)

    return d


def required_sample_size(effect_size, alpha=0.05, power=0.80, ratio=1.0):
    """
    Calculate required sample size per group to detect effect.

    effect_size: Cohen's d
    ratio: n_case/n_control ratio (1.0 = balanced)
    """
    z_alpha = stats.norm.ppf(1 - alpha/2)
    z_beta = stats.norm.ppf(power)

    # n per group for equal allocation
    n = 2 * ((z_alpha + z_beta) / effect_size) ** 2

    if ratio != 1.0:
        # Adjust for unequal allocation
        n = n * (1 + ratio) ** 2 / (4 * ratio)

    return int(np.ceil(n))


def effect_size_interpretation(d):
    """Interpret Cohen's d in fold-change terms."""
    # Assuming log-normal distribution, d maps approximately to fold change
    # d = log(FC) / pooled_sd, where pooled_sd ≈ 1 for normalized data
    fc = np.exp(d)
    if d < 0.2:
        return f"negligible (< 1.2x)"
    elif d < 0.5:
        return f"small (~{fc:.1f}x)"
    elif d < 0.8:
        return f"medium (~{fc:.1f}x)"
    else:
        return f"large (~{fc:.1f}x)"


def analyze_cohort_power(data_dir):
    """Analyze power for each cohort."""
    print("=" * 70)
    print("POWER ANALYSIS: CRC Cohorts")
    print("=" * 70)

    cohorts = ['zeller', 'zackular', 'xiang', 'zhao']
    results = []

    for cohort in cohorts:
        meta_file = os.path.join(data_dir, f'{cohort}_metadata.tsv')
        if not os.path.exists(meta_file):
            continue

        meta = pd.read_csv(meta_file, sep='\t')
        n_case = (meta['group'] == 'case').sum()
        n_ctrl = (meta['group'] == 'control').sum()

        mde = minimum_detectable_effect(n_case, n_ctrl)
        mde_fc = np.exp(mde)

        # What would we need to detect smaller effects?
        n_for_medium = required_sample_size(0.5)  # medium effect
        n_for_small = required_sample_size(0.3)   # small effect

        results.append({
            'cohort': cohort,
            'n_case': n_case,
            'n_control': n_ctrl,
            'n_total': n_case + n_ctrl,
            'mde_d': mde,
            'mde_fc': mde_fc,
            'interpretation': effect_size_interpretation(mde)
        })

        print(f"\n{cohort.upper()}")
        print(f"  Sample sizes: {n_case} cases, {n_ctrl} controls")
        print(f"  Total: {n_case + n_ctrl}")
        print(f"  Minimum detectable effect: d = {mde:.2f} ({effect_size_interpretation(mde)})")
        print(f"  Equivalent fold change: {mde_fc:.1f}x")

    print("\n" + "=" * 70)
    print("SUMMARY: Detection Capabilities")
    print("=" * 70)

    print("\n| Cohort    | N (case/ctrl) | Min Detectable | Fold Change |")
    print("|-----------|---------------|----------------|-------------|")
    for r in results:
        print(f"| {r['cohort']:9} | {r['n_case']:3}/{r['n_control']:3}        | d = {r['mde_d']:.2f}        | {r['mde_fc']:.1f}x        |")

    # What sample size would we need?
    print("\n" + "=" * 70)
    print("REQUIRED SAMPLE SIZES (balanced design)")
    print("=" * 70)

    for d, label in [(1.0, 'Large (2.7x)'), (0.8, 'Medium-Large (2.2x)'),
                      (0.5, 'Medium (1.6x)'), (0.3, 'Small-Medium (1.3x)')]:
        n = required_sample_size(d)
        print(f"  To detect {label} effect: {n} per group ({2*n} total)")

    return results


def analyze_fusobacterium_effect():
    """
    Estimate actual Fusobacterium effect size from data.
    """
    print("\n" + "=" * 70)
    print("FUSOBACTERIUM EFFECT SIZE ESTIMATION")
    print("=" * 70)

    # From our earlier analysis:
    # Zackular: Mean CRC=2465.2, Control=35.3, but Median CRC=1, Control=0
    # This suggests the effect is driven by outliers

    # Using log(mean + 1) for rough effect size:
    # log(2465 + 1) - log(35 + 1) ≈ 4.2 (huge!)
    # But this is misleading due to outliers

    # Using presence/absence:
    # CRC: 57% present, Control: 40% present
    # Odds ratio ≈ 2.0

    print("\nFrom Zackular data (n=58 CRC, n=30 control):")
    print("  Mean-based: 70x fold change (but driven by 2 outliers)")
    print("  Median-based: 1 vs 0 (minimal difference)")
    print("  Presence: 57% vs 40% (OR ≈ 2.0)")

    print("\nThe challenge:")
    print("  - The 'true' Fusobacterium effect is heterogeneous")
    print("  - A subset of CRC patients have very high loads")
    print("  - Most CRC patients have zero or low counts")
    print("  - Standard methods detect mean differences, not outlier enrichment")

    print("\nImplications:")
    print("  - Need MUCH larger samples to detect heterogeneous effects")
    print("  - Or need outlier-aware methods")
    print("  - Or need to define 'Fusobacterium-positive' subtype")


def analyze_observed_vs_detectable():
    """
    Compare observed effects to detectable thresholds.
    """
    print("\n" + "=" * 70)
    print("OBSERVED vs DETECTABLE EFFECTS")
    print("=" * 70)

    # From meta-analysis output:
    # Fusobacterium pooled effect = -0.26 (SE = 0.44)
    # This is in log-scale from hurdle model

    observed = [
        ('zeller', 0.42, 0.99, 53, 63),
        ('zackular', -1.11, 1.30, 58, 30),
        ('xiang', -0.52, 1.00, 21, 22),
        ('zhao', -0.22, 0.64, 46, 56),
    ]

    print("\n| Cohort    | Observed d | SE    | MDE d | Detectable? |")
    print("|-----------|------------|-------|-------|-------------|")

    for cohort, effect, se, n_case, n_ctrl in observed:
        mde = minimum_detectable_effect(n_case, n_ctrl)
        # Effect size in d units (roughly effect / pooled_sd)
        # For log-scale, effect IS approximately Cohen's d
        obs_d = abs(effect)
        detectable = "Yes" if obs_d > mde else "No"
        print(f"| {cohort:9} | {obs_d:10.2f} | {se:5.2f} | {mde:5.2f} | {detectable:11} |")

    print("\nConclusion:")
    print("  The observed Fusobacterium effects (|d| = 0.2 - 1.1)")
    print("  are SMALLER than the minimum detectable effects (d = 0.6 - 1.1)")
    print("  in most cohorts, explaining why it's not significant.")


if __name__ == '__main__':
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'diseases', 'crc', 'data_genus')
    data_dir = os.path.abspath(data_dir)

    analyze_cohort_power(data_dir)
    analyze_fusobacterium_effect()
    analyze_observed_vs_detectable()

    print("\n" + "=" * 70)
    print("KEY TAKEAWAY")
    print("=" * 70)
    print("""
These cohorts (n=40-100) can only detect LARGE effects (d > 0.7, ~2x fold change).

The Fusobacterium-CRC association, despite being one of the most replicated
findings in the literature, shows:
1. High heterogeneity (some patients have very high loads, most have zero)
2. Small median effect (medians are nearly equal)
3. Effect size below detection threshold for these sample sizes

To reliably detect Fusobacterium enrichment, we would need:
- n > 100 per group with current methods, OR
- Stratified analysis (Fusobacterium-high vs -low subtypes), OR
- Different statistical approach (e.g., presence/absence, outlier detection)
""")
