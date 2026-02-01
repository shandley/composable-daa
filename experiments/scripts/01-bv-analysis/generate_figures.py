#!/usr/bin/env python3
"""
Generate publication figures for BV Compositional Analysis

Figures:
1. Compositional closure demonstration (bar chart comparing sum of estimates)
2. Sensitivity analysis heatmap (load-corrected interpretations)
3. Method comparison (Hurdle vs LinDA effect sizes)
4. Top differential taxa (forest plot style)
5. Conceptual diagram: Load scenarios

Requirements:
    pip install matplotlib seaborn pandas numpy

Usage:
    python generate_figures.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# Paths
SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = SCRIPT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Colors
HEALTHY_COLOR = '#2ecc71'  # Green
BV_COLOR = '#e74c3c'       # Red
NEUTRAL_COLOR = '#95a5a6'  # Gray
HURDLE_COLOR = '#3498db'   # Blue
LINDA_COLOR = '#9b59b6'    # Purple


def load_data():
    """Load analysis results."""
    hurdle = pd.read_csv(RESULTS_DIR / "linda_results.tsv", sep='\t')
    linda = pd.read_csv(RESULTS_DIR / "linda_clr_results.tsv", sep='\t')
    sensitivity = pd.read_csv(RESULTS_DIR / "sensitivity_analysis.tsv", sep='\t')
    return hurdle, linda, sensitivity


def figure1_compositional_closure(hurdle, linda):
    """
    Figure 1: Compositional Closure Demonstration
    Bar chart comparing sum of estimates between methods.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Panel A: Distribution of estimates
    ax1 = axes[0]

    # Prepare data
    hurdle_est = hurdle['estimate'].dropna()
    linda_est = linda['estimate'].dropna()

    # Histogram
    bins = np.linspace(-5, 6, 23)
    ax1.hist(hurdle_est, bins=bins, alpha=0.6, label='Hurdle', color=HURDLE_COLOR, edgecolor='white')
    ax1.hist(linda_est, bins=bins, alpha=0.6, label='LinDA (CLR)', color=LINDA_COLOR, edgecolor='white')

    ax1.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.7)
    ax1.set_xlabel('Effect Size (log2 fold change)')
    ax1.set_ylabel('Number of Taxa')
    ax1.set_title('A. Distribution of Effect Sizes')
    ax1.legend()

    # Panel B: Sum of estimates
    ax2 = axes[1]

    hurdle_sum = hurdle_est.sum()
    linda_sum = linda_est.sum()

    bars = ax2.bar(['Hurdle\n(Count-based)', 'LinDA\n(CLR-based)'],
                   [hurdle_sum, linda_sum],
                   color=[HURDLE_COLOR, LINDA_COLOR],
                   edgecolor='black', linewidth=1.5)

    # Add value labels
    ax2.text(0, hurdle_sum + 1, f'{hurdle_sum:.2f}', ha='center', va='bottom', fontweight='bold', fontsize=12)
    ax2.text(1, 0.5, f'{linda_sum:.4f}', ha='center', va='bottom', fontweight='bold', fontsize=12)

    # Highlight the zero line
    ax2.axhline(y=0, color='red', linestyle='-', linewidth=2, alpha=0.8)
    ax2.annotate('Compositional\nClosure', xy=(1.3, 0), fontsize=10, color='red',
                 ha='left', va='center', fontweight='bold')

    ax2.set_ylabel('Sum of All Effect Sizes')
    ax2.set_title('B. Compositional Closure Test')
    ax2.set_ylim(-35, 5)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_compositional_closure.png")
    plt.savefig(FIGURES_DIR / "fig1_compositional_closure.pdf")
    plt.close()
    print("Generated: fig1_compositional_closure.png/pdf")


def figure2_sensitivity_heatmap(sensitivity):
    """
    Figure 2: Sensitivity Analysis Heatmap
    Shows how interpretation changes with load assumptions.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Prepare data for heatmap
    taxa = sensitivity['Taxon'].tolist()
    load_cols = ['Load_1x', 'Load_2x', 'Load_5x', 'Load_10x', 'Load_20x']

    data = sensitivity[load_cols].values

    # Create heatmap with diverging colormap
    vmax = np.abs(data).max()
    im = ax.imshow(data, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('log2 Fold Change (Healthy vs BV)', fontsize=10)

    # Set ticks
    ax.set_xticks(range(len(load_cols)))
    ax.set_xticklabels(['1x\n(equal)', '2x', '5x', '10x', '20x'])
    ax.set_yticks(range(len(taxa)))
    ax.set_yticklabels(taxa)

    # Add text annotations
    for i in range(len(taxa)):
        for j in range(len(load_cols)):
            value = data[i, j]
            color = 'white' if abs(value) > 3 else 'black'
            ax.text(j, i, f'{value:.1f}', ha='center', va='center',
                   fontsize=9, color=color, fontweight='bold')

    # Add vertical line at 10x (realistic scenario)
    ax.axvline(x=2.5, color='black', linestyle='--', linewidth=2, alpha=0.5)
    ax.annotate('Plausible\nBV load', xy=(3, -0.8), fontsize=9, ha='center',
               fontweight='bold', color='black')

    # Highlight the transition zone
    ax.add_patch(plt.Rectangle((2.5, 2.5), 2, 5.5, fill=False,
                               edgecolor='yellow', linewidth=3, linestyle='-'))

    ax.set_xlabel('BV Total Load (relative to Healthy)', fontsize=11)
    ax.set_title('Effect of Total Bacterial Load on Interpretation\n'
                '(Blue = Higher in Healthy, Red = Higher in BV)', fontsize=12)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_sensitivity_heatmap.png")
    plt.savefig(FIGURES_DIR / "fig2_sensitivity_heatmap.pdf")
    plt.close()
    print("Generated: fig2_sensitivity_heatmap.png/pdf")


def figure3_method_comparison(hurdle, linda):
    """
    Figure 3: Method Comparison (Hurdle vs LinDA)
    Scatter plot showing effect size correlation and key discrepancies.
    """
    # Merge on feature_id
    merged = hurdle.merge(linda, on='feature_id', suffixes=('_hurdle', '_linda'))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Scatter plot of effect sizes
    ax1 = axes[0]

    x = merged['estimate_hurdle']
    y = merged['estimate_linda']

    # Color by significance
    both_sig = (merged['q_value_hurdle'] < 0.05) & (merged['q_value_linda'] < 0.05)
    hurdle_only = (merged['q_value_hurdle'] < 0.05) & (merged['q_value_linda'] >= 0.05)
    linda_only = (merged['q_value_hurdle'] >= 0.05) & (merged['q_value_linda'] < 0.05)
    neither = ~(both_sig | hurdle_only | linda_only)

    ax1.scatter(x[neither], y[neither], c=NEUTRAL_COLOR, alpha=0.5, s=50, label='Neither sig')
    ax1.scatter(x[hurdle_only], y[hurdle_only], c=HURDLE_COLOR, alpha=0.7, s=80, label='Hurdle only')
    ax1.scatter(x[linda_only], y[linda_only], c=LINDA_COLOR, alpha=0.7, s=80, label='LinDA only')
    ax1.scatter(x[both_sig], y[both_sig], c='#27ae60', alpha=0.7, s=80, label='Both sig')

    # Add diagonal line
    lims = [-5, 5]
    ax1.plot(lims, lims, 'k--', alpha=0.5, label='1:1 line')

    # Highlight key taxa
    key_taxa = ['Lactobacillus crispatus', 'Prevotella', 'Gardnerella', 'Lactobacillus iners']
    for taxon in key_taxa:
        row = merged[merged['feature_id'] == taxon]
        if len(row) > 0:
            ax1.annotate(taxon.replace('Lactobacillus ', 'L. '),
                        xy=(row['estimate_hurdle'].values[0], row['estimate_linda'].values[0]),
                        xytext=(5, 5), textcoords='offset points', fontsize=8,
                        arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    ax1.set_xlabel('Hurdle Effect Size (log2 FC)')
    ax1.set_ylabel('LinDA Effect Size (log2 FC)')
    ax1.set_title('A. Effect Size Comparison')
    ax1.legend(loc='upper left', fontsize=8)
    ax1.set_xlim(-5, 5)
    ax1.set_ylim(-5, 5)
    ax1.set_aspect('equal')

    # Panel B: Significance comparison
    ax2 = axes[1]

    # Count categories
    categories = ['Both\nSignificant', 'Hurdle\nOnly', 'LinDA\nOnly', 'Neither']
    counts = [both_sig.sum(), hurdle_only.sum(), linda_only.sum(), neither.sum()]
    colors = ['#27ae60', HURDLE_COLOR, LINDA_COLOR, NEUTRAL_COLOR]

    bars = ax2.bar(categories, counts, color=colors, edgecolor='black', linewidth=1.5)

    # Add count labels
    for bar, count in zip(bars, counts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(count), ha='center', va='bottom', fontweight='bold', fontsize=12)

    ax2.set_ylabel('Number of Taxa')
    ax2.set_title('B. Significance Agreement (q < 0.05)')
    ax2.set_ylim(0, max(counts) * 1.15)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_method_comparison.png")
    plt.savefig(FIGURES_DIR / "fig3_method_comparison.pdf")
    plt.close()
    print("Generated: fig3_method_comparison.png/pdf")


def figure4_forest_plot(linda):
    """
    Figure 4: Forest Plot of Top Differential Taxa
    Shows effect sizes with confidence intervals.
    """
    # Select top taxa by absolute effect size
    linda_sorted = linda.copy()
    linda_sorted['abs_estimate'] = linda_sorted['estimate'].abs()
    linda_sorted = linda_sorted.sort_values('abs_estimate', ascending=True)

    # Top 20 taxa
    top_taxa = linda_sorted.tail(20).copy()

    fig, ax = plt.subplots(figsize=(10, 8))

    y_pos = range(len(top_taxa))

    # Calculate confidence intervals (estimate ± 1.96 * std_error)
    top_taxa['ci_lower'] = top_taxa['estimate'] - 1.96 * top_taxa['std_error']
    top_taxa['ci_upper'] = top_taxa['estimate'] + 1.96 * top_taxa['std_error']

    # Color by direction
    colors = [HEALTHY_COLOR if e > 0 else BV_COLOR for e in top_taxa['estimate']]

    # Plot points and error bars
    ax.errorbar(top_taxa['estimate'], y_pos,
                xerr=[top_taxa['estimate'] - top_taxa['ci_lower'],
                      top_taxa['ci_upper'] - top_taxa['estimate']],
                fmt='o', color='black', ecolor='gray', elinewidth=1.5,
                capsize=3, markersize=0)

    ax.scatter(top_taxa['estimate'], y_pos, c=colors, s=100, zorder=5, edgecolor='black')

    # Add vertical line at 0
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1)

    # Add shading for effect size magnitude
    ax.axvspan(-1, 1, alpha=0.1, color='gray', label='|log2FC| < 1')

    # Labels
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_taxa['feature_id'])
    ax.set_xlabel('log2 Fold Change (Healthy vs BV)')
    ax.set_title('Top 20 Differentially Abundant Taxa (LinDA/CLR)')

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor=HEALTHY_COLOR, edgecolor='black', label='Higher in Healthy'),
        mpatches.Patch(facecolor=BV_COLOR, edgecolor='black', label='Higher in BV'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig4_forest_plot.png")
    plt.savefig(FIGURES_DIR / "fig4_forest_plot.pdf")
    plt.close()
    print("Generated: fig4_forest_plot.png/pdf")


def figure5_conceptual_scenarios():
    """
    Figure 5: Conceptual Diagram - Load Scenarios
    Illustrates Scenario A (constant load) vs Scenario B (increased load in BV).
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    # Common parameters
    bar_width = 0.35
    categories = ['Lactobacillus', 'Gardnerella', 'Prevotella', 'Other']

    # Scenario A: Constant Load (Healthy)
    ax1 = axes[0]
    healthy_abs = [80, 5, 5, 10]  # Absolute counts (arbitrary units)
    ax1.bar(categories, healthy_abs, color=[HEALTHY_COLOR, BV_COLOR, BV_COLOR, NEUTRAL_COLOR],
            edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Absolute Abundance (×10^6)')
    ax1.set_title('Healthy Vagina\n(Total Load: 100 units)', fontsize=11, fontweight='bold')
    ax1.set_ylim(0, 100)
    ax1.tick_params(axis='x', rotation=45)

    # Scenario A: BV with constant load
    ax2 = axes[1]
    bv_constant = [10, 40, 30, 20]  # Lactobacillus died, pathogens replaced
    ax2.bar(categories, bv_constant, color=[HEALTHY_COLOR, BV_COLOR, BV_COLOR, NEUTRAL_COLOR],
            edgecolor='black', linewidth=1.5)
    ax2.set_ylabel('Absolute Abundance (×10^6)')
    ax2.set_title('BV - Scenario A\n(Constant Load: 100 units)\n"Lactobacillus dies, pathogens replace"',
                  fontsize=11, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2.tick_params(axis='x', rotation=45)

    # Add annotation
    ax2.annotate('', xy=(0, 10), xytext=(0, 80),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax2.text(0.3, 45, '↓ 8×', color='red', fontsize=12, fontweight='bold')

    # Scenario B: BV with increased load (10x)
    ax3 = axes[2]
    # If healthy has 100 units, BV has 1000 units
    # But we scale to show relative - Lactobacillus same absolute, others bloom
    bv_increased = [80, 400, 300, 220]  # Scale: Lactobacillus same, others bloomed
    # Normalize to percentage for comparison
    bv_pct = [8, 40, 30, 22]  # What the relative abundance looks like

    # Show stacked: original + bloom
    ax3.bar(categories, [80, 40, 30, 20], color=[HEALTHY_COLOR, BV_COLOR, BV_COLOR, NEUTRAL_COLOR],
            edgecolor='black', linewidth=1.5, alpha=0.5, label='Original abundance')
    ax3.bar(categories, [0, 360, 270, 200], bottom=[80, 40, 30, 20],
            color=[HEALTHY_COLOR, BV_COLOR, BV_COLOR, NEUTRAL_COLOR],
            edgecolor='black', linewidth=1.5, hatch='///', label='Biofilm bloom')

    ax3.set_ylabel('Absolute Abundance (×10^6)')
    ax3.set_title('BV - Scenario B\n(10× Load: 1000 units)\n"Lactobacillus unchanged, pathogens bloom"',
                  fontsize=11, fontweight='bold')
    ax3.set_ylim(0, 500)
    ax3.tick_params(axis='x', rotation=45)
    ax3.legend(loc='upper right', fontsize=8)

    # Add annotation showing Lactobacillus is same
    ax3.annotate('Same absolute\nabundance!', xy=(0, 80), xytext=(0.5, 150),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

    # Add bottom text
    fig.text(0.5, 0.02,
             'Both scenarios produce IDENTICAL relative abundance profiles.\n'
             'Standard 16S analysis cannot distinguish between them.',
             ha='center', fontsize=11, style='italic',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(FIGURES_DIR / "fig5_load_scenarios.png")
    plt.savefig(FIGURES_DIR / "fig5_load_scenarios.pdf")
    plt.close()
    print("Generated: fig5_load_scenarios.png/pdf")


def figure6_key_taxa_waterfall(sensitivity):
    """
    Figure 6: Waterfall chart showing load correction effect on key taxa.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    taxa = sensitivity['Taxon'].tolist()
    observed = sensitivity['Load_1x'].tolist()
    corrected_10x = sensitivity['Load_10x'].tolist()

    x = np.arange(len(taxa))
    width = 0.35

    # Bars
    bars1 = ax.bar(x - width/2, observed, width, label='Observed (CLR)',
                   color=LINDA_COLOR, edgecolor='black', linewidth=1)
    bars2 = ax.bar(x + width/2, corrected_10x, width, label='Load-corrected (10×)',
                   color='#f39c12', edgecolor='black', linewidth=1)

    # Add horizontal line at 0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)

    # Add shaded region for "no effect"
    ax.axhspan(-1, 1, alpha=0.15, color='gray')
    ax.text(7.5, 0, 'No meaningful\neffect zone', ha='right', va='center',
            fontsize=9, style='italic', color='gray')

    # Highlight taxa that flip interpretation
    for i, (obs, corr) in enumerate(zip(observed, corrected_10x)):
        if (obs < -1 and corr > -1) or (obs > 1 and corr < 1):
            ax.annotate('*', xy=(i, max(obs, corr) + 0.5), fontsize=20,
                       ha='center', color='red', fontweight='bold')

    # Labels
    ax.set_xlabel('Taxon')
    ax.set_ylabel('log2 Fold Change (Healthy vs BV)')
    ax.set_title('Effect of Load Correction on Interpretation\n'
                '(* = interpretation changes with 10× load assumption)')
    ax.set_xticks(x)
    ax.set_xticklabels(taxa, rotation=45, ha='right')
    ax.legend(loc='upper left')

    # Add arrows showing direction of change
    for i, (obs, corr) in enumerate(zip(observed, corrected_10x)):
        if abs(corr - obs) > 0.5:
            ax.annotate('', xy=(i + width/2, corr), xytext=(i - width/2, obs),
                       arrowprops=dict(arrowstyle='->', color='red', lw=1.5,
                                      connectionstyle='arc3,rad=0.3'))

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig6_load_correction_waterfall.png")
    plt.savefig(FIGURES_DIR / "fig6_load_correction_waterfall.pdf")
    plt.close()
    print("Generated: fig6_load_correction_waterfall.png/pdf")


def generate_summary_figure(hurdle, linda, sensitivity):
    """
    Generate a summary figure combining key elements.
    """
    fig = plt.figure(figsize=(14, 10))

    # Create grid
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # Panel A: Compositional closure
    ax1 = fig.add_subplot(gs[0, 0])
    hurdle_sum = hurdle['estimate'].dropna().sum()
    linda_sum = linda['estimate'].dropna().sum()
    bars = ax1.bar(['Hurdle', 'LinDA'], [hurdle_sum, linda_sum],
                   color=[HURDLE_COLOR, LINDA_COLOR], edgecolor='black')
    ax1.axhline(y=0, color='red', linestyle='-', linewidth=2)
    ax1.set_ylabel('Sum of Estimates')
    ax1.set_title('A. Compositional Closure', fontweight='bold')
    ax1.text(1, 1, f'{linda_sum:.4f}', ha='center', fontsize=10)

    # Panel B: Top taxa (simplified forest)
    ax2 = fig.add_subplot(gs[0, 1:])
    top_taxa_names = ['Lactobacillus crispatus', 'Prevotella', 'Megasphaera',
                      'Lactobacillus jensenii', 'Gardnerella']
    top_data = linda[linda['feature_id'].isin(top_taxa_names)].copy()
    top_data = top_data.set_index('feature_id').loc[top_taxa_names].reset_index()

    y_pos = range(len(top_data))
    colors = [HEALTHY_COLOR if e > 0 else BV_COLOR for e in top_data['estimate']]
    ax2.barh(y_pos, top_data['estimate'], color=colors, edgecolor='black')
    ax2.axvline(x=0, color='black', linestyle='-')
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([t.replace('Lactobacillus ', 'L. ') for t in top_data['feature_id']])
    ax2.set_xlabel('log2 FC (Healthy vs BV)')
    ax2.set_title('B. Key Differential Taxa (LinDA/CLR)', fontweight='bold')

    # Panel C: Sensitivity heatmap
    ax3 = fig.add_subplot(gs[1, :])
    taxa = sensitivity['Taxon'].tolist()
    load_cols = ['Load_1x', 'Load_2x', 'Load_5x', 'Load_10x', 'Load_20x']
    data = sensitivity[load_cols].values

    vmax = np.abs(data).max()
    im = ax3.imshow(data, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)

    ax3.set_xticks(range(len(load_cols)))
    ax3.set_xticklabels(['1×', '2×', '5×', '10×', '20×'])
    ax3.set_yticks(range(len(taxa)))
    ax3.set_yticklabels([t.replace('Lactobacillus ', 'L. ') for t in taxa])
    ax3.set_xlabel('BV Load Relative to Healthy')
    ax3.set_title('C. Load-Corrected Interpretation', fontweight='bold')

    # Add values
    for i in range(len(taxa)):
        for j in range(len(load_cols)):
            value = data[i, j]
            color = 'white' if abs(value) > 3 else 'black'
            ax3.text(j, i, f'{value:.1f}', ha='center', va='center',
                    fontsize=8, color=color)

    cbar = plt.colorbar(im, ax=ax3, shrink=0.6, pad=0.02)
    cbar.set_label('log2 FC')

    plt.savefig(FIGURES_DIR / "fig_summary.png")
    plt.savefig(FIGURES_DIR / "fig_summary.pdf")
    plt.close()
    print("Generated: fig_summary.png/pdf")


def main():
    print("Loading data...")
    hurdle, linda, sensitivity = load_data()

    print(f"Hurdle: {len(hurdle)} taxa")
    print(f"LinDA: {len(linda)} taxa")
    print(f"Sensitivity: {len(sensitivity)} taxa")

    print("\nGenerating figures...")

    figure1_compositional_closure(hurdle, linda)
    figure2_sensitivity_heatmap(sensitivity)
    figure3_method_comparison(hurdle, linda)
    figure4_forest_plot(linda)
    figure5_conceptual_scenarios()
    figure6_key_taxa_waterfall(sensitivity)
    generate_summary_figure(hurdle, linda, sensitivity)

    print(f"\nAll figures saved to: {FIGURES_DIR}")
    print("\nFigures generated:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  - {f.name}")


if __name__ == "__main__":
    main()
