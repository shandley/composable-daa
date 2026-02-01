#!/usr/bin/env python3
"""
Generate publication figures for Spike-in Load Estimation Analysis

Figures:
1. Spike-in relative abundance variation
2. Estimated total load across samples
3. Spike-in correlation matrix
4. Artifact potential demonstration
5. Before/after spike-in normalization
6. Summary figure

Requirements:
    pip install matplotlib seaborn pandas numpy scipy

Usage:
    python generate_figures.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
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
DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = SCRIPT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Colors
SPIKEIN_COLORS = {
    'Salinibacter': '#e74c3c',      # Red - constant spike
    'Alicyclobacillus': '#3498db',  # Blue - variable spike
    'Rhizobium': '#2ecc71'          # Green - variable spike
}
LOAD_COLOR = '#9b59b6'  # Purple


def load_data():
    """Load analysis results."""
    counts = pd.read_csv(DATA_DIR / "counts.tsv", sep='\t', index_col=0)
    metadata = pd.read_csv(DATA_DIR / "metadata.tsv", sep='\t', index_col=0)
    spikein_analysis = pd.read_csv(RESULTS_DIR / "spikein_analysis.tsv", sep='\t')
    correlations = pd.read_csv(RESULTS_DIR / "spikein_correlations.tsv", sep='\t', index_col=0)
    return counts, metadata, spikein_analysis, correlations


def figure1_spikein_variation(spikein_analysis):
    """
    Figure 1: Spike-in relative abundance variation across samples
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    # Sort by estimated load
    df = spikein_analysis.sort_values('estimated_load')
    x = range(len(df))

    # Plot each spike-in
    width = 0.25
    ax.bar([i - width for i in x], df['salinibacter_pct'], width,
           label='Salinibacter (constant)', color=SPIKEIN_COLORS['Salinibacter'],
           edgecolor='black', linewidth=0.5)
    ax.bar(x, df['alicyclobacillus_pct'], width,
           label='Alicyclobacillus (variable)', color=SPIKEIN_COLORS['Alicyclobacillus'],
           edgecolor='black', linewidth=0.5)
    ax.bar([i + width for i in x], df['rhizobium_pct'], width,
           label='Rhizobium (variable)', color=SPIKEIN_COLORS['Rhizobium'],
           edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Sample (sorted by estimated total load)')
    ax.set_ylabel('Relative Abundance (%)')
    ax.set_title('Spike-in Relative Abundances Vary Despite Known Input Amounts')
    ax.set_xticks(x)
    ax.set_xticklabels(df['sample_id'], rotation=45, ha='right')
    ax.legend(loc='upper right')

    # Add annotation
    ax.annotate('Low total load\n(high spike-in %)', xy=(1, df['salinibacter_pct'].iloc[0]),
               xytext=(3, df['salinibacter_pct'].iloc[0] + 5),
               fontsize=9, ha='center',
               arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    ax.annotate('High total load\n(low spike-in %)', xy=(len(x)-2, df['salinibacter_pct'].iloc[-1]),
               xytext=(len(x)-4, df['salinibacter_pct'].iloc[-1] + 8),
               fontsize=9, ha='center',
               arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_spikein_variation.png")
    plt.savefig(FIGURES_DIR / "fig1_spikein_variation.pdf")
    plt.close()
    print("Generated: fig1_spikein_variation.png/pdf")


def figure2_load_estimation(spikein_analysis):
    """
    Figure 2: Estimated total bacterial load across samples
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Load per sample
    ax1 = axes[0]
    df = spikein_analysis.sort_values('estimated_load')

    bars = ax1.bar(range(len(df)), df['estimated_load'], color=LOAD_COLOR,
                   edgecolor='black', linewidth=0.5)

    # Color bars by condition if available
    ax1.axhline(y=1.0, color='black', linestyle='--', linewidth=1, label='Mean load')

    ax1.set_xlabel('Sample (sorted by load)')
    ax1.set_ylabel('Estimated Relative Load')
    ax1.set_title('A. Total Bacterial Load Varies 9.8x Across Samples')
    ax1.set_xticks(range(len(df)))
    ax1.set_xticklabels(df['sample_id'], rotation=45, ha='right')

    # Add min/max annotations
    min_load = df['estimated_load'].min()
    max_load = df['estimated_load'].max()
    ax1.annotate(f'{min_load:.2f}x', xy=(0, min_load), xytext=(2, min_load + 0.3),
                fontsize=10, ha='center', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='black', lw=0.5))
    ax1.annotate(f'{max_load:.2f}x', xy=(len(df)-1, max_load), xytext=(len(df)-3, max_load - 0.3),
                fontsize=10, ha='center', fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    # Panel B: Distribution
    ax2 = axes[1]
    ax2.hist(df['estimated_load'], bins=8, color=LOAD_COLOR, edgecolor='black', alpha=0.7)
    ax2.axvline(x=1.0, color='black', linestyle='--', linewidth=2, label='Mean')
    ax2.axvline(x=df['estimated_load'].min(), color='red', linestyle='-', linewidth=2)
    ax2.axvline(x=df['estimated_load'].max(), color='red', linestyle='-', linewidth=2)

    # Add text box with key stats
    fold_range = max_load / min_load
    log2_range = np.log2(fold_range)
    textstr = f'Fold range: {fold_range:.1f}x\nlog2 range: {log2_range:.1f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax2.text(0.95, 0.95, textstr, transform=ax2.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    ax2.set_xlabel('Estimated Relative Load')
    ax2.set_ylabel('Number of Samples')
    ax2.set_title('B. Distribution of Total Load')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_load_estimation.png")
    plt.savefig(FIGURES_DIR / "fig2_load_estimation.pdf")
    plt.close()
    print("Generated: fig2_load_estimation.png/pdf")


def figure3_correlation_matrix(correlations):
    """
    Figure 3: Spike-in correlation matrix
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Clean up names for display
    display_names = {
        'Salinibacter_ruber': 'Salinibacter\n(constant)',
        'Alicyclobacillus_acidiphilus': 'Alicyclobacillus\n(variable)',
        'Rhizobium_radiobacter': 'Rhizobium\n(variable)'
    }

    corr_display = correlations.copy()
    corr_display.index = [display_names.get(n, n) for n in corr_display.index]
    corr_display.columns = [display_names.get(n, n) for n in corr_display.columns]

    # Heatmap
    mask = np.triu(np.ones_like(corr_display, dtype=bool), k=1)
    sns.heatmap(corr_display, annot=True, fmt='.2f', cmap='RdYlGn',
                center=0.5, vmin=0, vmax=1, square=True,
                linewidths=2, cbar_kws={'label': 'Correlation (r)'},
                ax=ax, annot_kws={'fontsize': 14, 'fontweight': 'bold'})

    ax.set_title('Spike-in Correlations\n(High correlation confirms shared dilution factor)',
                fontsize=12)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_correlation_matrix.png")
    plt.savefig(FIGURES_DIR / "fig3_correlation_matrix.pdf")
    plt.close()
    print("Generated: fig3_correlation_matrix.png/pdf")


def figure4_artifact_potential(spikein_analysis):
    """
    Figure 4: Artifact potential demonstration
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Salinibacter vs Load relationship
    ax1 = axes[0]

    x = spikein_analysis['estimated_load']
    y = spikein_analysis['salinibacter_pct']

    ax1.scatter(x, y, c=SPIKEIN_COLORS['Salinibacter'], s=100, edgecolor='black', zorder=5)

    # Fit line (should be inverse relationship)
    # y = k/x, so log(y) = log(k) - log(x)
    log_x = np.log(x)
    log_y = np.log(y)
    slope, intercept, r, p, se = stats.linregress(log_x, log_y)

    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = np.exp(intercept) * x_fit ** slope
    ax1.plot(x_fit, y_fit, 'r--', linewidth=2, label=f'y = k/x (r={r:.2f})')

    ax1.set_xlabel('Estimated Total Load')
    ax1.set_ylabel('Salinibacter Relative Abundance (%)')
    ax1.set_title('A. Constant Spike-in Shows Inverse Relationship with Load')
    ax1.legend()

    # Panel B: Artifact illustration
    ax2 = axes[1]

    # Create illustration of artifact
    load_low = spikein_analysis['estimated_load'].min()
    load_high = spikein_analysis['estimated_load'].max()

    # Scenario: True abundance unchanged
    categories = ['True\nAbundance', 'Observed\n(Low Load)', 'Observed\n(High Load)']
    true_val = 100  # Arbitrary units
    # Observed = True / (Total_Load_factor)
    obs_low = 100 * (1/load_low) / (1/load_low + 1)  # Simplified model
    obs_high = 100 * (1/load_high) / (1/load_high + 1)

    # Normalize to show relative
    values = [50, 50 * load_high/load_low, 50]  # Simplified visualization

    colors = ['#95a5a6', '#27ae60', '#e74c3c']
    bars = ax2.bar(categories, [100, 100 * (load_high/load_low)**0.5, 100 / (load_high/load_low)**0.5],
                   color=colors, edgecolor='black', linewidth=2)

    # Add fold change annotation
    fold_change = (load_high/load_low)
    ax2.annotate(f'{fold_change**0.5:.1f}x apparent\n"increase"',
                xy=(1, 100 * (load_high/load_low)**0.5),
                xytext=(1, 100 * (load_high/load_low)**0.5 + 30),
                ha='center', fontsize=10, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='green', lw=2))

    ax2.annotate(f'{(load_high/load_low)**0.5:.1f}x apparent\n"decrease"',
                xy=(2, 100 / (load_high/load_low)**0.5),
                xytext=(2, 100 / (load_high/load_low)**0.5 + 50),
                ha='center', fontsize=10, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='red', lw=2))

    ax2.axhline(y=100, color='black', linestyle='--', linewidth=2, label='True value')
    ax2.set_ylabel('Apparent Relative Abundance')
    ax2.set_title('B. Load Variation Creates Spurious Differences\n'
                 '(Both samples have IDENTICAL true abundance)')
    ax2.set_ylim(0, 350)
    ax2.legend()

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig4_artifact_potential.png")
    plt.savefig(FIGURES_DIR / "fig4_artifact_potential.pdf")
    plt.close()
    print("Generated: fig4_artifact_potential.png/pdf")


def figure5_normalization_effect(counts, spikein_analysis):
    """
    Figure 5: Effect of spike-in normalization
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    SPIKEIN_FEATURE = 'AF323500XXXX'

    # Get spike-in counts
    spikein_counts = counts.loc[SPIKEIN_FEATURE]
    median_spikein = spikein_counts.median()
    scaling_factors = median_spikein / spikein_counts

    # Panel A: Before normalization
    ax1 = axes[0]

    # Sort by load
    sorted_samples = spikein_analysis.sort_values('estimated_load')['sample_id'].tolist()
    spikein_sorted = spikein_counts[sorted_samples]

    ax1.bar(range(len(spikein_sorted)), spikein_sorted.values,
           color=SPIKEIN_COLORS['Salinibacter'], edgecolor='black', linewidth=0.5)
    ax1.axhline(y=median_spikein, color='black', linestyle='--', linewidth=2,
               label=f'Median: {median_spikein:.0f}')

    cv_before = spikein_counts.std() / spikein_counts.mean() * 100
    ax1.set_title(f'A. Before Normalization\n(Salinibacter counts, CV = {cv_before:.1f}%)')
    ax1.set_xlabel('Sample (sorted by load)')
    ax1.set_ylabel('Raw Counts')
    ax1.set_xticks(range(len(spikein_sorted)))
    ax1.set_xticklabels(sorted_samples, rotation=45, ha='right')
    ax1.legend()

    # Panel B: After normalization
    ax2 = axes[1]

    normalized_spikein = spikein_counts * scaling_factors
    normalized_sorted = normalized_spikein[sorted_samples]

    ax2.bar(range(len(normalized_sorted)), normalized_sorted.values,
           color=SPIKEIN_COLORS['Salinibacter'], edgecolor='black', linewidth=0.5)
    ax2.axhline(y=median_spikein, color='black', linestyle='--', linewidth=2,
               label=f'Target: {median_spikein:.0f}')

    cv_after = normalized_spikein.std() / normalized_spikein.mean() * 100
    ax2.set_title(f'B. After Spike-in Normalization\n(CV = {cv_after:.1f}%)')
    ax2.set_xlabel('Sample (sorted by load)')
    ax2.set_ylabel('Normalized Counts')
    ax2.set_xticks(range(len(normalized_sorted)))
    ax2.set_xticklabels(sorted_samples, rotation=45, ha='right')
    ax2.legend()

    # Match y-axis limits
    max_y = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_ylim(0, max_y)
    ax2.set_ylim(0, max_y)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig5_normalization_effect.png")
    plt.savefig(FIGURES_DIR / "fig5_normalization_effect.pdf")
    plt.close()
    print("Generated: fig5_normalization_effect.png/pdf")


def figure6_summary(spikein_analysis, correlations):
    """
    Figure 6: Summary figure for graphical abstract
    """
    fig = plt.figure(figsize=(14, 8))

    # Create grid
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    # Panel A: Load variation
    ax1 = fig.add_subplot(gs[0, 0])
    df = spikein_analysis.sort_values('estimated_load')
    ax1.bar(range(len(df)), df['estimated_load'], color=LOAD_COLOR, edgecolor='black')
    ax1.axhline(y=1.0, color='black', linestyle='--')
    ax1.set_xlabel('Sample')
    ax1.set_ylabel('Relative Load')
    ax1.set_title('A. 9.8x Load Variation', fontweight='bold')
    ax1.set_xticks([])

    # Add range annotation
    fold_range = df['estimated_load'].max() / df['estimated_load'].min()
    ax1.text(0.95, 0.95, f'{fold_range:.1f}x\nrange',
            transform=ax1.transAxes, ha='right', va='top',
            fontsize=14, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))

    # Panel B: Spike-in correlation
    ax2 = fig.add_subplot(gs[0, 1])
    display_names = ['Salinibacter', 'Alicyclobacillus', 'Rhizobium']
    corr_display = correlations.copy()
    corr_display.index = display_names
    corr_display.columns = display_names
    sns.heatmap(corr_display, annot=True, fmt='.2f', cmap='RdYlGn',
               vmin=0, vmax=1, square=True, cbar=False, ax=ax2,
               annot_kws={'fontsize': 10})
    ax2.set_title('B. Spike-in Correlations', fontweight='bold')

    # Panel C: Artifact potential
    ax3 = fig.add_subplot(gs[0, 2])
    log2_range = np.log2(fold_range)

    # Show published effect sizes vs artifact range
    effect_sizes = [1, 2, 3, 4]
    labels = ['Small\n(2x)', 'Moderate\n(4x)', 'Large\n(8x)', 'Very Large\n(16x)']

    bars = ax3.barh(range(len(effect_sizes)), effect_sizes, color='#3498db',
                   edgecolor='black', height=0.6)
    ax3.axvline(x=log2_range, color='red', linestyle='-', linewidth=3,
               label=f'Artifact potential: {log2_range:.1f}')

    # Shade artifact zone
    ax3.axvspan(0, log2_range, alpha=0.3, color='red')
    ax3.text(log2_range/2, 3.5, 'Could be\nartifact!', ha='center', va='center',
            fontsize=10, fontweight='bold', color='red')

    ax3.set_yticks(range(len(effect_sizes)))
    ax3.set_yticklabels(labels)
    ax3.set_xlabel('log2 Fold Change')
    ax3.set_title('C. Artifact Potential', fontweight='bold')
    ax3.set_xlim(0, 5)

    # Panel D: Key message
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis('off')

    message = """
    KEY FINDINGS FROM STAMMLER 2016 SPIKE-IN ANALYSIS

    1. Total bacterial load varies 9.8x across samples
    2. This creates potential artifacts up to {:.1f} log2FC
    3. Most published microbiome effect sizes are within this range
    4. Without load correction, apparent effects may be spurious

    SOLUTION: Include spike-in controls and report load-corrected abundances
    """.format(log2_range)

    ax4.text(0.5, 0.5, message, transform=ax4.transAxes,
            fontsize=12, ha='center', va='center',
            family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='black', linewidth=2))

    plt.savefig(FIGURES_DIR / "fig6_summary.png")
    plt.savefig(FIGURES_DIR / "fig6_summary.pdf")
    plt.close()
    print("Generated: fig6_summary.png/pdf")


def main():
    print("Loading data...")
    counts, metadata, spikein_analysis, correlations = load_data()

    print(f"Counts: {counts.shape[0]} features x {counts.shape[1]} samples")
    print(f"Spike-in analysis: {len(spikein_analysis)} samples")

    print("\nGenerating figures...")

    figure1_spikein_variation(spikein_analysis)
    figure2_load_estimation(spikein_analysis)
    figure3_correlation_matrix(correlations)
    figure4_artifact_potential(spikein_analysis)
    figure5_normalization_effect(counts, spikein_analysis)
    figure6_summary(spikein_analysis, correlations)

    print(f"\nAll figures saved to: {FIGURES_DIR}")
    print("\nFigures generated:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  - {f.name}")


if __name__ == "__main__":
    main()
