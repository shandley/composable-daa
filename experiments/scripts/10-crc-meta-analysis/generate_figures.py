#!/usr/bin/env python3
"""
Experiment 10: CRC Meta-Analysis Figures

Generates publication-quality figures demonstrating cross-study
consistency analysis for CRC microbiome associations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set publication style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 13
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'

# Paths
SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
FIGURES_DIR = SCRIPT_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def load_data():
    """Load all result files."""
    data = {}

    # Cohort summary
    summary_file = RESULTS_DIR / "cohort_summary.tsv"
    if summary_file.exists():
        data['cohort_summary'] = pd.read_csv(summary_file, sep='\t')

    # Consistency summary
    consistency_file = RESULTS_DIR / "consistency_summary.tsv"
    if consistency_file.exists():
        data['consistency'] = pd.read_csv(consistency_file, sep='\t')

    # Consistency matrix
    matrix_file = RESULTS_DIR / "consistency_matrix.tsv"
    if matrix_file.exists():
        data['matrix'] = pd.read_csv(matrix_file, sep='\t')

    return data


def fig1_cohort_comparison(data, save=True):
    """
    Figure 1: Per-Cohort Detection Rates
    """
    if 'cohort_summary' not in data:
        print("Skipping Figure 1: cohort summary not found")
        return None

    df = data['cohort_summary'].copy()

    # Focus on hurdle method
    df_hurdle = df[df['method'] == 'hurdle'].copy()

    if len(df_hurdle) == 0:
        print("Skipping Figure 1: no hurdle results")
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    cohorts = df_hurdle['cohort'].values
    x = np.arange(len(cohorts))
    width = 0.25

    n_sig_005 = df_hurdle['n_sig_005'].values
    n_sig_010 = df_hurdle['n_sig_010'].values
    n_robust = df_hurdle['n_robust'].values

    bars1 = ax.bar(x - width, n_sig_005, width, label='q < 0.05', color='#3498db', alpha=0.8)
    bars2 = ax.bar(x, n_sig_010, width, label='q < 0.10', color='#2ecc71', alpha=0.8)
    bars3 = ax.bar(x + width, n_robust, width, label='Robust (>3.3 log2FC)', color='#e74c3c', alpha=0.8)

    ax.set_xlabel('Cohort')
    ax.set_ylabel('Number of Significant Taxa')
    ax.set_title('Detection Rates by Cohort (Hurdle Model)', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(cohorts)
    ax.legend()

    # Add cohort characteristics
    characteristics = [
        'Large effect\nLow sparsity',
        'Large effect\nHigh sparsity',
        'Moderate effect\nLow sparsity',
        'Moderate effect\nHigh sparsity',
    ]
    for i, (xi, char) in enumerate(zip(x, characteristics)):
        ax.text(xi, -3, char, ha='center', fontsize=8, style='italic')

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig1_cohort_comparison.png")
        fig.savefig(FIGURES_DIR / "fig1_cohort_comparison.pdf")
        print("Saved: fig1_cohort_comparison.png/pdf")

    return fig


def fig2_consistency_distribution(data, save=True):
    """
    Figure 2: Cross-Cohort Consistency Distribution
    """
    if 'consistency' not in data:
        print("Skipping Figure 2: consistency data not found")
        return None

    df = data['consistency'].copy()

    fig, ax = plt.subplots(figsize=(10, 6))

    n_cohorts = df['n_cohorts'].values
    n_features = df['n_features'].values
    categories = df['category'].values

    colors = ['#95a5a6', '#f39c12', '#e67e22', '#e74c3c', '#2ecc71']
    bars = ax.bar(n_cohorts, n_features, color=colors, alpha=0.8, edgecolor='black')

    ax.set_xlabel('Number of Cohorts with Significant Detection')
    ax.set_ylabel('Number of Taxa')
    ax.set_title('Cross-Cohort Consistency of Findings', fontsize=14, fontweight='bold')
    ax.set_xticks(n_cohorts)

    # Add category labels
    for bar, cat in zip(bars, categories):
        height = bar.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width()/2, height),
                       ha='center', va='bottom', fontsize=10, fontweight='bold')

    # Add legend for interpretation
    ax.axhline(y=0, color='black', linestyle='-')
    ax.text(4, max(n_features) * 0.9, 'Most robust\n(4 cohorts)',
           ha='center', fontsize=9, color='#2ecc71', fontweight='bold')
    ax.text(0, max(n_features) * 0.9, 'Not significant\n(0 cohorts)',
           ha='center', fontsize=9, color='#95a5a6')

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig2_consistency_distribution.png")
        fig.savefig(FIGURES_DIR / "fig2_consistency_distribution.pdf")
        print("Saved: fig2_consistency_distribution.png/pdf")

    return fig


def fig3_consistency_heatmap(data, save=True):
    """
    Figure 3: Consistency Heatmap (Taxa x Cohorts)
    """
    if 'matrix' not in data:
        print("Skipping Figure 3: consistency matrix not found")
        return None

    df = data['matrix'].copy()

    # Filter to taxa significant in at least 1 cohort
    df['n_sig'] = df['n_cohorts_sig'].fillna(0).astype(int)
    df_sig = df[df['n_sig'] >= 1].head(30)  # Top 30 for visibility

    if len(df_sig) == 0:
        print("Skipping Figure 3: no significant taxa")
        return None

    fig, ax = plt.subplots(figsize=(10, 8))

    # Create binary matrix
    sig_cols = ['cohort1_sig', 'cohort2_sig', 'cohort3_sig', 'cohort4_sig']
    matrix = pd.DataFrame()
    for col in sig_cols:
        if col in df_sig.columns:
            matrix[col.replace('_sig', '')] = (df_sig[col] == 'SIG').astype(int)

    if matrix.empty:
        print("Skipping Figure 3: no data for heatmap")
        return None

    matrix.index = df_sig['feature'].values

    # Create heatmap
    cmap = sns.color_palette(['#ecf0f1', '#2ecc71'])
    sns.heatmap(matrix, cmap=cmap, cbar=False, linewidths=0.5, ax=ax)

    ax.set_xlabel('Cohort')
    ax.set_ylabel('Taxa')
    ax.set_title('Significance Pattern Across Cohorts', fontsize=14, fontweight='bold')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#ecf0f1', edgecolor='black', label='Not significant'),
                      Patch(facecolor='#2ecc71', edgecolor='black', label='Significant')]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1))

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig3_consistency_heatmap.png")
        fig.savefig(FIGURES_DIR / "fig3_consistency_heatmap.pdf")
        print("Saved: fig3_consistency_heatmap.png/pdf")

    return fig


def fig4_meta_analysis_framework(save=True):
    """
    Figure 4: Meta-Analysis Best Practices Diagram
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(6, 9.5, 'Meta-Analysis Framework for Microbiome Studies', fontsize=16, fontweight='bold', ha='center')

    # Step 1: Consistent methods
    ax.add_patch(plt.Rectangle((0.5, 7), 3.5, 2, fill=True,
                               facecolor='#3498db', alpha=0.3, edgecolor='#3498db', linewidth=2))
    ax.text(2.25, 8.5, '1. Consistent Methods', fontsize=11, fontweight='bold', ha='center')
    ax.text(2.25, 7.8, 'Same pipeline for', fontsize=9, ha='center')
    ax.text(2.25, 7.3, 'all cohorts', fontsize=9, ha='center')

    # Step 2: Per-cohort results
    ax.add_patch(plt.Rectangle((4.25, 7), 3.5, 2, fill=True,
                               facecolor='#2ecc71', alpha=0.3, edgecolor='#2ecc71', linewidth=2))
    ax.text(6, 8.5, '2. Per-Cohort Results', fontsize=11, fontweight='bold', ha='center')
    ax.text(6, 7.8, 'Effect sizes and', fontsize=9, ha='center')
    ax.text(6, 7.3, 'significance', fontsize=9, ha='center')

    # Step 3: Consistency assessment
    ax.add_patch(plt.Rectangle((8, 7), 3.5, 2, fill=True,
                               facecolor='#e74c3c', alpha=0.3, edgecolor='#e74c3c', linewidth=2))
    ax.text(9.75, 8.5, '3. Consistency Check', fontsize=11, fontweight='bold', ha='center')
    ax.text(9.75, 7.8, 'How many cohorts', fontsize=9, ha='center')
    ax.text(9.75, 7.3, 'agree?', fontsize=9, ha='center')

    # Arrows
    ax.annotate('', xy=(4.1, 8), xytext=(4, 8),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.annotate('', xy=(7.9, 8), xytext=(7.75, 8),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))

    # Outcomes
    ax.add_patch(plt.Rectangle((1, 4), 4.5, 2.5, fill=True,
                               facecolor='#2ecc71', alpha=0.2, edgecolor='#2ecc71', linewidth=2))
    ax.text(3.25, 6, 'ROBUST (3-4 cohorts)', fontsize=11, fontweight='bold', ha='center', color='#27ae60')
    ax.text(3.25, 5.3, '• Significant in most cohorts', fontsize=9, ha='center')
    ax.text(3.25, 4.8, '• Effect size > 3.3 log2FC', fontsize=9, ha='center')
    ax.text(3.25, 4.3, '• Consistent direction', fontsize=9, ha='center')

    ax.add_patch(plt.Rectangle((6.5, 4), 4.5, 2.5, fill=True,
                               facecolor='#e74c3c', alpha=0.2, edgecolor='#e74c3c', linewidth=2))
    ax.text(8.75, 6, 'AT-RISK (0-2 cohorts)', fontsize=11, fontweight='bold', ha='center', color='#c0392b')
    ax.text(8.75, 5.3, '• Significant in few cohorts', fontsize=9, ha='center')
    ax.text(8.75, 4.8, '• Small effect sizes', fontsize=9, ha='center')
    ax.text(8.75, 4.3, '• May be cohort-specific', fontsize=9, ha='center')

    # Arrows down
    ax.annotate('', xy=(3.25, 6.8), xytext=(6, 6.8),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.annotate('', xy=(8.75, 6.8), xytext=(6, 6.8),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))

    # Key message
    ax.add_patch(plt.Rectangle((2, 1), 8, 2, fill=True,
                facecolor='#f1c40f', alpha=0.3, edgecolor='#f39c12', linewidth=2))
    ax.text(6, 2.5, 'Cross-Study Consistency is the Gold Standard', fontsize=12, fontweight='bold', ha='center')
    ax.text(6, 1.7, 'True biological signals should replicate across independent cohorts', fontsize=10, ha='center', style='italic')
    ax.text(6, 1.2, 'Cohort-specific findings require independent validation', fontsize=10, ha='center', style='italic')

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig4_meta_framework.png")
        fig.savefig(FIGURES_DIR / "fig4_meta_framework.pdf")
        print("Saved: fig4_meta_framework.png/pdf")

    return fig


def fig5_forest_plot(data, save=True):
    """
    Figure 5: Forest Plot of Effect Sizes by Cohort
    """
    # Create simulated forest plot data
    fig, ax = plt.subplots(figsize=(10, 8))

    # Example taxa with varying consistency
    taxa = [
        'Taxon A (4 cohorts)',
        'Taxon B (3 cohorts)',
        'Taxon C (2 cohorts)',
        'Taxon D (1 cohort)',
        'Taxon E (1 cohort)',
    ]

    # Simulated effect sizes (mean, CI) for each cohort
    # Format: [cohort1, cohort2, cohort3, cohort4]
    effects = [
        [(3.5, 0.5), (4.0, 0.6), (3.2, 0.4), (3.8, 0.5)],  # Taxon A - consistent
        [(2.5, 0.4), (3.0, 0.5), (2.8, 0.6), (0, 1)],      # Taxon B - 3 cohorts
        [(2.0, 0.5), (0, 1), (1.8, 0.4), (0, 1)],           # Taxon C - 2 cohorts
        [(4.0, 0.3), (0, 1), (0, 1), (0, 1)],               # Taxon D - 1 cohort
        [(0, 1), (0, 1), (0, 1), (2.5, 0.4)],               # Taxon E - 1 cohort
    ]

    colors = ['#3498db', '#2ecc71', '#f39c12', '#e74c3c']
    cohort_labels = ['Cohort 1', 'Cohort 2', 'Cohort 3', 'Cohort 4']

    y_positions = np.arange(len(taxa))

    for i, taxon in enumerate(taxa):
        for j, (mean, se) in enumerate(effects[i]):
            if mean != 0:  # Only plot if significant
                offset = (j - 1.5) * 0.15
                ax.errorbar(mean, i + offset, xerr=1.96*se,
                           fmt='o', color=colors[j], capsize=3, markersize=6,
                           label=cohort_labels[j] if i == 0 else None)

    # Artifact threshold
    ax.axvline(x=3.3, color='red', linestyle='--', linewidth=2, label='Artifact threshold')
    ax.axvline(x=-3.3, color='red', linestyle='--', linewidth=2)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=1)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(taxa)
    ax.set_xlabel('Effect Size (log2FC)')
    ax.set_title('Forest Plot: Effect Sizes by Cohort', fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    ax.set_xlim(-1, 6)

    # Add shading for robust zone
    ax.axvspan(3.3, 6, alpha=0.1, color='green', label='Robust zone')

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig5_forest_plot.png")
        fig.savefig(FIGURES_DIR / "fig5_forest_plot.pdf")
        print("Saved: fig5_forest_plot.png/pdf")

    return fig


def main():
    """Generate all figures for Experiment 10."""
    print("=" * 60)
    print("Experiment 10: CRC Meta-Analysis Figures")
    print("=" * 60)

    # Load data
    data = load_data()
    print(f"\nLoaded data: {list(data.keys())}")

    # Generate figures
    print("\nGenerating figures...")

    fig1_cohort_comparison(data)
    fig2_consistency_distribution(data)
    fig3_consistency_heatmap(data)
    fig4_meta_analysis_framework()
    fig5_forest_plot(data)

    print("\n" + "=" * 60)
    print("Figure generation complete!")
    print(f"Output directory: {FIGURES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
