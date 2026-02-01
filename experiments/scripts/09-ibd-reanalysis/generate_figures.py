#!/usr/bin/env python3
"""
Experiment 09: IBD Cohort Reanalysis Figures

Generates publication-quality figures demonstrating artifact risk
assessment for IBD microbiome findings.
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

# Artifact threshold
ARTIFACT_THRESHOLD = 3.3


def load_data():
    """Load all result files."""
    data = {}

    # Artifact risk assessment
    risk_file = RESULTS_DIR / "artifact_risk.tsv"
    if risk_file.exists():
        data['risk'] = pd.read_csv(risk_file, sep='\t')

    # FPR validation
    fpr_file = RESULTS_DIR / "fpr_validation.tsv"
    if fpr_file.exists():
        data['fpr'] = pd.read_csv(fpr_file, sep='\t')

    # Load individual method results for effect size distribution
    for dataset in ['ibd_moderate', 'ibd_large']:
        for method in ['linda', 'hurdle']:
            filepath = RESULTS_DIR / f"{dataset}_{method}.tsv"
            if filepath.exists():
                data[f'{dataset}_{method}'] = pd.read_csv(filepath, sep='\t')

    return data


def fig1_artifact_risk_comparison(data, save=True):
    """
    Figure 1: Artifact Risk by Effect Size and Method
    """
    if 'risk' not in data:
        print("Skipping Figure 1: artifact risk data not found")
        return None

    df = data['risk'].copy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Moderate effect dataset
    ax = axes[0]
    df_mod = df[df['dataset'] == 'moderate_effect'].copy()
    if len(df_mod) > 0:
        methods = df_mod['method'].values
        x = np.arange(len(methods))
        width = 0.35

        robust = df_mod['n_robust'].values
        at_risk = df_mod['n_at_risk'].values

        ax.bar(x - width/2, robust, width, label='Robust (>3.3 log2FC)', color='#2ecc71', alpha=0.8)
        ax.bar(x + width/2, at_risk, width, label='At-Risk (<3.3 log2FC)', color='#e74c3c', alpha=0.8)

        ax.set_xlabel('Method')
        ax.set_ylabel('Number of Significant Taxa')
        ax.set_title('A. Moderate Effect (2.0 log2FC)')
        ax.set_xticks(x)
        ax.set_xticklabels([m.upper() for m in methods], rotation=45, ha='right')
        ax.legend()

    # Panel B: Large effect dataset
    ax = axes[1]
    df_large = df[df['dataset'] == 'large_effect'].copy()
    if len(df_large) > 0:
        methods = df_large['method'].values
        x = np.arange(len(methods))

        robust = df_large['n_robust'].values
        at_risk = df_large['n_at_risk'].values

        ax.bar(x - width/2, robust, width, label='Robust (>3.3 log2FC)', color='#2ecc71', alpha=0.8)
        ax.bar(x + width/2, at_risk, width, label='At-Risk (<3.3 log2FC)', color='#e74c3c', alpha=0.8)

        ax.set_xlabel('Method')
        ax.set_ylabel('Number of Significant Taxa')
        ax.set_title('B. Large Effect (4.0 log2FC)')
        ax.set_xticks(x)
        ax.set_xticklabels([m.upper() for m in methods], rotation=45, ha='right')
        ax.legend()

    plt.suptitle('Artifact Risk Assessment: IBD Data', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig1_artifact_risk.png")
        fig.savefig(FIGURES_DIR / "fig1_artifact_risk.pdf")
        print("Saved: fig1_artifact_risk.png/pdf")

    return fig


def fig2_effect_size_distribution(data, save=True):
    """
    Figure 2: Effect Size Distribution with Artifact Threshold
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Moderate effect
    ax = axes[0]
    if 'ibd_moderate_hurdle' in data:
        df = data['ibd_moderate_hurdle'].copy()
        effects = df['estimate'].dropna().values

        ax.hist(effects, bins=30, color='#3498db', alpha=0.7, edgecolor='black')
        ax.axvline(x=ARTIFACT_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Artifact threshold ({ARTIFACT_THRESHOLD})')
        ax.axvline(x=-ARTIFACT_THRESHOLD, color='red', linestyle='--', linewidth=2)

        ax.set_xlabel('Effect Size (log2FC)')
        ax.set_ylabel('Number of Taxa')
        ax.set_title('A. Moderate Effect Dataset')
        ax.legend()

    # Panel B: Large effect
    ax = axes[1]
    if 'ibd_large_hurdle' in data:
        df = data['ibd_large_hurdle'].copy()
        effects = df['estimate'].dropna().values

        ax.hist(effects, bins=30, color='#e74c3c', alpha=0.7, edgecolor='black')
        ax.axvline(x=ARTIFACT_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Artifact threshold ({ARTIFACT_THRESHOLD})')
        ax.axvline(x=-ARTIFACT_THRESHOLD, color='red', linestyle='--', linewidth=2)

        ax.set_xlabel('Effect Size (log2FC)')
        ax.set_ylabel('Number of Taxa')
        ax.set_title('B. Large Effect Dataset')
        ax.legend()

    plt.suptitle('Effect Size Distribution with Artifact Threshold', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig2_effect_distribution.png")
        fig.savefig(FIGURES_DIR / "fig2_effect_distribution.pdf")
        print("Saved: fig2_effect_distribution.png/pdf")

    return fig


def fig3_fpr_validation(data, save=True):
    """
    Figure 3: FPR Calibration on Null IBD Data
    """
    if 'fpr' not in data:
        print("Skipping Figure 3: FPR data not found")
        return None

    df = data['fpr'].copy()

    fig, ax = plt.subplots(figsize=(8, 6))

    methods = df['method'].values
    x = np.arange(len(methods))
    width = 0.35

    # Convert to numeric
    fpr_005 = pd.to_numeric(df['fpr_005'], errors='coerce').fillna(0).values
    fpr_010 = pd.to_numeric(df['fpr_010'], errors='coerce').fillna(0).values

    bars1 = ax.bar(x - width/2, fpr_005, width, label='alpha = 0.05', color='#e74c3c', alpha=0.8)
    bars2 = ax.bar(x + width/2, fpr_010, width, label='alpha = 0.10', color='#f39c12', alpha=0.8)

    # Nominal FPR lines
    ax.axhline(y=0.05, color='#e74c3c', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhline(y=0.10, color='#f39c12', linestyle='--', linewidth=1.5, alpha=0.7)

    ax.set_xlabel('Method')
    ax.set_ylabel('False Positive Rate')
    ax.set_title('FPR Calibration on Null IBD Data', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([m.upper() for m in methods], rotation=45, ha='right')
    ax.legend()
    ax.set_ylim(0, 0.20)

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig3_fpr_validation.png")
        fig.savefig(FIGURES_DIR / "fig3_fpr_validation.pdf")
        print("Saved: fig3_fpr_validation.png/pdf")

    return fig


def fig4_clinical_implications(save=True):
    """
    Figure 4: Clinical Implications Diagram
    """
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'Artifact Risk in IBD Microbiome Studies', fontsize=16, fontweight='bold', ha='center')

    # The problem
    ax.add_patch(plt.Rectangle((0.5, 6.5), 4, 2.5, fill=True,
                               facecolor='#e74c3c', alpha=0.2, edgecolor='#e74c3c', linewidth=2))
    ax.text(2.5, 8.5, 'IBD-Specific Confounders', fontsize=12, fontweight='bold', ha='center')
    ax.text(2.5, 7.8, '• Diarrhea → diluted samples', fontsize=10, ha='center')
    ax.text(2.5, 7.3, '• Inflammation → altered transit', fontsize=10, ha='center')
    ax.text(2.5, 6.8, '• Disease activity → load variation', fontsize=10, ha='center')

    # The consequence
    ax.add_patch(plt.Rectangle((5.5, 6.5), 4, 2.5, fill=True,
                               facecolor='#f39c12', alpha=0.2, edgecolor='#f39c12', linewidth=2))
    ax.text(7.5, 8.5, 'Artifact Risk', fontsize=12, fontweight='bold', ha='center')
    ax.text(7.5, 7.8, '• Load variation up to 10x', fontsize=10, ha='center')
    ax.text(7.5, 7.3, '• Creates ±3.3 log2FC artifacts', fontsize=10, ha='center')
    ax.text(7.5, 6.8, '• Mimics taxon-specific changes', fontsize=10, ha='center')

    # Arrow
    ax.annotate('', xy=(5.3, 7.75), xytext=(4.7, 7.75),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))

    # The solution
    ax.add_patch(plt.Rectangle((1, 3), 8, 2.5, fill=True,
                               facecolor='#2ecc71', alpha=0.2, edgecolor='#2ecc71', linewidth=2))
    ax.text(5, 5, 'Recommended Approach', fontsize=12, fontweight='bold', ha='center')
    ax.text(5, 4.3, '1. Report effect sizes, not just p-values', fontsize=10, ha='center')
    ax.text(5, 3.8, '2. Flag findings < 3.3 log2FC as "needs validation"', fontsize=10, ha='center')
    ax.text(5, 3.3, '3. Use spike-in controls when possible', fontsize=10, ha='center')

    # Arrows down
    ax.annotate('', xy=(2.5, 6.3), xytext=(2.5, 5.7),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))
    ax.annotate('', xy=(7.5, 6.3), xytext=(7.5, 5.7),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))

    # Outcome
    ax.add_patch(plt.Rectangle((2, 0.5), 6, 1.5, fill=True,
                               facecolor='#3498db', alpha=0.2, edgecolor='#3498db', linewidth=2))
    ax.text(5, 1.5, 'Robust, Reproducible Findings', fontsize=12, fontweight='bold', ha='center')
    ax.text(5, 0.9, 'Only large effects (>3.3 log2FC) or validated by spike-in', fontsize=10, ha='center')

    ax.annotate('', xy=(5, 2.8), xytext=(5, 2.2),
               arrowprops=dict(arrowstyle='->', color='gray', lw=2))

    plt.tight_layout()

    if save:
        fig.savefig(FIGURES_DIR / "fig4_clinical_implications.png")
        fig.savefig(FIGURES_DIR / "fig4_clinical_implications.pdf")
        print("Saved: fig4_clinical_implications.png/pdf")

    return fig


def main():
    """Generate all figures for Experiment 09."""
    print("=" * 60)
    print("Experiment 09: IBD Reanalysis Figures")
    print("=" * 60)

    # Load data
    data = load_data()
    print(f"\nLoaded data: {list(data.keys())}")

    # Generate figures
    print("\nGenerating figures...")

    fig1_artifact_risk_comparison(data)
    fig2_effect_size_distribution(data)
    fig3_fpr_validation(data)
    fig4_clinical_implications()

    print("\n" + "=" * 60)
    print("Figure generation complete!")
    print(f"Output directory: {FIGURES_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
