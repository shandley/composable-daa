#!/usr/bin/env python3
"""
Generate synthetic microbiome count data for testing DAA pipelines.

Creates datasets with known differential features for validation.
"""

import numpy as np
import argparse
from pathlib import Path


def generate_synthetic_counts(
    n_features: int = 200,
    n_samples: int = 40,
    n_differential: int = 20,
    fold_change: float = 2.0,
    sparsity: float = 0.3,
    base_depth: int = 10000,
    depth_variation: float = 0.3,
    seed: int = 42,
) -> tuple[np.ndarray, list[str], list[str], list[str], list[str]]:
    """
    Generate synthetic count matrix with known differential features.

    Returns:
        counts: n_features x n_samples count matrix
        feature_ids: list of feature names
        sample_ids: list of sample names
        groups: list of group assignments ('control' or 'treatment')
        differential_features: list of features with true differential abundance
    """
    rng = np.random.default_rng(seed)

    # Sample assignments (balanced groups)
    n_control = n_samples // 2
    n_treatment = n_samples - n_control
    groups = ['control'] * n_control + ['treatment'] * n_treatment

    # Generate base abundances (log-normal distribution, typical of microbiome)
    log_abundances = rng.normal(loc=2, scale=2, size=n_features)
    base_abundances = np.exp(log_abundances)
    base_abundances = base_abundances / base_abundances.sum()  # Normalize to proportions

    # Select differential features (from medium-abundance features for realism)
    abundance_order = np.argsort(base_abundances)
    # Pick from middle 60% of abundance range
    mid_start = int(n_features * 0.2)
    mid_end = int(n_features * 0.8)
    candidate_indices = abundance_order[mid_start:mid_end]
    differential_indices = rng.choice(candidate_indices, size=n_differential, replace=False)

    # Generate counts
    counts = np.zeros((n_features, n_samples), dtype=np.int64)

    for j in range(n_samples):
        # Vary library size
        depth = int(base_depth * (1 + rng.uniform(-depth_variation, depth_variation)))

        # Adjust abundances for differential features in treatment
        abundances = base_abundances.copy()
        if groups[j] == 'treatment':
            # Half up-regulated, half down-regulated
            for i, idx in enumerate(differential_indices):
                if i < n_differential // 2:
                    abundances[idx] *= fold_change
                else:
                    abundances[idx] /= fold_change

        # Renormalize
        abundances = abundances / abundances.sum()

        # Add sparsity (zero-inflation)
        zero_mask = rng.random(n_features) < sparsity
        abundances[zero_mask] = 0
        if abundances.sum() > 0:
            abundances = abundances / abundances.sum()

        # Generate counts from multinomial
        counts[:, j] = rng.multinomial(depth, abundances)

    # Generate IDs
    feature_ids = [f"Feature_{i:03d}" for i in range(n_features)]
    sample_ids = [f"Sample_{j:03d}" for j in range(n_samples)]
    differential_features = [feature_ids[i] for i in differential_indices]

    return counts, feature_ids, sample_ids, groups, differential_features


def write_count_matrix(counts: np.ndarray, feature_ids: list, sample_ids: list, path: Path):
    """Write count matrix to TSV."""
    with open(path, 'w') as f:
        # Header
        f.write("feature_id\t" + "\t".join(sample_ids) + "\n")
        # Data
        for i, fid in enumerate(feature_ids):
            row = "\t".join(str(c) for c in counts[i, :])
            f.write(f"{fid}\t{row}\n")


def write_metadata(sample_ids: list, groups: list, path: Path, add_covariates: bool = True):
    """Write metadata to TSV."""
    rng = np.random.default_rng(123)

    with open(path, 'w') as f:
        if add_covariates:
            f.write("sample_id\tgroup\tage\tbatch\n")
            batches = ['A', 'B']
            for sid, group in zip(sample_ids, groups):
                age = int(rng.normal(45, 15))
                age = max(18, min(80, age))  # Clamp to reasonable range
                batch = rng.choice(batches)
                f.write(f"{sid}\t{group}\t{age}\t{batch}\n")
        else:
            f.write("sample_id\tgroup\n")
            for sid, group in zip(sample_ids, groups):
                f.write(f"{sid}\t{group}\n")


def write_ground_truth(differential_features: list, fold_change: float, path: Path):
    """Write ground truth differential features."""
    n_diff = len(differential_features)
    with open(path, 'w') as f:
        f.write("feature_id\tdirection\tfold_change\n")
        for i, fid in enumerate(differential_features):
            if i < n_diff // 2:
                f.write(f"{fid}\tup\t{fold_change}\n")
            else:
                f.write(f"{fid}\tdown\t{1/fold_change:.4f}\n")


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic microbiome data")
    parser.add_argument("--output-dir", "-o", type=Path, default=Path("test_data"),
                        help="Output directory")
    parser.add_argument("--n-features", type=int, default=200,
                        help="Number of features/taxa")
    parser.add_argument("--n-samples", type=int, default=40,
                        help="Number of samples")
    parser.add_argument("--n-differential", type=int, default=20,
                        help="Number of truly differential features")
    parser.add_argument("--fold-change", type=float, default=2.0,
                        help="Fold change for differential features")
    parser.add_argument("--sparsity", type=float, default=0.3,
                        help="Proportion of zeros (0-1)")
    parser.add_argument("--base-depth", type=int, default=10000,
                        help="Mean sequencing depth per sample")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed")

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Generating synthetic data with:")
    print(f"  Features: {args.n_features}")
    print(f"  Samples: {args.n_samples}")
    print(f"  Differential features: {args.n_differential}")
    print(f"  Fold change: {args.fold_change}x")
    print(f"  Sparsity: {args.sparsity:.0%}")
    print(f"  Base depth: {args.base_depth}")
    print()

    # Generate data
    counts, feature_ids, sample_ids, groups, diff_features = generate_synthetic_counts(
        n_features=args.n_features,
        n_samples=args.n_samples,
        n_differential=args.n_differential,
        fold_change=args.fold_change,
        sparsity=args.sparsity,
        base_depth=args.base_depth,
        seed=args.seed,
    )

    # Write files
    counts_path = args.output_dir / "counts.tsv"
    metadata_path = args.output_dir / "metadata.tsv"
    truth_path = args.output_dir / "ground_truth.tsv"

    write_count_matrix(counts, feature_ids, sample_ids, counts_path)
    write_metadata(sample_ids, groups, metadata_path)
    write_ground_truth(diff_features, args.fold_change, truth_path)

    print(f"Written:")
    print(f"  {counts_path}")
    print(f"  {metadata_path}")
    print(f"  {truth_path}")
    print()

    # Summary statistics
    total_counts = counts.sum()
    zero_fraction = (counts == 0).sum() / counts.size
    lib_sizes = counts.sum(axis=0)

    print(f"Data summary:")
    print(f"  Total counts: {total_counts:,}")
    print(f"  Zero fraction: {zero_fraction:.1%}")
    print(f"  Library size range: {lib_sizes.min():,} - {lib_sizes.max():,}")
    print(f"  Differential features: {len(diff_features)}")
    print(f"    Up-regulated: {len(diff_features)//2}")
    print(f"    Down-regulated: {len(diff_features) - len(diff_features)//2}")


if __name__ == "__main__":
    main()
