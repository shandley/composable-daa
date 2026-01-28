//! Abundance spike-in: inject fold-change effects into count data.

use crate::data::{CountMatrix, Metadata, PrevalenceTier, Variable};
use crate::error::{DaaError, Result};
use crate::profile::profile_prevalence;
use crate::spike::types::{SpikeDiagnostics, SpikeMode, SpikeSelection, SpikeSpec, SpikeType, SpikedData};
use sprs::TriMat;
use std::collections::HashSet;

/// Simple LCG random number generator for reproducibility.
struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.state
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Shuffle a vector in place.
    fn shuffle<T>(&mut self, vec: &mut [T]) {
        for i in (1..vec.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            vec.swap(i, j);
        }
    }
}

/// Spike abundance: multiply counts by fold_change for selected features in target group.
///
/// This injects a known multiplicative effect into the data, simulating
/// differential abundance. Only non-zero counts are modified (zeros remain zeros).
/// Uses Raw mode (simple multiplication) by default.
///
/// # Arguments
/// * `counts` - Original count matrix
/// * `metadata` - Sample metadata
/// * `group_column` - Column name defining groups
/// * `n_spike` - Number of features to spike
/// * `fold_change` - Multiplicative factor (e.g., 2.0 for 2x increase)
/// * `target_group` - Group to receive the increased abundance
/// * `selection` - How to select features for spiking
/// * `seed` - Random seed for reproducibility
///
/// # Returns
/// SpikedData containing modified counts and spike specification.
pub fn spike_abundance(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    n_spike: usize,
    fold_change: f64,
    target_group: &str,
    selection: SpikeSelection,
    seed: u64,
) -> Result<SpikedData> {
    spike_abundance_with_mode(
        counts,
        metadata,
        group_column,
        n_spike,
        fold_change,
        target_group,
        selection,
        seed,
        SpikeMode::Raw,
    )
}

/// Spike abundance with explicit mode control.
///
/// This version allows specifying how the spike should handle compositional constraints:
///
/// - **Raw**: Simple multiplication of counts. Library size increases for spiked samples.
/// - **Compositional**: Spike then renormalize to original library size. Models relative
///   abundance changes while keeping total reads constant.
/// - **Absolute**: Models what true absolute changes look like after sequencing at fixed
///   depth. The most biologically realistic but results in smaller observed effects.
///
/// # Arguments
/// * `counts` - Original count matrix
/// * `metadata` - Sample metadata
/// * `group_column` - Column name defining groups
/// * `n_spike` - Number of features to spike
/// * `fold_change` - Multiplicative factor (e.g., 2.0 for 2x increase)
/// * `target_group` - Group to receive the increased abundance
/// * `selection` - How to select features for spiking
/// * `seed` - Random seed for reproducibility
/// * `mode` - How to handle compositional constraints
///
/// # Returns
/// SpikedData containing modified counts, spike specification, and diagnostics.
pub fn spike_abundance_with_mode(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    n_spike: usize,
    fold_change: f64,
    target_group: &str,
    selection: SpikeSelection,
    seed: u64,
    mode: SpikeMode,
) -> Result<SpikedData> {
    if fold_change <= 0.0 {
        return Err(DaaError::InvalidParameter(
            "fold_change must be positive".into(),
        ));
    }
    if n_spike == 0 {
        return Err(DaaError::InvalidParameter(
            "n_spike must be at least 1".into(),
        ));
    }

    let mut rng = SimpleRng::new(seed);

    // Identify samples in target group
    let target_samples = get_group_samples(counts, metadata, group_column, target_group)?;
    if target_samples.is_empty() {
        return Err(DaaError::InvalidParameter(format!(
            "No samples found in group '{}'",
            target_group
        )));
    }

    // Select features to spike
    let eligible_features = get_eligible_features(counts, metadata, group_column, &selection)?;
    if eligible_features.len() < n_spike {
        return Err(DaaError::InvalidParameter(format!(
            "Only {} eligible features, but {} requested for spiking",
            eligible_features.len(),
            n_spike
        )));
    }

    // Randomly select n_spike features from eligible
    let mut indices: Vec<usize> = eligible_features.clone();
    rng.shuffle(&mut indices);
    let selected_features: Vec<usize> = indices.into_iter().take(n_spike).collect();

    // Calculate original prevalence for selected features
    let prevalence_profile = profile_prevalence(counts);
    let original_prevalence: Vec<f64> = selected_features
        .iter()
        .map(|&i| prevalence_profile.feature_prevalence[i])
        .collect();

    // Create modified count matrix based on mode
    let (spiked_counts, spiked_feature_ids, diagnostics) = apply_abundance_spike_with_mode(
        counts,
        &selected_features,
        &target_samples,
        fold_change,
        mode,
    )?;

    // Build effect sizes (same fold_change for all)
    let effect_sizes = vec![fold_change; n_spike];

    let spec = SpikeSpec::with_diagnostics(
        SpikeType::Abundance,
        spiked_feature_ids,
        effect_sizes,
        target_group.to_string(),
        original_prevalence,
        seed,
        mode,
        diagnostics,
    );

    Ok(SpikedData::new(spiked_counts, spec, counts.clone()))
}

/// Get sample indices belonging to a specific group.
fn get_group_samples(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    target_group: &str,
) -> Result<Vec<usize>> {
    let mut indices = Vec::new();
    for (idx, sample_id) in counts.sample_ids().iter().enumerate() {
        if let Some(var) = metadata.get(sample_id, group_column) {
            if let Variable::Categorical(level) = var {
                if level == target_group {
                    indices.push(idx);
                }
            }
        }
    }
    Ok(indices)
}

/// Get feature indices eligible for spiking based on selection criteria.
fn get_eligible_features(
    counts: &CountMatrix,
    _metadata: &Metadata,
    _group_column: &str,
    selection: &SpikeSelection,
) -> Result<Vec<usize>> {
    let n_features = counts.n_features();
    let prevalence_profile = profile_prevalence(counts);

    match selection {
        SpikeSelection::Random => {
            // All features with non-zero prevalence
            Ok((0..n_features)
                .filter(|&i| prevalence_profile.feature_prevalence[i] > 0.0)
                .collect())
        }
        SpikeSelection::ByPrevalenceTier(tier) => {
            Ok((0..n_features)
                .filter(|&i| {
                    let prev = prevalence_profile.feature_prevalence[i];
                    PrevalenceTier::from_prevalence(prev) == *tier
                })
                .collect())
        }
        SpikeSelection::ByAbundance { min, max } => {
            // Calculate relative abundance per feature
            let total: u64 = counts.row_sums().iter().sum();
            if total == 0 {
                return Err(DaaError::EmptyData("All counts are zero".into()));
            }
            let row_sums = counts.row_sums();
            Ok((0..n_features)
                .filter(|&i| {
                    let rel_abund = row_sums[i] as f64 / total as f64;
                    rel_abund >= *min && rel_abund <= *max
                })
                .collect())
        }
        SpikeSelection::Specific(ids) => {
            let id_set: HashSet<&str> = ids.iter().map(|s| s.as_str()).collect();
            let mut indices = Vec::new();
            for (i, fid) in counts.feature_ids().iter().enumerate() {
                if id_set.contains(fid.as_str()) {
                    indices.push(i);
                }
            }
            if indices.len() != ids.len() {
                return Err(DaaError::InvalidParameter(
                    "Some specified feature IDs not found".into(),
                ));
            }
            Ok(indices)
        }
    }
}

/// Apply the abundance spike to selected features with mode-specific handling.
fn apply_abundance_spike_with_mode(
    counts: &CountMatrix,
    feature_indices: &[usize],
    target_samples: &[usize],
    fold_change: f64,
    mode: SpikeMode,
) -> Result<(CountMatrix, Vec<String>, SpikeDiagnostics)> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let feature_set: HashSet<usize> = feature_indices.iter().cloned().collect();

    // Calculate original library sizes and geometric means for target samples
    let col_sums = counts.col_sums();
    let original_lib_sizes: Vec<u64> = target_samples.iter().map(|&s| col_sums[s]).collect();

    // Calculate original geometric mean (using pseudocount of 1 to handle zeros)
    let orig_geo_means: Vec<f64> = target_samples
        .iter()
        .map(|&s| {
            let mut log_sum = 0.0;
            for row in 0..n_features {
                let val = counts.get(row, s);
                log_sum += (val as f64 + 1.0).ln();
            }
            (log_sum / n_features as f64).exp()
        })
        .collect();
    let avg_orig_geo_mean: f64 = orig_geo_means.iter().sum::<f64>() / orig_geo_means.len() as f64;

    // First pass: apply raw spike to get intermediate counts
    let mut spiked_values: Vec<Vec<u64>> = Vec::with_capacity(n_features);
    for row in 0..n_features {
        let row_data = counts.row_dense(row);
        let mut new_row: Vec<u64> = row_data.clone();
        for &col in target_samples {
            if row_data[col] > 0 && feature_set.contains(&row) {
                new_row[col] = ((row_data[col] as f64) * fold_change).round() as u64;
            }
        }
        spiked_values.push(new_row);
    }

    // Apply mode-specific adjustments
    let library_size_factor = match mode {
        SpikeMode::Raw => {
            // No adjustment - library size increases
            1.0
        }
        SpikeMode::Compositional => {
            // Renormalize each target sample to original library size
            for (i, &col) in target_samples.iter().enumerate() {
                let new_lib_size: u64 = spiked_values.iter().map(|row| row[col]).sum();
                if new_lib_size > 0 {
                    let scale = original_lib_sizes[i] as f64 / new_lib_size as f64;
                    for row in spiked_values.iter_mut() {
                        row[col] = ((row[col] as f64) * scale).round() as u64;
                    }
                }
            }
            1.0 // Library size preserved
        }
        SpikeMode::Absolute => {
            // Model absolute increase: spiked taxa increase, others stay same in absolute terms,
            // then renormalize to original library size (simulating fixed sequencing depth)
            //
            // This is equivalent to: increase spiked features by fold_change, then divide
            // everything by the new total to get back to original library size.
            // The effective fold change for spiked features becomes:
            // new_rel = (old_abs * FC) / new_total = (old_rel * FC) / (1 + delta)
            // where delta = sum of increases / original total

            for (i, &col) in target_samples.iter().enumerate() {
                // Calculate how much the spiked features increased
                let original_spiked_sum: u64 = feature_indices
                    .iter()
                    .map(|&row| counts.get(row, col))
                    .sum();
                let new_spiked_sum: u64 = feature_indices
                    .iter()
                    .map(|&row| spiked_values[row][col])
                    .sum();
                let increase = new_spiked_sum.saturating_sub(original_spiked_sum);

                // New total if we imagine absolute abundances
                let hypothetical_total = original_lib_sizes[i] + increase;

                if hypothetical_total > 0 {
                    // Scale everything to original library size
                    let scale = original_lib_sizes[i] as f64 / hypothetical_total as f64;
                    for row in spiked_values.iter_mut() {
                        row[col] = ((row[col] as f64) * scale).round() as u64;
                    }
                }
            }
            1.0 // Library size preserved
        }
    };

    // Build the count matrix
    let mut tri_mat = TriMat::new((n_features, n_samples));
    for row in 0..n_features {
        for col in 0..n_samples {
            let val = spiked_values[row][col];
            if val > 0 {
                tri_mat.add_triplet(row, col, val);
            }
        }
    }

    let spiked_feature_ids: Vec<String> = feature_indices
        .iter()
        .map(|&i| counts.feature_ids()[i].clone())
        .collect();

    let new_counts = CountMatrix::new(
        tri_mat.to_csr(),
        counts.feature_ids().to_vec(),
        counts.sample_ids().to_vec(),
    )?;

    // Calculate post-spike geometric means
    let spiked_geo_means: Vec<f64> = target_samples
        .iter()
        .map(|&s| {
            let mut log_sum = 0.0;
            for row in 0..n_features {
                let val = new_counts.get(row, s);
                log_sum += (val as f64 + 1.0).ln();
            }
            (log_sum / n_features as f64).exp()
        })
        .collect();
    let avg_spiked_geo_mean: f64 = spiked_geo_means.iter().sum::<f64>() / spiked_geo_means.len() as f64;

    // Calculate effective CLR effect
    // CLR effect ≈ log(fold_change) - log(geo_mean_ratio)
    let geo_mean_ratio = avg_spiked_geo_mean / avg_orig_geo_mean;
    let effective_clr_effect = fold_change.ln() - geo_mean_ratio.ln();

    // Generate warning if compositional effects may dominate
    let spike_fraction = feature_indices.len() as f64 / n_features as f64;
    let compositional_warning = if spike_fraction > 0.1 {
        Some(format!(
            "Spiking {:.1}% of features may cause substantial compositional effects. \
             Consider spiking fewer features for cleaner evaluation.",
            spike_fraction * 100.0
        ))
    } else {
        None
    };

    let diagnostics = SpikeDiagnostics {
        original_geometric_mean: avg_orig_geo_mean,
        spiked_geometric_mean: avg_spiked_geo_mean,
        geometric_mean_ratio: geo_mean_ratio,
        nominal_fold_change: fold_change,
        effective_clr_effect,
        library_size_factor,
        n_spiked: feature_indices.len(),
        n_total_features: n_features,
        compositional_warning,
    };

    Ok((new_counts, spiked_feature_ids, diagnostics))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        let mut tri_mat = TriMat::new((5, 8));
        // Feature 0: present in all samples
        for s in 0..8 {
            tri_mat.add_triplet(0, s, 100);
        }
        // Feature 1: present in all samples
        for s in 0..8 {
            tri_mat.add_triplet(1, s, 50);
        }
        // Feature 2: present in 6/8 samples
        for s in 0..6 {
            tri_mat.add_triplet(2, s, 30);
        }
        // Feature 3: present in 4/8 samples
        for s in 0..4 {
            tri_mat.add_triplet(3, s, 20);
        }
        // Feature 4: rare, only 2/8 samples
        tri_mat.add_triplet(4, 0, 10);
        tri_mat.add_triplet(4, 1, 10);

        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..8).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        // First 4 samples: control, last 4: treatment
        for i in 0..4 {
            writeln!(file, "S{}\tcontrol", i).unwrap();
        }
        for i in 4..8 {
            writeln!(file, "S{}\ttreatment", i).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_spike_abundance_basic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let spiked = spike_abundance(
            &counts,
            &metadata,
            "group",
            2,       // spike 2 features
            2.0,     // 2x fold change
            "treatment",
            SpikeSelection::Random,
            42,
        ).unwrap();

        assert_eq!(spiked.spec.n_spiked(), 2);
        assert_eq!(spiked.spec.spike_type, SpikeType::Abundance);
        assert_eq!(spiked.spec.affected_group, "treatment");
        assert!(spiked.spec.effect_sizes.iter().all(|&e| e == 2.0));
    }

    #[test]
    fn test_spike_abundance_effect() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Spike feature 0 specifically
        let spiked = spike_abundance(
            &counts,
            &metadata,
            "group",
            1,
            3.0, // 3x fold change
            "treatment",
            SpikeSelection::Specific(vec!["feat_0".into()]),
            42,
        ).unwrap();

        // Check that treatment samples have 3x counts for feat_0
        let orig_treatment = counts.get(0, 4); // treatment sample
        let spiked_treatment = spiked.counts.get(0, 4);
        assert_eq!(spiked_treatment, orig_treatment * 3);

        // Control samples should be unchanged
        let orig_control = counts.get(0, 0);
        let spiked_control = spiked.counts.get(0, 0);
        assert_eq!(spiked_control, orig_control);
    }

    #[test]
    fn test_spike_by_prevalence_tier() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Spike only high-prevalence features
        let spiked = spike_abundance(
            &counts,
            &metadata,
            "group",
            1,
            2.0,
            "treatment",
            SpikeSelection::ByPrevalenceTier(PrevalenceTier::VeryHigh),
            42,
        ).unwrap();

        // Should select from features with >75% prevalence
        // feat_0 and feat_1 have 100% prevalence
        let spiked_id = &spiked.spec.spiked_features[0];
        assert!(spiked_id == "feat_0" || spiked_id == "feat_1");
    }

    #[test]
    fn test_spike_invalid_params() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Invalid fold_change
        assert!(spike_abundance(
            &counts, &metadata, "group", 1, 0.0, "treatment",
            SpikeSelection::Random, 42
        ).is_err());

        // Invalid n_spike
        assert!(spike_abundance(
            &counts, &metadata, "group", 0, 2.0, "treatment",
            SpikeSelection::Random, 42
        ).is_err());

        // Too many features requested
        assert!(spike_abundance(
            &counts, &metadata, "group", 100, 2.0, "treatment",
            SpikeSelection::Random, 42
        ).is_err());
    }

    #[test]
    fn test_spike_mode_raw() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let spiked = spike_abundance_with_mode(
            &counts,
            &metadata,
            "group",
            1,
            2.0,
            "treatment",
            SpikeSelection::Specific(vec!["feat_0".into()]),
            42,
            SpikeMode::Raw,
        ).unwrap();

        // Raw mode: library size should increase for treatment samples
        let orig_lib_size: u64 = counts.col_sums()[4..8].iter().sum();
        let spiked_lib_size: u64 = spiked.counts.col_sums()[4..8].iter().sum();
        assert!(spiked_lib_size > orig_lib_size, "Raw mode should increase library size");

        // Diagnostics should be present
        assert!(spiked.spec.diagnostics.is_some());
        let diag = spiked.spec.diagnostics.as_ref().unwrap();
        assert_eq!(diag.nominal_fold_change, 2.0);
        assert!(diag.library_size_factor == 1.0);
    }

    #[test]
    fn test_spike_mode_compositional() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let spiked = spike_abundance_with_mode(
            &counts,
            &metadata,
            "group",
            1,
            2.0,
            "treatment",
            SpikeSelection::Specific(vec!["feat_0".into()]),
            42,
            SpikeMode::Compositional,
        ).unwrap();

        // Compositional mode: library size should be preserved (approximately)
        for sample_idx in 4..8 {
            let orig = counts.col_sums()[sample_idx];
            let spiked_val = spiked.counts.col_sums()[sample_idx];
            // Allow small rounding differences
            let diff = (orig as f64 - spiked_val as f64).abs();
            assert!(diff < 5.0, "Compositional mode should preserve library size, got diff={}", diff);
        }

        // Spiked feature should still have increased relative to others
        let orig_feat0_treatment: u64 = (4..8).map(|s| counts.get(0, s)).sum();
        let spiked_feat0_treatment: u64 = (4..8).map(|s| spiked.counts.get(0, s)).sum();
        let orig_total_treatment: u64 = counts.col_sums()[4..8].iter().sum();
        let spiked_total_treatment: u64 = spiked.counts.col_sums()[4..8].iter().sum();

        let orig_proportion = orig_feat0_treatment as f64 / orig_total_treatment as f64;
        let spiked_proportion = spiked_feat0_treatment as f64 / spiked_total_treatment as f64;
        assert!(spiked_proportion > orig_proportion, "Spiked feature should have higher relative abundance");
    }

    #[test]
    fn test_spike_mode_absolute() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let spiked = spike_abundance_with_mode(
            &counts,
            &metadata,
            "group",
            1,
            2.0,
            "treatment",
            SpikeSelection::Specific(vec!["feat_0".into()]),
            42,
            SpikeMode::Absolute,
        ).unwrap();

        // Absolute mode: library size should be preserved
        for sample_idx in 4..8 {
            let orig = counts.col_sums()[sample_idx];
            let spiked_val = spiked.counts.col_sums()[sample_idx];
            let diff = (orig as f64 - spiked_val as f64).abs();
            assert!(diff < 5.0, "Absolute mode should preserve library size");
        }

        // Non-spiked features should decrease (due to compositional closure)
        let orig_feat1_treatment: u64 = (4..8).map(|s| counts.get(1, s)).sum();
        let spiked_feat1_treatment: u64 = (4..8).map(|s| spiked.counts.get(1, s)).sum();
        assert!(spiked_feat1_treatment < orig_feat1_treatment,
            "Absolute mode: non-spiked features should decrease");
    }

    #[test]
    fn test_spike_diagnostics() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Use Raw mode so we can observe the geometric mean shift clearly
        let spiked = spike_abundance_with_mode(
            &counts,
            &metadata,
            "group",
            2,  // Spike 2 features
            3.0,
            "treatment",
            SpikeSelection::Random,
            42,
            SpikeMode::Raw,  // Raw mode to see geometric mean increase
        ).unwrap();

        let diag = spiked.spec.diagnostics.as_ref().unwrap();

        // Check diagnostic fields are populated
        assert!(diag.original_geometric_mean > 0.0);
        assert!(diag.spiked_geometric_mean > 0.0);
        assert!(diag.geometric_mean_ratio > 0.0);
        assert_eq!(diag.nominal_fold_change, 3.0);
        assert_eq!(diag.n_spiked, 2);
        assert_eq!(diag.n_total_features, 5);

        // In Raw mode with spiking, the geometric mean should increase
        // (spiked values go up, increasing the average)
        assert!(diag.geometric_mean_ratio >= 1.0,
            "Geometric mean should increase or stay same when spiking in Raw mode, got ratio={}",
            diag.geometric_mean_ratio);

        // The effective CLR effect accounts for the geometric mean shift
        // It should be positive (we're increasing features) but may differ from log(FC)
        assert!(diag.effective_clr_effect.is_finite(),
            "Effective CLR effect should be finite");
    }

    #[test]
    fn test_spike_compositional_warning() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Spike more than 10% of features (1 out of 5 = 20%)
        let spiked = spike_abundance_with_mode(
            &counts,
            &metadata,
            "group",
            1,
            2.0,
            "treatment",
            SpikeSelection::Random,
            42,
            SpikeMode::Raw,
        ).unwrap();

        let diag = spiked.spec.diagnostics.as_ref().unwrap();
        assert!(diag.compositional_warning.is_some(),
            "Should warn when spiking >10% of features");
    }
}
