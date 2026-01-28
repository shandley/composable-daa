//! Abundance spike-in: inject fold-change effects into count data.

use crate::data::{CountMatrix, Metadata, PrevalenceTier, Variable};
use crate::error::{DaaError, Result};
use crate::profile::profile_prevalence;
use crate::spike::types::{SpikeSelection, SpikeSpec, SpikeType, SpikedData};
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

    // Create modified count matrix
    let (spiked_counts, spiked_feature_ids) = apply_abundance_spike(
        counts,
        &selected_features,
        &target_samples,
        fold_change,
    )?;

    // Build effect sizes (same fold_change for all)
    let effect_sizes = vec![fold_change; n_spike];

    let spec = SpikeSpec::new(
        SpikeType::Abundance,
        spiked_feature_ids,
        effect_sizes,
        target_group.to_string(),
        original_prevalence,
        seed,
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

/// Apply the abundance spike to selected features.
fn apply_abundance_spike(
    counts: &CountMatrix,
    feature_indices: &[usize],
    target_samples: &[usize],
    fold_change: f64,
) -> Result<(CountMatrix, Vec<String>)> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let target_set: HashSet<usize> = target_samples.iter().cloned().collect();
    let feature_set: HashSet<usize> = feature_indices.iter().cloned().collect();

    let mut tri_mat = TriMat::new((n_features, n_samples));

    for row in 0..n_features {
        let row_data = counts.row_dense(row);
        for col in 0..n_samples {
            let mut val = row_data[col];
            if val > 0 && feature_set.contains(&row) && target_set.contains(&col) {
                // Apply fold change
                val = ((val as f64) * fold_change).round() as u64;
            }
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

    Ok((new_counts, spiked_feature_ids))
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
}
