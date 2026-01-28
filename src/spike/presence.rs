//! Presence spike-in: inject presence/absence effects into count data.

use crate::data::{CountMatrix, Metadata, PrevalenceTier, Variable};
use crate::error::{DaaError, Result};
use crate::profile::profile_prevalence;
use crate::spike::types::{AbundanceLevel, SpikeSelection, SpikeSpec, SpikeType, SpikedData};
use sprs::TriMat;
use std::collections::HashSet;

/// Simple LCG random number generator.
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

    fn shuffle<T>(&mut self, vec: &mut [T]) {
        for i in (1..vec.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            vec.swap(i, j);
        }
    }
}

/// Spike presence: convert zeros to non-zeros in target group.
///
/// This injects a known presence effect, simulating features that appear
/// more frequently in one group. Zeros in the target group are converted
/// to non-zero counts at the specified abundance level.
///
/// # Arguments
/// * `counts` - Original count matrix
/// * `metadata` - Sample metadata
/// * `group_column` - Column name defining groups
/// * `n_spike` - Number of features to spike
/// * `target_prevalence` - Desired prevalence in target group (0.0-1.0)
/// * `target_group` - Group to receive the presence increase
/// * `abundance_level` - What count value to use for new non-zeros
/// * `selection` - How to select features for spiking
/// * `seed` - Random seed for reproducibility
///
/// # Returns
/// SpikedData containing modified counts and spike specification.
pub fn spike_presence(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    n_spike: usize,
    target_prevalence: f64,
    target_group: &str,
    abundance_level: AbundanceLevel,
    selection: SpikeSelection,
    seed: u64,
) -> Result<SpikedData> {
    if !(0.0..=1.0).contains(&target_prevalence) {
        return Err(DaaError::InvalidParameter(
            "target_prevalence must be between 0 and 1".into(),
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

    // Select features to spike (prefer features with lower prevalence for presence spikes)
    let eligible_features = get_eligible_features_for_presence(counts, &target_samples, &selection)?;
    if eligible_features.len() < n_spike {
        return Err(DaaError::InvalidParameter(format!(
            "Only {} eligible features, but {} requested for spiking",
            eligible_features.len(),
            n_spike
        )));
    }

    // Randomly select n_spike features
    let mut indices: Vec<usize> = eligible_features.clone();
    rng.shuffle(&mut indices);
    let selected_features: Vec<usize> = indices.into_iter().take(n_spike).collect();

    // Calculate original prevalence (in target group) for selected features
    let original_prevalence: Vec<f64> = selected_features
        .iter()
        .map(|&feat_idx| {
            let row = counts.row_dense(feat_idx);
            let present = target_samples.iter().filter(|&&s| row[s] > 0).count();
            present as f64 / target_samples.len() as f64
        })
        .collect();

    // Determine abundance value for new non-zeros
    let fill_value = calculate_fill_value(counts, &abundance_level);

    // Create modified count matrix
    let (spiked_counts, spiked_feature_ids, prevalence_increases) = apply_presence_spike(
        counts,
        &selected_features,
        &target_samples,
        target_prevalence,
        fill_value,
        &original_prevalence,
        &mut rng,
    )?;

    let spec = SpikeSpec::new(
        SpikeType::Presence,
        spiked_feature_ids,
        prevalence_increases,
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

/// Get features eligible for presence spiking.
/// For presence spikes, we want features that have room to increase prevalence.
fn get_eligible_features_for_presence(
    counts: &CountMatrix,
    target_samples: &[usize],
    selection: &SpikeSelection,
) -> Result<Vec<usize>> {
    let n_features = counts.n_features();
    let prevalence_profile = profile_prevalence(counts);

    match selection {
        SpikeSelection::Random => {
            // Features that have at least one zero in target group
            Ok((0..n_features)
                .filter(|&feat| {
                    let row = counts.row_dense(feat);
                    target_samples.iter().any(|&s| row[s] == 0)
                })
                .collect())
        }
        SpikeSelection::ByPrevalenceTier(tier) => {
            Ok((0..n_features)
                .filter(|&i| {
                    let prev = prevalence_profile.feature_prevalence[i];
                    let tier_match = PrevalenceTier::from_prevalence(prev) == *tier;
                    let row = counts.row_dense(i);
                    let has_zeros = target_samples.iter().any(|&s| row[s] == 0);
                    tier_match && has_zeros
                })
                .collect())
        }
        SpikeSelection::ByAbundance { min, max } => {
            let total: u64 = counts.row_sums().iter().sum();
            if total == 0 {
                return Err(DaaError::EmptyData("All counts are zero".into()));
            }
            let row_sums = counts.row_sums();
            Ok((0..n_features)
                .filter(|&i| {
                    let rel_abund = row_sums[i] as f64 / total as f64;
                    let in_range = rel_abund >= *min && rel_abund <= *max;
                    let row = counts.row_dense(i);
                    let has_zeros = target_samples.iter().any(|&s| row[s] == 0);
                    in_range && has_zeros
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

/// Calculate the fill value for new non-zeros based on abundance level.
fn calculate_fill_value(counts: &CountMatrix, level: &AbundanceLevel) -> u64 {
    match level {
        AbundanceLevel::Fixed(v) => *v,
        AbundanceLevel::Low | AbundanceLevel::Median | AbundanceLevel::High => {
            // Collect all non-zero values
            let mut nonzeros: Vec<u64> = Vec::new();
            for row in 0..counts.n_features() {
                for col in 0..counts.n_samples() {
                    let val = counts.get(row, col);
                    if val > 0 {
                        nonzeros.push(val);
                    }
                }
            }
            if nonzeros.is_empty() {
                return 1;
            }
            nonzeros.sort();
            let n = nonzeros.len();
            match level {
                AbundanceLevel::Low => nonzeros[n / 10], // 10th percentile
                AbundanceLevel::Median => nonzeros[n / 2],
                AbundanceLevel::High => nonzeros[n * 9 / 10], // 90th percentile
                _ => unreachable!(),
            }
        }
    }
}

/// Apply the presence spike to selected features.
fn apply_presence_spike(
    counts: &CountMatrix,
    feature_indices: &[usize],
    target_samples: &[usize],
    target_prevalence: f64,
    fill_value: u64,
    original_prevalence: &[f64],
    rng: &mut SimpleRng,
) -> Result<(CountMatrix, Vec<String>, Vec<f64>)> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let n_target = target_samples.len();
    let target_set: HashSet<usize> = target_samples.iter().cloned().collect();
    let feature_set: HashSet<usize> = feature_indices.iter().cloned().collect();

    let mut tri_mat = TriMat::new((n_features, n_samples));
    let mut prevalence_increases = Vec::with_capacity(feature_indices.len());

    for (feat_list_idx, &feat_idx) in feature_indices.iter().enumerate() {
        // Calculate how many samples need to become non-zero
        let current_present = (original_prevalence[feat_list_idx] * n_target as f64).round() as usize;
        let target_present = (target_prevalence * n_target as f64).ceil() as usize;
        let need_to_fill = target_present.saturating_sub(current_present);

        // Find zeros to fill
        let row = counts.row_dense(feat_idx);
        let zero_samples: Vec<usize> = target_samples
            .iter()
            .filter(|&&s| row[s] == 0)
            .cloned()
            .collect();

        // Randomly select which zeros to fill
        let mut zeros_to_fill = zero_samples.clone();
        rng.shuffle(&mut zeros_to_fill);
        let fill_set: HashSet<usize> = zeros_to_fill.into_iter().take(need_to_fill).collect();

        // Calculate actual prevalence increase
        let new_present = current_present + fill_set.len();
        let new_prevalence = new_present as f64 / n_target as f64;
        prevalence_increases.push(new_prevalence - original_prevalence[feat_list_idx]);

        // Build row with filled values
        for col in 0..n_samples {
            let mut val = row[col];
            if feature_set.contains(&feat_idx) && target_set.contains(&col) && fill_set.contains(&col) {
                val = fill_value;
            }
            if val > 0 {
                tri_mat.add_triplet(feat_idx, col, val);
            }
        }
    }

    // Copy non-spiked features unchanged
    for feat_idx in 0..n_features {
        if feature_set.contains(&feat_idx) {
            continue;
        }
        let row = counts.row_dense(feat_idx);
        for col in 0..n_samples {
            if row[col] > 0 {
                tri_mat.add_triplet(feat_idx, col, row[col]);
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

    Ok((new_counts, spiked_feature_ids, prevalence_increases))
}

/// Spike both presence AND abundance (hurdle model validation).
///
/// This combines a presence spike with an abundance spike, simulating
/// features that both appear more frequently and have higher counts
/// when present in the target group.
pub fn spike_hurdle(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
    n_spike: usize,
    presence_effect: f64,  // target prevalence increase (additive)
    abundance_effect: f64, // fold change for non-zeros
    target_group: &str,
    selection: SpikeSelection,
    seed: u64,
) -> Result<SpikedData> {
    if !(0.0..=1.0).contains(&presence_effect) {
        return Err(DaaError::InvalidParameter(
            "presence_effect must be between 0 and 1".into(),
        ));
    }
    if abundance_effect <= 0.0 {
        return Err(DaaError::InvalidParameter(
            "abundance_effect must be positive".into(),
        ));
    }

    let mut rng = SimpleRng::new(seed);

    // Get target samples
    let target_samples = get_group_samples(counts, metadata, group_column, target_group)?;
    if target_samples.is_empty() {
        return Err(DaaError::InvalidParameter(format!(
            "No samples found in group '{}'",
            target_group
        )));
    }

    // Select eligible features
    let eligible_features = get_eligible_features_for_presence(counts, &target_samples, &selection)?;
    if eligible_features.len() < n_spike {
        return Err(DaaError::InvalidParameter(format!(
            "Only {} eligible features, but {} requested",
            eligible_features.len(),
            n_spike
        )));
    }

    let mut indices: Vec<usize> = eligible_features;
    rng.shuffle(&mut indices);
    let selected_features: Vec<usize> = indices.into_iter().take(n_spike).collect();

    // Calculate original prevalence
    let original_prevalence: Vec<f64> = selected_features
        .iter()
        .map(|&feat_idx| {
            let row = counts.row_dense(feat_idx);
            let present = target_samples.iter().filter(|&&s| row[s] > 0).count();
            present as f64 / target_samples.len() as f64
        })
        .collect();

    // Apply combined spike
    let fill_value = calculate_fill_value(counts, &AbundanceLevel::Median);
    let (spiked_counts, spiked_feature_ids, combined_effects) = apply_hurdle_spike(
        counts,
        &selected_features,
        &target_samples,
        presence_effect,
        abundance_effect,
        fill_value,
        &original_prevalence,
        &mut rng,
    )?;

    let spec = SpikeSpec::new(
        SpikeType::Hurdle,
        spiked_feature_ids,
        combined_effects,
        target_group.to_string(),
        original_prevalence,
        seed,
    );

    Ok(SpikedData::new(spiked_counts, spec, counts.clone()))
}

fn apply_hurdle_spike(
    counts: &CountMatrix,
    feature_indices: &[usize],
    target_samples: &[usize],
    presence_effect: f64,
    abundance_effect: f64,
    fill_value: u64,
    original_prevalence: &[f64],
    rng: &mut SimpleRng,
) -> Result<(CountMatrix, Vec<String>, Vec<f64>)> {
    let n_features = counts.n_features();
    let n_samples = counts.n_samples();
    let n_target = target_samples.len();
    let target_set: HashSet<usize> = target_samples.iter().cloned().collect();
    let feature_set: HashSet<usize> = feature_indices.iter().cloned().collect();

    let mut tri_mat = TriMat::new((n_features, n_samples));
    let mut combined_effects = Vec::with_capacity(feature_indices.len());

    for (feat_list_idx, &feat_idx) in feature_indices.iter().enumerate() {
        // Calculate presence spike
        let current_prev = original_prevalence[feat_list_idx];
        let target_prev = (current_prev + presence_effect).min(1.0);
        let current_present = (current_prev * n_target as f64).round() as usize;
        let target_present = (target_prev * n_target as f64).ceil() as usize;
        let need_to_fill = target_present.saturating_sub(current_present);

        let row = counts.row_dense(feat_idx);
        let zero_samples: Vec<usize> = target_samples
            .iter()
            .filter(|&&s| row[s] == 0)
            .cloned()
            .collect();

        let mut zeros_to_fill = zero_samples;
        rng.shuffle(&mut zeros_to_fill);
        let fill_set: HashSet<usize> = zeros_to_fill.into_iter().take(need_to_fill).collect();

        // Combined effect: presence increase * abundance fold change
        let prev_increase = fill_set.len() as f64 / n_target as f64;
        combined_effects.push(prev_increase + (abundance_effect - 1.0));

        // Build row with both presence and abundance effects
        for col in 0..n_samples {
            let mut val = row[col];
            if target_set.contains(&col) {
                if fill_set.contains(&col) {
                    val = fill_value;
                }
                if val > 0 {
                    val = ((val as f64) * abundance_effect).round() as u64;
                }
            }
            if val > 0 {
                tri_mat.add_triplet(feat_idx, col, val);
            }
        }
    }

    // Copy non-spiked features
    for feat_idx in 0..n_features {
        if feature_set.contains(&feat_idx) {
            continue;
        }
        let row = counts.row_dense(feat_idx);
        for col in 0..n_samples {
            if row[col] > 0 {
                tri_mat.add_triplet(feat_idx, col, row[col]);
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

    Ok((new_counts, spiked_feature_ids, combined_effects))
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
        // Feature 1: present in 4/8 (50%) - some in each group
        for s in 0..4 {
            tri_mat.add_triplet(1, s, 50);
        }
        // Feature 2: present only in control (0-3)
        for s in 0..4 {
            tri_mat.add_triplet(2, s, 30);
        }
        // Feature 3: present in 2/8 (25%)
        tri_mat.add_triplet(3, 0, 20);
        tri_mat.add_triplet(3, 1, 20);
        // Feature 4: rare
        tri_mat.add_triplet(4, 0, 10);

        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..8).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
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
    fn test_spike_presence_basic() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Spike feature 2 (only in control) to appear in treatment
        let spiked = spike_presence(
            &counts,
            &metadata,
            "group",
            1,
            0.75, // target 75% prevalence in treatment
            "treatment",
            AbundanceLevel::Median,
            SpikeSelection::Specific(vec!["feat_2".into()]),
            42,
        ).unwrap();

        assert_eq!(spiked.spec.n_spiked(), 1);
        assert_eq!(spiked.spec.spike_type, SpikeType::Presence);

        // Check that feature 2 now has non-zeros in treatment group
        let orig_treatment_nz: usize = (4..8).filter(|&s| counts.get(2, s) > 0).count();
        let spiked_treatment_nz: usize = (4..8).filter(|&s| spiked.counts.get(2, s) > 0).count();
        assert!(spiked_treatment_nz > orig_treatment_nz);
    }

    #[test]
    fn test_spike_hurdle() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let spiked = spike_hurdle(
            &counts,
            &metadata,
            "group",
            1,
            0.5,  // increase prevalence by 50%
            2.0,  // 2x abundance
            "treatment",
            SpikeSelection::Specific(vec!["feat_2".into()]),
            42,
        ).unwrap();

        assert_eq!(spiked.spec.spike_type, SpikeType::Hurdle);

        // Should have both more presence and higher abundance
        let spiked_treatment_nz: usize = (4..8).filter(|&s| spiked.counts.get(2, s) > 0).count();
        assert!(spiked_treatment_nz > 0);
    }
}
