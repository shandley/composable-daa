//! Integration tests for the LinDA-style pipeline.

use composable_daa::prelude::*;
use sprs::TriMat;
use std::io::Write;
use tempfile::NamedTempFile;

/// Create synthetic count data with known group effects.
fn create_synthetic_counts() -> CountMatrix {
    // 20 features × 40 samples (20 per group for better power)
    // - Features 0-4: strong treatment effect (5x increase)
    // - Features 5-9: moderate treatment effect (2x increase)
    // - Features 10-14: no effect (similar in both groups)
    // - Features 15-17: present in 50% of samples (no effect)
    // - Features 18-19: rare (present in 5% of samples)
    let n_features = 20;
    let n_samples = 40;
    let mut tri_mat = TriMat::new((n_features, n_samples));

    let mut rng_seed = 42u64;
    let simple_rand = |seed: &mut u64| -> f64 {
        *seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
        ((*seed >> 16) & 0x7FFF) as f64 / 32768.0
    };

    for feat in 0..n_features {
        for sample in 0..n_samples {
            let is_treatment = sample >= 20;

            let (base_count, effect_mult) = match feat {
                0..=4 => (100.0, if is_treatment { 5.0 } else { 1.0 }), // Strong 5x effect
                5..=9 => (150.0, if is_treatment { 2.0 } else { 1.0 }), // 2x effect
                10..=14 => (200.0, 1.0), // no effect
                15..=17 => {
                    // 50% prevalence
                    if sample % 2 == 0 {
                        (80.0, 1.0)
                    } else {
                        continue;
                    }
                }
                18..=19 => {
                    // Very rare: 5% prevalence (2 samples out of 40)
                    if sample < 2 {
                        (50.0, 1.0)
                    } else {
                        continue;
                    }
                }
                _ => unreachable!(),
            };

            // Add less noise for cleaner signal
            let noise = 0.9 + 0.2 * simple_rand(&mut rng_seed);
            let count = (base_count * effect_mult * noise).round() as u64;
            if count > 0 {
                tri_mat.add_triplet(feat, sample, count);
            }
        }
    }

    let feature_ids: Vec<String> = (0..n_features).map(|i| format!("taxon_{}", i)).collect();
    let sample_ids: Vec<String> = (0..n_samples).map(|i| format!("sample_{}", i)).collect();
    CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
}

/// Create metadata matching the synthetic counts.
fn create_synthetic_metadata() -> Metadata {
    let mut file = NamedTempFile::new().unwrap();
    writeln!(file, "sample_id\tgroup\tage\tbatch").unwrap();
    for i in 0..40 {
        let group = if i < 20 { "control" } else { "treatment" };
        let age = 25 + (i % 20) * 2;
        let batch = if i < 10 || (i >= 20 && i < 30) { "A" } else { "B" };
        writeln!(file, "sample_{}\t{}\t{}\t{}", i, group, age, batch).unwrap();
    }
    file.flush().unwrap();
    Metadata::from_tsv(file.path()).unwrap()
}

#[test]
fn test_full_linda_pipeline() {
    let counts = create_synthetic_counts();
    let metadata = create_synthetic_metadata();

    // Run the pipeline
    let results = Pipeline::new()
        .name("LinDA-test")
        .filter_prevalence(0.3)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh()
        .run(&counts, &metadata)
        .unwrap();

    // Basic checks
    assert!(!results.is_empty(), "Should have results");

    // Check that strongly affected features have lower q-values
    let significant = results.significant();
    assert!(
        significant.len() >= 3,
        "Should detect some significant features, got {}",
        significant.len()
    );

    // The strongly affected features (0-4) should have positive estimates
    // (treatment increases abundance)
    for r in &results.results {
        if r.feature_id.starts_with("taxon_") {
            let feat_num: usize = r.feature_id[6..].parse().unwrap_or(100);
            if feat_num < 5 && r.q_value < 0.1 {
                assert!(
                    r.estimate > 0.0,
                    "Feature {} should have positive effect, got {}",
                    r.feature_id,
                    r.estimate
                );
            }
        }
    }
}

#[test]
fn test_pipeline_with_covariate() {
    let counts = create_synthetic_counts();
    let metadata = create_synthetic_metadata();

    // Run with age as covariate
    let results = Pipeline::new()
        .filter_prevalence(0.3)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group + age")
        .test_wald("grouptreatment")
        .correct_bh()
        .run(&counts, &metadata)
        .unwrap();

    assert!(!results.is_empty());
}

#[test]
fn test_profile_before_pipeline() {
    let counts = create_synthetic_counts();

    // Profile the data
    let sparsity = profile_sparsity(&counts);
    let prevalence = profile_prevalence(&counts);
    let lib_size = profile_library_size(&counts);

    // Check profiling results
    assert_eq!(sparsity.total_entries, 20 * 40);
    assert!(sparsity.sparsity > 0.0, "Should have some zeros");
    assert!(sparsity.sparsity < 0.5, "Should not be too sparse");

    assert_eq!(prevalence.n_features, 20);
    assert!(prevalence.n_rare > 0, "Should have some rare features (<10%)");

    assert_eq!(lib_size.n_samples, 40);
    assert!(lib_size.mean > 0.0);
}

#[test]
fn test_groupwise_filtering() {
    let counts = create_synthetic_counts();
    let metadata = create_synthetic_metadata();

    // Filter requiring presence in both groups
    let results = Pipeline::new()
        .filter_prevalence_groupwise(0.3, "group", GroupwiseLogic::All)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh()
        .run(&counts, &metadata)
        .unwrap();

    assert!(!results.is_empty());
}

#[test]
fn test_run_linda_convenience() {
    let counts = create_synthetic_counts();
    let metadata = create_synthetic_metadata();

    let results = composable_daa::pipeline::run_linda(
        &counts,
        &metadata,
        "~ group",
        "grouptreatment",
        0.3,
        0.5,
    )
    .unwrap();

    assert_eq!(results.method, "LinDA");
    assert!(!results.is_empty());
}

#[test]
fn test_tsv_roundtrip() {
    let counts = create_synthetic_counts();

    // Write to TSV
    let temp = NamedTempFile::new().unwrap();
    counts.to_tsv(temp.path()).unwrap();

    // Read back
    let loaded = CountMatrix::from_tsv(temp.path()).unwrap();

    // Verify
    assert_eq!(loaded.n_features(), counts.n_features());
    assert_eq!(loaded.n_samples(), counts.n_samples());
    assert_eq!(loaded.feature_ids(), counts.feature_ids());
}

#[test]
fn test_result_output() {
    let counts = create_synthetic_counts();
    let metadata = create_synthetic_metadata();

    let results = Pipeline::new()
        .filter_prevalence(0.5)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh()
        .run(&counts, &metadata)
        .unwrap();

    // Write results to TSV
    let temp = NamedTempFile::new().unwrap();
    results.to_tsv(temp.path()).unwrap();

    // Verify file was written
    let content = std::fs::read_to_string(temp.path()).unwrap();
    assert!(content.contains("feature_id"));
    assert!(content.contains("q_value"));
}
