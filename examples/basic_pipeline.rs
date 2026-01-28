//! Basic example demonstrating the LinDA-style pipeline.
//!
//! This example shows how to:
//! 1. Create synthetic data
//! 2. Profile the data
//! 3. Run a differential abundance pipeline
//! 4. Examine results

use composable_daa::prelude::*;
use sprs::TriMat;

fn main() -> Result<()> {
    println!("=== Composable DAA Example ===\n");

    // Create synthetic data
    let (counts, metadata) = create_example_data();

    println!("Data dimensions:");
    println!("  Features: {}", counts.n_features());
    println!("  Samples:  {}", counts.n_samples());
    println!();

    // Profile the data
    println!("=== Data Profiling ===\n");

    let sparsity = profile_sparsity(&counts);
    println!("Sparsity: {:.1}%", sparsity.sparsity * 100.0);
    println!(
        "  {} total entries, {} non-zero",
        sparsity.total_entries, sparsity.nonzero_entries
    );
    println!();

    let prevalence = profile_prevalence(&counts);
    println!("Prevalence:");
    println!("  Mean: {:.1}%", prevalence.mean_prevalence * 100.0);
    println!("  Rare features (<10%): {}", prevalence.n_rare);
    println!("  Ubiquitous features: {}", prevalence.n_ubiquitous);
    println!();

    let lib_size = profile_library_size(&counts);
    println!("Library sizes:");
    println!("  Mean: {:.0}", lib_size.mean);
    println!("  Range: {} - {}", lib_size.min, lib_size.max);
    println!("  CV: {:.2}", lib_size.cv);
    println!();

    // Run the pipeline
    println!("=== Running LinDA Pipeline ===\n");

    let results = Pipeline::new()
        .name("LinDA-example")
        .filter_prevalence(0.2) // Keep features in ≥20% of samples
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh()
        .run(&counts, &metadata)?;

    println!("Pipeline complete!");
    println!("  Method: {}", results.method);
    println!("  Features tested: {}", results.len());
    println!();

    // Summary statistics
    let summary = results.summary();
    println!("=== Results Summary ===\n");
    println!("{}", summary);

    // Top significant results
    let sorted = results.sorted_by_qvalue();
    println!("=== Top 10 Results (by q-value) ===\n");
    println!(
        "{:<15} {:>10} {:>10} {:>12} {:>12}",
        "Feature", "Estimate", "Std.Err", "p-value", "q-value"
    );
    println!("{}", "-".repeat(65));

    for result in sorted.iter().take(10) {
        println!(
            "{:<15} {:>10.4} {:>10.4} {:>12.2e} {:>12.2e}",
            result.feature_id,
            result.estimate,
            result.std_error,
            result.p_value,
            result.q_value
        );
    }

    // Show significant features
    let significant = results.significant();
    println!("\n=== Significant Features (q < 0.05) ===\n");

    if significant.is_empty() {
        println!("No features significant at q < 0.05");
    } else {
        for result in &significant {
            let direction = if result.estimate > 0.0 {
                "increased"
            } else {
                "decreased"
            };
            println!(
                "  {} - {} in treatment (q = {:.4})",
                result.feature_id, direction, result.q_value
            );
        }
    }

    println!("\n=== Pipeline Configuration (YAML) ===\n");

    let config = Pipeline::new()
        .name("LinDA")
        .filter_prevalence(0.2)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh()
        .to_config(Some("LinDA-style differential abundance analysis"));

    println!("{}", config.to_yaml()?);

    Ok(())
}

/// Create example data with known effects.
fn create_example_data() -> (CountMatrix, Metadata) {
    let n_samples = 30; // 15 control + 15 treatment
    let n_features = 50;

    let mut tri_mat = TriMat::new((n_features, n_samples));
    let mut seed = 12345u64;

    let rand_uniform = |s: &mut u64| -> f64 {
        *s = s.wrapping_mul(1103515245).wrapping_add(12345);
        ((*s >> 16) & 0x7FFF) as f64 / 32768.0
    };

    for feat in 0..n_features {
        for sample in 0..n_samples {
            let is_treatment = sample >= 15;

            // Different feature categories:
            let (base, treatment_effect, prevalence) = match feat {
                // Strongly upregulated in treatment
                0..=4 => (100.0, 2.5, 1.0),
                // Moderately upregulated
                5..=9 => (150.0, 1.5, 1.0),
                // Strongly downregulated in treatment
                10..=14 => (120.0, 0.4, 1.0),
                // No effect
                15..=29 => (200.0, 1.0, 1.0),
                // Rare features (low prevalence)
                30..=39 => (50.0, 1.0, 0.3),
                // Very rare
                40..=49 => (30.0, 1.0, 0.15),
                _ => unreachable!(),
            };

            // Check prevalence
            if rand_uniform(&mut seed) > prevalence {
                continue;
            }

            // Calculate count with noise
            let effect = if is_treatment { treatment_effect } else { 1.0 };
            let noise = 0.7 + 0.6 * rand_uniform(&mut seed);
            let count = (base * effect * noise).round() as u64;

            if count > 0 {
                tri_mat.add_triplet(feat, sample, count);
            }
        }
    }

    let feature_ids: Vec<String> = (0..n_features)
        .map(|i| {
            match i {
                0..=4 => format!("up_strong_{}", i),
                5..=9 => format!("up_moderate_{}", i - 5),
                10..=14 => format!("down_strong_{}", i - 10),
                15..=29 => format!("no_effect_{}", i - 15),
                30..=39 => format!("rare_{}", i - 30),
                _ => format!("very_rare_{}", i - 40),
            }
        })
        .collect();

    let sample_ids: Vec<String> = (0..n_samples).map(|i| format!("S{:02}", i)).collect();

    let counts = CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids.clone()).unwrap();

    // Create metadata
    let metadata = create_metadata(&sample_ids);

    (counts, metadata)
}

fn create_metadata(sample_ids: &[String]) -> Metadata {
    use std::io::Write;

    let mut file = tempfile::NamedTempFile::new().unwrap();
    writeln!(file, "sample_id\tgroup\tage\tbatch").unwrap();

    for (i, sid) in sample_ids.iter().enumerate() {
        let group = if i < 15 { "control" } else { "treatment" };
        let age = 25 + (i % 15) * 2;
        let batch = if i % 3 == 0 { "A" } else { "B" };
        writeln!(file, "{}\t{}\t{}\t{}", sid, group, age, batch).unwrap();
    }
    file.flush().unwrap();

    Metadata::from_tsv(file.path()).unwrap()
}
