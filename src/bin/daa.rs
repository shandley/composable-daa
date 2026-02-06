//! DAA - Differential Abundance Analysis CLI
//!
//! Command-line interface for composable differential abundance analysis.

use clap::{Parser, Subcommand, ValueEnum};
use composable_daa::data::{CountMatrix, Metadata};
use composable_daa::error::Result;
use composable_daa::pipeline::{Pipeline, PipelineConfig};
use composable_daa::benchmark::{
    generate_synthetic, SyntheticConfig,
    fetch_dataset, list_datasets, clear_cache, BenchmarkDataset,
};
use composable_daa::profile::{profile_library_size, profile_prevalence, profile_sparsity, profile_for_llm};
use composable_daa::spike::{
    evaluate_spikes, spike_abundance_with_mode, SpikeMode, SpikeSelection,
    StressConfig, run_stress_test,
    optimize_prevalence_threshold, PrevalenceOptConfig, PrevalenceFilterLogic, OptimizationCriterion,
};
use std::path::PathBuf;

/// CLI-friendly spike mode enum
#[derive(Debug, Clone, Copy, ValueEnum)]
enum CliSpikeMode {
    /// Raw multiplication of counts (increases library size)
    Raw,
    /// Spike then renormalize to original library size
    Compositional,
    /// Model absolute changes then renormalize (most realistic)
    Absolute,
}

impl From<CliSpikeMode> for SpikeMode {
    fn from(mode: CliSpikeMode) -> Self {
        match mode {
            CliSpikeMode::Raw => SpikeMode::Raw,
            CliSpikeMode::Compositional => SpikeMode::Compositional,
            CliSpikeMode::Absolute => SpikeMode::Absolute,
        }
    }
}

/// Composable Differential Abundance Analysis
#[derive(Parser)]
#[command(name = "daa")]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a pipeline from a YAML configuration file
    Run {
        /// Path to pipeline configuration YAML
        #[arg(long)]
        config: PathBuf,

        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Output path for results TSV
        #[arg(short, long)]
        output: PathBuf,
    },

    /// Run analysis with permutation tests (non-parametric, distribution-free)
    Permutation {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Formula for the model (e.g., "~ group")
        #[arg(short, long)]
        formula: String,

        /// Coefficient to test (e.g., "grouptreatment")
        #[arg(short = 't', long)]
        test_coef: String,

        /// Output path for results TSV
        #[arg(short, long)]
        output: PathBuf,

        /// Prevalence threshold (default: 0.1)
        #[arg(long, default_value = "0.1")]
        prevalence: f64,

        /// Pseudocount value (default: 0.5)
        #[arg(long, default_value = "0.5")]
        pseudocount: f64,

        /// Number of permutations (default: 1000)
        #[arg(long, default_value = "1000")]
        n_permutations: usize,
    },

    /// Profile a count matrix
    Profile {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Output format: text, json, or yaml
        #[arg(short, long, default_value = "text")]
        format: String,
    },

    /// Validate a pipeline using spike-in analysis
    Validate {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Group column name
        #[arg(short, long)]
        group: String,

        /// Target group for spiking
        #[arg(short = 't', long)]
        target: String,

        /// Formula for the model
        #[arg(short, long)]
        formula: String,

        /// Coefficient to test
        #[arg(long)]
        test_coef: String,

        /// Number of features to spike (default: 20)
        #[arg(long, default_value = "20")]
        n_spike: usize,

        /// Fold change for spiked features (default: 2.0)
        #[arg(long, default_value = "2.0")]
        fold_change: f64,

        /// Random seed (default: 42)
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Spike mode: how to handle compositional constraints
        /// - raw: Simple multiplication (increases library size)
        /// - compositional: Renormalize to preserve library size
        /// - absolute: Model true absolute changes (most realistic)
        #[arg(long, value_enum, default_value = "compositional")]
        mode: CliSpikeMode,
    },

    /// Generate an example pipeline configuration
    Example {
        /// Output path for the example YAML
        #[arg(short, long, default_value = "pipeline.yaml")]
        output: PathBuf,
    },

    /// Run compositional stress testing to evaluate pipeline robustness
    Stress {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Group column name
        #[arg(short, long)]
        group: String,

        /// Target group for spiking
        #[arg(short = 't', long)]
        target: String,

        /// Formula for the model
        #[arg(short, long)]
        formula: String,

        /// Coefficient to test
        #[arg(long)]
        test_coef: String,

        /// Spike fractions (comma-separated, e.g., "0.01,0.05,0.10,0.25")
        #[arg(long, default_value = "0.01,0.05,0.10,0.25")]
        spike_fractions: String,

        /// Fold changes to test (comma-separated, e.g., "1.5,2.0,3.0,5.0")
        #[arg(long, default_value = "1.5,2.0,3.0,5.0")]
        fold_changes: String,

        /// Spike modes to test (comma-separated: raw,compositional,absolute)
        #[arg(long, default_value = "raw,compositional,absolute")]
        modes: String,

        /// Number of replicates per parameter combination
        #[arg(long, default_value = "5")]
        replicates: usize,

        /// Number of permutations for FDR calibration
        #[arg(long, default_value = "5")]
        permutations: usize,

        /// Random seed
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Output format: text, json, or csv
        #[arg(long, default_value = "text")]
        output_format: String,

        /// Quick mode: use fewer replicates and permutations
        #[arg(long)]
        quick: bool,
    },

    /// Find optimal prevalence threshold via spike-in validation
    OptimizePrevalence {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Group column name
        #[arg(short, long)]
        group: String,

        /// Target group for spiking
        #[arg(short = 't', long)]
        target: String,

        /// Formula for the model
        #[arg(short, long)]
        formula: String,

        /// Coefficient to test
        #[arg(long)]
        test_coef: String,

        /// Thresholds to evaluate (comma-separated, e.g., "0.01,0.05,0.10,0.20,0.30")
        #[arg(long, default_value = "0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30")]
        thresholds: String,

        /// Fold change for spike-ins
        #[arg(long, default_value = "2.0")]
        fold_change: f64,

        /// Number of features to spike per replicate
        #[arg(long, default_value = "20")]
        n_spike: usize,

        /// Number of replicates per threshold
        #[arg(long, default_value = "5")]
        replicates: usize,

        /// Optimization criterion: max-f1, max-sensitivity, max-efficiency, min-fdr
        #[arg(long, default_value = "max-f1")]
        criterion: String,

        /// Filter logic: overall, any-group, all-groups
        #[arg(long, default_value = "overall")]
        filter_logic: String,

        /// Random seed
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Output format: text, json, or csv
        #[arg(long, default_value = "text")]
        output_format: String,

        /// Quick mode: use fewer thresholds and replicates
        #[arg(long)]
        quick: bool,
    },

    /// Generate LLM-friendly data profile for AI-assisted pipeline selection
    ProfileLlm {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Group column name
        #[arg(short, long)]
        group: String,

        /// Output file (default: stdout for piping)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format: yaml (default) or markdown
        #[arg(long, default_value = "yaml")]
        format: String,
    },

    /// Generate synthetic benchmark data with known ground truth
    Generate {
        /// Preset: ideal, typical_16s, sparse_virome, extreme_sparse, group_specific, confounded, small_n
        #[arg(short, long, default_value = "typical_16s")]
        preset: String,

        /// Output directory for generated files
        #[arg(short, long)]
        output: PathBuf,

        /// Random seed for reproducibility
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Number of features (overrides preset)
        #[arg(long)]
        n_features: Option<usize>,

        /// Number of samples per group (overrides preset)
        #[arg(long)]
        n_samples: Option<usize>,

        /// Effect size in log2 fold change (overrides preset)
        #[arg(long)]
        effect_size: Option<f64>,

        /// Number of differential features (overrides preset)
        #[arg(long)]
        n_differential: Option<usize>,
    },

    /// Fetch classic benchmark datasets from Zenodo
    Fetch {
        /// Dataset: hmp_v13, hmp_v35, hmp_subset, hmp_wms, ravel, stammler
        #[arg(short, long)]
        dataset: Option<String>,

        /// Output directory for downloaded files
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// List available datasets
        #[arg(long)]
        list: bool,

        /// Clear the download cache
        #[arg(long)]
        clear_cache: bool,
    },

    /// Recommend and optionally run analysis based on data profile
    Recommend {
        /// Path to count matrix TSV
        #[arg(short = 'c', long)]
        counts: PathBuf,

        /// Path to metadata TSV
        #[arg(short, long)]
        metadata: PathBuf,

        /// Group column name
        #[arg(short, long)]
        group: String,

        /// Target level within group (e.g., "treatment", "disease")
        #[arg(short = 't', long)]
        target: String,

        /// Output file path (results TSV or pipeline YAML)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Run the recommended analysis (not just print recommendation)
        #[arg(long, conflicts_with = "yaml")]
        run: bool,

        /// Output an editable YAML pipeline config instead of running
        #[arg(long, conflicts_with = "run")]
        yaml: bool,

        /// Just print the command/config, no explanation
        #[arg(long)]
        quiet: bool,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Run {
            config,
            counts,
            metadata,
            output,
        } => cmd_run(&config, &counts, &metadata, &output),

        Commands::Permutation {
            counts,
            metadata,
            formula,
            test_coef,
            output,
            prevalence,
            pseudocount,
            n_permutations,
        } => cmd_permutation(
            &counts,
            &metadata,
            &formula,
            &test_coef,
            &output,
            prevalence,
            pseudocount,
            n_permutations,
        ),

        Commands::Profile { counts, format } => cmd_profile(&counts, &format),

        Commands::Validate {
            counts,
            metadata,
            group,
            target,
            formula,
            test_coef,
            n_spike,
            fold_change,
            seed,
            mode,
        } => cmd_validate(
            &counts, &metadata, &group, &target, &formula, &test_coef, n_spike, fold_change, seed, mode,
        ),

        Commands::Example { output } => cmd_example(&output),

        Commands::Stress {
            counts,
            metadata,
            group,
            target,
            formula,
            test_coef,
            spike_fractions,
            fold_changes,
            modes,
            replicates,
            permutations,
            seed,
            output_format,
            quick,
        } => cmd_stress(
            &counts,
            &metadata,
            &group,
            &target,
            &formula,
            &test_coef,
            &spike_fractions,
            &fold_changes,
            &modes,
            replicates,
            permutations,
            seed,
            &output_format,
            quick,
        ),

        Commands::OptimizePrevalence {
            counts,
            metadata,
            group,
            target,
            formula,
            test_coef,
            thresholds,
            fold_change,
            n_spike,
            replicates,
            criterion,
            filter_logic,
            seed,
            output_format,
            quick,
        } => cmd_optimize_prevalence(
            &counts,
            &metadata,
            &group,
            &target,
            &formula,
            &test_coef,
            &thresholds,
            fold_change,
            n_spike,
            replicates,
            &criterion,
            &filter_logic,
            seed,
            &output_format,
            quick,
        ),

        Commands::ProfileLlm {
            counts,
            metadata,
            group,
            output,
            format,
        } => cmd_profile_llm(&counts, &metadata, &group, output.as_ref(), &format),

        Commands::Generate {
            preset,
            output,
            seed,
            n_features,
            n_samples,
            effect_size,
            n_differential,
        } => cmd_generate(&preset, &output, seed, n_features, n_samples, effect_size, n_differential),

        Commands::Fetch {
            dataset,
            output,
            list,
            clear_cache: clear,
        } => cmd_fetch(dataset.as_deref(), output.as_ref(), list, clear),

        Commands::Recommend {
            counts,
            metadata,
            group,
            target,
            output,
            run,
            yaml,
            quiet,
        } => cmd_recommend(&counts, &metadata, &group, &target, output.as_ref(), run, yaml, quiet),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

/// Run a pipeline from configuration
fn cmd_run(
    config_path: &PathBuf,
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    output_path: &PathBuf,
) -> Result<()> {
    eprintln!("Loading pipeline configuration from {:?}...", config_path);
    let config_str = std::fs::read_to_string(config_path)?;
    let config = PipelineConfig::from_yaml(&config_str)?;

    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    eprintln!("Running pipeline '{}'...", config.name);
    let pipeline = Pipeline::from_config(&config);
    let results = pipeline.run(&counts, &metadata)?;

    eprintln!("Writing results to {:?}...", output_path);
    results.to_tsv(output_path)?;

    eprintln!("Done! {} features tested", results.len());
    let n_sig = results.significant().len();
    eprintln!("  {} significant at q < 0.05", n_sig);

    Ok(())
}

/// Truncate a string to max_len characters
fn truncate_str(s: &str, max_len: usize) -> String {
    if s.len() <= max_len {
        s.to_string()
    } else {
        format!("{}...", &s[..max_len-3])
    }
}

/// Print formatted results summary
fn print_results_summary(results: &composable_daa::data::DaResultSet, show_fold_change: bool) {
    let n_sig_05 = results.significant().len();
    let n_sig_10 = results.results.iter().filter(|r| r.q_value < 0.10).count();
    let n_up = results.results.iter().filter(|r| r.q_value < 0.05 && r.estimate > 0.0).count();
    let n_down = results.results.iter().filter(|r| r.q_value < 0.05 && r.estimate < 0.0).count();

    eprintln!();
    eprintln!("=== Results Summary ===");
    eprintln!("Features tested:     {}", results.len());
    eprintln!("Significant q<0.05:  {}", n_sig_05);
    eprintln!("Significant q<0.10:  {}", n_sig_10);
    if n_sig_05 > 0 {
        eprintln!("  Up in target:      {}", n_up);
        eprintln!("  Down in target:    {}", n_down);
    }

    // Print top hits
    let mut sorted = results.results.clone();
    sorted.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap());
    if !sorted.is_empty() {
        eprintln!();
        if show_fold_change {
            eprintln!("Top 5 hits:");
            eprintln!("  {:<20} {:>10} {:>8} {:>10}", "Feature", "Estimate", "FC", "q-value");
            eprintln!("  {}", "-".repeat(52));
            for r in sorted.iter().take(5) {
                let sig_marker = if r.q_value < 0.05 { "*" } else if r.q_value < 0.10 { "." } else { " " };
                let fc = r.estimate.exp();
                eprintln!(
                    "  {:<20} {:>10.3} {:>7.1}x {:>10.4}{}",
                    truncate_str(&r.feature_id, 20), r.estimate, fc, r.q_value, sig_marker
                );
            }
        } else {
            eprintln!("Top 5 hits:");
            eprintln!("  {:<20} {:>10} {:>10}", "Feature", "Estimate", "q-value");
            eprintln!("  {}", "-".repeat(42));
            for r in sorted.iter().take(5) {
                let sig_marker = if r.q_value < 0.05 { "*" } else if r.q_value < 0.10 { "." } else { " " };
                eprintln!(
                    "  {:<20} {:>10.3} {:>10.4}{}",
                    truncate_str(&r.feature_id, 20), r.estimate, r.q_value, sig_marker
                );
            }
        }
        eprintln!();
        eprintln!("  * = q < 0.05, . = q < 0.10");
    }
}
/// Run analysis with permutation tests
fn cmd_permutation(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    formula: &str,
    test_coef: &str,
    output_path: &PathBuf,
    prevalence: f64,
    pseudocount: f64,
    n_permutations: usize,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    eprintln!("Running permutation test analysis...");
    eprintln!("  Formula: {}", formula);
    eprintln!("  Testing: {}", test_coef);
    eprintln!("  Prevalence threshold: {:.1}%", prevalence * 100.0);
    eprintln!("  Permutations: {}", n_permutations);

    let results = Pipeline::new()
        .name("Permutation")
        .filter_prevalence(prevalence)
        .add_pseudocount(pseudocount)
        .normalize_clr()
        .model_lm(formula)
        .test_permutation(test_coef, n_permutations)
        .correct_bh()
        .run(&counts, &metadata)?;

    eprintln!("Writing results to {:?}...", output_path);
    results.to_tsv(output_path)?;

    print_results_summary(&results, false);

    Ok(())
}

/// Profile a count matrix
fn cmd_profile(counts_path: &PathBuf, format: &str) -> Result<()> {
    eprintln!("Loading count matrix...");
    let counts = CountMatrix::from_tsv(counts_path)?;

    let sparsity = profile_sparsity(&counts);
    let prevalence = profile_prevalence(&counts);
    let lib_size = profile_library_size(&counts);

    match format {
        "json" => {
            let profile = serde_json::json!({
                "dimensions": {
                    "n_features": counts.n_features(),
                    "n_samples": counts.n_samples()
                },
                "sparsity": {
                    "total_entries": sparsity.total_entries,
                    "non_zero": sparsity.nonzero_entries,
                    "zero": sparsity.zero_entries,
                    "sparsity": sparsity.sparsity
                },
                "prevalence": {
                    "n_features": prevalence.n_features,
                    "mean": prevalence.mean_prevalence,
                    "median": prevalence.median_prevalence,
                    "n_rare": prevalence.n_rare,
                    "n_ubiquitous": prevalence.n_ubiquitous
                },
                "library_size": {
                    "n_samples": lib_size.n_samples,
                    "mean": lib_size.mean,
                    "median": lib_size.median,
                    "min": lib_size.min,
                    "max": lib_size.max,
                    "cv": lib_size.cv
                }
            });
            println!("{}", serde_json::to_string_pretty(&profile).unwrap());
        }
        "yaml" => {
            let profile = serde_json::json!({
                "dimensions": {
                    "n_features": counts.n_features(),
                    "n_samples": counts.n_samples()
                },
                "sparsity": {
                    "total_entries": sparsity.total_entries,
                    "non_zero": sparsity.nonzero_entries,
                    "zero": sparsity.zero_entries,
                    "sparsity": sparsity.sparsity
                },
                "prevalence": {
                    "n_features": prevalence.n_features,
                    "mean": prevalence.mean_prevalence,
                    "median": prevalence.median_prevalence,
                    "n_rare": prevalence.n_rare,
                    "n_ubiquitous": prevalence.n_ubiquitous
                },
                "library_size": {
                    "n_samples": lib_size.n_samples,
                    "mean": lib_size.mean,
                    "median": lib_size.median,
                    "min": lib_size.min,
                    "max": lib_size.max,
                    "cv": lib_size.cv
                }
            });
            println!("{}", serde_yaml::to_string(&profile).unwrap());
        }
        _ => {
            // Text format
            println!("Data Profile");
            println!("============");
            println!();
            println!("Dimensions:");
            println!("  Features: {}", counts.n_features());
            println!("  Samples:  {}", counts.n_samples());
            println!();
            println!("Sparsity:");
            println!("  Total entries:    {}", sparsity.total_entries);
            println!("  Non-zero entries: {}", sparsity.nonzero_entries);
            println!("  Zero entries:     {}", sparsity.zero_entries);
            println!("  Sparsity:         {:.1}%", sparsity.sparsity * 100.0);
            println!();
            println!("Prevalence:");
            println!("  Mean prevalence:   {:.1}%", prevalence.mean_prevalence * 100.0);
            println!("  Median prevalence: {:.1}%", prevalence.median_prevalence * 100.0);
            println!(
                "  Rare features (<10%):    {} ({:.1}%)",
                prevalence.n_rare,
                prevalence.n_rare as f64 / prevalence.n_features as f64 * 100.0
            );
            println!(
                "  Ubiquitous features (100%): {} ({:.1}%)",
                prevalence.n_ubiquitous,
                prevalence.n_ubiquitous as f64 / prevalence.n_features as f64 * 100.0
            );
            println!();
            println!("Library Size:");
            println!("  Mean:   {:.0}", lib_size.mean);
            println!("  Median: {}", lib_size.median);
            println!("  Min:    {}", lib_size.min);
            println!("  Max:    {}", lib_size.max);
            println!("  CV:     {:.2}", lib_size.cv);
        }
    }

    Ok(())
}

/// Validate a pipeline using spike-in
fn cmd_validate(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    group: &str,
    target: &str,
    formula: &str,
    test_coef: &str,
    n_spike: usize,
    fold_change: f64,
    seed: u64,
    mode: CliSpikeMode,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    let mode_name = match mode {
        CliSpikeMode::Raw => "raw",
        CliSpikeMode::Compositional => "compositional",
        CliSpikeMode::Absolute => "absolute",
    };
    eprintln!("Spiking {} features with {}x fold change ({} mode)...", n_spike, fold_change, mode_name);

    let spiked = spike_abundance_with_mode(
        &counts,
        &metadata,
        group,
        n_spike,
        fold_change,
        target,
        SpikeSelection::Random,
        seed,
        mode.into(),
    )?;

    eprintln!("Running pipeline on spiked data...");
    let results = Pipeline::new()
        .name("validation")
        .filter_prevalence(0.1)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm(formula)
        .test_wald(test_coef)
        .correct_bh()
        .run(&spiked.counts, &metadata)?;

    eprintln!("Evaluating spike detection...");
    let eval = evaluate_spikes(&results, &spiked.spec, 0.05);

    println!("Spike-in Validation Results");
    println!("===========================");
    println!();
    println!("Configuration:");
    println!("  Features spiked: {}", n_spike);
    println!("  Fold change:     {}x", fold_change);
    println!("  Target group:    {}", target);
    println!("  Spike mode:      {}", mode_name);
    println!();

    // Print diagnostics if available
    if let Some(diag) = &spiked.spec.diagnostics {
        println!("Compositional Diagnostics:");
        println!("  Geometric mean ratio:  {:.3}", diag.geometric_mean_ratio);
        println!("  Effective CLR effect:  {:.3}", diag.effective_clr_effect);
        println!("  Nominal log(FC):       {:.3}", fold_change.ln());
        if let Some(warning) = &diag.compositional_warning {
            println!("  Warning: {}", warning);
        }
        println!();
    }

    println!("Detection (at FDR < 5%):");
    println!("  True positives:  {}", eval.true_positives);
    println!("  False positives: {}", eval.false_positives);
    println!("  False negatives: {}", eval.false_negatives);
    println!("  True negatives:  {}", eval.true_negatives);
    println!();
    println!("Metrics:");
    println!("  Sensitivity: {:.1}%", eval.sensitivity * 100.0);
    println!("  Specificity: {:.1}%", eval.specificity * 100.0);
    println!("  Precision:   {:.1}%", eval.precision * 100.0);
    println!("  FDR:         {:.1}%", eval.fdr * 100.0);
    println!("  F1 Score:    {:.3}", eval.f1_score);
    println!();

    if eval.effect_correlation.is_finite() {
        println!("Effect Size Recovery:");
        println!("  Correlation: {:.3}", eval.effect_correlation);
        println!("  Bias:        {:.3}", eval.effect_bias);
        println!("  MAE:         {:.3}", eval.effect_mae);
        println!();
    }

    if eval.is_good() {
        println!("Result: GOOD (sensitivity >= 70%, FDR <= 15%)");
    } else {
        println!("Result: NEEDS IMPROVEMENT");
        if eval.sensitivity < 0.7 {
            println!("  - Low sensitivity: consider larger effect sizes or more samples");
        }
        if eval.fdr > 0.15 {
            println!("  - High FDR: consider stricter filtering or more samples");
        }
    }

    Ok(())
}

/// Generate example pipeline configuration
fn cmd_example(output_path: &PathBuf) -> Result<()> {
    let pipeline = Pipeline::new()
        .name("example-linda")
        .filter_library_size(Some(1000), None)
        .filter_prevalence(0.1)
        .add_pseudocount(0.5)
        .normalize_clr()
        .model_lm("~ group")
        .test_wald("grouptreatment")
        .correct_bh();

    let config = pipeline.to_config(Some(
        "Example LinDA-style pipeline for differential abundance analysis",
    ));
    let yaml = config.to_yaml()?;

    std::fs::write(output_path, &yaml)?;
    eprintln!("Wrote example pipeline to {:?}", output_path);
    eprintln!();
    eprintln!("Contents:");
    println!("{}", yaml);

    Ok(())
}

/// Run compositional stress testing
#[allow(clippy::too_many_arguments)]
fn cmd_stress(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    group: &str,
    target: &str,
    formula: &str,
    test_coef: &str,
    spike_fractions_str: &str,
    fold_changes_str: &str,
    modes_str: &str,
    replicates: usize,
    permutations: usize,
    seed: u64,
    output_format: &str,
    quick: bool,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    // Parse spike fractions
    let spike_fractions: Vec<f64> = spike_fractions_str
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    // Parse fold changes
    let fold_changes: Vec<f64> = fold_changes_str
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    // Parse modes
    let modes: Vec<SpikeMode> = modes_str
        .split(',')
        .filter_map(|s| match s.trim().to_lowercase().as_str() {
            "raw" => Some(SpikeMode::Raw),
            "compositional" => Some(SpikeMode::Compositional),
            "absolute" => Some(SpikeMode::Absolute),
            _ => None,
        })
        .collect();

    // Build config
    let config = if quick {
        StressConfig::quick()
            .with_groups(group, target)
    } else {
        StressConfig::new("stress_test")
            .with_spike_fractions(spike_fractions)
            .with_fold_changes(fold_changes)
            .with_modes(modes)
            .with_groups(group, target)
            .with_replicates(replicates, permutations)
    };

    // Calculate total runs
    let grid_size = config.spike_fractions.len() * config.fold_changes.len() * config.spike_modes.len();
    let total_runs = grid_size * config.n_replicates;

    eprintln!("Running stress test...");
    eprintln!("  Spike fractions: {:?}", config.spike_fractions);
    eprintln!("  Fold changes: {:?}", config.fold_changes);
    eprintln!("  Modes: {:?}", config.spike_modes);
    eprintln!("  Replicates: {}", config.n_replicates);
    eprintln!("  Permutations: {}", config.n_permutations);
    eprintln!("  Total runs: {}", total_runs);
    eprintln!();

    // Create the pipeline closure
    let formula_owned = formula.to_string();
    let test_coef_owned = test_coef.to_string();

    let run_pipeline = move |c: &CountMatrix, m: &Metadata| {
        Pipeline::new()
            .name("stress_pipeline")
            .filter_prevalence(0.1)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm(&formula_owned)
            .test_wald(&test_coef_owned)
            .correct_bh()
            .run(c, m)
    };

    // Run stress test with seed
    let mut config_with_seed = config;
    config_with_seed.seed = seed;

    let summary = run_stress_test(&counts, &metadata, &config_with_seed, run_pipeline)?;

    // Output results
    match output_format {
        "json" => {
            let json = summary.to_json()?;
            println!("{}", json);
        }
        "csv" => {
            let csv = summary.to_csv();
            print!("{}", csv);
        }
        _ => {
            // Text format (default)
            println!("{}", summary);
        }
    }

    Ok(())
}

/// Optimize prevalence threshold via spike-in validation
#[allow(clippy::too_many_arguments)]
fn cmd_optimize_prevalence(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    group: &str,
    target: &str,
    formula: &str,
    test_coef: &str,
    thresholds_str: &str,
    fold_change: f64,
    n_spike: usize,
    replicates: usize,
    criterion_str: &str,
    filter_logic_str: &str,
    seed: u64,
    output_format: &str,
    quick: bool,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    // Parse thresholds
    let thresholds: Vec<f64> = thresholds_str
        .split(',')
        .filter_map(|s| s.trim().parse().ok())
        .collect();

    // Parse criterion
    let criterion = match criterion_str.to_lowercase().as_str() {
        "max-f1" => OptimizationCriterion::MaxF1,
        "max-sensitivity" => OptimizationCriterion::MaxSensitivityAtFdr(0.10),
        "max-efficiency" => OptimizationCriterion::MaxEfficiency,
        "min-fdr" => OptimizationCriterion::MinFdrAtSensitivity(0.70),
        _ => OptimizationCriterion::MaxF1,
    };

    // Parse filter logic
    let filter_logic = match filter_logic_str.to_lowercase().as_str() {
        "any-group" => PrevalenceFilterLogic::AnyGroup,
        "all-groups" => PrevalenceFilterLogic::AllGroups,
        _ => PrevalenceFilterLogic::Overall,
    };

    // Build config
    let config = if quick {
        PrevalenceOptConfig::quick()
    } else {
        PrevalenceOptConfig {
            thresholds,
            fold_change,
            n_spike,
            n_replicates: replicates,
            group_column: group.to_string(),
            target_group: target.to_string(),
            fdr_threshold: 0.05,
            selection_criterion: criterion,
            filter_logic,
            spike_mode: SpikeMode::Compositional,
            seed,
        }
    };

    let total_runs = config.thresholds.len() * config.n_replicates;

    eprintln!("Optimizing prevalence threshold...");
    eprintln!("  Thresholds: {:?}", config.thresholds);
    eprintln!("  Fold change: {}x", config.fold_change);
    eprintln!("  Features per spike: {}", config.n_spike);
    eprintln!("  Replicates: {}", config.n_replicates);
    eprintln!("  Total evaluations: {}", total_runs);
    eprintln!();

    // Create the pipeline closure
    let formula_owned = formula.to_string();
    let test_coef_owned = test_coef.to_string();

    let run_pipeline = move |c: &CountMatrix, m: &Metadata| {
        Pipeline::new()
            .name("optimize_prevalence")
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm(&formula_owned)
            .test_wald(&test_coef_owned)
            .correct_bh()
            .run(c, m)
    };

    let result = optimize_prevalence_threshold(&counts, &metadata, &config, run_pipeline)?;

    // Output results
    match output_format {
        "json" => {
            let json = serde_json::to_string_pretty(&result)?;
            println!("{}", json);
        }
        "csv" => {
            let csv = result.to_csv();
            print!("{}", csv);
        }
        _ => {
            // Text format (default)
            println!("{}", result);
        }
    }

    Ok(())
}

/// Generate LLM-friendly data profile
fn cmd_profile_llm(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    group: &str,
    output_path: Option<&PathBuf>,
    format: &str,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    eprintln!("Generating LLM-friendly profile...");
    let profile = profile_for_llm(&counts, &metadata, group)?;

    // Generate output in requested format
    let output = match format {
        "markdown" | "md" => profile.to_markdown(),
        _ => profile.to_yaml()?,
    };

    // Write to file or stdout
    match output_path {
        Some(path) => {
            std::fs::write(path, &output)?;
            eprintln!("Profile written to {:?}", path);
        }
        None => {
            println!("{}", output);
        }
    }

    Ok(())
}

/// Generate synthetic benchmark data
fn cmd_generate(
    preset: &str,
    output_dir: &PathBuf,
    seed: u64,
    n_features: Option<usize>,
    n_samples: Option<usize>,
    effect_size: Option<f64>,
    n_differential: Option<usize>,
) -> Result<()> {
    // Select preset config
    let mut config = match preset.to_lowercase().as_str() {
        "ideal" => SyntheticConfig::ideal(),
        "typical_16s" | "typical" => SyntheticConfig::typical_16s(),
        "sparse_virome" | "virome" | "sparse" => SyntheticConfig::sparse_virome(),
        "extreme_sparse" | "extreme" => SyntheticConfig::extreme_sparse(),
        "group_specific" | "groupspec" => SyntheticConfig::group_specific(),
        "confounded" => SyntheticConfig::confounded(),
        "small_n" | "small" => SyntheticConfig::small_n(),
        _ => {
            eprintln!("Unknown preset '{}'. Available: ideal, typical_16s, sparse_virome, extreme_sparse, group_specific, confounded, small_n", preset);
            return Err(composable_daa::error::DaaError::InvalidParameter(
                format!("Unknown preset: {}", preset)
            ));
        }
    };

    // Apply overrides
    config.seed = seed;
    if let Some(n) = n_features {
        config.n_features = n;
    }
    if let Some(n) = n_samples {
        config.n_samples_per_group = n;
    }
    if let Some(es) = effect_size {
        config.effect_size = es;
    }
    if let Some(n) = n_differential {
        config.n_differential = n;
    }

    eprintln!("Generating synthetic data...");
    eprintln!("  Preset: {}", config.name);
    eprintln!("  Features: {}", config.n_features);
    eprintln!("  Samples per group: {}", config.n_samples_per_group);
    eprintln!("  Target sparsity: {:.0}%", config.sparsity * 100.0);
    eprintln!("  Differential features: {}", config.n_differential);
    eprintln!("  Effect size: {:.1} log2FC", config.effect_size);
    eprintln!("  Seed: {}", config.seed);

    let data = generate_synthetic(&config)?;

    // Calculate actual sparsity
    let mut zeros = 0;
    let total = data.counts.n_features() * data.counts.n_samples();
    for f in 0..data.counts.n_features() {
        for s in 0..data.counts.n_samples() {
            if data.counts.get(f, s) == 0 {
                zeros += 1;
            }
        }
    }
    let actual_sparsity = zeros as f64 / total as f64;

    eprintln!();
    eprintln!("Generated:");
    eprintln!("  Actual sparsity: {:.1}%", actual_sparsity * 100.0);
    eprintln!("  Differential features: {}", data.ground_truth.differential_features.len());

    // Write files
    data.write_to_dir(output_dir)?;

    eprintln!();
    eprintln!("Files written to {:?}:", output_dir);
    eprintln!("  counts.tsv       - Count matrix");
    eprintln!("  metadata.tsv     - Sample metadata");
    eprintln!("  ground_truth.tsv - True differential features");
    eprintln!("  config.yaml      - Generation config");

    Ok(())
}

/// Fetch classic benchmark datasets from Zenodo
fn cmd_fetch(
    dataset: Option<&str>,
    output: Option<&PathBuf>,
    list: bool,
    clear: bool,
) -> Result<()> {
    // Handle cache clearing
    if clear {
        eprintln!("Clearing dataset cache...");
        clear_cache(None)?;
        eprintln!("Cache cleared.");
        return Ok(());
    }

    // Handle listing
    if list || dataset.is_none() {
        eprintln!("Available benchmark datasets:\n");
        eprintln!("  {:<15} {:<50} {}", "ID", "Description", "Cached");
        eprintln!("  {}", "-".repeat(75));

        for info in list_datasets(None) {
            let cached = if info.cached { "yes" } else { "no" };
            eprintln!(
                "  {:<15} {:<50} {}",
                match info.dataset {
                    BenchmarkDataset::HmpGingivalV13 => "hmp_v13",
                    BenchmarkDataset::HmpGingivalV35 => "hmp_v35",
                    BenchmarkDataset::HmpGingivalV35Subset => "hmp_subset",
                    BenchmarkDataset::HmpGingivalWms => "hmp_wms",
                    BenchmarkDataset::RavelBv => "ravel",
                    BenchmarkDataset::StammlerSpikein => "stammler",
                },
                info.dataset.description(),
                cached
            );
        }

        eprintln!();
        eprintln!("Usage: daa fetch -d <dataset> [-o <output_dir>]");
        eprintln!();
        eprintln!("Note: stammler has experimental spike-in ground truth!");
        return Ok(());
    }

    // Fetch the requested dataset
    let ds_name = dataset.unwrap();
    let ds = BenchmarkDataset::from_str(ds_name).ok_or_else(|| {
        composable_daa::error::DaaError::InvalidParameter(format!(
            "Unknown dataset: {}. Run 'daa fetch --list' to see available datasets.",
            ds_name
        ))
    })?;

    eprintln!("Fetching dataset: {}", ds.name());
    let fetched = fetch_dataset(ds, None)?;

    eprintln!("Loaded: {} features x {} samples",
        fetched.counts.n_features(),
        fetched.counts.n_samples()
    );

    // Copy to output directory if specified
    if let Some(out_dir) = output {
        std::fs::create_dir_all(out_dir)?;

        let counts_dest = out_dir.join("counts.tsv");
        let metadata_dest = out_dir.join("metadata.tsv");

        fetched.counts.to_tsv(&counts_dest)?;
        // Copy metadata from cache
        let metadata_src = fetched.cache_dir.join(format!("{}_sample_metadata.tsv", ds.name()));
        std::fs::copy(&metadata_src, &metadata_dest)?;

        eprintln!();
        eprintln!("Files copied to {:?}:", out_dir);
        eprintln!("  counts.tsv   - Count matrix");
        eprintln!("  metadata.tsv - Sample metadata");
    } else {
        eprintln!();
        eprintln!("Files cached at: {:?}", fetched.cache_dir);
        eprintln!("Use -o <dir> to copy to a specific location.");
    }

    if ds.has_ground_truth() {
        eprintln!();
        eprintln!("This dataset has experimental spike-in controls!");
    }

    Ok(())
}

/// Recommend analysis method based on data profile
fn cmd_recommend(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    group: &str,
    target: &str,
    output_path: Option<&PathBuf>,
    run: bool,
    yaml: bool,
    quiet: bool,
) -> Result<()> {
    // Load data
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    if !quiet {
        eprintln!("Analyzing data profile...");
    }

    // Get profile
    let profile = profile_for_llm(&counts, &metadata, group)?;

    // Extract key metrics from profile
    let sparsity = profile.sparsity.overall;
    let n_features = profile.data_summary.features;
    let n_samples = profile.data_summary.samples;

    // Get group sizes
    let group_sizes: Vec<usize> = profile.data_summary.groups.sizes.values().cloned().collect();
    let min_group_size = group_sizes.iter().min().copied().unwrap_or(0);

    // Detect potential covariates and study design from metadata columns
    let all_columns: Vec<&str> = metadata.column_names().iter().map(|s| s.as_str()).collect();

    // Patterns for different column types
    let subject_patterns = ["subject", "patient", "individual", "participant", "person", "donor"];
    let time_patterns = ["time", "timepoint", "visit", "day", "week", "month", "year", "tp"];
    let batch_patterns = ["batch", "run", "plate", "lane", "site", "location", "center", "sequencing"];
    let continuous_patterns = ["age", "bmi", "weight", "height"];

    // Detect subject column (for repeated measures / longitudinal)
    let subject_column: Option<&str> = all_columns
        .iter()
        .find(|col| {
            let col_lower = col.to_lowercase();
            col_lower != group.to_lowercase() &&
            subject_patterns.iter().any(|p| col_lower.contains(p))
        })
        .copied();

    // Detect time column (for longitudinal)
    let time_column: Option<&str> = all_columns
        .iter()
        .find(|col| {
            let col_lower = col.to_lowercase();
            col_lower != group.to_lowercase() &&
            time_patterns.iter().any(|p| col_lower.contains(p))
        })
        .copied();

    // Detect batch columns
    let batch_columns: Vec<&str> = all_columns
        .iter()
        .filter(|col| {
            let col_lower = col.to_lowercase();
            col_lower != group.to_lowercase() &&
            batch_patterns.iter().any(|p| col_lower.contains(p))
        })
        .copied()
        .collect();

    // Detect continuous covariates
    let continuous_covariates: Vec<&str> = all_columns
        .iter()
        .filter(|col| {
            let col_lower = col.to_lowercase();
            col_lower != group.to_lowercase() &&
            continuous_patterns.iter().any(|p| col_lower.contains(p))
        })
        .copied()
        .collect();

    // Determine study design
    let is_longitudinal = subject_column.is_some() && time_column.is_some();
    let is_repeated_measures = subject_column.is_some() && time_column.is_none();

    // Collect all detected covariates (excluding subject/time which are handled specially)
    let detected_covariates: Vec<&str> = batch_columns
        .iter()
        .chain(continuous_covariates.iter())
        .copied()
        .collect();

    // Determine recommended method based on sparsity (for cross-sectional)
    let (base_method, base_threshold, base_rationale) = if sparsity > 0.70 {
        ("hurdle", "0.05", "High sparsity (>70%) indicates structural zeros - Hurdle model separates presence from abundance")
    } else if sparsity > 0.50 {
        ("zinb", "0.05", "Moderate-high sparsity (50-70%) suggests zero-inflation - ZINB models excess zeros")
    } else if sparsity > 0.30 {
        ("zinb", "0.05", "Moderate sparsity (30-50%) - ZINB or NB work well")
    } else {
        ("linda", "0.10", "Low sparsity (<30%) - standard CLR + linear model works well. Use q<0.10 due to CLR attenuation")
    };

    // Override for longitudinal/repeated measures designs
    let (method, threshold, rationale, formula, study_design): (&str, &str, String, String, &str) = if is_longitudinal {
        let subj = subject_column.unwrap();
        let time = time_column.unwrap();
        let formula = format!("~ {} + {} + (1 | {})", group, time, subj);
        (
            "lmm",
            "0.05",
            format!("Longitudinal design detected ({} + {}) - LMM with random intercept per subject", time, subj),
            formula,
            "longitudinal",
        )
    } else if is_repeated_measures {
        let subj = subject_column.unwrap();
        let formula = format!("~ {} + (1 | {})", group, subj);
        (
            "lmm",
            "0.05",
            format!("Repeated measures detected ({}) - LMM with random intercept per subject", subj),
            formula,
            "repeated_measures",
        )
    } else {
        (
            base_method,
            base_threshold,
            base_rationale.to_string(),
            format!("~ {}", group),
            "cross_sectional",
        )
    };

    // Build the coefficient name and track if we need to flip interpretation
    // Statistical models use alphabetically-first level as reference by default.
    // The coefficient will be created for the non-reference level(s).
    let levels = metadata.levels(group)?;
    if levels.is_empty() {
        return Err(composable_daa::error::DaaError::InvalidParameter(
            format!("No levels found in group column '{}'", group)
        ));
    }
    if !levels.contains(&target.to_string()) {
        return Err(composable_daa::error::DaaError::InvalidParameter(
            format!("Target '{}' not found in group column '{}'. Available levels: {:?}",
                    target, group, levels)
        ));
    }
    if levels.len() != 2 {
        return Err(composable_daa::error::DaaError::InvalidParameter(
            format!("The recommend command currently only supports two-level comparisons. \
                    Group '{}' has {} levels: {:?}",
                    group, levels.len(), levels)
        ));
    }

    // The model will use the alphabetically-first level as reference (already sorted by levels())
    let reference_level = &levels[0];
    let other_level = &levels[1];

    // Determine which level should be the actual reference based on user's target
    // If user wants target as the comparison group, we want: target - reference
    // But model creates: other_level - reference_level (alphabetically)
    let (actual_reference, flip_sign) = if target == reference_level {
        // User wants: CELIAC - CONTROL (but CELIAC is alphabetically first, so it IS the reference)
        // Model creates coefficient for: CONTROL (which means CONTROL - CELIAC)
        // We need to flip the sign to get CELIAC - CONTROL
        (other_level.as_str(), true)
    } else {
        // User wants: CONTROL - CELIAC (and CELIAC is alphabetically first reference)
        // Model creates coefficient for: CONTROL (which means CONTROL - CELIAC)
        // This is what we want, no flip needed
        (reference_level.as_str(), false)
    };

    // The coefficient name will always be: group + non_reference_level
    let test_coef = format!("{}{}", group, other_level);

    // Handle YAML output mode
    if yaml {
        let yaml_file = output_path
            .map(|p| p.clone())
            .unwrap_or_else(|| PathBuf::from(format!("pipeline_{}.yaml", method)));

        let yaml_content = generate_pipeline_yaml(
            method,
            &formula,
            &test_coef,
            sparsity,
            &detected_covariates,
            group,
            min_group_size,
            &profile.library_size.by_group.as_ref().map(|g| g.ratio),
            study_design,
            &rationale,
            subject_column,
            time_column,
        );

        if quiet {
            // Just print the YAML to stdout
            println!("{}", yaml_content);
        } else {
            // Write to file with explanation
            std::fs::write(&yaml_file, &yaml_content)?;
            println!("Pipeline config written to: {}", yaml_file.display());
            println!();
            println!("To run this pipeline:");
            println!("  daa run -c {} -m {} --config {} -o results.tsv",
                     counts_path.display(), metadata_path.display(), yaml_file.display());
            println!();
            if is_longitudinal || is_repeated_measures {
                println!("Study design: {}", if is_longitudinal { "Longitudinal" } else { "Repeated measures" });
                if let Some(subj) = subject_column {
                    println!("  Subject column: {}", subj);
                }
                if let Some(time) = time_column {
                    println!("  Time column: {}", time);
                }
                println!();
            }
            if !detected_covariates.is_empty() {
                println!("Detected potential covariates: {}", detected_covariates.join(", "));
                println!("Edit the YAML to include them in the formula if needed.");
            }
        }
        return Ok(());
    }

    // Determine output file path for results
    let output_file = output_path
        .map(|p| p.clone())
        .unwrap_or_else(|| PathBuf::from(format!("results_{}.tsv", method)));

    // Build the command for display
    let cmd = format!(
        "daa recommend -c {} -m {} -g {} -t {} -o {} --run",
        counts_path.display(),
        metadata_path.display(),
        group,
        target,
        output_file.display()
    );

    if quiet && !run {
        // Just print the command
        println!("{}", cmd);
        return Ok(());
    }

    if !quiet {
        println!();
        println!("=== Data Summary ===");
        println!("Features:        {}", n_features);
        println!("Samples:         {}", n_samples);
        println!("Sparsity:        {:.1}%", sparsity * 100.0);
        println!("Min group size:  {}", min_group_size);
        if is_longitudinal || is_repeated_measures {
            println!("Study design:    {}", if is_longitudinal { "Longitudinal" } else { "Repeated measures" });
            if let Some(subj) = subject_column {
                println!("  Subject column: {}", subj);
            }
            if let Some(time) = time_column {
                println!("  Time column:    {}", time);
            }
        }
        println!();
        println!("=== Recommendation ===");
        println!("Method:    {}", method.to_uppercase());
        println!("Threshold: q < {}", threshold);
        println!();
        println!("Rationale: {}", rationale);
        println!();

        // Warnings
        if min_group_size < 10 {
            println!("WARNING: Very small group size (n={}). Only large effects will be detectable.", min_group_size);
            println!();
        } else if min_group_size < 20 {
            println!("NOTE: Limited power with n={}. Effects >4x fold change should be detectable.", min_group_size);
            println!();
        }

        // Check library size imbalance
        if let Some(ratio) = profile.library_size.by_group.as_ref().map(|g| g.ratio) {
            if ratio > 2.0 {
                println!("WARNING: Library size imbalance ({:.1}x). Results may be confounded.", ratio);
                println!();
            }
        }
    }

    // If run flag is set, execute the analysis
    if run {
        if !quiet {
            println!("=== Running {} Analysis ===", method.to_uppercase());
            eprintln!("  Formula: {}", formula);
            eprintln!("  Testing: {}", test_coef);
            eprintln!("  Comparison: {} vs {} (reference)", target, actual_reference);
            if flip_sign {
                eprintln!("  Note: Signs flipped to interpret as {} - {}", target, actual_reference);
            }
        }

        // Build and run the appropriate pipeline
        let mut results = match method {
            "hurdle" => {
                Pipeline::new()
                    .name("Hurdle")
                    .filter_prevalence(0.1)
                    .model_hurdle(&formula)
                    .test_wald(&test_coef)
                    .correct_bh()
                    .run(&counts, &metadata)?
            }
            "zinb" => {
                Pipeline::new()
                    .name("ZINB")
                    .filter_prevalence(0.1)
                    .model_zinb(&formula)
                    .test_wald(&test_coef)
                    .correct_bh()
                    .run(&counts, &metadata)?
            }
            "lmm" => {
                Pipeline::new()
                    .name("LMM")
                    .filter_prevalence(0.1)
                    .add_pseudocount(0.5)
                    .normalize_clr()
                    .model_lmm(&formula)
                    .test_wald(&test_coef)
                    .correct_bh()
                    .run(&counts, &metadata)?
            }
            "linda" | _ => {
                Pipeline::new()
                    .name("LinDA")
                    .filter_prevalence(0.1)
                    .add_pseudocount(0.5)
                    .normalize_clr()
                    .model_lm(&formula)
                    .test_wald(&test_coef)
                    .correct_bh()
                    .run(&counts, &metadata)?
            }
        };

        // Flip signs if needed to match user's interpretation
        // User specified target as -t, so results should be interpreted as "target - reference"
        // But if target is alphabetically first, the model uses it as reference
        // So we flip signs to get the correct interpretation
        if flip_sign {
            for result in &mut results.results {
                result.estimate = -result.estimate;
                result.statistic = -result.statistic;
            }
        }

        if !quiet {
            eprintln!("Writing results to {:?}...", output_file);
        }
        results.to_tsv(&output_file)?;

        // Print summary
        let is_linda = method == "linda";
        print_results_summary(&results, !is_linda);

        if !quiet {
            println!();
            println!("=== Interpretation Guide ===");
            if method == "linda" || method == "lmm" {
                println!("  - Use q < {} threshold", threshold);
                println!("  - Effect sizes are in CLR units; multiply by ~4 for approximate log2FC");
                if method == "lmm" {
                    println!("  - Random effects account for within-subject correlation");
                }
            } else {
                println!("  - Use q < 0.05 threshold");
                println!("  - Effect sizes are in log scale; exp(estimate) = fold change");
            }
        }
    } else {
        // Just print the recommendation
        if !quiet {
            println!("=== Ready-to-run Command ===");
            println!("{}", cmd);
            println!();
            println!("After running, interpret results:");
            if method == "linda" || method == "lmm" {
                println!("  - Use q < {} threshold", threshold);
                println!("  - Effect sizes are in CLR units; multiply by ~4 for approximate log2FC");
            } else {
                println!("  - Use q < 0.05 threshold");
                println!("  - Effect sizes are in log scale; exp(estimate) = fold change");
            }
        }
    }

    Ok(())
}

/// Generate a commented YAML pipeline configuration
fn generate_pipeline_yaml(
    method: &str,
    formula: &str,
    test_coef: &str,
    sparsity: f64,
    detected_covariates: &[&str],
    group: &str,
    min_group_size: usize,
    library_ratio: &Option<f64>,
    study_design: &str,
    rationale: &str,
    subject_column: Option<&str>,
    time_column: Option<&str>,
) -> String {
    let mut yaml = String::new();

    // Header comments
    yaml.push_str("# Pipeline Configuration\n");
    yaml.push_str("# Generated by: daa recommend\n");

    // Method and study design info
    if study_design == "longitudinal" || study_design == "repeated_measures" {
        yaml.push_str(&format!("# Method: {} (Linear Mixed Model)\n", method.to_uppercase()));
        yaml.push_str(&format!("# Study design: {}\n", study_design.replace("_", " ")));
        if let Some(subj) = subject_column {
            yaml.push_str(&format!("#   Subject column: {}\n", subj));
        }
        if let Some(time) = time_column {
            yaml.push_str(&format!("#   Time column: {}\n", time));
        }
    } else {
        yaml.push_str(&format!("# Method: {} (selected based on {:.1}% sparsity)\n", method.to_uppercase(), sparsity * 100.0));
    }
    yaml.push_str("#\n");

    // Method selection rationale
    yaml.push_str(&format!("# Rationale: {}\n", rationale));

    // Threshold guidance
    let threshold = if method == "linda" { "0.10" } else { "0.05" };
    yaml.push_str(&format!("# Recommended threshold: q < {}\n", threshold));
    if method == "linda" {
        yaml.push_str("#   (LinDA uses q<0.10 due to CLR effect attenuation)\n");
    }
    yaml.push_str("#\n");

    // Warnings
    if min_group_size < 10 {
        yaml.push_str(&format!("# WARNING: Small group size (n={}). Only large effects detectable.\n", min_group_size));
    } else if min_group_size < 20 {
        yaml.push_str(&format!("# NOTE: Limited power (n={}). Effects >4x fold change detectable.\n", min_group_size));
    }

    if let Some(ratio) = library_ratio {
        if *ratio > 2.0 {
            yaml.push_str(&format!("# WARNING: Library size imbalance ({:.1}x). Consider adding as covariate.\n", ratio));
        }
    }

    // Detected covariates
    if !detected_covariates.is_empty() {
        yaml.push_str("#\n");
        yaml.push_str("# Detected potential covariates (not included - edit formula to add):\n");
        for cov in detected_covariates {
            yaml.push_str(&format!("#   - {}\n", cov));
        }
        yaml.push_str(&format!("#   Example: formula: \"~ {} + {}\"\n", group, detected_covariates.join(" + ")));
    }

    yaml.push_str("#\n");
    yaml.push_str("# Edit this file to customize, then run with:\n");
    yaml.push_str("#   daa run -c counts.tsv -m metadata.tsv --config this_file.yaml -o results.tsv\n");
    yaml.push_str("\n");

    // Pipeline definition
    let method_upper = match method {
        "hurdle" => "Hurdle",
        "zinb" => "ZINB",
        "linda" => "LinDA",
        "lmm" => "LMM",
        _ => method,
    };

    let description = if study_design == "longitudinal" {
        format!("Auto-generated {} pipeline for longitudinal data", method_upper)
    } else if study_design == "repeated_measures" {
        format!("Auto-generated {} pipeline for repeated measures", method_upper)
    } else {
        format!("Auto-generated {} pipeline", method_upper)
    };

    yaml.push_str(&format!("name: {}\n", method_upper));
    yaml.push_str(&format!("description: {}\n", description));
    yaml.push_str("steps:\n");

    // Filter step
    yaml.push_str("- !FilterPrevalence\n");
    yaml.push_str("  threshold: 0.1\n");

    // Method-specific steps
    match method {
        "hurdle" => {
            yaml.push_str("- !ModelHurdle\n");
            yaml.push_str(&format!("  formula: \"{}\"\n", formula));
        }
        "zinb" => {
            yaml.push_str("- !ModelZINB\n");
            yaml.push_str(&format!("  formula: \"{}\"\n", formula));
        }
        "lmm" => {
            // LMM requires CLR transformation
            yaml.push_str("- !AddPseudocount\n");
            yaml.push_str("  value: 0.5\n");
            yaml.push_str("- NormalizeCLR\n");
            yaml.push_str("- !ModelLMM\n");
            yaml.push_str(&format!("  formula: \"{}\"\n", formula));
        }
        "linda" | _ => {
            yaml.push_str("- !AddPseudocount\n");
            yaml.push_str("  value: 0.5\n");
            yaml.push_str("- NormalizeCLR\n");
            yaml.push_str("- !ModelLM\n");
            yaml.push_str(&format!("  formula: \"{}\"\n", formula));
        }
    }

    // Test and correct
    yaml.push_str("- !TestWald\n");
    yaml.push_str(&format!("  coefficient: {}\n", test_coef));
    yaml.push_str("- CorrectBH\n");

    yaml
}
