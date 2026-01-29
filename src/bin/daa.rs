//! DAA - Differential Abundance Analysis CLI
//!
//! Command-line interface for composable differential abundance analysis.

use clap::{Parser, Subcommand, ValueEnum};
use composable_daa::data::{CountMatrix, Metadata};
use composable_daa::error::Result;
use composable_daa::pipeline::{Pipeline, PipelineConfig};
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
        #[arg(short, long)]
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

    /// Run a quick LinDA-style analysis
    Linda {
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

        Commands::Linda {
            counts,
            metadata,
            formula,
            test_coef,
            output,
            prevalence,
            pseudocount,
        } => cmd_linda(
            &counts,
            &metadata,
            &formula,
            &test_coef,
            &output,
            prevalence,
            pseudocount,
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

/// Run a quick LinDA analysis
fn cmd_linda(
    counts_path: &PathBuf,
    metadata_path: &PathBuf,
    formula: &str,
    test_coef: &str,
    output_path: &PathBuf,
    prevalence: f64,
    pseudocount: f64,
) -> Result<()> {
    eprintln!("Loading data...");
    let counts = CountMatrix::from_tsv(counts_path)?;
    let metadata = Metadata::from_tsv(metadata_path)?;

    eprintln!(
        "Loaded {} features x {} samples",
        counts.n_features(),
        counts.n_samples()
    );

    eprintln!("Running LinDA analysis...");
    eprintln!("  Formula: {}", formula);
    eprintln!("  Testing: {}", test_coef);
    eprintln!("  Prevalence threshold: {:.1}%", prevalence * 100.0);

    let results = Pipeline::new()
        .name("LinDA")
        .filter_prevalence(prevalence)
        .add_pseudocount(pseudocount)
        .normalize_clr()
        .model_lm(formula)
        .test_wald(test_coef)
        .correct_bh()
        .run(&counts, &metadata)?;

    eprintln!("Writing results to {:?}...", output_path);
    results.to_tsv(output_path)?;

    eprintln!("Done! {} features tested", results.len());
    let n_sig = results.significant().len();
    eprintln!("  {} significant at q < 0.05", n_sig);

    // Print top hits
    let mut sorted = results.results.clone();
    sorted.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap());
    if !sorted.is_empty() {
        eprintln!("\nTop 5 hits:");
        for r in sorted.iter().take(5) {
            eprintln!(
                "  {}: estimate={:.3}, q={:.4}",
                r.feature_id, r.estimate, r.q_value
            );
        }
    }

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
