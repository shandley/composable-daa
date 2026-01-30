# Composable Differential Abundance Analysis (DAA) Primitives Library

## Project Overview

A Rust library providing composable statistical primitives for differential abundance analysis of sparse count data (microbiome/virome). The core insight is that all existing DAA methods (DESeq2, LinDA, ANCOM-BC, ALDEx2, corncob, etc.) are simply different compositions of fundamental statistical operations. This library exposes those primitives directly, allowing data-driven pipeline construction and empirical validation via spike-in analysis.

## Core Goals

The library addresses four key problems in DAA:

1. **Method selection is arbitrary** - Researchers pick methods without knowing which works best for their data. We solve this with **empirical spike-in validation** that tests methods on the user's actual data structure.

2. **Prevalence filtering is a shot in the dark** - Thresholds like 0.1 or 0.05 are chosen arbitrarily. We solve this with **prevalence threshold optimization** via spike-in, finding the threshold that maximizes sensitivity while controlling FDR.

3. **Group-specific prevalence is ignored** - Healthy vs. disease groups may have different sparsity patterns. We support **group-aware prevalence filtering** with different thresholds per group.

4. **Pipeline selection requires expertise** - Choosing the right combination of primitives is complex. We enable **AI-assisted pipeline design** through structured data profiling that LLMs can interpret and use to recommend configurations.

## Design Philosophy

1. **Primitives over methods**: Expose fundamental operations, not opinionated workflows
2. **Sparse-first**: Virome/microbiome data is 70-90% zeros; design around this from the start
3. **Composable**: Primitives have consistent interfaces and can be chained arbitrarily
4. **Validatable**: Built-in spike-in framework to empirically evaluate pipeline performance on user's actual data structure
5. **Prevalence-aware**: Groupwise prevalence analysis as a first-class concept, not an afterthought
6. **Empirical over theoretical**: Validate on user's data, don't trust assumptions

---

## Implementation Status

### Fully Implemented (254 tests passing)

| Category | Components |
|----------|------------|
| **Data Structures** | CountMatrix, Metadata, Formula, DesignMatrix, DAResult, TransformedMatrix |
| **Profiling** | sparsity, prevalence (by group), library_size, LLM-friendly output |
| **Filtering** | prevalence (overall, groupwise, differential), abundance, library_size, stratified |
| **Zero Handling** | pseudocount (fixed, adaptive) |
| **Normalization** | CLR, TSS, TMM, spike-in normalization |
| **Models** | LM, NB-GLM, ZINB, model comparison (AIC/BIC), effect size shrinkage |
| **Testing** | Wald (LM, NB, ZINB), LRT (NB, ZINB), permutation tests |
| **Correction** | Benjamini-Hochberg |
| **Spike-in** | abundance spikes (3 modes), presence spikes, evaluation, validation, stress testing, threshold optimization |
| **Benchmarking** | synthetic data generation, classic dataset fetcher (Zenodo) |
| **Pipeline** | composition, YAML serialization, stratified presets, CLI |

### Completed (Originally Next Priorities)

1. **Prevalence threshold optimization** - Implemented via spike-in validation with multiple optimization criteria
2. **LLM-friendly profiling** - Structured YAML/markdown output for AI-assisted pipeline design
3. **TMM normalization** - edgeR-style trimmed mean of M-values
4. **Permutation tests** - Distribution-free hypothesis testing

### Future Work

- Linear mixed models (longitudinal data)
- Beta-binomial models
- Hurdle models
- ALR normalization with reference selection
- CSS normalization (metagenomeSeq-style)

---

## Core Data Structures

### CountMatrix

The fundamental input: taxa (rows) × samples (columns) of integer counts.

```rust
pub struct CountMatrix {
    // Sparse representation - CSR format recommended
    data: SparseMatrix<u64>,
    taxa_ids: Vec<String>,
    sample_ids: Vec<String>,
}
```

Requirements:
- Efficient row-wise and column-wise iteration
- Zero-copy views for subsetting
- Serialization support (for caching intermediate results)

### Metadata

Sample metadata with support for categorical and continuous variables.

```rust
pub struct Metadata {
    sample_ids: Vec<String>,
    variables: HashMap<String, Variable>,
}

pub enum Variable {
    Categorical(Vec<String>),      // group labels
    Continuous(Vec<f64>),          // numeric covariates
    Ordinal(Vec<i32>),             // ordered categories
}
```

### Formula

Model specification following R-style formula syntax.

```rust
pub struct Formula {
    response: Option<String>,      // usually None for DAA (taxa are implicit response)
    fixed_effects: Vec<Term>,
    random_effects: Vec<Term>,     // for mixed models
}

pub enum Term {
    Main(String),                  // single variable
    Interaction(Vec<String>),      // var1:var2
    Nested(String, String),        // var1/var2
}
```

Parse from string: `"~ group + age + (1|subject)"`

### TransformedMatrix

Result of normalization/transformation, preserving provenance.

```rust
pub struct TransformedMatrix {
    data: DenseMatrix<f64>,        // transformations typically densify
    taxa_ids: Vec<String>,
    sample_ids: Vec<String>,
    transformation: TransformationType,
    params: TransformParams,
}
```

### DAResult

Output from a differential abundance test.

```rust
pub struct DAResult {
    taxa_id: String,
    estimate: f64,                 // log fold change or coefficient
    std_error: f64,
    statistic: f64,                // t, z, or chi-squared
    p_value: f64,
    adjusted_p: Option<f64>,       // after multiple testing correction
    
    // Prevalence metadata
    prevalence_group1: f64,
    prevalence_group2: f64,
    prevalence_tier: PrevalenceTier,
    
    // Confidence annotation
    confidence: Confidence,
}

pub enum PrevalenceTier {
    Ubiquitous,      // >50% in all groups
    Moderate,        // 10-50% 
    GroupSpecific,   // >20% difference between groups
    Rare,            // <10% in all groups
}

pub enum Confidence {
    High,            // sufficient data, standard testing
    Medium,          // sparse but testable
    Low,             // very sparse, interpret cautiously
    Untestable,      // insufficient data
}
```

---

## Primitive Categories

### 1. Profiling Primitives

Read-only analysis of data characteristics to inform pipeline selection.

#### `profile_sparsity`

```rust
pub struct SparsityProfile {
    overall_sparsity: f64,                    // proportion of zeros
    per_sample_sparsity: Vec<f64>,
    per_taxon_sparsity: Vec<f64>,
    zero_pattern: ZeroPattern,
}

pub enum ZeroPattern {
    Random,          // zeros distributed randomly
    Structured,      // zeros cluster by sample or taxon
    BlockDiagonal,   // distinct communities with little overlap
}

fn profile_sparsity(counts: &CountMatrix) -> SparsityProfile;
```

#### `profile_prevalence`

```rust
pub struct PrevalenceProfile {
    overall: Vec<f64>,                        // per-taxon overall prevalence
    by_group: HashMap<String, Vec<f64>>,      // per-taxon, per-group prevalence
    differential: Vec<f64>,                   // max(group) - min(group)
    tiers: Vec<PrevalenceTier>,               // classification per taxon
    tier_counts: HashMap<PrevalenceTier, usize>,
}

fn profile_prevalence(
    counts: &CountMatrix, 
    groups: &[String]
) -> PrevalenceProfile;
```

#### `profile_library_size`

```rust
pub struct LibrarySizeProfile {
    sizes: Vec<u64>,
    mean: f64,
    median: f64,
    cv: f64,                                  // coefficient of variation
    size_by_group: HashMap<String, SummaryStats>,
    size_imbalance: f64,                      // between-group imbalance
}

fn profile_library_size(counts: &CountMatrix) -> LibrarySizeProfile;
```

#### `profile_compositionality`

```rust
pub struct CompositionalityProfile {
    dominance: f64,                           // proportion from top N taxa
    evenness: f64,                            // Shannon evenness
    dominant_taxa: Vec<(String, f64)>,        // taxa contributing >X%
    compositional_severity: CompositionalSeverity,
}

pub enum CompositionalSeverity {
    Low,             // relatively even, CLR stable
    Moderate,        // some dominance, CLR usable with caution
    High,            // few taxa dominate, CLR unstable
}

fn profile_compositionality(counts: &CountMatrix) -> CompositionalityProfile;
```

#### `profile_data` (aggregate)

```rust
pub struct DataProfile {
    sparsity: SparsityProfile,
    prevalence: PrevalenceProfile,
    library_size: LibrarySizeProfile,
    compositionality: CompositionalityProfile,
    n_samples: usize,
    n_taxa: usize,
    n_groups: usize,
    design_type: DesignType,
}

pub enum DesignType {
    Independent,
    Paired,
    Longitudinal,
    Nested,
}

fn profile_data(
    counts: &CountMatrix, 
    metadata: &Metadata, 
    formula: &Formula
) -> DataProfile;
```

---

### 2. Filtering Primitives

Taxa and sample filtering operations.

#### `filter_prevalence_overall`

```rust
pub struct FilterResult {
    kept: CountMatrix,
    removed: CountMatrix,
    kept_ids: Vec<String>,
    removed_ids: Vec<String>,
}

fn filter_prevalence_overall(
    counts: &CountMatrix,
    min_prevalence: f64,          // 0.0 to 1.0
) -> FilterResult;
```

#### `filter_prevalence_groupwise`

```rust
pub enum GroupwiseLogic {
    Any,             // keep if prevalent in ANY group
    All,             // keep if prevalent in ALL groups  
    Differential,    // keep if prevalence DIFFERS between groups
}

fn filter_prevalence_groupwise(
    counts: &CountMatrix,
    groups: &[String],
    min_prevalence: f64,
    logic: GroupwiseLogic,
    min_differential: Option<f64>,  // for Differential logic
) -> FilterResult;
```

#### `filter_abundance`

```rust
fn filter_abundance(
    counts: &CountMatrix,
    min_count: u64,               // minimum total count per taxon
    min_proportion: Option<f64>,  // minimum proportion of total reads
) -> FilterResult;
```

#### `filter_library_size`

```rust
fn filter_library_size(
    counts: &CountMatrix,
    min_size: u64,
) -> FilterResult;  // filters SAMPLES, not taxa
```

---

### 3. Stratification Primitives

Partition taxa for tiered analysis rather than filtering.

#### `stratify_by_prevalence`

```rust
pub struct StratifiedCounts {
    ubiquitous: CountMatrix,      // >50% all groups
    moderate: CountMatrix,        // 10-50%
    group_specific: CountMatrix,  // >20% differential
    rare: CountMatrix,            // <10% all groups
    assignments: HashMap<String, PrevalenceTier>,
}

fn stratify_by_prevalence(
    counts: &CountMatrix,
    groups: &[String],
    thresholds: PrevalenceThresholds,
) -> StratifiedCounts;

pub struct PrevalenceThresholds {
    ubiquitous_min: f64,          // default 0.5
    moderate_min: f64,            // default 0.1
    differential_min: f64,        // default 0.2
}
```

---

### 4. Zero Handling Primitives

Strategies for dealing with zeros before transformation/modeling.

#### `zeros_pseudocount`

```rust
fn zeros_pseudocount(
    counts: &CountMatrix,
    value: f64,                   // typically 0.5 or 1.0
) -> DenseMatrix<f64>;           // returns float matrix
```

#### `zeros_multiplicative_replacement`

```rust
// Replaces zeros with small value proportional to minimum non-zero
fn zeros_multiplicative_replacement(
    counts: &CountMatrix,
    delta: Option<f64>,           // replacement factor, default 0.65
) -> DenseMatrix<f64>;
```

#### `zeros_bayesian_impute`

```rust
// BMDD-style Bayesian imputation
fn zeros_bayesian_impute(
    counts: &CountMatrix,
    n_samples: usize,             // posterior samples
    seed: u64,
) -> Vec<DenseMatrix<f64>>;      // returns multiple imputed matrices
```

#### `zeros_detect_structural`

```rust
// ANCOM-BC style structural zero detection
pub struct StructuralZeroResult {
    is_structural: Vec<Vec<bool>>,  // taxa × samples
    structural_by_group: HashMap<String, Vec<bool>>,  // taxa that are structural zeros in entire groups
}

fn zeros_detect_structural(
    counts: &CountMatrix,
    groups: &[String],
    threshold: f64,               // prevalence threshold for structural
) -> StructuralZeroResult;
```

---

### 5. Normalization Primitives

Transform counts to account for compositionality and library size.

#### `norm_clr`

```rust
// Centered log-ratio transformation
fn norm_clr(
    counts: &DenseMatrix<f64>,    // after zero handling
) -> TransformedMatrix;
```

#### `norm_alr`

```rust
// Additive log-ratio (relative to reference taxon)
fn norm_alr(
    counts: &DenseMatrix<f64>,
    reference: ReferenceSelection,
) -> TransformedMatrix;

pub enum ReferenceSelection {
    Index(usize),
    TaxonId(String),
    LeastVariable,                // auto-select most stable taxon
}
```

#### `norm_tmm`

```rust
// Trimmed mean of M-values (edgeR style)
fn norm_tmm(
    counts: &CountMatrix,
    trim_m: f64,                  // default 0.3
    trim_a: f64,                  // default 0.05
) -> NormFactors;

pub struct NormFactors {
    factors: Vec<f64>,            // per-sample normalization factors
}
```

#### `norm_css`

```rust
// Cumulative sum scaling (metagenomeSeq style)
fn norm_css(
    counts: &CountMatrix,
    quantile: f64,                // default 0.5
) -> NormFactors;
```

#### `norm_median_ratios`

```rust
// DESeq2-style median of ratios
fn norm_median_ratios(
    counts: &CountMatrix,
) -> NormFactors;
```

#### `norm_rarefaction`

```rust
// Rarefaction (subsampling to even depth)
fn norm_rarefaction(
    counts: &CountMatrix,
    depth: RarefactionDepth,
    seed: u64,
) -> CountMatrix;                 // returns integer counts

pub enum RarefactionDepth {
    Fixed(u64),
    MinSample,
    Quantile(f64),                // e.g., 0.1 for 10th percentile
}
```

#### `norm_aldex_montecarlo`

```rust
// ALDEx2-style Monte Carlo sampling from Dirichlet
fn norm_aldex_montecarlo(
    counts: &CountMatrix,
    n_instances: usize,           // typically 128
    seed: u64,
) -> Vec<TransformedMatrix>;     // returns multiple CLR-transformed instances
```

---

### 6. Model Fitting Primitives

Statistical model fitting on transformed data.

#### `model_lm`

```rust
// Simple linear model
pub struct LmFit {
    coefficients: Vec<f64>,
    std_errors: Vec<f64>,
    residuals: Vec<f64>,
    r_squared: f64,
    df_residual: usize,
}

fn model_lm(
    y: &[f64],                    // response (one taxon's transformed abundances)
    x: &DesignMatrix,             // design matrix from formula + metadata
) -> LmFit;
```

#### `model_glm`

```rust
pub enum GlmFamily {
    Gaussian,
    Poisson,
    NegativeBinomial { theta: Option<f64> },
    Binomial,
    Gamma,
}

pub struct GlmFit {
    coefficients: Vec<f64>,
    std_errors: Vec<f64>,
    deviance: f64,
    null_deviance: f64,
    df_residual: usize,
    dispersion: f64,
    converged: bool,
}

fn model_glm(
    y: &[f64],
    x: &DesignMatrix,
    family: GlmFamily,
    max_iter: usize,
) -> Result<GlmFit, FitError>;
```

#### `model_lmm`

```rust
// Linear mixed model
pub struct LmmFit {
    fixed_effects: Vec<f64>,
    fixed_std_errors: Vec<f64>,
    random_effects: HashMap<String, Vec<f64>>,
    variance_components: Vec<f64>,
    log_likelihood: f64,
}

fn model_lmm(
    y: &[f64],
    fixed: &DesignMatrix,
    random: &RandomEffectSpec,
) -> Result<LmmFit, FitError>;
```

#### `model_betabinom`

```rust
// Beta-binomial regression (corncob-style)
pub struct BetaBinomFit {
    abundance_coef: Vec<f64>,     // coefficients for abundance
    dispersion_coef: Vec<f64>,    // coefficients for dispersion
    std_errors: Vec<f64>,
    log_likelihood: f64,
}

fn model_betabinom(
    counts: &[u64],               // raw counts for one taxon
    totals: &[u64],               // library sizes
    x_abundance: &DesignMatrix,
    x_dispersion: &DesignMatrix,  // can model heterogeneous dispersion
) -> Result<BetaBinomFit, FitError>;
```

#### `model_hurdle`

```rust
// Hurdle model: separate presence/absence and abundance
pub struct HurdleFit {
    presence_fit: GlmFit,         // binomial for zero/non-zero
    abundance_fit: GlmFit,        // for positive counts only
}

fn model_hurdle(
    counts: &[u64],
    x: &DesignMatrix,
    abundance_family: GlmFamily,  // typically NegativeBinomial or Gamma
) -> Result<HurdleFit, FitError>;
```

#### `model_zinb`

```rust
// Zero-inflated negative binomial
pub struct ZinfFit {
    count_coef: Vec<f64>,
    zero_coef: Vec<f64>,          // inflation model coefficients
    theta: f64,                   // dispersion
    log_likelihood: f64,
}

fn model_zinb(
    counts: &[u64],
    x_count: &DesignMatrix,
    x_zero: &DesignMatrix,
) -> Result<ZinfFit, FitError>;
```

---

### 7. Testing Primitives

Statistical tests on fitted models.

#### `test_wald`

```rust
pub struct TestResult {
    statistic: f64,
    df: Option<f64>,              // degrees of freedom if applicable
    p_value: f64,
}

fn test_wald(
    estimate: f64,
    std_error: f64,
    null_value: f64,              // typically 0.0
) -> TestResult;
```

#### `test_lrt`

```rust
// Likelihood ratio test
fn test_lrt(
    log_lik_full: f64,
    log_lik_reduced: f64,
    df_diff: usize,
) -> TestResult;
```

#### `test_permutation`

```rust
fn test_permutation(
    observed_statistic: f64,
    y: &[f64],
    x: &DesignMatrix,
    test_coef: usize,             // which coefficient to test
    n_perms: usize,
    seed: u64,
) -> TestResult;
```

#### `test_fisher`

```rust
// Fisher's exact test for presence/absence
fn test_fisher(
    present_group1: usize,
    absent_group1: usize,
    present_group2: usize,
    absent_group2: usize,
) -> TestResult;
```

#### `test_wilcoxon`

```rust
fn test_wilcoxon(
    group1: &[f64],
    group2: &[f64],
) -> TestResult;
```

---

### 8. Multiple Testing Correction Primitives

#### `correct_bh`

```rust
// Benjamini-Hochberg
fn correct_bh(p_values: &[f64]) -> Vec<f64>;  // returns adjusted p-values
```

#### `correct_qvalue`

```rust
// Storey's q-value
pub struct QvalueResult {
    q_values: Vec<f64>,
    pi0: f64,                     // estimated proportion of true nulls
}

fn correct_qvalue(p_values: &[f64]) -> QvalueResult;
```

#### `correct_locfdr`

```rust
// Local false discovery rate
pub struct LocfdrResult {
    lfdr: Vec<f64>,               // local fdr per test
    fdr: Vec<f64>,                // cumulative fdr
}

fn correct_locfdr(
    z_scores: &[f64],
    null_proportion: Option<f64>,
) -> LocfdrResult;
```

---

### 9. Spike-in Validation Primitives

Empirical validation by injecting known signal.

#### `spike_abundance`

```rust
pub struct SpikeSpec {
    spiked_taxa: Vec<String>,
    fold_changes: Vec<f64>,       // true effect sizes
    affected_group: String,
}

pub struct SpikedData {
    counts: CountMatrix,
    spec: SpikeSpec,
    original_counts: CountMatrix,
}

fn spike_abundance(
    counts: &CountMatrix,
    groups: &[String],
    n_spike: usize,
    fold_change: f64,
    target_group: &str,
    selection: SpikeSelection,
    seed: u64,
) -> SpikedData;

pub enum SpikeSelection {
    Random,
    ByPrevalenceTier(PrevalenceTier),
    ByAbundance { min: f64, max: f64 },
    Specific(Vec<String>),
}
```

#### `spike_presence`

```rust
// Spike presence/absence: convert zeros to non-zeros in one group
fn spike_presence(
    counts: &CountMatrix,
    groups: &[String],
    n_spike: usize,
    target_prevalence: f64,       // desired prevalence in target group
    target_group: &str,
    abundance_level: AbundanceLevel,
    seed: u64,
) -> SpikedData;

pub enum AbundanceLevel {
    Low,                          // 10th percentile of non-zeros
    Median,
    High,                         // 90th percentile
    Fixed(u64),
}
```

#### `spike_hurdle`

```rust
// Spike both presence AND abundance (for hurdle model validation)
fn spike_hurdle(
    counts: &CountMatrix,
    groups: &[String],
    n_spike: usize,
    presence_effect: f64,         // increase in prevalence
    abundance_effect: f64,        // fold change in non-zeros
    target_group: &str,
    seed: u64,
) -> SpikedData;
```

#### `shuffle_labels`

```rust
// Negative control: permute group labels
fn shuffle_labels(
    groups: &[String],
    seed: u64,
) -> Vec<String>;
```

#### `evaluate_spikes`

```rust
pub struct SpikeEvaluation {
    // Detection metrics
    true_positives: usize,
    false_positives: usize,
    false_negatives: usize,
    true_negatives: usize,
    
    // Rates
    sensitivity: f64,             // TP / (TP + FN)
    specificity: f64,             // TN / (TN + FP)
    precision: f64,               // TP / (TP + FP)
    fdr: f64,                     // FP / (TP + FP)
    
    // Effect size recovery
    effect_correlation: f64,      // correlation of estimated vs true effects
    effect_bias: f64,             // mean(estimated - true)
    
    // By prevalence tier
    by_tier: HashMap<PrevalenceTier, TierMetrics>,
}

fn evaluate_spikes(
    results: &[DAResult],
    spike_spec: &SpikeSpec,
    fdr_threshold: f64,
) -> SpikeEvaluation;
```

#### `run_spike_validation`

```rust
pub struct ValidationConfig {
    n_iterations: usize,          // number of spike-in rounds
    fold_changes: Vec<f64>,       // effect sizes to test
    n_spike_per_tier: usize,
    include_presence_spikes: bool,
    include_hurdle_spikes: bool,
    fdr_thresholds: Vec<f64>,     // e.g., [0.05, 0.1, 0.2]
    seed: u64,
}

pub struct ValidationResult {
    evaluations: Vec<SpikeEvaluation>,
    summary: ValidationSummary,
}

pub struct ValidationSummary {
    mean_sensitivity: f64,
    mean_fdr: f64,
    sensitivity_by_effect_size: HashMap<OrderedFloat<f64>, f64>,
    sensitivity_by_tier: HashMap<PrevalenceTier, f64>,
    recommended_fdr_threshold: f64,
}

fn run_spike_validation(
    counts: &CountMatrix,
    metadata: &Metadata,
    pipeline: &Pipeline,
    config: &ValidationConfig,
) -> ValidationResult;
```

---

## Pipeline Composition

### Pipeline Specification

```rust
pub struct Pipeline {
    name: String,
    steps: Vec<PipelineStep>,
}

pub enum PipelineStep {
    Filter(FilterSpec),
    Stratify(StratifySpec),
    ZeroHandle(ZeroHandleSpec),
    Normalize(NormalizeSpec),
    Model(ModelSpec),
    Test(TestSpec),
    Correct(CorrectionSpec),
}

// Example specs
pub struct FilterSpec {
    method: FilterMethod,
    params: HashMap<String, Value>,
}

pub enum FilterMethod {
    PrevalenceOverall,
    PrevalenceGroupwise,
    Abundance,
}
```

### Pipeline Execution

```rust
pub struct PipelineResult {
    results: Vec<DAResult>,
    intermediate: HashMap<String, IntermediateResult>,  // optional caching
    timing: HashMap<String, Duration>,
    warnings: Vec<String>,
}

fn run_pipeline(
    counts: &CountMatrix,
    metadata: &Metadata,
    formula: &Formula,
    pipeline: &Pipeline,
) -> Result<PipelineResult, PipelineError>;
```

### Tiered Pipeline

```rust
// Special pipeline type that runs different sub-pipelines per prevalence tier
pub struct TieredPipeline {
    ubiquitous: Pipeline,
    moderate: Pipeline,
    group_specific: Pipeline,
    rare: Option<Pipeline>,       // None to skip rare taxa
    stratification_thresholds: PrevalenceThresholds,
}

fn run_tiered_pipeline(
    counts: &CountMatrix,
    metadata: &Metadata,
    formula: &Formula,
    pipeline: &TieredPipeline,
) -> Result<PipelineResult, PipelineError>;
```

---

## Serialization

All pipeline specifications should be serializable for reproducibility.

```rust
// Support YAML, JSON, and RON formats
impl Pipeline {
    fn to_yaml(&self) -> String;
    fn from_yaml(s: &str) -> Result<Self, ParseError>;
    
    fn to_json(&self) -> String;
    fn from_json(s: &str) -> Result<Self, ParseError>;
}
```

Example YAML:

```yaml
name: "virome_sparse_pipeline"
steps:
  - filter:
      method: prevalence_groupwise
      params:
        min_prevalence: 0.05
        logic: any
  - zero_handle:
      method: pseudocount
      params:
        value: 0.5
  - normalize:
      method: clr
  - model:
      method: lm
  - test:
      method: wald
      params:
        coefficient: 1
  - correct:
      method: bh
```

---

## CLI Interface

```bash
# Profile data
daa profile --counts counts.tsv --metadata meta.tsv --formula "~ group"

# Run pipeline
daa run --counts counts.tsv --metadata meta.tsv --formula "~ group" \
        --pipeline pipeline.yaml --output results.tsv

# Run spike-in validation
daa validate --counts counts.tsv --metadata meta.tsv --formula "~ group" \
             --pipeline pipeline.yaml --iterations 100 --output validation.json

# Compare pipelines
daa compare --counts counts.tsv --metadata meta.tsv --formula "~ group" \
            --pipelines pipeline1.yaml pipeline2.yaml pipeline3.yaml \
            --output comparison.json
```

---

## Implementation Notes

### Dependencies (suggested)

- `nalgebra` or `ndarray` for dense matrices
- `sprs` for sparse matrices  
- `statrs` for statistical distributions
- `rayon` for parallelization
- `serde` for serialization
- `clap` for CLI

### Performance Considerations

1. **Sparse operations**: Keep data sparse as long as possible. Only densify after transformation.

2. **Parallelization**: 
   - Per-taxon model fitting is embarrassingly parallel
   - Spike-in iterations are independent
   - Use rayon for both

3. **Memory**: For large studies (>1000 samples), consider:
   - Memory-mapped files for count matrix
   - Streaming per-taxon rather than loading all

4. **Caching**: Model fitting is expensive. Cache intermediate results (design matrices, normalization factors) when running multiple pipelines on same data.

### Testing Strategy

1. **Unit tests**: Each primitive in isolation with known inputs/outputs

2. **Integration tests**: Standard pipelines reproduce published method results
   - LinDA pipeline matches LinDA R package
   - DESeq2 pipeline matches DESeq2 results (within tolerance)

3. **Spike-in tests**: Validate that spike-in framework correctly identifies known effects

4. **Property tests**: Use proptest for fuzzing with random count matrices

---

## Example Usage

```rust
use daa_primitives::*;

fn main() -> Result<(), Error> {
    // Load data
    let counts = CountMatrix::from_tsv("counts.tsv")?;
    let metadata = Metadata::from_tsv("meta.tsv")?;
    let formula = Formula::parse("~ disease + age")?;
    
    // Profile
    let profile = profile_data(&counts, &metadata, &formula);
    println!("Sparsity: {:.1}%", profile.sparsity.overall_sparsity * 100.0);
    println!("Prevalence tiers: {:?}", profile.prevalence.tier_counts);
    
    // Build tiered pipeline based on profile
    let pipeline = TieredPipeline {
        ubiquitous: Pipeline::from_yaml(include_str!("pipelines/abundant.yaml"))?,
        moderate: Pipeline::from_yaml(include_str!("pipelines/moderate.yaml"))?,
        group_specific: Pipeline::from_yaml(include_str!("pipelines/presence.yaml"))?,
        rare: None,
        stratification_thresholds: PrevalenceThresholds::default(),
    };
    
    // Validate with spike-ins
    let validation = run_spike_validation(
        &counts,
        &metadata,
        &pipeline,
        &ValidationConfig {
            n_iterations: 100,
            fold_changes: vec![2.0, 4.0, 8.0],
            n_spike_per_tier: 10,
            include_presence_spikes: true,
            include_hurdle_spikes: true,
            fdr_thresholds: vec![0.05, 0.1],
            seed: 42,
        },
    )?;
    
    println!("Validation: {:.1}% sensitivity, {:.1}% FDR", 
             validation.summary.mean_sensitivity * 100.0,
             validation.summary.mean_fdr * 100.0);
    
    // Run on real data
    let results = run_tiered_pipeline(&counts, &metadata, &formula, &pipeline)?;
    
    // Output
    results.to_tsv("results.tsv")?;
    
    Ok(())
}
```

---

## Future Extensions

### Near-term (Completed)

1. **Prevalence threshold optimization**: DONE - Implemented via spike-in validation with MaxF1, MaxSensitivity, MaxEfficiency, and MinFDR criteria.

2. **Group-specific prevalence optimization**: DONE - Supports overall, any-group, and all-groups filter logic.

3. **LLM-friendly data profiling**: DONE - Generates structured YAML/markdown summaries for AI-assisted pipeline recommendation.

4. **TMM normalization**: DONE - edgeR-style trimmed mean of M-values with configurable parameters.

5. **Permutation tests**: DONE - Distribution-free hypothesis testing with parallel execution.

### Medium-term (Model Extensions)

4. **Linear mixed models (LMM)**: Support longitudinal and repeated-measures designs with random effects for subjects.

5. **Beta-binomial models**: corncob-style modeling of both abundance and dispersion as functions of covariates.

6. **Hurdle models**: Explicit two-part models separating presence/absence from abundance among present features.

### Long-term (Advanced Features)

8. **Phylogenetic-aware methods**: Incorporate tree structure into testing (e.g., PGLMM)

9. **Multi-omics integration**: Joint modeling with metabolomics, host transcriptomics

10. **Network effects**: Test for differential co-occurrence, not just abundance

11. **Batch effects**: ComBat-style batch correction as a primitive

### Completed (Originally Future)

- ~~Effect size shrinkage~~: Implemented with empirical Bayes estimation (Normal/Adaptive methods)
- ~~TMM normalization~~: Implemented with configurable trimming and reference sample selection
- ~~Permutation tests~~: Implemented with parallel execution via rayon
- ~~Prevalence threshold optimization~~: Implemented with spike-in validation and multiple criteria
- ~~LLM-friendly profiling~~: Implemented with YAML and markdown output formats
