//! LLM-friendly data profiling for AI-assisted pipeline selection.
//!
//! This module generates comprehensive YAML profiles optimized for LLM consumption,
//! enabling AI-assisted workflow design based on data characteristics.

use crate::data::{CountMatrix, Metadata};
use crate::error::Result;
use crate::profile::{
    compositionality::{profile_compositionality, CompositionalityProfile, DominanceCategory},
    profile_library_size, profile_prevalence, profile_sparsity,
    PrevalenceProfile,
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Complete LLM-friendly profile for AI-assisted pipeline selection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LlmProfile {
    /// Generation metadata.
    pub meta: ProfileMeta,
    /// Basic data dimensions and grouping.
    pub data_summary: DataSummary,
    /// Sparsity analysis with implications.
    pub sparsity: SparsityAnalysis,
    /// Prevalence distribution analysis.
    pub prevalence: PrevalenceAnalysis,
    /// Library size variation analysis.
    pub library_size: LibrarySizeAnalysis,
    /// Compositionality and dominance analysis.
    pub compositionality: CompositionalityAnalysis,
    /// Pipeline component recommendations.
    pub recommendations: Recommendations,
    /// Validation command suggestions.
    pub validation: ValidationSuggestions,
    /// Warnings requiring attention.
    pub warnings: Vec<Warning>,
    /// Ready-to-use pipeline configurations.
    pub pipelines: SuggestedPipelines,
    /// Condensed summary for LLM context efficiency.
    pub llm_summary: String,
}

impl LlmProfile {
    /// Render as YAML string.
    pub fn to_yaml(&self) -> Result<String> {
        Ok(serde_yaml::to_string(self)?)
    }

    /// Render as markdown (alternative format).
    pub fn to_markdown(&self) -> String {
        let mut md = String::new();

        md.push_str("# Data Profile for LLM-Assisted Pipeline Selection\n\n");

        // Data summary
        md.push_str("## Data Summary\n\n");
        md.push_str(&format!("- **Features**: {}\n", self.data_summary.features));
        md.push_str(&format!("- **Samples**: {}\n", self.data_summary.samples));
        md.push_str(&format!("- **Group column**: {}\n", self.data_summary.groups.column));
        md.push_str(&format!(
            "- **Groups**: {} ({})\n",
            self.data_summary.groups.levels.join(", "),
            if self.data_summary.groups.balanced {
                "balanced"
            } else {
                "imbalanced"
            }
        ));
        md.push('\n');

        // Sparsity
        md.push_str("## Sparsity Analysis\n\n");
        md.push_str(&format!(
            "- **Overall sparsity**: {:.1}% ({})\n",
            self.sparsity.overall * 100.0,
            self.sparsity.category.as_str()
        ));
        md.push_str(&format!("\n{}\n\n", self.sparsity.context));

        // Prevalence
        md.push_str("## Prevalence Analysis\n\n");
        md.push_str(&format!("- **Mean prevalence**: {:.1}%\n", self.prevalence.mean * 100.0));
        md.push_str(&format!(
            "- **Group-specific features**: {} ({:.1}%)\n",
            self.prevalence.group_specific.count,
            self.prevalence.group_specific.percent
        ));
        md.push('\n');

        // Library size
        md.push_str("## Library Size Analysis\n\n");
        md.push_str(&format!(
            "- **CV**: {:.2} ({})\n",
            self.library_size.cv,
            self.library_size.category.as_str()
        ));
        if let Some(ref by_group) = self.library_size.by_group {
            md.push_str(&format!(
                "- **Group ratio**: {:.2}x\n",
                by_group.ratio
            ));
        }
        md.push('\n');

        // Compositionality
        md.push_str("## Compositionality Analysis\n\n");
        md.push_str(&format!(
            "- **Top 3 dominance**: {:.1}% ({})\n",
            self.compositionality.dominance.top_3_percent,
            self.compositionality.category.as_str()
        ));
        md.push_str(&format!("- **Evenness**: {:.2}\n", self.compositionality.evenness));
        md.push('\n');

        // Recommendations
        md.push_str("## Recommendations\n\n");
        md.push_str(&format!(
            "### Filtering\n{}\n\n",
            self.recommendations.filtering.prevalence.recommendation
        ));
        md.push_str(&format!(
            "### Zero Handling\n{}\n\n",
            self.recommendations.zero_handling.recommendation
        ));
        md.push_str(&format!(
            "### Normalization\n{}\n\n",
            self.recommendations.normalization.recommendation
        ));
        md.push_str(&format!(
            "### Model\n{}\n\n",
            self.recommendations.model.recommendation
        ));

        // Warnings
        if !self.warnings.is_empty() {
            md.push_str("## Warnings\n\n");
            for warning in &self.warnings {
                md.push_str(&format!(
                    "- **[{}] {}**: {}\n",
                    warning.severity.as_str(),
                    warning.code,
                    warning.message
                ));
            }
            md.push('\n');
        }

        // LLM summary
        md.push_str("## LLM Summary\n\n");
        md.push_str(&self.llm_summary);

        md
    }
}

/// Profile generation metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfileMeta {
    /// Generation timestamp.
    pub generated: String,
    /// Tool version.
    pub version: String,
    /// Tool name.
    pub tool: String,
    /// Source counts file (if known).
    pub source_counts: Option<String>,
    /// Source metadata file (if known).
    pub source_metadata: Option<String>,
}

/// Basic data dimensions and grouping info.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataSummary {
    /// Description of the data.
    pub description: String,
    /// Number of features.
    pub features: usize,
    /// Number of samples.
    pub samples: usize,
    /// Group information.
    pub groups: GroupInfo,
}

/// Group structure information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupInfo {
    /// Column name for grouping.
    pub column: String,
    /// Group levels.
    pub levels: Vec<String>,
    /// Sample counts per group.
    pub sizes: HashMap<String, usize>,
    /// Whether groups are balanced (within 1.5x).
    pub balanced: bool,
}

/// Sparsity analysis with context and implications.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparsityAnalysis {
    /// Overall sparsity (fraction of zeros).
    pub overall: f64,
    /// Sparsity category.
    pub category: SparsityCategory,
    /// Verbose context explaining what this means.
    pub context: String,
    /// Sparsity by group.
    pub by_group: HashMap<String, f64>,
    /// Interpretation of group differences.
    pub group_interpretation: String,
    /// Implications for analysis.
    pub implications: Vec<Implication>,
}

/// Sparsity category based on thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SparsityCategory {
    /// <50% zeros
    Low,
    /// 50-70% zeros
    Moderate,
    /// 70-90% zeros
    High,
    /// >90% zeros
    Extreme,
}

impl SparsityCategory {
    /// Categorize sparsity value.
    pub fn from_value(sparsity: f64) -> Self {
        if sparsity < 0.50 {
            SparsityCategory::Low
        } else if sparsity < 0.70 {
            SparsityCategory::Moderate
        } else if sparsity < 0.90 {
            SparsityCategory::High
        } else {
            SparsityCategory::Extreme
        }
    }

    /// Get string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            SparsityCategory::Low => "low",
            SparsityCategory::Moderate => "moderate",
            SparsityCategory::High => "high",
            SparsityCategory::Extreme => "extreme",
        }
    }

    /// Get verbose context for this category.
    pub fn context(&self) -> &'static str {
        match self {
            SparsityCategory::Low => {
                "Low sparsity (<50% zeros) indicates a well-sampled community where most \
                 features are detected in most samples. Standard statistical methods work \
                 well, and zero-handling is straightforward. Pseudocount addition is \
                 sufficient for CLR transformation."
            }
            SparsityCategory::Moderate => {
                "Moderate sparsity (50-70% zeros) is typical for many microbiome datasets. \
                 Zero handling with pseudocount is recommended for CLR. Consider filtering \
                 very rare features to reduce noise while preserving signal."
            }
            SparsityCategory::High => {
                "High sparsity (70-90% zeros) suggests either shallow sequencing, a diverse \
                 community, or both. Zero-inflated models (ZINB) may better capture the data \
                 structure. Prevalence filtering is important to remove noise features."
            }
            SparsityCategory::Extreme => {
                "Extreme sparsity (>90% zeros) presents significant analytical challenges. \
                 Most features are absent from most samples, making abundance-based testing \
                 difficult. Consider presence/absence testing, aggressive filtering, or \
                 specialized sparse-data methods."
            }
        }
    }
}

/// Prevalence distribution analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceAnalysis {
    /// Mean prevalence across features.
    pub mean: f64,
    /// Median prevalence.
    pub median: f64,
    /// Prevalence distribution by tier.
    pub distribution: PrevalenceDistribution,
    /// Group-specific feature analysis.
    pub group_specific: GroupSpecificFeatures,
    /// Implications for analysis.
    pub implications: Vec<Implication>,
}

/// Prevalence distribution across tiers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrevalenceDistribution {
    /// Very rare features (<5% prevalence).
    pub very_rare: TierInfo,
    /// Rare features (5-10% prevalence).
    pub rare: TierInfo,
    /// Low prevalence (10-25%).
    pub low: TierInfo,
    /// Moderate prevalence (25-50%).
    pub moderate: TierInfo,
    /// Common features (50-75%).
    pub common: TierInfo,
    /// Ubiquitous features (>75%).
    pub ubiquitous: TierInfo,
}

/// Information about a prevalence tier.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TierInfo {
    /// Upper threshold for this tier.
    pub threshold: f64,
    /// Count of features in this tier.
    pub count: usize,
    /// Percent of total features.
    pub percent: f64,
}

/// Group-specific feature analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupSpecificFeatures {
    /// Count of group-specific features.
    pub count: usize,
    /// Percent of total features.
    pub percent: f64,
    /// Definition used.
    pub definition: String,
    /// Verbose context.
    pub context: String,
}

/// Library size analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibrarySizeAnalysis {
    /// Mean library size.
    pub mean: f64,
    /// Median library size.
    pub median: f64,
    /// Minimum library size.
    pub min: u64,
    /// Maximum library size.
    pub max: u64,
    /// Coefficient of variation.
    pub cv: f64,
    /// Library size variation category.
    pub category: LibrarySizeCategory,
    /// Library size by group.
    pub by_group: Option<LibrarySizeByGroup>,
    /// Verbose context.
    pub context: String,
    /// Implications for analysis.
    pub implications: Vec<Implication>,
}

/// Library size statistics by group.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibrarySizeByGroup {
    /// Stats per group.
    pub stats: HashMap<String, GroupLibraryStats>,
    /// Ratio of largest to smallest group mean.
    pub ratio: f64,
    /// Interpretation of group differences.
    pub interpretation: String,
}

/// Library size stats for a group.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupLibraryStats {
    /// Mean library size.
    pub mean: f64,
    /// Median library size.
    pub median: f64,
}

/// Library size variation category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LibrarySizeCategory {
    /// CV < 0.5
    Low,
    /// CV 0.5-1.0
    Moderate,
    /// CV 1.0-2.0
    High,
    /// CV > 2.0
    Extreme,
}

impl LibrarySizeCategory {
    /// Categorize from CV value.
    pub fn from_cv(cv: f64) -> Self {
        if cv < 0.5 {
            LibrarySizeCategory::Low
        } else if cv < 1.0 {
            LibrarySizeCategory::Moderate
        } else if cv < 2.0 {
            LibrarySizeCategory::High
        } else {
            LibrarySizeCategory::Extreme
        }
    }

    /// Get string representation.
    pub fn as_str(&self) -> &'static str {
        match self {
            LibrarySizeCategory::Low => "low_variation",
            LibrarySizeCategory::Moderate => "moderate_variation",
            LibrarySizeCategory::High => "high_variation",
            LibrarySizeCategory::Extreme => "extreme_variation",
        }
    }

    /// Get verbose context.
    pub fn context(&self) -> &'static str {
        match self {
            LibrarySizeCategory::Low => {
                "Low library size variation (CV < 0.5) indicates consistent sequencing depth. \
                 Normalization is still recommended but has minimal impact on results."
            }
            LibrarySizeCategory::Moderate => {
                "Moderate library size variation (CV 0.5-1.0) is typical. Standard \
                 normalization methods (CLR, TSS) handle this well."
            }
            LibrarySizeCategory::High => {
                "High library size variation (CV 1.0-2.0) requires careful normalization. \
                 Check for outlier samples and consider whether variation is biological \
                 or technical."
            }
            LibrarySizeCategory::Extreme => {
                "Extreme library size variation (CV > 2.0) suggests possible technical \
                 issues or outlier samples. Consider filtering low-depth samples or \
                 investigating the source of variation before analysis."
            }
        }
    }
}

/// Compositionality analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CompositionalityAnalysis {
    /// Dominance metrics.
    pub dominance: DominanceMetrics,
    /// Shannon evenness (0-1).
    pub evenness: f64,
    /// Gini coefficient (0-1).
    pub gini: f64,
    /// Dominance category.
    pub category: String,
    /// Top dominant features.
    pub dominant_features: Vec<DominantFeatureInfo>,
    /// Verbose context.
    pub context: String,
    /// Implications for analysis.
    pub implications: Vec<Implication>,
}

/// Dominance metrics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DominanceMetrics {
    /// Percent from top 1 feature.
    pub top_1_percent: f64,
    /// Percent from top 3 features.
    pub top_3_percent: f64,
    /// Percent from top 10 features.
    pub top_10_percent: f64,
}

/// Information about a dominant feature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DominantFeatureInfo {
    /// Feature identifier.
    pub id: String,
    /// Mean proportion.
    pub mean_proportion: f64,
}

/// An implication with code and message.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Implication {
    /// Short code for the implication.
    pub code: String,
    /// Human-readable message.
    pub message: String,
    /// Severity (optional).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub severity: Option<Severity>,
}

/// Severity level for warnings and implications.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    Low,
    Moderate,
    High,
}

impl Severity {
    pub fn as_str(&self) -> &'static str {
        match self {
            Severity::Low => "low",
            Severity::Moderate => "moderate",
            Severity::High => "high",
        }
    }
}

/// Pipeline component recommendations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Recommendations {
    /// Filtering recommendation.
    pub filtering: FilteringRecommendation,
    /// Zero handling recommendation.
    pub zero_handling: ComponentRecommendation,
    /// Normalization recommendation.
    pub normalization: ComponentRecommendation,
    /// Model recommendation.
    pub model: ComponentRecommendation,
}

/// Filtering recommendation details.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FilteringRecommendation {
    /// Prevalence filtering recommendation.
    pub prevalence: ComponentRecommendation,
    /// Abundance filtering recommendation (optional).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub abundance: Option<ComponentRecommendation>,
}

/// A component recommendation with alternatives.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComponentRecommendation {
    /// Primary recommendation.
    pub recommendation: String,
    /// Confidence level.
    pub confidence: Confidence,
    /// Rationale for the recommendation.
    pub rationale: String,
    /// Alternative options.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub alternatives: Vec<Alternative>,
}

/// Confidence level for recommendations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Confidence {
    Low,
    Medium,
    High,
}

/// An alternative option.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alternative {
    /// The alternative option.
    pub option: String,
    /// When to use this alternative.
    pub when: String,
}

/// Validation suggestions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationSuggestions {
    /// Description of why validation matters.
    pub description: String,
    /// Essential validations.
    pub essential: Vec<ValidationStep>,
    /// Recommended additional validations.
    pub recommended: Vec<ValidationStep>,
}

/// A validation step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationStep {
    /// Name of the validation.
    pub name: String,
    /// Purpose of this validation.
    pub purpose: String,
    /// Command to run.
    pub command: String,
    /// Expected outcome.
    pub expected: String,
    /// How to interpret results.
    pub interpretation: String,
}

/// A warning about the data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Warning {
    /// Severity level.
    pub severity: Severity,
    /// Category of warning.
    pub category: String,
    /// Warning code.
    pub code: String,
    /// Short message.
    pub message: String,
    /// Detailed explanation.
    pub details: String,
    /// Suggested action.
    pub action: String,
}

/// Suggested pipeline configurations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuggestedPipelines {
    /// Recommended pipeline.
    pub recommended: PipelineConfig,
    /// Alternative pipelines.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub alternatives: Vec<PipelineConfig>,
}

/// A pipeline configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    /// Pipeline name.
    pub name: String,
    /// Description of when to use.
    pub description: String,
    /// YAML configuration.
    pub config: String,
}

/// Generate an LLM-friendly profile from count data.
pub fn profile_for_llm(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
) -> Result<LlmProfile> {
    // Get base profiles
    let sparsity_profile = profile_sparsity(counts);
    let prevalence_profile = profile_prevalence(counts);
    let library_profile = profile_library_size(counts);
    let comp_profile = profile_compositionality(counts);

    // Get group information
    let groups = extract_group_info(metadata, group_column)?;

    // Calculate group-specific metrics
    let sparsity_by_group = calculate_sparsity_by_group(counts, metadata, group_column)?;
    let library_by_group = calculate_library_size_by_group(counts, metadata, group_column)?;
    let group_specific_features =
        calculate_group_specific_features(counts, metadata, group_column)?;

    // Build categorized analyses
    let sparsity_cat = SparsityCategory::from_value(sparsity_profile.sparsity);
    let library_cat = LibrarySizeCategory::from_cv(library_profile.cv);

    // Build sparsity analysis
    let sparsity_analysis = SparsityAnalysis {
        overall: sparsity_profile.sparsity,
        category: sparsity_cat,
        context: sparsity_cat.context().to_string(),
        by_group: sparsity_by_group.clone(),
        group_interpretation: interpret_group_sparsity(&sparsity_by_group),
        implications: sparsity_implications(sparsity_cat),
    };

    // Build prevalence analysis
    let prevalence_analysis = build_prevalence_analysis(&prevalence_profile, &group_specific_features);

    // Build library size analysis
    let library_analysis = LibrarySizeAnalysis {
        mean: library_profile.mean,
        median: library_profile.median,
        min: library_profile.min,
        max: library_profile.max,
        cv: library_profile.cv,
        category: library_cat,
        by_group: Some(library_by_group),
        context: library_cat.context().to_string(),
        implications: library_size_implications(library_cat),
    };

    // Build compositionality analysis
    let comp_analysis = build_compositionality_analysis(&comp_profile);

    // Generate recommendations based on all profiles
    let recommendations = generate_recommendations(
        &sparsity_analysis,
        &prevalence_analysis,
        &library_analysis,
        &comp_analysis,
        &group_specific_features,
    );

    // Generate warnings
    let warnings = generate_warnings(
        &sparsity_analysis,
        &library_analysis,
        &comp_analysis,
        &groups,
    );

    // Generate validation suggestions
    let validation = generate_validation_suggestions(group_column);

    // Generate pipeline suggestions
    let pipelines = generate_pipeline_suggestions(
        &recommendations,
        group_column,
    );

    // Generate LLM summary
    let llm_summary = generate_llm_summary(
        &sparsity_analysis,
        &prevalence_analysis,
        &library_analysis,
        &comp_analysis,
        &recommendations,
        &warnings,
    );

    // Build meta
    let meta = ProfileMeta {
        generated: chrono::Utc::now().to_rfc3339(),
        version: env!("CARGO_PKG_VERSION").to_string(),
        tool: "composable-daa".to_string(),
        source_counts: None,
        source_metadata: None,
    };

    // Build data summary
    let data_summary = DataSummary {
        description: "Microbiome count matrix for differential abundance analysis".to_string(),
        features: counts.n_features(),
        samples: counts.n_samples(),
        groups,
    };

    Ok(LlmProfile {
        meta,
        data_summary,
        sparsity: sparsity_analysis,
        prevalence: prevalence_analysis,
        library_size: library_analysis,
        compositionality: comp_analysis,
        recommendations,
        validation,
        warnings,
        pipelines,
        llm_summary,
    })
}

// ============================================================================
// Helper functions
// ============================================================================

fn extract_group_info(metadata: &Metadata, group_column: &str) -> Result<GroupInfo> {
    let mut sizes: HashMap<String, usize> = HashMap::new();

    for sample_id in metadata.sample_ids() {
        if let Some(var) = metadata.get(sample_id, group_column) {
            if let Some(value) = var.as_categorical() {
                *sizes.entry(value.to_string()).or_insert(0) += 1;
            }
        }
    }

    let levels: Vec<String> = sizes.keys().cloned().collect();
    let max_size = sizes.values().max().copied().unwrap_or(1) as f64;
    let min_size = sizes.values().min().copied().unwrap_or(1) as f64;
    let balanced = max_size / min_size <= 1.5;

    Ok(GroupInfo {
        column: group_column.to_string(),
        levels,
        sizes,
        balanced,
    })
}

fn calculate_sparsity_by_group(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
) -> Result<HashMap<String, f64>> {
    let mut group_sparsity: HashMap<String, (usize, usize)> = HashMap::new();

    for (sample_idx, sample_id) in counts.sample_ids().iter().enumerate() {
        if let Some(var) = metadata.get(sample_id, group_column) {
            if let Some(group) = var.as_categorical() {
                let entry = group_sparsity.entry(group.to_string()).or_insert((0, 0));
                for feat_idx in 0..counts.n_features() {
                    let count = counts.get(feat_idx, sample_idx);
                    if count == 0 {
                        entry.0 += 1;
                    }
                    entry.1 += 1;
                }
            }
        }
    }

    Ok(group_sparsity
        .into_iter()
        .map(|(g, (zeros, total))| (g, zeros as f64 / total as f64))
        .collect())
}

fn calculate_library_size_by_group(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
) -> Result<LibrarySizeByGroup> {
    let mut group_sizes: HashMap<String, Vec<u64>> = HashMap::new();

    for (sample_idx, sample_id) in counts.sample_ids().iter().enumerate() {
        if let Some(var) = metadata.get(sample_id, group_column) {
            if let Some(group) = var.as_categorical() {
                let lib_size: u64 = (0..counts.n_features())
                    .map(|f| counts.get(f, sample_idx))
                    .sum();
                group_sizes
                    .entry(group.to_string())
                    .or_default()
                    .push(lib_size);
            }
        }
    }

    let mut stats: HashMap<String, GroupLibraryStats> = HashMap::new();
    let mut means: Vec<f64> = Vec::new();

    for (group, sizes) in &group_sizes {
        let mean = sizes.iter().sum::<u64>() as f64 / sizes.len() as f64;
        let mut sorted = sizes.clone();
        sorted.sort();
        let median = if sorted.len() % 2 == 0 {
            (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) as f64 / 2.0
        } else {
            sorted[sorted.len() / 2] as f64
        };
        stats.insert(group.clone(), GroupLibraryStats { mean, median });
        means.push(mean);
    }

    let ratio = if means.len() >= 2 {
        let max = means.iter().cloned().fold(f64::MIN, f64::max);
        let min = means.iter().cloned().fold(f64::MAX, f64::min);
        max / min
    } else {
        1.0
    };

    let interpretation = if ratio < 1.5 {
        "Balanced library sizes across groups".to_string()
    } else if ratio < 2.0 {
        "Moderate library size difference between groups".to_string()
    } else {
        format!("Substantial library size imbalance ({:.1}x) - potential confounding", ratio)
    };

    Ok(LibrarySizeByGroup {
        stats,
        ratio,
        interpretation,
    })
}

fn calculate_group_specific_features(
    counts: &CountMatrix,
    metadata: &Metadata,
    group_column: &str,
) -> Result<GroupSpecificFeatures> {
    // Get sample indices per group
    let mut group_samples: HashMap<String, Vec<usize>> = HashMap::new();
    for (sample_idx, sample_id) in counts.sample_ids().iter().enumerate() {
        if let Some(var) = metadata.get(sample_id, group_column) {
            if let Some(group) = var.as_categorical() {
                group_samples
                    .entry(group.to_string())
                    .or_default()
                    .push(sample_idx);
            }
        }
    }

    let groups: Vec<&String> = group_samples.keys().collect();
    if groups.len() < 2 {
        return Ok(GroupSpecificFeatures {
            count: 0,
            percent: 0.0,
            definition: ">20% prevalence difference between groups".to_string(),
            context: "Insufficient groups for comparison".to_string(),
        });
    }

    let mut group_specific_count = 0;

    for feat_idx in 0..counts.n_features() {
        let mut prevalences: Vec<f64> = Vec::new();

        for group in &groups {
            let samples = &group_samples[*group];
            let present = samples
                .iter()
                .filter(|&&s| counts.get(feat_idx, s) > 0)
                .count();
            prevalences.push(present as f64 / samples.len() as f64);
        }

        // Check if prevalence difference > 20%
        let max_prev = prevalences.iter().cloned().fold(f64::MIN, f64::max);
        let min_prev = prevalences.iter().cloned().fold(f64::MAX, f64::min);
        if max_prev - min_prev > 0.20 {
            group_specific_count += 1;
        }
    }

    let n_features = counts.n_features();
    let percent = group_specific_count as f64 / n_features as f64 * 100.0;

    let context = if percent < 5.0 {
        "Few group-specific features - overall prevalence filtering is appropriate".to_string()
    } else if percent < 15.0 {
        "Some group-specific features present - consider groupwise filtering with 'any' logic".to_string()
    } else if percent < 30.0 {
        "Many group-specific features - groupwise filtering with 'any' logic recommended".to_string()
    } else {
        "Dominant group-specific pattern - presence/absence testing may be valuable".to_string()
    };

    Ok(GroupSpecificFeatures {
        count: group_specific_count,
        percent,
        definition: ">20% prevalence difference between groups".to_string(),
        context,
    })
}

fn interpret_group_sparsity(by_group: &HashMap<String, f64>) -> String {
    if by_group.len() < 2 {
        return "Single group - no comparison available".to_string();
    }

    let values: Vec<f64> = by_group.values().cloned().collect();
    let max = values.iter().cloned().fold(f64::MIN, f64::max);
    let min = values.iter().cloned().fold(f64::MAX, f64::min);
    let diff = max - min;

    if diff < 0.05 {
        "Similar sparsity across groups".to_string()
    } else if diff < 0.10 {
        "Slight sparsity difference between groups".to_string()
    } else {
        format!(
            "Notable sparsity difference ({:.0}%) - one group may be more deeply sequenced",
            diff * 100.0
        )
    }
}

fn sparsity_implications(category: SparsityCategory) -> Vec<Implication> {
    match category {
        SparsityCategory::Low => vec![Implication {
            code: "STANDARD_METHODS_OK".to_string(),
            message: "Standard statistical methods work well with this sparsity level".to_string(),
            severity: None,
        }],
        SparsityCategory::Moderate => vec![
            Implication {
                code: "ZERO_HANDLING_REQUIRED".to_string(),
                message: "CLR requires pseudocount to avoid log(0)".to_string(),
                severity: None,
            },
            Implication {
                code: "FILTER_RARE".to_string(),
                message: "Consider filtering very rare features (<5% prevalence)".to_string(),
                severity: None,
            },
        ],
        SparsityCategory::High => vec![
            Implication {
                code: "ZERO_HANDLING_REQUIRED".to_string(),
                message: "CLR requires pseudocount to avoid log(0)".to_string(),
                severity: None,
            },
            Implication {
                code: "CONSIDER_ZINB".to_string(),
                message: "ZINB may better model structural zeros".to_string(),
                severity: None,
            },
            Implication {
                code: "FILTER_REQUIRED".to_string(),
                message: "Prevalence filtering important to reduce noise".to_string(),
                severity: Some(Severity::Moderate),
            },
        ],
        SparsityCategory::Extreme => vec![
            Implication {
                code: "SPARSE_DATA_CHALLENGES".to_string(),
                message: "Extreme sparsity presents significant analytical challenges".to_string(),
                severity: Some(Severity::High),
            },
            Implication {
                code: "CONSIDER_PRESENCE_ABSENCE".to_string(),
                message: "Presence/absence testing may be more appropriate".to_string(),
                severity: Some(Severity::Moderate),
            },
            Implication {
                code: "AGGRESSIVE_FILTERING".to_string(),
                message: "Aggressive prevalence filtering recommended".to_string(),
                severity: Some(Severity::Moderate),
            },
        ],
    }
}

fn build_prevalence_analysis(
    profile: &PrevalenceProfile,
    group_specific: &GroupSpecificFeatures,
) -> PrevalenceAnalysis {
    let n_features = profile.n_features;

    // Calculate distribution tiers (get counts from Vec lengths)
    let very_rare_count = profile.features_below(0.05).len();
    let rare_count = profile.features_in_range(0.05, 0.10).len();
    let low_count = profile.features_in_range(0.10, 0.25).len();
    let moderate_count = profile.features_in_range(0.25, 0.50).len();
    let common_count = profile.features_in_range(0.50, 0.75).len();
    let ubiquitous_count = profile.features_above(0.75).len();

    let distribution = PrevalenceDistribution {
        very_rare: TierInfo {
            threshold: 0.05,
            count: very_rare_count,
            percent: very_rare_count as f64 / n_features as f64 * 100.0,
        },
        rare: TierInfo {
            threshold: 0.10,
            count: rare_count,
            percent: rare_count as f64 / n_features as f64 * 100.0,
        },
        low: TierInfo {
            threshold: 0.25,
            count: low_count,
            percent: low_count as f64 / n_features as f64 * 100.0,
        },
        moderate: TierInfo {
            threshold: 0.50,
            count: moderate_count,
            percent: moderate_count as f64 / n_features as f64 * 100.0,
        },
        common: TierInfo {
            threshold: 0.75,
            count: common_count,
            percent: common_count as f64 / n_features as f64 * 100.0,
        },
        ubiquitous: TierInfo {
            threshold: 1.00,
            count: ubiquitous_count,
            percent: ubiquitous_count as f64 / n_features as f64 * 100.0,
        },
    };

    let mut implications = Vec::new();

    // Add implication based on group-specific percentage
    if group_specific.percent >= 15.0 {
        implications.push(Implication {
            code: "GROUPWISE_FILTERING".to_string(),
            message: "Use groupwise filtering to preserve group-specific signals".to_string(),
            severity: Some(Severity::Moderate),
        });
    }

    if very_rare_count as f64 / n_features as f64 > 0.30 {
        implications.push(Implication {
            code: "MANY_RARE_FEATURES".to_string(),
            message: "Many very rare features - consider filtering to reduce multiple testing burden"
                .to_string(),
            severity: None,
        });
    }

    PrevalenceAnalysis {
        mean: profile.mean_prevalence,
        median: profile.median_prevalence,
        distribution,
        group_specific: group_specific.clone(),
        implications,
    }
}

fn library_size_implications(category: LibrarySizeCategory) -> Vec<Implication> {
    match category {
        LibrarySizeCategory::Low => vec![Implication {
            code: "NORMALIZATION_MINIMAL".to_string(),
            message: "Low variation - normalization has minimal impact".to_string(),
            severity: None,
        }],
        LibrarySizeCategory::Moderate => vec![Implication {
            code: "NORMALIZATION_STANDARD".to_string(),
            message: "Standard normalization (CLR/TSS) is appropriate".to_string(),
            severity: None,
        }],
        LibrarySizeCategory::High => vec![
            Implication {
                code: "NORMALIZATION_CRITICAL".to_string(),
                message: "Normalization is critical for valid comparisons".to_string(),
                severity: Some(Severity::Moderate),
            },
            Implication {
                code: "CHECK_OUTLIERS".to_string(),
                message: "Check for outlier samples with very low or high depth".to_string(),
                severity: None,
            },
        ],
        LibrarySizeCategory::Extreme => vec![
            Implication {
                code: "NORMALIZATION_CRITICAL".to_string(),
                message: "Must normalize; consider investigating source of variation".to_string(),
                severity: Some(Severity::High),
            },
            Implication {
                code: "FILTER_LOW_DEPTH".to_string(),
                message: "Consider filtering very low-depth samples".to_string(),
                severity: Some(Severity::Moderate),
            },
        ],
    }
}

fn build_compositionality_analysis(profile: &CompositionalityProfile) -> CompositionalityAnalysis {
    let dominance = DominanceMetrics {
        top_1_percent: profile.dominance_top_1 * 100.0,
        top_3_percent: profile.dominance_top_3 * 100.0,
        top_10_percent: profile.dominance_top_10 * 100.0,
    };

    let dominant_features: Vec<DominantFeatureInfo> = profile
        .dominant_features
        .iter()
        .map(|f| DominantFeatureInfo {
            id: f.feature_id.clone(),
            mean_proportion: f.mean_proportion,
        })
        .collect();

    let category = profile.dominance_category();
    let category_str = category.as_str().to_string();

    let context = match category {
        DominanceCategory::Low => {
            "Low dominance means compositional effects are minimal. The community is \
             relatively even, so changes in one taxon don't strongly affect the apparent \
             abundance of others."
        }
        DominanceCategory::Moderate => {
            "Moderate dominance indicates some compositional bias is possible. When a \
             dominant taxon changes, it may create apparent changes in other taxa even \
             if their absolute abundance is unchanged."
        }
        DominanceCategory::High => {
            "High dominance means a few taxa comprise most of the reads. Changes in \
             dominant taxa can significantly affect relative abundances of other taxa, \
             potentially creating spurious associations or masking real effects."
        }
        DominanceCategory::Extreme => {
            "Extreme dominance is a serious concern. One or a few taxa dominate the \
             community, meaning most statistical signals may be driven by compositional \
             effects rather than true biological changes. Results require very careful \
             interpretation."
        }
    }
    .to_string();

    let mut implications = Vec::new();
    if category == DominanceCategory::High || category == DominanceCategory::Extreme {
        implications.push(Implication {
            code: "COMPOSITIONAL_ARTIFACTS".to_string(),
            message: "High dominance may attenuate true signals or create false associations"
                .to_string(),
            severity: Some(Severity::Moderate),
        });
    }

    CompositionalityAnalysis {
        dominance,
        evenness: profile.evenness,
        gini: profile.gini,
        category: category_str,
        dominant_features,
        context,
        implications,
    }
}

fn generate_recommendations(
    sparsity: &SparsityAnalysis,
    _prevalence: &PrevalenceAnalysis,
    _library: &LibrarySizeAnalysis,
    _comp: &CompositionalityAnalysis,
    group_specific: &GroupSpecificFeatures,
) -> Recommendations {
    // Filtering recommendation
    let (filter_rec, filter_conf, filter_rationale, filter_alts) = if group_specific.percent >= 15.0
    {
        (
            "5% threshold with groupwise 'any' logic".to_string(),
            Confidence::High,
            "With significant group-specific features, groupwise filtering with 'any' logic \
             preserves features present in only one group while still removing noise."
                .to_string(),
            vec![Alternative {
                option: "10% overall".to_string(),
                when: "Want simpler, more conservative analysis".to_string(),
            }],
        )
    } else {
        (
            "5-10% overall prevalence threshold".to_string(),
            Confidence::Medium,
            "Few group-specific features, so overall prevalence filtering is appropriate. \
             5-10% balances noise reduction with signal preservation."
                .to_string(),
            vec![
                Alternative {
                    option: "Groupwise filtering".to_string(),
                    when: "Suspect presence/absence differences are biologically important"
                        .to_string(),
                },
                Alternative {
                    option: "Run optimize-prevalence".to_string(),
                    when: "Want data-driven threshold selection".to_string(),
                },
            ],
        )
    };

    let filtering = FilteringRecommendation {
        prevalence: ComponentRecommendation {
            recommendation: filter_rec,
            confidence: filter_conf,
            rationale: filter_rationale,
            alternatives: filter_alts,
        },
        abundance: None,
    };

    // Zero handling recommendation
    let (zero_rec, zero_conf, zero_rationale) = match sparsity.category {
        SparsityCategory::Low | SparsityCategory::Moderate => (
            "Pseudocount 0.5".to_string(),
            Confidence::High,
            "Standard pseudocount is sufficient for this sparsity level.".to_string(),
        ),
        SparsityCategory::High => (
            "Pseudocount 0.5, or consider ZINB".to_string(),
            Confidence::Medium,
            "High sparsity may benefit from zero-inflated models that distinguish \
             structural from sampling zeros."
                .to_string(),
        ),
        SparsityCategory::Extreme => (
            "ZINB recommended, or presence/absence testing".to_string(),
            Confidence::Medium,
            "Extreme sparsity suggests many structural zeros. ZINB explicitly models \
             the zero-inflation process."
                .to_string(),
        ),
    };

    let zero_handling = ComponentRecommendation {
        recommendation: zero_rec,
        confidence: zero_conf,
        rationale: zero_rationale,
        alternatives: vec![Alternative {
            option: "Different pseudocount (0.1, 1.0)".to_string(),
            when: "Want to assess sensitivity to pseudocount choice".to_string(),
        }],
    };

    // Normalization recommendation
    let normalization = ComponentRecommendation {
        recommendation: "CLR (Centered Log-Ratio)".to_string(),
        confidence: Confidence::High,
        rationale: "CLR is the standard for compositional data, handling the unit-sum \
                    constraint while preserving relative relationships."
            .to_string(),
        alternatives: vec![
            Alternative {
                option: "TSS (Total Sum Scaling)".to_string(),
                when: "Need simple relative abundances without log transformation".to_string(),
            },
            Alternative {
                option: "No normalization with ZINB".to_string(),
                when: "Using count-based models that handle library size internally".to_string(),
            },
        ],
    };

    // Model recommendation
    let (model_rec, model_conf, model_rationale, model_alts) = match sparsity.category {
        SparsityCategory::Low | SparsityCategory::Moderate => (
            "Linear model on CLR-transformed data (LinDA approach)".to_string(),
            Confidence::High,
            "With moderate sparsity, linear models on CLR data perform well and are \
             computationally efficient."
                .to_string(),
            vec![Alternative {
                option: "ZINB".to_string(),
                when: "Strong biological reason to expect zero-inflation".to_string(),
            }],
        ),
        SparsityCategory::High | SparsityCategory::Extreme => (
            "ZINB or linear model with careful validation".to_string(),
            Confidence::Medium,
            "High sparsity may benefit from zero-inflated models, but linear models \
             can still work well with appropriate filtering."
                .to_string(),
            vec![
                Alternative {
                    option: "Linear model on CLR".to_string(),
                    when: "Prefer interpretability and speed".to_string(),
                },
                Alternative {
                    option: "Presence/absence logistic".to_string(),
                    when: "Binary presence is the main signal".to_string(),
                },
            ],
        ),
    };

    let model = ComponentRecommendation {
        recommendation: model_rec,
        confidence: model_conf,
        rationale: model_rationale,
        alternatives: model_alts,
    };

    Recommendations {
        filtering,
        zero_handling,
        normalization,
        model,
    }
}

fn generate_warnings(
    sparsity: &SparsityAnalysis,
    library: &LibrarySizeAnalysis,
    comp: &CompositionalityAnalysis,
    groups: &GroupInfo,
) -> Vec<Warning> {
    let mut warnings = Vec::new();

    // Group imbalance warning
    if !groups.balanced {
        let sizes: Vec<usize> = groups.sizes.values().cloned().collect();
        let max = sizes.iter().max().copied().unwrap_or(1);
        let min = sizes.iter().min().copied().unwrap_or(1);
        let ratio = max as f64 / min as f64;

        warnings.push(Warning {
            severity: if ratio > 2.0 {
                Severity::Moderate
            } else {
                Severity::Low
            },
            category: "design".to_string(),
            code: "GROUP_IMBALANCE".to_string(),
            message: format!("Unbalanced groups ({:.1}x ratio)", ratio),
            details: "Imbalanced group sizes can reduce statistical power and may affect \
                      some test assumptions."
                .to_string(),
            action: "Consider balanced sampling in future experiments".to_string(),
        });
    }

    // Library size imbalance between groups
    if let Some(ref by_group) = library.by_group {
        if by_group.ratio > 1.5 {
            warnings.push(Warning {
                severity: if by_group.ratio > 2.0 {
                    Severity::Moderate
                } else {
                    Severity::Low
                },
                category: "confounding".to_string(),
                code: "LIBRARY_SIZE_IMBALANCE".to_string(),
                message: format!(
                    "Library size differs between groups ({:.1}x)",
                    by_group.ratio
                ),
                details: "Systematic library size differences between groups can confound \
                          differential abundance results if not properly normalized."
                    .to_string(),
                action: "Ensure normalization is applied; consider including library size \
                         as a covariate."
                    .to_string(),
            });
        }
    }

    // Extreme sparsity
    if sparsity.category == SparsityCategory::Extreme {
        warnings.push(Warning {
            severity: Severity::Moderate,
            category: "data_quality".to_string(),
            code: "EXTREME_SPARSITY".to_string(),
            message: format!("Extreme sparsity ({:.0}% zeros)", sparsity.overall * 100.0),
            details: "Very high sparsity limits the statistical power for abundance-based \
                      testing and may indicate shallow sequencing."
                .to_string(),
            action: "Consider aggressive filtering or presence/absence analysis".to_string(),
        });
    }

    // High dominance
    if comp.category == "high" || comp.category == "extreme" {
        warnings.push(Warning {
            severity: Severity::Moderate,
            category: "compositional".to_string(),
            code: "HIGH_DOMINANCE".to_string(),
            message: format!(
                "Top 3 features comprise {:.0}% of reads",
                comp.dominance.top_3_percent
            ),
            details: "High dominance by a few features can cause compositional artifacts \
                      where true signals are attenuated or spurious signals appear."
                .to_string(),
            action: "Run spike-in validation to quantify effect attenuation".to_string(),
        });
    }

    warnings
}

fn generate_validation_suggestions(group_column: &str) -> ValidationSuggestions {
    ValidationSuggestions {
        description: "Validation ensures your pipeline performs as expected on your specific \
                      data. Spike-in tests inject known effects and measure detection accuracy."
            .to_string(),
        essential: vec![
            ValidationStep {
                name: "Spike-in validation".to_string(),
                purpose: "Verify detection of 2-fold changes at controlled FDR".to_string(),
                command: format!(
                    "daa validate -c counts.tsv -m metadata.tsv -g {} -t treatment \
                     -f \"~ {}\" --test-coef {}treatment --fold-change 2.0",
                    group_column, group_column, group_column
                ),
                expected: "Sensitivity >70%, FDR <15%".to_string(),
                interpretation: "If sensitivity is low, consider adjusting filtering or \
                                 trying a different model. If FDR is high, increase the \
                                 significance threshold."
                    .to_string(),
            },
            ValidationStep {
                name: "Prevalence optimization".to_string(),
                purpose: "Find optimal prevalence threshold for your data".to_string(),
                command: format!(
                    "daa optimize-prevalence -c counts.tsv -m metadata.tsv -g {} \
                     -t treatment -f \"~ {}\" --test-coef {}treatment",
                    group_column, group_column, group_column
                ),
                expected: "Recommended threshold with F1 score".to_string(),
                interpretation: "Use the recommended threshold, or compare results at \
                                 different thresholds to assess robustness."
                    .to_string(),
            },
        ],
        recommended: vec![ValidationStep {
            name: "Compositional stress test".to_string(),
            purpose: "Quantify compositional effect attenuation".to_string(),
            command: format!(
                "daa stress -c counts.tsv -m metadata.tsv -g {} -t treatment \
                 -f \"~ {}\" --test-coef {}treatment --quick",
                group_column, group_column, group_column
            ),
            expected: "Effect attenuation <30%".to_string(),
            interpretation: "High attenuation indicates compositional effects are \
                             substantial; interpret fold-changes cautiously."
                .to_string(),
        }],
    }
}

fn generate_pipeline_suggestions(
    recommendations: &Recommendations,
    group_column: &str,
) -> SuggestedPipelines {
    // Build recommended pipeline config
    let recommended_config = format!(
        r#"name: recommended_pipeline
steps:
  - filter_prevalence_groupwise:
      threshold: 0.05
      group_column: "{}"
      logic: any
  - add_pseudocount:
      value: 0.5
  - normalize_clr: {{}}
  - model_lm:
      formula: "~ {}"
  - test_wald:
      coefficient: "{}treatment"
  - correct_bh:
      alpha: 0.05"#,
        group_column, group_column, group_column
    );

    let recommended = PipelineConfig {
        name: "recommended_pipeline".to_string(),
        description: format!(
            "Based on your data profile: {}. This pipeline uses groupwise filtering \
             to preserve group-specific signals, CLR normalization for compositional \
             data, and linear modeling for efficiency.",
            recommendations.filtering.prevalence.rationale
        ),
        config: recommended_config,
    };

    // ZINB alternative
    let zinb_config = format!(
        r#"name: zinb_pipeline
steps:
  - filter_prevalence_groupwise:
      threshold: 0.05
      group_column: "{}"
      logic: any
  - model_zinb:
      formula: "~ {}"
      max_iter: 100
  - test_wald_zinb:
      coefficient: "{}treatment"
  - correct_bh:
      alpha: 0.05"#,
        group_column, group_column, group_column
    );

    let zinb_alt = PipelineConfig {
        name: "zinb_pipeline".to_string(),
        description: "Zero-inflated negative binomial model for data with structural zeros. \
                      More computationally intensive but may better capture zero-inflation."
            .to_string(),
        config: zinb_config,
    };

    SuggestedPipelines {
        recommended,
        alternatives: vec![zinb_alt],
    }
}

fn generate_llm_summary(
    sparsity: &SparsityAnalysis,
    prevalence: &PrevalenceAnalysis,
    library: &LibrarySizeAnalysis,
    comp: &CompositionalityAnalysis,
    recommendations: &Recommendations,
    warnings: &Vec<Warning>,
) -> String {
    let mut summary = String::new();

    summary.push_str("## Data Overview\n\n");
    summary.push_str(&format!(
        "This is a microbiome count matrix with {} sparsity ({:.0}% zeros), {} \
         library size variation (CV={:.2}), and {} dominance ({:.0}% from top 3 features).\n\n",
        sparsity.category.as_str(),
        sparsity.overall * 100.0,
        library.category.as_str(),
        library.cv,
        comp.category,
        comp.dominance.top_3_percent
    ));

    summary.push_str("## Key Characteristics\n\n");
    summary.push_str(&format!(
        "- Sparsity: {} - {}\n",
        sparsity.category.as_str(),
        sparsity.context.lines().next().unwrap_or("")
    ));
    summary.push_str(&format!(
        "- Group-specific features: {:.1}% - {}\n",
        prevalence.group_specific.percent,
        prevalence.group_specific.context.lines().next().unwrap_or("")
    ));
    summary.push_str(&format!(
        "- Library size: {} ({})\n",
        library.category.as_str(),
        library.by_group.as_ref().map_or(
            "".to_string(),
            |g| g.interpretation.clone()
        )
    ));
    summary.push_str(&format!(
        "- Compositionality: {} dominance, evenness={:.2}\n\n",
        comp.category, comp.evenness
    ));

    summary.push_str("## Recommended Approach\n\n");
    summary.push_str(&format!(
        "1. Filter: {}\n",
        recommendations.filtering.prevalence.recommendation
    ));
    summary.push_str(&format!(
        "2. Zeros: {}\n",
        recommendations.zero_handling.recommendation
    ));
    summary.push_str(&format!(
        "3. Normalize: {}\n",
        recommendations.normalization.recommendation
    ));
    summary.push_str(&format!("4. Model: {}\n\n", recommendations.model.recommendation));

    summary.push_str("## Critical Validations\n\n");
    summary.push_str("- Run spike-in validation to verify detection performance\n");
    summary.push_str("- Run optimize-prevalence to find data-driven threshold\n\n");

    if !warnings.is_empty() {
        summary.push_str("## Warnings to Communicate\n\n");
        for warning in warnings {
            summary.push_str(&format!(
                "- [{}] {}: {}\n",
                warning.severity.as_str(),
                warning.code,
                warning.message
            ));
        }
    }

    summary
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_counts() -> CountMatrix {
        let mut tri_mat = TriMat::new((20, 10));
        // Add some counts with varying sparsity
        for feat in 0..20 {
            for sample in 0..10 {
                if (feat + sample) % 3 != 0 {
                    tri_mat.add_triplet(feat, sample, ((feat + 1) * 10) as u64);
                }
            }
        }
        let feature_ids: Vec<String> = (0..20).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..10).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        use std::io::Write;
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup").unwrap();
        for i in 0..10 {
            let group = if i < 5 { "control" } else { "treatment" };
            writeln!(file, "S{}\t{}", i, group).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_sparsity_category() {
        assert_eq!(SparsityCategory::from_value(0.3), SparsityCategory::Low);
        assert_eq!(SparsityCategory::from_value(0.6), SparsityCategory::Moderate);
        assert_eq!(SparsityCategory::from_value(0.8), SparsityCategory::High);
        assert_eq!(SparsityCategory::from_value(0.95), SparsityCategory::Extreme);
    }

    #[test]
    fn test_library_size_category() {
        assert_eq!(LibrarySizeCategory::from_cv(0.3), LibrarySizeCategory::Low);
        assert_eq!(
            LibrarySizeCategory::from_cv(0.7),
            LibrarySizeCategory::Moderate
        );
        assert_eq!(LibrarySizeCategory::from_cv(1.5), LibrarySizeCategory::High);
        assert_eq!(
            LibrarySizeCategory::from_cv(2.5),
            LibrarySizeCategory::Extreme
        );
    }

    #[test]
    fn test_profile_generation() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let profile = profile_for_llm(&counts, &metadata, "group").unwrap();

        assert_eq!(profile.data_summary.features, 20);
        assert_eq!(profile.data_summary.samples, 10);
        assert_eq!(profile.data_summary.groups.levels.len(), 2);
        assert!(!profile.llm_summary.is_empty());
    }

    #[test]
    fn test_yaml_serialization() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let profile = profile_for_llm(&counts, &metadata, "group").unwrap();
        let yaml = profile.to_yaml().unwrap();

        assert!(yaml.contains("meta:"));
        assert!(yaml.contains("data_summary:"));
        assert!(yaml.contains("sparsity:"));
        assert!(yaml.contains("recommendations:"));
    }

    #[test]
    fn test_markdown_output() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let profile = profile_for_llm(&counts, &metadata, "group").unwrap();
        let md = profile.to_markdown();

        assert!(md.contains("# Data Profile"));
        assert!(md.contains("## Data Summary"));
        assert!(md.contains("## Recommendations"));
    }

    #[test]
    fn test_extract_group_info() {
        let metadata = create_test_metadata();
        let info = extract_group_info(&metadata, "group").unwrap();

        assert_eq!(info.levels.len(), 2);
        assert!(info.balanced);
        assert_eq!(info.sizes.get("control"), Some(&5));
        assert_eq!(info.sizes.get("treatment"), Some(&5));
    }
}
