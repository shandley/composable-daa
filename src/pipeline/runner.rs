//! Pipeline runner for composing and executing analysis steps.

use crate::correct::{bh::correct_bh_wald, bh::create_results};
use crate::data::{CountMatrix, DaResultSet, DesignMatrix, Formula, Metadata};
use crate::error::{DaaError, Result};
use crate::filter::{
    filter_abundance, filter_library_size, filter_library_size_quantile,
    filter_library_size_relative, filter_mean_abundance, filter_prevalence_groupwise,
    filter_prevalence_overall, filter_stratified, filter_total_count, GroupwiseLogic,
    TierThresholds,
};
use crate::model::model_lm;
use crate::normalize::{norm_clr, norm_tss, TransformedMatrix};
use crate::profile::{profile_prevalence, PrevalenceProfile};
use crate::test::test_wald;
use crate::zero::pseudocount::add_pseudocount;
use serde::{Deserialize, Serialize};

/// A step in the analysis pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PipelineStep {
    // === Feature Filtering ===
    /// Filter by overall prevalence threshold.
    FilterPrevalence { threshold: f64 },
    /// Filter by group-wise prevalence.
    FilterPrevalenceGroupwise {
        threshold: f64,
        group_column: String,
        logic: GroupwiseLogic,
    },
    /// Filter by relative abundance range.
    FilterAbundance {
        min_abundance: f64,
        max_abundance: Option<f64>,
    },
    /// Filter by mean abundance when present.
    FilterMeanAbundance { min_mean: f64 },
    /// Filter by total count across samples.
    FilterTotalCount { min_total: u64 },
    /// Filter using stratified thresholds by prevalence tier.
    FilterStratified { preset: StratifiedPreset },

    // === Sample Filtering ===
    /// Filter samples by library size (total reads).
    FilterLibrarySize {
        min_reads: Option<u64>,
        max_reads: Option<u64>,
    },
    /// Filter samples by library size quantile.
    FilterLibrarySizeQuantile {
        lower_quantile: f64,
        upper_quantile: f64,
    },
    /// Filter samples by library size relative to median.
    FilterLibrarySizeRelative { min_fraction_of_median: f64 },

    // === Zero Handling ===
    /// Add pseudocount for zero handling.
    AddPseudocount { value: f64 },

    // === Normalization ===
    /// Apply CLR normalization.
    NormalizeCLR,
    /// Apply TSS normalization (relative abundances).
    NormalizeTSS { scale_factor: f64 },

    // === Model Fitting ===
    /// Fit linear model.
    ModelLM { formula: String },

    // === Testing ===
    /// Wald test for a coefficient.
    TestWald { coefficient: String },

    // === Multiple Testing Correction ===
    /// Benjamini-Hochberg correction.
    CorrectBH,
}

/// Preset configurations for stratified filtering.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum StratifiedPreset {
    /// Very permissive, keeps most features.
    Permissive,
    /// Lenient on rare features (higher abundance requirement for rare).
    LenientRare,
    /// Strict filtering (removes most rare features).
    Strict,
    /// Custom thresholds (use with `filter_stratified_custom`).
    Custom,
}

/// Pipeline configuration for serialization.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    /// Name of the pipeline.
    pub name: String,
    /// Description.
    pub description: Option<String>,
    /// Steps to execute.
    pub steps: Vec<PipelineStep>,
}

impl PipelineConfig {
    /// Load from YAML string.
    pub fn from_yaml(yaml: &str) -> Result<Self> {
        serde_yaml::from_str(yaml).map_err(DaaError::from)
    }

    /// Save to YAML string.
    pub fn to_yaml(&self) -> Result<String> {
        serde_yaml::to_string(self).map_err(DaaError::from)
    }
}

/// Builder for constructing and running analysis pipelines.
#[derive(Debug, Clone)]
pub struct Pipeline {
    steps: Vec<PipelineStep>,
    name: String,
}

impl Default for Pipeline {
    fn default() -> Self {
        Self::new()
    }
}

impl Pipeline {
    /// Create a new empty pipeline.
    pub fn new() -> Self {
        Self {
            steps: Vec::new(),
            name: "unnamed".to_string(),
        }
    }

    /// Create from a config.
    pub fn from_config(config: &PipelineConfig) -> Self {
        Self {
            steps: config.steps.clone(),
            name: config.name.clone(),
        }
    }

    /// Set the pipeline name.
    pub fn name(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }

    /// Add overall prevalence filtering.
    pub fn filter_prevalence(mut self, threshold: f64) -> Self {
        self.steps.push(PipelineStep::FilterPrevalence { threshold });
        self
    }

    /// Add group-wise prevalence filtering.
    pub fn filter_prevalence_groupwise(
        mut self,
        threshold: f64,
        group_column: &str,
        logic: GroupwiseLogic,
    ) -> Self {
        self.steps.push(PipelineStep::FilterPrevalenceGroupwise {
            threshold,
            group_column: group_column.to_string(),
            logic,
        });
        self
    }

    /// Filter features by relative abundance range.
    ///
    /// Keeps features whose mean relative abundance falls within the specified range.
    pub fn filter_abundance(mut self, min_abundance: f64, max_abundance: Option<f64>) -> Self {
        self.steps.push(PipelineStep::FilterAbundance {
            min_abundance,
            max_abundance,
        });
        self
    }

    /// Filter features by mean abundance when present.
    ///
    /// Keeps features whose mean count (when non-zero) meets the threshold.
    pub fn filter_mean_abundance(mut self, min_mean: f64) -> Self {
        self.steps
            .push(PipelineStep::FilterMeanAbundance { min_mean });
        self
    }

    /// Filter features by total count across all samples.
    pub fn filter_total_count(mut self, min_total: u64) -> Self {
        self.steps
            .push(PipelineStep::FilterTotalCount { min_total });
        self
    }

    /// Filter features using stratified thresholds by prevalence tier.
    ///
    /// Uses preset configurations:
    /// - `Permissive`: Very lenient, keeps most features
    /// - `LenientRare`: Higher requirements for rare features
    /// - `Strict`: Removes most rare and low-abundance features
    pub fn filter_stratified(mut self, preset: StratifiedPreset) -> Self {
        self.steps.push(PipelineStep::FilterStratified { preset });
        self
    }

    /// Filter samples by library size (total reads).
    ///
    /// Removes samples with too few or too many total counts.
    pub fn filter_library_size(mut self, min_reads: Option<u64>, max_reads: Option<u64>) -> Self {
        self.steps.push(PipelineStep::FilterLibrarySize {
            min_reads,
            max_reads,
        });
        self
    }

    /// Filter samples by library size quantile.
    ///
    /// Removes samples in the lower or upper quantiles of the library size distribution.
    pub fn filter_library_size_quantile(
        mut self,
        lower_quantile: f64,
        upper_quantile: f64,
    ) -> Self {
        self.steps.push(PipelineStep::FilterLibrarySizeQuantile {
            lower_quantile,
            upper_quantile,
        });
        self
    }

    /// Filter samples by library size relative to median.
    ///
    /// Removes samples with library size below a fraction of the median.
    pub fn filter_library_size_relative(mut self, min_fraction_of_median: f64) -> Self {
        self.steps
            .push(PipelineStep::FilterLibrarySizeRelative {
                min_fraction_of_median,
            });
        self
    }

    /// Add pseudocount.
    pub fn add_pseudocount(mut self, value: f64) -> Self {
        self.steps.push(PipelineStep::AddPseudocount { value });
        self
    }

    /// Add CLR normalization.
    pub fn normalize_clr(mut self) -> Self {
        self.steps.push(PipelineStep::NormalizeCLR);
        self
    }

    /// Add TSS normalization (relative abundances).
    ///
    /// TSS divides each count by the sample total, converting to proportions.
    /// A scale factor can be applied (1.0 for proportions, 1e6 for CPM).
    pub fn normalize_tss(mut self, scale_factor: f64) -> Self {
        self.steps.push(PipelineStep::NormalizeTSS { scale_factor });
        self
    }

    /// Add TSS normalization as proportions (scale_factor = 1.0).
    pub fn normalize_tss_proportions(self) -> Self {
        self.normalize_tss(1.0)
    }

    /// Add TSS normalization as CPM (counts per million).
    pub fn normalize_tss_cpm(self) -> Self {
        self.normalize_tss(1_000_000.0)
    }

    /// Add linear model.
    pub fn model_lm(mut self, formula: &str) -> Self {
        self.steps.push(PipelineStep::ModelLM {
            formula: formula.to_string(),
        });
        self
    }

    /// Add Wald test.
    pub fn test_wald(mut self, coefficient: &str) -> Self {
        self.steps.push(PipelineStep::TestWald {
            coefficient: coefficient.to_string(),
        });
        self
    }

    /// Add BH correction.
    pub fn correct_bh(mut self) -> Self {
        self.steps.push(PipelineStep::CorrectBH);
        self
    }

    /// Convert to config for serialization.
    pub fn to_config(&self, description: Option<&str>) -> PipelineConfig {
        PipelineConfig {
            name: self.name.clone(),
            description: description.map(String::from),
            steps: self.steps.clone(),
        }
    }

    /// Run the pipeline on data.
    pub fn run(&self, counts: &CountMatrix, metadata: &Metadata) -> Result<DaResultSet> {
        let mut state = PipelineState::new(counts.clone(), metadata.clone());

        for (i, step) in self.steps.iter().enumerate() {
            state = state.apply(step).map_err(|e| {
                DaaError::Pipeline(format!("Step {} ({:?}) failed: {}", i + 1, step, e))
            })?;
        }

        state.finalize(&self.name)
    }
}

/// Internal state during pipeline execution.
struct PipelineState {
    counts: CountMatrix,
    metadata: Metadata,
    pseudocount_data: Option<nalgebra::DMatrix<f64>>,
    transformed: Option<TransformedMatrix>,
    design: Option<DesignMatrix>,
    lm_fit: Option<crate::model::LmFit>,
    wald_result: Option<crate::test::WaldResult>,
    bh_corrected: Option<crate::correct::bh::BhCorrected>,
    prevalence: Option<PrevalenceProfile>,
}

impl PipelineState {
    fn new(counts: CountMatrix, metadata: Metadata) -> Self {
        Self {
            counts,
            metadata,
            pseudocount_data: None,
            transformed: None,
            design: None,
            lm_fit: None,
            wald_result: None,
            bh_corrected: None,
            prevalence: None,
        }
    }

    fn apply(mut self, step: &PipelineStep) -> Result<Self> {
        match step {
            // === Feature Filtering ===
            PipelineStep::FilterPrevalence { threshold } => {
                self.counts = filter_prevalence_overall(&self.counts, *threshold)?;
            }
            PipelineStep::FilterPrevalenceGroupwise {
                threshold,
                group_column,
                logic,
            } => {
                self.counts = filter_prevalence_groupwise(
                    &self.counts,
                    &self.metadata,
                    group_column,
                    *threshold,
                    *logic,
                )?;
            }
            PipelineStep::FilterAbundance {
                min_abundance,
                max_abundance,
            } => {
                self.counts = filter_abundance(&self.counts, *min_abundance, *max_abundance)?;
            }
            PipelineStep::FilterMeanAbundance { min_mean } => {
                self.counts = filter_mean_abundance(&self.counts, *min_mean)?;
            }
            PipelineStep::FilterTotalCount { min_total } => {
                self.counts = filter_total_count(&self.counts, *min_total)?;
            }
            PipelineStep::FilterStratified { preset } => {
                let thresholds = match preset {
                    StratifiedPreset::Permissive => TierThresholds::new(),
                    StratifiedPreset::LenientRare => TierThresholds::lenient_rare(),
                    StratifiedPreset::Strict => TierThresholds::strict(),
                    StratifiedPreset::Custom => {
                        return Err(DaaError::Pipeline(
                            "Custom stratified filtering requires filter_stratified_custom()"
                                .to_string(),
                        ));
                    }
                };
                self.counts = filter_stratified(&self.counts, &thresholds)?;
            }

            // === Sample Filtering ===
            PipelineStep::FilterLibrarySize {
                min_reads,
                max_reads,
            } => {
                self.counts = filter_library_size(&self.counts, *min_reads, *max_reads)?;
                // Also filter metadata to keep in sync
                self.metadata = self.metadata.subset_samples(self.counts.sample_ids())?;
            }
            PipelineStep::FilterLibrarySizeQuantile {
                lower_quantile,
                upper_quantile,
            } => {
                self.counts =
                    filter_library_size_quantile(&self.counts, *lower_quantile, *upper_quantile)?;
                self.metadata = self.metadata.subset_samples(self.counts.sample_ids())?;
            }
            PipelineStep::FilterLibrarySizeRelative {
                min_fraction_of_median,
            } => {
                self.counts =
                    filter_library_size_relative(&self.counts, *min_fraction_of_median)?;
                self.metadata = self.metadata.subset_samples(self.counts.sample_ids())?;
            }

            // === Zero Handling ===
            PipelineStep::AddPseudocount { value } => {
                self.pseudocount_data = Some(add_pseudocount(&self.counts, *value)?);
            }

            // === Normalization ===
            PipelineStep::NormalizeCLR => {
                let data = self.pseudocount_data.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must add pseudocount before CLR".to_string())
                })?;
                self.transformed = Some(norm_clr(
                    data,
                    self.counts.feature_ids().to_vec(),
                    self.counts.sample_ids().to_vec(),
                )?);
                // Store prevalence profile for results
                self.prevalence = Some(profile_prevalence(&self.counts));
            }

            PipelineStep::NormalizeTSS { scale_factor } => {
                // TSS works directly on counts (no pseudocount required)
                let tss_result = norm_tss(&self.counts, *scale_factor)?;
                self.transformed = Some(tss_result.to_transformed());
                // Store prevalence profile for results
                self.prevalence = Some(profile_prevalence(&self.counts));
            }

            // === Model Fitting ===
            PipelineStep::ModelLM { formula } => {
                let parsed_formula = Formula::parse(formula)?;
                self.design = Some(DesignMatrix::from_formula(&self.metadata, &parsed_formula)?);
                let transformed = self.transformed.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must normalize before fitting model".to_string())
                })?;
                let design = self.design.as_ref().unwrap();
                self.lm_fit = Some(model_lm(transformed, design)?);
            }

            // === Testing ===
            PipelineStep::TestWald { coefficient } => {
                let fit = self.lm_fit.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must fit model before Wald test".to_string())
                })?;
                self.wald_result = Some(test_wald(fit, coefficient)?);
            }

            // === Correction ===
            PipelineStep::CorrectBH => {
                let wald = self.wald_result.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must perform Wald test before BH correction".to_string())
                })?;
                self.bh_corrected = Some(correct_bh_wald(wald));
            }
        }
        Ok(self)
    }

    fn finalize(self, method_name: &str) -> Result<DaResultSet> {
        let wald = self.wald_result.ok_or_else(|| {
            DaaError::Pipeline("Pipeline must include a test step".to_string())
        })?;
        let bh = self.bh_corrected.ok_or_else(|| {
            DaaError::Pipeline("Pipeline must include correction step".to_string())
        })?;
        let prevalence = self.prevalence.ok_or_else(|| {
            DaaError::Pipeline("Prevalence not computed".to_string())
        })?;

        // Calculate mean abundances from transformed data
        let transformed = self.transformed.ok_or_else(|| {
            DaaError::Pipeline("Transformed data not available".to_string())
        })?;
        let mean_abundances: Vec<f64> = (0..transformed.n_features())
            .map(|i| {
                let row = transformed.row(i);
                row.iter().sum::<f64>() / row.len() as f64
            })
            .collect();

        Ok(create_results(
            &wald,
            &bh,
            &prevalence,
            &mean_abundances,
            method_name,
        ))
    }
}

/// Convenience function to run a LinDA-style pipeline.
pub fn run_linda(
    counts: &CountMatrix,
    metadata: &Metadata,
    formula: &str,
    coefficient: &str,
    prevalence_threshold: f64,
    pseudocount: f64,
) -> Result<DaResultSet> {
    Pipeline::new()
        .name("LinDA")
        .filter_prevalence(prevalence_threshold)
        .add_pseudocount(pseudocount)
        .normalize_clr()
        .model_lm(formula)
        .test_wald(coefficient)
        .correct_bh()
        .run(counts, metadata)
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_counts() -> CountMatrix {
        // 10 features × 8 samples
        let mut tri_mat = TriMat::new((10, 8));

        // Features 0-4: present in all samples, varying by group
        for feat in 0usize..5 {
            for sample in 0usize..8 {
                let is_treatment = sample % 2 == 1;
                let base = 100u64 + feat as u64 * 10;
                let effect = if feat == 0 { 0u64 } else { feat as u64 * 20 };
                let value = if is_treatment {
                    base + effect + (sample as u64 * 5)
                } else {
                    base + (sample as u64 * 3)
                };
                tri_mat.add_triplet(feat, sample, value);
            }
        }

        // Features 5-7: present in 50% of samples
        for feat in 5usize..8 {
            for sample in 0usize..4 {
                tri_mat.add_triplet(feat, sample, 50 + (feat as u64 * 10));
            }
        }

        // Features 8-9: rare (present in 2 samples only)
        tri_mat.add_triplet(8, 0, 20);
        tri_mat.add_triplet(8, 1, 25);
        tri_mat.add_triplet(9, 2, 15);
        tri_mat.add_triplet(9, 3, 18);

        let feature_ids: Vec<String> = (0..10).map(|i| format!("feat_{}", i)).collect();
        let sample_ids: Vec<String> = (0..8).map(|i| format!("S{}", i)).collect();
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        for i in 0..8 {
            let group = if i % 2 == 0 { "control" } else { "treatment" };
            let age = 25 + i * 2;
            writeln!(file, "S{}\t{}\t{}", i, group, age).unwrap();
        }
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_pipeline_builder() {
        let pipeline = Pipeline::new()
            .name("test")
            .filter_prevalence(0.1)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh();

        let config = pipeline.to_config(Some("Test pipeline"));
        assert_eq!(config.steps.len(), 6);
        assert_eq!(config.name, "test");
    }

    #[test]
    fn test_pipeline_run() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let results = Pipeline::new()
            .name("test")
            .filter_prevalence(0.3) // Keep features present in at least 30%
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh()
            .run(&counts, &metadata)
            .unwrap();

        // Should have results for filtered features
        assert!(results.len() > 0);
        assert!(results.len() < 10); // Some features should be filtered

        // Check that results have valid values
        for r in results.iter() {
            assert!(r.p_value >= 0.0 && r.p_value <= 1.0);
            assert!(r.q_value >= 0.0 && r.q_value <= 1.0);
        }
    }

    #[test]
    fn test_run_linda() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let results = run_linda(
            &counts,
            &metadata,
            "~ group",
            "grouptreatment",
            0.3,
            0.5,
        )
        .unwrap();

        assert!(results.len() > 0);
        assert_eq!(results.method, "LinDA");
    }

    #[test]
    fn test_pipeline_config_yaml() {
        let pipeline = Pipeline::new()
            .name("example")
            .filter_prevalence(0.1)
            .add_pseudocount(1.0)
            .normalize_clr()
            .model_lm("~ treatment")
            .test_wald("treatmentYes")
            .correct_bh();

        let config = pipeline.to_config(Some("Example LinDA pipeline"));
        let yaml = config.to_yaml().unwrap();

        // Verify it can be parsed back
        let parsed = PipelineConfig::from_yaml(&yaml).unwrap();
        assert_eq!(parsed.name, "example");
        assert_eq!(parsed.steps.len(), 6);
    }

    #[test]
    fn test_pipeline_error_handling() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Try to run model without normalization
        let result = Pipeline::new()
            .add_pseudocount(0.5)
            .model_lm("~ group")
            .run(&counts, &metadata);

        assert!(result.is_err());
    }

    #[test]
    fn test_pipeline_with_abundance_filter() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let results = Pipeline::new()
            .name("test-abundance")
            .filter_abundance(0.001, None) // Keep features with >= 0.1% abundance
            .filter_prevalence(0.3)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh()
            .run(&counts, &metadata)
            .unwrap();

        assert!(results.len() > 0);
    }

    #[test]
    fn test_pipeline_with_library_size_filter() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let results = Pipeline::new()
            .name("test-libsize")
            .filter_library_size(Some(100), None) // Minimum 100 reads per sample
            .filter_prevalence(0.3)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh()
            .run(&counts, &metadata)
            .unwrap();

        assert!(results.len() > 0);
    }

    #[test]
    fn test_pipeline_with_stratified_filter() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        let results = Pipeline::new()
            .name("test-stratified")
            .filter_stratified(StratifiedPreset::LenientRare)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh()
            .run(&counts, &metadata)
            .unwrap();

        assert!(results.len() > 0);
    }

    #[test]
    fn test_pipeline_with_multiple_filters() {
        let counts = create_test_counts();
        let metadata = create_test_metadata();

        // Apply multiple filters in sequence
        let results = Pipeline::new()
            .name("test-multi-filter")
            .filter_library_size(Some(50), None)
            .filter_prevalence(0.2)
            .filter_mean_abundance(5.0)
            .add_pseudocount(0.5)
            .normalize_clr()
            .model_lm("~ group")
            .test_wald("grouptreatment")
            .correct_bh()
            .run(&counts, &metadata)
            .unwrap();

        assert!(results.len() > 0);
    }

    #[test]
    fn test_pipeline_config_with_new_filters() {
        let pipeline = Pipeline::new()
            .name("advanced")
            .filter_library_size(Some(1000), None)
            .filter_abundance(0.001, Some(0.5))
            .filter_stratified(StratifiedPreset::Strict)
            .add_pseudocount(1.0)
            .normalize_clr()
            .model_lm("~ treatment")
            .test_wald("treatmentYes")
            .correct_bh();

        let config = pipeline.to_config(Some("Advanced filtering pipeline"));
        let yaml = config.to_yaml().unwrap();

        // Verify it can be parsed back
        let parsed = PipelineConfig::from_yaml(&yaml).unwrap();
        assert_eq!(parsed.name, "advanced");
        assert_eq!(parsed.steps.len(), 8);
    }
}
