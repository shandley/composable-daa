//! Pipeline runner for composing and executing analysis steps.

use crate::correct::{bh::correct_bh_wald, bh::create_results};
use crate::data::{CountMatrix, DaResultSet, DesignMatrix, Formula, Metadata};
use crate::error::{DaaError, Result};
use crate::filter::{filter_prevalence_groupwise, filter_prevalence_overall, GroupwiseLogic};
use crate::model::model_lm;
use crate::normalize::{norm_clr, TransformedMatrix};
use crate::profile::{profile_prevalence, PrevalenceProfile};
use crate::test::test_wald;
use crate::zero::pseudocount::add_pseudocount;
use serde::{Deserialize, Serialize};

/// A step in the analysis pipeline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PipelineStep {
    /// Filter by overall prevalence threshold.
    FilterPrevalence { threshold: f64 },
    /// Filter by group-wise prevalence.
    FilterPrevalenceGroupwise {
        threshold: f64,
        group_column: String,
        logic: GroupwiseLogic,
    },
    /// Add pseudocount for zero handling.
    AddPseudocount { value: f64 },
    /// Apply CLR normalization.
    NormalizeCLR,
    /// Fit linear model.
    ModelLM { formula: String },
    /// Wald test for a coefficient.
    TestWald { coefficient: String },
    /// Benjamini-Hochberg correction.
    CorrectBH,
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
            PipelineStep::AddPseudocount { value } => {
                self.pseudocount_data = Some(add_pseudocount(&self.counts, *value)?);
            }
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
            PipelineStep::ModelLM { formula } => {
                let parsed_formula = Formula::parse(formula)?;
                self.design = Some(DesignMatrix::from_formula(&self.metadata, &parsed_formula)?);
                let transformed = self.transformed.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must normalize before fitting model".to_string())
                })?;
                let design = self.design.as_ref().unwrap();
                self.lm_fit = Some(model_lm(transformed, design)?);
            }
            PipelineStep::TestWald { coefficient } => {
                let fit = self.lm_fit.as_ref().ok_or_else(|| {
                    DaaError::Pipeline("Must fit model before Wald test".to_string())
                })?;
                self.wald_result = Some(test_wald(fit, coefficient)?);
            }
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
}
