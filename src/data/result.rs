//! Result types for differential abundance analysis.

use crate::error::Result;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Prevalence tier classification for features.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum PrevalenceTier {
    /// Very rare: prevalence < 5%
    VeryRare,
    /// Rare: prevalence 5-10%
    Rare,
    /// Low: prevalence 10-25%
    Low,
    /// Medium: prevalence 25-50%
    Medium,
    /// High: prevalence 50-75%
    High,
    /// VeryHigh: prevalence > 75%
    VeryHigh,
}

impl PrevalenceTier {
    /// Classify a prevalence value (0.0-1.0) into a tier.
    pub fn from_prevalence(prevalence: f64) -> Self {
        if prevalence < 0.05 {
            Self::VeryRare
        } else if prevalence < 0.10 {
            Self::Rare
        } else if prevalence < 0.25 {
            Self::Low
        } else if prevalence < 0.50 {
            Self::Medium
        } else if prevalence < 0.75 {
            Self::High
        } else {
            Self::VeryHigh
        }
    }

    /// Get the descriptive name.
    pub fn name(&self) -> &'static str {
        match self {
            Self::VeryRare => "very_rare",
            Self::Rare => "rare",
            Self::Low => "low",
            Self::Medium => "medium",
            Self::High => "high",
            Self::VeryHigh => "very_high",
        }
    }
}

/// Confidence level based on statistical evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Confidence {
    /// Very high confidence (q < 0.001)
    VeryHigh,
    /// High confidence (q < 0.01)
    High,
    /// Moderate confidence (q < 0.05)
    Moderate,
    /// Low confidence (q < 0.10)
    Low,
    /// Suggestive (q < 0.20)
    Suggestive,
    /// Not significant
    NotSignificant,
}

impl Confidence {
    /// Classify based on adjusted p-value (q-value).
    pub fn from_qvalue(q: f64) -> Self {
        if q < 0.001 {
            Self::VeryHigh
        } else if q < 0.01 {
            Self::High
        } else if q < 0.05 {
            Self::Moderate
        } else if q < 0.10 {
            Self::Low
        } else if q < 0.20 {
            Self::Suggestive
        } else {
            Self::NotSignificant
        }
    }

    /// Get the descriptive name.
    pub fn name(&self) -> &'static str {
        match self {
            Self::VeryHigh => "very_high",
            Self::High => "high",
            Self::Moderate => "moderate",
            Self::Low => "low",
            Self::Suggestive => "suggestive",
            Self::NotSignificant => "not_significant",
        }
    }

    /// Check if significant at alpha = 0.05.
    pub fn is_significant(&self) -> bool {
        matches!(self, Self::VeryHigh | Self::High | Self::Moderate)
    }
}

/// Result for a single feature from differential abundance analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DaResult {
    /// Feature identifier.
    pub feature_id: String,
    /// Name of the coefficient/contrast being tested.
    pub coefficient: String,
    /// Estimated effect size (log-fold change or similar).
    pub estimate: f64,
    /// Standard error of the estimate.
    pub std_error: f64,
    /// Test statistic (t, z, or Wald).
    pub statistic: f64,
    /// Raw p-value.
    pub p_value: f64,
    /// Adjusted p-value (q-value) after multiple testing correction.
    pub q_value: f64,
    /// Prevalence of this feature (proportion of samples with non-zero counts).
    pub prevalence: f64,
    /// Mean abundance (after any transformation).
    pub mean_abundance: f64,
    /// Prevalence tier classification.
    pub prevalence_tier: PrevalenceTier,
    /// Confidence level classification.
    pub confidence: Confidence,
}

impl DaResult {
    /// Create a new DA result.
    pub fn new(
        feature_id: String,
        coefficient: String,
        estimate: f64,
        std_error: f64,
        statistic: f64,
        p_value: f64,
        q_value: f64,
        prevalence: f64,
        mean_abundance: f64,
    ) -> Self {
        let prevalence_tier = PrevalenceTier::from_prevalence(prevalence);
        let confidence = Confidence::from_qvalue(q_value);
        Self {
            feature_id,
            coefficient,
            estimate,
            std_error,
            statistic,
            p_value,
            q_value,
            prevalence,
            mean_abundance,
            prevalence_tier,
            confidence,
        }
    }

    /// Check if this result is significant at the default threshold (q < 0.05).
    pub fn is_significant(&self) -> bool {
        self.confidence.is_significant()
    }

    /// Check if this result is significant at a custom threshold.
    pub fn is_significant_at(&self, alpha: f64) -> bool {
        self.q_value < alpha
    }
}

/// Collection of DA results.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DaResultSet {
    /// Method name used to generate these results.
    pub method: String,
    /// Individual results for each feature.
    pub results: Vec<DaResult>,
}

impl DaResultSet {
    /// Create a new result set.
    pub fn new(method: String, results: Vec<DaResult>) -> Self {
        Self { method, results }
    }

    /// Number of results.
    pub fn len(&self) -> usize {
        self.results.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.results.is_empty()
    }

    /// Get results sorted by p-value (ascending).
    pub fn sorted_by_pvalue(&self) -> Vec<&DaResult> {
        let mut sorted: Vec<_> = self.results.iter().collect();
        sorted.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());
        sorted
    }

    /// Get results sorted by q-value (ascending).
    pub fn sorted_by_qvalue(&self) -> Vec<&DaResult> {
        let mut sorted: Vec<_> = self.results.iter().collect();
        sorted.sort_by(|a, b| a.q_value.partial_cmp(&b.q_value).unwrap());
        sorted
    }

    /// Get significant results (q < 0.05).
    pub fn significant(&self) -> Vec<&DaResult> {
        self.results.iter().filter(|r| r.is_significant()).collect()
    }

    /// Get significant results at a custom threshold.
    pub fn significant_at(&self, alpha: f64) -> Vec<&DaResult> {
        self.results
            .iter()
            .filter(|r| r.is_significant_at(alpha))
            .collect()
    }

    /// Count significant results at various thresholds.
    pub fn summary(&self) -> ResultSummary {
        ResultSummary {
            total: self.len(),
            significant_001: self.results.iter().filter(|r| r.q_value < 0.001).count(),
            significant_01: self.results.iter().filter(|r| r.q_value < 0.01).count(),
            significant_05: self.results.iter().filter(|r| r.q_value < 0.05).count(),
            significant_10: self.results.iter().filter(|r| r.q_value < 0.10).count(),
        }
    }

    /// Write results to TSV file.
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Header
        writeln!(
            writer,
            "feature_id\tcoefficient\testimate\tstd_error\tstatistic\tp_value\tq_value\tprevalence\tmean_abundance\tprevalence_tier\tconfidence"
        )?;

        // Data rows
        for r in &self.results {
            writeln!(
                writer,
                "{}\t{}\t{:.6}\t{:.6}\t{:.4}\t{:.2e}\t{:.2e}\t{:.4}\t{:.4}\t{}\t{}",
                r.feature_id,
                r.coefficient,
                r.estimate,
                r.std_error,
                r.statistic,
                r.p_value,
                r.q_value,
                r.prevalence,
                r.mean_abundance,
                r.prevalence_tier.name(),
                r.confidence.name()
            )?;
        }

        Ok(())
    }

    /// Iterate over results.
    pub fn iter(&self) -> impl Iterator<Item = &DaResult> {
        self.results.iter()
    }
}

/// Summary statistics for a result set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResultSummary {
    pub total: usize,
    pub significant_001: usize,
    pub significant_01: usize,
    pub significant_05: usize,
    pub significant_10: usize,
}

impl std::fmt::Display for ResultSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Total features tested: {}", self.total)?;
        writeln!(f, "Significant at q < 0.001: {}", self.significant_001)?;
        writeln!(f, "Significant at q < 0.01:  {}", self.significant_01)?;
        writeln!(f, "Significant at q < 0.05:  {}", self.significant_05)?;
        writeln!(f, "Significant at q < 0.10:  {}", self.significant_10)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prevalence_tier() {
        assert_eq!(PrevalenceTier::from_prevalence(0.02), PrevalenceTier::VeryRare);
        assert_eq!(PrevalenceTier::from_prevalence(0.07), PrevalenceTier::Rare);
        assert_eq!(PrevalenceTier::from_prevalence(0.15), PrevalenceTier::Low);
        assert_eq!(PrevalenceTier::from_prevalence(0.35), PrevalenceTier::Medium);
        assert_eq!(PrevalenceTier::from_prevalence(0.60), PrevalenceTier::High);
        assert_eq!(PrevalenceTier::from_prevalence(0.90), PrevalenceTier::VeryHigh);
    }

    #[test]
    fn test_confidence() {
        assert_eq!(Confidence::from_qvalue(0.0005), Confidence::VeryHigh);
        assert_eq!(Confidence::from_qvalue(0.005), Confidence::High);
        assert_eq!(Confidence::from_qvalue(0.03), Confidence::Moderate);
        assert_eq!(Confidence::from_qvalue(0.08), Confidence::Low);
        assert_eq!(Confidence::from_qvalue(0.15), Confidence::Suggestive);
        assert_eq!(Confidence::from_qvalue(0.30), Confidence::NotSignificant);
    }

    #[test]
    fn test_da_result() {
        let result = DaResult::new(
            "taxon_1".to_string(),
            "groupB".to_string(),
            1.5,
            0.3,
            5.0,
            0.001,
            0.02,
            0.75,
            100.0,
        );

        assert!(result.is_significant());
        assert_eq!(result.confidence, Confidence::Moderate);
        assert_eq!(result.prevalence_tier, PrevalenceTier::VeryHigh); // 0.75 >= 0.75 threshold
    }

    #[test]
    fn test_result_set_summary() {
        let results = vec![
            DaResult::new("t1".into(), "g".into(), 1.0, 0.1, 10.0, 0.0001, 0.0005, 0.5, 100.0),
            DaResult::new("t2".into(), "g".into(), 0.5, 0.1, 5.0, 0.01, 0.02, 0.6, 80.0),
            DaResult::new("t3".into(), "g".into(), 0.2, 0.1, 2.0, 0.1, 0.15, 0.4, 50.0),
            DaResult::new("t4".into(), "g".into(), 0.1, 0.1, 1.0, 0.5, 0.6, 0.3, 30.0),
        ];

        let set = DaResultSet::new("test".to_string(), results);
        let summary = set.summary();

        assert_eq!(summary.total, 4);
        assert_eq!(summary.significant_001, 1);
        assert_eq!(summary.significant_05, 2);
        assert_eq!(summary.significant_10, 2);
    }
}
