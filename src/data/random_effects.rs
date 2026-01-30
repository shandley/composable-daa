//! Random effects specification for mixed models.
//!
//! Supports lme4-style syntax for specifying random effects:
//! - `(1 | subject)` - random intercept per subject
//! - `(1 + time | subject)` - random intercept and slope
//! - `(0 + time | subject)` - random slope only (no intercept)

use crate::data::{Formula, Metadata, Variable};
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A single random effect term.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RandomEffect {
    /// Terms in the random effect (e.g., ["1"] for intercept, ["1", "time"] for intercept + slope).
    pub terms: Vec<String>,
    /// Grouping variable (e.g., "subject_id").
    pub grouping: String,
    /// Whether to include a random intercept.
    pub has_intercept: bool,
}

impl RandomEffect {
    /// Parse a random effect specification string.
    ///
    /// Supports lme4-style syntax:
    /// - `(1 | subject)` -> intercept only
    /// - `(1 + time | subject)` -> intercept + slope
    /// - `(0 + time | subject)` -> slope only
    /// - `(time | subject)` -> intercept + slope (implicit intercept)
    pub fn parse(spec: &str) -> Result<Self> {
        let spec = spec.trim();

        // Must be wrapped in parentheses
        if !spec.starts_with('(') || !spec.ends_with(')') {
            return Err(DaaError::FormulaParse(format!(
                "Random effect must be wrapped in parentheses: {}",
                spec
            )));
        }

        let inner = &spec[1..spec.len() - 1].trim();

        // Split on |
        let parts: Vec<&str> = inner.split('|').collect();
        if parts.len() != 2 {
            return Err(DaaError::FormulaParse(format!(
                "Random effect must have exactly one '|': {}",
                spec
            )));
        }

        let terms_str = parts[0].trim();
        let grouping = parts[1].trim().to_string();

        if grouping.is_empty() {
            return Err(DaaError::FormulaParse(
                "Random effect grouping variable cannot be empty".to_string(),
            ));
        }

        // Parse terms
        let mut terms = Vec::new();
        let mut has_intercept = true;

        // Check for explicit no-intercept
        let terms_str = if terms_str.starts_with("0 +") || terms_str.starts_with("0+") {
            has_intercept = false;
            terms_str.trim_start_matches("0").trim_start_matches('+').trim()
        } else if terms_str.starts_with("-1 +") || terms_str.starts_with("-1+") {
            has_intercept = false;
            terms_str.trim_start_matches("-1").trim_start_matches('+').trim()
        } else {
            terms_str
        };

        // Split remaining terms by +
        for term in terms_str.split('+') {
            let term = term.trim();
            if term.is_empty() {
                continue;
            }
            if term == "1" {
                // Explicit intercept - already handled by has_intercept
                continue;
            }
            if term == "0" || term == "-1" {
                has_intercept = false;
                continue;
            }
            terms.push(term.to_string());
        }

        // If has_intercept is true, add "1" to terms
        if has_intercept {
            terms.insert(0, "1".to_string());
        }

        if terms.is_empty() {
            return Err(DaaError::FormulaParse(
                "Random effect must have at least one term".to_string(),
            ));
        }

        Ok(Self {
            terms,
            grouping,
            has_intercept,
        })
    }

    /// Check if this is a random intercept only.
    pub fn is_intercept_only(&self) -> bool {
        self.terms.len() == 1 && self.terms[0] == "1"
    }

    /// Number of random effect terms per group.
    pub fn n_terms(&self) -> usize {
        self.terms.len()
    }
}

impl std::fmt::Display for RandomEffect {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({} | {})", self.terms.join(" + "), self.grouping)
    }
}

/// A formula with both fixed and random effects.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MixedFormula {
    /// Fixed effects formula.
    pub fixed: Formula,
    /// Random effects.
    pub random: Vec<RandomEffect>,
    /// Original formula string.
    pub formula_str: String,
}

impl MixedFormula {
    /// Parse a mixed model formula.
    ///
    /// Supports lme4-style syntax:
    /// - `~ group + time + (1 | subject)`
    /// - `~ group + (1 + time | subject)`
    /// - `~ group + (1 | subject) + (1 | batch)`
    ///
    /// # Examples
    /// ```
    /// use composable_daa::data::MixedFormula;
    /// let f = MixedFormula::parse("~ group + (1 | subject)").unwrap();
    /// assert!(f.fixed.intercept);
    /// assert_eq!(f.random.len(), 1);
    /// assert!(f.random[0].is_intercept_only());
    /// ```
    pub fn parse(formula: &str) -> Result<Self> {
        let formula_str = formula.to_string();
        let formula = formula.trim();

        // Must start with ~
        if !formula.starts_with('~') {
            return Err(DaaError::FormulaParse(
                "Formula must start with '~'".to_string(),
            ));
        }

        let rhs = formula[1..].trim();
        if rhs.is_empty() {
            return Err(DaaError::FormulaParse(
                "Formula right-hand side is empty".to_string(),
            ));
        }

        // Extract random effects using regex
        let re = Regex::new(r"\([^)]+\|[^)]+\)").unwrap();
        let random_strs: Vec<&str> = re.find_iter(rhs).map(|m: regex::Match| m.as_str()).collect();

        // Parse random effects
        let random: Vec<RandomEffect> = random_strs
            .iter()
            .map(|s| RandomEffect::parse(s))
            .collect::<Result<Vec<_>>>()?;

        // Remove random effects from formula to get fixed effects
        let mut fixed_str = rhs.to_string();
        for re_str in &random_strs {
            fixed_str = fixed_str.replace(re_str, "");
        }

        // Clean up: remove double +, trailing +, etc.
        let fixed_str = fixed_str
            .split('+')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect::<Vec<_>>()
            .join(" + ");

        // If no fixed effects remain, default to intercept only
        let fixed_formula_str = if fixed_str.is_empty() {
            "~ 1".to_string()
        } else {
            format!("~ {}", fixed_str)
        };

        let fixed = Formula::parse(&fixed_formula_str)?;

        Ok(Self {
            fixed,
            random,
            formula_str: formula_str.to_string(),
        })
    }

    /// Check if the formula has any random effects.
    pub fn has_random_effects(&self) -> bool {
        !self.random.is_empty()
    }

    /// Get all grouping variables.
    pub fn grouping_variables(&self) -> Vec<&str> {
        self.random.iter().map(|r| r.grouping.as_str()).collect()
    }

    /// Get all variables used (fixed + random grouping + random slope terms).
    pub fn all_variables(&self) -> Vec<&str> {
        let mut vars: Vec<&str> = self.fixed.variables();

        for re in &self.random {
            vars.push(&re.grouping);
            for term in &re.terms {
                if term != "1" && term != "0" {
                    vars.push(term);
                }
            }
        }

        vars.sort();
        vars.dedup();
        vars
    }

    /// Total number of random effects (sum of terms across all random effect specifications).
    pub fn total_random_effects(&self) -> usize {
        self.random.iter().map(|r| r.n_terms()).sum()
    }
}

impl std::fmt::Display for MixedFormula {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.formula_str)
    }
}

/// Design matrix for random effects (Z matrix).
#[derive(Debug, Clone)]
pub struct RandomDesignMatrix {
    /// The Z matrix (samples × total random effects).
    /// For random intercepts: n_samples × n_groups
    /// For random intercepts + slopes: n_samples × (n_groups * n_terms_per_group)
    pub matrix: DMatrix<f64>,
    /// Group indices for each sample (which group each sample belongs to).
    pub group_indices: Vec<usize>,
    /// Unique group IDs.
    pub group_ids: Vec<String>,
    /// Number of groups.
    pub n_groups: usize,
    /// Number of random effect terms per group.
    pub n_random_per_group: usize,
    /// Sample IDs.
    pub sample_ids: Vec<String>,
}

impl RandomDesignMatrix {
    /// Build the random effects design matrix from metadata and a random effect.
    ///
    /// Currently supports:
    /// - Random intercepts: `(1 | group)`
    ///
    /// Future: Random slopes
    pub fn from_random_effect(
        metadata: &Metadata,
        random_effect: &RandomEffect,
    ) -> Result<Self> {
        let sample_ids = metadata.sample_ids().to_vec();
        let n_samples = sample_ids.len();

        // Get group assignments
        let group_column = metadata.column(&random_effect.grouping)?;

        // Build mapping from group value to index
        let mut group_map: HashMap<String, usize> = HashMap::new();
        let mut group_ids: Vec<String> = Vec::new();
        let mut group_indices: Vec<usize> = Vec::with_capacity(n_samples);

        for val in &group_column {
            let group_str = match val {
                Variable::Categorical(s) => s.clone(),
                Variable::Continuous(x) => x.to_string(),
                Variable::Ordinal(x) => x.to_string(),
                Variable::Missing => {
                    return Err(DaaError::InvalidParameter(
                        "Missing values not allowed in grouping variable".to_string(),
                    ));
                }
            };

            let idx = if let Some(&idx) = group_map.get(&group_str) {
                idx
            } else {
                let idx = group_ids.len();
                group_map.insert(group_str.clone(), idx);
                group_ids.push(group_str);
                idx
            };
            group_indices.push(idx);
        }

        let n_groups = group_ids.len();

        // For now, only support random intercepts
        if !random_effect.is_intercept_only() {
            return Err(DaaError::NotImplemented(
                "Random slopes not yet implemented. Use (1 | group) for random intercepts only.".to_string(),
            ));
        }

        let n_random_per_group = 1; // Only intercepts for now
        let n_random_total = n_groups * n_random_per_group;

        // Build Z matrix
        // For random intercepts: Z[i, j] = 1 if sample i belongs to group j, else 0
        let mut z = DMatrix::zeros(n_samples, n_random_total);

        for (sample_idx, &group_idx) in group_indices.iter().enumerate() {
            z[(sample_idx, group_idx)] = 1.0;
        }

        Ok(Self {
            matrix: z,
            group_indices,
            group_ids,
            n_groups,
            n_random_per_group,
            sample_ids,
        })
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.matrix.nrows()
    }

    /// Total number of random effects (columns in Z).
    pub fn n_random_effects(&self) -> usize {
        self.matrix.ncols()
    }

    /// Get the Z matrix.
    pub fn matrix(&self) -> &DMatrix<f64> {
        &self.matrix
    }

    /// Get number of observations per group.
    pub fn observations_per_group(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.n_groups];
        for &idx in &self.group_indices {
            counts[idx] += 1;
        }
        counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\ttime\tsubject").unwrap();
        writeln!(file, "S1\tcontrol\t0\tA").unwrap();
        writeln!(file, "S2\tcontrol\t1\tA").unwrap();
        writeln!(file, "S3\ttreatment\t0\tB").unwrap();
        writeln!(file, "S4\ttreatment\t1\tB").unwrap();
        writeln!(file, "S5\tcontrol\t0\tC").unwrap();
        writeln!(file, "S6\tcontrol\t1\tC").unwrap();
        writeln!(file, "S7\ttreatment\t0\tD").unwrap();
        writeln!(file, "S8\ttreatment\t1\tD").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_parse_random_intercept() {
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        assert_eq!(re.grouping, "subject");
        assert!(re.has_intercept);
        assert!(re.is_intercept_only());
        assert_eq!(re.terms, vec!["1"]);
    }

    #[test]
    fn test_parse_random_slope() {
        let re = RandomEffect::parse("(1 + time | subject)").unwrap();
        assert_eq!(re.grouping, "subject");
        assert!(re.has_intercept);
        assert!(!re.is_intercept_only());
        assert_eq!(re.terms, vec!["1", "time"]);
    }

    #[test]
    fn test_parse_random_slope_no_intercept() {
        let re = RandomEffect::parse("(0 + time | subject)").unwrap();
        assert_eq!(re.grouping, "subject");
        assert!(!re.has_intercept);
        assert_eq!(re.terms, vec!["time"]);
    }

    #[test]
    fn test_parse_implicit_intercept() {
        let re = RandomEffect::parse("(time | subject)").unwrap();
        assert_eq!(re.grouping, "subject");
        assert!(re.has_intercept);
        assert_eq!(re.terms, vec!["1", "time"]);
    }

    #[test]
    fn test_random_effect_display() {
        let re = RandomEffect::parse("(1 + time | subject)").unwrap();
        assert_eq!(re.to_string(), "(1 + time | subject)");
    }

    #[test]
    fn test_mixed_formula_basic() {
        let f = MixedFormula::parse("~ group + (1 | subject)").unwrap();
        assert!(f.fixed.intercept);
        assert_eq!(f.fixed.terms.len(), 1);
        assert_eq!(f.random.len(), 1);
        assert!(f.random[0].is_intercept_only());
        assert_eq!(f.random[0].grouping, "subject");
    }

    #[test]
    fn test_mixed_formula_multiple_fixed() {
        let f = MixedFormula::parse("~ group + time + (1 | subject)").unwrap();
        assert_eq!(f.fixed.terms.len(), 2);
        assert_eq!(f.random.len(), 1);
    }

    #[test]
    fn test_mixed_formula_multiple_random() {
        let f = MixedFormula::parse("~ group + (1 | subject) + (1 | batch)").unwrap();
        assert_eq!(f.fixed.terms.len(), 1);
        assert_eq!(f.random.len(), 2);
        assert_eq!(f.random[0].grouping, "subject");
        assert_eq!(f.random[1].grouping, "batch");
    }

    #[test]
    fn test_mixed_formula_random_slope() {
        let f = MixedFormula::parse("~ group + (1 + time | subject)").unwrap();
        assert_eq!(f.random.len(), 1);
        assert!(!f.random[0].is_intercept_only());
        assert_eq!(f.random[0].terms, vec!["1", "time"]);
    }

    #[test]
    fn test_mixed_formula_grouping_variables() {
        let f = MixedFormula::parse("~ group + (1 | subject) + (1 | batch)").unwrap();
        let vars = f.grouping_variables();
        assert_eq!(vars, vec!["subject", "batch"]);
    }

    #[test]
    fn test_mixed_formula_all_variables() {
        let f = MixedFormula::parse("~ group + time + (1 + time | subject)").unwrap();
        let vars = f.all_variables();
        assert!(vars.contains(&"group"));
        assert!(vars.contains(&"time"));
        assert!(vars.contains(&"subject"));
    }

    #[test]
    fn test_invalid_random_effect() {
        assert!(RandomEffect::parse("1 | subject").is_err()); // Missing parentheses
        assert!(RandomEffect::parse("(1 subject)").is_err()); // Missing |
        assert!(RandomEffect::parse("(1 | )").is_err()); // Empty grouping
    }

    #[test]
    fn test_random_design_matrix_basic() {
        let metadata = create_test_metadata();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let z = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();

        assert_eq!(z.n_samples(), 8);
        assert_eq!(z.n_groups, 4); // A, B, C, D
        assert_eq!(z.n_random_per_group, 1);
        assert_eq!(z.n_random_effects(), 4);

        // Check structure: each row should have exactly one 1
        for i in 0..z.n_samples() {
            let row_sum: f64 = (0..z.n_random_effects()).map(|j| z.matrix[(i, j)]).sum();
            assert_eq!(row_sum, 1.0);
        }
    }

    #[test]
    fn test_random_design_matrix_structure() {
        let metadata = create_test_metadata();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let z = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();

        // Samples S1, S2 belong to subject A (should have same group index)
        assert_eq!(z.group_indices[0], z.group_indices[1]);
        // Samples S3, S4 belong to subject B
        assert_eq!(z.group_indices[2], z.group_indices[3]);
        // Different subjects have different indices
        assert_ne!(z.group_indices[0], z.group_indices[2]);
    }

    #[test]
    fn test_observations_per_group() {
        let metadata = create_test_metadata();
        let re = RandomEffect::parse("(1 | subject)").unwrap();
        let z = RandomDesignMatrix::from_random_effect(&metadata, &re).unwrap();

        let obs = z.observations_per_group();
        // Each subject has 2 observations
        assert!(obs.iter().all(|&n| n == 2));
    }
}
