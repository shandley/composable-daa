//! Design matrix construction from metadata and formula.

use crate::data::{Formula, Metadata, Term, Variable, VariableType};
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use std::collections::HashMap;

/// A design matrix for linear modeling.
#[derive(Debug, Clone)]
pub struct DesignMatrix {
    /// The design matrix (samples × coefficients).
    matrix: DMatrix<f64>,
    /// Names of the coefficients (columns).
    coefficient_names: Vec<String>,
    /// Sample IDs (rows).
    sample_ids: Vec<String>,
    /// Reference levels for categorical variables.
    reference_levels: HashMap<String, String>,
}

impl DesignMatrix {
    /// Create a design matrix directly from components.
    ///
    /// This is useful for creating reduced models for LRT or other custom designs.
    pub fn from_matrix(
        matrix: DMatrix<f64>,
        coefficient_names: Vec<String>,
        sample_ids: Vec<String>,
    ) -> Self {
        Self {
            matrix,
            coefficient_names,
            sample_ids,
            reference_levels: HashMap::new(),
        }
    }

    /// Build a design matrix from metadata and formula.
    pub fn from_formula(metadata: &Metadata, formula: &Formula) -> Result<Self> {
        let sample_ids = metadata.sample_ids().to_vec();
        let n_samples = sample_ids.len();

        // Validate that all formula variables exist in metadata
        for var in formula.variables() {
            if !metadata.has_column(var) {
                return Err(DaaError::MissingColumn(var.to_string()));
            }
        }

        // Determine reference levels for categorical variables (alphabetically first)
        let mut reference_levels = HashMap::new();
        for var in formula.variables() {
            if metadata.column_type(var) == Some(VariableType::Categorical) {
                let levels = metadata.levels(var)?;
                if !levels.is_empty() {
                    reference_levels.insert(var.to_string(), levels[0].clone());
                }
            }
        }

        // Build coefficient names and matrix columns
        let mut coefficient_names = Vec::new();
        let mut columns: Vec<Vec<f64>> = Vec::new();

        // Intercept
        if formula.intercept {
            coefficient_names.push("(Intercept)".to_string());
            columns.push(vec![1.0; n_samples]);
        }

        // Process each term
        for term in &formula.terms {
            match term {
                Term::Intercept => {
                    // Already handled above
                }
                Term::Main(var_name) => {
                    let var_type = metadata.column_type(var_name);
                    let values = metadata.column(var_name)?;

                    match var_type {
                        Some(VariableType::Continuous) | Some(VariableType::Ordinal) => {
                            // Single column for numeric variable
                            coefficient_names.push(var_name.clone());
                            let col: Vec<f64> = values
                                .iter()
                                .map(|v| match v {
                                    Variable::Continuous(x) => *x,
                                    Variable::Ordinal(x) => *x as f64,
                                    _ => 0.0, // Handle missing as 0 (or could error)
                                })
                                .collect();
                            columns.push(col);
                        }
                        Some(VariableType::Categorical) | None => {
                            // Dummy coding for categorical variable
                            let levels = metadata.levels(var_name)?;
                            let ref_level = reference_levels.get(var_name);

                            for level in &levels {
                                // Skip reference level only when there's an intercept
                                if formula.intercept && Some(level) == ref_level {
                                    continue;
                                }
                                coefficient_names.push(format!("{}{}", var_name, level));
                                let col: Vec<f64> = values
                                    .iter()
                                    .map(|v| {
                                        if let Variable::Categorical(s) = v {
                                            if s == level { 1.0 } else { 0.0 }
                                        } else {
                                            0.0
                                        }
                                    })
                                    .collect();
                                columns.push(col);
                            }
                        }
                    }
                }
                Term::Interaction(var1, var2) => {
                    // Build interaction columns
                    let cols1 = Self::get_term_columns(metadata, var1, &reference_levels)?;
                    let cols2 = Self::get_term_columns(metadata, var2, &reference_levels)?;

                    for (name1, col1) in &cols1 {
                        for (name2, col2) in &cols2 {
                            coefficient_names.push(format!("{}:{}", name1, name2));
                            let col: Vec<f64> = col1
                                .iter()
                                .zip(col2.iter())
                                .map(|(a, b)| a * b)
                                .collect();
                            columns.push(col);
                        }
                    }
                }
            }
        }

        // Build matrix
        let n_coef = columns.len();
        let mut matrix = DMatrix::zeros(n_samples, n_coef);
        for (col_idx, col) in columns.iter().enumerate() {
            for (row_idx, &val) in col.iter().enumerate() {
                matrix[(row_idx, col_idx)] = val;
            }
        }

        Ok(Self {
            matrix,
            coefficient_names,
            sample_ids,
            reference_levels,
        })
    }

    /// Get columns for a single variable (helper for interactions).
    fn get_term_columns(
        metadata: &Metadata,
        var_name: &str,
        reference_levels: &HashMap<String, String>,
    ) -> Result<Vec<(String, Vec<f64>)>> {
        let var_type = metadata.column_type(var_name);
        let values = metadata.column(var_name)?;
        let mut result = Vec::new();

        match var_type {
            Some(VariableType::Continuous) | Some(VariableType::Ordinal) => {
                let col: Vec<f64> = values
                    .iter()
                    .map(|v| match v {
                        Variable::Continuous(x) => *x,
                        Variable::Ordinal(x) => *x as f64,
                        _ => 0.0,
                    })
                    .collect();
                result.push((var_name.to_string(), col));
            }
            Some(VariableType::Categorical) | None => {
                let levels = metadata.levels(var_name)?;
                let ref_level = reference_levels.get(var_name);

                for level in &levels {
                    if Some(level) == ref_level {
                        continue;
                    }
                    let col: Vec<f64> = values
                        .iter()
                        .map(|v| {
                            if let Variable::Categorical(s) = v {
                                if s == level { 1.0 } else { 0.0 }
                            } else {
                                0.0
                            }
                        })
                        .collect();
                    result.push((format!("{}{}", var_name, level), col));
                }
            }
        }

        Ok(result)
    }

    /// Get the design matrix.
    pub fn matrix(&self) -> &DMatrix<f64> {
        &self.matrix
    }

    /// Get coefficient names.
    pub fn coefficient_names(&self) -> &[String] {
        &self.coefficient_names
    }

    /// Get sample IDs.
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Number of samples (rows).
    pub fn n_samples(&self) -> usize {
        self.matrix.nrows()
    }

    /// Number of coefficients (columns).
    pub fn n_coefficients(&self) -> usize {
        self.matrix.ncols()
    }

    /// Get the reference level for a categorical variable.
    pub fn reference_level(&self, variable: &str) -> Option<&str> {
        self.reference_levels.get(variable).map(|s| s.as_str())
    }

    /// Get the index of a coefficient by name.
    pub fn coefficient_index(&self, name: &str) -> Option<usize> {
        self.coefficient_names.iter().position(|n| n == name)
    }

    /// Check if the matrix has an intercept.
    pub fn has_intercept(&self) -> bool {
        self.coefficient_names.first() == Some(&"(Intercept)".to_string())
    }

    /// Set a custom reference level for a categorical variable.
    /// This rebuilds the design matrix with the new reference.
    pub fn with_reference_level(
        mut self,
        variable: &str,
        level: &str,
        metadata: &Metadata,
        formula: &Formula,
    ) -> Result<Self> {
        // Validate the level exists
        let levels = metadata.levels(variable)?;
        if !levels.contains(&level.to_string()) {
            return Err(DaaError::InvalidParameter(format!(
                "Level '{}' not found for variable '{}'",
                level, variable
            )));
        }

        self.reference_levels.insert(variable.to_string(), level.to_string());

        // Rebuild with new reference
        DesignMatrix::from_formula(metadata, formula)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_metadata() -> Metadata {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        writeln!(file, "S1\tcontrol\t25").unwrap();
        writeln!(file, "S2\ttreatment\t30").unwrap();
        writeln!(file, "S3\tcontrol\t35").unwrap();
        writeln!(file, "S4\ttreatment\t28").unwrap();
        file.flush().unwrap();
        Metadata::from_tsv(file.path()).unwrap()
    }

    #[test]
    fn test_intercept_only() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ 1").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert_eq!(dm.n_samples(), 4);
        assert_eq!(dm.n_coefficients(), 1);
        assert_eq!(dm.coefficient_names(), &["(Intercept)"]);
        assert!(dm.matrix().iter().all(|&v| v == 1.0));
    }

    #[test]
    fn test_continuous_variable() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ age").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert_eq!(dm.n_coefficients(), 2); // intercept + age
        assert_eq!(dm.coefficient_names(), &["(Intercept)", "age"]);

        // Check age values
        let age_col: Vec<f64> = (0..4).map(|i| dm.matrix()[(i, 1)]).collect();
        assert_eq!(age_col, vec![25.0, 30.0, 35.0, 28.0]);
    }

    #[test]
    fn test_categorical_variable() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ group").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert_eq!(dm.n_coefficients(), 2); // intercept + grouptreatment
        assert_eq!(dm.coefficient_names(), &["(Intercept)", "grouptreatment"]);
        assert_eq!(dm.reference_level("group"), Some("control"));

        // Check dummy coding: S1=control(0), S2=treatment(1), S3=control(0), S4=treatment(1)
        let group_col: Vec<f64> = (0..4).map(|i| dm.matrix()[(i, 1)]).collect();
        assert_eq!(group_col, vec![0.0, 1.0, 0.0, 1.0]);
    }

    #[test]
    fn test_multiple_variables() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ group + age").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert_eq!(dm.n_coefficients(), 3);
        assert_eq!(dm.coefficient_names(), &["(Intercept)", "grouptreatment", "age"]);
    }

    #[test]
    fn test_no_intercept() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ 0 + group").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert!(!dm.has_intercept());
        // Without intercept, both levels should be included
        assert_eq!(dm.coefficient_names(), &["groupcontrol", "grouptreatment"]);
    }

    #[test]
    fn test_interaction() {
        let meta = create_test_metadata();
        let formula = Formula::parse("~ group * age").unwrap();
        let dm = DesignMatrix::from_formula(&meta, &formula).unwrap();

        assert_eq!(dm.n_coefficients(), 4);
        assert_eq!(
            dm.coefficient_names(),
            &["(Intercept)", "grouptreatment", "age", "grouptreatment:age"]
        );

        // Check interaction: grouptreatment * age
        // S1: 0 * 25 = 0, S2: 1 * 30 = 30, S3: 0 * 35 = 0, S4: 1 * 28 = 28
        let interaction_col: Vec<f64> = (0..4).map(|i| dm.matrix()[(i, 3)]).collect();
        assert_eq!(interaction_col, vec![0.0, 30.0, 0.0, 28.0]);
    }
}
