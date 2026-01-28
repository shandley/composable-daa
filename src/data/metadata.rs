//! Sample metadata handling for differential abundance analysis.

use crate::error::{DaaError, Result};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A variable value that can be categorical, continuous, or ordinal.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Variable {
    /// Categorical variable with string levels.
    Categorical(String),
    /// Continuous numeric variable.
    Continuous(f64),
    /// Ordinal variable with integer rank.
    Ordinal(i64),
    /// Missing value.
    Missing,
}

impl Variable {
    /// Check if this is a missing value.
    pub fn is_missing(&self) -> bool {
        matches!(self, Variable::Missing)
    }

    /// Try to get as categorical string.
    pub fn as_categorical(&self) -> Option<&str> {
        match self {
            Variable::Categorical(s) => Some(s),
            _ => None,
        }
    }

    /// Try to get as continuous f64.
    pub fn as_continuous(&self) -> Option<f64> {
        match self {
            Variable::Continuous(v) => Some(*v),
            _ => None,
        }
    }

    /// Try to get as ordinal i64.
    pub fn as_ordinal(&self) -> Option<i64> {
        match self {
            Variable::Ordinal(v) => Some(*v),
            _ => None,
        }
    }
}

/// Type hint for columns when loading metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VariableType {
    Categorical,
    Continuous,
    Ordinal,
}

/// Sample metadata containing variables for each sample.
#[derive(Debug, Clone)]
pub struct Metadata {
    /// Sample IDs in order.
    sample_ids: Vec<String>,
    /// Column names.
    column_names: Vec<String>,
    /// Data stored as sample_id -> column_name -> Variable.
    data: HashMap<String, HashMap<String, Variable>>,
    /// Type hints for each column.
    column_types: HashMap<String, VariableType>,
}

impl Metadata {
    /// Create empty metadata.
    pub fn new() -> Self {
        Self {
            sample_ids: Vec::new(),
            column_names: Vec::new(),
            data: HashMap::new(),
            column_types: HashMap::new(),
        }
    }

    /// Load metadata from a TSV file.
    ///
    /// Expected format:
    /// - First row: header with column names (first column is sample ID)
    /// - Subsequent rows: sample ID followed by variable values
    ///
    /// By default, columns are inferred as continuous if all values parse as numbers,
    /// otherwise categorical. Use `with_column_types` to override.
    pub fn from_tsv<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Parse header
        let header_line = lines
            .next()
            .ok_or_else(|| DaaError::EmptyData("Empty metadata file".to_string()))??;
        let header: Vec<&str> = header_line.split('\t').collect();
        if header.len() < 2 {
            return Err(DaaError::EmptyData(
                "Metadata must have at least one variable column".to_string(),
            ));
        }
        let column_names: Vec<String> = header[1..].iter().map(|s| s.to_string()).collect();

        // First pass: collect all values to infer types
        let mut raw_data: Vec<(String, Vec<String>)> = Vec::new();
        for line_result in lines {
            let line = line_result?;
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }
            let sample_id = fields[0].to_string();
            let values: Vec<String> = fields[1..]
                .iter()
                .map(|s| s.to_string())
                .collect();
            raw_data.push((sample_id, values));
        }

        if raw_data.is_empty() {
            return Err(DaaError::EmptyData("No samples in metadata".to_string()));
        }

        // Infer column types
        let mut column_types = HashMap::new();
        for (col_idx, col_name) in column_names.iter().enumerate() {
            let all_numeric = raw_data.iter().all(|(_, values)| {
                if col_idx >= values.len() {
                    return true; // missing, skip
                }
                let v = values[col_idx].trim();
                v.is_empty() || v == "NA" || v == "na" || v.parse::<f64>().is_ok()
            });
            let var_type = if all_numeric {
                VariableType::Continuous
            } else {
                VariableType::Categorical
            };
            column_types.insert(col_name.clone(), var_type);
        }

        // Build metadata
        let mut sample_ids = Vec::new();
        let mut data = HashMap::new();

        for (sample_id, values) in raw_data {
            sample_ids.push(sample_id.clone());
            let mut sample_data = HashMap::new();

            for (col_idx, col_name) in column_names.iter().enumerate() {
                let var = if col_idx >= values.len() {
                    Variable::Missing
                } else {
                    let raw = values[col_idx].trim();
                    if raw.is_empty() || raw == "NA" || raw == "na" {
                        Variable::Missing
                    } else {
                        match column_types.get(col_name) {
                            Some(VariableType::Continuous) => {
                                match raw.parse::<f64>() {
                                    Ok(v) => Variable::Continuous(v),
                                    Err(_) => Variable::Missing,
                                }
                            }
                            Some(VariableType::Ordinal) => {
                                match raw.parse::<i64>() {
                                    Ok(v) => Variable::Ordinal(v),
                                    Err(_) => Variable::Missing,
                                }
                            }
                            Some(VariableType::Categorical) | None => {
                                Variable::Categorical(raw.to_string())
                            }
                        }
                    }
                };
                sample_data.insert(col_name.clone(), var);
            }
            data.insert(sample_id, sample_data);
        }

        Ok(Self {
            sample_ids,
            column_names,
            data,
            column_types,
        })
    }

    /// Set type hints for specific columns.
    pub fn with_column_types(mut self, types: HashMap<String, VariableType>) -> Self {
        // Re-parse values according to new types
        for (col_name, var_type) in &types {
            self.column_types.insert(col_name.clone(), *var_type);

            // Re-interpret values
            for sample_data in self.data.values_mut() {
                if let Some(var) = sample_data.get_mut(col_name) {
                    *var = match var {
                        Variable::Categorical(s) => {
                            let trimmed = s.trim();
                            match var_type {
                                VariableType::Continuous => {
                                    trimmed.parse::<f64>()
                                        .map(Variable::Continuous)
                                        .unwrap_or(Variable::Missing)
                                }
                                VariableType::Ordinal => {
                                    trimmed.parse::<i64>()
                                        .map(Variable::Ordinal)
                                        .unwrap_or(Variable::Missing)
                                }
                                VariableType::Categorical => Variable::Categorical(s.clone()),
                            }
                        }
                        Variable::Continuous(v) => match var_type {
                            VariableType::Continuous => Variable::Continuous(*v),
                            VariableType::Ordinal => Variable::Ordinal(*v as i64),
                            VariableType::Categorical => Variable::Categorical(v.to_string()),
                        },
                        Variable::Ordinal(v) => match var_type {
                            VariableType::Continuous => Variable::Continuous(*v as f64),
                            VariableType::Ordinal => Variable::Ordinal(*v),
                            VariableType::Categorical => Variable::Categorical(v.to_string()),
                        },
                        Variable::Missing => Variable::Missing,
                    };
                }
            }
        }
        self
    }

    /// Sample IDs in order.
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Column names.
    pub fn column_names(&self) -> &[String] {
        &self.column_names
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.sample_ids.len()
    }

    /// Number of columns (variables).
    pub fn n_columns(&self) -> usize {
        self.column_names.len()
    }

    /// Get a variable value for a specific sample and column.
    pub fn get(&self, sample_id: &str, column: &str) -> Option<&Variable> {
        self.data.get(sample_id).and_then(|m| m.get(column))
    }

    /// Get all values for a column.
    pub fn column(&self, column: &str) -> Result<Vec<&Variable>> {
        if !self.column_names.contains(&column.to_string()) {
            return Err(DaaError::MissingColumn(column.to_string()));
        }
        Ok(self
            .sample_ids
            .iter()
            .map(|sid| {
                self.data
                    .get(sid)
                    .and_then(|m| m.get(column))
                    .unwrap_or(&Variable::Missing)
            })
            .collect())
    }

    /// Get the type of a column.
    pub fn column_type(&self, column: &str) -> Option<VariableType> {
        self.column_types.get(column).copied()
    }

    /// Get unique levels for a categorical column.
    pub fn levels(&self, column: &str) -> Result<Vec<String>> {
        let values = self.column(column)?;
        let mut levels: Vec<String> = values
            .iter()
            .filter_map(|v| v.as_categorical().map(String::from))
            .collect::<HashSet<_>>()
            .into_iter()
            .collect();
        levels.sort();
        Ok(levels)
    }

    /// Subset metadata to only include specified samples.
    pub fn subset_samples(&self, sample_ids: &[String]) -> Result<Self> {
        let mut new_data = HashMap::new();
        let mut new_sample_ids = Vec::new();

        for sid in sample_ids {
            if let Some(sample_data) = self.data.get(sid) {
                new_data.insert(sid.clone(), sample_data.clone());
                new_sample_ids.push(sid.clone());
            } else {
                return Err(DaaError::SampleMismatch(format!(
                    "Sample '{}' not found in metadata",
                    sid
                )));
            }
        }

        Ok(Self {
            sample_ids: new_sample_ids,
            column_names: self.column_names.clone(),
            data: new_data,
            column_types: self.column_types.clone(),
        })
    }

    /// Align metadata to match the sample order in a count matrix.
    pub fn align_to(&self, sample_ids: &[String]) -> Result<Self> {
        self.subset_samples(sample_ids)
    }

    /// Check if a sample exists.
    pub fn has_sample(&self, sample_id: &str) -> bool {
        self.data.contains_key(sample_id)
    }

    /// Check if a column exists.
    pub fn has_column(&self, column: &str) -> bool {
        self.column_names.contains(&column.to_string())
    }
}

impl Default for Metadata {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_tsv() -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage\tseverity").unwrap();
        writeln!(file, "S1\tcontrol\t25\t1").unwrap();
        writeln!(file, "S2\ttreatment\t30\t2").unwrap();
        writeln!(file, "S3\tcontrol\t35\t1").unwrap();
        writeln!(file, "S4\ttreatment\t28\t3").unwrap();
        file.flush().unwrap();
        file
    }

    #[test]
    fn test_load_metadata() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        assert_eq!(meta.n_samples(), 4);
        assert_eq!(meta.n_columns(), 3);
        assert_eq!(meta.sample_ids(), &["S1", "S2", "S3", "S4"]);
        assert_eq!(meta.column_names(), &["group", "age", "severity"]);
    }

    #[test]
    fn test_get_value() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        let val = meta.get("S1", "group").unwrap();
        assert_eq!(val.as_categorical(), Some("control"));

        let val = meta.get("S2", "age").unwrap();
        assert_eq!(val.as_continuous(), Some(30.0));
    }

    #[test]
    fn test_column_type_inference() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        // "group" should be categorical (strings)
        assert_eq!(meta.column_type("group"), Some(VariableType::Categorical));
        // "age" should be continuous (numbers)
        assert_eq!(meta.column_type("age"), Some(VariableType::Continuous));
    }

    #[test]
    fn test_with_column_types() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        // Force severity to be ordinal
        let mut types = HashMap::new();
        types.insert("severity".to_string(), VariableType::Ordinal);
        let meta = meta.with_column_types(types);

        let val = meta.get("S2", "severity").unwrap();
        assert_eq!(val.as_ordinal(), Some(2));
    }

    #[test]
    fn test_levels() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        let levels = meta.levels("group").unwrap();
        assert_eq!(levels, vec!["control", "treatment"]);
    }

    #[test]
    fn test_subset_samples() {
        let file = create_test_tsv();
        let meta = Metadata::from_tsv(file.path()).unwrap();

        let subset = meta.subset_samples(&["S1".to_string(), "S3".to_string()]).unwrap();
        assert_eq!(subset.n_samples(), 2);
        assert_eq!(subset.sample_ids(), &["S1", "S3"]);
    }

    #[test]
    fn test_missing_values() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, "sample_id\tgroup\tage").unwrap();
        writeln!(file, "S1\tcontrol\t25").unwrap();
        writeln!(file, "S2\ttreatment\tNA").unwrap();
        writeln!(file, "S3\t\t30").unwrap();
        file.flush().unwrap();

        let meta = Metadata::from_tsv(file.path()).unwrap();

        assert!(meta.get("S2", "age").unwrap().is_missing());
        assert!(meta.get("S3", "group").unwrap().is_missing());
    }
}
