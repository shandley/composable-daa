//! Error types for the composable-daa library.

use thiserror::Error;

/// Main error type for the library.
#[derive(Error, Debug)]
pub enum DaaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("CSV parsing error: {0}")]
    Csv(#[from] csv::Error),

    #[error("Invalid count value '{value}' at row {row}, column {col}")]
    InvalidCount {
        value: String,
        row: usize,
        col: usize,
    },

    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch { expected: usize, actual: usize },

    #[error("Sample ID mismatch: {0}")]
    SampleMismatch(String),

    #[error("Missing column '{0}' in metadata")]
    MissingColumn(String),

    #[error("Invalid variable type for column '{column}': {reason}")]
    InvalidVariableType { column: String, reason: String },

    #[error("Formula parse error: {0}")]
    FormulaParse(String),

    #[error("Empty data: {0}")]
    EmptyData(String),

    #[error("Numerical error: {0}")]
    Numerical(String),

    #[error("Invalid parameter: {0}")]
    InvalidParameter(String),

    #[error("Pipeline error: {0}")]
    Pipeline(String),

    #[error("YAML serialization error: {0}")]
    Yaml(#[from] serde_yaml::Error),

    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
}

/// Result type alias for library operations.
pub type Result<T> = std::result::Result<T, DaaError>;
