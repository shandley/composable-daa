//! Additive Log-Ratio (ALR) transformation for compositional data.
//!
//! ALR transforms compositional data by taking the log-ratio of each component
//! relative to a reference component. Unlike CLR, ALR produces an unconstrained
//! coordinate system but requires choosing a reference taxon.
//!
//! # Reference Selection
//!
//! The choice of reference taxon affects interpretation:
//! - Results are relative to the reference (reference taxon is excluded from output)
//! - Choosing a stable/invariant reference improves interpretability
//! - `LeastVariable` automatically selects the most stable taxon
//!
//! # Comparison with CLR
//!
//! | Property | CLR | ALR |
//! |----------|-----|-----|
//! | Output dimensions | n | n-1 |
//! | Sum constraint | Sums to zero | None |
//! | Reference | Geometric mean | Single taxon |
//! | Interpretation | Relative to average | Relative to reference |

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use crate::normalize::clr::TransformedMatrix;
use crate::zero::pseudocount::add_pseudocount;
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Strategy for selecting the reference taxon in ALR transformation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ReferenceSelection {
    /// Use a specific taxon index as reference.
    Index(usize),
    /// Use a specific taxon ID as reference.
    TaxonId(String),
    /// Automatically select the least variable taxon (most stable).
    /// Variability is measured as coefficient of variation (CV).
    LeastVariable,
    /// Automatically select the most abundant taxon (highest mean).
    MostAbundant,
    /// Automatically select the most prevalent taxon (present in most samples).
    MostPrevalent,
}

impl Default for ReferenceSelection {
    fn default() -> Self {
        ReferenceSelection::LeastVariable
    }
}

/// Result of ALR transformation with reference information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlrMatrix {
    /// The transformed data (features × samples).
    /// Note: Has one fewer row than input (reference excluded).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers (excludes reference).
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// Name of the transformation.
    pub transformation: String,
    /// Index of the reference taxon in the original matrix.
    pub reference_index: usize,
    /// ID of the reference taxon.
    pub reference_id: String,
    /// Selection method used.
    pub selection_method: ReferenceSelection,
    /// Reference taxon values (for back-transformation if needed).
    pub reference_values: Vec<f64>,
}

impl AlrMatrix {
    /// Get the transformed value for a feature and sample.
    pub fn get(&self, feature: usize, sample: usize) -> f64 {
        self.data[(feature, sample)]
    }

    /// Number of features (one less than input).
    pub fn n_features(&self) -> usize {
        self.data.nrows()
    }

    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.data.ncols()
    }

    /// Get a row (feature) as a vector.
    pub fn row(&self, feature: usize) -> Vec<f64> {
        self.data.row(feature).iter().cloned().collect()
    }

    /// Get a column (sample) as a vector.
    pub fn col(&self, sample: usize) -> Vec<f64> {
        self.data.column(sample).iter().cloned().collect()
    }

    /// Get reference to the underlying matrix.
    pub fn matrix(&self) -> &DMatrix<f64> {
        &self.data
    }

    /// Convert to a TransformedMatrix for pipeline compatibility.
    ///
    /// Note: This loses some ALR-specific information (reference details).
    pub fn to_transformed_matrix(&self) -> TransformedMatrix {
        TransformedMatrix {
            data: self.data.clone(),
            feature_ids: self.feature_ids.clone(),
            sample_ids: self.sample_ids.clone(),
            transformation: format!("ALR(ref={})", self.reference_id),
            geometric_means: self.reference_values.clone(), // Store reference values here
        }
    }
}

/// Apply Additive Log-Ratio (ALR) transformation.
///
/// ALR transforms compositional data by computing log-ratios relative to
/// a reference taxon. The output has n-1 features (reference excluded).
///
/// # Formula
/// ALR(x_ij) = log(x_ij / x_ref,j)
///
/// Where x_ref,j is the reference taxon value for sample j.
///
/// # Arguments
/// * `data` - Dense matrix with pseudocounts already added (no zeros!)
/// * `feature_ids` - Feature identifiers
/// * `sample_ids` - Sample identifiers
/// * `reference` - Strategy for selecting the reference taxon
///
/// # Returns
/// An AlrMatrix containing ALR-transformed values (n-1 features).
///
/// # Example
/// ```ignore
/// use composable_daa::normalize::alr::{norm_alr, ReferenceSelection};
///
/// // Use least variable taxon as reference
/// let result = norm_alr(&data, feature_ids, sample_ids, ReferenceSelection::LeastVariable)?;
///
/// // Use specific taxon as reference
/// let result = norm_alr(&data, feature_ids, sample_ids,
///     ReferenceSelection::TaxonId("Escherichia".into()))?;
/// ```
pub fn norm_alr(
    data: &DMatrix<f64>,
    feature_ids: Vec<String>,
    sample_ids: Vec<String>,
    reference: ReferenceSelection,
) -> Result<AlrMatrix> {
    let (n_features, n_samples) = data.shape();

    if n_features < 2 {
        return Err(DaaError::InvalidParameter(
            "ALR requires at least 2 features".to_string(),
        ));
    }

    if n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply ALR to empty matrix".to_string(),
        ));
    }

    // Check for zeros or negative values
    for i in 0..n_features {
        for j in 0..n_samples {
            let val = data[(i, j)];
            if val <= 0.0 {
                return Err(DaaError::Numerical(format!(
                    "ALR requires positive values; found {} at ({}, {})",
                    val, i, j
                )));
            }
        }
    }

    // Select reference taxon
    let ref_idx = select_reference(data, &feature_ids, &reference)?;
    let ref_id = feature_ids[ref_idx].clone();

    // Extract reference values
    let reference_values: Vec<f64> = (0..n_samples)
        .map(|j| data[(ref_idx, j)])
        .collect();

    // Calculate log of reference values
    let log_ref: Vec<f64> = reference_values.iter().map(|x| x.ln()).collect();

    // Build output feature list (excluding reference)
    let output_feature_ids: Vec<String> = feature_ids
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != ref_idx)
        .map(|(_, id)| id.clone())
        .collect();

    // Apply ALR transformation
    let n_output_features = n_features - 1;
    let mut alr_data = DMatrix::zeros(n_output_features, n_samples);

    let mut output_row = 0;
    for i in 0..n_features {
        if i == ref_idx {
            continue;
        }
        for j in 0..n_samples {
            // ALR = log(x_i) - log(x_ref)
            alr_data[(output_row, j)] = data[(i, j)].ln() - log_ref[j];
        }
        output_row += 1;
    }

    Ok(AlrMatrix {
        data: alr_data,
        feature_ids: output_feature_ids,
        sample_ids,
        transformation: "ALR".to_string(),
        reference_index: ref_idx,
        reference_id: ref_id,
        selection_method: reference,
        reference_values,
    })
}

/// Select reference taxon based on the specified strategy.
fn select_reference(
    data: &DMatrix<f64>,
    feature_ids: &[String],
    reference: &ReferenceSelection,
) -> Result<usize> {
    let (n_features, n_samples) = data.shape();

    match reference {
        ReferenceSelection::Index(idx) => {
            if *idx >= n_features {
                return Err(DaaError::InvalidParameter(format!(
                    "Reference index {} out of bounds (n_features = {})",
                    idx, n_features
                )));
            }
            Ok(*idx)
        }

        ReferenceSelection::TaxonId(id) => {
            feature_ids
                .iter()
                .position(|x| x == id)
                .ok_or_else(|| DaaError::InvalidParameter(format!(
                    "Reference taxon '{}' not found in feature IDs",
                    id
                )))
        }

        ReferenceSelection::LeastVariable => {
            // Select taxon with lowest coefficient of variation
            let cvs: Vec<(usize, f64)> = (0..n_features)
                .into_par_iter()
                .map(|i| {
                    let values: Vec<f64> = (0..n_samples).map(|j| data[(i, j)]).collect();
                    let mean = values.iter().sum::<f64>() / n_samples as f64;
                    let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                        / n_samples as f64;
                    let cv = if mean > 0.0 {
                        variance.sqrt() / mean
                    } else {
                        f64::INFINITY
                    };
                    (i, cv)
                })
                .collect();

            cvs.into_iter()
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|(idx, _)| idx)
                .ok_or_else(|| DaaError::EmptyData("No features to select from".to_string()))
        }

        ReferenceSelection::MostAbundant => {
            // Select taxon with highest mean abundance
            let means: Vec<(usize, f64)> = (0..n_features)
                .into_par_iter()
                .map(|i| {
                    let mean = (0..n_samples).map(|j| data[(i, j)]).sum::<f64>() / n_samples as f64;
                    (i, mean)
                })
                .collect();

            means
                .into_iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .map(|(idx, _)| idx)
                .ok_or_else(|| DaaError::EmptyData("No features to select from".to_string()))
        }

        ReferenceSelection::MostPrevalent => {
            // Select taxon present in the most samples
            // Since we require positive values (pseudocount added), use a threshold
            let threshold = 1.0; // Values > 1.0 are considered "present" beyond pseudocount

            let prevalences: Vec<(usize, usize)> = (0..n_features)
                .into_par_iter()
                .map(|i| {
                    let count = (0..n_samples).filter(|&j| data[(i, j)] > threshold).count();
                    (i, count)
                })
                .collect();

            prevalences
                .into_iter()
                .max_by_key(|(_, count)| *count)
                .map(|(idx, _)| idx)
                .ok_or_else(|| DaaError::EmptyData("No features to select from".to_string()))
        }
    }
}

/// Convenience function to apply ALR directly to a CountMatrix with pseudocount.
pub fn norm_alr_with_pseudocount(
    counts: &CountMatrix,
    pseudocount: f64,
    reference: ReferenceSelection,
) -> Result<AlrMatrix> {
    let data = add_pseudocount(counts, pseudocount)?;
    norm_alr(
        &data,
        counts.feature_ids().to_vec(),
        counts.sample_ids().to_vec(),
        reference,
    )
}

/// Convenience function using default reference selection (least variable).
pub fn norm_alr_default(
    data: &DMatrix<f64>,
    feature_ids: Vec<String>,
    sample_ids: Vec<String>,
) -> Result<AlrMatrix> {
    norm_alr(data, feature_ids, sample_ids, ReferenceSelection::default())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn create_test_data() -> (DMatrix<f64>, Vec<String>, Vec<String>) {
        // 4 features x 3 samples (already with pseudocount added)
        // Feature C is designed to be the least variable
        let data = DMatrix::from_row_slice(
            4,
            3,
            &[
                10.0, 20.0, 15.0, // feature A - variable
                30.0, 40.0, 35.0, // feature B - variable
                10.0, 10.0, 10.0, // feature C - constant (least variable)
                5.0, 15.0, 10.0,  // feature D - variable
            ],
        );
        let feature_ids = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into()];
        (data, feature_ids, sample_ids)
    }

    #[test]
    fn test_alr_basic() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::Index(2), // Use C as reference
        )
        .unwrap();

        // Output should have 3 features (4 - 1)
        assert_eq!(result.n_features(), 3);
        assert_eq!(result.n_samples(), 3);
        assert_eq!(result.reference_index, 2);
        assert_eq!(result.reference_id, "C");

        // Feature C should not be in output
        assert!(!result.feature_ids.contains(&"C".to_string()));
        assert_eq!(result.feature_ids, vec!["A", "B", "D"]);
    }

    #[test]
    fn test_alr_manual_calculation() {
        // Simple case: 3 features x 2 samples
        let data = DMatrix::from_row_slice(
            3,
            2,
            &[
                2.0, 8.0, // feature A
                4.0, 4.0, // feature B (reference)
                8.0, 2.0, // feature C
            ],
        );
        let feature_ids = vec!["A".into(), "B".into(), "C".into()];
        let sample_ids = vec!["S1".into(), "S2".into()];

        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::Index(1), // B as reference
        )
        .unwrap();

        // ALR[A, S1] = log(2) - log(4) = log(0.5)
        assert_relative_eq!(result.get(0, 0), (0.5_f64).ln(), epsilon = 1e-10);

        // ALR[A, S2] = log(8) - log(4) = log(2)
        assert_relative_eq!(result.get(0, 1), (2.0_f64).ln(), epsilon = 1e-10);

        // ALR[C, S1] = log(8) - log(4) = log(2)
        assert_relative_eq!(result.get(1, 0), (2.0_f64).ln(), epsilon = 1e-10);

        // ALR[C, S2] = log(2) - log(4) = log(0.5)
        assert_relative_eq!(result.get(1, 1), (0.5_f64).ln(), epsilon = 1e-10);
    }

    #[test]
    fn test_alr_by_taxon_id() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::TaxonId("B".into()),
        )
        .unwrap();

        assert_eq!(result.reference_id, "B");
        assert_eq!(result.reference_index, 1);
        assert!(!result.feature_ids.contains(&"B".to_string()));
    }

    #[test]
    fn test_alr_least_variable() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::LeastVariable,
        )
        .unwrap();

        // Feature C (constant at 10.0) should be selected as least variable
        assert_eq!(result.reference_id, "C");
        assert_eq!(result.reference_index, 2);
    }

    #[test]
    fn test_alr_most_abundant() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::MostAbundant,
        )
        .unwrap();

        // Feature B (mean = 35) should be selected as most abundant
        assert_eq!(result.reference_id, "B");
        assert_eq!(result.reference_index, 1);
    }

    #[test]
    fn test_alr_invalid_index() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::Index(10), // Out of bounds
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_alr_invalid_taxon_id() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::TaxonId("NonExistent".into()),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_alr_rejects_zeros() {
        let data = DMatrix::from_row_slice(
            2,
            2,
            &[
                1.0, 0.0, // has zero!
                4.0, 1.0,
            ],
        );
        let result = norm_alr(
            &data,
            vec!["A".into(), "B".into()],
            vec!["S1".into(), "S2".into()],
            ReferenceSelection::Index(0),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_alr_requires_two_features() {
        let data = DMatrix::from_row_slice(1, 2, &[1.0, 2.0]);
        let result = norm_alr(
            &data,
            vec!["A".into()],
            vec!["S1".into(), "S2".into()],
            ReferenceSelection::Index(0),
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_alr_to_transformed_matrix() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let alr = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::Index(2),
        )
        .unwrap();

        let transformed = alr.to_transformed_matrix();
        assert_eq!(transformed.n_features(), 3);
        assert_eq!(transformed.n_samples(), 3);
        assert!(transformed.transformation.contains("ALR"));
        assert!(transformed.transformation.contains("C")); // Reference ID
    }

    #[test]
    fn test_alr_reference_values_preserved() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_alr(
            &data,
            feature_ids,
            sample_ids,
            ReferenceSelection::Index(2), // C has values [10, 10, 10]
        )
        .unwrap();

        assert_eq!(result.reference_values, vec![10.0, 10.0, 10.0]);
    }
}
