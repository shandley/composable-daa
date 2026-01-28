//! Centered Log-Ratio (CLR) transformation for compositional data.

use crate::error::{DaaError, Result};
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// A transformed matrix with metadata about the transformation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransformedMatrix {
    /// The transformed data (features × samples).
    #[serde(skip)]
    pub data: DMatrix<f64>,
    /// Feature identifiers.
    pub feature_ids: Vec<String>,
    /// Sample identifiers.
    pub sample_ids: Vec<String>,
    /// Name of the transformation applied.
    pub transformation: String,
    /// Geometric means per sample (used in CLR).
    pub geometric_means: Vec<f64>,
}

impl TransformedMatrix {
    /// Get the transformed value for a feature and sample.
    pub fn get(&self, feature: usize, sample: usize) -> f64 {
        self.data[(feature, sample)]
    }

    /// Number of features.
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
}

/// Apply Centered Log-Ratio (CLR) transformation.
///
/// CLR transforms compositional data by taking the log of each value
/// divided by the geometric mean of that sample. This removes the
/// simplex constraint and allows standard statistical methods.
///
/// # Formula
/// For sample j: CLR(x_ij) = log(x_ij) - mean(log(x_j))
///
/// Where mean(log(x_j)) is the arithmetic mean of log values for sample j,
/// which equals log(geometric_mean(x_j)).
///
/// # Arguments
/// * `data` - Dense matrix with pseudocounts already added (no zeros!)
/// * `feature_ids` - Feature identifiers
/// * `sample_ids` - Sample identifiers
///
/// # Returns
/// A TransformedMatrix containing CLR-transformed values.
///
/// # Note
/// Input must have no zeros (add pseudocounts first).
pub fn norm_clr(
    data: &DMatrix<f64>,
    feature_ids: Vec<String>,
    sample_ids: Vec<String>,
) -> Result<TransformedMatrix> {
    let (n_features, n_samples) = data.shape();

    if n_features == 0 || n_samples == 0 {
        return Err(DaaError::EmptyData(
            "Cannot apply CLR to empty matrix".to_string(),
        ));
    }

    // Check for zeros or negative values
    for i in 0..n_features {
        for j in 0..n_samples {
            let val = data[(i, j)];
            if val <= 0.0 {
                return Err(DaaError::Numerical(format!(
                    "CLR requires positive values; found {} at ({}, {})",
                    val, i, j
                )));
            }
        }
    }

    // Calculate log-transformed data
    let log_data: DMatrix<f64> = data.map(|x| x.ln());

    // Calculate geometric mean for each sample (column)
    // Geometric mean = exp(mean(log(x)))
    let geometric_means: Vec<f64> = (0..n_samples)
        .into_par_iter()
        .map(|j| {
            let col_sum: f64 = (0..n_features).map(|i| log_data[(i, j)]).sum();
            let mean_log = col_sum / n_features as f64;
            mean_log.exp()
        })
        .collect();

    // Apply CLR: log(x) - mean(log(x)) = log(x / geometric_mean)
    let mut clr_data = DMatrix::zeros(n_features, n_samples);
    for j in 0..n_samples {
        let log_geom_mean = geometric_means[j].ln();
        for i in 0..n_features {
            clr_data[(i, j)] = log_data[(i, j)] - log_geom_mean;
        }
    }

    Ok(TransformedMatrix {
        data: clr_data,
        feature_ids,
        sample_ids,
        transformation: "CLR".to_string(),
        geometric_means,
    })
}

/// Convenience function to apply CLR directly to a CountMatrix with pseudocount.
pub fn norm_clr_with_pseudocount(
    counts: &crate::data::CountMatrix,
    pseudocount: f64,
) -> Result<TransformedMatrix> {
    let data = crate::zero::pseudocount::add_pseudocount(counts, pseudocount)?;
    norm_clr(
        &data,
        counts.feature_ids().to_vec(),
        counts.sample_ids().to_vec(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn create_test_data() -> (DMatrix<f64>, Vec<String>, Vec<String>) {
        // 3 features × 4 samples (already with pseudocount added)
        let data = DMatrix::from_row_slice(3, 4, &[
            10.5, 20.5, 15.5, 5.5,   // feature 0
            30.5, 40.5, 35.5, 25.5,  // feature 1
            5.5,  10.5, 8.5,  3.5,   // feature 2
        ]);
        let feature_ids = vec!["A".into(), "B".into(), "C".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into(), "S4".into()];
        (data, feature_ids, sample_ids)
    }

    #[test]
    fn test_clr_basic() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_clr(&data, feature_ids, sample_ids).unwrap();

        assert_eq!(result.n_features(), 3);
        assert_eq!(result.n_samples(), 4);
        assert_eq!(result.transformation, "CLR");
    }

    #[test]
    fn test_clr_column_sums_zero() {
        // CLR values within each sample should sum to zero
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_clr(&data, feature_ids, sample_ids).unwrap();

        for j in 0..result.n_samples() {
            let col_sum: f64 = (0..result.n_features()).map(|i| result.get(i, j)).sum();
            assert_relative_eq!(col_sum, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_clr_geometric_mean() {
        let (data, feature_ids, sample_ids) = create_test_data();
        let result = norm_clr(&data, feature_ids, sample_ids).unwrap();

        // Verify geometric mean calculation for first sample
        // geom_mean = (10.5 * 30.5 * 5.5)^(1/3)
        let expected_geom_mean = (10.5_f64 * 30.5 * 5.5).powf(1.0 / 3.0);
        assert_relative_eq!(result.geometric_means[0], expected_geom_mean, epsilon = 1e-10);
    }

    #[test]
    fn test_clr_manual_calculation() {
        // Simple case: 2 features × 2 samples
        let data = DMatrix::from_row_slice(2, 2, &[
            1.0, 4.0,  // feature 0
            4.0, 1.0,  // feature 1
        ]);
        let feature_ids = vec!["A".into(), "B".into()];
        let sample_ids = vec!["S1".into(), "S2".into()];
        let result = norm_clr(&data, feature_ids, sample_ids).unwrap();

        // Sample 0: geom_mean = sqrt(1*4) = 2
        // CLR[0,0] = log(1) - log(2) = -log(2)
        // CLR[1,0] = log(4) - log(2) = log(2)
        assert_relative_eq!(result.get(0, 0), -2.0_f64.ln(), epsilon = 1e-10);
        assert_relative_eq!(result.get(1, 0), 2.0_f64.ln(), epsilon = 1e-10);

        // Sample 1: geom_mean = sqrt(4*1) = 2
        // CLR[0,1] = log(4) - log(2) = log(2)
        // CLR[1,1] = log(1) - log(2) = -log(2)
        assert_relative_eq!(result.get(0, 1), 2.0_f64.ln(), epsilon = 1e-10);
        assert_relative_eq!(result.get(1, 1), -2.0_f64.ln(), epsilon = 1e-10);
    }

    #[test]
    fn test_clr_rejects_zeros() {
        let data = DMatrix::from_row_slice(2, 2, &[
            1.0, 0.0,  // has zero!
            4.0, 1.0,
        ]);
        let result = norm_clr(&data, vec!["A".into(), "B".into()], vec!["S1".into(), "S2".into()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_clr_rejects_negative() {
        let data = DMatrix::from_row_slice(2, 2, &[
            1.0, -1.0,  // has negative!
            4.0, 1.0,
        ]);
        let result = norm_clr(&data, vec!["A".into(), "B".into()], vec!["S1".into(), "S2".into()]);
        assert!(result.is_err());
    }
}
