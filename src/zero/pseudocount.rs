//! Pseudocount addition for zero handling.

use crate::data::CountMatrix;
use crate::error::{DaaError, Result};
use nalgebra::DMatrix;

/// Add a pseudocount to all entries, converting to dense matrix.
///
/// This is a simple but effective approach for handling zeros before
/// log transformations like CLR. The pseudocount is added to all entries
/// (both zero and non-zero) to maintain relative relationships.
///
/// # Arguments
/// * `counts` - The count matrix
/// * `pseudocount` - Value to add (typically 0.5 or 1.0)
///
/// # Returns
/// A dense matrix with pseudocounts added.
///
/// # Common choices
/// - 0.5: Often used as a compromise
/// - 1.0: Laplace smoothing, simple and interpretable
/// - Minimum non-zero / 2: Data-adaptive approach
pub fn add_pseudocount(counts: &CountMatrix, pseudocount: f64) -> Result<DMatrix<f64>> {
    if pseudocount <= 0.0 {
        return Err(DaaError::InvalidParameter(
            "Pseudocount must be positive".to_string(),
        ));
    }

    let dense = counts.to_dense();
    let result = dense.add_scalar(pseudocount);
    Ok(result)
}

/// Calculate a data-adaptive pseudocount based on minimum non-zero value.
///
/// Returns min(non-zero values) / 2, which is a commonly used heuristic.
pub fn adaptive_pseudocount(counts: &CountMatrix) -> f64 {
    let mut min_nonzero = f64::INFINITY;

    for row_vec in counts.data().outer_iterator() {
        for (_, &val) in row_vec.iter() {
            if val > 0 {
                let fval = val as f64;
                if fval < min_nonzero {
                    min_nonzero = fval;
                }
            }
        }
    }

    if min_nonzero.is_infinite() {
        0.5 // fallback if all zeros (shouldn't happen in practice)
    } else {
        min_nonzero / 2.0
    }
}

/// Add an adaptive pseudocount (min non-zero / 2).
pub fn add_pseudocount_adaptive(counts: &CountMatrix) -> Result<DMatrix<f64>> {
    let pc = adaptive_pseudocount(counts);
    add_pseudocount(counts, pc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use sprs::TriMat;

    fn create_test_matrix() -> CountMatrix {
        let mut tri_mat = TriMat::new((2, 3));
        tri_mat.add_triplet(0, 0, 10);
        tri_mat.add_triplet(0, 1, 20);
        // (0, 2) is zero
        tri_mat.add_triplet(1, 0, 5);
        tri_mat.add_triplet(1, 1, 0); // explicit zero
        tri_mat.add_triplet(1, 2, 15);

        let feature_ids = vec!["A".into(), "B".into()];
        let sample_ids = vec!["S1".into(), "S2".into(), "S3".into()];
        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_add_pseudocount() {
        let counts = create_test_matrix();
        let result = add_pseudocount(&counts, 0.5).unwrap();

        assert_eq!(result.nrows(), 2);
        assert_eq!(result.ncols(), 3);
        assert!((result[(0, 0)] - 10.5).abs() < 1e-10);
        assert!((result[(0, 2)] - 0.5).abs() < 1e-10); // was zero
        assert!((result[(1, 1)] - 0.5).abs() < 1e-10); // was zero
    }

    #[test]
    fn test_adaptive_pseudocount() {
        let counts = create_test_matrix();
        let pc = adaptive_pseudocount(&counts);

        // Minimum non-zero is 5, so adaptive should be 2.5
        assert!((pc - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_pseudocount() {
        let counts = create_test_matrix();
        assert!(add_pseudocount(&counts, 0.0).is_err());
        assert!(add_pseudocount(&counts, -1.0).is_err());
    }
}
