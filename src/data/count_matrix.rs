//! Count matrix with sparse storage for microbiome abundance data.

use crate::error::{DaaError, Result};
use rayon::prelude::*;
use sprs::{CsMat, TriMat};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// A sparse count matrix storing feature abundances across samples.
///
/// Rows represent features (taxa/genes), columns represent samples.
/// Uses CSR (Compressed Sparse Row) format for efficient row-wise operations.
#[derive(Debug, Clone)]
pub struct CountMatrix {
    /// Sparse matrix in CSR format (features × samples)
    data: CsMat<u64>,
    /// Feature identifiers (row names)
    feature_ids: Vec<String>,
    /// Sample identifiers (column names)
    sample_ids: Vec<String>,
}

impl CountMatrix {
    /// Create a new CountMatrix from a sparse matrix and identifiers.
    pub fn new(
        data: CsMat<u64>,
        feature_ids: Vec<String>,
        sample_ids: Vec<String>,
    ) -> Result<Self> {
        let (nrows, ncols) = data.shape();
        if nrows != feature_ids.len() {
            return Err(DaaError::DimensionMismatch {
                expected: nrows,
                actual: feature_ids.len(),
            });
        }
        if ncols != sample_ids.len() {
            return Err(DaaError::DimensionMismatch {
                expected: ncols,
                actual: sample_ids.len(),
            });
        }
        Ok(Self {
            data,
            feature_ids,
            sample_ids,
        })
    }

    /// Load a count matrix from a TSV file.
    ///
    /// Expected format:
    /// - First row: header with sample IDs (first column is feature ID header)
    /// - Subsequent rows: feature ID followed by counts
    pub fn from_tsv<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        // Parse header
        let header_line = lines
            .next()
            .ok_or_else(|| DaaError::EmptyData("Empty TSV file".to_string()))??;
        let header: Vec<&str> = header_line.split('\t').collect();
        if header.len() < 2 {
            return Err(DaaError::EmptyData(
                "TSV must have at least one sample".to_string(),
            ));
        }
        let sample_ids: Vec<String> = header[1..].iter().map(|s| s.to_string()).collect();
        let n_samples = sample_ids.len();

        // Parse data rows into triplets for sparse matrix construction
        let mut triplets: Vec<(usize, usize, u64)> = Vec::new();
        let mut feature_ids: Vec<String> = Vec::new();

        for (row_idx, line_result) in lines.enumerate() {
            let line = line_result?;
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                continue;
            }

            let feature_id = fields[0].to_string();
            feature_ids.push(feature_id);

            for (col_idx, value_str) in fields[1..].iter().enumerate() {
                if col_idx >= n_samples {
                    break;
                }
                let value: u64 = value_str.trim().parse().map_err(|_| DaaError::InvalidCount {
                    value: value_str.to_string(),
                    row: row_idx,
                    col: col_idx,
                })?;
                if value > 0 {
                    triplets.push((row_idx, col_idx, value));
                }
            }
        }

        let n_features = feature_ids.len();
        if n_features == 0 {
            return Err(DaaError::EmptyData("No features in TSV".to_string()));
        }

        // Build sparse matrix from triplets
        let mut tri_mat = TriMat::new((n_features, n_samples));
        for (row, col, val) in triplets {
            tri_mat.add_triplet(row, col, val);
        }
        let data: CsMat<u64> = tri_mat.to_csr();

        Self::new(data, feature_ids, sample_ids)
    }

    /// Write the count matrix to a TSV file.
    pub fn to_tsv<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Write header
        write!(writer, "feature_id")?;
        for sample_id in &self.sample_ids {
            write!(writer, "\t{}", sample_id)?;
        }
        writeln!(writer)?;

        // Write data rows
        for (row_idx, feature_id) in self.feature_ids.iter().enumerate() {
            write!(writer, "{}", feature_id)?;
            for col_idx in 0..self.n_samples() {
                let value = self.get(row_idx, col_idx);
                write!(writer, "\t{}", value)?;
            }
            writeln!(writer)?;
        }

        Ok(())
    }

    /// Get the value at (row, col), returning 0 for missing entries.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> u64 {
        self.data.get(row, col).copied().unwrap_or(0)
    }

    /// Number of features (rows).
    #[inline]
    pub fn n_features(&self) -> usize {
        self.data.rows()
    }

    /// Number of samples (columns).
    #[inline]
    pub fn n_samples(&self) -> usize {
        self.data.cols()
    }

    /// Total number of non-zero entries.
    #[inline]
    pub fn nnz(&self) -> usize {
        self.data.nnz()
    }

    /// Feature identifiers.
    #[inline]
    pub fn feature_ids(&self) -> &[String] {
        &self.feature_ids
    }

    /// Sample identifiers.
    #[inline]
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Get the underlying sparse matrix.
    #[inline]
    pub fn data(&self) -> &CsMat<u64> {
        &self.data
    }

    /// Iterate over rows (features) as sparse vectors.
    pub fn row_iter(&self) -> impl Iterator<Item = sprs::CsVecViewI<'_, u64, usize>> + '_ {
        self.data.outer_iterator()
    }

    /// Get a dense vector for a specific row (feature).
    pub fn row_dense(&self, row: usize) -> Vec<u64> {
        let mut dense = vec![0u64; self.n_samples()];
        if let Some(row_vec) = self.data.outer_view(row) {
            for (col, &val) in row_vec.iter() {
                dense[col] = val;
            }
        }
        dense
    }

    /// Get a dense vector for a specific column (sample).
    pub fn col_dense(&self, col: usize) -> Vec<u64> {
        (0..self.n_features())
            .map(|row| self.get(row, col))
            .collect()
    }

    /// Compute row sums (total counts per feature).
    pub fn row_sums(&self) -> Vec<u64> {
        (0..self.n_features())
            .into_par_iter()
            .map(|row| {
                self.data
                    .outer_view(row)
                    .map(|v| v.iter().map(|(_, &val)| val).sum())
                    .unwrap_or(0)
            })
            .collect()
    }

    /// Compute column sums (library sizes per sample).
    pub fn col_sums(&self) -> Vec<u64> {
        let mut sums = vec![0u64; self.n_samples()];
        for row_vec in self.data.outer_iterator() {
            for (col, &val) in row_vec.iter() {
                sums[col] += val;
            }
        }
        sums
    }

    /// Subset the matrix to include only specified features (by index).
    pub fn subset_features(&self, indices: &[usize]) -> Result<Self> {
        let n_features = indices.len();
        let n_samples = self.n_samples();

        let mut triplets = Vec::new();
        let mut new_feature_ids = Vec::with_capacity(n_features);

        for (new_row, &old_row) in indices.iter().enumerate() {
            if old_row >= self.n_features() {
                return Err(DaaError::InvalidParameter(format!(
                    "Feature index {} out of bounds",
                    old_row
                )));
            }
            new_feature_ids.push(self.feature_ids[old_row].clone());

            if let Some(row_vec) = self.data.outer_view(old_row) {
                for (col, &val) in row_vec.iter() {
                    triplets.push((new_row, col, val));
                }
            }
        }

        let mut tri_mat = TriMat::new((n_features, n_samples));
        for (row, col, val) in triplets {
            tri_mat.add_triplet(row, col, val);
        }

        Self::new(tri_mat.to_csr(), new_feature_ids, self.sample_ids.clone())
    }

    /// Subset the matrix to include only specified samples (by index).
    pub fn subset_samples(&self, indices: &[usize]) -> Result<Self> {
        let n_features = self.n_features();
        let n_samples = indices.len();

        // Build column index mapping
        let col_map: HashMap<usize, usize> = indices
            .iter()
            .enumerate()
            .map(|(new_idx, &old_idx)| (old_idx, new_idx))
            .collect();

        let mut triplets = Vec::new();
        let mut new_sample_ids = Vec::with_capacity(n_samples);

        for &old_col in indices {
            if old_col >= self.n_samples() {
                return Err(DaaError::InvalidParameter(format!(
                    "Sample index {} out of bounds",
                    old_col
                )));
            }
            new_sample_ids.push(self.sample_ids[old_col].clone());
        }

        for (row, row_vec) in self.data.outer_iterator().enumerate() {
            for (old_col, &val) in row_vec.iter() {
                if let Some(&new_col) = col_map.get(&old_col) {
                    triplets.push((row, new_col, val));
                }
            }
        }

        let mut tri_mat = TriMat::new((n_features, n_samples));
        for (row, col, val) in triplets {
            tri_mat.add_triplet(row, col, val);
        }

        Self::new(tri_mat.to_csr(), self.feature_ids.clone(), new_sample_ids)
    }

    /// Convert to a dense matrix (f64).
    pub fn to_dense(&self) -> nalgebra::DMatrix<f64> {
        let mut dense = nalgebra::DMatrix::zeros(self.n_features(), self.n_samples());
        for (row, row_vec) in self.data.outer_iterator().enumerate() {
            for (col, &val) in row_vec.iter() {
                dense[(row, col)] = val as f64;
            }
        }
        dense
    }

    /// Create from a dense matrix.
    pub fn from_dense(
        data: &nalgebra::DMatrix<f64>,
        feature_ids: Vec<String>,
        sample_ids: Vec<String>,
    ) -> Result<Self> {
        let (nrows, ncols) = data.shape();
        let mut tri_mat = TriMat::new((nrows, ncols));

        for row in 0..nrows {
            for col in 0..ncols {
                let val = data[(row, col)];
                if val > 0.0 {
                    tri_mat.add_triplet(row, col, val.round() as u64);
                }
            }
        }

        Self::new(tri_mat.to_csr(), feature_ids, sample_ids)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_matrix() -> CountMatrix {
        // 3 features × 4 samples
        let mut tri_mat = TriMat::new((3, 4));
        tri_mat.add_triplet(0, 0, 10);
        tri_mat.add_triplet(0, 1, 20);
        tri_mat.add_triplet(0, 2, 0); // explicit zero, won't be stored
        tri_mat.add_triplet(0, 3, 5);
        tri_mat.add_triplet(1, 0, 100);
        tri_mat.add_triplet(1, 1, 200);
        tri_mat.add_triplet(1, 2, 150);
        tri_mat.add_triplet(1, 3, 175);
        tri_mat.add_triplet(2, 0, 1);
        // feature 2 is sparse - only present in sample 0

        let feature_ids = vec!["feat_A".to_string(), "feat_B".to_string(), "feat_C".to_string()];
        let sample_ids = vec![
            "sample1".to_string(),
            "sample2".to_string(),
            "sample3".to_string(),
            "sample4".to_string(),
        ];

        CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids).unwrap()
    }

    #[test]
    fn test_dimensions() {
        let mat = create_test_matrix();
        assert_eq!(mat.n_features(), 3);
        assert_eq!(mat.n_samples(), 4);
    }

    #[test]
    fn test_get_values() {
        let mat = create_test_matrix();
        assert_eq!(mat.get(0, 0), 10);
        assert_eq!(mat.get(0, 2), 0); // zero entry
        assert_eq!(mat.get(2, 0), 1);
        assert_eq!(mat.get(2, 1), 0); // sparse entry
    }

    #[test]
    fn test_row_dense() {
        let mat = create_test_matrix();
        let row0 = mat.row_dense(0);
        assert_eq!(row0, vec![10, 20, 0, 5]);
    }

    #[test]
    fn test_col_sums() {
        let mat = create_test_matrix();
        let sums = mat.col_sums();
        assert_eq!(sums, vec![111, 220, 150, 180]);
    }

    #[test]
    fn test_row_sums() {
        let mat = create_test_matrix();
        let sums = mat.row_sums();
        assert_eq!(sums, vec![35, 625, 1]);
    }

    #[test]
    fn test_tsv_roundtrip() {
        let mat = create_test_matrix();

        let mut temp_file = NamedTempFile::new().unwrap();
        mat.to_tsv(temp_file.path()).unwrap();

        let loaded = CountMatrix::from_tsv(temp_file.path()).unwrap();
        assert_eq!(loaded.n_features(), mat.n_features());
        assert_eq!(loaded.n_samples(), mat.n_samples());
        assert_eq!(loaded.feature_ids(), mat.feature_ids());
        assert_eq!(loaded.sample_ids(), mat.sample_ids());

        for row in 0..mat.n_features() {
            for col in 0..mat.n_samples() {
                assert_eq!(loaded.get(row, col), mat.get(row, col));
            }
        }
    }

    #[test]
    fn test_subset_features() {
        let mat = create_test_matrix();
        let subset = mat.subset_features(&[0, 2]).unwrap();

        assert_eq!(subset.n_features(), 2);
        assert_eq!(subset.n_samples(), 4);
        assert_eq!(subset.feature_ids(), &["feat_A", "feat_C"]);
        assert_eq!(subset.get(0, 0), 10);
        assert_eq!(subset.get(1, 0), 1);
    }

    #[test]
    fn test_subset_samples() {
        let mat = create_test_matrix();
        let subset = mat.subset_samples(&[1, 3]).unwrap();

        assert_eq!(subset.n_features(), 3);
        assert_eq!(subset.n_samples(), 2);
        assert_eq!(subset.sample_ids(), &["sample2", "sample4"]);
        assert_eq!(subset.get(0, 0), 20);
        assert_eq!(subset.get(0, 1), 5);
    }
}
