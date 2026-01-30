//! Benjamini-Hochberg false discovery rate correction.

use crate::data::{DaResult, DaResultSet};
use crate::profile::PrevalenceProfile;
use crate::test::{PermutationResults, WaldResult};
use serde::{Deserialize, Serialize};

/// Result of BH correction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BhCorrected {
    /// Feature IDs in original order.
    pub feature_ids: Vec<String>,
    /// Original p-values.
    pub p_values: Vec<f64>,
    /// Adjusted p-values (q-values).
    pub q_values: Vec<f64>,
    /// Number of tests.
    pub n_tests: usize,
}

impl BhCorrected {
    /// Get q-value for a specific feature.
    pub fn get_qvalue(&self, feature_id: &str) -> Option<f64> {
        let idx = self.feature_ids.iter().position(|f| f == feature_id)?;
        self.q_values.get(idx).copied()
    }

    /// Count significant results at a threshold.
    pub fn n_significant(&self, alpha: f64) -> usize {
        self.q_values.iter().filter(|&&q| q < alpha).count()
    }

    /// Get indices of significant results.
    pub fn significant_indices(&self, alpha: f64) -> Vec<usize> {
        self.q_values
            .iter()
            .enumerate()
            .filter(|(_, &q)| q < alpha)
            .map(|(i, _)| i)
            .collect()
    }
}

/// Apply Benjamini-Hochberg FDR correction.
///
/// The BH procedure controls the false discovery rate (FDR) at level α.
/// For each p-value, the adjusted p-value (q-value) is calculated as:
/// q[i] = min(p[i] * n / rank[i], q[i+1])
///
/// # Arguments
/// * `p_values` - Raw p-values
/// * `feature_ids` - Feature identifiers (same order as p_values)
///
/// # Returns
/// BhCorrected containing q-values.
pub fn correct_bh(p_values: &[f64], feature_ids: &[String]) -> BhCorrected {
    let n = p_values.len();
    if n == 0 {
        return BhCorrected {
            feature_ids: vec![],
            p_values: vec![],
            q_values: vec![],
            n_tests: 0,
        };
    }

    // Create sorted index
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        p_values[a]
            .partial_cmp(&p_values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Calculate adjusted p-values
    let mut q_sorted = vec![0.0; n];
    let n_f64 = n as f64;

    // Start from largest p-value
    q_sorted[n - 1] = p_values[indices[n - 1]].min(1.0);

    // Work backwards
    for i in (0..n - 1).rev() {
        let rank = i + 1;
        let adjusted = p_values[indices[i]] * n_f64 / rank as f64;
        q_sorted[i] = adjusted.min(q_sorted[i + 1]).min(1.0);
    }

    // Restore original order
    let mut q_values = vec![0.0; n];
    for (i, &orig_idx) in indices.iter().enumerate() {
        q_values[orig_idx] = q_sorted[i];
    }

    BhCorrected {
        feature_ids: feature_ids.to_vec(),
        p_values: p_values.to_vec(),
        q_values,
        n_tests: n,
    }
}

/// Apply BH correction to Wald test results.
pub fn correct_bh_wald(wald: &WaldResult) -> BhCorrected {
    let p_values = wald.p_values();
    let feature_ids: Vec<String> = wald.feature_ids().iter().map(|s| s.to_string()).collect();
    correct_bh(&p_values, &feature_ids)
}

/// Apply BH correction to permutation test results.
pub fn correct_bh_permutation(perm: &PermutationResults) -> BhCorrected {
    let p_values: Vec<f64> = perm.results.iter().map(|r| r.p_value).collect();
    let feature_ids: Vec<String> = perm.results.iter().map(|r| r.feature_id.clone()).collect();
    correct_bh(&p_values, &feature_ids)
}

/// Create full DA results from Wald test and BH correction.
///
/// Combines Wald test statistics with corrected p-values and prevalence
/// information to create a complete result set.
pub fn create_results(
    wald: &WaldResult,
    bh: &BhCorrected,
    prevalence: &PrevalenceProfile,
    mean_abundances: &[f64],
    method: &str,
) -> DaResultSet {
    let results: Vec<DaResult> = wald
        .results
        .iter()
        .enumerate()
        .map(|(i, w)| {
            let q_value = bh.q_values.get(i).copied().unwrap_or(f64::NAN);
            let prev = prevalence.feature_prevalence.get(i).copied().unwrap_or(0.0);
            let mean_abund = mean_abundances.get(i).copied().unwrap_or(0.0);

            DaResult::new(
                w.feature_id.clone(),
                w.coefficient.clone(),
                w.estimate,
                w.std_error,
                w.statistic,
                w.p_value,
                q_value,
                prev,
                mean_abund,
            )
        })
        .collect();

    DaResultSet::new(method.to_string(), results)
}

/// Create full DA results from permutation test and BH correction.
///
/// Combines permutation test statistics with corrected p-values and prevalence
/// information to create a complete result set.
pub fn create_results_permutation(
    perm: &PermutationResults,
    bh: &BhCorrected,
    prevalence: &PrevalenceProfile,
    mean_abundances: &[f64],
    method: &str,
) -> DaResultSet {
    let results: Vec<DaResult> = perm
        .results
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let q_value = bh.q_values.get(i).copied().unwrap_or(f64::NAN);
            let prev = prevalence.feature_prevalence.get(i).copied().unwrap_or(0.0);
            let mean_abund = mean_abundances.get(i).copied().unwrap_or(0.0);

            // Use observed_stat as the test statistic
            DaResult::new(
                p.feature_id.clone(),
                perm.coefficient.clone(),
                p.estimate,
                p.std_error,
                p.observed_stat, // The observed t-statistic
                p.p_value,
                q_value,
                prev,
                mean_abund,
            )
        })
        .collect();

    DaResultSet::new(method.to_string(), results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bh_basic() {
        let p_values = vec![0.01, 0.04, 0.03, 0.005];
        let feature_ids: Vec<String> = (0..4).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        assert_eq!(corrected.n_tests, 4);
        assert_eq!(corrected.p_values, p_values);
    }

    #[test]
    fn test_bh_ordering() {
        // P-values in non-sorted order
        let p_values = vec![0.04, 0.01, 0.03, 0.005];
        let feature_ids: Vec<String> = (0..4).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        // Smallest p-value (0.005 at index 3) should give smallest q-value
        // q = 0.005 * 4 / 1 = 0.02
        assert_relative_eq!(corrected.q_values[3], 0.02, epsilon = 1e-10);

        // Second smallest (0.01 at index 1)
        // q = min(0.01 * 4 / 2, q[next]) = min(0.02, 0.02) = 0.02
        assert_relative_eq!(corrected.q_values[1], 0.02, epsilon = 1e-10);
    }

    #[test]
    fn test_bh_monotonicity() {
        let p_values = vec![0.001, 0.01, 0.02, 0.05, 0.1, 0.5];
        let feature_ids: Vec<String> = (0..6).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        // Q-values should be monotonically non-decreasing when p-values are sorted
        let mut sorted_q: Vec<_> = p_values
            .iter()
            .enumerate()
            .map(|(i, _)| corrected.q_values[i])
            .collect();
        let mut prev = 0.0;
        for q in &sorted_q {
            assert!(*q >= prev || (*q - prev).abs() < 1e-10);
            prev = *q;
        }
    }

    #[test]
    fn test_bh_bounded() {
        let p_values = vec![0.5, 0.6, 0.7, 0.8, 0.9];
        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        // All q-values should be <= 1.0
        for q in &corrected.q_values {
            assert!(*q <= 1.0);
        }
    }

    #[test]
    fn test_bh_empty() {
        let corrected = correct_bh(&[], &[]);
        assert_eq!(corrected.n_tests, 0);
        assert!(corrected.q_values.is_empty());
    }

    #[test]
    fn test_bh_single() {
        let p_values = vec![0.05];
        let feature_ids = vec!["feat_0".to_string()];

        let corrected = correct_bh(&p_values, &feature_ids);

        assert_eq!(corrected.n_tests, 1);
        assert_relative_eq!(corrected.q_values[0], 0.05, epsilon = 1e-10);
    }

    #[test]
    fn test_n_significant() {
        let p_values = vec![0.001, 0.01, 0.03, 0.1, 0.5];
        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        // Count significant at various thresholds
        let n_sig_01 = corrected.n_significant(0.01);
        let n_sig_05 = corrected.n_significant(0.05);
        let n_sig_10 = corrected.n_significant(0.10);

        assert!(n_sig_01 <= n_sig_05);
        assert!(n_sig_05 <= n_sig_10);
    }

    #[test]
    fn test_bh_known_values() {
        // Example from Benjamini & Hochberg (1995) style calculation
        // 5 tests, p = [0.005, 0.01, 0.02, 0.04, 0.1]
        let p_values = vec![0.005, 0.01, 0.02, 0.04, 0.1];
        let feature_ids: Vec<String> = (0..5).map(|i| format!("feat_{}", i)).collect();

        let corrected = correct_bh(&p_values, &feature_ids);

        // Manual calculation:
        // Rank 1: 0.005 * 5/1 = 0.025
        // Rank 2: 0.01 * 5/2 = 0.025, min with next = 0.025
        // Rank 3: 0.02 * 5/3 = 0.0333, min with next = 0.0333
        // Rank 4: 0.04 * 5/4 = 0.05, min with next = 0.05
        // Rank 5: 0.1 * 5/5 = 0.1

        assert_relative_eq!(corrected.q_values[0], 0.025, epsilon = 1e-10);
        assert_relative_eq!(corrected.q_values[1], 0.025, epsilon = 1e-10);
        assert_relative_eq!(corrected.q_values[2], 1.0/30.0, epsilon = 1e-10);
        assert_relative_eq!(corrected.q_values[3], 0.05, epsilon = 1e-10);
        assert_relative_eq!(corrected.q_values[4], 0.1, epsilon = 1e-10);
    }
}
