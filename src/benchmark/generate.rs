//! Synthetic data generation for benchmarking and testing.
//!
//! Generates count matrices with known ground truth for validating
//! differential abundance methods.

use crate::data::{CountMatrix, Metadata};
use crate::error::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for synthetic data generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyntheticConfig {
    /// Name/identifier for this dataset.
    pub name: String,
    /// Number of features (taxa/genes).
    pub n_features: usize,
    /// Number of samples per group.
    pub n_samples_per_group: usize,
    /// Target sparsity (proportion of zeros, 0.0-1.0).
    pub sparsity: f64,
    /// Dominance level (proportion of reads from top 3 features, 0.0-1.0).
    pub dominance: f64,
    /// Coefficient of variation for library sizes.
    pub library_size_cv: f64,
    /// Mean library size per sample.
    pub mean_library_size: u64,
    /// Library size ratio between groups (>1 means group2 has more).
    pub group_library_ratio: f64,
    /// Number of truly differential features.
    pub n_differential: usize,
    /// Log2 fold change for differential features.
    pub effect_size: f64,
    /// Proportion of differential features that are group-specific (presence/absence).
    pub group_specific_fraction: f64,
    /// Random seed for reproducibility.
    pub seed: u64,
}

impl Default for SyntheticConfig {
    fn default() -> Self {
        Self {
            name: "synthetic".to_string(),
            n_features: 200,
            n_samples_per_group: 20,
            sparsity: 0.6,
            dominance: 0.4,
            library_size_cv: 0.3,
            mean_library_size: 50000,
            group_library_ratio: 1.0,
            n_differential: 20,
            effect_size: 1.0, // 2-fold change
            group_specific_fraction: 0.0,
            seed: 42,
        }
    }
}

impl SyntheticConfig {
    /// Create a new config with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// Set sparsity level.
    pub fn with_sparsity(mut self, sparsity: f64) -> Self {
        self.sparsity = sparsity.clamp(0.0, 0.99);
        self
    }

    /// Set dominance level.
    pub fn with_dominance(mut self, dominance: f64) -> Self {
        self.dominance = dominance.clamp(0.0, 0.99);
        self
    }

    /// Set library size parameters.
    pub fn with_library_size(mut self, mean: u64, cv: f64) -> Self {
        self.mean_library_size = mean;
        self.library_size_cv = cv;
        self
    }

    /// Set group library size imbalance.
    pub fn with_group_imbalance(mut self, ratio: f64) -> Self {
        self.group_library_ratio = ratio;
        self
    }

    /// Set differential features.
    pub fn with_differential(mut self, n: usize, effect_size: f64) -> Self {
        self.n_differential = n;
        self.effect_size = effect_size;
        self
    }

    /// Set group-specific (presence/absence) fraction.
    pub fn with_group_specific(mut self, fraction: f64) -> Self {
        self.group_specific_fraction = fraction.clamp(0.0, 1.0);
        self
    }

    /// Set random seed.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Set dimensions.
    pub fn with_dimensions(mut self, n_features: usize, n_samples_per_group: usize) -> Self {
        self.n_features = n_features;
        self.n_samples_per_group = n_samples_per_group;
        self
    }

    // Preset configurations

    /// Ideal case: low sparsity, balanced, easy to analyze.
    pub fn ideal() -> Self {
        Self::new("ideal")
            .with_sparsity(0.2)
            .with_dominance(0.2)
            .with_library_size(50000, 0.2)
            .with_differential(20, 1.0)
    }

    /// Typical 16S microbiome data.
    pub fn typical_16s() -> Self {
        Self::new("typical_16s")
            .with_sparsity(0.6)
            .with_dominance(0.4)
            .with_library_size(30000, 0.4)
            .with_differential(15, 1.0)
    }

    /// Sparse virome-like data.
    pub fn sparse_virome() -> Self {
        Self::new("sparse_virome")
            .with_sparsity(0.85)
            .with_dominance(0.6)
            .with_library_size(20000, 0.6)
            .with_differential(10, 1.5)
    }

    /// Extreme sparsity edge case.
    pub fn extreme_sparse() -> Self {
        Self::new("extreme_sparse")
            .with_sparsity(0.95)
            .with_dominance(0.8)
            .with_library_size(10000, 0.8)
            .with_differential(5, 2.0)
    }

    /// Dataset with group-specific features.
    pub fn group_specific() -> Self {
        Self::new("group_specific")
            .with_sparsity(0.5)
            .with_dominance(0.3)
            .with_group_specific(0.5)
            .with_differential(20, 1.0)
    }

    /// Confounded by library size imbalance.
    pub fn confounded() -> Self {
        Self::new("confounded")
            .with_sparsity(0.5)
            .with_dominance(0.4)
            .with_group_imbalance(3.0)
            .with_differential(15, 1.0)
    }

    /// Small sample size.
    pub fn small_n() -> Self {
        Self::new("small_n")
            .with_dimensions(200, 5)
            .with_sparsity(0.5)
            .with_dominance(0.4)
            .with_differential(10, 1.5)
    }
}

/// Ground truth for a synthetic dataset.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroundTruth {
    /// Feature IDs that are truly differential.
    pub differential_features: Vec<String>,
    /// True log2 fold changes (positive = higher in group2).
    pub true_log2fc: HashMap<String, f64>,
    /// Effect type for each differential feature.
    pub effect_types: HashMap<String, EffectType>,
    /// Direction of effect.
    pub directions: HashMap<String, Direction>,
    /// Null features (not differential).
    pub null_features: Vec<String>,
}

/// Type of differential effect.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EffectType {
    /// Abundance difference (present in both groups).
    Abundance,
    /// Presence/absence difference (group-specific).
    Presence,
}

/// Direction of effect.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Direction {
    /// Higher in group 2 (e.g., treatment/disease).
    Up,
    /// Lower in group 2.
    Down,
}

impl GroundTruth {
    /// Check if a feature is truly differential.
    pub fn is_differential(&self, feature_id: &str) -> bool {
        self.differential_features.contains(&feature_id.to_string())
    }

    /// Get the true log2 fold change for a feature (0.0 for null features).
    pub fn get_log2fc(&self, feature_id: &str) -> f64 {
        self.true_log2fc.get(feature_id).copied().unwrap_or(0.0)
    }

    /// Write ground truth to TSV format.
    pub fn to_tsv(&self) -> String {
        let mut lines = vec!["feature_id\ttrue_log2fc\tdirection\teffect_type\tis_differential".to_string()];

        // Differential features
        for feat in &self.differential_features {
            let log2fc = self.true_log2fc.get(feat).copied().unwrap_or(0.0);
            let dir = self.directions.get(feat).map(|d| match d {
                Direction::Up => "up",
                Direction::Down => "down",
            }).unwrap_or("none");
            let effect = self.effect_types.get(feat).map(|e| match e {
                EffectType::Abundance => "abundance",
                EffectType::Presence => "presence",
            }).unwrap_or("unknown");
            lines.push(format!("{}\t{:.4}\t{}\t{}\ttrue", feat, log2fc, dir, effect));
        }

        // Null features
        for feat in &self.null_features {
            lines.push(format!("{}\t0.0\tnone\tnull\tfalse", feat));
        }

        lines.join("\n")
    }
}

/// Result of synthetic data generation.
#[derive(Debug, Clone)]
pub struct SyntheticData {
    /// Generated count matrix.
    pub counts: CountMatrix,
    /// Sample metadata.
    pub metadata: Metadata,
    /// Ground truth.
    pub ground_truth: GroundTruth,
    /// Configuration used.
    pub config: SyntheticConfig,
}

impl SyntheticData {
    /// Write all files to a directory.
    pub fn write_to_dir(&self, dir: &std::path::Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;

        self.counts.to_tsv(&dir.join("counts.tsv"))?;

        // Write metadata manually
        let metadata_content = self.metadata_to_tsv();
        std::fs::write(dir.join("metadata.tsv"), metadata_content)?;

        let gt_content = self.ground_truth.to_tsv();
        std::fs::write(dir.join("ground_truth.tsv"), gt_content)?;

        let config_yaml = serde_yaml::to_string(&self.config)
            .map_err(|e| crate::error::DaaError::Pipeline(e.to_string()))?;
        std::fs::write(dir.join("config.yaml"), config_yaml)?;

        Ok(())
    }

    /// Convert metadata to TSV string.
    fn metadata_to_tsv(&self) -> String {
        let mut lines = vec!["sample_id\tgroup\tlibrary_size".to_string()];
        for sample_id in self.metadata.sample_ids() {
            let group = self.metadata.get(sample_id, "group")
                .and_then(|v| v.as_categorical())
                .unwrap_or("unknown");
            let lib_size = self.metadata.get(sample_id, "library_size")
                .and_then(|v| v.as_continuous())
                .unwrap_or(0.0) as u64;
            lines.push(format!("{}\t{}\t{}", sample_id, group, lib_size));
        }
        lines.join("\n")
    }
}

/// Simple deterministic RNG (xorshift64).
struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self { state: seed.max(1) }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() as f64) / (u64::MAX as f64)
    }

    /// Generate a value from approximate normal distribution (Box-Muller).
    fn next_normal(&mut self, mean: f64, std: f64) -> f64 {
        let u1 = self.next_f64().max(1e-10);
        let u2 = self.next_f64();
        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
        mean + std * z
    }

    /// Generate from negative binomial (using Poisson-Gamma mixture).
    fn next_negbin(&mut self, mean: f64, dispersion: f64) -> u64 {
        if mean <= 0.0 {
            return 0;
        }
        // Gamma-Poisson mixture
        // Shape = 1/dispersion, Scale = mean * dispersion
        let shape = 1.0 / dispersion;
        let scale = mean * dispersion;

        // Generate gamma using Marsaglia and Tsang's method (simplified)
        let gamma_sample = self.next_gamma(shape, scale);

        // Generate Poisson with gamma rate
        self.next_poisson(gamma_sample)
    }

    /// Simplified gamma sampling.
    fn next_gamma(&mut self, shape: f64, scale: f64) -> f64 {
        if shape < 1.0 {
            // Use Ahrens-Dieter method for shape < 1
            let u = self.next_f64();
            return self.next_gamma(shape + 1.0, scale) * u.powf(1.0 / shape);
        }

        // Marsaglia and Tsang's method for shape >= 1
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();

        loop {
            let x = self.next_normal(0.0, 1.0);
            let v = (1.0 + c * x).powi(3);

            if v > 0.0 {
                let u = self.next_f64();
                if u < 1.0 - 0.0331 * x.powi(4) {
                    return d * v * scale;
                }
                if u.ln() < 0.5 * x.powi(2) + d * (1.0 - v + v.ln()) {
                    return d * v * scale;
                }
            }
        }
    }

    /// Generate Poisson random variable.
    fn next_poisson(&mut self, lambda: f64) -> u64 {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda < 30.0 {
            // Direct method for small lambda
            let l = (-lambda).exp();
            let mut k = 0u64;
            let mut p = 1.0;
            loop {
                k += 1;
                p *= self.next_f64();
                if p <= l {
                    return k - 1;
                }
            }
        } else {
            // Normal approximation for large lambda
            let sample = self.next_normal(lambda, lambda.sqrt());
            sample.max(0.0).round() as u64
        }
    }

    /// Shuffle a vector in place.
    fn shuffle<T>(&mut self, vec: &mut [T]) {
        for i in (1..vec.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            vec.swap(i, j);
        }
    }
}

/// Generate synthetic data with known ground truth.
pub fn generate_synthetic(config: &SyntheticConfig) -> Result<SyntheticData> {
    let mut rng = Rng::new(config.seed);

    let n_features = config.n_features;
    let n_samples = config.n_samples_per_group * 2;
    let n_diff = config.n_differential.min(n_features);
    let n_group_specific = ((n_diff as f64) * config.group_specific_fraction).round() as usize;
    let n_abundance_diff = n_diff - n_group_specific;

    // Generate feature IDs
    let feature_ids: Vec<String> = (0..n_features)
        .map(|i| format!("Feature_{:04}", i))
        .collect();

    // Generate sample IDs
    let sample_ids: Vec<String> = (0..n_samples)
        .map(|i| format!("S{:03}", i))
        .collect();

    // Assign groups
    let groups: Vec<String> = (0..n_samples)
        .map(|i| {
            if i < config.n_samples_per_group {
                "control".to_string()
            } else {
                "treatment".to_string()
            }
        })
        .collect();

    // Generate base abundances with dominance structure
    let base_abundances = generate_abundance_profile(n_features, config.dominance, &mut rng);

    // Select differential features
    let mut feature_indices: Vec<usize> = (0..n_features).collect();
    rng.shuffle(&mut feature_indices);

    let diff_indices: Vec<usize> = feature_indices[..n_diff].to_vec();
    let abundance_diff_indices: Vec<usize> = diff_indices[..n_abundance_diff].to_vec();
    let presence_diff_indices: Vec<usize> = diff_indices[n_abundance_diff..].to_vec();

    // Assign directions (half up, half down)
    let mut directions: HashMap<usize, Direction> = HashMap::new();
    for (i, &idx) in diff_indices.iter().enumerate() {
        directions.insert(idx, if i % 2 == 0 { Direction::Up } else { Direction::Down });
    }

    // Generate library sizes per sample
    let library_sizes: Vec<f64> = (0..n_samples)
        .map(|i| {
            let base = config.mean_library_size as f64;
            let group_factor = if i < config.n_samples_per_group {
                1.0 / config.group_library_ratio.sqrt()
            } else {
                config.group_library_ratio.sqrt()
            };
            let sample_factor = 1.0 + rng.next_normal(0.0, config.library_size_cv);
            (base * group_factor * sample_factor.max(0.1)).max(1000.0)
        })
        .collect();

    // Generate count matrix
    let dispersion = 0.5; // Typical overdispersion
    let mut counts: Vec<Vec<u64>> = vec![vec![0; n_samples]; n_features];

    for feat_idx in 0..n_features {
        let base_prop = base_abundances[feat_idx];

        // Use sparsity directly as probability of zero per entry
        let zero_prob = config.sparsity;

        for sample_idx in 0..n_samples {
            let is_treatment = sample_idx >= config.n_samples_per_group;

            // Handle group-specific (presence/absence) features
            if presence_diff_indices.contains(&feat_idx) {
                let dir = directions.get(&feat_idx).copied().unwrap_or(Direction::Up);
                let present = match (dir, is_treatment) {
                    (Direction::Up, true) => true,   // Present in treatment
                    (Direction::Up, false) => false, // Absent in control
                    (Direction::Down, true) => false, // Absent in treatment
                    (Direction::Down, false) => true, // Present in control
                };
                if !present {
                    counts[feat_idx][sample_idx] = 0;
                    continue;
                }
            }

            // Apply zero probability
            if rng.next_f64() < zero_prob {
                counts[feat_idx][sample_idx] = 0;
                continue;
            }

            // Calculate expected count
            let lib_size = library_sizes[sample_idx];
            let mut expected = base_prop * lib_size;

            // Apply fold change for abundance-differential features
            if abundance_diff_indices.contains(&feat_idx) {
                let dir = directions.get(&feat_idx).copied().unwrap_or(Direction::Up);
                let fc = 2.0_f64.powf(config.effect_size);
                expected *= match (dir, is_treatment) {
                    (Direction::Up, true) => fc,
                    (Direction::Down, true) => 1.0 / fc,
                    _ => 1.0,
                };
            }

            // Generate count from negative binomial
            counts[feat_idx][sample_idx] = rng.next_negbin(expected, dispersion);
        }
    }

    // Build ground truth
    let differential_features: Vec<String> = diff_indices
        .iter()
        .map(|&i| feature_ids[i].clone())
        .collect();

    let null_features: Vec<String> = (0..n_features)
        .filter(|i| !diff_indices.contains(i))
        .map(|i| feature_ids[i].clone())
        .collect();

    let mut true_log2fc: HashMap<String, f64> = HashMap::new();
    let mut effect_types: HashMap<String, EffectType> = HashMap::new();
    let mut dir_map: HashMap<String, Direction> = HashMap::new();

    for &idx in &abundance_diff_indices {
        let feat_id = feature_ids[idx].clone();
        let dir = directions.get(&idx).copied().unwrap_or(Direction::Up);
        let log2fc = match dir {
            Direction::Up => config.effect_size,
            Direction::Down => -config.effect_size,
        };
        true_log2fc.insert(feat_id.clone(), log2fc);
        effect_types.insert(feat_id.clone(), EffectType::Abundance);
        dir_map.insert(feat_id, dir);
    }

    for &idx in &presence_diff_indices {
        let feat_id = feature_ids[idx].clone();
        let dir = directions.get(&idx).copied().unwrap_or(Direction::Up);
        // Presence effects have "infinite" log2fc, we use a large value
        let log2fc = match dir {
            Direction::Up => 10.0,
            Direction::Down => -10.0,
        };
        true_log2fc.insert(feat_id.clone(), log2fc);
        effect_types.insert(feat_id.clone(), EffectType::Presence);
        dir_map.insert(feat_id, dir);
    }

    let ground_truth = GroundTruth {
        differential_features,
        true_log2fc,
        effect_types,
        directions: dir_map,
        null_features,
    };

    // Build count matrix
    let mut tri_mat = sprs::TriMat::new((n_features, n_samples));
    for feat_idx in 0..n_features {
        for sample_idx in 0..n_samples {
            let count = counts[feat_idx][sample_idx];
            if count > 0 {
                tri_mat.add_triplet(feat_idx, sample_idx, count);
            }
        }
    }
    let count_matrix = CountMatrix::new(tri_mat.to_csr(), feature_ids, sample_ids.clone())?;

    // Build metadata
    let metadata = build_metadata(&sample_ids, &groups, &library_sizes)?;

    Ok(SyntheticData {
        counts: count_matrix,
        metadata,
        ground_truth,
        config: config.clone(),
    })
}

/// Generate abundance profile with specified dominance.
fn generate_abundance_profile(n_features: usize, dominance: f64, rng: &mut Rng) -> Vec<f64> {
    // Use log-normal distribution, then scale top features for dominance
    let mut abundances: Vec<f64> = (0..n_features)
        .map(|_| rng.next_normal(0.0, 2.0).exp())
        .collect();

    // Sort descending
    abundances.sort_by(|a, b| b.partial_cmp(a).unwrap());

    // Adjust top 3 to achieve target dominance
    let total: f64 = abundances.iter().sum();
    let current_top3: f64 = abundances[..3.min(n_features)].iter().sum();
    let current_dominance = current_top3 / total;

    if dominance > current_dominance && dominance < 1.0 {
        // Boost top 3
        let boost_factor = (dominance * total - current_top3 * (1.0 - dominance / current_dominance))
            / (current_top3 * (1.0 - dominance / current_dominance + 1.0));
        for i in 0..3.min(n_features) {
            abundances[i] *= 1.0 + boost_factor.max(0.0);
        }
    }

    // Normalize to proportions
    let total: f64 = abundances.iter().sum();
    abundances.iter().map(|&x| x / total).collect()
}

/// Build metadata from vectors.
fn build_metadata(
    sample_ids: &[String],
    groups: &[String],
    library_sizes: &[f64],
) -> Result<Metadata> {
    let mut content = String::from("sample_id\tgroup\tlibrary_size\n");
    for (i, sample_id) in sample_ids.iter().enumerate() {
        content.push_str(&format!(
            "{}\t{}\t{}\n",
            sample_id,
            groups[i],
            library_sizes[i] as u64
        ));
    }

    // Write to a unique temporary file
    use std::time::{SystemTime, UNIX_EPOCH};
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let temp_path = std::env::temp_dir().join(format!(
        "daa_meta_{}_{}.tsv",
        std::process::id(),
        timestamp
    ));
    std::fs::write(&temp_path, &content)?;
    let metadata = Metadata::from_tsv(&temp_path)?;
    let _ = std::fs::remove_file(&temp_path); // Clean up

    Ok(metadata)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = SyntheticConfig::default();
        assert_eq!(config.n_features, 200);
        assert_eq!(config.n_samples_per_group, 20);
        assert_eq!(config.seed, 42);
    }

    #[test]
    fn test_preset_configs() {
        let ideal = SyntheticConfig::ideal();
        assert!(ideal.sparsity < 0.3);
        assert!(ideal.dominance < 0.3);

        let sparse = SyntheticConfig::sparse_virome();
        assert!(sparse.sparsity > 0.8);

        let extreme = SyntheticConfig::extreme_sparse();
        assert!(extreme.sparsity > 0.9);
    }

    #[test]
    fn test_generate_synthetic() {
        let config = SyntheticConfig::new("test")
            .with_dimensions(50, 10)
            .with_differential(5, 1.0)
            .with_seed(12345);

        let data = generate_synthetic(&config).unwrap();

        assert_eq!(data.counts.n_features(), 50);
        assert_eq!(data.counts.n_samples(), 20);
        assert_eq!(data.ground_truth.differential_features.len(), 5);
        assert_eq!(data.ground_truth.null_features.len(), 45);
    }

    #[test]
    fn test_deterministic_generation() {
        let config = SyntheticConfig::new("test")
            .with_dimensions(50, 10)
            .with_seed(999);

        let data1 = generate_synthetic(&config).unwrap();
        let data2 = generate_synthetic(&config).unwrap();

        // Same seed should produce identical results
        assert_eq!(
            data1.ground_truth.differential_features,
            data2.ground_truth.differential_features
        );

        // Check a few counts are identical
        for feat in 0..5 {
            for sample in 0..5 {
                assert_eq!(
                    data1.counts.get(feat, sample),
                    data2.counts.get(feat, sample)
                );
            }
        }
    }

    #[test]
    fn test_sparsity_control() {
        let low_sparse = generate_synthetic(
            &SyntheticConfig::new("low")
                .with_dimensions(100, 10)
                .with_sparsity(0.0) // No forced zeros
                .with_seed(1)
        ).unwrap();

        let high_sparse = generate_synthetic(
            &SyntheticConfig::new("high")
                .with_dimensions(100, 10)
                .with_sparsity(0.9) // 90% forced zeros
                .with_seed(1)
        ).unwrap();

        // Count zeros
        let count_zeros = |data: &SyntheticData| -> f64 {
            let mut zeros = 0;
            let mut total = 0;
            for f in 0..data.counts.n_features() {
                for s in 0..data.counts.n_samples() {
                    if data.counts.get(f, s) == 0 {
                        zeros += 1;
                    }
                    total += 1;
                }
            }
            zeros as f64 / total as f64
        };

        let low_sparsity = count_zeros(&low_sparse);
        let high_sparsity = count_zeros(&high_sparse);

        // Main assertion: high sparsity config should produce significantly more zeros
        assert!(
            high_sparsity > low_sparsity + 0.2,
            "High sparsity ({:.2}) should be much greater than low sparsity ({:.2})",
            high_sparsity, low_sparsity
        );

        // High sparsity config should produce >90% zeros (since we set 90% zero prob)
        assert!(
            high_sparsity > 0.85,
            "High sparsity config should produce >85% zeros, got {:.2}",
            high_sparsity
        );
    }

    #[test]
    fn test_ground_truth_tsv() {
        let config = SyntheticConfig::new("test")
            .with_dimensions(20, 5)
            .with_differential(3, 1.0)
            .with_seed(42);

        let data = generate_synthetic(&config).unwrap();
        let tsv = data.ground_truth.to_tsv();

        assert!(tsv.contains("feature_id\ttrue_log2fc"));
        assert!(tsv.contains("true"));
        assert!(tsv.contains("false"));
    }

    #[test]
    fn test_group_specific_features() {
        let config = SyntheticConfig::new("test")
            .with_dimensions(50, 10)
            .with_differential(10, 1.0)
            .with_group_specific(0.5)
            .with_seed(42);

        let data = generate_synthetic(&config).unwrap();

        let presence_count = data.ground_truth.effect_types
            .values()
            .filter(|&&e| e == EffectType::Presence)
            .count();

        assert!(presence_count >= 3); // Roughly half of 10
        assert!(presence_count <= 7);
    }
}
