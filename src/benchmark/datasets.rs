//! Fetch classic benchmark datasets from Zenodo.
//!
//! Available datasets from MicrobiomeBenchmarkData (Zenodo 6911027):
//! - HMP gingival (16S V13/V35) - body site comparison
//! - Ravel 2011 BV - bacterial vaginosis case-control
//! - Stammler 2016 spike-in - experimental ground truth

use crate::data::{CountMatrix, Metadata};
use crate::error::{DaaError, Result};
use std::fs;
use std::path::{Path, PathBuf};

const BASE_URL: &str = "https://zenodo.org/records/6911027/files";

/// Available benchmark datasets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BenchmarkDataset {
    /// HMP 2012 gingival 16S V13 region (larger, ~300 samples)
    HmpGingivalV13,
    /// HMP 2012 gingival 16S V35 region
    HmpGingivalV35,
    /// HMP 2012 gingival V35 subset (small, ~50 samples, good for quick tests)
    HmpGingivalV35Subset,
    /// HMP 2012 whole metagenome gingival
    HmpGingivalWms,
    /// Ravel 2011 bacterial vaginosis study
    RavelBv,
    /// Stammler 2016 spike-in experiment (has ground truth!)
    StammlerSpikein,
}

impl BenchmarkDataset {
    /// Get the dataset name for display.
    pub fn name(&self) -> &'static str {
        match self {
            Self::HmpGingivalV13 => "HMP_2012_16S_gingival_V13",
            Self::HmpGingivalV35 => "HMP_2012_16S_gingival_V35",
            Self::HmpGingivalV35Subset => "HMP_2012_16S_gingival_V35_subset",
            Self::HmpGingivalWms => "HMP_2012_WMS_gingival",
            Self::RavelBv => "Ravel_2011_16S_BV",
            Self::StammlerSpikein => "Stammler_2016_16S_spikein",
        }
    }

    /// Get a short description.
    pub fn description(&self) -> &'static str {
        match self {
            Self::HmpGingivalV13 => "HMP gingival microbiome, 16S V1-V3 region, ~300 samples",
            Self::HmpGingivalV35 => "HMP gingival microbiome, 16S V3-V5 region",
            Self::HmpGingivalV35Subset => "HMP gingival subset, 50 samples (quick tests)",
            Self::HmpGingivalWms => "HMP gingival whole metagenome shotgun",
            Self::RavelBv => "Bacterial vaginosis case-control study",
            Self::StammlerSpikein => "Spike-in experiment with known ground truth",
        }
    }

    /// Whether this dataset has experimental ground truth.
    pub fn has_ground_truth(&self) -> bool {
        matches!(self, Self::StammlerSpikein)
    }

    /// Get all available datasets.
    pub fn all() -> Vec<Self> {
        vec![
            Self::HmpGingivalV13,
            Self::HmpGingivalV35,
            Self::HmpGingivalV35Subset,
            Self::HmpGingivalWms,
            Self::RavelBv,
            Self::StammlerSpikein,
        ]
    }

    /// Parse from string.
    pub fn from_str(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "hmp_v13" | "hmp-v13" | "hmpv13" => Some(Self::HmpGingivalV13),
            "hmp_v35" | "hmp-v35" | "hmpv35" => Some(Self::HmpGingivalV35),
            "hmp_subset" | "hmp-subset" | "hmpsubset" => Some(Self::HmpGingivalV35Subset),
            "hmp_wms" | "hmp-wms" | "hmpwms" => Some(Self::HmpGingivalWms),
            "ravel" | "ravel_bv" | "bv" => Some(Self::RavelBv),
            "stammler" | "spikein" | "spike-in" => Some(Self::StammlerSpikein),
            _ => None,
        }
    }

    fn count_matrix_filename(&self) -> String {
        format!("{}_count_matrix.tsv", self.name())
    }

    fn metadata_filename(&self) -> String {
        format!("{}_sample_metadata.tsv", self.name())
    }
}

/// Fetched benchmark data.
pub struct FetchedDataset {
    pub dataset: BenchmarkDataset,
    pub counts: CountMatrix,
    pub metadata: Metadata,
    pub cache_dir: PathBuf,
}

/// Fetch a benchmark dataset, using cache if available.
pub fn fetch_dataset(dataset: BenchmarkDataset, cache_dir: Option<&Path>) -> Result<FetchedDataset> {
    let cache = cache_dir
        .map(PathBuf::from)
        .unwrap_or_else(default_cache_dir);

    fs::create_dir_all(&cache)?;

    let counts_file = cache.join(dataset.count_matrix_filename());
    let metadata_file = cache.join(dataset.metadata_filename());

    // Download if not cached
    if !counts_file.exists() {
        eprintln!("Downloading {}...", dataset.count_matrix_filename());
        download_file(
            &format!("{}/{}?download=1", BASE_URL, dataset.count_matrix_filename()),
            &counts_file,
        )?;
    }

    if !metadata_file.exists() {
        eprintln!("Downloading {}...", dataset.metadata_filename());
        download_file(
            &format!("{}/{}?download=1", BASE_URL, dataset.metadata_filename()),
            &metadata_file,
        )?;
    }

    // Load data
    let counts = CountMatrix::from_tsv(&counts_file)?;
    let metadata = Metadata::from_tsv(&metadata_file)?;

    Ok(FetchedDataset {
        dataset,
        counts,
        metadata,
        cache_dir: cache,
    })
}

/// List available datasets with their cache status.
pub fn list_datasets(cache_dir: Option<&Path>) -> Vec<DatasetInfo> {
    let cache = cache_dir
        .map(PathBuf::from)
        .unwrap_or_else(default_cache_dir);

    BenchmarkDataset::all()
        .into_iter()
        .map(|ds| {
            let counts_cached = cache.join(ds.count_matrix_filename()).exists();
            let metadata_cached = cache.join(ds.metadata_filename()).exists();

            DatasetInfo {
                dataset: ds,
                cached: counts_cached && metadata_cached,
            }
        })
        .collect()
}

/// Clear the dataset cache.
pub fn clear_cache(cache_dir: Option<&Path>) -> Result<()> {
    let cache = cache_dir
        .map(PathBuf::from)
        .unwrap_or_else(default_cache_dir);

    if cache.exists() {
        fs::remove_dir_all(&cache)?;
    }
    Ok(())
}

/// Dataset info for listing.
pub struct DatasetInfo {
    pub dataset: BenchmarkDataset,
    pub cached: bool,
}

fn default_cache_dir() -> PathBuf {
    dirs::cache_dir()
        .unwrap_or_else(|| PathBuf::from("."))
        .join("composable-daa")
        .join("datasets")
}

fn download_file(url: &str, dest: &Path) -> Result<()> {
    // Use curl or reqwest-like approach
    // For simplicity, shell out to curl (available on most systems)
    let output = std::process::Command::new("curl")
        .args(["-sL", "-o", dest.to_str().unwrap(), url])
        .output()
        .map_err(|e| DaaError::Pipeline(format!("Failed to run curl: {}", e)))?;

    if !output.status.success() {
        return Err(DaaError::Pipeline(format!(
            "Download failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }

    // Verify file exists and has content
    let meta = fs::metadata(dest)?;
    if meta.len() == 0 {
        return Err(DaaError::Pipeline("Downloaded file is empty".into()));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dataset_parsing() {
        assert_eq!(
            BenchmarkDataset::from_str("hmp_v13"),
            Some(BenchmarkDataset::HmpGingivalV13)
        );
        assert_eq!(
            BenchmarkDataset::from_str("stammler"),
            Some(BenchmarkDataset::StammlerSpikein)
        );
        assert_eq!(
            BenchmarkDataset::from_str("ravel"),
            Some(BenchmarkDataset::RavelBv)
        );
    }

    #[test]
    fn test_all_datasets() {
        let all = BenchmarkDataset::all();
        assert_eq!(all.len(), 6);
    }
}
