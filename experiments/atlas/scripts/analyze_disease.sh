#!/bin/bash
#
# Microbiome Reproducibility Atlas - Disease Analysis Script
#
# Usage: ./analyze_disease.sh <disease> <cohort1> <cohort2> ...
# Example: ./analyze_disease.sh crc baxter zeller zackular xiang zhao zhu
#
# This script:
# 1. Downloads cohort data from MicrobiomeHD
# 2. Converts to standard format
# 3. Runs LinDA, Hurdle, and Permutation tests
# 4. Calculates cross-cohort replication scores

set -e

DISEASE=$1
shift
COHORTS=("$@")

if [ -z "$DISEASE" ] || [ ${#COHORTS[@]} -eq 0 ]; then
    echo "Usage: $0 <disease> <cohort1> <cohort2> ..."
    echo "Example: $0 crc baxter zeller zackular"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ATLAS_DIR="$(dirname "$SCRIPT_DIR")"
DISEASE_DIR="$ATLAS_DIR/diseases/$DISEASE"
DATA_DIR="$DISEASE_DIR/data"
RESULTS_DIR="$DISEASE_DIR/results"
CONFIGS_DIR="$DISEASE_DIR/configs"
DAA_BIN="${HOME}/.cargo/target/release/daa"

echo "=============================================="
echo "Microbiome Reproducibility Atlas"
echo "Disease: $DISEASE"
echo "Cohorts: ${COHORTS[*]}"
echo "=============================================="

# Create directories
mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$CONFIGS_DIR"

# Ensure daa binary exists
if [ ! -f "$DAA_BIN" ]; then
    echo "Building daa binary..."
    (cd "$ATLAS_DIR/../.." && cargo build --release --bin daa)
fi

# ============================================
# PHASE 1: Download Datasets
# ============================================

echo ""
echo "=== Phase 1: Downloading Datasets ==="

BASE_URL="https://zenodo.org/records/840333/files"

download_cohort() {
    local disease=$1
    local cohort=$2
    local filename="${disease}_${cohort}_results.tar.gz"

    if [ ! -d "$DATA_DIR/${disease}_${cohort}" ]; then
        echo "Downloading ${cohort}..."
        curl -sL "${BASE_URL}/${filename}?download=1" -o "$DATA_DIR/${cohort}.tar.gz" 2>/dev/null || {
            echo "  Warning: Could not download $filename"
            return 1
        }
        mkdir -p "$DATA_DIR/${disease}_${cohort}"
        tar -xzf "$DATA_DIR/${cohort}.tar.gz" -C "$DATA_DIR/${disease}_${cohort}" --strip-components=1 2>/dev/null || {
            echo "  Warning: Could not extract $filename"
            rm -f "$DATA_DIR/${cohort}.tar.gz"
            return 1
        }
        rm "$DATA_DIR/${cohort}.tar.gz"
        echo "  Downloaded ${cohort}"
    else
        echo "  ${cohort} already present"
    fi
}

for cohort in "${COHORTS[@]}"; do
    download_cohort "$DISEASE" "$cohort"
done

# ============================================
# PHASE 2: Convert to Standard Format
# ============================================

echo ""
echo "=== Phase 2: Converting to Standard Format ==="

python3 << PYTHON_SCRIPT
import pandas as pd
import numpy as np
import os
import glob
import sys

data_dir = "$DATA_DIR"
disease = "$DISEASE"
cohorts = "${COHORTS[*]}".split()

def process_microbiomehd_dataset(disease, cohort):
    """Process a MicrobiomeHD dataset to standard format."""
    name = f"{disease}_{cohort}"
    print(f"\nProcessing {name}...")

    base_path = os.path.join(data_dir, name)

    if not os.path.exists(base_path):
        print(f"  Warning: Directory not found: {base_path}")
        return None, None

    # Find OTU table
    otu_files = glob.glob(os.path.join(base_path, 'RDP', '*.otu_table.*'))
    if not otu_files:
        otu_files = glob.glob(os.path.join(base_path, '*.otu_table.*'))
    if not otu_files:
        print(f"  Warning: No OTU table found")
        return None, None

    otu_file = otu_files[0]
    print(f"  OTU file: {os.path.basename(otu_file)}")

    # Find metadata
    meta_files = glob.glob(os.path.join(base_path, '*.metadata.txt'))
    if not meta_files:
        print(f"  Warning: No metadata found")
        return None, None

    meta_file = meta_files[0]
    print(f"  Metadata: {os.path.basename(meta_file)}")

    # Read OTU table
    try:
        otu_df = pd.read_csv(otu_file, sep='\t', index_col=0)
    except Exception as e:
        print(f"  Error reading OTU table: {e}")
        return None, None

    # Read metadata with encoding fallback
    try:
        meta_df = pd.read_csv(meta_file, sep='\t')
    except UnicodeDecodeError:
        meta_df = pd.read_csv(meta_file, sep='\t', encoding='latin-1')

    # Find sample ID column
    sample_col = None
    for col in ['#SampleID', 'SampleID', 'sample_id', 'Sample', 'Subject ID']:
        if col in meta_df.columns:
            sample_col = col
            break
    if sample_col is None:
        sample_col = meta_df.columns[0]

    # Find disease state column
    disease_col = None
    for col in ['DiseaseState', 'disease_state', 'Disease', 'disease', 'Group', 'group', 'Status']:
        if col in meta_df.columns:
            disease_col = col
            break

    if disease_col is None:
        for col in meta_df.columns:
            vals = meta_df[col].astype(str).str.lower().unique()
            if any(k in str(vals) for k in ['case', 'control', 'healthy', 'disease', 'h', 'patient']):
                disease_col = col
                break

    if disease_col is None:
        print(f"  Warning: No disease state column found")
        print(f"  Columns: {list(meta_df.columns)}")
        return None, None

    print(f"  Sample column: {sample_col}")
    print(f"  Disease column: {disease_col}")
    print(f"  Values: {meta_df[disease_col].unique()[:5]}")

    # Standardize disease labels based on disease type
    meta_df['group_raw'] = meta_df[disease_col].astype(str).str.lower().str.strip()

    # Generic disease/control mapping
    def classify_group(x):
        x = str(x).lower().strip()
        # Control patterns
        if any(k in x for k in ['control', 'healthy', 'normal', 'non-', 'negative']):
            return 'control'
        if x == 'h' or x == 'n':
            return 'control'
        # Disease patterns (be generous - anything not control is disease)
        if any(k in x for k in ['case', 'patient', 'disease', 'positive', 'affected']):
            return 'case'
        # Disease-specific patterns
        if disease == 'crc' and any(k in x for k in ['crc', 'cancer', 'carcinoma', 'tumor']):
            return 'case'
        if disease == 'ibd' and any(k in x for k in ['ibd', 'crohn', 'uc', 'colitis', 'cd']):
            return 'case'
        if disease == 'obesity' or disease == 'ob':
            if any(k in x for k in ['obese', 'overweight', 'ob']):
                return 'case'
            if any(k in x for k in ['lean', 'normal']):
                return 'control'
        if disease == 'cdi' and any(k in x for k in ['cdi', 'cdiff', 'c. diff', 'positive']):
            return 'case'
        if disease == 'hiv' and any(k in x for k in ['hiv', 'positive', 'infected']):
            return 'case'
        if disease == 't1d' and any(k in x for k in ['t1d', 'diabetic', 'diabetes']):
            return 'case'
        # Exclude ambiguous categories
        if any(k in x for k in ['adenoma', 'intermediate', 'unknown', 'other']):
            return 'exclude'
        return 'unknown'

    meta_df['group'] = meta_df['group_raw'].apply(classify_group)

    # Keep only case and control
    meta_df = meta_df[meta_df['group'].isin(['case', 'control'])].copy()

    if len(meta_df) < 20:
        print(f"  Warning: Only {len(meta_df)} case/control samples, skipping")
        return None, None

    # Create standardized metadata
    meta_out = pd.DataFrame({
        'sample_id': meta_df[sample_col].values,
        'group': meta_df['group'].values,
        'cohort': cohort
    })

    # Match samples between OTU and metadata
    otu_samples = set(otu_df.columns)
    meta_samples = set(meta_out['sample_id'])
    common_samples = list(otu_samples & meta_samples)

    if len(common_samples) < 20:
        # Try using Subject ID if Sample ID didn't work
        if 'Subject ID' in meta_df.columns and sample_col != 'Subject ID':
            meta_out['sample_id'] = meta_df['Subject ID'].values
            meta_samples = set(meta_out['sample_id'])
            common_samples = list(otu_samples & meta_samples)

    if len(common_samples) < 20:
        print(f"  Warning: Only {len(common_samples)} matching samples")
        return None, None

    otu_filtered = otu_df[common_samples]
    meta_out = meta_out[meta_out['sample_id'].isin(common_samples)]

    # Convert to integers
    otu_filtered = otu_filtered.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    otu_filtered.index.name = 'feature_id'

    # Clean taxon names
    def clean_taxon(name):
        if ';' in str(name):
            parts = str(name).split(';')
            for p in reversed(parts):
                if p.strip() and p.strip() != 'unclassified':
                    return p.strip()
        return str(name)

    otu_filtered.index = [clean_taxon(n) for n in otu_filtered.index]

    print(f"  Final: {otu_filtered.shape[0]} features x {otu_filtered.shape[1]} samples")
    print(f"  Groups: {meta_out['group'].value_counts().to_dict()}")

    # Save
    otu_filtered.to_csv(os.path.join(data_dir, f'{cohort}_counts.tsv'), sep='\t')
    meta_out.to_csv(os.path.join(data_dir, f'{cohort}_metadata.tsv'), sep='\t', index=False)

    return otu_filtered, meta_out

# Process all cohorts
successful = []
for cohort in cohorts:
    otu, meta = process_microbiomehd_dataset(disease, cohort)
    if otu is not None:
        successful.append(cohort)

print(f"\nSuccessfully processed: {len(successful)}/{len(cohorts)} cohorts")
print(f"Cohorts: {successful}")

# Save list of successful cohorts
with open(os.path.join(data_dir, 'successful_cohorts.txt'), 'w') as f:
    for c in successful:
        f.write(f"{c}\n")

PYTHON_SCRIPT

# ============================================
# PHASE 3: Create Pipeline Configs
# ============================================

echo ""
echo "=== Phase 3: Creating Pipeline Configs ==="

# LinDA pipeline
cat > "$CONFIGS_DIR/linda.yaml" << 'EOF'
name: linda
description: LinDA-style pipeline (CLR + Linear Model)
steps:
- !FilterPrevalence
  threshold: 0.10
- !AddPseudocount
  value: 0.5
- NormalizeCLR
- !ModelLM
  formula: "~ group"
- !TestWald
  coefficient: groupcontrol
- CorrectBH
EOF

# Hurdle pipeline
cat > "$CONFIGS_DIR/hurdle.yaml" << 'EOF'
name: hurdle
description: Hurdle model pipeline
steps:
- !FilterPrevalence
  threshold: 0.10
- !ModelHurdle
  formula: "~ group"
- !TestWald
  coefficient: groupcontrol
- CorrectBH
EOF

echo "Created pipeline configs"

# ============================================
# PHASE 4: Run Analysis on Each Cohort
# ============================================

echo ""
echo "=== Phase 4: Running DAA on Each Cohort ==="

# Read successful cohorts
if [ ! -f "$DATA_DIR/successful_cohorts.txt" ]; then
    echo "No cohorts to analyze"
    exit 1
fi

while IFS= read -r cohort; do
    counts="$DATA_DIR/${cohort}_counts.tsv"
    metadata="$DATA_DIR/${cohort}_metadata.tsv"

    if [ ! -f "$counts" ]; then
        echo "Skipping $cohort (no counts file)"
        continue
    fi

    echo ""
    echo "Analyzing $cohort..."

    # LinDA
    echo "  Running LinDA..."
    "$DAA_BIN" run \
        -c "$counts" \
        -m "$metadata" \
        --config "$CONFIGS_DIR/linda.yaml" \
        -o "$RESULTS_DIR/${cohort}_linda.tsv" 2>&1 | head -5 || echo "  (LinDA completed with warnings)"

    # Hurdle
    echo "  Running Hurdle..."
    "$DAA_BIN" run \
        -c "$counts" \
        -m "$metadata" \
        --config "$CONFIGS_DIR/hurdle.yaml" \
        -o "$RESULTS_DIR/${cohort}_hurdle.tsv" 2>&1 | head -5 || echo "  (Hurdle completed with warnings)"

    # Permutation
    echo "  Running Permutation..."
    "$DAA_BIN" permutation \
        -c "$counts" \
        -m "$metadata" \
        -f "~ group" \
        -t groupcontrol \
        --prevalence 0.10 \
        --n-permutations 999 \
        -o "$RESULTS_DIR/${cohort}_permutation.tsv" 2>&1 | head -5 || echo "  (Permutation completed with warnings)"

done < "$DATA_DIR/successful_cohorts.txt"

# ============================================
# PHASE 5: Calculate Replication Scores
# ============================================

echo ""
echo "=== Phase 5: Calculating Replication Scores ==="

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import numpy as np
import os
from collections import Counter

results_dir = "$RESULTS_DIR"
data_dir = "$DATA_DIR"
disease = "$DISEASE"

print(f"Calculating replication scores for {disease}")
print("=" * 60)

# Read successful cohorts
cohorts_file = os.path.join(data_dir, 'successful_cohorts.txt')
if not os.path.exists(cohorts_file):
    print("No cohorts file found")
    exit(1)

with open(cohorts_file) as f:
    cohorts = [line.strip() for line in f if line.strip()]

print(f"Cohorts: {cohorts}")

# Load results
LINDA_THRESHOLD = 0.10
HURDLE_THRESHOLD = 0.05

results = {}
for cohort in cohorts:
    results[cohort] = {}
    for method in ['linda', 'hurdle']:
        filepath = os.path.join(results_dir, f'{cohort}_{method}.tsv')
        if os.path.exists(filepath):
            try:
                df = pd.read_csv(filepath, sep='\t')
                if 'feature_id' in df.columns:
                    df = df.drop_duplicates(subset='feature_id', keep='first')
                    df = df.set_index('feature_id')
                elif 'feature' in df.columns:
                    df = df.drop_duplicates(subset='feature', keep='first')
                    df = df.set_index('feature')
                results[cohort][method] = df
            except Exception as e:
                print(f"Warning: Could not load {cohort}_{method}: {e}")

# Count significant taxa per cohort
print("\n=== Per-Cohort Summary ===\n")
for cohort in cohorts:
    for method in ['linda', 'hurdle']:
        if method in results.get(cohort, {}):
            df = results[cohort][method]
            threshold = LINDA_THRESHOLD if method == 'linda' else HURDLE_THRESHOLD
            n_sig = (df['q_value'] < threshold).sum()
            print(f"{cohort} {method} (q<{threshold}): {n_sig} significant")

# Calculate replication scores
print("\n=== Replication Scores ===\n")

replication_data = []

for method in ['linda', 'hurdle']:
    threshold = LINDA_THRESHOLD if method == 'linda' else HURDLE_THRESHOLD

    # Get all significant taxa from all cohorts
    all_sig_taxa = {}
    for cohort in cohorts:
        if method in results.get(cohort, {}):
            df = results[cohort][method]
            sig = df[df['q_value'] < threshold].index.tolist()
            all_sig_taxa[cohort] = set(sig)

    if len(all_sig_taxa) < 2:
        print(f"{method}: Not enough cohorts")
        continue

    # Count how many cohorts each taxon is significant in
    all_taxa = set()
    for taxa_set in all_sig_taxa.values():
        all_taxa.update(taxa_set)

    for taxon in all_taxa:
        n_sig = sum(1 for c in cohorts if taxon in all_sig_taxa.get(c, set()))
        n_tested = len([c for c in cohorts if method in results.get(c, {})])
        score = n_sig / n_tested if n_tested > 0 else 0

        # Get effect directions
        effects = []
        for c in cohorts:
            if method in results.get(c, {}) and taxon in results[c][method].index:
                est = results[c][method].loc[taxon, 'estimate']
                effects.append(est)

        direction_consistent = all(e > 0 for e in effects) or all(e < 0 for e in effects) if effects else False

        replication_data.append({
            'taxon': taxon,
            'method': method,
            'n_significant': n_sig,
            'n_tested': n_tested,
            'replication_score': score,
            'direction_consistent': direction_consistent,
            'mean_effect': np.mean(effects) if effects else np.nan,
            'effect_range': max(effects) - min(effects) if len(effects) > 1 else 0
        })

# Create summary DataFrame
rep_df = pd.DataFrame(replication_data)

if len(rep_df) > 0:
    # Summary statistics
    for method in ['linda', 'hurdle']:
        method_df = rep_df[rep_df['method'] == method]
        if len(method_df) == 0:
            continue

        n_cohorts = method_df['n_tested'].max()

        print(f"\n{method.upper()} (n={n_cohorts} cohorts):")

        for n in range(1, n_cohorts + 1):
            count = len(method_df[method_df['n_significant'] == n])
            pct = 100 * count / len(method_df) if len(method_df) > 0 else 0
            label = "ROBUST" if n == n_cohorts else ("moderate" if n >= n_cohorts/2 else "study-specific")
            print(f"  {n}/{n_cohorts} cohorts: {count} taxa ({pct:.1f}%) - {label}")

        # Show top robust taxa
        robust = method_df[method_df['n_significant'] == n_cohorts].sort_values('mean_effect', key=abs, ascending=False)
        if len(robust) > 0:
            print(f"\n  Top robust taxa (all {n_cohorts} cohorts):")
            for _, row in robust.head(10).iterrows():
                direction = "+" if row['mean_effect'] > 0 else "-"
                print(f"    {row['taxon'][:40]}: {direction} (score={row['replication_score']:.2f})")

    # Save results
    rep_df.to_csv(os.path.join(results_dir, 'replication_scores.tsv'), sep='\t', index=False)
    print(f"\nSaved replication scores to {results_dir}/replication_scores.tsv")

    # Summary by replication level
    summary = []
    for method in ['linda', 'hurdle']:
        method_df = rep_df[rep_df['method'] == method]
        if len(method_df) == 0:
            continue
        n_cohorts = method_df['n_tested'].max()

        summary.append({
            'disease': disease,
            'method': method,
            'n_cohorts': n_cohorts,
            'total_taxa_sig': len(method_df),
            'taxa_all_cohorts': len(method_df[method_df['n_significant'] == n_cohorts]),
            'taxa_majority': len(method_df[method_df['n_significant'] >= n_cohorts/2]),
            'taxa_one_only': len(method_df[method_df['n_significant'] == 1]),
            'replication_rate': len(method_df[method_df['n_significant'] == n_cohorts]) / len(method_df) if len(method_df) > 0 else 0
        })

    pd.DataFrame(summary).to_csv(os.path.join(results_dir, 'summary.tsv'), sep='\t', index=False)

    print("\n=== SUMMARY ===")
    for s in summary:
        rate = s['replication_rate'] * 100
        print(f"{s['method'].upper()}: {s['taxa_all_cohorts']}/{s['total_taxa_sig']} taxa replicate in all {s['n_cohorts']} cohorts ({rate:.1f}%)")

else:
    print("No replication data generated")

PYTHON_SCRIPT

echo ""
echo "=============================================="
echo "Analysis Complete for $DISEASE"
echo "=============================================="
echo ""
echo "Results saved to: $RESULTS_DIR"
ls -la "$RESULTS_DIR"/*.tsv 2>/dev/null || echo "(no results files)"
