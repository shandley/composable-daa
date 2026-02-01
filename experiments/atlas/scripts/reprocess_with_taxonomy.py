#!/usr/bin/env python3
"""
Reprocess CRC data with proper taxonomic handling.

This script:
1. Reads raw OTU tables with full taxonomy strings
2. Aggregates to genus level
3. Re-runs DAA analysis
4. Validates with known positive (Fusobacterium)
"""

import pandas as pd
import numpy as np
import os
import glob
import subprocess
import sys

# Add parent directory to path for taxonomy module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taxonomy import (
    parse_taxonomy_string, get_taxonomy_at_level, aggregate_to_genus,
    validate_known_positives, KNOWN_POSITIVES
)


def process_cohort_with_taxonomy(disease_dir, cohort_name, raw_folder):
    """
    Process a cohort with proper taxonomic handling.
    """
    print(f"\n{'='*60}")
    print(f"Processing {cohort_name} with taxonomy")
    print('='*60)

    base_path = os.path.join(disease_dir, 'data', raw_folder)

    # Find OTU table (RDP classified version has full taxonomy)
    otu_files = glob.glob(os.path.join(base_path, 'RDP', '*.rdp_assigned'))
    if not otu_files:
        otu_files = glob.glob(os.path.join(base_path, '*.otu_table.*'))
    if not otu_files:
        print(f"  No OTU table found in {base_path}")
        return None

    otu_file = otu_files[0]
    print(f"  OTU file: {os.path.basename(otu_file)}")

    # Find metadata
    meta_files = glob.glob(os.path.join(base_path, '*.metadata.txt'))
    if not meta_files:
        print(f"  No metadata found")
        return None

    meta_file = meta_files[0]

    # Read OTU table - keep full taxonomy in index
    otu_df = pd.read_csv(otu_file, sep='\t', index_col=0)
    print(f"  Raw OTUs: {otu_df.shape[0]} features x {otu_df.shape[1]} samples")

    # Check for Fusobacterium in raw data
    fuso_otus = [idx for idx in otu_df.index if 'Fusobacterium' in str(idx)]
    print(f"  Fusobacterium OTUs in raw data: {len(fuso_otus)}")
    if fuso_otus:
        for f in fuso_otus[:3]:
            total_reads = otu_df.loc[f].sum()
            print(f"    {f[:60]}... ({total_reads} reads)")

    # Read metadata
    try:
        meta_df = pd.read_csv(meta_file, sep='\t')
    except:
        meta_df = pd.read_csv(meta_file, sep='\t', encoding='latin-1')

    # Find sample and disease columns
    sample_col = None
    for col in ['#SampleID', 'SampleID', 'sample_id', 'Sample', 'sample', 'Subject ID']:
        if col in meta_df.columns:
            sample_col = col
            break
    if not sample_col:
        sample_col = meta_df.columns[0]

    disease_col = 'DiseaseState' if 'DiseaseState' in meta_df.columns else None
    if not disease_col:
        for col in meta_df.columns:
            if 'disease' in col.lower():
                disease_col = col
                break

    if not disease_col:
        print(f"  No disease column found")
        return None

    print(f"  Sample col: {sample_col}, Disease col: {disease_col}")
    print(f"  Disease values: {meta_df[disease_col].unique()}")

    # Classify groups for CRC
    def classify_crc(x):
        x = str(x).lower().strip()
        if x in ['crc', 'cancer', 'carcinoma', 'tumor']:
            return 'case'
        if any(k in x for k in ['crc', 'cancer', 'carcinoma']):
            return 'case'
        if x in ['h', 'healthy', 'control', 'normal', 'noncrc', 'non-crc']:
            return 'control'
        if 'control' in x or 'healthy' in x or 'normal' in x or 'non' in x:
            return 'control'
        return 'exclude'

    meta_df['group'] = meta_df[disease_col].astype(str).apply(classify_crc)
    print(f"  Group counts: {meta_df['group'].value_counts().to_dict()}")

    meta_df = meta_df[meta_df['group'].isin(['case', 'control'])]

    if len(meta_df) < 20:
        print(f"  Not enough samples")
        return None

    # Match samples
    meta_out = pd.DataFrame({
        'sample_id': meta_df[sample_col].values,
        'group': meta_df['group'].values,
        'cohort': cohort_name
    })

    common = list(set(otu_df.columns) & set(meta_out['sample_id']))
    if len(common) < 20:
        # Try alternate sample ID column
        if 'Subject ID' in meta_df.columns and sample_col != 'Subject ID':
            meta_out['sample_id'] = meta_df['Subject ID'].values
            common = list(set(otu_df.columns) & set(meta_out['sample_id']))

    if len(common) < 20:
        print(f"  Only {len(common)} matching samples")
        return None

    otu_matched = otu_df[common]
    meta_out = meta_out[meta_out['sample_id'].isin(common)]

    print(f"  Matched: {otu_matched.shape[0]} OTUs x {otu_matched.shape[1]} samples")

    # Aggregate to genus level
    print(f"\n  Aggregating to genus level...")
    genus_df = aggregate_to_genus(otu_matched)
    print(f"  Genus-level: {genus_df.shape[0]} genera x {genus_df.shape[1]} samples")

    # Check for Fusobacterium at genus level
    if 'Fusobacterium' in genus_df.index:
        fuso_counts = genus_df.loc['Fusobacterium']
        fuso_total = fuso_counts.sum()
        fuso_prev = (fuso_counts > 0).mean()
        print(f"  Fusobacterium: {fuso_total} total reads, {fuso_prev:.1%} prevalence")

        # Check case vs control
        case_samples = meta_out[meta_out['group'] == 'case']['sample_id'].values
        ctrl_samples = meta_out[meta_out['group'] == 'control']['sample_id'].values

        case_mean = genus_df.loc['Fusobacterium', [s for s in case_samples if s in genus_df.columns]].mean()
        ctrl_mean = genus_df.loc['Fusobacterium', [s for s in ctrl_samples if s in genus_df.columns]].mean()
        print(f"  Fusobacterium mean: CRC={case_mean:.1f}, Control={ctrl_mean:.1f}")
    else:
        print(f"  Warning: Fusobacterium not found at genus level!")

    # Convert to integers
    genus_df = genus_df.apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)

    # Save genus-level data
    output_dir = os.path.join(disease_dir, 'data_genus')
    os.makedirs(output_dir, exist_ok=True)

    genus_df.to_csv(os.path.join(output_dir, f'{cohort_name}_counts.tsv'), sep='\t')
    meta_out.to_csv(os.path.join(output_dir, f'{cohort_name}_metadata.tsv'), sep='\t', index=False)

    print(f"  Saved to {output_dir}/")

    return genus_df, meta_out


def run_daa_on_genus_data(disease_dir, cohort_name):
    """Run DAA on genus-level data using daa run with YAML config."""
    print(f"\n  Running DAA on {cohort_name}...")

    data_dir = os.path.join(disease_dir, 'data_genus')
    results_dir = os.path.join(disease_dir, 'results_genus')
    os.makedirs(results_dir, exist_ok=True)

    daa_bin = os.path.expanduser('~/.cargo/target/release/daa')

    counts_file = os.path.join(data_dir, f'{cohort_name}_counts.tsv')
    meta_file = os.path.join(data_dir, f'{cohort_name}_metadata.tsv')

    if not os.path.exists(counts_file):
        print(f"    No counts file: {counts_file}")
        return None

    # Create YAML config for Hurdle model
    config_file = os.path.join(results_dir, f'{cohort_name}_hurdle_config.yaml')
    config_content = """name: hurdle_genus
steps:
  - !FilterPrevalence
    threshold: 0.05
  - !ModelHurdle
    formula: "~ group"
  - !TestWald
    coefficient: groupcontrol
  - CorrectBH
"""
    with open(config_file, 'w') as f:
        f.write(config_content)

    result_file = os.path.join(results_dir, f'{cohort_name}_hurdle.tsv')

    cmd = [
        daa_bin, 'run',
        '-c', counts_file,
        '-m', meta_file,
        '--config', config_file,
        '-o', result_file
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            print(f"    Saved to {result_file}")
            return result_file
        else:
            print(f"    Error: {result.stderr[:200]}")
            return None
    except Exception as e:
        print(f"    Exception: {e}")
        return None


def main():
    disease_dir = os.path.join(os.path.dirname(__file__), '..', 'diseases', 'crc')
    disease_dir = os.path.abspath(disease_dir)

    print("="*70)
    print("REPROCESSING CRC DATA WITH PROPER TAXONOMY")
    print("="*70)

    # Define cohorts
    cohorts = [
        ('zeller', 'crc_zeller'),
        ('zackular', 'crc_zackular'),
        ('xiang', 'crc_xiang'),
        ('zhao', 'crc_zhao'),
    ]

    # Process each cohort
    successful = []
    for cohort_name, raw_folder in cohorts:
        result = process_cohort_with_taxonomy(disease_dir, cohort_name, raw_folder)
        if result is not None:
            successful.append(cohort_name)

    print(f"\n\nSuccessfully processed: {successful}")

    # Run DAA on genus-level data
    print("\n" + "="*70)
    print("RUNNING DAA ON GENUS-LEVEL DATA")
    print("="*70)

    for cohort_name in successful:
        run_daa_on_genus_data(disease_dir, cohort_name)

    # Validate known positives
    print("\n" + "="*70)
    print("VALIDATING KNOWN POSITIVES (Fusobacterium)")
    print("="*70)

    results_dir = os.path.join(disease_dir, 'results_genus')

    for cohort_name in successful:
        result_file = os.path.join(results_dir, f'{cohort_name}_hurdle.tsv')
        if os.path.exists(result_file):
            print(f"\n{cohort_name}:")
            df = pd.read_csv(result_file, sep='\t')
            if 'feature_id' in df.columns:
                df = df.set_index('feature_id')
            elif 'feature' in df.columns:
                df = df.set_index('feature')

            validate_known_positives(df, 'crc', threshold=0.10)

    # Cross-cohort check for Fusobacterium
    print("\n" + "="*70)
    print("FUSOBACTERIUM ACROSS ALL COHORTS")
    print("="*70)

    fuso_results = []
    for cohort_name in successful:
        result_file = os.path.join(results_dir, f'{cohort_name}_hurdle.tsv')
        if os.path.exists(result_file):
            df = pd.read_csv(result_file, sep='\t')
            if 'feature_id' in df.columns:
                df = df.set_index('feature_id')
            elif 'feature' in df.columns:
                df = df.set_index('feature')

            # Search for Fusobacterium
            fuso_rows = df[df.index.str.contains('Fusobacterium', case=False, na=False)]
            if len(fuso_rows) > 0:
                best = fuso_rows.loc[fuso_rows['q_value'].idxmin()]
                fuso_results.append({
                    'cohort': cohort_name,
                    'taxon': fuso_rows.index[0],
                    'estimate': best['estimate'],
                    'q_value': best['q_value'],
                    'significant': best['q_value'] < 0.10
                })

    if fuso_results:
        print("\nFusobacterium results:")
        for r in fuso_results:
            sig = "✓ SIG" if r['significant'] else "  ns"
            print(f"  {r['cohort']}: effect={r['estimate']:.2f}, q={r['q_value']:.3e} {sig}")

        n_sig = sum(1 for r in fuso_results if r['significant'])
        print(f"\nFusobacterium significant in {n_sig}/{len(fuso_results)} cohorts")

        # Direction consistency
        directions = ['up' if r['estimate'] < 0 else 'down' for r in fuso_results]
        consistent = len(set(directions)) == 1
        print(f"Direction consistent: {consistent} (all {directions[0] if consistent else 'mixed'})")
    else:
        print("Fusobacterium not found in any cohort!")


if __name__ == '__main__':
    main()
