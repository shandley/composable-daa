#!/usr/bin/env python3
"""
Taxonomic Harmonization Module for Microbiome Reproducibility Atlas

This module provides functions to:
1. Parse various taxonomy string formats
2. Extract taxonomy at specified level (genus, family, etc.)
3. Aggregate OTU-level data to genus level
4. Validate known positive associations
"""

import pandas as pd
import numpy as np
import re
from collections import defaultdict


def parse_taxonomy_string(name):
    """
    Parse various taxonomy formats into a dictionary.

    Handles:
    - Greengenes: k__Bacteria;p__Firmicutes;c__Clostridia;...;g__Genus;s__;d__denovo123
    - SILVA: Bacteria;Firmicutes;Clostridia;...;Genus
    - Plain: Genus species
    - With denovo ID: Genus_denovo123

    Returns:
        dict with keys: kingdom, phylum, class, order, family, genus, species, otu_id
    """
    result = {
        'kingdom': None,
        'phylum': None,
        'class': None,
        'order': None,
        'family': None,
        'genus': None,
        'species': None,
        'otu_id': None,
        'raw': name
    }

    if not name or pd.isna(name):
        return result

    name = str(name).strip()

    # Extract denovo/OTU ID if present
    denovo_match = re.search(r'd__denovo\d+|denovo\d+|OTU_?\d+|ASV_?\d+', name)
    if denovo_match:
        result['otu_id'] = denovo_match.group()

    # Greengenes format: k__X;p__Y;c__Z;...
    if '__' in name and ';' in name:
        parts = name.split(';')
        level_map = {
            'k': 'kingdom', 'd': 'kingdom',  # domain/kingdom
            'p': 'phylum',
            'c': 'class',
            'o': 'order',
            'f': 'family',
            'g': 'genus',
            's': 'species'
        }
        for part in parts:
            if '__' in part:
                prefix, value = part.split('__', 1)
                prefix = prefix.lower()
                if prefix in level_map and value and value not in ['', 'unclassified', 'Unknown']:
                    # Remove denovo ID from value if present
                    value = re.sub(r';?d__denovo\d+', '', value)
                    if value:
                        result[level_map[prefix]] = value

    # SILVA-like format: Level1;Level2;Level3;... (no prefixes)
    elif ';' in name and '__' not in name:
        parts = [p.strip() for p in name.split(';') if p.strip()]
        levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        for i, part in enumerate(parts):
            if i < len(levels) and part not in ['', 'unclassified', 'Unknown']:
                result[levels[i]] = part

    # Plain name: "Genus species" or just "Genus"
    elif ' ' in name or name[0].isupper():
        parts = name.split()
        if parts:
            result['genus'] = parts[0]
            if len(parts) > 1 and not parts[1].startswith('denovo'):
                result['species'] = parts[1]

    return result


def get_taxonomy_at_level(name, level='genus'):
    """
    Extract taxonomy at specified level.

    Args:
        name: Taxonomy string in any format
        level: 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'

    Returns:
        str: Taxonomy name at level, or 'Unknown' if not available
    """
    parsed = parse_taxonomy_string(name)
    value = parsed.get(level)

    if value:
        return value

    # If genus not found, try to get lowest available level
    if level == 'genus':
        for l in ['species', 'family', 'order', 'class', 'phylum', 'kingdom']:
            if parsed.get(l):
                return f"Unknown_{l[:3]}_{parsed[l]}"

    return 'Unknown'


def aggregate_to_genus(counts_df, method='sum'):
    """
    Aggregate OTU-level counts to genus level.

    Args:
        counts_df: DataFrame with OTUs as index, samples as columns
        method: 'sum' or 'mean'

    Returns:
        DataFrame with genera as index, samples as columns
    """
    # Parse all taxonomy strings
    genus_mapping = {}
    for otu in counts_df.index:
        genus = get_taxonomy_at_level(otu, 'genus')
        genus_mapping[otu] = genus

    # Add genus column and aggregate
    counts_df = counts_df.copy()
    counts_df['genus'] = counts_df.index.map(genus_mapping)

    if method == 'sum':
        genus_df = counts_df.groupby('genus').sum()
    else:
        genus_df = counts_df.groupby('genus').mean()

    # Remove Unknown if it dominates
    if 'Unknown' in genus_df.index:
        unknown_pct = genus_df.loc['Unknown'].sum() / genus_df.values.sum()
        if unknown_pct > 0.5:
            print(f"Warning: {unknown_pct:.1%} of reads map to 'Unknown' genus")

    return genus_df


def search_taxon(results_df, taxon_name, level='genus', case_insensitive=True):
    """
    Search for a taxon in results DataFrame.

    Args:
        results_df: DataFrame with taxonomy as index
        taxon_name: Name to search for (e.g., 'Fusobacterium')
        level: Taxonomy level to search at
        case_insensitive: Whether to ignore case

    Returns:
        DataFrame subset matching the taxon
    """
    if case_insensitive:
        pattern = taxon_name.lower()
        matches = [idx for idx in results_df.index
                   if pattern in get_taxonomy_at_level(idx, level).lower()]
    else:
        pattern = taxon_name
        matches = [idx for idx in results_df.index
                   if pattern in get_taxonomy_at_level(idx, level)]

    return results_df.loc[matches] if matches else pd.DataFrame()


# Known positive associations for validation
KNOWN_POSITIVES = {
    'crc': [
        {'taxon': 'Fusobacterium', 'direction': 'up', 'confidence': 'high',
         'evidence': 'Castellarin 2012, Kostic 2012, multiple meta-analyses'},
        {'taxon': 'Parvimonas', 'direction': 'up', 'confidence': 'high',
         'evidence': 'Multiple CRC metagenomics studies'},
        {'taxon': 'Porphyromonas', 'direction': 'up', 'confidence': 'medium',
         'evidence': 'Oral bacteria enriched in CRC'},
        {'taxon': 'Peptostreptococcus', 'direction': 'up', 'confidence': 'medium',
         'evidence': 'Wirbel 2019, Thomas 2019'},
        {'taxon': 'Solobacterium', 'direction': 'up', 'confidence': 'medium',
         'evidence': 'Thomas 2019 meta-analysis'},
    ],
    'ibd': [
        {'taxon': 'Faecalibacterium', 'direction': 'down', 'confidence': 'high',
         'evidence': 'F. prausnitzii depletion is hallmark of IBD'},
        {'taxon': 'Roseburia', 'direction': 'down', 'confidence': 'high',
         'evidence': 'Butyrate producer, depleted in IBD'},
        {'taxon': 'Escherichia', 'direction': 'up', 'confidence': 'high',
         'evidence': 'Enterobacteriaceae bloom in inflammation'},
        {'taxon': 'Ruminococcus', 'direction': 'down', 'confidence': 'medium',
         'evidence': 'R. gnavus complex findings'},
    ],
    'cdi': [
        {'taxon': 'Clostridioides', 'direction': 'up', 'confidence': 'high',
         'evidence': 'Causative agent'},
        {'taxon': 'Clostridium', 'direction': 'up', 'confidence': 'medium',
         'evidence': 'May include C. difficile reads'},
    ],
}


def validate_known_positives(results_df, disease, threshold=0.1, verbose=True):
    """
    Check if known positive associations are detected.

    Args:
        results_df: DataFrame with taxonomy index and q_value, estimate columns
        disease: Disease name ('crc', 'ibd', 'cdi')
        threshold: q-value threshold for significance
        verbose: Print detailed output

    Returns:
        dict with validation results
    """
    expected = KNOWN_POSITIVES.get(disease.lower(), [])

    if not expected:
        print(f"No known positives defined for {disease}")
        return {}

    validation = {
        'disease': disease,
        'found': 0,
        'correct_direction': 0,
        'significant': 0,
        'total_expected': len(expected),
        'details': []
    }

    for kp in expected:
        taxon = kp['taxon']
        expected_dir = kp['direction']

        # Search for taxon
        matches = search_taxon(results_df, taxon)

        if len(matches) == 0:
            status = 'NOT_FOUND'
            if verbose:
                print(f"✗ {taxon}: Not found in results")
        else:
            validation['found'] += 1

            # Get best (most significant) match
            best = matches.loc[matches['q_value'].idxmin()]
            q = best['q_value']
            effect = best['estimate']

            # Check direction (positive effect = higher in case vs control)
            # Note: direction depends on coefficient coding
            observed_dir = 'up' if effect < 0 else 'down'  # groupcontrol coefficient

            dir_match = (observed_dir == expected_dir)
            if dir_match:
                validation['correct_direction'] += 1

            is_sig = q < threshold
            if is_sig:
                validation['significant'] += 1

            status = '✓' if (is_sig and dir_match) else '~' if dir_match else '✗'

            if verbose:
                dir_symbol = '↑' if observed_dir == 'up' else '↓'
                exp_symbol = '↑' if expected_dir == 'up' else '↓'
                print(f"{status} {taxon}: q={q:.3e}, effect={effect:.2f} "
                      f"(observed: {dir_symbol}, expected: {exp_symbol})")

        validation['details'].append({
            'taxon': taxon,
            'status': status,
            'expected_direction': expected_dir
        })

    if verbose:
        print(f"\nSummary: {validation['found']}/{validation['total_expected']} found, "
              f"{validation['significant']} significant, "
              f"{validation['correct_direction']} correct direction")

    return validation


if __name__ == '__main__':
    # Test parsing
    test_cases = [
        "k__Bacteria;p__Fusobacteria;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Fusobacterium;s__;d__denovo58195",
        "Fusobacterium nucleatum",
        "d__denovo12345",
        "Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Roseburia;R. intestinalis",
    ]

    print("=== Taxonomy Parsing Tests ===\n")
    for tc in test_cases:
        parsed = parse_taxonomy_string(tc)
        genus = get_taxonomy_at_level(tc, 'genus')
        print(f"Input: {tc[:60]}...")
        print(f"Genus: {genus}")
        print(f"Parsed: {parsed}")
        print()
