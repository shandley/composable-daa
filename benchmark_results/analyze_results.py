#!/usr/bin/env python3
"""Analyze benchmark results against ground truth."""
import csv
from pathlib import Path

def load_ground_truth(path):
    """Load true differential features."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        return {row['feature_id']: float(row['true_log2fc']) for row in reader
                if row.get('is_differential', 'true') == 'true'}

def load_results(path, q_threshold=0.05):
    """Load significant features from results file."""
    significant = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                q = float(row.get('q_value', 1.0))
                if q < q_threshold:
                    est = float(row.get('estimate', 0))
                    significant[row['feature_id']] = est
            except (ValueError, KeyError):
                continue
    return significant

def compute_metrics(true_positives, significant):
    """Compute sensitivity and FDR."""
    true_set = set(true_positives.keys())
    sig_set = set(significant.keys())
    
    tp = len(true_set & sig_set)
    fp = len(sig_set - true_set)
    fn = len(true_set - sig_set)
    
    sensitivity = tp / len(true_set) if true_set else 0
    fdr = fp / len(sig_set) if sig_set else 0
    precision = tp / len(sig_set) if sig_set else 0
    
    return {
        'tp': tp,
        'fp': fp,
        'fn': fn,
        'total_true': len(true_set),
        'total_sig': len(sig_set),
        'sensitivity': sensitivity,
        'fdr': fdr,
        'precision': precision
    }

def main():
    datasets = ['typical_16s', 'sparse_virome', 'extreme_sparse', 'confounded', 'group_specific']
    models = ['linda', 'bb', 'hurdle']
    
    print(f"{'Dataset':<20} {'Model':<10} {'True':<6} {'Sig':<6} {'TP':<6} {'FP':<6} {'Sens':<8} {'FDR':<8} {'Prec':<8}")
    print("-" * 90)
    
    for dataset in datasets:
        gt_path = Path(f'benchmark_results/{dataset}/ground_truth.tsv')
        if not gt_path.exists():
            continue
        
        ground_truth = load_ground_truth(gt_path)
        
        for model in models:
            results_path = Path(f'benchmark_results/{dataset}/{model}_results.tsv')
            if not results_path.exists():
                continue
            
            significant = load_results(results_path)
            metrics = compute_metrics(ground_truth, significant)
            
            print(f"{dataset:<20} {model:<10} {metrics['total_true']:<6} {metrics['total_sig']:<6} "
                  f"{metrics['tp']:<6} {metrics['fp']:<6} "
                  f"{metrics['sensitivity']:.4f}   {metrics['fdr']:.4f}   {metrics['precision']:.4f}")
        print()

if __name__ == '__main__':
    main()
