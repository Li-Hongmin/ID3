#!/usr/bin/env python3
"""


"""

import json
import sys
import time
import numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import multiprocessing as mp
from typing import Dict, List, Any
import logging


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)

def load_single_experiment(file_path: Path) -> Dict[str, Any]:

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        

        result = {
            'protein': data.get('protein_name'),
            'constraint': data.get('constraint_type'),
            'variant': data.get('variant'),
            'best_accessibility': data.get('best_accessibility'),
            'final_accessibility': data.get('final_accessibility'),
            'aa_match_rate': data.get('aa_match_rate'),
            'iterations': data.get('iterations'),
            'filename': file_path.stem
        }
        

        if 'trajectory' in data and 'discrete_sequences' in data['trajectory']:
            sequences = data['trajectory']['discrete_sequences']
            unique_sequences = len(set(sequences))
            result['uniqueness_rate'] = unique_sequences / len(sequences) if sequences else 0
            result['total_sequences'] = len(sequences)
            result['unique_sequences'] = unique_sequences
        

        if 'final_ecai' in data:
            result['final_ecai'] = data['final_ecai']
            result['cai_target_achieved'] = data.get('cai_target_achieved', False)
        
        return result
    except Exception as e:
        logger.error(f"Error loading {file_path}: {e}")
        return None

def parallel_load_experiments(directory: Path, max_workers: int = None) -> List[Dict]:
    """并行加载所有实验文件"""
    if max_workers is None:
        max_workers = min(mp.cpu_count(), 8)
    
    experiment_files = list(directory.glob("*.json"))

    experiment_files = [f for f in experiment_files 
                       if f.name not in ["config.json", "summary.json"]]
    
    logger.info(f"Loading {len(experiment_files)} files using {max_workers} workers...")
    
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        future_to_file = {executor.submit(load_single_experiment, f): f 
                         for f in experiment_files}
        

        completed = 0
        for future in as_completed(future_to_file):
            result = future.result()
            if result:
                results.append(result)
            completed += 1
            if completed % 20 == 0:
                logger.info(f"  Loaded {completed}/{len(experiment_files)} files...")
    
    return results

def analyze_performance(results: List[Dict]) -> Dict[str, Any]:

    if not results:
        return {}
    

    best_accs = [r['best_accessibility'] for r in results if 'best_accessibility' in r]
    final_accs = [r['final_accessibility'] for r in results if 'final_accessibility' in r]
    

    by_constraint = defaultdict(list)
    for r in results:
        if 'best_accessibility' in r:
            by_constraint[r['constraint']].append(r['best_accessibility'])
    

    constraint_stats = {}
    for constraint, values in by_constraint.items():
        constraint_stats[constraint] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'count': len(values)
        }
    
    return {
        'total_experiments': len(results),
        'best_overall': np.min(best_accs) if best_accs else None,
        'worst_overall': np.max(best_accs) if best_accs else None,
        'mean_best': np.mean(best_accs) if best_accs else None,
        'mean_final': np.mean(final_accs) if final_accs else None,
        'by_constraint': constraint_stats
    }

def analyze_uniqueness(results: List[Dict]) -> Dict[str, Any]:
    """快速唯一率分析"""
    uniqueness_data = [r.get('uniqueness_rate', 0) for r in results 
                       if 'uniqueness_rate' in r]
    
    if not uniqueness_data:
        return {'error': 'No uniqueness data found'}
    

    by_constraint = defaultdict(list)
    for r in results:
        if 'uniqueness_rate' in r:
            by_constraint[r['constraint']].append(r['uniqueness_rate'])
    

    constraint_stats = {}
    for constraint, values in by_constraint.items():
        constraint_stats[constraint] = {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values)
        }
    
    return {
        'global_uniqueness_rate': np.mean(uniqueness_data),
        'uniqueness_std': np.std(uniqueness_data),
        'min_uniqueness': np.min(uniqueness_data),
        'max_uniqueness': np.max(uniqueness_data),
        'by_constraint': constraint_stats
    }

def analyze_cai(results: List[Dict]) -> Dict[str, Any]:

    cai_data = [r for r in results if 'final_ecai' in r]
    
    if not cai_data:
        return {'has_cai': False, 'message': 'No CAI data found (non-CAI experiment)'}
    
    ecai_values = [r['final_ecai'] for r in cai_data]
    achieved = sum(1 for r in cai_data if r.get('cai_target_achieved', False))
    
    return {
        'has_cai': True,
        'total_experiments': len(cai_data),
        'mean_ecai': np.mean(ecai_values),
        'std_ecai': np.std(ecai_values),
        'target_achieved': achieved,
        'achievement_rate': achieved / len(cai_data) if cai_data else 0
    }

def print_results(perf: Dict, unique: Dict, cai: Dict, elapsed_time: float):
    """打印分析结果"""
    print("\n" + "=" * 80)
    print("📊 FAST EXPERIMENT ANALYSIS RESULTS")
    print("=" * 80)
    

    print("\n🎯 Performance Analysis:")
    print(f"  Total Experiments: {perf.get('total_experiments', 0)}")
    print(f"  Best Accessibility: {perf.get('best_overall', 'N/A'):.4f}")
    print(f"  Mean Best: {perf.get('mean_best', 'N/A'):.4f}")
    
    print("\n  By Constraint Type:")
    for constraint, stats in perf.get('by_constraint', {}).items():
        print(f"    {constraint:12s}: {stats['mean']:.4f} ± {stats['std']:.4f} "
              f"(min: {stats['min']:.4f}, max: {stats['max']:.4f})")
    

    print("\n🔍 Uniqueness Analysis:")
    if 'error' not in unique:
        print(f"  Global Uniqueness Rate: {unique.get('global_uniqueness_rate', 0):.1%}")
        print(f"  Std Dev: {unique.get('uniqueness_std', 0):.1%}")
        print(f"  Range: {unique.get('min_uniqueness', 0):.1%} - {unique.get('max_uniqueness', 0):.1%}")
        
        print("\n  By Constraint Type:")
        for constraint, stats in unique.get('by_constraint', {}).items():
            print(f"    {constraint:12s}: {stats['mean']:.1%} ± {stats['std']:.1%}")
    else:
        print(f"  {unique['error']}")
    

    print("\n🧬 CAI Analysis:")
    if cai.get('has_cai'):
        print(f"  Mean ECAI: {cai.get('mean_ecai', 0):.4f}")
        print(f"  Target Achievement: {cai.get('achievement_rate', 0):.1%}")
    else:
        print(f"  {cai.get('message', 'No CAI data')}")
    
    print(f"\n⏱️  Analysis completed in {elapsed_time:.1f} seconds")
    print("=" * 80)

def compare_experiments(dir1: Path, dir2: Path):

    print("\n" + "=" * 80)
    print("📊 COMPARATIVE ANALYSIS")
    print("=" * 80)
    

    print(f"\n1️⃣  Loading {dir1.name}...")
    results1 = parallel_load_experiments(dir1)
    perf1 = analyze_performance(results1)
    unique1 = analyze_uniqueness(results1)
    cai1 = analyze_cai(results1)
    
    print(f"\n2️⃣  Loading {dir2.name}...")
    results2 = parallel_load_experiments(dir2)
    perf2 = analyze_performance(results2)
    unique2 = analyze_uniqueness(results2)
    cai2 = analyze_cai(results2)
    

    print("\n" + "-" * 80)
    print("COMPARISON RESULTS")
    print("-" * 80)
    
    print(f"\n{'Metric':<30} {'Exp 1':<20} {'Exp 2':<20} {'Difference'}")
    print("-" * 80)
    

    if perf1.get('best_overall') and perf2.get('best_overall'):
        diff = perf2['best_overall'] - perf1['best_overall']
        print(f"{'Best Accessibility':<30} {perf1['best_overall']:<20.4f} {perf2['best_overall']:<20.4f} {diff:+.4f}")
    
    if perf1.get('mean_best') and perf2.get('mean_best'):
        diff = perf2['mean_best'] - perf1['mean_best']
        print(f"{'Mean Best Accessibility':<30} {perf1['mean_best']:<20.4f} {perf2['mean_best']:<20.4f} {diff:+.4f}")
    

    if unique1.get('global_uniqueness_rate') and unique2.get('global_uniqueness_rate'):
        diff = unique2['global_uniqueness_rate'] - unique1['global_uniqueness_rate']
        print(f"{'Uniqueness Rate':<30} {unique1['global_uniqueness_rate']:<20.1%} {unique2['global_uniqueness_rate']:<20.1%} {diff:+.1%}")
    

    print(f"{'Has CAI':<30} {str(cai1.get('has_cai', False)):<20} {str(cai2.get('has_cai', False)):<20}")
    
    print("\n" + "=" * 80)

def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("Usage: python fast_analyze.py <experiment_dir> [--compare <experiment_dir2>]")
        sys.exit(1)
    
    dir1 = Path(sys.argv[1])
    if not dir1.exists():
        print(f"Error: Directory {dir1} does not exist")
        sys.exit(1)
    
    start_time = time.time()
    

    if len(sys.argv) > 2 and sys.argv[2] == '--compare' and len(sys.argv) > 3:
        dir2 = Path(sys.argv[3])
        if not dir2.exists():
            print(f"Error: Directory {dir2} does not exist")
            sys.exit(1)
        compare_experiments(dir1, dir2)
    else:

        logger.info(f"Analyzing {dir1}...")
        

        results = parallel_load_experiments(dir1)
        

        perf = analyze_performance(results)
        unique = analyze_uniqueness(results)
        cai = analyze_cai(results)
        

        elapsed_time = time.time() - start_time
        print_results(perf, unique, cai, elapsed_time)
        

        output_file = dir1 / "fast_analysis_results.json"
        with open(output_file, 'w') as f:
            json.dump({
                'performance': perf,
                'uniqueness': unique,
                'cai': cai,
                'analysis_time': elapsed_time
            }, f, indent=2)
        
        print(f"\n💾 Results saved to: {output_file}")

if __name__ == "__main__":
    main()