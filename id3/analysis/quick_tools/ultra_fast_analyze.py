#!/usr/bin/env python3
"""


"""

import json
import sys
import time
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging


logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def fast_json_extract(file_path: Path, need_sequences: bool = False) -> Optional[Dict]:
    """


    """
    try:
        with open(file_path, 'r') as f:

            content = f.read()
            

            result = {'filename': file_path.stem}
            


            import re
            

            patterns = {
                'protein_name': r'"protein_name":\s*"([^"]+)"',
                'constraint_type': r'"constraint_type":\s*"([^"]+)"',
                'variant': r'"variant":\s*"([^"]+)"',
                'best_accessibility': r'"best_accessibility":\s*([\d.]+)',
                'final_accessibility': r'"final_accessibility":\s*([\d.]+)',
                'aa_match_rate': r'"aa_match_rate":\s*([\d.]+)',
                'iterations': r'"iterations":\s*(\d+)',
                'final_ecai': r'"final_ecai":\s*([\d.]+)',
                'cai_target_achieved': r'"cai_target_achieved":\s*(true|false)',
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    value = match.group(1)

                    if key in ['best_accessibility', 'final_accessibility', 'aa_match_rate', 'final_ecai']:
                        result[key] = float(value)
                    elif key == 'iterations':
                        result[key] = int(value)
                    elif key == 'cai_target_achieved':
                        result[key] = value == 'true'
                    else:
                        result[key] = value
            

            if need_sequences:

                data = json.loads(content)
                if 'trajectory' in data and 'discrete_sequences' in data['trajectory']:
                    sequences = data['trajectory']['discrete_sequences']
                    unique_sequences = len(set(sequences))
                    result['uniqueness_rate'] = unique_sequences / len(sequences) if sequences else 0
                    result['total_sequences'] = len(sequences)
                    result['unique_sequences'] = unique_sequences
            
            return result
            
    except Exception as e:
        logger.error(f"Error loading {file_path}: {e}")
        return None

def batch_analyze(directory: Path, need_sequences: bool = False) -> List[Dict]:

    experiment_files = list(directory.glob("*.json"))

    experiment_files = [f for f in experiment_files 
                       if f.name not in ["config.json", "summary.json"]]
    
    logger.info(f"⚡ Fast analyzing {len(experiment_files)} files...")
    
    results = []
    for i, file_path in enumerate(experiment_files):
        result = fast_json_extract(file_path, need_sequences)
        if result:
            results.append(result)
        

        if (i + 1) % 30 == 0:
            logger.info(f"  Processed {i + 1}/{len(experiment_files)} files...")
    
    return results

def quick_stats(results: List[Dict]) -> None:
    """快速统计并打印结果"""
    if not results:
        print("No results to analyze")
        return
    

    best_accs = [r['best_accessibility'] for r in results if 'best_accessibility' in r]
    

    by_constraint = {}
    for r in results:
        if 'best_accessibility' in r:
            constraint = r.get('constraint_type', 'unknown')
            if constraint not in by_constraint:
                by_constraint[constraint] = []
            by_constraint[constraint].append(r['best_accessibility'])
    

    has_cai = any('final_ecai' in r for r in results)
    

    print("\n" + "=" * 80)
    print("⚡ ULTRA FAST ANALYSIS RESULTS")
    print("=" * 80)
    
    print(f"\n📊 Dataset Overview:")
    print(f"  Total Experiments: {len(results)}")
    print(f"  Experiment Type: {'CAI Optimization' if has_cai else 'Pure Accessibility Optimization'}")
    
    if best_accs:
        print(f"\n🎯 Performance Summary:")
        print(f"  Best Accessibility: {np.min(best_accs):.4f}")
        print(f"  Worst Accessibility: {np.max(best_accs):.4f}")
        print(f"  Mean Accessibility: {np.mean(best_accs):.4f} ± {np.std(best_accs):.4f}")
        
        print(f"\n📈 By Constraint Type:")
        for constraint in sorted(by_constraint.keys()):
            values = by_constraint[constraint]
            print(f"  {constraint:12s}: {np.mean(values):6.4f} ± {np.std(values):.4f} "
                  f"(best: {np.min(values):.4f}, n={len(values)})")
    

    if has_cai:
        ecai_values = [r['final_ecai'] for r in results if 'final_ecai' in r]
        achieved = sum(1 for r in results if r.get('cai_target_achieved', False))
        
        print(f"\n🧬 CAI Optimization:")
        print(f"  Mean ECAI: {np.mean(ecai_values):.4f} ± {np.std(ecai_values):.4f}")
        print(f"  Target Achievement: {achieved}/{len(results)} ({achieved/len(results)*100:.1f}%)")
    

    if best_accs:
        best_idx = np.argmin([r.get('best_accessibility', float('inf')) for r in results])
        best_exp = results[best_idx]
        
        print(f"\n🏆 Best Configuration:")
        print(f"  Protein: {best_exp.get('protein_name', 'N/A')}")
        print(f"  Constraint: {best_exp.get('constraint_type', 'N/A')}")
        print(f"  Variant: {best_exp.get('variant', 'N/A')}")
        print(f"  Best Accessibility: {best_exp.get('best_accessibility', 'N/A'):.4f}")
        if has_cai and 'final_ecai' in best_exp:
            print(f"  Final ECAI: {best_exp.get('final_ecai', 'N/A'):.4f}")

def analyze_with_sequences(directory: Path) -> None:

    results = batch_analyze(directory, need_sequences=True)
    

    quick_stats(results)
    

    uniqueness_data = [r.get('uniqueness_rate', 0) for r in results 
                       if 'uniqueness_rate' in r]
    
    if uniqueness_data:
        print(f"\n🔍 Sequence Uniqueness:")
        print(f"  Mean Uniqueness Rate: {np.mean(uniqueness_data):.1%}")
        print(f"  Std Dev: {np.std(uniqueness_data):.1%}")
        print(f"  Range: {np.min(uniqueness_data):.1%} - {np.max(uniqueness_data):.1%}")
        

        by_constraint = {}
        for r in results:
            if 'uniqueness_rate' in r:
                constraint = r.get('constraint_type', 'unknown')
                if constraint not in by_constraint:
                    by_constraint[constraint] = []
                by_constraint[constraint].append(r['uniqueness_rate'])
        
        if by_constraint:
            print(f"\n  By Constraint Type:")
            for constraint in sorted(by_constraint.keys()):
                values = by_constraint[constraint]
                print(f"    {constraint:12s}: {np.mean(values):.1%} ± {np.std(values):.1%}")

def compare_two_experiments(dir1: Path, dir2: Path) -> None:
    """快速比较两个实验"""
    print("\n" + "=" * 80)
    print("⚡ ULTRA FAST COMPARISON")
    print("=" * 80)
    

    print(f"\nLoading experiments...")
    results1 = batch_analyze(dir1, need_sequences=False)
    results2 = batch_analyze(dir2, need_sequences=False)
    

    has_cai1 = any('final_ecai' in r for r in results1)
    has_cai2 = any('final_ecai' in r for r in results2)
    
    print(f"\n📁 Experiment 1: {dir1.name}")
    print(f"   Type: {'CAI Optimization' if has_cai1 else 'Pure Accessibility'}")
    print(f"   Count: {len(results1)} experiments")
    
    print(f"\n📁 Experiment 2: {dir2.name}")
    print(f"   Type: {'CAI Optimization' if has_cai2 else 'Pure Accessibility'}")
    print(f"   Count: {len(results2)} experiments")
    

    best1 = [r['best_accessibility'] for r in results1 if 'best_accessibility' in r]
    best2 = [r['best_accessibility'] for r in results2 if 'best_accessibility' in r]
    
    if best1 and best2:
        print(f"\n📊 Performance Comparison:")
        print(f"{'Metric':<25} {'Exp 1':>12} {'Exp 2':>12} {'Difference':>12}")
        print("-" * 62)
        
        min1, min2 = np.min(best1), np.min(best2)
        print(f"{'Best Accessibility':<25} {min1:>12.4f} {min2:>12.4f} {min2-min1:>+12.4f}")
        
        mean1, mean2 = np.mean(best1), np.mean(best2)
        print(f"{'Mean Accessibility':<25} {mean1:>12.4f} {mean2:>12.4f} {mean2-mean1:>+12.4f}")
        
        std1, std2 = np.std(best1), np.std(best2)
        print(f"{'Std Dev':<25} {std1:>12.4f} {std2:>12.4f} {std2-std1:>+12.4f}")
        

        print(f"\n🏆 Winner: {'Experiment 1' if min1 < min2 else 'Experiment 2'} "
              f"(better by {abs(min2-min1):.4f})")

def main():

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python ultra_fast_analyze.py <experiment_dir>              # Quick analysis")
        print("  python ultra_fast_analyze.py <experiment_dir> --sequences  # Include uniqueness")
        print("  python ultra_fast_analyze.py <dir1> --compare <dir2>       # Compare two experiments")
        sys.exit(1)
    
    dir1 = Path(sys.argv[1])
    if not dir1.exists():
        print(f"Error: Directory {dir1} does not exist")
        sys.exit(1)
    
    start_time = time.time()
    

    if len(sys.argv) > 2:
        if sys.argv[2] == '--sequences':

            analyze_with_sequences(dir1)
        elif sys.argv[2] == '--compare' and len(sys.argv) > 3:

            dir2 = Path(sys.argv[3])
            if not dir2.exists():
                print(f"Error: Directory {dir2} does not exist")
                sys.exit(1)
            compare_two_experiments(dir1, dir2)
        else:

            results = batch_analyze(dir1, need_sequences=False)
            quick_stats(results)
    else:

        results = batch_analyze(dir1, need_sequences=False)
        quick_stats(results)
    
    elapsed_time = time.time() - start_time
    print(f"\n⏱️  Analysis completed in {elapsed_time:.1f} seconds")
    print("=" * 80)

if __name__ == "__main__":
    main()