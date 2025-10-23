#!/usr/bin/env python3
"""


"""

import json
import os
import sys
from pathlib import Path
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

def quick_analyze_method(batch_dir: Path, method: str, max_samples: int = 5) -> dict:

    method_files = list(batch_dir.glob(f"*{method}*.json"))[:max_samples]
    
    accessibility_values = []
    final_accessibility_values = []
    
    for file_path in method_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            best_acc = data.get('best_accessibility')
            final_acc = data.get('final_accessibility')
            
            if best_acc is not None:
                accessibility_values.append(best_acc)
            if final_acc is not None:
                final_accessibility_values.append(final_acc)
                
        except Exception as e:

            continue
    
    result = {
        'total_files': len(list(batch_dir.glob(f"*{method}*.json"))),
        'analyzed_files': len(method_files),
        'best_accessibility': {
            'values': accessibility_values,
            'min': min(accessibility_values) if accessibility_values else None,
            'mean': np.mean(accessibility_values) if accessibility_values else None
        },
        'final_accessibility': {
            'values': final_accessibility_values,
            'min': min(final_accessibility_values) if final_accessibility_values else None,
            'mean': np.mean(final_accessibility_values) if final_accessibility_values else None
        }
    }
    
    return result

def main():

    

    print("=" * 80)
    
    batches = [



    ]
    
    all_results = {}
    
    for batch_name, batch_desc in batches:
        batch_dir = reliable_dir / batch_name
        
        print(f"\n📊 {batch_desc} ({batch_name})")
        print("-" * 60)
        
        batch_results = {}
        
        for method in ['lagrangian', 'ams', 'cpc']:
            result = quick_analyze_method(batch_dir, method, max_samples=5)
            batch_results[method] = result
            
            best_min = result['best_accessibility']['min']
            best_mean = result['best_accessibility']['mean']
            final_min = result['final_accessibility']['min']
            

            if best_min and best_min < 1.0:

            elif best_min and best_min < 1.5:

            elif best_min and best_min < 2.0:

            else:

            

            if best_min:


            if final_min:

        
        all_results[batch_name] = batch_results
    


    print("=" * 80)

    print("-" * 90)
    
    for batch_name, batch_desc in batches:
        short_name = batch_name.split('_')[0]  # 20250909
        results = all_results[batch_name]
        
        for method in ['lagrangian', 'ams', 'cpc']:
            result = results[method]
            best_min = result['best_accessibility']['min']
            best_mean = result['best_accessibility']['mean']
            final_min = result['final_accessibility']['min']
            

            if best_min and best_min < 1.0:

            elif best_min and best_min < 1.5:

            elif best_min and best_min < 2.0:

            else:

            
            best_min_str = f"{best_min:.3f}" if best_min else "N/A"
            best_mean_str = f"{best_mean:.3f}" if best_mean else "N/A"
            final_min_str = f"{final_min:.3f}" if final_min else "N/A"
            
            print(f"{short_name:<15} {method.upper():<12} {result['total_files']:<8} {best_min_str:<12} {best_mean_str:<12} {final_min_str:<12} {level:<12}")
    


    print("=" * 80)
    
    access1_results = all_results["20250910_004126_unified_access_experiments"]
    access2_results = all_results["20250910_022355_unified_access_experiments"]
    
    for method in ['lagrangian', 'ams', 'cpc']:
        mean1 = access1_results[method]['best_accessibility']['mean']
        mean2 = access2_results[method]['best_accessibility']['mean']
        
        if mean1 and mean2:
            diff = abs(mean1 - mean2)
            relative_diff = diff / mean1 * 100
            

            

        else:

    

    print("=" * 80)






if __name__ == "__main__":
    main()