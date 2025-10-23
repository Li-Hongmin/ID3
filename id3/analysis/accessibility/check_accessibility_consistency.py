#!/usr/bin/env python3
"""


"""

import json
import numpy as np
from pathlib import Path

def check_accessibility_consistency(file_path):

    with open(file_path, 'r') as f:
        data = json.load(f)
    

    best_accessibility = data.get('best_accessibility')
    final_accessibility = data.get('final_accessibility')

    trajectory = data.get('trajectory', {})
    accessibility_values = trajectory.get('accessibility', [])
    
    if not accessibility_values:
        return {
            'file': file_path.name,
            'error': 'No accessibility_values found'
        }
    

    min_accessibility = min(accessibility_values)
    max_accessibility = max(accessibility_values)
    final_value = accessibility_values[-1] if accessibility_values else None
    

    best_matches_min = abs(best_accessibility - min_accessibility) < 0.0001
    final_matches_last = abs(final_accessibility - final_value) < 0.0001
    

    min_positions = [i for i, v in enumerate(accessibility_values) if abs(v - min_accessibility) < 0.0001]
    
    return {
        'file': file_path.name,
        'best_accessibility': best_accessibility,
        'final_accessibility': final_accessibility,
        'min_accessibility_values': min_accessibility,
        'final_accessibility_values': final_value,
        'total_iterations': len(accessibility_values),
        
        'consistency_checks': {
            'best_matches_min': best_matches_min,
            'final_matches_last': final_matches_last,
            'best_vs_min_diff': abs(best_accessibility - min_accessibility),
            'final_vs_last_diff': abs(final_accessibility - final_value) if final_value else None
        },
        
        'min_positions': min_positions,
        'min_position_percentage': min_positions[0] / len(accessibility_values) * 100 if min_positions else None
    }

def main():

    import glob
    pattern = "results_full/20250909_105804_unified_access_experiments/*lagrangian*.json"

    

    print("=" * 80)
    
    for file_path_str in test_files:
        file_path = Path(file_path_str)
        if not file_path.exists():

            continue
        
        try:
            result = check_accessibility_consistency(file_path)
            

            print(f"   Best accessibility: {result['best_accessibility']:.6f}")
            print(f"   Min accessibility_values: {result['min_accessibility_values']:.6f}")
            print(f"   Final accessibility: {result['final_accessibility']:.6f}")
            print(f"   Final accessibility_values: {result['final_accessibility_values']:.6f}")
            
            consistency = result['consistency_checks']
            print(f"   ✅ Best matches min: {consistency['best_matches_min']} (diff: {consistency['best_vs_min_diff']:.8f})")
            print(f"   ✅ Final matches last: {consistency['final_matches_last']} (diff: {consistency['final_vs_last_diff']:.8f})")
            
            if result['min_position_percentage']:
                print(f"   📊 Min value at {result['min_position_percentage']:.1f}% of iterations")
            
        except Exception as e:


if __name__ == "__main__":
    main()