#!/usr/bin/env python3
"""

"""

import json
import numpy as np

def compare_accessibility_arrays(file_path):

    with open(file_path, 'r') as f:
        data = json.load(f)
    

    trajectory = data.get('trajectory', {})
    accessibility1 = trajectory.get('accessibility', [])
    accessibility_values = data.get('accessibility_values', [])
    



    
    if len(accessibility1) == 0:

        return
    
    if len(accessibility_values) == 0:

        return
    

    min_len = min(len(accessibility1), len(accessibility_values))
    


    print("-" * 65)
    
    differences = []
    for i in range(min(10, min_len)):
        val1 = accessibility1[i]
        val2 = accessibility_values[i]
        diff = abs(val1 - val2)
        differences.append(diff)
        print(f"{i:4d} | {val1:20.10f} | {val2:17.10f} | {diff:.2e}")
    

    if len(accessibility1) == len(accessibility_values):
        all_same = all(abs(a - b) < 1e-10 for a, b in zip(accessibility1, accessibility_values))
        if all_same:

        else:
            max_diff = max(abs(a - b) for a, b in zip(accessibility1, accessibility_values))

    else:

    

    min1 = min(accessibility1)
    min2 = min(accessibility_values)




    

    best_accessibility = data.get('best_accessibility')
    if best_accessibility:
        print(f"\nbest_accessibility: {best_accessibility:.10f}")



def main():
    file_path = "results_full/20250909_105804_unified_access_experiments/20250909_105838_O15263_lagrangian_01_seed42.json"
    

    print("=" * 80)
    
    try:
        compare_accessibility_arrays(file_path)
    except Exception as e:


if __name__ == "__main__":
    main()