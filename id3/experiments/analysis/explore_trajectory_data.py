#!/usr/bin/env python3
"""

"""

import json
from pathlib import Path

def explore_trajectory_structure():


    sample_file = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/merged_cai_experiments/20250911_165648_P42212_lagrangian_11_seed42.json")
    
    print("="*60)

    print("="*60)
    
    with open(sample_file, 'r') as f:
        data = json.load(f)
    


    






    

    trajectory = data.get('trajectory', {})


    
    if isinstance(trajectory, dict):


        keys = list(trajectory.keys())[:10]

        

        if keys:
            first_key = keys[0]
            first_point = trajectory[first_key]


            if isinstance(first_point, dict):

                for key, value in first_point.items():
                    if isinstance(value, (list, tuple)):
                        print(f"    {key}: {type(value).__name__} with {len(value)} items")
                    else:
                        print(f"    {key}: {value}")
            


            for key in keys:
                value = trajectory[key]
                if isinstance(value, list):

                else:
                    print(f"  {key}: {type(value)} = {value}")
                    

            if 'accessibility_values' in trajectory:
                acc_values = trajectory['accessibility_values']
                iterations = trajectory.get('iterations', list(range(len(acc_values))))



                if len(acc_values) >= 10:





                    
    elif isinstance(trajectory, list):

        if len(trajectory) > 0:

            if isinstance(trajectory[0], dict):

    


    for key in data.keys():
        if 'trajectory' in key.lower() or 'history' in key.lower() or 'step' in key.lower():
            value = data[key]
            if isinstance(value, (list, dict)):
                print(f"  {key}: {type(value).__name__} with {len(value)} items")
            else:
                print(f"  {key}: {value}")

if __name__ == "__main__":
    explore_trajectory_structure()