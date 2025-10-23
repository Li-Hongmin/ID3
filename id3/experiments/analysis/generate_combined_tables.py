#!/usr/bin/env python3
"""

"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def process_single_file(file_path: Path) -> Optional[Dict]:

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            

        result = {
            'file_name': file_path.name,

            'constraint_type': data.get('constraint_type'),

            'best_accessibility': data.get('best_accessibility', np.nan),
        }
        

        if 'best_seq_design' in data and isinstance(data['best_seq_design'], dict):
            result['discrete_cai'] = data['best_seq_design'].get('discrete_cai', np.nan)
        else:
            result['discrete_cai'] = np.nan
            
        return result
    except Exception as e:
        return None

def load_experiments_parallel(directory: Path, num_workers: int = 40) -> Dict[str, Dict]:
    """并行加载实验数据"""
    json_files = [f for f in directory.glob("*.json") 
                  if not f.name.startswith(('config', 'progress', 'summary'))]
    
    print(f"   找到 {len(json_files)} 个实验文件")
    
    results = {}
    with ProcessPoolExecutor(max_workers=num_workers) as executor:

        future_to_file = {
            executor.submit(process_single_file, f): f 
            for f in json_files
        }
        

        with tqdm(total=len(json_files), desc="   处理文件", leave=False) as pbar:
            for future in as_completed(future_to_file):
                data = future.result()
                if data:

                    key = f"{data['protein_id']}_{data['constraint_type']}_{data['discretization_variant']}"
                    results[key] = data
                pbar.update(1)
    
    print(f"   成功处理 {len(results)}/{len(json_files)} 个文件")
    return results

def generate_combined_table(penalty_dir: Path, no_penalty_dir: Path, output_dir: Path):

    
    proteins = ['O15263', 'P00004', 'P01308', 'P01825', 'P04637', 'P0CG48', 
                'P0DTC2', 'P0DTC9', 'P31417', 'P42212', 'P61626', 'P99999']
    

    variant_map = {
        ('lagrangian', '00'): 'L00', ('lagrangian', '01'): 'L01',
        ('lagrangian', '10'): 'L10', ('lagrangian', '11'): 'L11',
        ('ams', '00'): 'A00', ('ams', '01'): 'A01',
        ('ams', '10'): 'A10', ('ams', '11'): 'A11',
        ('cpc', '00'): 'C00', ('cpc', '01'): 'C01',
        ('cpc', '10'): 'C10', ('cpc', '11'): 'C11',
    }
    


    penalty_data = load_experiments_parallel(penalty_dir)
    

    no_penalty_data = load_experiments_parallel(no_penalty_dir)
    


    penalty_rows = []
    for protein in proteins:
        row = {'Protein': protein}
        for (method, variant), short_name in variant_map.items():
            key = f"{protein}_{method}_{variant}"
            if key in penalty_data:
                data = penalty_data[key]
                row[f'{short_name}_Access'] = data['best_accessibility']
                row[f'{short_name}_CAI'] = data['discrete_cai']
            else:
                row[f'{short_name}_Access'] = np.nan
                row[f'{short_name}_CAI'] = np.nan
        
        penalty_rows.append(row)
    

    mean_row = {'Protein': 'Mean'}
    for (method, variant), short_name in variant_map.items():
        access_vals = [row[f'{short_name}_Access'] for row in penalty_rows if not pd.isna(row[f'{short_name}_Access'])]
        cai_vals = [row[f'{short_name}_CAI'] for row in penalty_rows if not pd.isna(row[f'{short_name}_CAI'])]
        mean_row[f'{short_name}_Access'] = np.mean(access_vals) if access_vals else np.nan
        mean_row[f'{short_name}_CAI'] = np.mean(cai_vals) if cai_vals else np.nan
    penalty_rows.append(mean_row)
    
    df_penalty = pd.DataFrame(penalty_rows)
    df_penalty.to_csv(output_dir / 'table2_combined_with_penalty.csv', index=False)

    


    no_penalty_rows = []
    for protein in proteins:
        row = {'Protein': protein}
        for (method, variant), short_name in variant_map.items():
            key = f"{protein}_{method}_{variant}"
            if key in no_penalty_data:
                data = no_penalty_data[key]
                row[f'{short_name}_Access'] = data['best_accessibility']
                row[f'{short_name}_CAI'] = data['discrete_cai']
            else:
                row[f'{short_name}_Access'] = np.nan
                row[f'{short_name}_CAI'] = np.nan
        
        no_penalty_rows.append(row)
    

    mean_row = {'Protein': 'Mean'}
    for (method, variant), short_name in variant_map.items():
        access_vals = [row[f'{short_name}_Access'] for row in no_penalty_rows if not pd.isna(row[f'{short_name}_Access'])]
        cai_vals = [row[f'{short_name}_CAI'] for row in no_penalty_rows if not pd.isna(row[f'{short_name}_CAI'])]
        mean_row[f'{short_name}_Access'] = np.mean(access_vals) if access_vals else np.nan
        mean_row[f'{short_name}_CAI'] = np.mean(cai_vals) if cai_vals else np.nan
    no_penalty_rows.append(mean_row)
    
    df_no_penalty = pd.DataFrame(no_penalty_rows)
    df_no_penalty.to_csv(output_dir / 'table2_combined_no_penalty.csv', index=False)

    
    return df_penalty, df_no_penalty

if __name__ == '__main__':

    penalty_dir = Path('paper_experiment_results/cai_with_penalty')
    no_penalty_dir = Path('paper_experiment_results/cai_no_penalty')
    output_dir = Path('paper_experiment_results/tables')
    output_dir.mkdir(parents=True, exist_ok=True)
    


    df_penalty, df_no_penalty = generate_combined_table(penalty_dir, no_penalty_dir, output_dir)
    

