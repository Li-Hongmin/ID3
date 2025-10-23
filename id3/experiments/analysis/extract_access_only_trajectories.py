#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from functools import partial
import multiprocessing as mp
import time
import sys

# Add project path for imports
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

def process_single_experiment(json_file: Path, max_steps: int = 1000):

    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        

        protein = data.get('protein_name', 'unknown')
        constraint_type = data.get('constraint_type', 'unknown')
        variant = data.get('variant', 'unknown')
        seed = data.get('seed', 0)
        

        constraint_map = {
            'lagrangian': 'L',
            'ams': 'A', 
            'cpc': 'C'
        }
        constraint_code = constraint_map.get(constraint_type, 'X')
        variant_id = f"ID3-{constraint_code}{variant}"
        

        trajectory = data.get('trajectory', {})
        


        if not accessibility_values:

            accessibility_values = trajectory.get('accessibility_values', [])
            
        iterations = trajectory.get('iterations', list(range(len(accessibility_values))))
        cai_values = trajectory.get('ecai_values', [])
        discrete_cai_values = trajectory.get('discrete_cai_values', [])
        

        if len(accessibility_values) > max_steps:
            accessibility_values = accessibility_values[:max_steps]
            iterations = iterations[:max_steps]
            if cai_values:
                cai_values = cai_values[:max_steps]
            if discrete_cai_values:
                discrete_cai_values = discrete_cai_values[:max_steps]
        

        experiment_data = {
            'protein': protein,
            'constraint_type': constraint_type,
            'variant': variant,
            'variant_id': variant_id,
            'seed': seed,
            'file_name': json_file.name,
            'accessibility_values': accessibility_values,
            'iterations': iterations,
            'cai_values': cai_values,
            'discrete_cai_values': discrete_cai_values,
            'steps': len(accessibility_values)
        }
        
        return experiment_data
        
    except Exception as e:

        return None

def extract_trajectories_parallel(data_dir: str, max_workers: int = None):
    """并行提取所有实验的轨迹数据"""
    data_path = Path(data_dir)
    json_files = list(data_path.glob("*.json"))
    
    print(f"发现 {len(json_files)} 个实验文件")
    print(f"使用 {max_workers or mp.cpu_count()} 个进程并行处理")
    
    start_time = time.time()
    

    with ProcessPoolExecutor(max_workers=max_workers) as executor:

        process_func = partial(process_single_experiment, max_steps=1000)
        

        results = list(executor.map(process_func, json_files))
    

    valid_results = [r for r in results if r is not None]
    
    processing_time = time.time() - start_time
    print(f"✅ 成功处理 {len(valid_results)}/{len(json_files)} 个文件")
    print(f"⚡ 处理时间: {processing_time:.2f} 秒")
    
    return valid_results

def convert_to_detailed_dataframe(trajectory_data):



    
    if not trajectory_data:

        return pd.DataFrame()
    

    if trajectory_data:
        first_exp = trajectory_data[0]


    
    all_rows = []
    
    for exp_idx, exp_data in enumerate(trajectory_data):
        if exp_data is None:
            continue
            
        protein = exp_data.get('protein', 'unknown')
        variant_id = exp_data.get('variant_id', 'unknown') 
        seed = exp_data.get('seed', 0)
        
        accessibility_values = exp_data.get('accessibility_values', [])
        iterations = exp_data.get('iterations', [])
        cai_values = exp_data.get('cai_values', [])
        
        if not accessibility_values:

            continue
            

        if not iterations:
            iterations = list(range(len(accessibility_values)))
        

        for i, (step, acc_val) in enumerate(zip(iterations, accessibility_values)):
            row = {
                'protein': protein,
                'variant': variant_id,
                'seed': seed,
                'step': step,
                'true_accessibility': acc_val,



            }
            

            if cai_values and i < len(cai_values):
                row['cai_loss'] = cai_values[i]
            
            all_rows.append(row)
        


    
    df = pd.DataFrame(all_rows)

    
    return df

def main():
    """主函数"""

    access_data_dir = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/access_only"
    
    print("="*80)
    print("🚀 并行提取Access-only实验轨迹数据")
    print("="*80)
    

    if not Path(access_data_dir).exists():
        print(f"❌ 数据目录不存在: {access_data_dir}")
        return
    

    trajectory_data = extract_trajectories_parallel(access_data_dir)
    
    if not trajectory_data:
        print("❌ 没有成功提取到轨迹数据")
        return
    

    detailed_df = convert_to_detailed_dataframe(trajectory_data)
    

    output_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures")
    output_dir.mkdir(exist_ok=True)
    

    trajectory_csv_path = output_dir / "trajectory_access_only_1000steps.csv"
    detailed_df.to_csv(trajectory_csv_path, index=False)
    
    print(f"✅ Access-only轨迹数据已保存: {trajectory_csv_path}")
    print(f"📊 数据维度: {detailed_df.shape}")
    

    if len(detailed_df) > 0:
        print("\n📈 数据统计:")
        print(f"  变体数量: {detailed_df['variant'].nunique()}")
        print(f"  蛋白质数量: {detailed_df['protein'].nunique()}")
        print(f"  总步骤数: {detailed_df['step'].nunique()}")
        print(f"  实验种子数: {detailed_df['seed'].nunique()}")
        

        print("\n🧬 变体分布:")
        variant_counts = detailed_df['variant'].value_counts()
        for variant, count in variant_counts.items():
            print(f"  {variant}: {count:,} 数据点")
    

    trajectory_json_path = output_dir / "trajectory_access_only_data.json"
    

    json_data = []
    for exp_data in trajectory_data:
        json_exp = exp_data.copy()

        for key, value in json_exp.items():
            if isinstance(value, np.ndarray):
                json_exp[key] = value.tolist()
            elif isinstance(value, (np.float64, np.int64)):
                json_exp[key] = float(value) if 'float' in str(type(value)) else int(value)
        json_data.append(json_exp)
    
    with open(trajectory_json_path, 'w') as f:
        json.dump(json_data, f, indent=2)
    
    print(f"✅ 原始轨迹数据已保存: {trajectory_json_path}")
    print("\n🎯 数据准备完成，可以开始生成收敛分析图表!")

if __name__ == "__main__":
    main()