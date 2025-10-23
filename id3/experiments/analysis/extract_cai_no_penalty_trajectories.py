#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
import os
from typing import Dict, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

def load_experiment_file(file_path: str) -> Dict:

    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

def extract_trajectory_from_experiment(file_path: str) -> pd.DataFrame:
    """从单个实验文件提取轨迹数据"""
    data = load_experiment_file(file_path)
    if not data:
        return None


    protein = data.get('protein_name', 'Unknown')
    constraint_type = data.get('constraint_type', 'Unknown')
    variant_str = data.get('variant', 'Unknown')
    seed = data.get('seed', 42)


    constraint_map = {
        'lagrangian': 'L',
        'ams': 'A',
        'cpc': 'C'
    }
    variant_prefix = constraint_map.get(constraint_type.lower(), 'U')
    variant = f"ID3-{variant_prefix}{variant_str}"


    trajectory = []


    if 'trajectory' in data:
        traj_data = data['trajectory']


        access_values = traj_data.get('accessibility', [])
        if not access_values:
            access_values = traj_data.get('accessibility_values', [])

        loss_values = traj_data.get('unified_loss', [])
        if not loss_values:
            loss_values = traj_data.get('loss_values', [])

        cai_losses = traj_data.get('cai_loss', [0.0] * len(access_values))
        discrete_cai_values = traj_data.get('discrete_cai_values', [])


        for step in range(min(1000, len(access_values))):
            trajectory.append({
                'protein': protein,
                'variant': variant,
                'seed': seed,
                'step': step,
                'true_accessibility': access_values[step] if step < len(access_values) else np.nan,
                'accessibility_loss': loss_values[step] if step < len(loss_values) else (access_values[step] if step < len(access_values) else np.nan),
                'constraint_loss': 0.0,
                'temperature': 1.0,
                'cai_loss': cai_losses[step] if step < len(cai_losses) else 0.0,
                'discrete_cai': discrete_cai_values[step] if step < len(discrete_cai_values) else np.nan
            })

    if trajectory:
        return pd.DataFrame(trajectory)
    return None

def process_file(file_path: Path) -> pd.DataFrame:

    return extract_trajectory_from_experiment(str(file_path))

def extract_all_trajectories(data_dir: str, output_dir: str = None):
    """提取所有实验的轨迹数据"""

    data_path = Path(data_dir)
    if not data_path.exists():
        print(f"数据目录不存在: {data_dir}")
        return


    json_files = list(data_path.glob("*.json"))
    print(f"找到 {len(json_files)} 个实验文件")

    if not json_files:
        print("没有找到任何JSON文件")
        return


    all_trajectories = []

    with ProcessPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(process_file, f): f for f in json_files}

        for future in tqdm(as_completed(futures), total=len(futures), desc="处理实验文件"):
            try:
                df = future.result()
                if df is not None and not df.empty:
                    all_trajectories.append(df)
            except Exception as e:
                file_path = futures[future]
                print(f"处理文件 {file_path} 时出错: {e}")

    if not all_trajectories:
        print("没有成功提取任何轨迹数据")
        return


    combined_df = pd.concat(all_trajectories, ignore_index=True)
    print(f"\n成功提取 {len(combined_df)} 行轨迹数据")
    print(f"覆盖 {combined_df['protein'].nunique()} 个蛋白质")
    print(f"覆盖 {combined_df['variant'].nunique()} 个变体")


    if output_dir is None:
        output_dir = "paper_experiment_results/figures"

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)


    csv_path = output_path / "trajectory_cai_no_penalty_1000steps.csv"
    combined_df.to_csv(csv_path, index=False)
    print(f"轨迹数据已保存到: {csv_path}")


    json_path = output_path / "trajectory_cai_no_penalty_data.json"
    combined_df.to_json(json_path, orient='records', indent=2)
    print(f"JSON数据已保存到: {json_path}")


    print("\n=== 统计信息 ===")
    print(f"总步数: {combined_df['step'].max() + 1}")
    print(f"平均最终可及性: {combined_df[combined_df['step'] == combined_df['step'].max()]['true_accessibility'].mean():.3f}")


    print("\n各变体最终性能:")
    final_step = combined_df['step'].max()
    final_data = combined_df[combined_df['step'] == final_step]

    for variant in sorted(final_data['variant'].unique()):
        variant_data = final_data[final_data['variant'] == variant]
        mean_access = variant_data['true_accessibility'].mean()
        std_access = variant_data['true_accessibility'].std()
        print(f"  {variant}: {mean_access:.3f} ± {std_access:.3f}")

    return combined_df

if __name__ == "__main__":

    cai_no_penalty_dir = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/cai_no_penalty"

    print("开始提取CAI无惩罚实验轨迹...")
    print("=" * 50)

    df = extract_all_trajectories(cai_no_penalty_dir)

    print("\n提取完成！")