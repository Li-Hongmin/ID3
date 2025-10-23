#!/usr/bin/env python3
"""
CAI No Penalty Performance Table Generator







"""

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import sys
import time

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from id3.cai.unified_calculator import UnifiedCAICalculator


CONSTRAINT_NAMES = {
    'lagrangian': 'Lag',
    'ams': 'AMS',
    'cpc': 'CPC'
}


VARIANT_NAMES = {
    '00': 'DS',  # Deterministic Soft
    '01': 'DH',  # Deterministic Hard  
    '10': 'SS',  # Stochastic Soft
    '11': 'SH'   # Stochastic Hard
}

def validate_sequence_constraint(rna_sequence: str, amino_sequence: str) -> bool:

    if not rna_sequence or not amino_sequence:
        return False
    
    try:

        rna_seq = Seq(rna_sequence.replace('T', 'U'))
        

        translated = str(rna_seq.translate())
        

        if translated.endswith('*'):
            translated = translated[:-1]
        
        return translated == amino_sequence
    except Exception:
        return False

def process_single_file(file_path: Path) -> Optional[Dict[str, Any]]:
    """处理单个实验文件"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        

        constraint_type = data.get('constraint_type', '')
        variant = data.get('variant', '')
        protein_name = data.get('protein_name', '')
        

        best_design = data.get('best_seq_design', {})
        if not best_design:
            return None
        
        accessibility = best_design.get('accessibility', None)
        discrete_cai = best_design.get('discrete_cai', None)
        discrete_sequence = best_design.get('discrete_sequence', '')
        

        if accessibility is None:
            return None
        

        expected_amino = data.get('expected_amino_acids', '')
        constraint_satisfied = validate_sequence_constraint(discrete_sequence, expected_amino)
        

        final_ecai = data.get('final_ecai', discrete_cai)
        
        return {
            'constraint_type': constraint_type,
            'variant': variant,
            'protein_name': protein_name,
            'accessibility': accessibility,
            'discrete_cai': discrete_cai if discrete_cai is not None else final_ecai,
            'final_ecai': final_ecai,
            'constraint_satisfied': constraint_satisfied,
            'file_path': str(file_path)
        }
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

def process_files_parallel(directory: Path, num_processes: int = 20) -> List[Dict[str, Any]]:


    json_files = list(directory.glob("*.json"))
    

    json_files = [f for f in json_files if f.name not in ['config.json', 'progress.json', 'summary.json']]
    

    
    results = []
    

    with ProcessPoolExecutor(max_workers=num_processes) as executor:

        future_to_file = {executor.submit(process_single_file, file_path): file_path 
                         for file_path in json_files}
        


            for future in as_completed(future_to_file):
                result = future.result()
                if result:
                    results.append(result)
                pbar.update(1)
    
    return results

def organize_results_matrix(results: List[Dict[str, Any]]) -> pd.DataFrame:
    """组织结果为12x12矩阵格式的DataFrame"""

    proteins = sorted(set(r['protein_name'] for r in results if r['protein_name'] != 'unknown'))
    

    columns = []
    for constraint in ['lagrangian', 'ams', 'cpc']:
        for variant in ['00', '01', '10', '11']:
            columns.append(f"{CONSTRAINT_NAMES[constraint]}_{VARIANT_NAMES[variant]}")
    

    df_access = pd.DataFrame(index=proteins, columns=columns)
    df_cai = pd.DataFrame(index=proteins, columns=columns)
    

    for result in results:
        if result['protein_name'] == 'unknown':
            continue
        protein = result['protein_name']
        col_name = f"{CONSTRAINT_NAMES[result['constraint_type']]}_{VARIANT_NAMES[result['variant']]}"
        
        df_access.loc[protein, col_name] = result['accessibility']
        df_cai.loc[protein, col_name] = result['discrete_cai']
    
    return df_access, df_cai

def generate_latex_table(df_access: pd.DataFrame, df_cai: pd.DataFrame, output_path: Path) -> None:

    latex_lines = []
    latex_lines.append(r"\begin{table}[htbp]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{CAI No Penalty Ablation Study: Accessibility Performance Without CAI Penalty}")
    latex_lines.append(r"\label{tab:cai_no_penalty}")
    latex_lines.append(r"\scriptsize")
    latex_lines.append(r"\begin{adjustbox}{width=\textwidth}")
    latex_lines.append(r"\begin{tabular}{l|cccc|cccc|cccc}")
    latex_lines.append(r"\toprule")
    

    latex_lines.append(r"\multirow{2}{*}{\textbf{Protein}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Lagrangian}} & " +
                      r"\multicolumn{4}{c|}{\textbf{AMS}} & " +
                      r"\multicolumn{4}{c}{\textbf{CPC}} \\")
    
    header_row = r" & " + " & ".join(["DS", "DH", "SS", "SH"] * 3) + r" \\"
    latex_lines.append(header_row)
    latex_lines.append(r"\midrule")
    

    for protein in df_access.index:
        row_values = []
        for col in df_access.columns:
            access_val = df_access.loc[protein, col]
            cai_val = df_cai.loc[protein, col]
            
            if pd.notna(access_val):

                if access_val < 1.0:
                    cell = f"\\textbf{{{access_val:.3f}}}"
                else:
                    cell = f"{access_val:.3f}"
                

                if pd.notna(cai_val):
                    cell = f"\\makecell{{{cell}\\\\\\scriptsize({cai_val:.3f})}}"
                    
                row_values.append(cell)
            else:
                row_values.append("--")
        
        latex_lines.append(f"{protein} & " + " & ".join(row_values) + r" \\")
    

    latex_lines.append(r"\midrule")
    

    avg_row = ["\\textbf{Avg}"]
    for col in df_access.columns:
        avg_access = df_access[col].mean()
        avg_cai = df_cai[col].mean()
        if pd.notna(avg_access):
            avg_row.append(f"\\makecell{{{avg_access:.3f}\\\\\\scriptsize({avg_cai:.3f})}}")
        else:
            avg_row.append("--")
    latex_lines.append(" & ".join(avg_row) + r" \\")
    

    min_row = ["\\textbf{Min}"]
    for col in df_access.columns:
        min_access = df_access[col].min()
        if pd.notna(min_access):
            min_row.append(f"{min_access:.3f}")
        else:
            min_row.append("--")
    latex_lines.append(" & ".join(min_row) + r" \\")
    

    latex_lines.append(r"\bottomrule")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{adjustbox}")
    latex_lines.append(r"\begin{tablenotes}")
    latex_lines.append(r"\scriptsize")
    latex_lines.append(r"\item Each cell shows Accessibility (kcal/mol) with CAI value in parentheses")
    latex_lines.append(r"\item Bold values indicate excellent performance (<1.0 kcal/mol)")
    latex_lines.append(r"\item CAI penalty disabled ($\lambda_{CAI}=0$), CAI values shown for reference only")
    latex_lines.append(r"\item DS: Deterministic Soft, DH: Deterministic Hard, SS: Stochastic Soft, SH: Stochastic Hard")
    latex_lines.append(r"\end{tablenotes}")
    latex_lines.append(r"\end{table}")
    

    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))
    


def generate_csv_table(df_access: pd.DataFrame, df_cai: pd.DataFrame, output_path: Path) -> None:
    """生成CSV格式的表格"""

    combined_df = pd.DataFrame()
    
    for col in df_access.columns:
        combined_df[f"{col}_Access"] = df_access[col]
        combined_df[f"{col}_CAI"] = df_cai[col]
    

    combined_df.to_csv(output_path)
    print(f"CSV表格已保存到: {output_path}")

def generate_summary_stats(results: List[Dict[str, Any]]) -> Dict[str, Any]:

    stats = {
        'total_experiments': len(results),
        'constraint_satisfied': sum(1 for r in results if r['constraint_satisfied']),
        'avg_accessibility': np.mean([r['accessibility'] for r in results]),
        'min_accessibility': min(r['accessibility'] for r in results),
        'max_accessibility': max(r['accessibility'] for r in results),
        'avg_cai': np.mean([r['discrete_cai'] for r in results if r['discrete_cai'] is not None]),
        'excellent_count': sum(1 for r in results if r['accessibility'] < 1.0),
        'good_count': sum(1 for r in results if 1.0 <= r['accessibility'] < 1.5),
    }
    

    stats['by_constraint'] = {}
    for constraint in ['lagrangian', 'ams', 'cpc']:
        constraint_results = [r for r in results if r['constraint_type'] == constraint]
        if constraint_results:
            stats['by_constraint'][constraint] = {
                'count': len(constraint_results),
                'avg': np.mean([r['accessibility'] for r in constraint_results]),
                'min': min(r['accessibility'] for r in constraint_results),
                'max': max(r['accessibility'] for r in constraint_results),
            }
    
    return stats

def main(experiment_dir: str, output_dir: str):
    """主函数"""
    start_time = time.time()
    

    exp_dir = Path(experiment_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"从 {exp_dir} 读取实验数据")
    print(f"输出到 {out_dir}")
    

    results = process_files_parallel(exp_dir, num_processes=20)
    
    print(f"\n成功处理 {len(results)} 个实验文件")
    

    stats = generate_summary_stats(results)
    
    print("\n=== CAI无惩罚实验统计 ===")
    print(f"总实验数: {stats['total_experiments']}")
    print(f"约束满足: {stats['constraint_satisfied']}/{stats['total_experiments']} " +
          f"({stats['constraint_satisfied']/stats['total_experiments']*100:.1f}%)")
    print(f"\nAccessibility性能:")
    print(f"  平均值: {stats['avg_accessibility']:.3f} kcal/mol")
    print(f"  最小值: {stats['min_accessibility']:.3f} kcal/mol")
    print(f"  最大值: {stats['max_accessibility']:.3f} kcal/mol")
    print(f"  优秀(<1.0): {stats['excellent_count']} ({stats['excellent_count']/stats['total_experiments']*100:.1f}%)")
    print(f"  良好(1.0-1.5): {stats['good_count']} ({stats['good_count']/stats['total_experiments']*100:.1f}%)")
    print(f"\n平均CAI值: {stats['avg_cai']:.3f}")
    
    print("\n按约束类型统计:")
    for constraint, data in stats['by_constraint'].items():
        print(f"  {CONSTRAINT_NAMES[constraint]}:")
        print(f"    平均: {data['avg']:.3f} kcal/mol")
        print(f"    最佳: {data['min']:.3f} kcal/mol")
    

    df_access, df_cai = organize_results_matrix(results)
    

    latex_path = out_dir / "table_cai_no_penalty.tex"
    generate_latex_table(df_access, df_cai, latex_path)
    

    csv_path = out_dir / "table_cai_no_penalty.csv"
    generate_csv_table(df_access, df_cai, csv_path)
    

    json_path = out_dir / "data_cai_no_penalty.json"
    with open(json_path, 'w') as f:
        json.dump({
            'results': results,
            'stats': stats,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S')
        }, f, indent=2, default=str)
    print(f"原始数据已保存到: {json_path}")
    
    elapsed_time = time.time() - start_time
    print(f"\n总耗时: {elapsed_time:.2f} 秒")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="生成CAI无惩罚实验性能表格")
    parser.add_argument('--exp-dir', type=str, 
                       default='paper_experiment_results/ablation_cai_no_penalty/20250912_013943_unified_cai_experiments',
                       help='实验结果目录')
    parser.add_argument('--output-dir', type=str,
                       default='paper_experiment_results/tables',
                       help='输出目录')
    
    args = parser.parse_args()
    main(args.exp_dir, args.output_dir)