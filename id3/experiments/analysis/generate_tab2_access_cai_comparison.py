#!/usr/bin/env python3
"""
Tab2: Access+CAI Comparison Table Generator











"""

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import sys
import time

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent.parent.parent))

from id3.cai.unified_calculator import UnifiedCAICalculator


CONSTRAINT_FULL_NAMES = {
    'lagrangian': 'Lagrangian Multiplier',
    'ams': 'Amino Matching Softmax',
    'cpc': 'Codon Profile Constraint'
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

def validate_cai_calculation(rna_sequence: str, expected_cai: float, 
                            species: str = 'ecoli_bl21de3', tolerance: float = 0.01) -> bool:
    """验证CAI计算是否正确"""
    try:

        calculator = UnifiedCAICalculator(species=species)
        

        calculated_cai = calculator.compute_cai(rna_sequence, method='standard')
        

        return abs(calculated_cai - expected_cai) < tolerance
    except Exception:
        return False

def process_single_file(file_path: Path) -> Optional[Dict[str, Any]]:

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
        

        if accessibility is None or discrete_cai is None:
            return None
        

        expected_amino = data.get('expected_amino_acids', '')
        constraint_satisfied = validate_sequence_constraint(discrete_sequence, expected_amino)
        

        cai_valid = validate_cai_calculation(discrete_sequence, discrete_cai)
        
        return {
            'constraint_type': constraint_type,
            'variant': variant,
            'protein_name': protein_name,
            'accessibility': accessibility,
            'discrete_cai': discrete_cai,
            'constraint_satisfied': constraint_satisfied,
            'cai_valid': cai_valid,
            'file_path': str(file_path)
        }
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return None

def process_files_parallel(directory: Path, num_processes: int = 40) -> List[Dict[str, Any]]:
    """并行处理所有实验文件"""

    json_files = list(directory.glob("*.json"))
    

    json_files = [f for f in json_files if f.name not in ['config.json', 'progress.json', 'summary.json']]
    
    print(f"找到 {len(json_files)} 个实验文件")
    
    results = []
    

    with ProcessPoolExecutor(max_workers=num_processes) as executor:

        future_to_file = {executor.submit(process_single_file, file_path): file_path 
                         for file_path in json_files}
        

        with tqdm(total=len(json_files), desc="处理文件") as pbar:
            for future in as_completed(future_to_file):
                result = future.result()
                if result:
                    results.append(result)
                pbar.update(1)
    
    return results

def organize_results_by_protein(results: List[Dict[str, Any]]) -> Dict[str, Dict[Tuple[str, str], Dict[str, Any]]]:


    protein_matrix = {}
    
    for result in results:
        protein = result['protein_name']
        key = (result['constraint_type'], result['variant'])
        
        if protein not in protein_matrix:
            protein_matrix[protein] = {}
        
        protein_matrix[protein][key] = result
    
    return protein_matrix

def generate_full_latex_table(protein_matrix: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]], 
                             output_path: Path) -> None:
    """生成完整的144个实验结果的LaTeX表格"""
    

    constraint_order = ['lagrangian', 'ams', 'cpc']
    variant_order = ['00', '01', '10', '11']
    

    proteins = sorted(protein_matrix.keys())
    

    latex_lines = []
    latex_lines.append(r"\begin{table}[htbp]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{Comprehensive CAI Integration Performance: All 144 Individual Experiments}")
    latex_lines.append(r"\label{tab:access_cai_comparison_full}")
    latex_lines.append(r"\scriptsize")
    latex_lines.append(r"\begin{adjustbox}{width=\textwidth}")
    latex_lines.append(r"\begin{tabular}{l|cccc|cccc|cccc}")
    latex_lines.append(r"\toprule")
    

    latex_lines.append(r"\multirow{2}{*}{\textbf{Protein}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Lagrangian Multiplier}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Amino Matching Softmax}} & " +
                      r"\multicolumn{4}{c}{\textbf{Codon Profile Constraint}} \\")
    
    header_row = r" & " + " & ".join(["DS", "DH", "SS", "SH"] * 3) + r" \\"
    latex_lines.append(header_row)
    latex_lines.append(r"\midrule")
    

    for protein in proteins:
        if protein not in protein_matrix or protein == 'unknown':
            continue
            
        protein_data = protein_matrix[protein]
        

        access_row = f"\\multirow{{2}}{{*}}{{{protein}}} & "
        access_values = []
        
        for constraint in constraint_order:
            for variant in variant_order:
                key = (constraint, variant)
                if key in protein_data:
                    result = protein_data[key]
                    access_val = result['accessibility']
                    

                    if not result['constraint_satisfied']:
                        access_values.append(f"\\cellcolor{{red!20}}{access_val:.3f}")
                    else:
                        access_values.append(f"{access_val:.3f}")
                else:
                    access_values.append("--")
        
        access_row += " & ".join(access_values) + r" \\"
        latex_lines.append(access_row)
        

        cai_row = r" & "
        cai_values = []
        
        for constraint in constraint_order:
            for variant in variant_order:
                key = (constraint, variant)
                if key in protein_data:
                    result = protein_data[key]
                    cai_val = result['discrete_cai']
                    

                    if not result['cai_valid']:
                        cai_values.append(f"\\cellcolor{{yellow!20}}({cai_val:.3f})")
                    else:
                        cai_values.append(f"({cai_val:.3f})")
                else:
                    cai_values.append("--")
        
        cai_row += " & ".join(cai_values) + r" \\"
        latex_lines.append(cai_row)
        

        if protein != proteins[-1]:
            latex_lines.append(r"\hdashline")
    

    latex_lines.append(r"\bottomrule")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{adjustbox}")
    latex_lines.append(r"\begin{tablenotes}")
    latex_lines.append(r"\scriptsize")
    latex_lines.append(r"\item Each cell shows Accessibility (kcal/mol) on top and (CAI) below in parentheses")
    latex_lines.append(r"\item DS: Deterministic Soft, DH: Deterministic Hard, SS: Stochastic Soft, SH: Stochastic Hard")
    latex_lines.append(r"\item Red cells indicate constraint violation; Yellow cells indicate CAI validation failure")
    latex_lines.append(r"\item All experiments use CAI target=0.8, $\lambda_{CAI}$=0.1, E. coli BL21(DE3) codon table")
    latex_lines.append(r"\end{tablenotes}")
    latex_lines.append(r"\end{table}")
    

    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))
    
    print(f"完整LaTeX表格已保存到: {output_path}")

def generate_simplified_latex_table(protein_matrix: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]], 
                                   output_path: Path) -> None:

    

    constraint_order = ['lagrangian', 'ams', 'cpc']
    variant_order = ['00', '01', '10', '11']
    

    proteins = sorted([p for p in protein_matrix.keys() if p != 'unknown'])
    

    latex_lines = []
    latex_lines.append(r"\begin{table}[htbp]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{CAI Integration Performance: Complete 144-Experiment Results}")
    latex_lines.append(r"\label{tab:access_cai_comparison}")

    latex_lines.append(r"\begin{adjustbox}{width=\textwidth}")
    latex_lines.append(r"\begin{tabular}{l|cccc|cccc|cccc}")
    latex_lines.append(r"\toprule")
    

    latex_lines.append(r"\multirow{2}{*}{\textbf{Protein}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Lagrangian Multiplier}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Amino Matching Softmax}} & " +
                      r"\multicolumn{4}{c}{\textbf{Codon Profile Constraint}} \\")
    
    header_row = r" & " + " & ".join(["DS", "DH", "SS", "SH"] * 3) + r" \\"
    latex_lines.append(header_row)
    latex_lines.append(r"\midrule")
    

    for protein in proteins:
        protein_data = protein_matrix[protein]
        

        row = f"{protein} & "
        cell_values = []
        
        for constraint in constraint_order:
            for variant in variant_order:
                key = (constraint, variant)
                if key in protein_data:
                    result = protein_data[key]
                    access_val = result['accessibility']
                    cai_val = result['discrete_cai']
                    

                    cell_content = f"\\makecell{{{access_val:.2f}\\\\\\footnotesize({cai_val:.3f})}}"
                    cell_values.append(cell_content)
                else:
                    cell_values.append("--")
        
        row += " & ".join(cell_values) + r" \\"
        latex_lines.append(row)
    

    latex_lines.append(r"\bottomrule")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{adjustbox}")
    latex_lines.append(r"\begin{tablenotes}")
    latex_lines.append(r"\scriptsize")
    latex_lines.append(r"\item Each cell shows Accessibility (kcal/mol) and CAI value in parentheses")
    latex_lines.append(r"\item DS: Deterministic Soft, DH: Deterministic Hard, SS: Stochastic Soft, SH: Stochastic Hard")
    latex_lines.append(r"\item All experiments use CAI target=0.8, $\lambda_{CAI}$=0.1, E. coli BL21(DE3) codon table")
    latex_lines.append(r"\end{tablenotes}")
    latex_lines.append(r"\end{table}")
    

    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))
    


def generate_summary_stats(protein_matrix: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]]) -> None:
    """生成汇总统计信息"""
    print("\n=== 完整实验统计 ===\n")
    

    total_proteins = len([p for p in protein_matrix.keys() if p != 'unknown'])
    total_experiments = sum(len(data) for data in protein_matrix.values())
    
    print(f"蛋白质数量: {total_proteins}")
    print(f"总实验数: {total_experiments}")
    

    print("\n每个蛋白质的实验完整性:")
    for protein in sorted(protein_matrix.keys()):
        if protein == 'unknown':
            continue
        count = len(protein_matrix[protein])
        print(f"  {protein}: {count}/12 实验")
    

    best_access = float('inf')
    best_config = None
    
    for protein, data in protein_matrix.items():
        if protein == 'unknown':
            continue
        for key, result in data.items():
            if result['accessibility'] < best_access:
                best_access = result['accessibility']
                best_config = (protein, key[0], key[1])
    
    if best_config:
        print(f"\n最佳Accessibility性能:")
        print(f"  蛋白质: {best_config[0]}")
        print(f"  配置: {CONSTRAINT_FULL_NAMES[best_config[1]]} - {VARIANT_NAMES[best_config[2]]}")
        print(f"  值: {best_access:.3f} kcal/mol")

def main(experiment_dir: str, output_dir: str):

    start_time = time.time()
    

    exp_dir = Path(experiment_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    


    

    results = process_files_parallel(exp_dir, num_processes=40)
    

    

    protein_matrix = organize_results_by_protein(results)
    

    generate_summary_stats(protein_matrix)
    

    simplified_path = out_dir / "tab2_access_cai_comparison.tex"
    generate_simplified_latex_table(protein_matrix, simplified_path)
    

    full_path = out_dir / "tab2_access_cai_comparison_full.tex"
    generate_full_latex_table(protein_matrix, full_path)
    

    json_path = out_dir / "data_tab2_access_cai_comparison.json"
    with open(json_path, 'w') as f:

        json_data = {}
        for protein, data in protein_matrix.items():
            json_data[protein] = {f"{k[0]}_{k[1]}": v for k, v in data.items()}
        json.dump(json_data, f, indent=2)

    
    elapsed_time = time.time() - start_time


if __name__ == "__main__":
    import argparse
    

    parser.add_argument('--exp-dir', type=str,
                       default='paper_experiment_results/cai_with_penalty',

    parser.add_argument('--output-dir', type=str,
                       default='paper_experiment_results/tables',

    
    args = parser.parse_args()
    main(args.exp_dir, args.output_dir)