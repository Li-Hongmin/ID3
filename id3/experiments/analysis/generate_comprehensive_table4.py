#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path

def load_cai_data():


    base_path = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables")
    

    with_penalty = pd.read_csv(base_path / "table2_combined_with_penalty_transposed.csv")
    no_penalty = pd.read_csv(base_path / "table2_combined_no_penalty_transposed.csv")
    
    return with_penalty, no_penalty

def parse_cell_value(cell):
    """解析单元格值，提取可及性和CAI"""
    if pd.isna(cell) or cell == '':
        return None, None
    
    try:

        parts = str(cell).split(' (')
        if len(parts) == 2:
            accessibility = float(parts[0])
            cai = float(parts[1].rstrip(')'))
            return accessibility, cai
    except:
        pass
    
    return None, None

def generate_comprehensive_table():

    

    with_penalty, no_penalty = load_cai_data()
    

    latex_lines = []
    

    latex_lines.append(r"\begin{table*}[!ht]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{Comprehensive CAI optimization performance across all ID3 variants\label{tab:cai_comprehensive}}")
    latex_lines.append(r"\small")
    latex_lines.append(r"\setlength{\tabcolsep}{3pt}")
    latex_lines.append(r"\begin{tabular}{l|cc|cc|cc|cc}")
    latex_lines.append(r"\hline")
    

    latex_lines.append(r"\multirow{3}{*}{\textbf{Variant}} & \multicolumn{4}{c|}{\textbf{With Constraint Penalty}} & \multicolumn{4}{c}{\textbf{Without Constraint Penalty}} \\\\")
    latex_lines.append(r"\cline{2-9}")
    latex_lines.append(r" & \multicolumn{2}{c|}{\textbf{Accessibility}} & \multicolumn{2}{c|}{\textbf{CAI}} & \multicolumn{2}{c|}{\textbf{Accessibility}} & \multicolumn{2}{c}{\textbf{CAI}} \\\\")
    latex_lines.append(r"\cline{2-9}")
    latex_lines.append(r" & Mean & Std & Mean & Std & Mean & Std & Mean & Std \\\\")
    latex_lines.append(r"\hline")
    

    variant_map = {
        'L00': 'ID3-L00', 'L01': 'ID3-L01', 'L10': 'ID3-L10', 'L11': 'ID3-L11',
        'A00': 'ID3-A00', 'A01': 'ID3-A01', 'A10': 'ID3-A10', 'A11': 'ID3-A11',
        'C00': 'ID3-C00', 'C01': 'ID3-C01', 'C10': 'ID3-C10', 'C11': 'ID3-C11'
    }
    

    proteins = ['O15263', 'P00004', 'P01308', 'P01825', 'P04637', 'P0CG48', 
                'P0DTC2', 'P0DTC9', 'P31417', 'P42212', 'P61626', 'P99999']
    

    all_data = {}


    for _, row_with in with_penalty.iterrows():
        variant = row_with['Variant']
        if variant not in variant_map:
            continue


        row_without = no_penalty[no_penalty['Variant'] == variant].iloc[0] if not no_penalty[no_penalty['Variant'] == variant].empty else None


        access_with = []
        cai_with = []
        for protein in proteins:
            if protein in row_with:
                acc, cai = parse_cell_value(row_with[protein])
                if acc is not None:
                    access_with.append(acc)
                    cai_with.append(cai)


        access_without = []
        cai_without = []
        if row_without is not None:
            for protein in proteins:
                if protein in row_without:
                    acc, cai = parse_cell_value(row_without[protein])
                    if acc is not None:
                        access_without.append(acc)
                        cai_without.append(cai)


        all_data[variant] = {
            'access_with': access_with,
            'cai_with': cai_with,
            'access_without': access_without,
            'cai_without': cai_without
        }


    best_access_with = float('inf')
    best_access_without = float('inf')
    best_cai_with = 0
    best_cai_without = 0
    best_variant_access_with = ''
    best_variant_access_without = ''
    best_variant_cai_with = ''
    best_variant_cai_without = ''

    for variant, data in all_data.items():
        if data['access_with']:
            mean_access_with = np.mean(data['access_with'])
            mean_cai_with = np.mean(data['cai_with'])

            if mean_access_with < best_access_with:
                best_access_with = mean_access_with
                best_variant_access_with = variant
            if mean_cai_with > best_cai_with:
                best_cai_with = mean_cai_with
                best_variant_cai_with = variant

        if data['access_without']:
            mean_access_without = np.mean(data['access_without'])
            mean_cai_without = np.mean(data['cai_without'])

            if mean_access_without < best_access_without:
                best_access_without = mean_access_without
                best_variant_access_without = variant
            if mean_cai_without > best_cai_without:
                best_cai_without = mean_cai_without
                best_variant_cai_without = variant


    variant_order = ['L00', 'L01', 'L10', 'L11', 'A00', 'A01', 'A10', 'A11', 'C00', 'C01', 'C10', 'C11']

    for variant in variant_order:
        if variant not in all_data:
            continue

        data = all_data[variant]
        line = f"{variant_map[variant]} & "


        if data['access_with']:
            mean_access_with = np.mean(data['access_with'])
            std_access_with = np.std(data['access_with'])
            mean_cai_with = np.mean(data['cai_with'])
            std_cai_with = np.std(data['cai_with'])


            if variant == best_variant_access_with:
                line += f"\\textbf{{{mean_access_with:.3f}}} & {std_access_with:.3f} & "
            else:
                line += f"{mean_access_with:.3f} & {std_access_with:.3f} & "

            if variant == best_variant_cai_with:
                line += f"\\textbf{{{mean_cai_with:.3f}}} & {std_cai_with:.3f} & "
            else:
                line += f"{mean_cai_with:.3f} & {std_cai_with:.3f} & "
        else:
            line += "- & - & - & - & "


        if data['access_without']:
            mean_access_without = np.mean(data['access_without'])
            std_access_without = np.std(data['access_without'])
            mean_cai_without = np.mean(data['cai_without'])
            std_cai_without = np.std(data['cai_without'])


            if variant == best_variant_access_without:
                line += f"\\textbf{{{mean_access_without:.3f}}} & {std_access_without:.3f} & "
            else:
                line += f"{mean_access_without:.3f} & {std_access_without:.3f} & "

            if variant == best_variant_cai_without:
                line += f"\\textbf{{{mean_cai_without:.3f}}} & {std_cai_without:.3f}"
            else:
                line += f"{mean_cai_without:.3f} & {std_cai_without:.3f}"
        else:
            line += "- & - & - & -"

        line += " \\\\"
        latex_lines.append(line)


        if variant in ['L11', 'A11']:
            latex_lines.append(r"\hline")
    

    latex_lines.append(r"\hline")
    latex_lines.append(r"\multicolumn{9}{l}{\\textit{Note: Bold values indicate best performance in each category. Accessibility in kcal/mol (lower is better), CAI (higher is better).}} \\\\")
    latex_lines.append(r"\hline")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{table*}")
    

    output_path = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables")
    output_path.mkdir(parents=True, exist_ok=True)
    
    with open(output_path / "table4_cai_comprehensive.tex", 'w') as f:
        f.write('\n'.join(latex_lines))
    
    print(f"Table 4 (comprehensive) saved to {output_path}/table4_cai_comprehensive.tex")
    

    print("\n=== CAI Performance Summary ===")
    print(f"Best Accessibility (with penalty): {best_variant_access_with} = {best_access_with:.3f} kcal/mol")
    print(f"Best Accessibility (without penalty): {best_variant_access_without} = {best_access_without:.3f} kcal/mol")
    print(f"Best CAI (with penalty): {best_variant_cai_with} = {best_cai_with:.3f}")
    print(f"Best CAI (without penalty): {best_variant_cai_without} = {best_cai_without:.3f}")
    

    generate_detailed_protein_table(with_penalty, no_penalty)

def generate_detailed_protein_table(with_penalty, no_penalty):
    """生成包含所有蛋白质详细数据的表格"""
    
    latex_lines = []
    

    latex_lines.append(r"\begin{table*}[!ht]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{Detailed CAI optimization results for representative proteins\label{tab:cai_detailed}}")
    latex_lines.append(r"\tiny")
    latex_lines.append(r"\setlength{\tabcolsep}{2pt}")
    

    selected_proteins = ['P00004', 'P01308', 'P04637', 'P0DTC2', 'P42212', 'P99999']
    
    latex_lines.append(r"\begin{tabular}{l|" + 'cc|' * len(selected_proteins) + "}")
    latex_lines.append(r"\hline")
    

    header = r"\textbf{Variant} & "
    for protein in selected_proteins:
        header += f"\\multicolumn{{2}}{{c|}}{{\\textbf{{{protein}}}}} & "
    header = header.rstrip(' & ') + " \\\\"
    latex_lines.append(header)
    
    latex_lines.append(r"\cline{2-" + str(len(selected_proteins) * 2 + 1) + "}")
    
    subheader = r" & "
    for _ in selected_proteins:
        subheader += "Access & CAI & "
    subheader = subheader.rstrip(' & ') + " \\\\"
    latex_lines.append(subheader)
    latex_lines.append(r"\hline")
    

    for _, row in with_penalty.iterrows():
        variant = row['Variant']
        if pd.isna(variant):
            continue
        
        line = f"ID3-{variant} & "
        
        for protein in selected_proteins:
            if protein in row:
                acc, cai = parse_cell_value(row[protein])
                if acc is not None:
                    line += f"{acc:.2f} & {cai:.3f} & "
                else:
                    line += "- & - & "
            else:
                line += "- & - & "
        
        line = line.rstrip(' & ') + " \\\\"
        latex_lines.append(line)
        

        if variant in ['L11', 'A11']:
            latex_lines.append(r"\hline")
    

    latex_lines.append(r"\hline")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{table*}")
    

    output_path = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables")
    with open(output_path / "table4_cai_detailed.tex", 'w') as f:
        f.write('\n'.join(latex_lines))
    
    print(f"Detailed table saved to {output_path}/table4_cai_detailed.tex")

if __name__ == "__main__":
    generate_comprehensive_table()