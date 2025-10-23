#!/usr/bin/env python3
"""


"""

import pandas as pd
from pathlib import Path
import json


VARIANT_MAPPING = {
    # Lagrangian variants
    'Lag_DS': 'ID3-L00',  # Lagrangian + Deterministic Soft
    'Lag_DH': 'ID3-L01',  # Lagrangian + Deterministic Hard  
    'Lag_SS': 'ID3-L10',  # Lagrangian + Stochastic Soft
    'Lag_SH': 'ID3-L11',  # Lagrangian + Stochastic Hard
    # AMS variants
    'AMS_DS': 'ID3-A00',  # Amino Matching Softmax + Deterministic Soft
    'AMS_DH': 'ID3-A01',  # Amino Matching Softmax + Deterministic Hard
    'AMS_SS': 'ID3-A10',
    'AMS_SH': 'ID3-A11',  # Amino Matching Softmax + Stochastic Hard
    # CPC variants  
    'CPC_DS': 'ID3-C00',  # Codon Profile Constraint + Deterministic Soft
    'CPC_DH': 'ID3-C01',  # Codon Profile Constraint + Deterministic Hard
    'CPC_SS': 'ID3-C10',  # Codon Profile Constraint + Stochastic Soft
    'CPC_SH': 'ID3-C11',  # Codon Profile Constraint + Stochastic Hard
}

def rename_csv_columns(csv_path: Path, output_path: Path):

    df = pd.read_csv(csv_path, index_col=0)
    

    new_columns = []
    for col in df.columns:
        renamed = col
        for old_name, new_name in VARIANT_MAPPING.items():
            renamed = renamed.replace(old_name, new_name)
        new_columns.append(renamed)
    
    df.columns = new_columns
    


    ordered_columns = []
    for prefix in ['ID3-L', 'ID3-A', 'ID3-C']:
        for suffix in ['00', '01', '10', '11']:
            variant = f"{prefix}{suffix}"

            for col in df.columns:
                if variant in col:
                    if col not in ordered_columns:
                        ordered_columns.append(col)
    

    for col in df.columns:
        if col not in ordered_columns:
            ordered_columns.append(col)
    
    df = df[ordered_columns]
    df.to_csv(output_path)

    return df

def generate_latex_table_cai_no_penalty(df_access: pd.DataFrame, df_cai: pd.DataFrame, output_path: Path):
    """生成论文格式的LaTeX表格（CAI无惩罚）"""
    latex_lines = []
    latex_lines.append(r"\begin{table}[htbp]")
    latex_lines.append(r"\centering")
    latex_lines.append(r"\caption{CAI No-Penalty Ablation: RNA Accessibility Performance without CAI Optimization\label{tab:cai_no_penalty}}")
    latex_lines.append(r"\tiny")
    latex_lines.append(r"\begin{tabular}{l|cccc|cccc|cccc|c}")
    latex_lines.append(r"\toprule")
    

    latex_lines.append(r"\multirow{2}{*}{\textbf{Protein}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Lagrangian}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Amino Matching}} & " +
                      r"\multicolumn{4}{c|}{\textbf{Codon Profile}} & " +
                      r"\multirow{2}{*}{\textbf{Avg}} \\")
    latex_lines.append(r" & L00 & L01 & L10 & L11 & A00 & A01 & A10 & A11 & C00 & C01 & C10 & C11 & \\")
    latex_lines.append(r"\midrule")
    

    for protein in df_access.index:
        row_values = []
        row_sum = 0
        row_count = 0
        
        for variant in ['ID3-L00', 'ID3-L01', 'ID3-L10', 'ID3-L11',
                       'ID3-A00', 'ID3-A01', 'ID3-A10', 'ID3-A11',
                       'ID3-C00', 'ID3-C01', 'ID3-C10', 'ID3-C11']:

            access_col = f"{variant}_Access"
            cai_col = f"{variant}_CAI"
            
            if access_col in df_access.columns:
                access_val = df_access.loc[protein, access_col]
                if pd.notna(access_val):

                    if access_val < 1.0:
                        cell = f"\\textbf{{{access_val:.3f}}}"
                    else:
                        cell = f"{access_val:.3f}"
                    row_values.append(cell)
                    row_sum += access_val
                    row_count += 1
                else:
                    row_values.append("--")
            else:
                row_values.append("--")
        

        row_avg = row_sum / row_count if row_count > 0 else 0
        row_values.append(f"{row_avg:.3f}")
        

        short_name = protein[:6] if len(protein) > 6 else protein
        latex_lines.append(f"{short_name} & " + " & ".join(row_values) + r" \\")
    

    latex_lines.append(r"\midrule")
    avg_row = ["\\textbf{Avg}"]
    for variant in ['ID3-L00', 'ID3-L01', 'ID3-L10', 'ID3-L11',
                   'ID3-A00', 'ID3-A01', 'ID3-A10', 'ID3-A11',
                   'ID3-C00', 'ID3-C01', 'ID3-C10', 'ID3-C11']:
        access_col = f"{variant}_Access"
        if access_col in df_access.columns:
            avg_val = df_access[access_col].mean()
            if pd.notna(avg_val):
                if avg_val < 1.0:
                    avg_row.append(f"\\textbf{{{avg_val:.3f}}}")
                else:
                    avg_row.append(f"{avg_val:.3f}")
            else:
                avg_row.append("--")
        else:
            avg_row.append("--")
    

    total_avg = df_access.mean().mean()
    avg_row.append(f"\\textbf{{{total_avg:.3f}}}")
    latex_lines.append(" & ".join(avg_row) + r" \\")
    

    latex_lines.append(r"\bottomrule")
    latex_lines.append(r"\end{tabular}")
    latex_lines.append(r"\end{table}")
    

    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))
    
    print(f"LaTeX表格已保存到: {output_path}")

def main():

    tables_dir = Path("paper_experiment_results/tables")
    


    csv_path = tables_dir / "table_cai_no_penalty.csv"
    if csv_path.exists():

        renamed_csv = tables_dir / "table_cai_no_penalty_paper_format.csv"
        df = rename_csv_columns(csv_path, renamed_csv)
        

        access_cols = [col for col in df.columns if '_Access' in col]
        cai_cols = [col for col in df.columns if '_CAI' in col]
        
        df_access = df[access_cols] if access_cols else df
        df_cai = df[cai_cols] if cai_cols else pd.DataFrame()
        

        latex_path = tables_dir / "table_cai_no_penalty_paper_format.tex"
        generate_latex_table_cai_no_penalty(df_access, df_cai, latex_path)
    


    csv_path = tables_dir / "data_tab1_access_comparison.csv"
    if csv_path.exists():
        renamed_csv = tables_dir / "table_access_only_paper_format.csv"
        df = rename_csv_columns(csv_path, renamed_csv)
    



    


if __name__ == "__main__":
    main()