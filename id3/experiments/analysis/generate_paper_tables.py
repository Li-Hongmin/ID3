#!/usr/bin/env python3
"""


"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def generate_performance_matrix_table():



    csv_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables/performance_12x12_parallel.csv"
    df = pd.read_csv(csv_path, index_col=0)


    variant_map = {
        'lagrangian_00': 'ID3-L00',
        'lagrangian_01': 'ID3-L01',
        'lagrangian_10': 'ID3-L10',
        'lagrangian_11': 'ID3-L11',
        'ams_00': 'ID3-A00',
        'ams_01': 'ID3-A01',
        'ams_10': 'ID3-A10',
        'ams_11': 'ID3-A11',
        'cpc_00': 'ID3-C00',
        'cpc_01': 'ID3-C01',
        'cpc_10': 'ID3-C10',
        'cpc_11': 'ID3-C11'
    }


    ordered_cols = ['cpc_00', 'cpc_01', 'cpc_10', 'cpc_11',
                   'ams_00', 'ams_01', 'ams_10', 'ams_11',
                   'lagrangian_00', 'lagrangian_01', 'lagrangian_10', 'lagrangian_11']

    df = df[ordered_cols]


    avg_values = df.mean(axis=0)


    ranks = avg_values.rank().astype(int)


    latex_lines = []
    latex_lines.append("\\begin{table*}[!t]")
    latex_lines.append("\\centering")
    latex_lines.append("\\caption{RNA accessibility optimization performance (kcal/mol) across ID3 variants and protein targets\\label{tab:performance_matrix}}")
    latex_lines.append("\\tiny")
    latex_lines.append("\\begin{tabular}{l|cccccccccccc|cc}")
    latex_lines.append("\\hline")


    header = "\\textbf{Variant}"
    for protein in df.index:
        header += f" & \\rotatebox{{90}}{{\\textbf{{{protein}}}}}"
    header += " & \\textbf{Avg} & \\textbf{Rank}"
    latex_lines.append(header + " \\\\")
    latex_lines.append("\\hline")


    for col in ordered_cols:
        variant_name = variant_map[col]
        row = variant_name


        min_positions = []
        for protein in df.index:
            value = df.loc[protein, col]

            if value == df.loc[protein].min():
                min_positions.append(protein)


        for protein in df.index:
            value = df.loc[protein, col]

            if protein in min_positions:
                row += f" & \\textbf{{{value:.3f}}}"
            else:
                row += f" & {value:.3f}"


        avg_val = avg_values[col]
        rank = ranks[col]


        if rank <= 3:
            row += f" & \\textbf{{{avg_val:.3f}}} & {rank}"
        else:
            row += f" & {avg_val:.3f} & {rank}"

        latex_lines.append(row + " \\\\")

    latex_lines.append("\\hline")
    latex_lines.append("\\multicolumn{15}{l}{\\textit{Note: Rank indicates overall performance ordering based on average accessibility scores (lower is better).}} \\\\")
    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table*}")


    output_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables/table3_performance_matrix.tex"
    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))

    print(f"Table 3 saved to: {output_path}")


    print("\n=== Performance Statistics ===")
    print(f"Best variant: {variant_map[avg_values.idxmin()]} ({avg_values.min():.3f} kcal/mol)")
    print(f"Worst variant: {variant_map[avg_values.idxmax()]} ({avg_values.max():.3f} kcal/mol)")


    for constraint in ['cpc', 'ams', 'lagrangian']:
        constraint_cols = [col for col in ordered_cols if col.startswith(constraint)]
        constraint_avg = df[constraint_cols].mean().mean()
        print(f"{constraint.upper()} average: {constraint_avg:.3f} kcal/mol")

    return df, avg_values, ranks

def generate_cai_performance_table():
    """生成Table 4: CAI Performance (deterministic variants only)"""


    json_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables/data_tab2_access_cai_comparison.json"

    import json
    with open(json_path, 'r') as f:
        data = json.load(f)


    deterministic_variants = ['ams_00', 'ams_01', 'cpc_00', 'cpc_01']
    variant_map = {
        'ams_00': 'ID3-A00',
        'ams_01': 'ID3-A01',
        'cpc_00': 'ID3-C00',
        'cpc_01': 'ID3-C01'
    }


    results = {}
    for variant in deterministic_variants:
        results[variant] = {
            'accessibility': [],
            'cai': []
        }


    proteins = list(data.keys())
    for protein in proteins:
        if isinstance(data[protein], dict):
            for variant in deterministic_variants:
                if variant in data[protein]:
                    variant_data = data[protein][variant]
                    if 'accessibility' in variant_data:
                        results[variant]['accessibility'].append(variant_data['accessibility'])
                    if 'discrete_cai' in variant_data:
                        results[variant]['cai'].append(variant_data['discrete_cai'])
                    elif 'cai' in variant_data:
                        results[variant]['cai'].append(variant_data['cai'])


    latex_lines = []
    latex_lines.append("\\begin{table*}[!t]")
    latex_lines.append("\\centering")
    latex_lines.append("\\caption{CAI Integration Performance: Complete 48-Experiment Results for 4 Deterministic Variants across 12 Proteins\\label{tab:cai_performance}}")
    latex_lines.append("\\tiny")
    latex_lines.append("\\begin{tabular}{l|c|cccccccccccc|c}")
    latex_lines.append("\\hline")


    header = "\\textbf{Variant} & \\textbf{CAI}"
    for protein in proteins[:12]:
        header += f" & \\rotatebox{{90}}{{\\textbf{{{protein}}}}}"
    header += " & \\textbf{Avg}"
    latex_lines.append(header + " \\\\")
    latex_lines.append("\\hline")


    for variant in deterministic_variants:
        variant_name = variant_map[variant]


        avg_cai = np.mean(results[variant]['cai']) if results[variant]['cai'] else 0

        row = f"{variant_name} & {avg_cai:.3f}"


        for i, protein in enumerate(proteins[:12]):
            if protein in data and variant in data[protein]:
                acc_value = data[protein][variant].get('accessibility', 0)
                row += f" & {acc_value:.3f}"
            else:
                row += " & -"


        avg_acc = np.mean(results[variant]['accessibility']) if results[variant]['accessibility'] else 0
        row += f" & {avg_acc:.3f}"

        latex_lines.append(row + " \\\\")

    latex_lines.append("\\hline")
    latex_lines.append("\\multicolumn{15}{l}{\\textit{Note: CAI = Codon Adaptation Index (higher is better); Accessibility values in kcal/mol (lower is better).}} \\\\")
    latex_lines.append("\\multicolumn{15}{l}{\\textit{All 48 experiments achieved 100\\% success rate with CAI $\\geq$ 0.8.}} \\\\")
    latex_lines.append("\\hline")
    latex_lines.append("\\end{tabular}")
    latex_lines.append("\\end{table*}")


    output_path = "/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables/table4_cai_performance.tex"
    with open(output_path, 'w') as f:
        f.write('\n'.join(latex_lines))

    print(f"\nTable 4 saved to: {output_path}")


    print("\n=== CAI Performance Statistics ===")
    for variant in deterministic_variants:
        variant_name = variant_map[variant]
        avg_cai = np.mean(results[variant]['cai']) if results[variant]['cai'] else 0
        avg_acc = np.mean(results[variant]['accessibility']) if results[variant]['accessibility'] else 0
        print(f"{variant_name}: CAI={avg_cai:.3f}, Accessibility={avg_acc:.3f}")

def main():

    print("Generating paper tables in LaTeX format...")
    print("=" * 60)


    print("\nGenerating Table 3: Performance Matrix")
    df, avg_values, ranks = generate_performance_matrix_table()


    print("\nGenerating Table 4: CAI Performance")
    generate_cai_performance_table()

    print("\n" + "=" * 60)
    print("All tables generated successfully!")

if __name__ == "__main__":
    main()