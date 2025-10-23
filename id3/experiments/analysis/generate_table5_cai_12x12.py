#!/usr/bin/env python3
"""


"""

import pandas as pd
import os

def parse_cell_value(cell):

    if pd.isna(cell) or cell == '':
        return None, None


    try:
        parts = str(cell).split(' (')
        if len(parts) == 2:
            accessibility = float(parts[0])
            cai = float(parts[1].rstrip(')'))
            return accessibility, cai
        else:
            return None, None
    except:
        return None, None

def generate_cai_12x12_table():
    """生成12×12 CAI表格"""


    csv_path = "../../../paper_experiment_results/tables/table2_combined_with_penalty_transposed.csv"

    if not os.path.exists(csv_path):
        print(f"错误：找不到数据文件 {csv_path}")
        return

    df = pd.read_csv(csv_path)


    latex_content = []


    latex_content.append("\\begin{table*}[!ht]")
    latex_content.append("\\centering")
    latex_content.append("\\caption{CAI optimization performance across all ID3 variants and protein targets\\label{tab:cai_12x12}}")
    latex_content.append("\\tiny")
    latex_content.append("\\setlength{\\tabcolsep}{2.5pt}")


    proteins = [col for col in df.columns if col not in ['Variant', 'Mean']]


    col_def = "l|" + "c" * len(proteins) + "|cc"
    latex_content.append(f"\\begin{{tabular}}{{{col_def}}}")

    latex_content.append("\\hline")


    header_row = "\\textbf{Variant}"
    for protein in proteins:
        header_row += f" & \\rotatebox{{90}}{{\\textbf{{{protein}}}}}"
    header_row += " & \\textbf{Avg} & \\textbf{Rank} \\\\"
    latex_content.append(header_row)

    latex_content.append("\\hline")


    variant_avg_cai = {}

    for _, row in df.iterrows():
        variant = row['Variant']
        cai_values = []

        for protein in proteins:
            cell_value = row[protein]
            _, cai = parse_cell_value(cell_value)
            if cai is not None:
                cai_values.append(cai)

        if cai_values:
            avg_cai = sum(cai_values) / len(cai_values)
            variant_avg_cai[variant] = avg_cai


    sorted_variants = sorted(variant_avg_cai.items(), key=lambda x: x[1], reverse=True)
    variant_ranks = {variant: rank+1 for rank, (variant, _) in enumerate(sorted_variants)}


    lagrangian_variants = ['L00', 'L01', 'L10', 'L11']
    ams_variants = ['A00', 'A01', 'A10', 'A11']
    cpc_variants = ['C00', 'C01', 'C10', 'C11']


    for variant_group, group_name in [(cpc_variants, "CPC"), (ams_variants, "AMS"), (lagrangian_variants, "Lagrangian")]:
        for variant in variant_group:
            if variant in df['Variant'].values:
                row_data = df[df['Variant'] == variant].iloc[0]


                data_row = f"ID3-{variant}"


                cai_values = []
                for protein in proteins:
                    cell_value = row_data[protein]
                    _, cai = parse_cell_value(cell_value)

                    if cai is not None:

                        if cai == max([parse_cell_value(df[df['Variant'] == v].iloc[0][protein])[1]
                                      for v in df['Variant'] if parse_cell_value(df[df['Variant'] == v].iloc[0][protein])[1] is not None]):
                            data_row += f" & \\textbf{{{cai:.3f}}}"
                        else:
                            data_row += f" & {cai:.3f}"
                        cai_values.append(cai)
                    else:
                        data_row += " & --"


                if cai_values:
                    avg_cai = sum(cai_values) / len(cai_values)
                    rank = variant_ranks.get(variant, len(variant_ranks))


                    if rank == 1:
                        data_row += f" & \\textbf{{{avg_cai:.3f}}} & \\textbf{{{rank}}}"
                    else:
                        data_row += f" & {avg_cai:.3f} & {rank}"
                else:
                    data_row += " & -- & --"

                data_row += " \\\\"
                latex_content.append(data_row)


        if variant_group != lagrangian_variants:
            latex_content.append("\\hline")

    latex_content.append("\\hline")
    latex_content.append("\\multicolumn{" + str(len(proteins) + 3) + "}{l}{\\textit{Note: Bold values indicate best performance in each protein. CAI values range from 0 to 1 (higher is better).}} \\\\")
    latex_content.append("\\hline")
    latex_content.append("\\end{tabular}")
    latex_content.append("\\end{table*}")


    output_path = "../../../paper_experiment_results/tables/table5_cai_12x12.tex"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(latex_content))

    print(f"✅ 生成完整12×12 CAI表格: {output_path}")
    print(f"📊 表格包含 {len(df)} 个变体 × {len(proteins)} 个蛋白质的CAI数据")


    print("\n🏆 CAI性能排名（前5位）:")
    for i, (variant, avg_cai) in enumerate(sorted_variants[:5]):
        print(f"{i+1:2d}. ID3-{variant}: {avg_cai:.3f}")

if __name__ == "__main__":
    generate_cai_12x12_table()