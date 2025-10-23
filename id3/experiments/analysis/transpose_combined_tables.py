#!/usr/bin/env python3
"""

"""

import pandas as pd
import numpy as np
from pathlib import Path

def transpose_combined_table(input_file: Path, output_file: Path):

    

    df = pd.read_csv(input_file)
    

    mean_row = df[df['Protein'] == 'Mean'].iloc[0]
    df_data = df[df['Protein'] != 'Mean']
    

    df_data = df_data.set_index('Protein')
    

    transposed_data = []
    

    variants = []
    for col in df_data.columns:
        if col.endswith('_Access'):
            variant = col.replace('_Access', '')
            if variant not in variants:
                variants.append(variant)
    

    for variant in variants:
        row = {'Variant': variant}
        

        for protein in df_data.index:
            access_col = f'{variant}_Access'
            cai_col = f'{variant}_CAI'
            
            if access_col in df_data.columns and cai_col in df_data.columns:
                access_val = df_data.loc[protein, access_col]
                cai_val = df_data.loc[protein, cai_col]
                

                if pd.notna(access_val) and pd.notna(cai_val):
                    row[protein] = f"{access_val:.2f} ({cai_val:.3f})"
                else:
                    row[protein] = ""
        

        access_mean = mean_row.get(f'{variant}_Access', np.nan)
        cai_mean = mean_row.get(f'{variant}_CAI', np.nan)
        if pd.notna(access_mean) and pd.notna(cai_mean):
            row['Mean'] = f"{access_mean:.2f} ({cai_mean:.3f})"
        else:
            row['Mean'] = ""
            
        transposed_data.append(row)
    

    df_transposed = pd.DataFrame(transposed_data)
    

    proteins = ['O15263', 'P00004', 'P01308', 'P01825', 'P04637', 'P0CG48', 
                'P0DTC2', 'P0DTC9', 'P31417', 'P42212', 'P61626', 'P99999']
    column_order = ['Variant'] + proteins + ['Mean']
    df_transposed = df_transposed[column_order]
    

    variant_order = ['L00', 'L01', 'L10', 'L11', 
                     'A00', 'A01', 'A10', 'A11', 
                     'C00', 'C01', 'C10', 'C11']
    df_transposed['Variant'] = pd.Categorical(df_transposed['Variant'], 
                                              categories=variant_order, 
                                              ordered=True)
    df_transposed = df_transposed.sort_values('Variant')
    

    df_transposed.to_csv(output_file, index=False)

    
    return df_transposed

def main():
    """主函数"""
    
    tables_dir = Path('paper_experiment_results/tables')
    
    print("📊 转置组合表格...")
    

    print("\n转置With Penalty表格...")
    penalty_input = tables_dir / 'table2_combined_with_penalty.csv'
    penalty_output = tables_dir / 'table2_combined_with_penalty_transposed.csv'
    df_penalty = transpose_combined_table(penalty_input, penalty_output)
    

    print("\n转置No Penalty表格...")
    no_penalty_input = tables_dir / 'table2_combined_no_penalty.csv'
    no_penalty_output = tables_dir / 'table2_combined_no_penalty_transposed.csv'
    df_no_penalty = transpose_combined_table(no_penalty_input, no_penalty_output)
    

    print("\n📋 With Penalty表格预览（前4行）:")
    print(df_penalty.head(4).to_string(index=False))
    
    print("\n📋 No Penalty表格预览（前4行）:")
    print(df_no_penalty.head(4).to_string(index=False))
    

    print("\n📝 生成Markdown格式...")
    
    # With Penalty Markdown
    with open(tables_dir / 'table2_combined_with_penalty_transposed.md', 'w') as f:
        f.write("## Combined Access and CAI Results - With Penalty (λ=0.1, Target=0.8)\n\n")
        f.write("| Variant | " + " | ".join(df_penalty.columns[1:]) + " |\n")
        f.write("|" + "---------|" * (len(df_penalty.columns)) + "\n")
        
        for _, row in df_penalty.iterrows():
            f.write("| " + " | ".join(str(row[col]) for col in df_penalty.columns) + " |\n")
    
    print(f"✅ Markdown表格已保存: {tables_dir / 'table2_combined_with_penalty_transposed.md'}")
    
    # No Penalty Markdown
    with open(tables_dir / 'table2_combined_no_penalty_transposed.md', 'w') as f:
        f.write("## Combined Access and CAI Results - No Penalty (Baseline)\n\n")
        f.write("| Variant | " + " | ".join(df_no_penalty.columns[1:]) + " |\n")
        f.write("|" + "---------|" * (len(df_no_penalty.columns)) + "\n")
        
        for _, row in df_no_penalty.iterrows():
            f.write("| " + " | ".join(str(row[col]) for col in df_no_penalty.columns) + " |\n")
    
    print(f"✅ Markdown表格已保存: {tables_dir / 'table2_combined_no_penalty_transposed.md'}")

if __name__ == '__main__':
    main()