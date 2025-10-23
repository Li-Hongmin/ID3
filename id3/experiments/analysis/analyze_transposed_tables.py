#!/usr/bin/env python3
"""

"""

import pandas as pd
import numpy as np
from pathlib import Path

def analyze_combined_table(csv_file: Path):

    

    df = pd.read_csv(csv_file)
    

    results = []
    
    for _, row in df.iterrows():
        variant = row['Variant']
        

        access_values = []

            if pd.notna(row[col]) and row[col]:

                try:
                    parts = str(row[col]).split(' (')
                    if parts:
                        access_val = float(parts[0])
                        access_values.append(access_val)
                except:
                    continue
        

        if access_values:
            mean_access = np.mean(access_values)
            std_access = np.std(access_values)
            min_access = np.min(access_values)
            max_access = np.max(access_values)
        else:
            mean_access = np.nan
            std_access = np.nan
            min_access = np.nan
            max_access = np.nan
        
        results.append({
            'Variant': variant,
            'Mean_Accessibility': mean_access,
            'Std': std_access,
            'Min': min_access,
            'Max': max_access,
            'Count': len(access_values)
        })
    

    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values('Mean_Accessibility')
    df_results['Rank'] = range(1, len(df_results) + 1)
    
    return df_results

def main():
    """主函数"""
    
    tables_dir = Path('paper_experiment_results/tables')
    
    print("📊 分析转置后的组合表格...")
    print("=" * 80)
    

    print("\n### With Penalty (λ=0.1, Target=0.8) - Accessibility Analysis")
    print("-" * 60)
    
    penalty_file = tables_dir / 'table2_combined_with_penalty_transposed.csv'
    df_penalty = analyze_combined_table(penalty_file)
    

    print("\n**按Accessibility性能排名（越低越好）:**\n")
    print("| Rank | Variant | Mean Access | Std | Min | Max | Count |")
    print("|------|---------|-------------|-----|-----|-----|-------|")
    
    for _, row in df_penalty.iterrows():
        print(f"| {row['Rank']:2d} | {row['Variant']:7s} | {row['Mean_Accessibility']:.4f} | "
              f"{row['Std']:.4f} | {row['Min']:.2f} | {row['Max']:.2f} | {row['Count']:2d} |")
    

    print("\n**按约束机制分组统计:**\n")
    

    lag_variants = ['L00', 'L01', 'L10', 'L11']
    lag_means = df_penalty[df_penalty['Variant'].isin(lag_variants)]['Mean_Accessibility'].values
    print(f"• **Lagrangian (L)**: {np.mean(lag_means):.4f} ± {np.std(lag_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_penalty[df_penalty['Variant'].isin(lag_variants)].iloc[0]['Variant']} ({df_penalty[df_penalty['Variant'].isin(lag_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    

    ams_variants = ['A00', 'A01', 'A10', 'A11']
    ams_means = df_penalty[df_penalty['Variant'].isin(ams_variants)]['Mean_Accessibility'].values
    print(f"• **AMS (A)**: {np.mean(ams_means):.4f} ± {np.std(ams_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_penalty[df_penalty['Variant'].isin(ams_variants)].iloc[0]['Variant']} ({df_penalty[df_penalty['Variant'].isin(ams_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    

    cpc_variants = ['C00', 'C01', 'C10', 'C11']
    cpc_means = df_penalty[df_penalty['Variant'].isin(cpc_variants)]['Mean_Accessibility'].values
    print(f"• **CPC (C)**: {np.mean(cpc_means):.4f} ± {np.std(cpc_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_penalty[df_penalty['Variant'].isin(cpc_variants)].iloc[0]['Variant']} ({df_penalty[df_penalty['Variant'].isin(cpc_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    
    print("\n" + "=" * 80)
    

    print("\n### No Penalty (Baseline) - Accessibility Analysis")
    print("-" * 60)
    
    no_penalty_file = tables_dir / 'table2_combined_no_penalty_transposed.csv'
    df_no_penalty = analyze_combined_table(no_penalty_file)
    

    print("\n**按Accessibility性能排名（越低越好）:**\n")
    print("| Rank | Variant | Mean Access | Std | Min | Max | Count |")
    print("|------|---------|-------------|-----|-----|-----|-------|")
    
    for _, row in df_no_penalty.iterrows():
        print(f"| {row['Rank']:2d} | {row['Variant']:7s} | {row['Mean_Accessibility']:.4f} | "
              f"{row['Std']:.4f} | {row['Min']:.2f} | {row['Max']:.2f} | {row['Count']:2d} |")
    

    print("\n**按约束机制分组统计:**\n")
    

    lag_means = df_no_penalty[df_no_penalty['Variant'].isin(lag_variants)]['Mean_Accessibility'].values
    print(f"• **Lagrangian (L)**: {np.mean(lag_means):.4f} ± {np.std(lag_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_no_penalty[df_no_penalty['Variant'].isin(lag_variants)].iloc[0]['Variant']} ({df_no_penalty[df_no_penalty['Variant'].isin(lag_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    

    ams_means = df_no_penalty[df_no_penalty['Variant'].isin(ams_variants)]['Mean_Accessibility'].values
    print(f"• **AMS (A)**: {np.mean(ams_means):.4f} ± {np.std(ams_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_no_penalty[df_no_penalty['Variant'].isin(ams_variants)].iloc[0]['Variant']} ({df_no_penalty[df_no_penalty['Variant'].isin(ams_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    

    cpc_means = df_no_penalty[df_no_penalty['Variant'].isin(cpc_variants)]['Mean_Accessibility'].values
    print(f"• **CPC (C)**: {np.mean(cpc_means):.4f} ± {np.std(cpc_means):.4f} kcal/mol")
    print(f"  - 最佳变体: {df_no_penalty[df_no_penalty['Variant'].isin(cpc_variants)].iloc[0]['Variant']} ({df_no_penalty[df_no_penalty['Variant'].isin(cpc_variants)].iloc[0]['Mean_Accessibility']:.4f})")
    
    print("\n" + "=" * 80)
    

    print("\n
    print("-" * 60)
    

    df_comparison = pd.merge(
        df_penalty[['Variant', 'Mean_Accessibility', 'Rank']],
        df_no_penalty[['Variant', 'Mean_Accessibility', 'Rank']],
        on='Variant',
        suffixes=('_Penalty', '_NoPenalty')
    )
    
    df_comparison['Diff'] = df_comparison['Mean_Accessibility_Penalty'] - df_comparison['Mean_Accessibility_NoPenalty']
    df_comparison['Rank_Change'] = df_comparison['Rank_NoPenalty'] - df_comparison['Rank_Penalty']
    
    print("\n**性能变化（Penalty - No Penalty）:**\n")
    print("| Variant | With Penalty | No Penalty | Difference | Rank Change |")
    print("|---------|--------------|------------|------------|-------------|")
    
    for _, row in df_comparison.sort_values('Variant').iterrows():
        rank_symbol = "↑" if row['Rank_Change'] > 0 else ("↓" if row['Rank_Change'] < 0 else "=")
        print(f"| {row['Variant']:7s} | {row['Mean_Accessibility_Penalty']:.4f} (#{row['Rank_Penalty']:2d}) | "
              f"{row['Mean_Accessibility_NoPenalty']:.4f} (#{row['Rank_NoPenalty']:2d}) | "
              f"{row['Diff']:+.4f} | {rank_symbol} {abs(row['Rank_Change'])} |")
    

    df_penalty.to_csv(tables_dir / 'accessibility_analysis_with_penalty.csv', index=False)
    df_no_penalty.to_csv(tables_dir / 'accessibility_analysis_no_penalty.csv', index=False)
    df_comparison.to_csv(tables_dir / 'accessibility_comparison.csv', index=False)
    
    print("\n✅ 分析结果已保存到:")
    print(f"  - {tables_dir / 'accessibility_analysis_with_penalty.csv'}")
    print(f"  - {tables_dir / 'accessibility_analysis_no_penalty.csv'}")
    print(f"  - {tables_dir / 'accessibility_comparison.csv'}")

if __name__ == '__main__':
    main()