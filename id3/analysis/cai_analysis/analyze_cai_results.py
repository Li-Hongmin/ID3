#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List

def load_experiment_results(results_dir: Path) -> List[Dict]:

    results = []
    
    for file_path in results_dir.glob("*.json"):

        if file_path.name in ["config.json", "summary.json"]:
            continue
            
        with open(file_path) as f:
            data = json.load(f)
            results.append(data)
    
    return results

def analyze_performance_by_constraint(results: List[Dict]) -> pd.DataFrame:
    """按约束类型分析性能"""
    data = []
    
    for result in results:

        final_cai = None
        if 'trajectory' in result and 'cai_loss' in result['trajectory']:
            cai_losses = result['trajectory']['cai_loss']
            if cai_losses:
                final_cai = cai_losses[-1] if isinstance(cai_losses[-1], (int, float)) else None
        
        data.append({
            'protein': result['protein_name'],
            'constraint_type': result['constraint_type'],
            'variant': result['variant'],
            'initial_accessibility': result['initial_accessibility'],
            'final_accessibility': result['final_accessibility'],
            'improvement': result['improvement'],
            'best_accessibility': result['best_accessibility'],
            'amino_acids_match': result['amino_acids_match'],
            'final_cai_loss': final_cai,
            'cai_enabled': result.get('cai_enabled', False)
        })
    
    return pd.DataFrame(data)

def print_summary_statistics(df: pd.DataFrame):

    print("\n" + "="*80)

    print("="*80)
    




    


    print("-"*60)
    
    constraint_stats = df.groupby('constraint_type').agg({
        'final_accessibility': ['mean', 'std', 'min', 'max'],
        'improvement': 'mean',
        'amino_acids_match': 'mean',
        'final_cai_loss': 'mean'
    }).round(4)
    
    print(constraint_stats)
    


    print("-"*60)
    
    variant_stats = df.groupby('variant').agg({
        'final_accessibility': ['mean', 'std', 'min', 'max'],
        'improvement': 'mean',
        'amino_acids_match': 'mean',
        'final_cai_loss': 'mean'
    }).round(4)
    
    print(variant_stats)
    


    print("-"*60)
    
    best_by_protein = df.groupby('protein')['final_accessibility'].min().sort_values()
    for protein, accessibility in best_by_protein.items():
        best_config = df[df['protein'] == protein].nsmallest(1, 'final_accessibility').iloc[0]
        print(f"{protein:8s}: {accessibility:.4f} kcal/mol "

    


    print("-"*60)
    

    df_with_cai = df[df['final_cai_loss'].notna()]
    
    if not df_with_cai.empty:



        

        cai_by_constraint = df_with_cai.groupby('constraint_type')['final_cai_loss'].agg(['mean', 'std'])

        print(cai_by_constraint)
    else:


def create_visualizations(df: pd.DataFrame, output_dir: Path):
    """创建可视化图表"""
    output_dir.mkdir(parents=True, exist_ok=True)
    

    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    

    sns.boxplot(data=df, x='constraint_type', y='final_accessibility', ax=axes[0, 0])
    axes[0, 0].set_title('Final Accessibility by Constraint Type')
    axes[0, 0].set_ylabel('Accessibility (kcal/mol)')
    

    sns.boxplot(data=df, x='constraint_type', y='improvement', ax=axes[0, 1])
    axes[0, 1].set_title('Improvement by Constraint Type')
    axes[0, 1].set_ylabel('Improvement (kcal/mol)')
    

    sns.boxplot(data=df, x='variant', y='final_accessibility', hue='constraint_type', ax=axes[1, 0])
    axes[1, 0].set_title('Final Accessibility by Variant and Constraint')
    axes[1, 0].set_ylabel('Accessibility (kcal/mol)')
    axes[1, 0].legend(title='Constraint', bbox_to_anchor=(1.05, 1), loc='upper left')
    

    df_with_cai = df[df['final_cai_loss'].notna()]
    if not df_with_cai.empty:
        sns.boxplot(data=df_with_cai, x='constraint_type', y='final_cai_loss', ax=axes[1, 1])
        axes[1, 1].set_title('CAI Loss by Constraint Type')
        axes[1, 1].set_ylabel('CAI Loss')
    else:
        axes[1, 1].text(0.5, 0.5, 'No CAI Loss Data Available', 
                       ha='center', va='center', transform=axes[1, 1].transAxes)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'constraint_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()
    

    fig, ax = plt.subplots(figsize=(14, 10))
    

    pivot = df.pivot_table(
        values='final_accessibility',
        index='protein',
        columns=['constraint_type', 'variant'],
        aggfunc='mean'
    )
    
    sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn_r', 
                cbar_kws={'label': 'Accessibility (kcal/mol)'}, ax=ax)
    ax.set_title('Final Accessibility Heatmap by Protein, Constraint Type and Variant')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'protein_performance_heatmap.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n可视化图表已保存到: {output_dir}")

def main():

    results_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_003751_unified_cai_experiments")
    output_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_003751_analysis")
    

    print("正在加载实验结果...")
    results = load_experiment_results(results_dir)
    print(f"成功加载 {len(results)} 个实验结果")
    

    df = analyze_performance_by_constraint(results)
    

    print_summary_statistics(df)
    

    print("\n正在生成可视化图表...")
    create_visualizations(df, output_dir)
    

    analysis_results = {
        'total_experiments': len(df),
        'cai_enabled': df['cai_enabled'].all(),
        'amino_acid_match_rate': df['amino_acids_match'].mean(),
        'best_accessibility': df['final_accessibility'].min(),
        'mean_accessibility': df['final_accessibility'].mean(),
        'mean_improvement': df['improvement'].mean(),
        'best_configs': {}
    }
    

    for protein in df['protein'].unique():
        best = df[df['protein'] == protein].nsmallest(1, 'final_accessibility').iloc[0]
        analysis_results['best_configs'][protein] = {
            'constraint_type': best['constraint_type'],
            'variant': best['variant'],
            'final_accessibility': best['final_accessibility'],
            'improvement': best['improvement']
        }
    

    with open(output_dir / 'analysis_results.json', 'w') as f:
        json.dump(analysis_results, f, indent=2, default=float)
    
    print(f"\n分析结果已保存到: {output_dir / 'analysis_results.json'}")
    print("\n分析完成！")

if __name__ == "__main__":
    main()