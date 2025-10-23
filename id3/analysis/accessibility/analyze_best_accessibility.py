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

def analyze_best_accessibility(results: List[Dict]) -> pd.DataFrame:
    """分析best_accessibility性能"""
    data = []
    
    for result in results:
        data.append({
            'protein': result['protein_name'],
            'constraint_type': result['constraint_type'],
            'variant': result['variant'],
            'best_accessibility': result['best_accessibility'],
            'final_accessibility': result['final_accessibility'],
            'initial_accessibility': result['initial_accessibility'],
            'improvement': result['best_accessibility'] - result['initial_accessibility'],
            'amino_acids_match': result['amino_acids_match'],
            'iterations': result.get('iterations', 1000)
        })
    
    return pd.DataFrame(data)

def print_best_accessibility_analysis(df: pd.DataFrame):

    print("\n" + "="*80)

    print("="*80)
    






    


    print("-"*70)
    
    constraint_stats = df.groupby('constraint_type')['best_accessibility'].agg([





    ]).round(4)
    
    print(constraint_stats.to_string())
    


    print("-"*70)
    
    variant_stats = df.groupby('variant')['best_accessibility'].agg([





    ]).round(4)
    
    print(variant_stats.to_string())
    


    print("-"*70)
    
    combined_stats = df.groupby(['constraint_type', 'variant'])['best_accessibility'].agg([



    ]).round(4)
    
    print(combined_stats.to_string())
    


    print("-"*70)
    
    top10 = df.nsmallest(10, 'best_accessibility')[
        ['protein', 'constraint_type', 'variant', 'best_accessibility', 'improvement']
    ]
    
    for idx, row in top10.iterrows():
        print(f"{row['protein']:8s} | {row['constraint_type']:10s} | {row['variant']:2s} | "

    


    print("-"*70)
    
    best_by_protein = []
    for protein in df['protein'].unique():
        protein_data = df[df['protein'] == protein]
        best_idx = protein_data['best_accessibility'].idxmin()
        best_config = protein_data.loc[best_idx]
        best_by_protein.append({
            'protein': protein,
            'constraint': best_config['constraint_type'],
            'variant': best_config['variant'],
            'best_accessibility': best_config['best_accessibility']
        })
    
    best_by_protein_df = pd.DataFrame(best_by_protein).sort_values('best_accessibility')
    
    for _, row in best_by_protein_df.iterrows():
        print(f"{row['protein']:8s}: {row['best_accessibility']:.4f} kcal/mol "

    


    print("-"*70)
    

    best_constraint_counts = best_by_protein_df['constraint'].value_counts()

    for constraint, count in best_constraint_counts.items():

    

    best_variant_counts = best_by_protein_df['variant'].value_counts()

    for variant, count in best_variant_counts.items():


def create_best_accessibility_visualizations(df: pd.DataFrame, output_dir: Path):
    """创建基于best_accessibility的可视化"""
    output_dir.mkdir(parents=True, exist_ok=True)
    

    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_palette("husl")
    

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    

    sns.boxplot(data=df, x='constraint_type', y='best_accessibility', ax=axes[0, 0])
    axes[0, 0].set_title('Best Accessibility by Constraint Type', fontsize=14, fontweight='bold')
    axes[0, 0].set_ylabel('Best Accessibility (kcal/mol)')
    axes[0, 0].set_xlabel('Constraint Type')
    

    sns.boxplot(data=df, x='variant', y='best_accessibility', ax=axes[0, 1])
    axes[0, 1].set_title('Best Accessibility by Variant', fontsize=14, fontweight='bold')
    axes[0, 1].set_ylabel('Best Accessibility (kcal/mol)')
    axes[0, 1].set_xlabel('Variant')
    

    sns.boxplot(data=df, x='variant', y='best_accessibility', hue='constraint_type', ax=axes[0, 2])
    axes[0, 2].set_title('Best Accessibility by Variant and Constraint', fontsize=14, fontweight='bold')
    axes[0, 2].set_ylabel('Best Accessibility (kcal/mol)')
    axes[0, 2].legend(title='Constraint', bbox_to_anchor=(1.05, 1), loc='upper left')
    

    sns.scatterplot(data=df, x='final_accessibility', y='best_accessibility', 
                    hue='constraint_type', style='variant', s=100, ax=axes[1, 0])
    axes[1, 0].set_title('Best vs Final Accessibility', fontsize=14, fontweight='bold')
    axes[1, 0].set_xlabel('Final Accessibility (kcal/mol)')
    axes[1, 0].set_ylabel('Best Accessibility (kcal/mol)')
    axes[1, 0].plot([0, 7], [0, 7], 'k--', alpha=0.3)
    

    df['best_improvement'] = df['initial_accessibility'] - df['best_accessibility']
    sns.boxplot(data=df, x='constraint_type', y='best_improvement', ax=axes[1, 1])
    axes[1, 1].set_title('Improvement (from Initial to Best)', fontsize=14, fontweight='bold')
    axes[1, 1].set_ylabel('Improvement (kcal/mol)')
    axes[1, 1].set_xlabel('Constraint Type')
    

    best_per_protein = df.groupby('protein')['best_accessibility'].min().sort_values()
    best_per_protein.plot(kind='barh', ax=axes[1, 2])
    axes[1, 2].set_title('Best Accessibility by Protein (Minimum across all configs)', fontsize=14, fontweight='bold')
    axes[1, 2].set_xlabel('Best Accessibility (kcal/mol)')
    axes[1, 2].set_ylabel('Protein')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'best_accessibility_analysis.png', dpi=150, bbox_inches='tight')
    plt.close()
    

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    

    pivot_best = df.pivot_table(
        values='best_accessibility',
        index='protein',
        columns=['constraint_type', 'variant'],
        aggfunc='mean'
    )
    
    sns.heatmap(pivot_best, annot=True, fmt='.3f', cmap='RdYlGn_r', 
                cbar_kws={'label': 'Best Accessibility (kcal/mol)'}, ax=ax1)
    ax1.set_title('Best Accessibility Heatmap', fontsize=14, fontweight='bold')
    

    pivot_improvement = df.pivot_table(
        values='best_improvement',
        index='protein',
        columns=['constraint_type', 'variant'],
        aggfunc='mean'
    )
    
    sns.heatmap(pivot_improvement, annot=True, fmt='.3f', cmap='YlOrRd', 
                cbar_kws={'label': 'Improvement (kcal/mol)'}, ax=ax2)
    ax2.set_title('Improvement Heatmap (Initial - Best)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'performance_heatmaps.png', dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\n📁 可视化图表已保存到: {output_dir}")

def save_best_accessibility_summary(df: pd.DataFrame, output_dir: Path):

    

    summary = {
        'global_metrics': {
            'total_experiments': len(df),
            'global_best_accessibility': float(df['best_accessibility'].min()),
            'mean_best_accessibility': float(df['best_accessibility'].mean()),
            'std_best_accessibility': float(df['best_accessibility'].std())
        },
        'by_constraint': {},
        'by_variant': {},
        'by_combination': {},
        'best_configs_per_protein': {},
        'top_10_configs': []
    }
    

    for constraint in df['constraint_type'].unique():
        constraint_data = df[df['constraint_type'] == constraint]
        summary['by_constraint'][constraint] = {
            'mean': float(constraint_data['best_accessibility'].mean()),
            'min': float(constraint_data['best_accessibility'].min()),
            'max': float(constraint_data['best_accessibility'].max()),
            'std': float(constraint_data['best_accessibility'].std())
        }
    

    for variant in df['variant'].unique():
        variant_data = df[df['variant'] == variant]
        summary['by_variant'][variant] = {
            'mean': float(variant_data['best_accessibility'].mean()),
            'min': float(variant_data['best_accessibility'].min()),
            'max': float(variant_data['best_accessibility'].max()),
            'std': float(variant_data['best_accessibility'].std())
        }
    

    for (constraint, variant), group in df.groupby(['constraint_type', 'variant']):
        key = f"{constraint}_{variant}"
        summary['by_combination'][key] = {
            'mean': float(group['best_accessibility'].mean()),
            'min': float(group['best_accessibility'].min()),
            'count': len(group)
        }
    

    for protein in df['protein'].unique():
        protein_data = df[df['protein'] == protein]
        best_idx = protein_data['best_accessibility'].idxmin()
        best_config = protein_data.loc[best_idx]
        summary['best_configs_per_protein'][protein] = {
            'constraint_type': best_config['constraint_type'],
            'variant': best_config['variant'],
            'best_accessibility': float(best_config['best_accessibility']),
            'improvement': float(best_config['improvement'])
        }
    

    top10 = df.nsmallest(10, 'best_accessibility')
    for _, row in top10.iterrows():
        summary['top_10_configs'].append({
            'protein': row['protein'],
            'constraint_type': row['constraint_type'],
            'variant': row['variant'],
            'best_accessibility': float(row['best_accessibility'])
        })
    

    with open(output_dir / 'best_accessibility_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    


def main():

    results_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_003751_unified_cai_experiments")
    output_dir = Path("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250909_003751_analysis")
    


    results = load_experiment_results(results_dir)

    

    df = analyze_best_accessibility(results)
    

    print_best_accessibility_analysis(df)
    


    create_best_accessibility_visualizations(df, output_dir)
    

    save_best_accessibility_summary(df, output_dir)
    


if __name__ == "__main__":
    main()