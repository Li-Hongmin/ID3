#!/usr/bin/env python3
"""
Generate fig:convergence_performance - ID3 Framework Performance Analysis




"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')


plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("tab20")


AMS_COLOR = '#2E86AB'  # Blue for Amino Matching Softmax
CPC_COLOR = '#A23B72'  # Purple for Codon Profile Constraint  
LAG_COLOR = '#F18F01'  # Orange for Lagrangian Multiplier


STRATEGY_COLORS = {
    '00': '
    '01': '
    '10': '
    '11': '
}


STRATEGY_LABELS = {
    '00': 'Det. Soft',
    '01': 'Det. Hard', 
    '10': 'Stoch. Soft',
    '11': 'Stoch. Hard'
}

def load_table_data():

    tables_dir = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/tables')
    

    tab1_path = tables_dir / 'data_tab1_access_comparison.csv'
    access_df = pd.read_csv(tab1_path)
    

    tab2_path = tables_dir / 'data_tab2_access_cai_comparison.json'
    with open(tab2_path, 'r') as f:
        tab2_data = json.load(f)
    

    cai_results = []
    for protein, experiments in tab2_data.items():
        for exp_key, exp_data in experiments.items():
            constraint, variant = exp_key.split('_')
            cai_results.append({
                'protein': protein,
                'constraint': constraint,
                'variant': variant,
                'accessibility': exp_data['accessibility'],
                'cai': exp_data['discrete_cai']
            })
    
    cai_df = pd.DataFrame(cai_results)
    


    access_results = []
    proteins = access_df['Protein'].values
    
    for idx, protein in enumerate(proteins):
        row = access_df.iloc[idx]

            if '_' in col:
                constraint, variant = col.rsplit('_', 1)
                value = row[col]
                if pd.notna(value):
                    access_results.append({
                        'protein': protein,
                        'constraint': constraint,
                        'variant': variant,
                        'accessibility': float(value)
                    })
    
    access_df_processed = pd.DataFrame(access_results)
    
    return access_df_processed, cai_df

def create_performance_figure(access_df: pd.DataFrame, cai_df: pd.DataFrame, output_dir: Path):
    """创建性能分析图表"""
    

    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(1, 2, hspace=0.3, wspace=0.3)
    

    print("生成Panel A: 约束机制对比...")
    ax1 = fig.add_subplot(gs[0, 0])
    

    constraints = ['lagrangian', 'ams', 'cpc']
    constraint_labels = {
        'lagrangian': 'Lagrangian\nMultiplier',
        'ams': 'Amino Matching\nSoftmax',
        'cpc': 'Codon Profile\nConstraint'
    }
    constraint_colors = {
        'lagrangian': LAG_COLOR,
        'ams': AMS_COLOR,
        'cpc': CPC_COLOR
    }
    

    access_means = []
    access_stds = []
    cai_means = []
    cai_stds = []
    
    for constraint in constraints:

        access_data = access_df[access_df['constraint'] == constraint]['accessibility']
        access_means.append(access_data.mean())
        access_stds.append(access_data.std())
        

        cai_data = cai_df[cai_df['constraint'] == constraint]['accessibility']
        cai_means.append(cai_data.mean())
        cai_stds.append(cai_data.std())
    

    x = np.arange(len(constraints))
    width = 0.35
    

    bars1 = ax1.bar(x - width/2, access_means, width, yerr=access_stds,
                    label='Access-only', capsize=5, alpha=0.8, 
                    edgecolor='black', linewidth=1.5)
    bars2 = ax1.bar(x + width/2, cai_means, width, yerr=cai_stds,
                    label='Access+CAI', capsize=5, alpha=0.8,
                    edgecolor='black', linewidth=1.5, hatch='//')
    

    for i, constraint in enumerate(constraints):
        bars1[i].set_facecolor(constraint_colors[constraint])
        bars2[i].set_facecolor(constraint_colors[constraint])
    

    ax1.set_xlabel('Constraint Mechanism', fontsize=12, fontweight='bold')
    ax1.set_ylabel('RNA Accessibility (kcal/mol)', fontsize=12, fontweight='bold')
    ax1.set_title('A) Performance by Constraint Mechanism', fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels([constraint_labels[c] for c in constraints])
    ax1.legend(loc='upper right', fontsize=10)
    ax1.grid(True, alpha=0.3, axis='y')
    

    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=8)
    

    print("生成Panel B: 策略性能分析...")
    ax2 = fig.add_subplot(gs[0, 1])
    

    variants = ['00', '01', '10', '11']
    

    access_strategy_means = []
    access_strategy_stds = []
    cai_strategy_means = []
    cai_strategy_stds = []
    
    for variant in variants:

        access_data = access_df[access_df['variant'] == variant]['accessibility']
        access_strategy_means.append(access_data.mean())
        access_strategy_stds.append(access_data.std())
        

        cai_data = cai_df[cai_df['variant'] == variant]['accessibility']
        cai_strategy_means.append(cai_data.mean())
        cai_strategy_stds.append(cai_data.std())
    

    x = np.arange(len(variants))
    

    bars3 = ax2.bar(x - width/2, access_strategy_means, width, yerr=access_strategy_stds,
                    label='Access-only', capsize=5, alpha=0.8,
                    edgecolor='black', linewidth=1.5)
    bars4 = ax2.bar(x + width/2, cai_strategy_means, width, yerr=cai_strategy_stds,
                    label='Access+CAI', capsize=5, alpha=0.8,
                    edgecolor='black', linewidth=1.5, hatch='//')
    

    for i, variant in enumerate(variants):
        bars3[i].set_facecolor(STRATEGY_COLORS[variant])
        bars4[i].set_facecolor(STRATEGY_COLORS[variant])
    

    ax2.set_xlabel('Optimization Strategy', fontsize=12, fontweight='bold')
    ax2.set_ylabel('RNA Accessibility (kcal/mol)', fontsize=12, fontweight='bold')
    ax2.set_title('B) Performance by Strategy', fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels([STRATEGY_LABELS[v] for v in variants])
    ax2.legend(loc='upper right', fontsize=10)
    ax2.grid(True, alpha=0.3, axis='y')
    

    for bars in [bars3, bars4]:
        for bar in bars:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.3f}', ha='center', va='bottom', fontsize=8)
    

    fig.suptitle('ID3 Framework Performance Analysis: Access-only vs Access+CAI', 
                 fontsize=16, fontweight='bold', y=1.02)
    

    plt.tight_layout()
    

    output_dir.mkdir(parents=True, exist_ok=True)
    
    png_path = output_dir / 'fig_convergence_performance.png'
    pdf_path = output_dir / 'fig_convergence_performance.pdf'
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
    
    print(f"✅ 图表已保存到: {png_path}")
    print(f"✅ 图表已保存到: {pdf_path}")
    

    # plt.show()
    plt.close()
    
    return fig

def generate_performance_summary(access_df: pd.DataFrame, cai_df: pd.DataFrame) -> dict:

    summary = {}
    

    summary['overall'] = {
        'access_only': {
            'mean': access_df['accessibility'].mean(),
            'std': access_df['accessibility'].std(),
            'min': access_df['accessibility'].min(),
            'max': access_df['accessibility'].max(),
            'count': len(access_df)
        },
        'access_cai': {
            'mean': cai_df['accessibility'].mean(),
            'std': cai_df['accessibility'].std(),
            'min': cai_df['accessibility'].min(),
            'max': cai_df['accessibility'].max(),
            'mean_cai': cai_df['cai'].mean(),
            'count': len(cai_df)
        }
    }
    

    summary['by_constraint'] = {}
    for constraint in ['lagrangian', 'ams', 'cpc']:
        summary['by_constraint'][constraint] = {
            'access_only': {
                'mean': access_df[access_df['constraint'] == constraint]['accessibility'].mean(),
                'std': access_df[access_df['constraint'] == constraint]['accessibility'].std(),
                'count': len(access_df[access_df['constraint'] == constraint])
            },
            'access_cai': {
                'mean': cai_df[cai_df['constraint'] == constraint]['accessibility'].mean(),
                'std': cai_df[cai_df['constraint'] == constraint]['accessibility'].std(),
                'mean_cai': cai_df[cai_df['constraint'] == constraint]['cai'].mean(),
                'count': len(cai_df[cai_df['constraint'] == constraint])
            }
        }
    

    summary['by_strategy'] = {}
    for variant in ['00', '01', '10', '11']:
        summary['by_strategy'][variant] = {
            'access_only': {
                'mean': access_df[access_df['variant'] == variant]['accessibility'].mean(),
                'std': access_df[access_df['variant'] == variant]['accessibility'].std(),
                'count': len(access_df[access_df['variant'] == variant])
            },
            'access_cai': {
                'mean': cai_df[cai_df['variant'] == variant]['accessibility'].mean(),
                'std': cai_df[cai_df['variant'] == variant]['accessibility'].std(),
                'mean_cai': cai_df[cai_df['variant'] == variant]['cai'].mean(),
                'count': len(cai_df[cai_df['variant'] == variant])
            }
        }
    

    access_mean = summary['overall']['access_only']['mean']
    cai_mean = summary['overall']['access_cai']['mean']
    summary['improvement'] = {
        'absolute': access_mean - cai_mean,
        'percentage': ((access_mean - cai_mean) / access_mean) * 100
    }
    
    return summary

def main():
    """主函数"""
    output_dir = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/paper_experiment_results/figures')
    
    print("="*60)
    print("ID3 Framework Performance Analysis (Fast Version)")
    print("="*60)
    

    print("\n📊 加载表格数据...")
    access_df, cai_df = load_table_data()
    
    print(f"Access-only数据: {len(access_df)} 条记录")
    print(f"Access+CAI数据: {len(cai_df)} 条记录")
    

    print("\n🎨 生成性能分析图表...")
    create_performance_figure(access_df, cai_df, output_dir)
    

    print("\n📈 性能统计摘要:")
    summary = generate_performance_summary(access_df, cai_df)
    
    print(f"\n总体性能:")
    print(f"  Access-only: {summary['overall']['access_only']['count']} 个实验")
    print(f"    平均: {summary['overall']['access_only']['mean']:.4f} kcal/mol")
    print(f"    标准差: {summary['overall']['access_only']['std']:.4f}")
    print(f"  Access+CAI: {summary['overall']['access_cai']['count']} 个实验")
    print(f"    平均: {summary['overall']['access_cai']['mean']:.4f} kcal/mol")
    print(f"    标准差: {summary['overall']['access_cai']['std']:.4f}")
    print(f"    平均CAI: {summary['overall']['access_cai']['mean_cai']:.4f}")
    print(f"  改进: {summary['improvement']['absolute']:.4f} kcal/mol ({summary['improvement']['percentage']:.1f}%)")
    
    print(f"\n按约束类型:")
    for constraint in ['lagrangian', 'ams', 'cpc']:
        stats = summary['by_constraint'][constraint]
        print(f"  {constraint.upper()}:")
        print(f"    Access-only: {stats['access_only']['mean']:.4f} ± {stats['access_only']['std']:.4f}")
        print(f"    Access+CAI: {stats['access_cai']['mean']:.4f} ± {stats['access_cai']['std']:.4f} (CAI: {stats['access_cai']['mean_cai']:.3f})")
    
    print(f"\n按策略类型:")
    for variant in ['00', '01', '10', '11']:
        stats = summary['by_strategy'][variant]
        label = STRATEGY_LABELS[variant]
        print(f"  {label} ({variant}):")
        print(f"    Access-only: {stats['access_only']['mean']:.4f} ± {stats['access_only']['std']:.4f}")
        print(f"    Access+CAI: {stats['access_cai']['mean']:.4f} ± {stats['access_cai']['std']:.4f} (CAI: {stats['access_cai']['mean_cai']:.3f})")
    

    summary_path = output_dir / 'performance_summary.json'
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n✅ 摘要已保存到: {summary_path}")
    
    print("\n✨ 分析完成!")

if __name__ == "__main__":
    main()