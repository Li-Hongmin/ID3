#!/usr/bin/env python3
"""


"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
from typing import Dict, List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')


plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

class PerformanceVisualizer:

    
    def __init__(self):
        """初始化可视化器"""
        sns.set_style("whitegrid")
        self.colors = {
            'lagrangian': '
            'ams': '
            'cpc': '
        }
        self.variant_colors = {
            '00': '
            '01': '
            '10': '
            '11': '
        }
    
    def load_experiment_data(self, dir_path: Path) -> pd.DataFrame:

        results = []
        

        json_files = list(dir_path.glob("*.json"))
        for json_file in json_files:
            if json_file.name.startswith(('config', 'progress', 'summary')):
                continue
            
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    results.append({
                        'protein': data.get('protein_name', 'unknown'),
                        'constraint': data.get('constraint_type', 'unknown'),
                        'variant': data.get('variant', 'unknown'),
                        'accessibility': data.get('best_accessibility', None),
                        'initial_acc': data.get('initial_accessibility', None),
                        'improvement': data.get('improvement', None),
                        'aa_correct': data.get('amino_acids_correct', 0),
                        'time': data.get('optimization_time', 0),
                        'iterations': data.get('iterations', 0),
                        'cai_enabled': data.get('cai_enabled', False),
                        'final_cai': data.get('final_ecai', None)
                    })
            except Exception as e:
                print(f"Error loading {json_file.name}: {e}")
        
        return pd.DataFrame(results)
    
    def plot_constraint_comparison(self, df: pd.DataFrame, save_path: Optional[Path] = None):
        """绘制约束类型对比图"""
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        

        ax1 = axes[0]
        constraint_order = ['lagrangian', 'ams', 'cpc']
        df_sorted = df[df['constraint'].isin(constraint_order)]
        
        box_plot = ax1.boxplot([df_sorted[df_sorted['constraint'] == c]['accessibility'].values 
                                for c in constraint_order],
                               labels=constraint_order,
                               patch_artist=True)
        
        for patch, constraint in zip(box_plot['boxes'], constraint_order):
            patch.set_facecolor(self.colors[constraint])
            patch.set_alpha(0.7)
        
        ax1.set_xlabel('Constraint Type')
        ax1.set_ylabel('Accessibility (kcal/mol)')
        ax1.set_title('Performance by Constraint Type')
        ax1.grid(True, alpha=0.3)
        

        means = [df_sorted[df_sorted['constraint'] == c]['accessibility'].mean() 
                for c in constraint_order]
        ax1.plot(range(1, len(constraint_order) + 1), means, 'r--', label='Mean')
        ax1.legend()
        

        ax2 = axes[1]
        variant_order = ['00', '01', '10', '11']
        df_variant = df[df['variant'].isin(variant_order)]
        
        positions = range(len(variant_order))
        violin_parts = ax2.violinplot(
            [df_variant[df_variant['variant'] == v]['accessibility'].values 
             for v in variant_order],
            positions=positions,
            showmeans=True,
            showmedians=True
        )
        
        ax2.set_xticks(positions)
        ax2.set_xticklabels(variant_order)
        ax2.set_xlabel('Variant')
        ax2.set_ylabel('Accessibility (kcal/mol)')
        ax2.set_title('Performance by Variant')
        ax2.grid(True, alpha=0.3)
        

        ax3 = axes[2]
        pivot_table = df.pivot_table(
            values='accessibility',
            index='constraint',
            columns='variant',
            aggfunc='mean'
        )
        
        im = ax3.imshow(pivot_table.values, cmap='RdYlGn_r', aspect='auto')
        ax3.set_xticks(range(len(variant_order)))
        ax3.set_xticklabels(variant_order)
        ax3.set_yticks(range(len(constraint_order)))
        ax3.set_yticklabels(constraint_order)
        ax3.set_xlabel('Variant')
        ax3.set_ylabel('Constraint Type')
        ax3.set_title('Mean Accessibility Heatmap')
        

        for i in range(len(constraint_order)):
            for j in range(len(variant_order)):
                value = pivot_table.iloc[i, j]
                if not np.isnan(value):
                    ax3.text(j, i, f'{value:.2f}', ha='center', va='center', color='white')
        
        plt.colorbar(im, ax=ax3, label='Accessibility (kcal/mol)')
        
        plt.suptitle('Experimental Performance Analysis', fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def plot_performance_distribution(self, df: pd.DataFrame, save_path: Optional[Path] = None):

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        

        ax1 = axes[0, 0]
        ax1.hist(df['accessibility'], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
        ax1.axvline(df['accessibility'].mean(), color='red', linestyle='--', label=f'Mean: {df["accessibility"].mean():.3f}')
        ax1.axvline(df['accessibility'].median(), color='green', linestyle='--', label=f'Median: {df["accessibility"].median():.3f}')
        ax1.set_xlabel('Accessibility (kcal/mol)')
        ax1.set_ylabel('Count')
        ax1.set_title('Overall Performance Distribution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        

        ax2 = axes[0, 1]
        for constraint in ['lagrangian', 'ams', 'cpc']:
            data = df[df['constraint'] == constraint]['accessibility']
            if len(data) > 0:
                data.plot.kde(ax=ax2, label=constraint.upper(), color=self.colors[constraint], linewidth=2)
        ax2.set_xlabel('Accessibility (kcal/mol)')
        ax2.set_ylabel('Density')
        ax2.set_title('Performance Density by Constraint')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        

        ax3 = axes[1, 0]
        for constraint in ['lagrangian', 'ams', 'cpc']:
            data = df[df['constraint'] == constraint]['accessibility'].sort_values()
            if len(data) > 0:
                y = np.arange(1, len(data) + 1) / len(data)
                ax3.plot(data, y, label=constraint.upper(), color=self.colors[constraint], linewidth=2)
        

        ax3.axvline(1.0, color='gray', linestyle=':', alpha=0.5, label='Excellent (<1.0)')
        ax3.axvline(1.5, color='gray', linestyle='--', alpha=0.5, label='Good (<1.5)')
        ax3.set_xlabel('Accessibility (kcal/mol)')
        ax3.set_ylabel('Cumulative Probability')
        ax3.set_title('Cumulative Distribution Function')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        

        ax4 = axes[1, 1]
        if 'improvement' in df.columns and df['improvement'].notna().any():
            for constraint in ['lagrangian', 'ams', 'cpc']:
                data = df[df['constraint'] == constraint]
                if len(data) > 0:
                    ax4.scatter(data['initial_acc'], data['accessibility'], 
                              label=constraint.upper(), color=self.colors[constraint], alpha=0.6)
            

            min_val = min(df['initial_acc'].min(), df['accessibility'].min())
            max_val = max(df['initial_acc'].max(), df['accessibility'].max())
            ax4.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.3, label='No improvement')
            
            ax4.set_xlabel('Initial Accessibility (kcal/mol)')
            ax4.set_ylabel('Final Accessibility (kcal/mol)')
            ax4.set_title('Optimization Improvement')
            ax4.legend()
            ax4.grid(True, alpha=0.3)
        
        plt.suptitle('Performance Distribution Analysis', fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def plot_protein_performance(self, df: pd.DataFrame, save_path: Optional[Path] = None):
        """绘制蛋白质性能对比"""
        fig, ax = plt.subplots(figsize=(14, 6))
        

        protein_means = df.groupby('protein')['accessibility'].mean().sort_values()
        

        bars = ax.bar(range(len(protein_means)), protein_means.values, color='steelblue', alpha=0.7)
        

        protein_stds = df.groupby('protein')['accessibility'].std()
        protein_stds = protein_stds.reindex(protein_means.index)
        ax.errorbar(range(len(protein_means)), protein_means.values, 
                   yerr=protein_stds.values, fmt='none', color='black', capsize=3)
        

        ax.set_xticks(range(len(protein_means)))
        ax.set_xticklabels(protein_means.index, rotation=45, ha='right')
        ax.set_xlabel('Protein')
        ax.set_ylabel('Mean Accessibility (kcal/mol)')
        ax.set_title('Performance by Protein')
        

        ax.axhline(1.0, color='green', linestyle='--', alpha=0.5, label='Excellent')
        ax.axhline(1.5, color='orange', linestyle='--', alpha=0.5, label='Good')
        ax.axhline(2.0, color='red', linestyle='--', alpha=0.5, label='Acceptable')
        
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def generate_performance_report(self, df: pd.DataFrame) -> Dict:

        report = {
            'overall': {
                'mean': df['accessibility'].mean(),
                'median': df['accessibility'].median(),
                'std': df['accessibility'].std(),
                'min': df['accessibility'].min(),
                'max': df['accessibility'].max(),
                'excellent_rate': (df['accessibility'] < 1.0).sum() / len(df) * 100,
                'good_rate': (df['accessibility'] < 1.5).sum() / len(df) * 100
            },
            'by_constraint': {},
            'by_variant': {},
            'by_protein': {}
        }
        

        for constraint in df['constraint'].unique():
            data = df[df['constraint'] == constraint]['accessibility']
            report['by_constraint'][constraint] = {
                'mean': data.mean(),
                'median': data.median(),
                'std': data.std(),
                'min': data.min(),
                'max': data.max(),
                'count': len(data)
            }
        

        for variant in df['variant'].unique():
            data = df[df['variant'] == variant]['accessibility']
            report['by_variant'][variant] = {
                'mean': data.mean(),
                'median': data.median(),
                'count': len(data)
            }
        

        protein_means = df.groupby('protein')['accessibility'].mean().sort_values()
        report['by_protein']['best_5'] = protein_means.head(5).to_dict()
        report['by_protein']['worst_5'] = protein_means.tail(5).to_dict()
        
        return report
    
    def compare_experiments(self, experiments: List[Tuple[str, Path]], save_path: Optional[Path] = None):
        """对比多个实验批次"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        all_data = []
        for exp_name, exp_path in experiments:
            df = self.load_experiment_data(exp_path)
            df['experiment'] = exp_name
            all_data.append(df)
        
        combined_df = pd.concat(all_data, ignore_index=True)
        

        ax1 = axes[0, 0]
        exp_names = [name for name, _ in experiments]
        box_data = [combined_df[combined_df['experiment'] == name]['accessibility'].values 
                   for name in exp_names]
        
        bp = ax1.boxplot(box_data, labels=exp_names, patch_artist=True)
        for patch, color in zip(bp['boxes'], plt.cm.Set3.colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax1.set_ylabel('Accessibility (kcal/mol)')
        ax1.set_title('Overall Performance Comparison')
        ax1.grid(True, alpha=0.3)
        

        ax2 = axes[0, 1]
        for i, (exp_name, _) in enumerate(experiments):
            exp_data = combined_df[combined_df['experiment'] == exp_name]
            constraint_means = exp_data.groupby('constraint')['accessibility'].mean()
            x = np.arange(len(constraint_means)) + i * 0.25
            ax2.bar(x, constraint_means.values, width=0.2, label=exp_name, alpha=0.7)
        
        ax2.set_xticks(np.arange(len(constraint_means)) + 0.25)
        ax2.set_xticklabels(constraint_means.index)
        ax2.set_ylabel('Mean Accessibility (kcal/mol)')
        ax2.set_title('Performance by Constraint Type')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        

        ax3 = axes[1, 0]
        summary_data = []
        for exp_name, _ in experiments:
            exp_data = combined_df[combined_df['experiment'] == exp_name]
            summary_data.append({
                'experiment': exp_name,
                'mean': exp_data['accessibility'].mean(),
                'std': exp_data['accessibility'].std()
            })
        
        summary_df = pd.DataFrame(summary_data)
        ax3.errorbar(range(len(summary_df)), summary_df['mean'], 
                    yerr=summary_df['std'], marker='o', capsize=5, linewidth=2)
        ax3.set_xticks(range(len(summary_df)))
        ax3.set_xticklabels(summary_df['experiment'], rotation=45, ha='right')
        ax3.set_ylabel('Mean Accessibility (kcal/mol)')
        ax3.set_title('Performance Trend')
        ax3.grid(True, alpha=0.3)
        

        ax4 = axes[1, 1]
        ax4.axis('tight')
        ax4.axis('off')
        
        table_data = []
        for exp_name, _ in experiments:
            exp_data = combined_df[combined_df['experiment'] == exp_name]
            table_data.append([
                exp_name[:15],
                f"{exp_data['accessibility'].mean():.3f}",
                f"{exp_data['accessibility'].median():.3f}",
                f"{exp_data['accessibility'].min():.3f}",
                f"{(exp_data['accessibility'] < 1.0).sum() / len(exp_data) * 100:.1f}%"
            ])
        
        table = ax4.table(cellText=table_data,
                         colLabels=['Experiment', 'Mean', 'Median', 'Best', 'Excellent%'],
                         cellLoc='center',
                         loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1.2, 1.5)
        
        plt.suptitle('Multi-Experiment Comparison', fontsize=14, y=1.02)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig


def main():

    from pathlib import Path
    

    visualizer = PerformanceVisualizer()
    

    exp_dir = Path('results_full/20250909_105804_unified_access_experiments')
    if exp_dir.exists():
        df = visualizer.load_experiment_data(exp_dir)
        

        visualizer.plot_constraint_comparison(df)
        visualizer.plot_performance_distribution(df)
        visualizer.plot_protein_performance(df)
        

        report = visualizer.generate_performance_report(df)





if __name__ == "__main__":
    main()