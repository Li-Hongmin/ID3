"""
Uniqueness Visualizer Module

Creates visualizations for sequence uniqueness analysis.
"""

import logging
from pathlib import Path
from typing import Dict, Any, Optional
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

# Set matplotlib style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")


class UniquenessVisualizer:
    """
    Create visualizations for sequence uniqueness analysis.
    """
    
    def __init__(self, output_dir: Path):
        """
        Initialize visualizer.
        
        Args:
            output_dir: Directory to save visualizations
        """
        self.output_dir = Path(output_dir)
        self.figures_dir = self.output_dir / 'figures'
        self.figures_dir.mkdir(parents=True, exist_ok=True)
    
    def create_all_visualizations(self, df: pd.DataFrame, summary: dict, patterns: dict) -> Dict[str, Path]:
        """
        Create all visualizations.
        
        Args:
            df: DataFrame with analysis results
            summary: Summary statistics
            patterns: Identified patterns
            
        Returns:
            Dictionary mapping visualization names to file paths
        """
        files = {}
        
        if df.empty:
            logger.warning("Empty dataframe, skipping visualizations")
            return files
        
        # Create individual visualizations
        try:
            files['distribution'] = self.plot_uniqueness_distribution(df)
        except Exception as e:
            logger.warning(f"Failed to create distribution plot: {e}")
        
        try:
            files['by_constraint'] = self.plot_by_constraint(df)
        except Exception as e:
            logger.warning(f"Failed to create constraint plot: {e}")
        
        try:
            files['by_variant'] = self.plot_by_variant(df)
        except Exception as e:
            logger.warning(f"Failed to create variant plot: {e}")
        
        try:
            files['heatmap'] = self.plot_configuration_heatmap(df)
        except Exception as e:
            logger.warning(f"Failed to create heatmap: {e}")
        
        try:
            files['correlation'] = self.plot_correlation_scatter(df)
        except Exception as e:
            logger.warning(f"Failed to create correlation plot: {e}")
        
        try:
            files['temporal'] = self.plot_temporal_analysis(df)
        except Exception as e:
            logger.warning(f"Failed to create temporal plot: {e}")
        
        try:
            files['comprehensive'] = self.create_comprehensive_figure(df, summary, patterns)
        except Exception as e:
            logger.warning(f"Failed to create comprehensive figure: {e}")
        
        return files
    
    def plot_uniqueness_distribution(self, df: pd.DataFrame) -> Path:
        """
        Plot distribution of uniqueness rates.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Histogram
        axes[0].hist(df['uniqueness_rate'] * 100, bins=30, edgecolor='black', alpha=0.7)
        axes[0].set_xlabel('Uniqueness Rate (%)')
        axes[0].set_ylabel('Count')
        axes[0].set_title('Distribution of Sequence Uniqueness Rates')
        axes[0].axvline(df['uniqueness_rate'].mean() * 100, color='red', 
                       linestyle='--', label=f'Mean: {df["uniqueness_rate"].mean():.1%}')
        axes[0].legend()
        
        # Violin plot by constraint type
        if 'constraint_type' in df.columns:
            df_plot = df.copy()
            df_plot['uniqueness_percent'] = df_plot['uniqueness_rate'] * 100
            sns.violinplot(data=df_plot, x='constraint_type', y='uniqueness_percent', ax=axes[1])
            axes[1].set_xlabel('Constraint Type')
            axes[1].set_ylabel('Uniqueness Rate (%)')
            axes[1].set_title('Uniqueness Distribution by Constraint Type')
        
        plt.tight_layout()
        filepath = self.figures_dir / 'uniqueness_distribution.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def plot_by_constraint(self, df: pd.DataFrame) -> Path:
        """
        Plot uniqueness by constraint type.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # Box plot
        df_plot = df.copy()
        df_plot['uniqueness_percent'] = df_plot['uniqueness_rate'] * 100
        
        sns.boxplot(data=df_plot, x='constraint_type', y='uniqueness_percent', ax=axes[0, 0])
        axes[0, 0].set_title('Uniqueness Rate by Constraint Type')
        axes[0, 0].set_ylabel('Uniqueness Rate (%)')
        
        # Shannon entropy comparison
        if 'shannon_entropy' in df.columns:
            sns.boxplot(data=df_plot, x='constraint_type', y='shannon_entropy', ax=axes[0, 1])
            axes[0, 1].set_title('Shannon Entropy by Constraint Type')
            axes[0, 1].set_ylabel('Shannon Entropy')
        
        # Repetition factor
        if 'repetition_factor' in df.columns:
            sns.boxplot(data=df_plot, x='constraint_type', y='repetition_factor', ax=axes[1, 0])
            axes[1, 0].set_title('Repetition Factor by Constraint Type')
            axes[1, 0].set_ylabel('Repetition Factor')
            axes[1, 0].set_ylim(bottom=0)
        
        # Mean comparison with error bars
        means = df_plot.groupby('constraint_type')['uniqueness_percent'].agg(['mean', 'std', 'count'])
        x_pos = np.arange(len(means))
        axes[1, 1].bar(x_pos, means['mean'], yerr=means['std'], capsize=5, alpha=0.7)
        axes[1, 1].set_xticks(x_pos)
        axes[1, 1].set_xticklabels(means.index)
        axes[1, 1].set_xlabel('Constraint Type')
        axes[1, 1].set_ylabel('Mean Uniqueness Rate (%)')
        axes[1, 1].set_title('Mean Uniqueness Rate Comparison')
        
        # Add sample sizes
        for i, (idx, row) in enumerate(means.iterrows()):
            axes[1, 1].text(i, row['mean'] + row['std'] + 0.5, f'n={row["count"]}', 
                          ha='center', va='bottom')
        
        plt.tight_layout()
        filepath = self.figures_dir / 'uniqueness_by_constraint.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def plot_by_variant(self, df: pd.DataFrame) -> Path:
        """
        Plot uniqueness by variant.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        df_plot = df.copy()
        df_plot['uniqueness_percent'] = df_plot['uniqueness_rate'] * 100
        
        # Box plot by variant
        sns.boxplot(data=df_plot, x='variant', y='uniqueness_percent', ax=axes[0])
        axes[0].set_title('Uniqueness Rate by Variant')
        axes[0].set_xlabel('Variant')
        axes[0].set_ylabel('Uniqueness Rate (%)')
        
        # Grouped by constraint and variant
        if 'constraint_type' in df.columns:
            sns.boxplot(data=df_plot, x='variant', y='uniqueness_percent', 
                       hue='constraint_type', ax=axes[1])
            axes[1].set_title('Uniqueness by Variant and Constraint Type')
            axes[1].set_xlabel('Variant')
            axes[1].set_ylabel('Uniqueness Rate (%)')
            axes[1].legend(title='Constraint', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.tight_layout()
        filepath = self.figures_dir / 'uniqueness_by_variant.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def plot_configuration_heatmap(self, df: pd.DataFrame) -> Path:
        """
        Plot heatmap of uniqueness by configuration.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        # Create pivot table
        pivot = df.pivot_table(
            values='uniqueness_rate',
            index='protein',
            columns=['constraint_type', 'variant'],
            aggfunc='mean'
        )
        
        # Convert to percentage
        pivot = pivot * 100
        
        # Create figure
        fig, ax = plt.subplots(figsize=(16, 10))
        
        sns.heatmap(pivot, annot=True, fmt='.1f', cmap='RdYlGn', 
                   cbar_kws={'label': 'Uniqueness Rate (%)'},
                   vmin=0, vmax=pivot.max().max(), ax=ax)
        
        ax.set_title('Sequence Uniqueness Rate Heatmap\nby Protein, Constraint Type, and Variant')
        ax.set_xlabel('Constraint Type / Variant')
        ax.set_ylabel('Protein')
        
        plt.tight_layout()
        filepath = self.figures_dir / 'uniqueness_heatmap.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def plot_correlation_scatter(self, df: pd.DataFrame) -> Path:
        """
        Plot correlation between uniqueness and performance.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        df_plot = df.copy()
        df_plot['uniqueness_percent'] = df_plot['uniqueness_rate'] * 100
        
        # Uniqueness vs Best Accessibility
        if 'best_accessibility' in df.columns:
            axes[0, 0].scatter(df_plot['uniqueness_percent'], df_plot['best_accessibility'], alpha=0.6)
            axes[0, 0].set_xlabel('Uniqueness Rate (%)')
            axes[0, 0].set_ylabel('Best Accessibility (kcal/mol)')
            axes[0, 0].set_title('Uniqueness vs Best Accessibility')
            
            # Add trend line
            z = np.polyfit(df_plot['uniqueness_percent'], df_plot['best_accessibility'], 1)
            p = np.poly1d(z)
            axes[0, 0].plot(df_plot['uniqueness_percent'].sort_values(), 
                          p(df_plot['uniqueness_percent'].sort_values()),
                          "r--", alpha=0.5, label=f'Trend (r={df_plot["uniqueness_percent"].corr(df_plot["best_accessibility"]):.3f})')
            axes[0, 0].legend()
        
        # Uniqueness vs Sequence Length
        if 'sequence_length' in df.columns:
            axes[0, 1].scatter(df_plot['sequence_length'], df_plot['uniqueness_percent'], alpha=0.6)
            axes[0, 1].set_xlabel('Sequence Length (nt)')
            axes[0, 1].set_ylabel('Uniqueness Rate (%)')
            axes[0, 1].set_title('Sequence Length vs Uniqueness')
        
        # Shannon Entropy vs Simpson Diversity
        if 'shannon_entropy' in df.columns and 'simpson_diversity' in df.columns:
            axes[1, 0].scatter(df_plot['shannon_entropy'], df_plot['simpson_diversity'], alpha=0.6)
            axes[1, 0].set_xlabel('Shannon Entropy')
            axes[1, 0].set_ylabel('Simpson Diversity Index')
            axes[1, 0].set_title('Entropy vs Diversity Index')
        
        # Repetition Factor Distribution
        if 'repetition_factor' in df.columns:
            axes[1, 1].hist(df_plot['repetition_factor'], bins=30, edgecolor='black', alpha=0.7)
            axes[1, 1].set_xlabel('Repetition Factor')
            axes[1, 1].set_ylabel('Count')
            axes[1, 1].set_title('Distribution of Repetition Factors')
            axes[1, 1].axvline(df_plot['repetition_factor'].mean(), color='red', 
                             linestyle='--', label=f'Mean: {df_plot["repetition_factor"].mean():.1f}')
            axes[1, 1].legend()
        
        plt.tight_layout()
        filepath = self.figures_dir / 'uniqueness_correlations.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def plot_temporal_analysis(self, df: pd.DataFrame) -> Path:
        """
        Plot temporal analysis of uniqueness.
        
        Args:
            df: DataFrame with results
            
        Returns:
            Path to saved figure
        """
        # Check for temporal columns
        temporal_cols = [col for col in df.columns if col.startswith('uniqueness_at_')]
        
        if not temporal_cols:
            logger.warning("No temporal data available")
            return None
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Extract iteration numbers
        iterations = []
        mean_values = []
        std_values = []
        
        for col in sorted(temporal_cols):
            iter_num = int(col.split('_')[-1])
            iterations.append(iter_num)
            mean_values.append(df[col].mean() * 100)
            std_values.append(df[col].std() * 100)
        
        # Plot mean uniqueness over iterations
        axes[0].errorbar(iterations, mean_values, yerr=std_values, 
                        marker='o', capsize=5, capthick=2)
        axes[0].set_xlabel('Iteration')
        axes[0].set_ylabel('Mean Uniqueness Rate (%)')
        axes[0].set_title('Uniqueness Rate Evolution During Optimization')
        axes[0].grid(True, alpha=0.3)
        
        # Plot by constraint type
        for constraint in df['constraint_type'].unique():
            constraint_df = df[df['constraint_type'] == constraint]
            means = [constraint_df[col].mean() * 100 for col in sorted(temporal_cols)]
            axes[1].plot(iterations, means, marker='o', label=constraint)
        
        axes[1].set_xlabel('Iteration')
        axes[1].set_ylabel('Mean Uniqueness Rate (%)')
        axes[1].set_title('Uniqueness Evolution by Constraint Type')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        filepath = self.figures_dir / 'temporal_uniqueness.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath
    
    def create_comprehensive_figure(self, df: pd.DataFrame, summary: dict, patterns: dict) -> Path:
        """
        Create a comprehensive summary figure.
        
        Args:
            df: DataFrame with results
            summary: Summary statistics
            patterns: Identified patterns
            
        Returns:
            Path to saved figure
        """
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)
        
        df_plot = df.copy()
        df_plot['uniqueness_percent'] = df_plot['uniqueness_rate'] * 100
        
        # 1. Overall distribution
        ax1 = fig.add_subplot(gs[0, :2])
        ax1.hist(df_plot['uniqueness_percent'], bins=30, edgecolor='black', alpha=0.7, color='skyblue')
        ax1.axvline(df_plot['uniqueness_percent'].mean(), color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {df_plot["uniqueness_percent"].mean():.1f}%')
        ax1.set_xlabel('Uniqueness Rate (%)', fontsize=12)
        ax1.set_ylabel('Count', fontsize=12)
        ax1.set_title('Distribution of Sequence Uniqueness Rates', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Summary statistics box
        ax2 = fig.add_subplot(gs[0, 2])
        ax2.axis('off')
        summary_text = f"""
Global Statistics
─────────────────
Mean: {summary.get('global_uniqueness_rate', 0):.1%}
Std: {summary.get('global_uniqueness_std', 0):.1%}
Min: {summary.get('min_uniqueness_rate', 0):.1%}
Max: {summary.get('max_uniqueness_rate', 0):.1%}

Total Experiments: {summary.get('total_experiments', 0)}
Total Sequences: {summary.get('total_sequences', 0):,}
Unique Sequences: {summary.get('total_unique_sequences', 0):,}

Mean Entropy: {summary.get('mean_shannon_entropy', 0):.2f}
Mean Diversity: {summary.get('mean_simpson_diversity', 0):.3f}
        """
        ax2.text(0.1, 0.9, summary_text, transform=ax2.transAxes, 
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # 3. By constraint type
        ax3 = fig.add_subplot(gs[1, 0])
        sns.boxplot(data=df_plot, x='constraint_type', y='uniqueness_percent', ax=ax3)
        ax3.set_xlabel('Constraint Type', fontsize=12)
        ax3.set_ylabel('Uniqueness Rate (%)', fontsize=12)
        ax3.set_title('By Constraint Type', fontsize=12, fontweight='bold')
        
        # 4. By variant
        ax4 = fig.add_subplot(gs[1, 1])
        sns.boxplot(data=df_plot, x='variant', y='uniqueness_percent', ax=ax4)
        ax4.set_xlabel('Variant', fontsize=12)
        ax4.set_ylabel('Uniqueness Rate (%)', fontsize=12)
        ax4.set_title('By Variant', fontsize=12, fontweight='bold')
        
        # 5. Correlation with performance
        ax5 = fig.add_subplot(gs[1, 2])
        if 'best_accessibility' in df.columns:
            ax5.scatter(df_plot['uniqueness_percent'], df_plot['best_accessibility'], 
                       alpha=0.5, s=20)
            ax5.set_xlabel('Uniqueness Rate (%)', fontsize=12)
            ax5.set_ylabel('Best Accessibility', fontsize=12)
            ax5.set_title('Uniqueness vs Performance', fontsize=12, fontweight='bold')
            
            # Add correlation coefficient
            corr = df_plot['uniqueness_percent'].corr(df_plot['best_accessibility'])
            ax5.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax5.transAxes,
                    fontsize=10, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # 6. Heatmap
        ax6 = fig.add_subplot(gs[2:, :])
        pivot = df.pivot_table(
            values='uniqueness_rate',
            index='protein',
            columns=['constraint_type', 'variant'],
            aggfunc='mean'
        ) * 100
        
        sns.heatmap(pivot, annot=True, fmt='.1f', cmap='RdYlGn', 
                   cbar_kws={'label': 'Uniqueness Rate (%)'},
                   ax=ax6, vmin=0)
        ax6.set_title('Uniqueness Rate Heatmap by Configuration', fontsize=14, fontweight='bold')
        ax6.set_xlabel('Constraint Type / Variant', fontsize=12)
        ax6.set_ylabel('Protein', fontsize=12)
        
        # Overall title
        fig.suptitle('Sequence Uniqueness Analysis - Comprehensive Summary', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        plt.tight_layout()
        filepath = self.figures_dir / 'comprehensive_uniqueness_summary.png'
        plt.savefig(filepath, dpi=150, bbox_inches='tight')
        plt.close()
        
        return filepath