#!/usr/bin/env python3
"""
Visualization Module for Unified Experiment Analysis
Creates various plots and figures for experimental results
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import logging
from collections import defaultdict

# Set style for publication-ready figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ExperimentVisualizer:
    """Visualizer for unified experiment results"""
    
    def __init__(self, output_dir: str = None):
        """
        Initialize visualizer
        
        Args:
            output_dir: Directory to save figures
        """
        if output_dir:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.output_dir = Path.cwd()
            
        # Set default figure properties
        self.figure_dpi = 300
        self.figure_format = 'png'
        
    def plot_convergence_curves(self, results: List[Dict], protein: str = None, 
                              save_path: str = None) -> plt.Figure:
        """
        Plot convergence curves for optimization
        
        Args:
            results: List of experiment results
            protein: Specific protein to plot (if None, plot all)
            save_path: Path to save figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        axes = axes.flatten()
        
        # Group by constraint type
        constraint_groups = defaultdict(list)
        for result in results:
            if protein and result['protein_name'] != protein:
                continue
            constraint_groups[result['constraint_type']].append(result)
        
        # Plot for each constraint type
        for idx, (constraint, group) in enumerate(constraint_groups.items()):
            if idx >= 4:
                break
                
            ax = axes[idx]
            
            for result in group:
                trajectory = result.get('trajectory', {})
                iterations = trajectory.get('iterations', [])
                accessibility = trajectory.get('accessibility', [])
                
                if iterations and accessibility:
                    label = f"{result['protein_name']}_{result['variant']}"
                    ax.plot(iterations, accessibility, label=label, alpha=0.7)
            
            ax.set_xlabel('Iteration')
            ax.set_ylabel('Accessibility Score')
            ax.set_title(f'{constraint.upper()} Convergence')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(len(constraint_groups), 4):
            axes[idx].set_visible(False)
        
        fig.suptitle(f'Optimization Convergence Curves{" - " + protein if protein else ""}', 
                    fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved convergence plot to {save_path}")
        
        return fig
    
    def plot_performance_comparison(self, results: List[Dict], 
                                  save_path: str = None) -> plt.Figure:
        """
        Plot performance comparison across constraint types and variants
        
        Args:
            results: List of experiment results
            save_path: Path to save figure
        """
        # Prepare data
        constraints = ['ams', 'cpc', 'lagrangian']
        variants = ['00', '01', '10', '11']
        
        # Calculate average final accessibility for each combination
        performance_data = defaultdict(list)
        
        for constraint in constraints:
            for variant in variants:
                values = []
                for result in results:
                    if (result['constraint_type'] == constraint and 
                        result['variant'] == variant):
                        values.append(result.get('final_accessibility', np.nan))
                
                if values:
                    performance_data[constraint].append(np.nanmean(values))
                else:
                    performance_data[constraint].append(np.nan)
        
        # Create plot
        fig, ax = plt.subplots(figsize=(12, 6))
        
        x = np.arange(len(variants))
        width = 0.25
        
        colors = {'ams': '#1f77b4', 'cpc': '#ff7f0e', 'lagrangian': '#2ca02c'}
        
        for i, constraint in enumerate(constraints):
            values = performance_data[constraint]
            ax.bar(x + i * width, values, width, label=constraint.upper(), 
                  color=colors[constraint], alpha=0.8)
        
        ax.set_xlabel('Variant', fontsize=12)
        ax.set_ylabel('Average Final Accessibility', fontsize=12)
        ax.set_title('Performance Comparison Across Constraints and Variants', 
                    fontsize=14, fontweight='bold')
        ax.set_xticks(x + width)
        ax.set_xticklabels(variants)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved performance comparison to {save_path}")
        
        return fig
    
    def plot_cai_effect_analysis(self, results: List[Dict], 
                                save_path: str = None) -> plt.Figure:
        """
        Plot CAI optimization effect analysis
        
        Args:
            results: List of experiment results with CAI enabled
            save_path: Path to save figure
        """
        # Filter CAI-enabled results
        cai_results = [r for r in results if r.get('cai_enabled', False)]
        
        if not cai_results:
            logger.warning("No CAI-enabled results found")
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # 1. CAI improvement over iterations (top-left)
        ax1 = axes[0, 0]
        for result in cai_results[:5]:  # Limit to 5 for clarity
            trajectory = result.get('trajectory', {})
            iterations = trajectory.get('iterations', [])
            ecai_values = trajectory.get('ecai_values', [])
            
            if iterations and ecai_values:
                label = f"{result['protein_name']}_{result['constraint_type']}_{result['variant']}"
                ax1.plot(iterations, ecai_values, label=label, alpha=0.7)
        
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('ECAI Value')
        ax1.set_title('CAI Optimization Progress')
        ax1.legend(fontsize=8, loc='best')
        ax1.grid(True, alpha=0.3)
        
        # 2. Initial vs Final CAI comparison (top-right)
        ax2 = axes[0, 1]
        initial_cais = []
        final_cais = []
        
        for result in cai_results:
            trajectory = result.get('trajectory', {})
            ecai_values = trajectory.get('ecai_values', [])
            if ecai_values:
                initial_cais.append(ecai_values[0])
                final_cais.append(ecai_values[-1])
        
        if initial_cais and final_cais:
            ax2.scatter(initial_cais, final_cais, alpha=0.6)
            ax2.plot([min(initial_cais), max(initial_cais)], 
                    [min(initial_cais), max(initial_cais)], 
                    'r--', alpha=0.5, label='No change line')
            ax2.set_xlabel('Initial ECAI')
            ax2.set_ylabel('Final ECAI')
            ax2.set_title('CAI Improvement Scatter Plot')
            ax2.legend()
            ax2.grid(True, alpha=0.3)
        
        # 3. CAI vs Accessibility trade-off (bottom-left)
        ax3 = axes[1, 0]
        final_accessibilities = []
        final_ecais = []
        
        for result in cai_results:
            final_acc = result.get('final_accessibility')
            trajectory = result.get('trajectory', {})
            ecai_values = trajectory.get('ecai_values', [])
            
            if final_acc and ecai_values:
                final_accessibilities.append(final_acc)
                final_ecais.append(ecai_values[-1])
        
        if final_accessibilities and final_ecais:
            ax3.scatter(final_ecais, final_accessibilities, alpha=0.6)
            ax3.set_xlabel('Final ECAI')
            ax3.set_ylabel('Final Accessibility')
            ax3.set_title('CAI vs Accessibility Trade-off')
            ax3.grid(True, alpha=0.3)
        
        # 4. CAI distribution by constraint type (bottom-right)
        ax4 = axes[1, 1]
        cai_by_constraint = defaultdict(list)
        
        for result in cai_results:
            constraint = result['constraint_type']
            trajectory = result.get('trajectory', {})
            ecai_values = trajectory.get('ecai_values', [])
            if ecai_values:
                cai_by_constraint[constraint].append(ecai_values[-1])
        
        if cai_by_constraint:
            data_to_plot = [cai_by_constraint[c] for c in ['ams', 'cpc', 'lagrangian'] 
                          if c in cai_by_constraint]
            labels = [c.upper() for c in ['ams', 'cpc', 'lagrangian'] 
                     if c in cai_by_constraint]
            
            bp = ax4.boxplot(data_to_plot, labels=labels)
            ax4.set_ylabel('Final ECAI')
            ax4.set_title('CAI Distribution by Constraint Type')
            ax4.grid(True, alpha=0.3, axis='y')
        
        fig.suptitle('CAI Optimization Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved CAI analysis to {save_path}")
        
        return fig
    
    def plot_protein_comparison(self, results: List[Dict], 
                              save_path: str = None) -> plt.Figure:
        """
        Plot comparison across different proteins
        
        Args:
            results: List of experiment results
            save_path: Path to save figure
        """
        # Group results by protein
        protein_groups = defaultdict(list)
        for result in results:
            protein_groups[result['protein_name']].append(result)
        
        # Calculate statistics for each protein
        protein_stats = {}
        for protein, group in protein_groups.items():
            improvements = [r['improvement'] for r in group if r.get('improvement')]
            final_accs = [r['final_accessibility'] for r in group if r.get('final_accessibility')]
            
            protein_stats[protein] = {
                'avg_improvement': np.mean(improvements) if improvements else 0,
                'std_improvement': np.std(improvements) if improvements else 0,
                'avg_final': np.mean(final_accs) if final_accs else 0,
                'std_final': np.std(final_accs) if final_accs else 0,
                'best_final': min(final_accs) if final_accs else np.inf,
                'count': len(group)
            }
        
        # Create plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Sort proteins by average final accessibility
        sorted_proteins = sorted(protein_stats.keys(), 
                               key=lambda x: protein_stats[x]['avg_final'])
        
        # 1. Average improvement by protein
        ax1 = axes[0]
        improvements = [protein_stats[p]['avg_improvement'] for p in sorted_proteins]
        errors = [protein_stats[p]['std_improvement'] for p in sorted_proteins]
        
        bars1 = ax1.bar(range(len(sorted_proteins)), improvements, yerr=errors, 
                       capsize=5, alpha=0.7, color='skyblue')
        ax1.set_xlabel('Protein')
        ax1.set_ylabel('Average Improvement')
        ax1.set_title('Average Accessibility Improvement by Protein')
        ax1.set_xticks(range(len(sorted_proteins)))
        ax1.set_xticklabels(sorted_proteins, rotation=45, ha='right')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # 2. Final accessibility by protein
        ax2 = axes[1]
        finals = [protein_stats[p]['avg_final'] for p in sorted_proteins]
        errors = [protein_stats[p]['std_final'] for p in sorted_proteins]
        bests = [protein_stats[p]['best_final'] for p in sorted_proteins]
        
        bars2 = ax2.bar(range(len(sorted_proteins)), finals, yerr=errors, 
                       capsize=5, alpha=0.7, color='lightcoral')
        ax2.scatter(range(len(sorted_proteins)), bests, color='red', 
                   marker='*', s=100, zorder=5, label='Best result')
        ax2.set_xlabel('Protein')
        ax2.set_ylabel('Final Accessibility')
        ax2.set_title('Final Accessibility by Protein')
        ax2.set_xticks(range(len(sorted_proteins)))
        ax2.set_xticklabels(sorted_proteins, rotation=45, ha='right')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.suptitle('Protein Performance Comparison', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved protein comparison to {save_path}")
        
        return fig
    
    def plot_performance_heatmap(self, matrix: np.ndarray, proteins: List[str], 
                                combinations: List[str], save_path: str = None) -> plt.Figure:
        """
        Plot performance heatmap for all experiments
        
        Args:
            matrix: Performance matrix (proteins x constraint_variant)
            proteins: List of protein names
            combinations: List of constraint_variant combinations
            save_path: Path to save figure
        """
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Create heatmap
        im = ax.imshow(matrix, cmap='RdYlGn_r', aspect='auto')
        
        # Set ticks and labels
        ax.set_xticks(np.arange(len(combinations)))
        ax.set_yticks(np.arange(len(proteins)))
        ax.set_xticklabels(combinations, rotation=45, ha='right')
        ax.set_yticklabels(proteins)
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Final Accessibility', rotation=270, labelpad=20)
        
        # Add text annotations
        for i in range(len(proteins)):
            for j in range(len(combinations)):
                if not np.isnan(matrix[i, j]):
                    text = ax.text(j, i, f'{matrix[i, j]:.2f}',
                                 ha="center", va="center", color="black", fontsize=8)
        
        ax.set_xlabel('Constraint_Variant')
        ax.set_ylabel('Protein')
        ax.set_title('Performance Heatmap: Final Accessibility Across All Experiments', 
                    fontsize=14, fontweight='bold')
        
        # Add grid
        ax.set_xticks(np.arange(len(combinations) + 1) - 0.5, minor=True)
        ax.set_yticks(np.arange(len(proteins) + 1) - 0.5, minor=True)
        ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved performance heatmap to {save_path}")
        
        return fig
    
    def plot_constraint_variant_analysis(self, results: List[Dict], 
                                        save_path: str = None) -> plt.Figure:
        """
        Detailed analysis of constraint types and variants
        
        Args:
            results: List of experiment results
            save_path: Path to save figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Prepare data
        constraint_variant_data = defaultdict(lambda: defaultdict(list))
        
        for result in results:
            constraint = result['constraint_type']
            variant = result['variant']
            final_acc = result.get('final_accessibility')
            if final_acc:
                constraint_variant_data[constraint][variant].append(final_acc)
        
        # 1. Box plot by constraint type (top-left)
        ax1 = axes[0, 0]
        constraint_data = []
        constraint_labels = []
        
        for constraint in ['ams', 'cpc', 'lagrangian']:
            all_values = []
            for variant_values in constraint_variant_data[constraint].values():
                all_values.extend(variant_values)
            if all_values:
                constraint_data.append(all_values)
                constraint_labels.append(constraint.upper())
        
        if constraint_data:
            bp1 = ax1.boxplot(constraint_data, labels=constraint_labels)
            ax1.set_ylabel('Final Accessibility')
            ax1.set_title('Performance Distribution by Constraint Type')
            ax1.grid(True, alpha=0.3, axis='y')
        
        # 2. Box plot by variant (top-right)
        ax2 = axes[0, 1]
        variant_data = []
        variant_labels = []
        
        for variant in ['00', '01', '10', '11']:
            all_values = []
            for constraint in constraint_variant_data:
                all_values.extend(constraint_variant_data[constraint][variant])
            if all_values:
                variant_data.append(all_values)
                variant_labels.append(variant)
        
        if variant_data:
            bp2 = ax2.boxplot(variant_data, labels=variant_labels)
            ax2.set_ylabel('Final Accessibility')
            ax2.set_title('Performance Distribution by Variant')
            ax2.grid(True, alpha=0.3, axis='y')
        
        # 3. Interaction plot (bottom-left)
        ax3 = axes[1, 0]
        variants = ['00', '01', '10', '11']
        colors = {'ams': 'blue', 'cpc': 'orange', 'lagrangian': 'green'}
        
        for constraint in ['ams', 'cpc', 'lagrangian']:
            means = []
            for variant in variants:
                values = constraint_variant_data[constraint][variant]
                means.append(np.mean(values) if values else np.nan)
            
            ax3.plot(variants, means, marker='o', label=constraint.upper(), 
                    color=colors[constraint], linewidth=2)
        
        ax3.set_xlabel('Variant')
        ax3.set_ylabel('Mean Final Accessibility')
        ax3.set_title('Constraint-Variant Interaction')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Best performers table (bottom-right)
        ax4 = axes[1, 1]
        ax4.axis('tight')
        ax4.axis('off')
        
        # Find best performer for each constraint
        best_performers = []
        for constraint in ['ams', 'cpc', 'lagrangian']:
            best_variant = None
            best_score = float('inf')
            
            for variant in variants:
                values = constraint_variant_data[constraint][variant]
                if values:
                    mean_score = np.mean(values)
                    if mean_score < best_score:
                        best_score = mean_score
                        best_variant = variant
            
            if best_variant:
                best_performers.append([constraint.upper(), best_variant, f'{best_score:.3f}'])
        
        if best_performers:
            table = ax4.table(cellText=best_performers,
                            colLabels=['Constraint', 'Best Variant', 'Avg Score'],
                            cellLoc='center',
                            loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.2, 1.5)
            ax4.set_title('Best Performing Variants', fontweight='bold', pad=20)
        
        fig.suptitle('Constraint Type and Variant Analysis', fontsize=16, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved constraint-variant analysis to {save_path}")
        
        return fig
    
    def plot_constraint_mechanism_boxplot(self, df_access: pd.DataFrame, 
                                         df_penalty: pd.DataFrame = None,
                                         df_no_penalty: pd.DataFrame = None,
                                         save_path: str = None) -> plt.Figure:
        """
        Plot boxplots grouped by constraint mechanism (L, A, C)
        
        Args:
            df_access: Access-only results DataFrame
            df_penalty: With penalty results DataFrame (optional)
            df_no_penalty: No penalty results DataFrame (optional)
            save_path: Path to save figure
        """
        fig, axes = plt.subplots(1, 3, figsize=(15, 6))
        
        # Color scheme
        colors = {
            'Access-only': '#1f77b4',
            'With Penalty': '#ff7f0e', 
            'No Penalty': '#2ca02c'
        }
        
        mechanisms = ['L', 'A', 'C']
        mechanism_names = {
            'L': 'Lagrangian',
            'A': 'Amino Matching Softmax',
            'C': 'Codon Profile Constraint'
        }
        
        for idx, mechanism in enumerate(mechanisms):
            ax = axes[idx]
            data_to_plot = []
            labels = []
            colors_list = []
            
            # Filter variants for this mechanism
            variants = [v for v in df_access.index if mechanism in v]
            
            # Access-only data
            if df_access is not None and len(variants) > 0:
                access_values = []
                for variant in variants:
                    if variant in df_access.index:
                        # Get mean value across proteins
                        row_values = df_access.loc[variant].values
                        # Filter out non-numeric values
                        numeric_values = [float(v) for v in row_values if isinstance(v, (int, float)) or (isinstance(v, str) and v.replace('.','').isdigit())]
                        if numeric_values:
                            access_values.extend(numeric_values)
                
                if access_values:
                    data_to_plot.append(access_values)
                    labels.append('Access-only')
                    colors_list.append(colors['Access-only'])
            
            # With Penalty data
            if df_penalty is not None:
                penalty_values = []
                for variant in variants:
                    if variant in df_penalty.index:
                        row_values = df_penalty.loc[variant].values
                        # Extract accessibility values from "Access (CAI)" format
                        for val in row_values:
                            if isinstance(val, str) and '(' in val:
                                try:
                                    access_val = float(val.split('(')[0].strip())
                                    penalty_values.append(access_val)
                                except:
                                    pass
                
                if penalty_values:
                    data_to_plot.append(penalty_values)
                    labels.append('With Penalty')
                    colors_list.append(colors['With Penalty'])
            
            # No Penalty data
            if df_no_penalty is not None:
                no_penalty_values = []
                for variant in variants:
                    if variant in df_no_penalty.index:
                        row_values = df_no_penalty.loc[variant].values
                        # Extract accessibility values
                        for val in row_values:
                            if isinstance(val, str) and '(' in val:
                                try:
                                    access_val = float(val.split('(')[0].strip())
                                    no_penalty_values.append(access_val)
                                except:
                                    pass
                
                if no_penalty_values:
                    data_to_plot.append(no_penalty_values)
                    labels.append('No Penalty')
                    colors_list.append(colors['No Penalty'])
            
            # Create boxplot
            if data_to_plot:
                bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True)
                
                # Color the boxes
                for patch, color in zip(bp['boxes'], colors_list):
                    patch.set_facecolor(color)
                    patch.set_alpha(0.7)
                
                # Customize boxplot appearance
                for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
                    if element in bp:
                        plt.setp(bp[element], color='black')
            
            ax.set_title(mechanism_names[mechanism], fontsize=12, fontweight='bold')
            ax.set_ylabel('Accessibility (kcal/mol)' if idx == 0 else '')
            ax.grid(True, alpha=0.3, axis='y')
            
            # Add variant count info
            ax.text(0.02, 0.98, f'{len(variants)} variants',
                   transform=ax.transAxes, fontsize=9,
                   verticalalignment='top')
        
        fig.suptitle('Constraint Mechanism Performance Distribution', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved constraint mechanism boxplot to {save_path}")
        
        return fig
    
    def plot_accessibility_heatmap(self, df: pd.DataFrame, title: str = "Accessibility Heatmap",
                                  save_path: str = None) -> plt.Figure:
        """
        Plot heatmap of accessibility values
        
        Args:
            df: DataFrame with variants as rows and proteins as columns
            title: Title for the heatmap
            save_path: Path to save figure
        """
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Prepare data - extract numeric values only
        numeric_df = pd.DataFrame()
        for col in df.columns:
            if col not in ['Mean', 'Rank']:  # Skip summary columns
                numeric_col = []
                for val in df[col]:
                    if isinstance(val, (int, float)):
                        numeric_col.append(float(val))
                    elif isinstance(val, str):
                        # Try to extract number from "Access (CAI)" format
                        if '(' in val:
                            try:
                                numeric_col.append(float(val.split('(')[0].strip()))
                            except:
                                numeric_col.append(np.nan)
                        else:
                            try:
                                numeric_col.append(float(val))
                            except:
                                numeric_col.append(np.nan)
                    else:
                        numeric_col.append(np.nan)
                
                if numeric_col and not all(np.isnan(numeric_col)):
                    numeric_df[col] = numeric_col
        
        if not numeric_df.empty:
            numeric_df.index = df.index
            
            # Create heatmap
            im = ax.imshow(numeric_df.values, cmap='RdYlGn_r', aspect='auto')
            
            # Set ticks
            ax.set_xticks(np.arange(len(numeric_df.columns)))
            ax.set_yticks(np.arange(len(numeric_df.index)))
            ax.set_xticklabels(numeric_df.columns, rotation=45, ha='right')
            ax.set_yticklabels(numeric_df.index)
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Accessibility (kcal/mol)', rotation=270, labelpad=20)
            
            # Add text annotations
            for i in range(len(numeric_df.index)):
                for j in range(len(numeric_df.columns)):
                    val = numeric_df.iloc[i, j]
                    if not np.isnan(val):
                        text = ax.text(j, i, f'{val:.2f}',
                                     ha="center", va="center",
                                     color="white" if val > np.nanmean(numeric_df.values) else "black",
                                     fontsize=8)
        
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Protein', fontsize=12)
        ax.set_ylabel('Variant', fontsize=12)
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved heatmap to {save_path}")
        
        return fig
    
    def plot_ranking_change(self, df_comparison: pd.DataFrame, save_path: str = None) -> plt.Figure:
        """
        Plot ranking change between with penalty and no penalty
        
        Args:
            df_comparison: DataFrame with ranking comparison data
            save_path: Path to save figure
        """
        fig, ax = plt.subplots(figsize=(10, 12))
        
        # Color scheme by constraint type
        colors = {
            'L': '#F18F01',  # Orange for Lagrangian
            'A': '#2E86AB',  # Blue for AMS
            'C': '#A23B72'   # Purple for CPC
        }
        
        # Plot connections
        for idx, row in df_comparison.iterrows():
            variant = row['Variant']
            rank_penalty = row['Rank_Penalty']
            rank_no_penalty = row['Rank_NoPenalty']
            
            # Determine color based on constraint type
            constraint_type = variant[0]  # First letter determines type
            color = colors.get(constraint_type, '#333333')
            
            # Draw line
            ax.plot([0, 1], [rank_penalty, rank_no_penalty], 
                   color=color, alpha=0.6, linewidth=1.5)
            
            # Add points
            ax.scatter(0, rank_penalty, color=color, s=50, zorder=5)
            ax.scatter(1, rank_no_penalty, color=color, s=50, zorder=5)
            
            # Add labels
            ax.text(-0.05, rank_penalty, variant, ha='right', va='center', fontsize=9)
            ax.text(1.05, rank_no_penalty, variant, ha='left', va='center', fontsize=9)
        
        # Customize plot
        ax.set_xlim(-0.2, 1.2)
        ax.set_ylim(max(df_comparison['Rank_Penalty'].max(), 
                       df_comparison['Rank_NoPenalty'].max()) + 0.5, 0.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['With Penalty\n(λ=0.1)', 'No Penalty'], fontsize=12)
        ax.set_ylabel('Rank (1 = Best)', fontsize=12)
        ax.set_title('Ranking Change: With vs Without CAI Penalty', fontsize=14, fontweight='bold')
        
        # Add legend
        legend_elements = [mpatches.Patch(color=colors['L'], label='Lagrangian', alpha=0.7),
                          mpatches.Patch(color=colors['A'], label='AMS', alpha=0.7),
                          mpatches.Patch(color=colors['C'], label='CPC', alpha=0.7)]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        
        if save_path:
            fig.savefig(save_path, dpi=self.figure_dpi, bbox_inches='tight')
            logger.info(f"Saved ranking change plot to {save_path}")
        
        return fig


if __name__ == "__main__":
    # Test visualization
    from result_parser import ExperimentResultParser
    
    # Load results
    parser = ExperimentResultParser("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250908_122909_unified_cai_experiments")
    results = parser.load_all_results()
    
    # Create visualizer
    viz = ExperimentVisualizer("./test_figures")
    
    # Test plots
    viz.plot_convergence_curves(results, save_path="./test_figures/convergence.png")
    viz.plot_performance_comparison(results, save_path="./test_figures/performance.png")
    viz.plot_protein_comparison(results, save_path="./test_figures/proteins.png")
    
    print("Test visualizations created in ./test_figures/")