"""
Sequence Uniqueness Analyzer

Core module for analyzing sequence diversity and uniqueness in CAI optimization experiments.
"""

import logging
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from collections import defaultdict, Counter

from .base_analyzer import BaseAnalyzer
from .diversity_metrics import DiversityMetrics

logger = logging.getLogger(__name__)


class SequenceUniquenessAnalyzer(BaseAnalyzer):
    """
    Analyzer for sequence uniqueness and diversity metrics.
    
    Focuses on understanding the exploration efficiency of the optimization process.
    """
    
    def get_analyzer_name(self) -> str:
        """Get the name of this analyzer."""
        return "uniqueness"
    
    def _perform_analysis(self) -> Dict[str, Any]:
        """
        Perform sequence uniqueness analysis.
        
        Returns:
            Dictionary containing analysis results
        """
        logger.info("Analyzing sequence uniqueness and diversity...")
        
        # Extract sequences from all experiments
        all_experiment_metrics = []
        grouped_metrics = defaultdict(list)
        
        for result in self.filtered_results:
            # Extract discrete sequences from trajectory
            if 'trajectory' not in result or 'discrete_sequences' not in result['trajectory']:
                logger.warning(f"No discrete sequences found in {result.get('_filename', 'unknown')}")
                continue
            
            sequences = result['trajectory']['discrete_sequences']
            
            # Calculate metrics for this experiment
            metrics = DiversityMetrics.calculate_all_metrics(sequences)
            
            # Add experiment metadata
            metrics.update({
                'protein': result.get('protein_name'),
                'constraint_type': result.get('constraint_type'),
                'variant': result.get('variant'),
                'seed': result.get('seed'),
                'best_accessibility': result.get('best_accessibility'),
                'final_accessibility': result.get('final_accessibility'),
                'sequence_length': len(sequences[0]) if sequences else 0,
                'num_iterations': len(sequences),
                'filename': result.get('_filename'),
            })
            
            # Analyze temporal diversity if requested
            if self.config.uniqueness_params.get('time_windows'):
                time_windows = self.config.uniqueness_params['time_windows']
                for window in time_windows:
                    if window <= len(sequences):
                        window_seqs = sequences[:window]
                        window_uniqueness = DiversityMetrics.uniqueness_rate(window_seqs)
                        metrics[f'uniqueness_at_{window}'] = window_uniqueness
            
            # Find most repeated sequences
            if self.config.uniqueness_params.get('top_k_sequences', 0) > 0:
                top_k = self.config.uniqueness_params['top_k_sequences']
                most_repeated = DiversityMetrics.most_repeated_sequences(sequences, top_k)
                metrics['most_repeated_sequences'] = most_repeated
            
            # Calculate diversity decay curve
            decay_curve = DiversityMetrics.temporal_diversity_decay(sequences, window_size=100)
            if decay_curve:
                metrics['diversity_decay_mean'] = np.mean(decay_curve)
                metrics['diversity_decay_final'] = decay_curve[-1] if decay_curve else None
            
            all_experiment_metrics.append(metrics)
            
            # Group by constraint type and variant
            group_key = f"{metrics['constraint_type']}_{metrics['variant']}"
            grouped_metrics[group_key].append(metrics)
        
        # Convert to DataFrame for easier analysis
        df = pd.DataFrame(all_experiment_metrics)
        
        # Calculate summary statistics
        summary = self._calculate_summary_statistics(df, grouped_metrics)
        
        # Analyze correlations
        correlations = self._analyze_correlations(df)
        
        # Identify patterns
        patterns = self._identify_patterns(df, grouped_metrics)
        
        return {
            'summary': summary,
            'correlations': correlations,
            'patterns': patterns,
            'dataframe': df,
            'detailed_metrics': all_experiment_metrics,
            'grouped_metrics': dict(grouped_metrics),
        }
    
    def _calculate_summary_statistics(self, df: pd.DataFrame, grouped_metrics: dict) -> dict:
        """
        Calculate summary statistics for uniqueness metrics.
        
        Args:
            df: DataFrame with all metrics
            grouped_metrics: Metrics grouped by configuration
            
        Returns:
            Dictionary of summary statistics
        """
        if df.empty:
            return {}
        
        summary = {
            'global_uniqueness_rate': df['uniqueness_rate'].mean(),
            'global_uniqueness_std': df['uniqueness_rate'].std(),
            'min_uniqueness_rate': df['uniqueness_rate'].min(),
            'max_uniqueness_rate': df['uniqueness_rate'].max(),
            'median_uniqueness_rate': df['uniqueness_rate'].median(),
            
            'total_experiments': len(df),
            'total_sequences': df['total_sequences'].sum(),
            'total_unique_sequences': df['unique_sequences'].sum(),
            
            'mean_shannon_entropy': df['shannon_entropy'].mean(),
            'mean_simpson_diversity': df['simpson_diversity'].mean(),
            'mean_repetition_factor': df['repetition_factor'].mean(),
        }
        
        # Statistics by constraint type
        constraint_stats = {}
        for constraint in df['constraint_type'].unique():
            constraint_df = df[df['constraint_type'] == constraint]
            constraint_stats[constraint] = {
                'mean_uniqueness': constraint_df['uniqueness_rate'].mean(),
                'std_uniqueness': constraint_df['uniqueness_rate'].std(),
                'min_uniqueness': constraint_df['uniqueness_rate'].min(),
                'max_uniqueness': constraint_df['uniqueness_rate'].max(),
                'count': len(constraint_df),
            }
        summary['by_constraint'] = constraint_stats
        
        # Statistics by variant
        variant_stats = {}
        for variant in df['variant'].unique():
            variant_df = df[df['variant'] == variant]
            variant_stats[variant] = {
                'mean_uniqueness': variant_df['uniqueness_rate'].mean(),
                'std_uniqueness': variant_df['uniqueness_rate'].std(),
                'min_uniqueness': variant_df['uniqueness_rate'].min(),
                'max_uniqueness': variant_df['uniqueness_rate'].max(),
                'count': len(variant_df),
            }
        summary['by_variant'] = variant_stats
        
        # Statistics by protein length bins
        if 'sequence_length' in df.columns:
            df['length_bin'] = pd.cut(df['sequence_length'], bins=5)
            length_stats = df.groupby('length_bin', observed=False)['uniqueness_rate'].agg(['mean', 'std', 'count'])
            # Convert to JSON-serializable format
            summary['by_sequence_length'] = {
                'bins': [str(bin) for bin in length_stats.index],
                'mean': length_stats['mean'].tolist(),
                'std': length_stats['std'].tolist(),
                'count': length_stats['count'].tolist(),
            }
        
        # Temporal analysis
        if 'uniqueness_at_100' in df.columns:
            summary['early_exploration'] = {
                'uniqueness_at_100': df['uniqueness_at_100'].mean() if 'uniqueness_at_100' in df else None,
                'uniqueness_at_500': df['uniqueness_at_500'].mean() if 'uniqueness_at_500' in df else None,
                'uniqueness_at_1000': df['uniqueness_at_1000'].mean() if 'uniqueness_at_1000' in df else None,
            }
        
        # Convergence analysis
        if 'convergence_iteration' in df.columns:
            converged = df['convergence_iteration'].notna()
            summary['convergence'] = {
                'experiments_converged': converged.sum(),
                'convergence_rate': converged.mean(),
                'mean_convergence_iteration': df.loc[converged, 'convergence_iteration'].mean() if converged.any() else None,
            }
        
        return summary
    
    def _analyze_correlations(self, df: pd.DataFrame) -> dict:
        """
        Analyze correlations between uniqueness and performance metrics.
        
        Args:
            df: DataFrame with metrics
            
        Returns:
            Dictionary of correlation results
        """
        if df.empty:
            return {}
        
        correlations = {}
        
        # Correlation with performance metrics
        if 'best_accessibility' in df.columns:
            correlations['uniqueness_vs_best_accessibility'] = df['uniqueness_rate'].corr(df['best_accessibility'])
            correlations['uniqueness_vs_final_accessibility'] = df['uniqueness_rate'].corr(df['final_accessibility'])
        
        # Correlation with sequence length
        if 'sequence_length' in df.columns:
            correlations['uniqueness_vs_sequence_length'] = df['uniqueness_rate'].corr(df['sequence_length'])
        
        # Correlation with diversity metrics
        correlations['uniqueness_vs_shannon_entropy'] = df['uniqueness_rate'].corr(df['shannon_entropy'])
        correlations['uniqueness_vs_simpson_diversity'] = df['uniqueness_rate'].corr(df['simpson_diversity'])
        
        # Correlation matrix for key metrics
        key_metrics = ['uniqueness_rate', 'shannon_entropy', 'simpson_diversity', 
                      'best_accessibility', 'repetition_factor']
        available_metrics = [m for m in key_metrics if m in df.columns]
        
        if len(available_metrics) > 1:
            correlation_matrix = df[available_metrics].corr()
            correlations['correlation_matrix'] = correlation_matrix.to_dict()
        
        return correlations
    
    def _identify_patterns(self, df: pd.DataFrame, grouped_metrics: dict) -> dict:
        """
        Identify patterns in sequence diversity.
        
        Args:
            df: DataFrame with metrics
            grouped_metrics: Metrics grouped by configuration
            
        Returns:
            Dictionary of identified patterns
        """
        patterns = {}
        
        # Find configurations with lowest/highest diversity
        if not df.empty:
            worst_config = df.nsmallest(1, 'uniqueness_rate').iloc[0]
            best_config = df.nlargest(1, 'uniqueness_rate').iloc[0]
            
            patterns['worst_diversity'] = {
                'protein': worst_config['protein'],
                'constraint': worst_config['constraint_type'],
                'variant': worst_config['variant'],
                'uniqueness_rate': worst_config['uniqueness_rate'],
                'filename': worst_config.get('filename', 'unknown'),
            }
            
            patterns['best_diversity'] = {
                'protein': best_config['protein'],
                'constraint': best_config['constraint_type'],
                'variant': best_config['variant'],
                'uniqueness_rate': best_config['uniqueness_rate'],
                'filename': best_config.get('filename', 'unknown'),
            }
        
        # Identify experiments with extreme repetition
        if 'max_repetition_streak' in df.columns:
            extreme_repetition = df.nlargest(5, 'max_repetition_streak')
            patterns['extreme_repetition'] = [
                {
                    'experiment': row['filename'],
                    'streak_length': row['max_repetition_streak'],
                    'uniqueness_rate': row['uniqueness_rate'],
                }
                for _, row in extreme_repetition.iterrows()
            ]
        
        # Analyze diversity decay patterns
        if 'diversity_decay_mean' in df.columns:
            patterns['diversity_decay'] = {
                'mean_decay_rate': df['diversity_decay_mean'].mean(),
                'fastest_decay': df.nsmallest(1, 'diversity_decay_mean')['filename'].iloc[0] if not df.empty else None,
                'slowest_decay': df.nlargest(1, 'diversity_decay_mean')['filename'].iloc[0] if not df.empty else None,
            }
        
        # Configuration rankings
        config_rankings = []
        for config_key in grouped_metrics:
            metrics_list = grouped_metrics[config_key]
            if metrics_list:
                mean_uniqueness = np.mean([m['uniqueness_rate'] for m in metrics_list])
                config_rankings.append((config_key, mean_uniqueness, len(metrics_list)))
        
        config_rankings.sort(key=lambda x: x[1], reverse=True)
        patterns['configuration_rankings'] = [
            {'configuration': cfg, 'mean_uniqueness': uniqueness, 'count': count}
            for cfg, uniqueness, count in config_rankings[:10]
        ]
        
        return patterns
    
    def _generate_visualizations(self) -> Dict[str, Path]:
        """
        Generate visualizations for uniqueness analysis.
        
        Returns:
            Dictionary mapping visualization names to file paths
        """
        from id3.experiments.analysis.visualizers.uniqueness_visualizer import UniquenessVisualizer
        
        visualizer = UniquenessVisualizer(self.output_dir)
        visualization_files = {}
        
        if 'dataframe' in self.analysis_results:
            df = self.analysis_results['dataframe']
            
            # Generate various plots
            files = visualizer.create_all_visualizations(
                df,
                self.analysis_results.get('summary', {}),
                self.analysis_results.get('patterns', {})
            )
            
            visualization_files.update(files)
        
        return visualization_files
    
    def _write_specific_summary(self, file_handle) -> None:
        """
        Write uniqueness-specific summary content.
        
        Args:
            file_handle: Open file handle for writing
        """
        summary = self.analysis_results.get('summary', {})
        patterns = self.analysis_results.get('patterns', {})
        
        # Global statistics
        file_handle.write("GLOBAL STATISTICS\n")
        file_handle.write("-" * 40 + "\n")
        file_handle.write(f"Average Uniqueness Rate: {summary.get('global_uniqueness_rate', 0):.1%}\n")
        file_handle.write(f"Uniqueness Std Dev: {summary.get('global_uniqueness_std', 0):.1%}\n")
        file_handle.write(f"Min Uniqueness: {summary.get('min_uniqueness_rate', 0):.1%}\n")
        file_handle.write(f"Max Uniqueness: {summary.get('max_uniqueness_rate', 0):.1%}\n")
        file_handle.write(f"Total Sequences Analyzed: {summary.get('total_sequences', 0):,}\n")
        file_handle.write(f"Total Unique Sequences: {summary.get('total_unique_sequences', 0):,}\n\n")
        
        # By constraint type
        if 'by_constraint' in summary:
            file_handle.write("BY CONSTRAINT TYPE\n")
            file_handle.write("-" * 40 + "\n")
            for constraint, stats in summary['by_constraint'].items():
                file_handle.write(f"{constraint:12s}: {stats['mean_uniqueness']:.1%} ± {stats['std_uniqueness']:.1%}\n")
            file_handle.write("\n")
        
        # By variant
        if 'by_variant' in summary:
            file_handle.write("BY VARIANT\n")
            file_handle.write("-" * 40 + "\n")
            for variant, stats in summary['by_variant'].items():
                file_handle.write(f"Variant {variant}: {stats['mean_uniqueness']:.1%} ± {stats['std_uniqueness']:.1%}\n")
            file_handle.write("\n")
        
        # Key findings
        file_handle.write("KEY FINDINGS\n")
        file_handle.write("-" * 40 + "\n")
        
        if 'worst_diversity' in patterns:
            worst = patterns['worst_diversity']
            file_handle.write(f"Worst Diversity: {worst['uniqueness_rate']:.1%} "
                            f"({worst['protein']}_{worst['constraint']}_{worst['variant']})\n")
        
        if 'best_diversity' in patterns:
            best = patterns['best_diversity']
            file_handle.write(f"Best Diversity: {best['uniqueness_rate']:.1%} "
                            f"({best['protein']}_{best['constraint']}_{best['variant']})\n")
        
        # Recommendations
        file_handle.write("\nRECOMMENDATIONS\n")
        file_handle.write("-" * 40 + "\n")
        
        avg_uniqueness = summary.get('global_uniqueness_rate', 0)
        if avg_uniqueness < 0.05:
            file_handle.write("⚠️  CRITICAL: Extremely low sequence diversity detected!\n")
            file_handle.write("   Average uniqueness rate is only {:.1%}\n".format(avg_uniqueness))
            file_handle.write("   Recommendations:\n")
            file_handle.write("   1. Implement diversity enhancement mechanisms (e.g., SADO)\n")
            file_handle.write("   2. Adjust temperature parameters to increase exploration\n")
            file_handle.write("   3. Add diversity rewards to the optimization objective\n")
            file_handle.write("   4. Consider using ensemble methods with different initializations\n")
        elif avg_uniqueness < 0.10:
            file_handle.write("⚠️  WARNING: Low sequence diversity\n")
            file_handle.write("   Consider implementing diversity enhancement strategies\n")
        else:
            file_handle.write("✓  Sequence diversity is reasonable\n")