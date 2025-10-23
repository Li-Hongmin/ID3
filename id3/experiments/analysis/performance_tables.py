#!/usr/bin/env python3
"""
Performance Tables Module for Unified Experiment Analysis
Generates various performance tables in multiple formats
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PerformanceTableGenerator:
    """Generator for performance tables"""
    
    def __init__(self, output_dir: str = None):
        """
        Initialize table generator
        
        Args:
            output_dir: Directory to save tables
        """
        if output_dir:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.output_dir = Path.cwd()
    
    def create_comprehensive_performance_table(self, results: List[Dict]) -> pd.DataFrame:
        """
        Create comprehensive performance table with all key metrics
        
        Args:
            results: List of experiment results
            
        Returns:
            DataFrame with performance metrics
        """
        table_data = []
        
        for result in results:
            row = {
                'Protein': result.get('protein_name'),
                'Constraint': result.get('constraint_type'),
                'Variant': result.get('variant'),
                'Initial_Acc': result.get('initial_accessibility'),
                'Final_Acc': result.get('final_accessibility'),
                'Best_Acc': result.get('best_accessibility'),
                'Improvement': result.get('improvement'),
                'Iterations': result.get('iterations'),
                'Time(s)': result.get('optimization_time'),
                'AA_Match': result.get('amino_acids_match'),
                'CAI_Enabled': result.get('cai_enabled', False),
            }
            
            # Add CAI metrics if available
            if result.get('cai_enabled'):
                trajectory = result.get('trajectory', {})
                if 'ecai_values' in trajectory and trajectory['ecai_values']:
                    row['Initial_ECAI'] = trajectory['ecai_values'][0]
                    row['Final_ECAI'] = trajectory['ecai_values'][-1]
                    row['ECAI_Improvement'] = trajectory['ecai_values'][-1] - trajectory['ecai_values'][0]
            
            table_data.append(row)
        
        df = pd.DataFrame(table_data)
        
        # Sort by protein, constraint, variant
        df = df.sort_values(['Protein', 'Constraint', 'Variant'])
        
        return df
    
    def create_best_configuration_table(self, results: List[Dict]) -> pd.DataFrame:
        """
        Create table showing best configuration for each protein
        
        Args:
            results: List of experiment results
            
        Returns:
            DataFrame with best configurations
        """
        # Group by protein and find best
        protein_best = {}
        
        for result in results:
            protein = result['protein_name']
            final_acc = result.get('final_accessibility', float('inf'))
            
            if protein not in protein_best or final_acc < protein_best[protein]['Final_Acc']:
                protein_best[protein] = {
                    'Protein': protein,
                    'Best_Constraint': result['constraint_type'],
                    'Best_Variant': result['variant'],
                    'Final_Acc': final_acc,
                    'Best_Acc': result.get('best_accessibility'),
                    'Improvement': result.get('improvement'),
                    'CAI_Enabled': result.get('cai_enabled', False),
                    'Configuration': f"{result['constraint_type']}_{result['variant']}"
                }
        
        df = pd.DataFrame(list(protein_best.values()))
        df = df.sort_values('Final_Acc')
        
        return df
    
    def create_cai_comparison_table(self, results: List[Dict]) -> pd.DataFrame:
        """
        Create table comparing performance with and without CAI
        
        Args:
            results: List of experiment results
            
        Returns:
            DataFrame with CAI comparison
        """
        # Separate CAI and non-CAI results
        cai_enabled = [r for r in results if r.get('cai_enabled', False)]
        cai_disabled = [r for r in results if not r.get('cai_enabled', False)]
        
        comparison_data = []
        
        # Group by protein and constraint
        proteins = set(r['protein_name'] for r in results)
        constraints = set(r['constraint_type'] for r in results)
        
        for protein in sorted(proteins):
            for constraint in sorted(constraints):
                # Find matching results
                cai_results = [r for r in cai_enabled 
                             if r['protein_name'] == protein and r['constraint_type'] == constraint]
                no_cai_results = [r for r in cai_disabled 
                                if r['protein_name'] == protein and r['constraint_type'] == constraint]
                
                if cai_results or no_cai_results:
                    row = {
                        'Protein': protein,
                        'Constraint': constraint,
                    }
                    
                    if cai_results:
                        cai_finals = [r['final_accessibility'] for r in cai_results]
                        cai_improvements = [r['improvement'] for r in cai_results]
                        row['CAI_Avg_Final'] = np.mean(cai_finals)
                        row['CAI_Std_Final'] = np.std(cai_finals)
                        row['CAI_Avg_Improvement'] = np.mean(cai_improvements)
                        
                        # Add ECAI metrics
                        ecai_finals = []
                        for r in cai_results:
                            trajectory = r.get('trajectory', {})
                            if 'ecai_values' in trajectory and trajectory['ecai_values']:
                                ecai_finals.append(trajectory['ecai_values'][-1])
                        if ecai_finals:
                            row['CAI_Avg_ECAI'] = np.mean(ecai_finals)
                    
                    if no_cai_results:
                        no_cai_finals = [r['final_accessibility'] for r in no_cai_results]
                        no_cai_improvements = [r['improvement'] for r in no_cai_results]
                        row['NoCAI_Avg_Final'] = np.mean(no_cai_finals)
                        row['NoCAI_Std_Final'] = np.std(no_cai_finals)
                        row['NoCAI_Avg_Improvement'] = np.mean(no_cai_improvements)
                    
                    # Calculate differences if both exist
                    if 'CAI_Avg_Final' in row and 'NoCAI_Avg_Final' in row:
                        row['Diff_Final'] = row['CAI_Avg_Final'] - row['NoCAI_Avg_Final']
                        row['Diff_Improvement'] = row.get('CAI_Avg_Improvement', 0) - row.get('NoCAI_Avg_Improvement', 0)
                    
                    comparison_data.append(row)
        
        df = pd.DataFrame(comparison_data)
        return df
    
    def create_statistical_summary_table(self, results: List[Dict]) -> pd.DataFrame:
        """
        Create statistical summary table
        
        Args:
            results: List of experiment results
            
        Returns:
            DataFrame with statistical summary
        """
        stats_data = []
        
        # Overall statistics
        all_finals = [r['final_accessibility'] for r in results if r.get('final_accessibility')]
        all_improvements = [r['improvement'] for r in results if r.get('improvement')]
        
        stats_data.append({
            'Category': 'Overall',
            'Metric': 'All Experiments',
            'Count': len(results),
            'Mean': np.mean(all_finals) if all_finals else np.nan,
            'Std': np.std(all_finals) if all_finals else np.nan,
            'Min': np.min(all_finals) if all_finals else np.nan,
            'Max': np.max(all_finals) if all_finals else np.nan,
            'Median': np.median(all_finals) if all_finals else np.nan,
        })
        
        # By constraint type
        for constraint in ['ams', 'cpc', 'lagrangian']:
            constraint_results = [r for r in results if r['constraint_type'] == constraint]
            finals = [r['final_accessibility'] for r in constraint_results if r.get('final_accessibility')]
            
            if finals:
                stats_data.append({
                    'Category': 'Constraint',
                    'Metric': constraint.upper(),
                    'Count': len(constraint_results),
                    'Mean': np.mean(finals),
                    'Std': np.std(finals),
                    'Min': np.min(finals),
                    'Max': np.max(finals),
                    'Median': np.median(finals),
                })
        
        # By variant
        for variant in ['00', '01', '10', '11']:
            variant_results = [r for r in results if r['variant'] == variant]
            finals = [r['final_accessibility'] for r in variant_results if r.get('final_accessibility')]
            
            if finals:
                stats_data.append({
                    'Category': 'Variant',
                    'Metric': variant,
                    'Count': len(variant_results),
                    'Mean': np.mean(finals),
                    'Std': np.std(finals),
                    'Min': np.min(finals),
                    'Max': np.max(finals),
                    'Median': np.median(finals),
                })
        
        # By protein (top 5 and bottom 5)
        protein_means = {}
        for protein in set(r['protein_name'] for r in results):
            protein_results = [r for r in results if r['protein_name'] == protein]
            finals = [r['final_accessibility'] for r in protein_results if r.get('final_accessibility')]
            if finals:
                protein_means[protein] = np.mean(finals)
        
        sorted_proteins = sorted(protein_means.items(), key=lambda x: x[1])
        
        # Add top 5 performers
        for i, (protein, mean_val) in enumerate(sorted_proteins[:5]):
            protein_results = [r for r in results if r['protein_name'] == protein]
            finals = [r['final_accessibility'] for r in protein_results if r.get('final_accessibility')]
            
            stats_data.append({
                'Category': f'Top Protein #{i+1}',
                'Metric': protein,
                'Count': len(protein_results),
                'Mean': mean_val,
                'Std': np.std(finals),
                'Min': np.min(finals),
                'Max': np.max(finals),
                'Median': np.median(finals),
            })
        
        df = pd.DataFrame(stats_data)
        return df
    
    def create_variant_performance_matrix(self, results: List[Dict]) -> pd.DataFrame:
        """
        Create matrix showing performance across variants
        
        Args:
            results: List of experiment results
            
        Returns:
            DataFrame matrix
        """
        # Create pivot table
        data_for_pivot = []
        
        for result in results:
            data_for_pivot.append({
                'Protein': result['protein_name'],
                'Config': f"{result['constraint_type']}_{result['variant']}",
                'Final_Acc': result.get('final_accessibility', np.nan)
            })
        
        df = pd.DataFrame(data_for_pivot)
        
        # Create pivot table
        pivot = df.pivot_table(
            index='Protein',
            columns='Config',
            values='Final_Acc',
            aggfunc='mean'
        )
        
        # Sort columns in a logical order
        col_order = []
        for constraint in ['ams', 'cpc', 'lagrangian']:
            for variant in ['00', '01', '10', '11']:
                col_name = f"{constraint}_{variant}"
                if col_name in pivot.columns:
                    col_order.append(col_name)
        
        pivot = pivot[col_order]
        
        # Add best configuration column
        numeric_cols = [col for col in pivot.columns if col != 'Best_Config']
        pivot['Best_Config'] = pivot[numeric_cols].idxmin(axis=1)
        pivot['Best_Score'] = pivot[numeric_cols].min(axis=1)
        
        # Sort by best score
        pivot = pivot.sort_values('Best_Score')
        
        return pivot
    
    def save_tables(self, results: List[Dict], format: str = 'all'):
        """
        Save all tables in specified format(s)
        
        Args:
            results: List of experiment results
            format: 'csv', 'latex', 'markdown', or 'all'
        """
        # Generate all tables
        tables = {
            'comprehensive_performance': self.create_comprehensive_performance_table(results),
            'best_configurations': self.create_best_configuration_table(results),
            'cai_comparison': self.create_cai_comparison_table(results),
            'statistical_summary': self.create_statistical_summary_table(results),
            'variant_matrix': self.create_variant_performance_matrix(results),
        }
        
        # Save in requested formats
        formats_to_save = ['csv', 'latex', 'markdown'] if format == 'all' else [format]
        
        for table_name, df in tables.items():
            for fmt in formats_to_save:
                if fmt == 'csv':
                    filepath = self.output_dir / f"{table_name}.csv"
                    df.to_csv(filepath, index=True)
                    logger.info(f"Saved {table_name}.csv")
                    
                elif fmt == 'latex':
                    filepath = self.output_dir / f"{table_name}.tex"
                    latex_str = df.to_latex(index=True, float_format="%.3f")
                    with open(filepath, 'w') as f:
                        f.write(latex_str)
                    logger.info(f"Saved {table_name}.tex")
                    
                elif fmt == 'markdown':
                    filepath = self.output_dir / f"{table_name}.md"
                    md_str = df.to_markdown(index=True, floatfmt=".3f")
                    with open(filepath, 'w') as f:
                        f.write(md_str)
                    logger.info(f"Saved {table_name}.md")
    
    def print_summary(self, results: List[Dict]):
        """
        Print summary of results to console
        
        Args:
            results: List of experiment results
        """
        print("\n" + "="*60)
        print("EXPERIMENT RESULTS SUMMARY")
        print("="*60)
        
        # Overall stats
        total = len(results)
        cai_enabled = sum(1 for r in results if r.get('cai_enabled', False))
        
        print(f"\nTotal Experiments: {total}")
        print(f"CAI-Enabled: {cai_enabled}")
        print(f"CAI-Disabled: {total - cai_enabled}")
        
        # Best overall result
        best_result = min(results, key=lambda x: x.get('final_accessibility', float('inf')))
        print(f"\nBest Overall Result:")
        print(f"  Protein: {best_result['protein_name']}")
        print(f"  Configuration: {best_result['constraint_type']}_{best_result['variant']}")
        print(f"  Final Accessibility: {best_result['final_accessibility']:.4f}")
        print(f"  Improvement: {best_result['improvement']:.4f}")
        
        # Average performance by constraint
        print("\nAverage Performance by Constraint Type:")
        for constraint in ['ams', 'cpc', 'lagrangian']:
            constraint_results = [r for r in results if r['constraint_type'] == constraint]
            if constraint_results:
                avg_final = np.mean([r['final_accessibility'] for r in constraint_results 
                                    if r.get('final_accessibility')])
                print(f"  {constraint.upper()}: {avg_final:.4f}")
        
        print("\n" + "="*60)


if __name__ == "__main__":
    # Test table generation
    from result_parser import ExperimentResultParser
    
    # Load results
    parser = ExperimentResultParser("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250908_122909_unified_cai_experiments")
    results = parser.load_all_results()
    
    # Create table generator
    gen = PerformanceTableGenerator("./test_tables")
    
    # Generate and save tables
    gen.save_tables(results, format='all')
    gen.print_summary(results)
    
    print("\nTest tables created in ./test_tables/")