#!/usr/bin/env python3
"""
Result Parser Module for Unified Experiment Analysis
Handles parsing and aggregation of experimental results from JSON files
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Any, Optional
import numpy as np
from collections import defaultdict
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class ExperimentResultParser:
    """Parser for unified experiment results"""
    
    def __init__(self, results_dir: str):
        """
        Initialize parser with results directory
        
        Args:
            results_dir: Path to directory containing JSON result files
        """
        self.results_dir = Path(results_dir)
        if not self.results_dir.exists():
            raise ValueError(f"Results directory does not exist: {results_dir}")
        
        self.raw_results = []
        self.aggregated_results = {}
        
    def load_all_results(self) -> List[Dict]:
        """Load all JSON result files from the directory"""
        json_files = list(self.results_dir.glob("*.json"))
        logger.info(f"Found {len(json_files)} result files")
        
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    # Add filename for reference
                    data['filename'] = json_file.name
                    self.raw_results.append(data)
            except Exception as e:
                logger.error(f"Error loading {json_file}: {e}")
                
        logger.info(f"Successfully loaded {len(self.raw_results)} results")
        return self.raw_results
    
    def aggregate_by_protein(self) -> Dict[str, List[Dict]]:
        """Aggregate results by protein name"""
        protein_results = defaultdict(list)
        
        for result in self.raw_results:
            protein_name = result.get('protein_name', 'Unknown')
            protein_results[protein_name].append(result)
            
        return dict(protein_results)
    
    def aggregate_by_constraint(self) -> Dict[str, List[Dict]]:
        """Aggregate results by constraint type"""
        constraint_results = defaultdict(list)
        
        for result in self.raw_results:
            constraint_type = result.get('constraint_type', 'Unknown')
            constraint_results[constraint_type].append(result)
            
        return dict(constraint_results)
    
    def aggregate_by_variant(self) -> Dict[str, List[Dict]]:
        """Aggregate results by variant (00, 01, 10, 11)"""
        variant_results = defaultdict(list)
        
        for result in self.raw_results:
            variant = result.get('variant', 'Unknown')
            variant_results[variant].append(result)
            
        return dict(variant_results)
    
    def aggregate_by_protein_constraint(self) -> Dict[str, Dict[str, List[Dict]]]:
        """Aggregate results by protein and constraint type"""
        pc_results = defaultdict(lambda: defaultdict(list))
        
        for result in self.raw_results:
            protein = result.get('protein_name', 'Unknown')
            constraint = result.get('constraint_type', 'Unknown')
            pc_results[protein][constraint].append(result)
            
        # Convert to regular dict
        return {p: dict(c) for p, c in pc_results.items()}
    
    def extract_key_metrics(self, result: Dict) -> Dict[str, Any]:
        """
        Extract key metrics from a single result
        
        Args:
            result: Single experiment result dictionary
            
        Returns:
            Dictionary of key metrics
        """
        metrics = {
            'protein_name': result.get('protein_name'),
            'constraint_type': result.get('constraint_type'),
            'variant': result.get('variant'),
            'initial_accessibility': result.get('initial_accessibility'),
            'final_accessibility': result.get('final_accessibility'),
            'best_accessibility': result.get('best_accessibility'),
            'improvement': result.get('improvement'),
            'iterations': result.get('iterations'),
            'optimization_time': result.get('optimization_time'),
            'amino_acids_match': result.get('amino_acids_match'),
            'cai_enabled': result.get('cai_enabled', False),
        }
        
        # Add CAI metrics if available
        if result.get('cai_enabled'):
            trajectory = result.get('trajectory', {})
            if 'ecai_values' in trajectory and trajectory['ecai_values']:
                metrics['initial_ecai'] = trajectory['ecai_values'][0]
                metrics['final_ecai'] = trajectory['ecai_values'][-1]
            if 'discrete_cai_values' in trajectory and trajectory['discrete_cai_values']:
                metrics['initial_discrete_cai'] = trajectory['discrete_cai_values'][0]
                metrics['final_discrete_cai'] = trajectory['discrete_cai_values'][-1]
                
        return metrics
    
    def get_summary_statistics(self) -> Dict[str, Any]:
        """Calculate summary statistics across all results"""
        if not self.raw_results:
            logger.warning("No results loaded")
            return {}
        
        all_metrics = [self.extract_key_metrics(r) for r in self.raw_results]
        
        # Calculate statistics
        improvements = [m['improvement'] for m in all_metrics if m['improvement'] is not None]
        final_accessibilities = [m['final_accessibility'] for m in all_metrics if m['final_accessibility'] is not None]
        best_accessibilities = [m['best_accessibility'] for m in all_metrics if m['best_accessibility'] is not None]
        
        stats = {
            'total_experiments': len(self.raw_results),
            'unique_proteins': len(set(m['protein_name'] for m in all_metrics)),
            'unique_constraints': len(set(m['constraint_type'] for m in all_metrics)),
            'unique_variants': len(set(m['variant'] for m in all_metrics)),
            'avg_improvement': np.mean(improvements) if improvements else 0,
            'std_improvement': np.std(improvements) if improvements else 0,
            'max_improvement': max(improvements) if improvements else 0,
            'min_improvement': min(improvements) if improvements else 0,
            'avg_final_accessibility': np.mean(final_accessibilities) if final_accessibilities else 0,
            'std_final_accessibility': np.std(final_accessibilities) if final_accessibilities else 0,
            'avg_best_accessibility': np.mean(best_accessibilities) if best_accessibilities else 0,
            'std_best_accessibility': np.std(best_accessibilities) if best_accessibilities else 0,
        }
        
        # CAI statistics if available
        cai_metrics = [m for m in all_metrics if m.get('cai_enabled')]
        if cai_metrics:
            final_ecais = [m.get('final_ecai', 0) for m in cai_metrics if m.get('final_ecai')]
            stats['cai_experiments'] = len(cai_metrics)
            stats['avg_final_ecai'] = np.mean(final_ecais) if final_ecais else 0
            stats['std_final_ecai'] = np.std(final_ecais) if final_ecais else 0
            
        return stats
    
    def find_best_configurations(self) -> Dict[str, Dict]:
        """Find best configuration for each protein"""
        protein_best = {}
        
        by_protein = self.aggregate_by_protein()
        
        for protein, results in by_protein.items():
            # Find result with lowest final accessibility
            best_result = min(results, key=lambda x: x.get('final_accessibility', float('inf')))
            
            protein_best[protein] = {
                'constraint_type': best_result.get('constraint_type'),
                'variant': best_result.get('variant'),
                'final_accessibility': best_result.get('final_accessibility'),
                'best_accessibility': best_result.get('best_accessibility'),
                'improvement': best_result.get('improvement'),
                'cai_enabled': best_result.get('cai_enabled', False),
            }
            
        return protein_best
    
    def get_performance_matrix(self) -> np.ndarray:
        """
        Create performance matrix for all experiments
        Rows: proteins, Columns: constraint_variant combinations
        """
        proteins = sorted(set(r['protein_name'] for r in self.raw_results))
        
        # Create constraint-variant combinations
        constraints = ['ams', 'cpc', 'lagrangian']
        variants = ['00', '01', '10', '11']
        combinations = [f"{c}_{v}" for c in constraints for v in variants]
        
        # Initialize matrix
        matrix = np.full((len(proteins), len(combinations)), np.nan)
        
        # Fill matrix
        for i, protein in enumerate(proteins):
            for j, combo in enumerate(combinations):
                constraint, variant = combo.split('_')
                
                # Find matching result
                for result in self.raw_results:
                    if (result['protein_name'] == protein and 
                        result['constraint_type'] == constraint and 
                        result['variant'] == variant):
                        matrix[i, j] = result.get('final_accessibility', np.nan)
                        break
                        
        return matrix, proteins, combinations
    
    def compare_cai_effect(self) -> Dict[str, Any]:
        """Compare results with and without CAI optimization"""
        cai_enabled = [r for r in self.raw_results if r.get('cai_enabled', False)]
        cai_disabled = [r for r in self.raw_results if not r.get('cai_enabled', False)]
        
        comparison = {
            'cai_enabled_count': len(cai_enabled),
            'cai_disabled_count': len(cai_disabled),
        }
        
        if cai_enabled:
            cai_improvements = [r['improvement'] for r in cai_enabled if r.get('improvement')]
            cai_finals = [r['final_accessibility'] for r in cai_enabled if r.get('final_accessibility')]
            comparison['cai_avg_improvement'] = np.mean(cai_improvements) if cai_improvements else 0
            comparison['cai_avg_final'] = np.mean(cai_finals) if cai_finals else 0
            
        if cai_disabled:
            no_cai_improvements = [r['improvement'] for r in cai_disabled if r.get('improvement')]
            no_cai_finals = [r['final_accessibility'] for r in cai_disabled if r.get('final_accessibility')]
            comparison['no_cai_avg_improvement'] = np.mean(no_cai_improvements) if no_cai_improvements else 0
            comparison['no_cai_avg_final'] = np.mean(no_cai_finals) if no_cai_finals else 0
            
        return comparison


if __name__ == "__main__":
    # Test the parser
    parser = ExperimentResultParser("/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/results/20250908_122909_unified_cai_experiments")
    
    # Load results
    results = parser.load_all_results()
    print(f"Loaded {len(results)} results")
    
    # Get summary statistics
    stats = parser.get_summary_statistics()
    print("\nSummary Statistics:")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")
    
    # Find best configurations
    best = parser.find_best_configurations()
    print("\nBest Configurations by Protein:")
    for protein, config in best.items():
        print(f"  {protein}: {config['constraint_type']}_{config['variant']} (accessibility: {config['final_accessibility']:.4f})")