"""
Analysis Configuration Module

Provides configuration classes and presets for experiment analysis.
Follows the pattern of unified_experiment_config.py.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)


@dataclass
class AnalysisConfig:
    """
    Configuration for experiment analysis.
    
    Attributes:
        input_dir: Directory containing experiment results
        output_subdir: Subdirectory name for analysis results (under input_dir/analysis/)
        analysis_types: List of analysis types to perform
        visualize: Whether to generate visualizations
        output_format: Output format for data files ('json', 'csv', 'both')
        verbose: Show detailed progress
        parallel: Number of parallel processes
        filter_proteins: Only analyze specific proteins
        filter_constraints: Only analyze specific constraint types
        filter_variants: Only analyze specific variants
    """
    
    input_dir: Optional[Path] = None
    output_subdir: Optional[str] = None
    analysis_types: List[str] = field(default_factory=lambda: ['uniqueness'])
    visualize: bool = True
    output_format: str = 'both'
    verbose: bool = False
    parallel: int = 1
    
    # Filters
    filter_proteins: Optional[List[str]] = None
    filter_constraints: Optional[List[str]] = None
    filter_variants: Optional[List[str]] = None
    
    # Analysis-specific parameters
    uniqueness_params: Dict[str, Any] = field(default_factory=lambda: {
        'compute_entropy': True,
        'compute_hamming_distance': True,
        'time_windows': [100, 500, 1000],  # Analyze uniqueness at these iteration points
        'top_k_sequences': 10,  # Show top K most repeated sequences
    })
    
    performance_params: Dict[str, Any] = field(default_factory=lambda: {
        'metrics': ['best_accessibility', 'final_accessibility', 'improvement'],
        'group_by': ['constraint_type', 'variant'],
        'compute_statistics': True,
    })
    
    convergence_params: Dict[str, Any] = field(default_factory=lambda: {
        'window_size': 50,
        'convergence_threshold': 0.001,
        'plot_trajectories': True,
    })
    
    def validate(self) -> None:
        """
        Validate configuration.
        
        Raises:
            ValueError: If configuration is invalid
        """
        if not self.input_dir:
            raise ValueError("input_dir is required")
        
        if not self.input_dir.exists():
            raise FileNotFoundError(f"Input directory does not exist: {self.input_dir}")
        
        if not self.analysis_types:
            raise ValueError("At least one analysis type must be specified")
        
        # Check for valid analysis types
        valid_types = {'uniqueness', 'performance', 'convergence', 'comprehensive'}
        invalid_types = set(self.analysis_types) - valid_types
        if invalid_types:
            raise ValueError(f"Invalid analysis types: {invalid_types}")
        
        # Expand 'comprehensive' to include all types
        if 'comprehensive' in self.analysis_types:
            self.analysis_types = ['uniqueness', 'performance', 'convergence']
    
    def get_output_dir(self) -> Path:
        """
        Get output directory path.
        Creates the directory if it doesn't exist.
        
        Returns:
            Path to output directory
        """
        if self.output_subdir:
            output_dir = self.input_dir / 'analysis' / self.output_subdir
        else:
            # Auto-generate subdirectory name based on analysis types
            subdir_name = '_'.join(sorted(self.analysis_types))
            output_dir = self.input_dir / 'analysis' / subdir_name
        
        output_dir.mkdir(parents=True, exist_ok=True)
        return output_dir
    
    def print_summary(self) -> None:
        """Print configuration summary."""
        print("\n" + "=" * 60)
        print("Analysis Configuration")
        print("=" * 60)
        print(f"Input directory: {self.input_dir}")
        print(f"Output directory: {self.get_output_dir()}")
        print(f"Analysis types: {', '.join(self.analysis_types)}")
        print(f"Visualization: {self.visualize}")
        print(f"Output format: {self.output_format}")
        
        if self.filter_proteins:
            print(f"Filter proteins: {', '.join(self.filter_proteins)}")
        if self.filter_constraints:
            print(f"Filter constraints: {', '.join(self.filter_constraints)}")
        if self.filter_variants:
            print(f"Filter variants: {', '.join(self.filter_variants)}")
        
        print("-" * 60)
    
    def to_dict(self) -> dict:
        """Convert configuration to dictionary."""
        return {
            'input_dir': str(self.input_dir) if self.input_dir else None,
            'output_dir': str(self.get_output_dir()),
            'analysis_types': self.analysis_types,
            'visualize': self.visualize,
            'output_format': self.output_format,
            'verbose': self.verbose,
            'parallel': self.parallel,
            'filters': {
                'proteins': self.filter_proteins,
                'constraints': self.filter_constraints,
                'variants': self.filter_variants,
            },
            'parameters': {
                'uniqueness': self.uniqueness_params,
                'performance': self.performance_params,
                'convergence': self.convergence_params,
            }
        }


class AnalysisPresets:
    """
    Preset configurations for common analysis scenarios.
    """
    
    @staticmethod
    def quick_uniqueness() -> AnalysisConfig:
        """
        Quick uniqueness analysis preset.
        Focuses on sequence diversity metrics.
        """
        config = AnalysisConfig()
        config.analysis_types = ['uniqueness']
        config.visualize = True
        config.uniqueness_params = {
            'compute_entropy': True,
            'compute_hamming_distance': False,  # Skip for speed
            'time_windows': [100, 500, 1000],
            'top_k_sequences': 5,
        }
        return config
    
    @staticmethod
    def comprehensive() -> AnalysisConfig:
        """
        Comprehensive analysis preset.
        Includes all analysis types with full parameters.
        """
        config = AnalysisConfig()
        config.analysis_types = ['uniqueness', 'performance', 'convergence']
        config.visualize = True
        config.output_format = 'both'
        config.uniqueness_params = {
            'compute_entropy': True,
            'compute_hamming_distance': True,
            'time_windows': [100, 200, 500, 800, 1000],
            'top_k_sequences': 20,
        }
        config.performance_params = {
            'metrics': ['best_accessibility', 'final_accessibility', 'improvement', 'optimization_time'],
            'group_by': ['constraint_type', 'variant', 'protein_name'],
            'compute_statistics': True,
            'compute_correlations': True,
        }
        config.convergence_params = {
            'window_size': 50,
            'convergence_threshold': 0.001,
            'plot_trajectories': True,
            'analyze_plateaus': True,
        }
        return config
    
    @staticmethod
    def performance_only() -> AnalysisConfig:
        """
        Performance-focused analysis preset.
        Analyzes optimization performance metrics.
        """
        config = AnalysisConfig()
        config.analysis_types = ['performance']
        config.visualize = True
        config.performance_params = {
            'metrics': ['best_accessibility', 'final_accessibility', 'improvement'],
            'group_by': ['constraint_type', 'variant'],
            'compute_statistics': True,
            'generate_tables': True,
        }
        return config
    
    @staticmethod
    def from_experiment_config(experiment_config_path: Path) -> AnalysisConfig:
        """
        Create analysis config based on experiment configuration.
        
        Args:
            experiment_config_path: Path to experiment config.json
            
        Returns:
            AnalysisConfig instance
        """
        import json
        
        with open(experiment_config_path) as f:
            exp_config = json.load(f)
        
        config = AnalysisConfig()
        config.input_dir = experiment_config_path.parent
        
        # Set analysis based on experiment type
        if exp_config.get('enable_cai'):
            # CAI experiments need uniqueness analysis
            config.analysis_types = ['uniqueness', 'performance']
        else:
            config.analysis_types = ['performance', 'convergence']
        
        return config