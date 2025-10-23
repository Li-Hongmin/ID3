"""
CAI Experiment Configuration - Incremental Extension

This module provides CAI-enhanced experiment configurations that extend
the existing 12x12 experiment framework without modifying the original code.

Design Principles:
- Incremental: Extends, not replaces existing configurations
- Modular: Self-contained CAI configuration module  
- Unix Philosophy: Do one thing well - add CAI to experiments
- Non-intrusive: Original experiments continue to work unchanged
"""

from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field


@dataclass
class CAIExperimentConfig:
    """
    CAI-specific experiment configuration.
    
    This configuration can be composed with existing experiment configs
    to add CAI functionality without modification.
    """
    
    # CAI optimization parameters
    cai_enabled: bool = False
    cai_target: float = 0.8
    cai_method: str = 'retreat'  # 'retreat', 'advance', 'dp_optimal'
    cai_weight: float = 0.1
    species: str = 'ecoli_bl21de3'
    
    # Experiment variants with CAI
    cai_variants: List[str] = field(default_factory=lambda: [
        'lagrangian-cai-00', 'lagrangian-cai-01', 'lagrangian-cai-10', 'lagrangian-cai-11',
        'cpc-cai-00', 'cpc-cai-01', 'cpc-cai-10', 'cpc-cai-11',
        'ams-cai-00', 'ams-cai-01', 'ams-cai-10', 'ams-cai-11'
    ])
    
    # Performance tuning
    enable_caching: bool = True
    batch_optimization: bool = True
    fallback_to_rams: bool = True
    
    # Results tracking
    track_cai_metrics: bool = True
    save_cai_sequences: bool = True
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            'cai_enabled': self.cai_enabled,
            'cai_target': self.cai_target,
            'cai_method': self.cai_method,
            'cai_weight': self.cai_weight,
            'species': self.species,
            'cai_variants': self.cai_variants,
            'enable_caching': self.enable_caching,
            'batch_optimization': self.batch_optimization,
            'fallback_to_rams': self.fallback_to_rams,
            'track_cai_metrics': self.track_cai_metrics,
            'save_cai_sequences': self.save_cai_sequences
        }
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'CAIExperimentConfig':
        """Create configuration from dictionary."""
        return cls(**config_dict)


class CAIVariantBuilder:
    """
    Builder for CAI-enhanced experiment variants.
    
    This builder creates CAI variants that can be added to the existing
    12x12 experiment matrix without modifying the original variants.
    """
    
    @staticmethod
    def build_cai_variants(base_constraints: List[str] = None,
                          alpha_beta_combinations: List[tuple] = None,
                          cai_config: CAIExperimentConfig = None) -> Dict[str, Dict]:
        """
        Build CAI-enhanced variants for experiments.
        
        Args:
            base_constraints: Base constraint types (lagrangian, cpc, ams)
            alpha_beta_combinations: (alpha, beta) tuples for variants
            cai_config: CAI configuration
            
        Returns:
            Dictionary of CAI variant configurations
        """
        if base_constraints is None:
            base_constraints = ['lagrangian', 'cpc', 'ams']
        
        if alpha_beta_combinations is None:
            alpha_beta_combinations = [
                (0, 0), (0, 1), (1, 0), (1, 1)
            ]
        
        if cai_config is None:
            cai_config = CAIExperimentConfig(cai_enabled=True)
        
        variants = {}
        
        for constraint in base_constraints:
            for alpha, beta in alpha_beta_combinations:
                variant_name = f"{constraint}-cai-{alpha}{beta}"
                
                variants[variant_name] = {
                    'constraint_type': constraint,
                    'alpha': alpha,
                    'beta': beta,
                    'cai_enabled': True,
                    'cai_config': cai_config.to_dict(),
                    'description': f"CAI-enhanced {constraint} with α={alpha}, β={beta}"
                }
        
        return variants


class CAIExperimentExtension:
    """
    Extension module for adding CAI to existing experiments.
    
    This class provides methods to enhance existing experiment configurations
    with CAI functionality without modifying the original experiment code.
    """
    
    def __init__(self, cai_config: Optional[CAIExperimentConfig] = None):
        """Initialize CAI extension."""
        self.cai_config = cai_config or CAIExperimentConfig()
        self.variant_builder = CAIVariantBuilder()
    
    def extend_experiment_config(self, 
                                base_config: Dict[str, Any],
                                add_cai_variants: bool = True) -> Dict[str, Any]:
        """
        Extend existing experiment configuration with CAI.
        
        Args:
            base_config: Original experiment configuration
            add_cai_variants: Whether to add CAI variants to the matrix
            
        Returns:
            Extended configuration with CAI support
        """
        # Create extended config as a copy (non-destructive)
        extended_config = base_config.copy()
        
        # Add CAI configuration
        extended_config['cai_extension'] = self.cai_config.to_dict()
        
        # Optionally add CAI variants
        if add_cai_variants:
            cai_variants = self.variant_builder.build_cai_variants(
                cai_config=self.cai_config
            )
            
            # Add to variants list if it exists
            if 'variants' in extended_config:
                # Preserve original variants
                original_variants = extended_config.get('variants', {})
                if isinstance(original_variants, list):
                    # Convert list to dict format
                    original_dict = {v: {} for v in original_variants}
                    extended_config['variants'] = {**original_dict, **cai_variants}
                else:
                    extended_config['variants'] = {**original_variants, **cai_variants}
            else:
                extended_config['cai_variants'] = cai_variants
        
        return extended_config
    
    def get_cai_constraint_factory(self, constraint_type: str):
        """
        Get CAI-enhanced constraint factory for a given type.
        
        Args:
            constraint_type: Type of constraint (lagrangian, cpc, ams)
            
        Returns:
            CAI-enhanced constraint factory function
        """
        from id3.cai.integration import (
            create_cai_enhanced_lagrangian,
            create_cai_enhanced_codon_profile
        )
        
        factories = {
            'lagrangian': lambda **kwargs: create_cai_enhanced_lagrangian(
                cai_target=self.cai_config.cai_target,
                cai_weight=self.cai_config.cai_weight,
                cai_method=self.cai_config.cai_method,
                **kwargs
            ),
            'cpc': lambda **kwargs: create_cai_enhanced_codon_profile(
                cai_target=self.cai_config.cai_target,
                cai_weight=self.cai_config.cai_weight,
                cai_method=self.cai_config.cai_method,
                **kwargs
            ),
            # AMS not yet implemented, use fallback
            'ams': lambda **kwargs: create_cai_enhanced_lagrangian(
                cai_target=self.cai_config.cai_target,
                cai_weight=self.cai_config.cai_weight,
                cai_method=self.cai_config.cai_method,
                **kwargs
            )
        }
        
        return factories.get(constraint_type, factories['lagrangian'])


# Convenience functions for quick setup
def get_default_cai_config() -> CAIExperimentConfig:
    """Get default CAI experiment configuration."""
    return CAIExperimentConfig(
        cai_enabled=True,
        cai_target=0.8,
        cai_method='retreat',
        cai_weight=0.1
    )


def extend_12x12_with_cai(proteins: List[str] = None,
                         cai_target: float = 0.8,
                         cai_method: str = 'retreat') -> Dict[str, Any]:
    """
    Create extended 12x12 experiment configuration with CAI.
    
    This creates a complete experiment configuration that includes
    both the original 12 variants and 12 new CAI-enhanced variants.
    
    Args:
        proteins: List of protein IDs (default: standard 12 proteins)
        cai_target: Target CAI value
        cai_method: CAI optimization method
        
    Returns:
        Complete experiment configuration with 24 variants
    """
    if proteins is None:
        # Use standard 12 proteins from original experiments
        proteins = [
            'P00004', 'P01308', 'P01825', 'P04637', 'P07527',
            'P0CG48', 'P24385', 'P42212', 'Q9Y6K9', 'Q15165',
            'Q9H9Y6', 'P01137'
        ]
    
    # Create CAI configuration
    cai_config = CAIExperimentConfig(
        cai_enabled=True,
        cai_target=cai_target,
        cai_method=cai_method
    )
    
    # Build complete configuration
    config = {
        'experiment_type': '12x24_with_cai',  # 12 proteins x 24 variants
        'proteins': proteins,
        'original_variants': [
            'lagrangian-00', 'lagrangian-01', 'lagrangian-10', 'lagrangian-11',
            'cpc-00', 'cpc-01', 'cpc-10', 'cpc-11',
            'ams-00', 'ams-01', 'ams-10', 'ams-11'
        ],
        'cai_extension': cai_config.to_dict(),
        'cai_variants': CAIVariantBuilder.build_cai_variants(cai_config=cai_config),
        'total_experiments': 12 * 24,  # 288 total experiments
        'description': 'Extended 12x12 experiment with CAI optimization'
    }
    
    return config


# Export key components
__all__ = [
    'CAIExperimentConfig',
    'CAIVariantBuilder',
    'CAIExperimentExtension',
    'get_default_cai_config',
    'extend_12x12_with_cai'
]