"""
Constraint Configuration Module

This module provides configuration management for all three constraint types,
including adaptive lambda_cai settings and backward compatibility support.

Provides centralized configuration for:
- Fixed vs adaptive lambda_cai modes
- Hyperparameter defaults for each constraint type
- Experimental settings and research parameters
"""

from typing import Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class AdaptiveLambdaCAIConfig:
    """Configuration for adaptive lambda_cai adjustment mechanism."""
    
    # Core adaptive settings
    adaptive_lambda_cai: bool = False  # Default: disabled for backward compatibility
    lambda_cai_lr: float = 0.1  # Learning rate for adaptation
    lambda_cai_max: float = 2.0  # Maximum allowed lambda_cai value
    lambda_cai_min: float = 0.01  # Minimum allowed lambda_cai value
    cai_tolerance: float = 0.05  # Tolerance for CAI target satisfaction
    
    # Advanced control parameters
    smoothing_factor: float = 0.9  # Exponential smoothing for CAI measurements
    convergence_window: int = 20  # Window size for convergence analysis
    update_frequency: int = 1  # How often to update (every N iterations)
    
    # Safety and stability
    max_lambda_change_per_step: float = 0.5  # Maximum change per update
    stability_threshold: float = 0.01  # Stability detection threshold
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary for easy passing to constraint constructors."""
        return {
            'adaptive_lambda_cai': self.adaptive_lambda_cai,
            'lambda_cai_lr': self.lambda_cai_lr,
            'lambda_cai_max': self.lambda_cai_max,
            'cai_tolerance': self.cai_tolerance,
            'smoothing_factor': self.smoothing_factor,
        }


@dataclass 
class ConstraintConfig:
    """Master configuration for all constraint types."""
    
    # Basic CAI settings
    enable_cai: bool = False
    cai_target: float = 0.8
    cai_weight: float = 0.1  # Initial/fixed lambda_cai value
    species: str = 'ecoli_bl21de3'
    
    # Adaptive lambda_cai configuration
    adaptive_config: AdaptiveLambdaCAIConfig = field(default_factory=AdaptiveLambdaCAIConfig)
    
    # Constraint-specific settings
    constraint_specific: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    
    # Experimental features
    experimental_features: Dict[str, bool] = field(default_factory=dict)
    
    def __post_init__(self):
        """Initialize constraint-specific defaults."""
        if not self.constraint_specific:
            self.constraint_specific = {
                'ams': {
                    'batch_size': 1,
                    'device': 'cuda'
                },
                'cpc': {
                    'batch_size': 1,
                    'device': 'cuda'
                },
                'lagrangian': {
                    'batch_size': 1,
                    'initial_lambda': 0.01,
                    'adaptive_lambda': True,
                    'lambda_lr': 0.01,
                    'lambda_max': 10.0,
                    'device': None
                }
            }
        
        if not self.experimental_features:
            self.experimental_features = {
                'enhanced_convergence_tracking': False,
                'advanced_stability_detection': False,
                'multi_objective_adaptation': False
            }
    
    def get_constraint_kwargs(self, constraint_type: str) -> Dict[str, Any]:
        """
        Get complete kwargs for initializing a specific constraint type.
        
        Args:
            constraint_type: One of 'ams', 'cpc', 'lagrangian'
            
        Returns:
            Dictionary of keyword arguments for constraint constructor
        """
        if constraint_type not in ['ams', 'cpc', 'lagrangian']:
            raise ValueError(f"Unknown constraint type: {constraint_type}")
        
        # Base CAI configuration
        kwargs = {
            'enable_cai': self.enable_cai,
            'cai_target': self.cai_target,
            'cai_weight': self.cai_weight,
        }
        
        # Add adaptive lambda_cai configuration if CAI is enabled
        if self.enable_cai:
            kwargs.update(self.adaptive_config.to_dict())
        
        # Add constraint-specific settings
        kwargs.update(self.constraint_specific.get(constraint_type, {}))
        
        # Add experimental features
        for feature, enabled in self.experimental_features.items():
            if enabled:
                kwargs[feature] = enabled
        
        return kwargs


# Predefined configurations for common use cases
class PresetConfigs:
    """Predefined configuration presets for common scenarios."""
    
    @staticmethod
    def default_fixed() -> ConstraintConfig:
        """Default configuration with fixed lambda_cai (backward compatible)."""
        return ConstraintConfig()
    
    @staticmethod
    def adaptive_conservative() -> ConstraintConfig:
        """Conservative adaptive configuration with slow adaptation."""
        adaptive_config = AdaptiveLambdaCAIConfig(
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.05,  # Slow learning
            lambda_cai_max=1.5,  # Conservative maximum
            cai_tolerance=0.03,  # Tight tolerance
            smoothing_factor=0.95  # Strong smoothing
        )
        
        return ConstraintConfig(
            enable_cai=True,
            adaptive_config=adaptive_config
        )
    
    @staticmethod
    def adaptive_aggressive() -> ConstraintConfig:
        """Aggressive adaptive configuration with fast adaptation."""
        adaptive_config = AdaptiveLambdaCAIConfig(
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.2,  # Fast learning
            lambda_cai_max=3.0,  # High maximum
            cai_tolerance=0.08,  # Relaxed tolerance
            smoothing_factor=0.8  # Less smoothing
        )
        
        return ConstraintConfig(
            enable_cai=True,
            adaptive_config=adaptive_config
        )
    
    @staticmethod
    def research_mode() -> ConstraintConfig:
        """Research configuration with all experimental features enabled."""
        adaptive_config = AdaptiveLambdaCAIConfig(
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.1,
            lambda_cai_max=2.0,
            cai_tolerance=0.05,
        )
        
        experimental_features = {
            'enhanced_convergence_tracking': True,
            'advanced_stability_detection': True,
            'multi_objective_adaptation': True
        }
        
        return ConstraintConfig(
            enable_cai=True,
            adaptive_config=adaptive_config,
            experimental_features=experimental_features
        )


def load_config_from_dict(config_dict: Dict[str, Any]) -> ConstraintConfig:
    """
    Load configuration from a dictionary (e.g., from JSON/YAML file).
    
    Args:
        config_dict: Configuration dictionary
        
    Returns:
        ConstraintConfig object
    """
    # Extract adaptive config if present
    adaptive_dict = config_dict.pop('adaptive_config', {})
    adaptive_config = AdaptiveLambdaCAIConfig(**adaptive_dict)
    
    # Create main config
    config = ConstraintConfig(adaptive_config=adaptive_config, **config_dict)
    
    return config


def get_backward_compatible_kwargs(**legacy_kwargs) -> Dict[str, Any]:
    """
    Convert legacy keyword arguments to new format for backward compatibility.
    
    This function ensures existing code continues to work while providing
    access to new adaptive features.
    
    Args:
        **legacy_kwargs: Legacy keyword arguments
        
    Returns:
        Updated kwargs with new parameter names and defaults
    """
    # Map old parameter names to new ones
    param_mapping = {
        'cai_lambda': 'cai_weight',  # Old name for lambda_cai
        'adaptive_cai': 'adaptive_lambda_cai',  # Old adaptive flag
        'cai_lr': 'lambda_cai_lr',  # Old learning rate name
        'max_cai_weight': 'lambda_cai_max',  # Old maximum name
    }
    
    # Apply mapping
    updated_kwargs = {}
    for old_name, new_name in param_mapping.items():
        if old_name in legacy_kwargs:
            updated_kwargs[new_name] = legacy_kwargs.pop(old_name)
    
    # Keep remaining kwargs as-is
    updated_kwargs.update(legacy_kwargs)
    
    # Apply defaults for backward compatibility
    if 'adaptive_lambda_cai' not in updated_kwargs:
        updated_kwargs['adaptive_lambda_cai'] = False  # Default: disabled
    
    return updated_kwargs


# Example usage configurations
EXAMPLE_CONFIGS = {
    'paper_reproduction': ConstraintConfig(
        enable_cai=True,
        cai_target=0.8,
        cai_weight=0.1,
        adaptive_config=AdaptiveLambdaCAIConfig(adaptive_lambda_cai=False)
    ),
    
    'adaptive_optimization': ConstraintConfig(
        enable_cai=True,
        cai_target=0.8,
        cai_weight=0.1,
        adaptive_config=AdaptiveLambdaCAIConfig(
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.1,
            lambda_cai_max=2.0,
            cai_tolerance=0.05
        )
    ),
    
    'high_cai_target': ConstraintConfig(
        enable_cai=True,
        cai_target=0.9,  # Higher target
        cai_weight=0.15,  # Higher initial weight
        adaptive_config=AdaptiveLambdaCAIConfig(
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.15,  # Faster adaptation for high targets
            lambda_cai_max=3.0,  # Higher maximum
            cai_tolerance=0.03   # Tighter tolerance
        )
    )
}