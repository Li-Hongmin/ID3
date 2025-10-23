"""
CAI Optimization Configuration

This module provides centralized configuration for CAI optimization methods.
Default method is now Discrete + Incremental due to 22x performance improvement.
"""

from dataclasses import dataclass
from typing import Optional


@dataclass
class CAIOptimizationConfig:
    """
    Configuration for CAI optimization methods.
    
    Attributes:
        use_discrete_search: Use discrete binary search (default: True)
        use_incremental_calc: Use incremental CAI calculation (default: True)
        cai_target: Target CAI value (default: 0.8)
        verbose: Print debug information (default: False)
        cache_size: Cache size for incremental calculator (default: 100)
        use_symmetry_breaking: Apply perturbation for simultaneous switches (default: True)
    """
    
    # Optimization method selection
    use_discrete_search: bool = True  # NEW DEFAULT: Discrete search
    use_incremental_calc: bool = True  # NEW DEFAULT: Incremental calculation
    
    # CAI parameters
    cai_target: float = 0.8
    max_iterations: int = 50
    precision: float = 1e-6
    
    # Debug and monitoring
    verbose: bool = False
    track_performance: bool = False
    
    # Advanced settings
    cache_size: int = 100
    use_symmetry_breaking: bool = True
    perturbation_magnitude: float = 1e-10
    
    @classmethod
    def default(cls):
        """Get default configuration with optimal settings."""
        return cls()
    
    @classmethod
    def fast(cls):
        """Get configuration optimized for speed (discrete + incremental)."""
        return cls(
            use_discrete_search=True,
            use_incremental_calc=True,
            verbose=False
        )
    
    @classmethod
    def accurate(cls):
        """Get configuration optimized for accuracy (continuous search)."""
        return cls(
            use_discrete_search=False,
            use_incremental_calc=True,  # Still use incremental for speed
            precision=1e-8,
            max_iterations=100
        )
    
    @classmethod
    def debug(cls):
        """Get configuration for debugging with verbose output."""
        return cls(
            use_discrete_search=True,
            use_incremental_calc=True,
            verbose=True,
            track_performance=True
        )
    
    @classmethod
    def legacy(cls):
        """Get legacy configuration (original continuous method)."""
        return cls(
            use_discrete_search=False,
            use_incremental_calc=False,
            verbose=False
        )
    
    def get_method_description(self) -> str:
        """Get human-readable description of the optimization method."""
        if self.use_discrete_search and self.use_incremental_calc:
            return "Discrete Binary Search + Incremental CAI (Default, ~22x faster)"
        elif self.use_discrete_search:
            return "Discrete Binary Search (~7x faster)"
        elif self.use_incremental_calc:
            return "Continuous Search + Incremental CAI (~2x faster)"
        else:
            return "Original Continuous Search (Legacy)"
    
    def estimated_speedup(self) -> float:
        """Estimate speedup compared to original method."""
        if self.use_discrete_search and self.use_incremental_calc:
            return 22.0
        elif self.use_discrete_search:
            return 7.0
        elif self.use_incremental_calc:
            return 2.0
        else:
            return 1.0


# Global default configuration
DEFAULT_CAI_CONFIG = CAIOptimizationConfig.default()


def set_global_cai_method(method: str = "default"):
    """
    Set global CAI optimization method.
    
    Args:
        method: One of "default", "fast", "accurate", "debug", "legacy"
    """
    global DEFAULT_CAI_CONFIG
    
    method_map = {
        "default": CAIOptimizationConfig.default,
        "fast": CAIOptimizationConfig.fast,
        "accurate": CAIOptimizationConfig.accurate,
        "debug": CAIOptimizationConfig.debug,
        "legacy": CAIOptimizationConfig.legacy
    }
    
    if method not in method_map:
        raise ValueError(f"Unknown method: {method}. Choose from {list(method_map.keys())}")
    
    DEFAULT_CAI_CONFIG = method_map[method]()
    print(f"CAI optimization method set to: {DEFAULT_CAI_CONFIG.get_method_description()}")