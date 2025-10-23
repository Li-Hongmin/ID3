"""
CAI Optimizers

Provides discrete optimization algorithms for CAI enhancement while
maintaining amino acid sequence constraints.

Available optimizers:
- BinarySearchCAIOptimizer: Binary search for optimal CAI factor (default)
- IncrementalCAIOptimizer: Incremental optimization with caching
"""

from .binary_search import BinarySearchCAIOptimizer

# Default optimizer
DefaultCAIOptimizer = BinarySearchCAIOptimizer

__all__ = ['BinarySearchCAIOptimizer', 'DefaultCAIOptimizer']
