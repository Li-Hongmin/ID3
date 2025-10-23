"""
ID3-CAI Integration Module (KISS Simplified)

Simplified CAI integration for ID3 framework following KISS principles.
Removed complex binary search optimizer that broke amino acid constraints.

Key Components:
- CAIModule: Basic CAI calculation
- Enhanced constraints with post-processing CAI optimization
- Simple validation and computation functions

Usage Examples:
    # Basic CAI calculation
    from id3.cai import CAIModule
    cai_calculator = CAIModule()
    cai_score = cai_calculator.compute_cai(dna_sequence)
    
    # CAI-enhanced constraints (built-in)
    constraint = CodonProfileConstraint(protein_seq, enable_cai=True, cai_target=0.8)
    result = constraint.forward()
    sequence = constraint.get_discrete_sequence()  # CAI-optimized but amino-acid preserving
"""

from .module import CAIModule

# Import available core components
from .probability import rna_to_codon_probabilities
from .validator import compute_cai_from_sequence
from .unified_enhancer import UnifiedCAIEnhancer

# All components are now simplified and available
CAI_INTEGRATION_AVAILABLE = True

__all__ = [
    'CAIModule',
    'rna_to_codon_probabilities',
    'compute_cai_from_sequence',
    'UnifiedCAIEnhancer'
]