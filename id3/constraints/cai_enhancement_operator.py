"""
CAI Enhancement Operator - Router Implementation







Mathematical Formulation from Paper:
    Ψ_{π→τ}(π, τ_CAI) = arg max_S P(S|π) s.t. CAI(S) ≥ τ_CAI
"""

import torch
import logging
from typing import Dict, Tuple, Optional, Any
from pathlib import Path


from id3.optimizers.cai import BinarySearchCAIOptimizer

logger = logging.getLogger(__name__)


class CAIEnhancementOperator:
    """

    





    


        operator = CAIEnhancementOperator()
        

        operator = CAIEnhancementOperator(method='binary_search')
        

        operator = CAIEnhancementOperator(method='binary_search')
    """
    

    OPTIMIZERS = {
        'binary_search': BinarySearchCAIOptimizer,
        'incremental': None,
    }
    

    _amino_acid_sequence_cache = {}
    _amino_acid_weights_cache = {}
    
    def __init__(self, 
                 method: str = 'incremental',
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        """

        
        Args:








        """
        if method not in self.OPTIMIZERS:
            raise ValueError(
                f"Unknown optimization method: {method}. "
                f"Available methods: {list(self.OPTIMIZERS.keys())}"
            )
        
        self.method = method
        self.species = species
        self.device = device if device is not None else torch.device('cpu')
        self.amino_acid_sequence = amino_acid_sequence
        

        if method == 'incremental':

            from id3.optimizers.cai.incremental import IncrementalCAIOptimizer
            self.optimizer = IncrementalCAIOptimizer(
                species=species,
                device=self.device,
                amino_acid_sequence=amino_acid_sequence
            )
        else:
            optimizer_class = self.OPTIMIZERS[method]
            self.optimizer = optimizer_class(
                species=species,
                device=self.device,
                amino_acid_sequence=amino_acid_sequence
            )
        

        if hasattr(self.optimizer, 'wi_table'):
            self.wi_table = self.optimizer.wi_table
        if hasattr(self.optimizer, 'weights_tensor'):
            self.weights_tensor = self.optimizer.weights_tensor
        

        self._sync_cache_state()
        
        logger.info(f"CAIEnhancementOperator initialized with method: {method}")
    
    def _sync_cache_state(self):

        if hasattr(self.optimizer, 'max_achievable_cai'):
            self.max_achievable_cai = self.optimizer.max_achievable_cai
        else:
            self.max_achievable_cai = None
            
        if hasattr(self.optimizer, 'cai_optimal_distribution'):
            self.cai_optimal_distribution = self.optimizer.cai_optimal_distribution
        else:
            self.cai_optimal_distribution = None
            
        if hasattr(self.optimizer, 'cached_amino_acid_sequence'):
            self.cached_amino_acid_sequence = self.optimizer.cached_amino_acid_sequence
        else:
            self.cached_amino_acid_sequence = None
    
    def apply_cai_enhancement(self,
                             pi_accessibility: torch.Tensor,
                             amino_acid_sequence: str,
                             valid_codon_mask: torch.Tensor,
                             codon_indices: torch.Tensor,
                             target_cai: float,
                             **kwargs) -> Tuple[torch.Tensor, Dict]:
        """
        Apply CAI enhancement operator Ψ_{π→τ}

        This is the main entry point for the Ψ_{π→τ} operator from the paper.
        Routes to the appropriate optimizer based on the method selected during initialization.

        Args:
            pi_accessibility: Accessibility-optimized distribution π
            amino_acid_sequence: Target amino acid sequence
            valid_codon_mask: Valid codon mask
            codon_indices: Codon indices (may be ignored by some optimizers)
            target_cai: Target CAI value τ_CAI
            **kwargs: Other method-specific parameters

        Returns:
            Tuple of (discrete distribution, metadata dictionary)
        """

        # All optimizers use the same interface
        discrete_distribution, metadata = self.optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            amino_acid_sequence=amino_acid_sequence,
            valid_codon_mask=valid_codon_mask,
            **kwargs
            )
        

        self._sync_cache_state()
        

        if 'method' not in metadata:
            metadata['method'] = self.method
        
        return discrete_distribution, metadata
    
    def enhance(self,
                pi_accessibility: torch.Tensor,
                target_cai: float = 0.8,
                **kwargs) -> Tuple[torch.Tensor, Dict]:
        """
        Simplified enhancement interface (backward compatible)

        Args:
            pi_accessibility: Accessibility-optimized distribution
            target_cai: Target CAI value
            **kwargs: Other parameters

        Returns:
            Tuple of (optimized distribution, metadata)
        """
        return self.optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            **kwargs
        )
    
    def reset(self):
        """
        Reset optimizer state

        Some optimizers maintain internal state that needs to be reset
        when starting a new optimization sequence.
        """
        if hasattr(self.optimizer, 'reset'):
            self.optimizer.reset()
            logger.debug(f"Optimizer {self.method} state reset")
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        Get optimizer performance statistics

        Returns:
            Dictionary containing performance metrics
        """
        stats = {
            'method': self.method,
            'species': self.species,
        }
        

        if hasattr(self.optimizer, 'get_statistics'):
            optimizer_stats = self.optimizer.get_statistics()
            stats.update(optimizer_stats)
        

        if hasattr(self.optimizer, 'get_diversity_stats'):
            diversity_stats = self.optimizer.get_diversity_stats()
            stats['diversity'] = diversity_stats
        
        return stats
    
    def switch_method(self, new_method: str):
        """
        Dynamically switch optimization method

        Args:
            new_method: New optimization method name
        """
        if new_method not in self.OPTIMIZERS:
            raise ValueError(f"Unknown method: {new_method}")
        
        if new_method != self.method:
            logger.info(f"Switching optimization method from {self.method} to {new_method}")
            

            optimizer_class = self.OPTIMIZERS[new_method]
            self.optimizer = optimizer_class(
                species=self.species,
                device=self.device,
                amino_acid_sequence=self.amino_acid_sequence
            )
            
            self.method = new_method
            self._sync_cache_state()
    

    
    def _load_or_compute_amino_acid_cache(self, amino_acid_sequence: str):
        """Backward compatible: load or compute amino acid sequence cache"""
        if hasattr(self.optimizer, '_load_or_compute_amino_acid_cache'):
            self.optimizer._load_or_compute_amino_acid_cache(amino_acid_sequence)
            self._sync_cache_state()
    
    def _precompute_amino_acid_weights(self):

        if hasattr(self.optimizer, '_precompute_amino_acid_weights'):
            self.optimizer._precompute_amino_acid_weights()
    
    def discretize_distribution(self, distribution: torch.Tensor, valid_mask: torch.Tensor) -> torch.Tensor:
        """Backward compatible: discretize distribution"""
        from id3.optimizers.cai.utils import discretize_distribution
        return discretize_distribution(distribution, valid_mask)
    
    def compute_discrete_cai(self, 
                            discrete_dist: torch.Tensor,
                            amino_acid_sequence: str,
                            valid_codon_mask: torch.Tensor,
                            codon_indices: torch.Tensor) -> float:

        if hasattr(self.optimizer, '_compute_cai_from_indices'):
            indices = discrete_dist.argmax(dim=-1)
            return self.optimizer._compute_cai_from_indices(indices, amino_acid_sequence)
        return 0.0
    
    def interpolate_distributions(self,
                                 dist1: torch.Tensor,
                                 dist2: torch.Tensor,
                                 gamma: float) -> torch.Tensor:
        """Backward compatible: distribution interpolation"""
        from id3.optimizers.cai.utils import interpolate_distributions
        return interpolate_distributions(dist1, dist2, gamma)
    
    @classmethod
    def available_methods(cls) -> list:

        return list(cls.OPTIMIZERS.keys())
    
    def __repr__(self) -> str:
        return (
            f"CAIEnhancementOperator(method='{self.method}', "
            f"species='{self.species}', device={self.device})"
        )