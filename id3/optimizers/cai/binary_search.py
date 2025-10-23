"""
Binary Search CAI Optimizer

Implements a binary search-based approach for finding the optimal interpolation parameter
to satisfy CAI constraints while maximizing RNA accessibility probability.
"""

import torch
import logging
from typing import Dict, Tuple, Optional, Any
from pathlib import Path
import json

from ..base import BaseCAIOptimizer
from .utils import (
    load_cai_weights, 
    compute_cai,
    discretize_distribution,
    interpolate_distributions
)
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


class BinarySearchCAIOptimizer(BaseCAIOptimizer):
    """

    


    



    """
    

    _amino_acid_sequence_cache = {}
    _amino_acid_weights_cache = {}
    
    def __init__(self, 
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        """
        Initialize the Binary Search CAI Optimizer.

        Args:
            species: Species name for loading CAI weights
            device: PyTorch device for tensor operations
            amino_acid_sequence: Amino acid sequence for optimization (optional)
        """
        super().__init__(species, device, amino_acid_sequence)
        

        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)
        

        self._precompute_amino_acid_weights()
        

        self.max_achievable_cai = None
        self.cai_optimal_distribution = None
        self.cached_amino_acid_sequence = None
        

        if amino_acid_sequence is not None:
            self._load_or_compute_amino_acid_cache(amino_acid_sequence)
    
    def _precompute_amino_acid_weights(self):
        """Precompute and normalize weights for each amino acid's codons."""
        for aa in amino_acids_to_codons:
            cache_key = f"{self.species}_{aa}"
            
            if cache_key not in self._amino_acid_weights_cache:
                codons = amino_acids_to_codons[aa]
                weights = torch.zeros(len(codons))
                
                for i, codon in enumerate(codons):
                    weights[i] = self.wi_table.get(codon, 0.0)

                # Normalize weights (divide by maximum weight)
                max_weight = weights.max()
                if max_weight > 0:
                    weights = weights / max_weight
                
                self._amino_acid_weights_cache[cache_key] = weights.to(self.device)
    
    def _load_or_compute_amino_acid_cache(self, amino_acid_sequence: str):
        """Load or compute cached data for amino acid sequence."""
        cache_key = f"{self.species}_{amino_acid_sequence}"
        
        if cache_key in self._amino_acid_sequence_cache:
            # Load from cache
            cache_data = self._amino_acid_sequence_cache[cache_key]
            self.max_achievable_cai = cache_data['max_achievable_cai']
            self.cai_optimal_distribution = cache_data['cai_optimal_distribution']
            self.cached_amino_acid_sequence = amino_acid_sequence
            logger.debug(f"Loaded cached data for sequence (max CAI: {self.max_achievable_cai:.4f})")
        else:
            # Compute and cache
            self._compute_and_cache_sequence_data(amino_acid_sequence)
    
    def _compute_and_cache_sequence_data(self, amino_acid_sequence: str):
        """Compute and cache maximum achievable CAI and optimal distribution for a sequence."""
        # Build CAI-optimal distribution (select highest weight codon for each position)
        cai_optimal_dist = []
        # Get maximum codon count across all amino acids
        max_codons = max(len(amino_acids_to_codons[aa]) for aa in amino_acids_to_codons)

        for aa in amino_acid_sequence:
            if aa in amino_acids_to_codons:
                cache_key = f"{self.species}_{aa}"
                weights = self._amino_acid_weights_cache[cache_key]
                # Pad weights to match max_codons dimension
                if len(weights) < max_codons:
                    padded = torch.zeros(max_codons, device=self.device)
                    padded[:len(weights)] = weights
                    cai_optimal_dist.append(padded)
                else:
                    cai_optimal_dist.append(weights[:max_codons])
        
        if cai_optimal_dist:
            self.cai_optimal_distribution = torch.stack(cai_optimal_dist)

            # Compute maximum achievable CAI (using codons with highest weights)
            max_indices = self.cai_optimal_distribution.argmax(dim=-1)
            self.max_achievable_cai = self._compute_cai_from_indices(
                max_indices, amino_acid_sequence
            )

            # Cache the computed data
            cache_key = f"{self.species}_{amino_acid_sequence}"
            self._amino_acid_sequence_cache[cache_key] = {
                'max_achievable_cai': self.max_achievable_cai,
                'cai_optimal_distribution': self.cai_optimal_distribution
            }
            
            self.cached_amino_acid_sequence = amino_acid_sequence
            logger.debug(f"Computed and cached data for sequence (max CAI: {self.max_achievable_cai:.4f})")
    
    def _compute_cai_from_indices(self, indices: torch.Tensor, amino_acid_sequence: str) -> float:
        """Compute CAI value from codon indices."""
        cai_product = 1.0
        
        for pos, (idx, aa) in enumerate(zip(indices, amino_acid_sequence)):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                if idx < len(codons):
                    codon = codons[idx.item() if torch.is_tensor(idx) else idx]
                    weight = self.wi_table.get(codon, 0.0)
                    cai_product *= weight

        # CAI is the geometric mean of weights
        cai = cai_product ** (1.0 / len(amino_acid_sequence))
        return cai
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float = 0.8,
                 amino_acid_sequence: str = None,
                 valid_codon_mask: torch.Tensor = None,
                 max_iterations: int = 20,
                 tolerance: float = 1e-4,
                 **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Execute binary search optimization to find optimal gamma parameter.

        Args:
            pi_accessibility: RNA accessibility probability distribution
            target_cai: Target CAI value to achieve
            amino_acid_sequence: Amino acid sequence to optimize
            valid_codon_mask: Mask indicating valid codons at each position
            max_iterations: Maximum number of binary search iterations
            tolerance: Convergence tolerance for binary search

        Returns:
            Tuple of (discrete_distribution, metadata)
        """
        # Handle 3D tensors by squeezing batch dimension
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)

        # Move tensor to the correct device
        pi_accessibility = pi_accessibility.to(self.device)

        # Load or compute amino acid cache if sequence changed
        if amino_acid_sequence and amino_acid_sequence != self.cached_amino_acid_sequence:
            self._load_or_compute_amino_acid_cache(amino_acid_sequence)

        # Ensure we have CAI optimal distribution
        if self.cai_optimal_distribution is None:
            raise ValueError("No amino acid sequence provided for CAI optimization")
        
        w_cai_optimal = self.cai_optimal_distribution

        # Pad CAI optimal distribution if necessary
        if w_cai_optimal.shape[1] < pi_accessibility.shape[1]:
            padding = torch.zeros(
                w_cai_optimal.shape[0], 
                pi_accessibility.shape[1] - w_cai_optimal.shape[1],
                device=self.device
            )
            w_cai_optimal = torch.cat([w_cai_optimal, padding], dim=1)

        # Move mask to device if provided
        if valid_codon_mask is not None:
            valid_codon_mask = valid_codon_mask.to(self.device)

        # Execute binary search to find optimal gamma
        optimal_gamma, discrete_dist = self._binary_search_gamma(
            pi_accessibility, w_cai_optimal, amino_acid_sequence,
            valid_codon_mask, target_cai, max_iterations, tolerance
        )

        # Compute final CAI from discrete distribution
        final_indices = discrete_dist.argmax(dim=-1)
        final_cai = self._compute_cai_from_indices(final_indices, amino_acid_sequence)

        # Build metadata dictionary
        metadata = {
            'optimal_gamma': optimal_gamma,
            'final_cai': final_cai,
            'target_cai': target_cai,
            'max_achievable_cai': self.max_achievable_cai,
            'constraint_satisfied': final_cai >= target_cai,
            'method': 'binary_search'
        }
        
        return discrete_dist, metadata
    
    def _binary_search_gamma(self,
                            pi_accessibility: torch.Tensor,
                            w_cai_optimal: torch.Tensor,
                            amino_acid_sequence: str,
                            valid_codon_mask: torch.Tensor,
                            target_cai: float,
                            max_iterations: int,
                            tolerance: float) -> Tuple[float, torch.Tensor]:
        """
        Execute binary search to find optimal gamma parameter.

        Returns optimal gamma and corresponding discrete distribution.
        """
        # Test gamma=0 (pure accessibility-optimal)
        discrete_dist_gamma0 = discretize_distribution(pi_accessibility, valid_codon_mask)
        cai_gamma0 = self._compute_cai_from_indices(
            discrete_dist_gamma0.argmax(dim=-1), amino_acid_sequence
        )

        # Early stop if gamma=0 already satisfies CAI constraint
        if cai_gamma0 >= target_cai:
            logger.debug(f"Early stop: gamma=0 satisfies CAI ({cai_gamma0:.4f} >= {target_cai:.4f})")
            return 0.0, discrete_dist_gamma0

        # Check if target CAI is achievable
        if self.max_achievable_cai < target_cai:
            logger.warning(f"Target CAI ({target_cai:.4f}) exceeds maximum ({self.max_achievable_cai:.4f})")
            return 1.0, discretize_distribution(w_cai_optimal, valid_codon_mask)

        # Initialize binary search bounds
        gamma_low = 0.0
        gamma_high = 1.0
        best_gamma = 1.0
        best_discrete_dist = discretize_distribution(w_cai_optimal, valid_codon_mask)

        # Linear interpolation-based initial guess
        if self.max_achievable_cai > cai_gamma0:
            estimated_gamma = (target_cai - cai_gamma0) / (self.max_achievable_cai - cai_gamma0)
            estimated_gamma = max(0.0, min(1.0, estimated_gamma))

            # Test estimated gamma
            interpolated = interpolate_distributions(pi_accessibility, w_cai_optimal, estimated_gamma)
            discrete_dist = discretize_distribution(interpolated, valid_codon_mask)
            estimated_cai = self._compute_cai_from_indices(
                discrete_dist.argmax(dim=-1), amino_acid_sequence
            )
            
            if estimated_cai >= target_cai:
                gamma_high = estimated_gamma
                best_gamma = estimated_gamma
                best_discrete_dist = discrete_dist
            else:
                gamma_low = estimated_gamma

        # Binary search loop
        for iteration in range(max_iterations):
            if abs(gamma_high - gamma_low) < tolerance:
                break
            
            gamma_mid = (gamma_low + gamma_high) / 2.0

            # Interpolate distributions at gamma_mid
            interpolated = interpolate_distributions(pi_accessibility, w_cai_optimal, gamma_mid)
            discrete_dist = discretize_distribution(interpolated, valid_codon_mask)

            # Compute CAI from discretized distribution
            discrete_cai = self._compute_cai_from_indices(
                discrete_dist.argmax(dim=-1), amino_acid_sequence
            )

            # Update search bounds based on CAI constraint
            if discrete_cai >= target_cai:
                gamma_high = gamma_mid
                best_gamma = gamma_mid
                best_discrete_dist = discrete_dist

                # Early termination if we're close to lower bound
                if gamma_mid - gamma_low < tolerance * 2:
                    break
            else:
                gamma_low = gamma_mid
        
        logger.debug(f"Binary search completed (gamma={best_gamma:.4f})")
        return best_gamma, best_discrete_dist