"""




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

        
        Args:



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

        for aa in amino_acids_to_codons:
            cache_key = f"{self.species}_{aa}"
            
            if cache_key not in self._amino_acid_weights_cache:
                codons = amino_acids_to_codons[aa]
                weights = torch.zeros(len(codons))
                
                for i, codon in enumerate(codons):
                    weights[i] = self.wi_table.get(codon, 0.0)
                

                max_weight = weights.max()
                if max_weight > 0:
                    weights = weights / max_weight
                
                self._amino_acid_weights_cache[cache_key] = weights.to(self.device)
    
    def _load_or_compute_amino_acid_cache(self, amino_acid_sequence: str):
        """加载或计算氨基酸序列的缓存数据"""
        cache_key = f"{self.species}_{amino_acid_sequence}"
        
        if cache_key in self._amino_acid_sequence_cache:

            cache_data = self._amino_acid_sequence_cache[cache_key]
            self.max_achievable_cai = cache_data['max_achievable_cai']
            self.cai_optimal_distribution = cache_data['cai_optimal_distribution']
            self.cached_amino_acid_sequence = amino_acid_sequence
            logger.debug(f"Loaded cached data for sequence (max CAI: {self.max_achievable_cai:.4f})")
        else:

            self._compute_and_cache_sequence_data(amino_acid_sequence)
    
    def _compute_and_cache_sequence_data(self, amino_acid_sequence: str):


        cai_optimal_dist = []

        
        for aa in amino_acid_sequence:
            if aa in amino_acids_to_codons:
                cache_key = f"{self.species}_{aa}"
                weights = self._amino_acid_weights_cache[cache_key]

                if len(weights) < max_codons:
                    padded = torch.zeros(max_codons, device=self.device)
                    padded[:len(weights)] = weights
                    cai_optimal_dist.append(padded)
                else:
                    cai_optimal_dist.append(weights[:max_codons])
        
        if cai_optimal_dist:
            self.cai_optimal_distribution = torch.stack(cai_optimal_dist)
            

            max_indices = self.cai_optimal_distribution.argmax(dim=-1)
            self.max_achievable_cai = self._compute_cai_from_indices(
                max_indices, amino_acid_sequence
            )
            

            cache_key = f"{self.species}_{amino_acid_sequence}"
            self._amino_acid_sequence_cache[cache_key] = {
                'max_achievable_cai': self.max_achievable_cai,
                'cai_optimal_distribution': self.cai_optimal_distribution
            }
            
            self.cached_amino_acid_sequence = amino_acid_sequence
            logger.debug(f"Computed and cached data for sequence (max CAI: {self.max_achievable_cai:.4f})")
    
    def _compute_cai_from_indices(self, indices: torch.Tensor, amino_acid_sequence: str) -> float:
        """从密码子索引计算CAI值"""
        cai_product = 1.0
        
        for pos, (idx, aa) in enumerate(zip(indices, amino_acid_sequence)):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                if idx < len(codons):
                    codon = codons[idx.item() if torch.is_tensor(idx) else idx]
                    weight = self.wi_table.get(codon, 0.0)
                    cai_product *= weight
        

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

        
        Args:






            
        Returns:

        """

        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)
            

        pi_accessibility = pi_accessibility.to(self.device)
        

        if amino_acid_sequence and amino_acid_sequence != self.cached_amino_acid_sequence:
            self._load_or_compute_amino_acid_cache(amino_acid_sequence)
        

        if self.cai_optimal_distribution is None:
            raise ValueError("No amino acid sequence provided for CAI optimization")
        
        w_cai_optimal = self.cai_optimal_distribution
        

        if w_cai_optimal.shape[1] < pi_accessibility.shape[1]:
            padding = torch.zeros(
                w_cai_optimal.shape[0], 
                pi_accessibility.shape[1] - w_cai_optimal.shape[1],
                device=self.device
            )
            w_cai_optimal = torch.cat([w_cai_optimal, padding], dim=1)
        

        if valid_codon_mask is not None:
            valid_codon_mask = valid_codon_mask.to(self.device)
        

        optimal_gamma, discrete_dist = self._binary_search_gamma(
            pi_accessibility, w_cai_optimal, amino_acid_sequence,
            valid_codon_mask, target_cai, max_iterations, tolerance
        )
        

        final_indices = discrete_dist.argmax(dim=-1)
        final_cai = self._compute_cai_from_indices(final_indices, amino_acid_sequence)
        

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

        

        """

        discrete_dist_gamma0 = discretize_distribution(pi_accessibility, valid_codon_mask)
        cai_gamma0 = self._compute_cai_from_indices(
            discrete_dist_gamma0.argmax(dim=-1), amino_acid_sequence
        )
        

        if cai_gamma0 >= target_cai:
            logger.debug(f"Early stop: gamma=0 satisfies CAI ({cai_gamma0:.4f} >= {target_cai:.4f})")
            return 0.0, discrete_dist_gamma0
        

        if self.max_achievable_cai < target_cai:
            logger.warning(f"Target CAI ({target_cai:.4f}) exceeds maximum ({self.max_achievable_cai:.4f})")
            return 1.0, discretize_distribution(w_cai_optimal, valid_codon_mask)
        

        gamma_low = 0.0
        gamma_high = 1.0
        best_gamma = 1.0
        best_discrete_dist = discretize_distribution(w_cai_optimal, valid_codon_mask)
        

        if self.max_achievable_cai > cai_gamma0:
            estimated_gamma = (target_cai - cai_gamma0) / (self.max_achievable_cai - cai_gamma0)
            estimated_gamma = max(0.0, min(1.0, estimated_gamma))
            

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
        

        for iteration in range(max_iterations):
            if abs(gamma_high - gamma_low) < tolerance:
                break
            
            gamma_mid = (gamma_low + gamma_high) / 2.0
            

            interpolated = interpolate_distributions(pi_accessibility, w_cai_optimal, gamma_mid)
            discrete_dist = discretize_distribution(interpolated, valid_codon_mask)
            

            discrete_cai = self._compute_cai_from_indices(
                discrete_dist.argmax(dim=-1), amino_acid_sequence
            )
            

            if discrete_cai >= target_cai:
                gamma_high = gamma_mid
                best_gamma = gamma_mid
                best_discrete_dist = discrete_dist
                

                if gamma_mid - gamma_low < tolerance * 2:
                    break
            else:
                gamma_low = gamma_mid
        
        logger.debug(f"Binary search completed (gamma={best_gamma:.4f})")
        return best_gamma, best_discrete_dist