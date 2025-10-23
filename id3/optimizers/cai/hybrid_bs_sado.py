"""
Hybrid Binary Search-SADO Optimizer

Combines Binary Search and SADO algorithms with intelligent switching based on
sequence repetition detection to avoid local optima.
"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, Tuple, Optional, Set
from ..base import BaseCAIOptimizer
from .binary_search import BinarySearchCAIOptimizer
from .sado import SADOOptimizer

logger = logging.getLogger(__name__)


class HybridBSSADOOptimizer(BaseCAIOptimizer):
    """

    





    """
    
    def __init__(self, 
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None,
                 repetition_threshold: int = 3):
        """
        Initialize the Hybrid optimizer.

        Args:
            species: Species name for loading CAI weights
            device: PyTorch device for tensor operations
            amino_acid_sequence: Amino acid sequence for optimization
            repetition_threshold: Number of consecutive repetitions before switching to SADO
        """
        super().__init__(species, device, amino_acid_sequence)

        # Initialize Binary Search optimizer
        self.bs_optimizer = BinarySearchCAIOptimizer(
            species=species,
            device=device,
            amino_acid_sequence=amino_acid_sequence
        )
        # Initialize SADO optimizer
        self.sado_optimizer = SADOOptimizer(
            species=species,
            device=device,
            amino_acid_sequence=amino_acid_sequence
        )

        # Repetition detection state
        self.sequence_history: Set[str] = set()
        self.recent_sequences = []
        self.repetition_threshold = repetition_threshold
        self.repetition_count = 0
        self.last_sequence_hash = None

        # Statistics tracking
        self.bs_calls = 0
        self.sado_calls = 0
        self.switch_count = 0
        
        logger.info(f"HybridBSSADOOptimizer initialized with threshold={repetition_threshold}")
    
    def _compute_sequence_hash(self, indices: torch.Tensor) -> str:
        """Compute hash of codon sequence for duplicate detection."""
        if indices.dim() > 1:
            indices = indices.argmax(dim=-1)
        indices_np = indices.cpu().numpy()
        return hashlib.md5(indices_np.tobytes()).hexdigest()
    
    def _check_repetition(self, sequence_hash: str) -> bool:
        """Check if sequence is repetitive."""
        # Check consecutive repetition
        if sequence_hash == self.last_sequence_hash:
            self.repetition_count += 1
        else:
            self.repetition_count = 0
            self.last_sequence_hash = sequence_hash

        # Check if hash appears in history
        is_duplicate = sequence_hash in self.sequence_history

        # Add to history
        self.sequence_history.add(sequence_hash)
        self.recent_sequences.append(sequence_hash)

        # Keep only recent sequences (sliding window)
        if len(self.recent_sequences) > 10:
            self.recent_sequences.pop(0)
        
        return is_duplicate or (self.repetition_count >= self.repetition_threshold)
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float = 0.8,
                 amino_acid_sequence: Optional[str] = None,
                 valid_codon_mask: Optional[torch.Tensor] = None,
                 **kwargs) -> Tuple[torch.Tensor, Dict]:
        """
        Execute optimization with automatic Binary Search to SADO switching.

        Args:
            pi_accessibility: RNA accessibility probability distribution
            target_cai: Target CAI value
            amino_acid_sequence: Amino acid sequence to optimize
            valid_codon_mask: Mask indicating valid codons
            **kwargs: Additional parameters

        Returns:
            Tuple of (optimized_distribution, metadata)
        """
        if amino_acid_sequence is None:
            amino_acid_sequence = self.amino_acid_sequence

        # Track Binary Search usage
        self.bs_calls += 1

        try:
            # First attempt: Binary Search
            bs_result, bs_metadata = self.bs_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )

            # Get gamma value from Binary Search
            bs_gamma = bs_metadata.get('optimal_gamma', 0.5)

            # Convert to indices for hash computation
            if bs_result.dim() > 1:
                indices = bs_result.argmax(dim=-1)
            else:
                indices = bs_result
                
            seq_hash = self._compute_sequence_hash(indices)

            # Check for repetition and switch to SADO if detected
            if self._check_repetition(seq_hash):
                logger.info(f"Binary Search produced repetitive sequence, switching to SADO (gamma={bs_gamma:.3f})")
                self.switch_count += 1
                self.sado_calls += 1

                # Prepare SADO kwargs (remove valid_codon_mask as SADO doesn't use it)
                # Use gamma from Binary Search as initialization
                sado_kwargs = {k: v for k, v in kwargs.items() if k != 'valid_codon_mask'}
                sado_kwargs['gamma'] = bs_gamma
                
                sado_result, sado_metadata = self.sado_optimizer.optimize(
                    pi_accessibility=pi_accessibility,
                    target_cai=target_cai,
                    amino_acid_sequence=amino_acid_sequence,
                    **sado_kwargs
                )

                # Add hybrid method information to metadata
                sado_metadata['method'] = 'hybrid_bs_sado'
                sado_metadata['switched_to_sado'] = True
                sado_metadata['bs_gamma'] = bs_gamma
                sado_metadata['switch_reason'] = 'repetition_detected'
                sado_metadata['bs_calls'] = self.bs_calls
                sado_metadata['sado_calls'] = self.sado_calls
                sado_metadata['switch_count'] = self.switch_count
                
                return sado_result, sado_metadata

            else:
                # Binary Search successful, no switching needed
                bs_metadata['method'] = 'hybrid_bs_sado'
                bs_metadata['switched_to_sado'] = False
                bs_metadata['bs_calls'] = self.bs_calls
                bs_metadata['sado_calls'] = self.sado_calls
                bs_metadata['switch_count'] = self.switch_count
                
                return bs_result, bs_metadata
                
        except Exception as e:
            logger.error(f"Binary Search failed: {e}, switching to SADO")
            self.switch_count += 1
            self.sado_calls += 1

            # Fallback to SADO with default gamma
            sado_kwargs = {k: v for k, v in kwargs.items() if k != 'valid_codon_mask'}
            sado_kwargs['gamma'] = 0.5
            
            sado_result, sado_metadata = self.sado_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                **sado_kwargs
            )
            
            sado_metadata['method'] = 'hybrid_bs_sado'
            sado_metadata['switched_to_sado'] = True
            sado_metadata['switch_reason'] = 'bs_error'
            sado_metadata['error'] = str(e)
            sado_metadata['bs_calls'] = self.bs_calls
            sado_metadata['sado_calls'] = self.sado_calls
            sado_metadata['switch_count'] = self.switch_count
            
            return sado_result, sado_metadata
    
    def reset(self):
        """Reset optimizer state."""
        self.sequence_history.clear()
        self.recent_sequences.clear()
        self.repetition_count = 0
        self.last_sequence_hash = None
        self.bs_calls = 0
        self.sado_calls = 0
        self.switch_count = 0

        # Reset sub-optimizers
        if hasattr(self.bs_optimizer, 'reset'):
            self.bs_optimizer.reset()
        if hasattr(self.sado_optimizer, 'reset'):
            self.sado_optimizer.reset()
    
    def get_statistics(self) -> Dict:
        """Get optimizer statistics."""
        return {
            'method': 'hybrid_bs_sado',
            'bs_calls': self.bs_calls,
            'sado_calls': self.sado_calls,
            'switch_count': self.switch_count,
            'switch_rate': self.switch_count / max(1, self.bs_calls),
            'unique_sequences': len(self.sequence_history),
            'repetition_threshold': self.repetition_threshold
        }