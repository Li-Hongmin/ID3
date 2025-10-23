#!/usr/bin/env python3
"""
Dimension Converter Module

Handles conversions between codon space (64-dimensional) and nucleotide space (4-dimensional)
for ID3 framework tensor operations.
"""

import torch
import torch.nn.functional as F
from typing import Tuple, Optional, Union
import numpy as np

class DimensionConverter:
    """
    Dimension converter between codon and nucleotide spaces.

    Converts between:
    - Codon space: 64-dimensional (all possible triplet combinations)
    - Nucleotide space: 4-dimensional (A, C, G, U/T)

    Maintains gradient flow for optimization and supports batch operations.
    """
    
    def __init__(self, device: torch.device = None):
        """
        Initialize dimension converter.

        Args:
            device: PyTorch device (defaults to CUDA if available, else CPU)
        """
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Create mapping matrices for codon ↔ nucleotide conversion
        self._create_codon_nucleotide_mapping()
        
    def _create_codon_nucleotide_mapping(self):
        """Create mapping matrix from 64 codons to nucleotide positions."""


        bases = ['A', 'C', 'G', 'T']
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        

        codons = []
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    codon = bases[i] + bases[j] + bases[k]
                    codons.append(codon)
        

        codon_to_nucleotide_matrix = torch.zeros(64, 3, 4, device=self.device)
        
        for codon_idx, codon in enumerate(codons):
            for pos, base in enumerate(codon):
                base_idx = base_to_idx[base]
                codon_to_nucleotide_matrix[codon_idx, pos, base_idx] = 1.0
        
        self.codon_to_nucleotide_matrix = codon_to_nucleotide_matrix
        self.codons = codons
        




    
    def detect_tensor_space(self, tensor: torch.Tensor) -> str:
        """
        Detect the tensor space type (nucleotide or codon).

        Args:
            tensor: Input tensor

        Returns:
            'nucleotide' (4-dimensional) or 'codon' (64-dimensional)
        """
        last_dim = tensor.shape[-1]
        
        if last_dim == 4:
            return 'nucleotide'
        elif last_dim == 64:
            return 'codon'
        else:
            raise ValueError(f"Cannot detect tensor space: last dimension is {last_dim}, expected 4 or 64")

    def codon_to_nucleotide(self,
                           codon_probs: torch.Tensor,
                           preserve_gradients: bool = True) -> torch.Tensor:
        """
        Convert codon probabilities to nucleotide probabilities.

        Args:
            codon_probs: Codon probability tensor [..., seq_len, 64]
            preserve_gradients: Whether to preserve gradient flow

        Returns:
            nucleotide_probs: Nucleotide probability tensor [..., seq_len*3, 4]
        """

        # Validate input dimensions
        if codon_probs.shape[-1] != 64:
            raise ValueError(f"Expected codon_probs last dimension to be 64, got {codon_probs.shape[-1]}")

        
        original_shape = codon_probs.shape
        batch_dims = original_shape[:-2]
        seq_len = original_shape[-2]
        
        # Reshape for batch processing
        codon_probs_flat = codon_probs.view(-1, seq_len, 64)
        batch_size = codon_probs_flat.shape[0]
        


        nucleotide_positions = torch.einsum('bsc,cij->bsij', codon_probs_flat, self.codon_to_nucleotide_matrix)
        

        nucleotide_probs = nucleotide_positions.view(batch_size, seq_len * 3, 4)
        

        final_shape = batch_dims + (seq_len * 3, 4)
        nucleotide_probs = nucleotide_probs.view(*final_shape)
        

        if preserve_gradients:
            nucleotide_probs = self._ensure_normalization(nucleotide_probs)
        
        return nucleotide_probs
    
    def nucleotide_to_codon_approximate(self,
                                      nucleotide_probs: torch.Tensor,
                                      amino_acid_sequence: str = None) -> torch.Tensor:
        """
        Convert nucleotide probabilities to codon probabilities (approximate).

        This is a more complex inverse transformation because multiple nucleotide
        combinations correspond to a single codon.

        Args:
            nucleotide_probs: Nucleotide probability tensor [..., seq_len*3, 4]
            amino_acid_sequence: Amino acid sequence to constrain valid codons

        Returns:
            codon_probs: Codon probability tensor [..., seq_len, 64]
        """

        if nucleotide_probs.shape[-1] != 4:
            raise ValueError(f"Expected nucleotide_probs last dimension to be 4, got {nucleotide_probs.shape[-1]}")

        original_shape = nucleotide_probs.shape
        batch_dims = original_shape[:-2]
        nucleotide_len = original_shape[-2]

        if nucleotide_len % 3 != 0:
            raise ValueError(f"Nucleotide length must be divisible by 3, got {nucleotide_len}")

        
        seq_len = nucleotide_len // 3
        
        # Reshape: [..., seq_len*3, 4] → [..., seq_len, 3, 4]
        nucleotide_reshaped = nucleotide_probs.view(*batch_dims, seq_len, 3, 4)
        batch_size = np.prod(batch_dims) if batch_dims else 1
        nucleotide_flat = nucleotide_reshaped.view(batch_size, seq_len, 3, 4)
        


        codon_probs = torch.zeros(batch_size, seq_len, 64, device=self.device)
        
        for codon_idx in range(64):

            codon_pattern = self.codon_to_nucleotide_matrix[codon_idx]  # [3, 4]
            

            codon_prob = torch.ones(batch_size, seq_len, device=self.device)
            
            for pos in range(3):

                position_prob = (nucleotide_flat[:, :, pos, :] * codon_pattern[pos].unsqueeze(0).unsqueeze(0)).sum(dim=-1)
                codon_prob = codon_prob * position_prob
            
            codon_probs[:, :, codon_idx] = codon_prob
        

        codon_probs = F.softmax(codon_probs, dim=-1)
        

        final_shape = batch_dims + (seq_len, 64)
        codon_probs = codon_probs.view(*final_shape)
        
        return codon_probs
    
    def _ensure_normalization(self, probs: torch.Tensor) -> torch.Tensor:
        """Ensure probability distribution is normalized."""
        prob_sums = probs.sum(dim=-1, keepdim=True)
        normalized_probs = probs / (prob_sums + 1e-8)
        return normalized_probs
    
    def auto_convert_to_nucleotide(self,
                                  tensor: torch.Tensor,
                                  preserve_gradients: bool = True) -> torch.Tensor:
        """
        Automatically convert tensor to nucleotide space if needed.

        Args:
            tensor: Input tensor (either codon or nucleotide space)
            preserve_gradients: Whether to preserve gradient flow

        Returns:
            nucleotide_tensor: Tensor in nucleotide space
        """
        space = self.detect_tensor_space(tensor)

        if space == 'nucleotide':
            # Already in nucleotide space, no conversion needed
            return tensor
        elif space == 'codon':
            # Convert from codon to nucleotide space
            return self.codon_to_nucleotide(tensor, preserve_gradients)
        else:
            raise ValueError(f"Unknown tensor space: {space}")

    def get_conversion_info(self, tensor: torch.Tensor) -> dict:
        """Get information about tensor space and conversion requirements."""
        space = self.detect_tensor_space(tensor)
        
        info = {
            'current_space': space,
            'tensor_shape': list(tensor.shape),
            'last_dim': tensor.shape[-1],
            'requires_conversion': space == 'codon',
            'device': tensor.device
        }

        if space == 'codon':
            # Calculate what the nucleotide shape would be after conversion
            original_shape = tensor.shape
            nucleotide_shape = original_shape[:-2] + (original_shape[-2] * 3, 4)
            info['converted_shape'] = list(nucleotide_shape)
        
        return info



_global_converter = None

def get_dimension_converter(device: torch.device = None) -> DimensionConverter:
    """Get global dimension converter instance."""
    global _global_converter
    
    if _global_converter is None:
        _global_converter = DimensionConverter(device)
    
    return _global_converter



def convert_to_nucleotide(tensor: torch.Tensor, device: torch.device = None) -> torch.Tensor:
    """Convenience function: convert tensor to nucleotide space."""
    converter = get_dimension_converter(device or tensor.device)
    return converter.auto_convert_to_nucleotide(tensor)

def detect_tensor_space(tensor: torch.Tensor) -> str:
    """Convenience function: detect tensor space type."""
    converter = get_dimension_converter(tensor.device)
    return converter.detect_tensor_space(tensor)