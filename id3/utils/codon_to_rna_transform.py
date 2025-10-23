#!/usr/bin/env python3
"""




"""

import torch
from typing import Optional


def reconstruct_rna_from_codon_probs(
    codon_probs: torch.Tensor,
    codon_encodings: torch.Tensor,
    batch_first: bool = True
) -> torch.Tensor:
    """

    


    
    Args:

            - 2D: [num_positions, max_codons]
            - 3D: [batch_size, num_positions, max_codons]


        
    Returns:




    """

    if codon_probs.dim() == 3:
        batch_size = codon_probs.shape[0]
        num_positions = codon_probs.shape[1]
    else:  # dim == 2
        batch_size = 1
        num_positions = codon_probs.shape[0]
        codon_probs = codon_probs.unsqueeze(0)
    

    # codon_probs: [batch, positions, max_codons] -> [batch, positions, max_codons, 1, 1]
    probs_expanded = codon_probs.unsqueeze(-1).unsqueeze(-1)
    

    # codon_encodings: [positions, max_codons, 3, 4] -> [batch, positions, max_codons, 3, 4]
    encodings_expanded = codon_encodings.unsqueeze(0).expand(batch_size, -1, -1, -1, -1)
    

    weighted_encodings = probs_expanded * encodings_expanded
    rna_per_position = weighted_encodings.sum(dim=2)
    

    rna_sequence = rna_per_position.view(batch_size, -1, 4)
    

    if batch_size == 1 and not batch_first:
        rna_sequence = rna_sequence.squeeze(0)
    
    return rna_sequence


def pi_codon_to_amino(
    codon_probs: torch.Tensor,
    codon_encodings: torch.Tensor
) -> torch.Tensor:
    """

    

    
    Args:


        
    Returns:

    """
    return reconstruct_rna_from_codon_probs(
        codon_probs, 
        codon_encodings, 
        batch_first=True
    )