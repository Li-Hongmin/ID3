"""



"""

import torch
import numpy as np
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import json
import logging

logger = logging.getLogger(__name__)


def load_cai_weights(species: str = 'ecoli_bl21de3') -> Tuple[Dict[str, float], torch.Tensor]:
    """

    
    Args:

        
    Returns:

    """
    weights_file = Path(__file__).parent.parent.parent.parent / 'data' / 'codon_references' / f'{species}_wi_weights_comparison.json'
    
    if not weights_file.exists():
        raise FileNotFoundError(f"CAI weights file not found: {weights_file}")
    
    with open(weights_file, 'r') as f:
        data = json.load(f)
    

    wi_table = data['wi_table']
    

    rna_wi_table = {}
    for dna_codon, weight in wi_table.items():
        rna_codon = dna_codon.replace('T', 'U')
        rna_wi_table[rna_codon] = weight
    

    weights_tensor = create_standard_weights_tensor(rna_wi_table)
    
    return rna_wi_table, weights_tensor


def create_standard_weights_tensor(rna_wi_table: Dict[str, float]) -> torch.Tensor:
    """

    
    Args:

        
    Returns:

    """
    bases = ['U', 'C', 'A', 'G']
    weights = []
    
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon = base1 + base2 + base3
                weight = rna_wi_table.get(codon, 0.0)
                weights.append(weight)
    
    return torch.tensor(weights, dtype=torch.float32)


def compute_cai(sequence_indices: torch.Tensor, 
                weights_tensor: torch.Tensor,
                amino_acid_sequence: str) -> float:
    """

    
    Args:



        
    Returns:

    """

    if sequence_indices.dim() == 1:
        weights = weights_tensor[sequence_indices]
    else:
        weights = torch.gather(weights_tensor.unsqueeze(0).expand_as(sequence_indices), 
                              1, sequence_indices)
    

    log_weights = torch.log(weights + 1e-10)
    cai = torch.exp(log_weights.mean()).item()
    
    return cai


def discretize_distribution(distribution: torch.Tensor, 
                           valid_mask: torch.Tensor) -> torch.Tensor:
    """

    
    Args:


        
    Returns:

    """

    masked_dist = distribution * valid_mask
    

    masked_dist = masked_dist / (masked_dist.sum(dim=-1, keepdim=True) + 1e-10)
    

    max_indices = masked_dist.argmax(dim=-1)
    

    discrete = torch.zeros_like(distribution)
    discrete.scatter_(1, max_indices.unsqueeze(-1), 1.0)
    

    discrete = discrete * valid_mask
    
    return discrete


def compute_distribution_similarity(dist1: torch.Tensor, 
                                   dist2: torch.Tensor) -> float:
    """

    
    Args:


        
    Returns:

    """

    dist1 = dist1 + 1e-10
    dist2 = dist2 + 1e-10
    

    dist1 = dist1 / dist1.sum(dim=-1, keepdim=True)
    dist2 = dist2 / dist2.sum(dim=-1, keepdim=True)
    

    kl_div = (dist1 * (dist1.log() - dist2.log())).sum()
    
    return kl_div.item()


def interpolate_distributions(dist1: torch.Tensor,
                             dist2: torch.Tensor,
                             gamma: float) -> torch.Tensor:
    """

    
    Args:



        
    Returns:

    """
    return (1 - gamma) * dist1 + gamma * dist2