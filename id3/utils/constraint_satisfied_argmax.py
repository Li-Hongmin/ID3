#!/usr/bin/env python3
"""


"""

import torch
import torch.nn.functional as F
from .constants import amino_acid_token_map, amino_acid_to_codon_matrix, codon_to_rna_matrix


def get_constraint_satisfied_argmax(rna_probs: torch.Tensor,
                                  amino_acid_sequence: str) -> torch.Tensor:
    """




    Args:



    Returns:

    """
    batch_size, seq_len, num_bases = rna_probs.shape
    num_positions = len(amino_acid_sequence)
    device = rna_probs.device
    dtype = rna_probs.dtype


    assert seq_len == num_positions * 3, f"RNA长度{seq_len} != 氨基酸长度{num_positions} * 3"


    rna_reshaped = rna_probs.view(batch_size, num_positions, 3, 4)


    amino_acid_indices = torch.tensor([amino_acid_token_map[aa] for aa in amino_acid_sequence], device=device)


    # amino_acid_to_codon_matrix: [20, 64]
    # amino_acid_indices: [num_positions]
    aa_to_codon = amino_acid_to_codon_matrix.to(device)
    valid_codon_matrix = aa_to_codon[amino_acid_indices]  # [num_positions, 64]



    max_codons = valid_codon_matrix.sum(dim=-1).max().int().item()



    valid_positions = valid_codon_matrix.nonzero(as_tuple=False)  # [total_valid, 2]
    pos_indices = valid_positions[:, 0]
    codon_indices_flat = valid_positions[:, 1]


    counts = torch.bincount(pos_indices, minlength=num_positions)


    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=device)
    mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=device)


    start_idx = 0
    for pos in range(num_positions):
        n_valid = counts[pos].item()
        if n_valid > 0:
            codon_indices[pos, :n_valid] = codon_indices_flat[start_idx:start_idx + n_valid]
            mask[pos, :n_valid] = True
            start_idx += n_valid


    codon_to_rna_device = codon_to_rna_matrix.to(device)
    valid_codons_padded = codon_to_rna_device[codon_indices].to(dtype)
    valid_codons_padded[~mask] = 0


    eps = torch.finfo(dtype).eps
    log_probs = torch.log(rna_reshaped + eps)


    codon_scores = torch.einsum('bpjk,pmjk->bpm', log_probs, valid_codons_padded)


    mask_expanded = mask.unsqueeze(0).expand(batch_size, -1, -1)
    codon_scores[~mask_expanded] = -float('inf')


    best_indices = codon_scores.argmax(dim=-1)  # [batch_size, num_positions]


    valid_codons_expanded = valid_codons_padded.unsqueeze(0).expand(batch_size, -1, -1, -1, -1)
    indices_expanded = best_indices.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, 3, 4).unsqueeze(2)


    result = torch.gather(valid_codons_expanded, 2, indices_expanded).squeeze(2)

    return result.view(batch_size, seq_len, num_bases)


if __name__ == "__main__":

    amino_acid_sequence = "MGK"
    batch_size = 2
    seq_len = len(amino_acid_sequence) * 3
    
    torch.manual_seed(42)
    rna_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)
    
    result = get_constraint_satisfied_argmax(rna_probs, amino_acid_sequence)
    print(f"测试通过: {result.shape}")