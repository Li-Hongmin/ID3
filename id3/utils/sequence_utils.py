"""
Utility functions for sequence conversion and validation.
"""

# Copyright (c) 2025 University of Tokyo
# Licensed under CC BY-NC-SA 4.0
# For commercial use, contact: lihongmin@edu.k.u-tokyo.ac.jp

import torch
from typing import Union

# Genetic code
GENETIC_CODE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def rna_string_to_amino_acids(rna_sequence: str) -> str:
    """
    Convert RNA sequence string to amino acid sequence.
    
    Args:
        rna_sequence: RNA sequence string (e.g., "AUGAAACUG")
        
    Returns:
        Amino acid sequence string
    """
    # Make sure it's uppercase and replace T with U if needed
    rna_sequence = rna_sequence.upper().replace('T', 'U')
    
    amino_acids = []
    for i in range(0, len(rna_sequence), 3):
        if i + 3 <= len(rna_sequence):
            codon = rna_sequence[i:i+3]
            if codon in GENETIC_CODE:
                amino_acids.append(GENETIC_CODE[codon])
            else:
                # Handle unknown codons
                amino_acids.append('X')
    
    return ''.join(amino_acids)


def rna_to_amino_acids(rna_input: Union[str, torch.Tensor]) -> str:
    """
    Convert RNA input (string or tensor) to amino acid sequence.
    
    Args:
        rna_input: RNA sequence as string or tensor
        
    Returns:
        Amino acid sequence string
    """
    if isinstance(rna_input, str):
        return rna_string_to_amino_acids(rna_input)
    
    # Handle tensor input
    if isinstance(rna_input, torch.Tensor):
        # Convert tensor to string first
        nucleotides = ['A', 'C', 'G', 'U']
        
        if rna_input.dim() == 3:
            # Take first batch
            rna_input = rna_input[0]
        
        if rna_input.dim() == 2:
            # Convert probabilities to indices
            indices = torch.argmax(rna_input, dim=-1)
            rna_string = ''.join([nucleotides[idx] for idx in indices.cpu().numpy()])
        else:
            # Already indices
            rna_string = ''.join([nucleotides[idx] for idx in rna_input.cpu().numpy()])
        
        return rna_string_to_amino_acids(rna_string)
    
    raise ValueError(f"Unsupported input type: {type(rna_input)}")


def validate_amino_acid_constraint(rna_sequence: Union[str, torch.Tensor],
                                  target_amino_acids: str) -> bool:
    """
    Check if RNA sequence translates to target amino acids.

    Args:
        rna_sequence: RNA sequence (string or tensor)
        target_amino_acids: Target amino acid sequence

    Returns:
        True if translation matches target
    """
    translated = rna_to_amino_acids(rna_sequence)
    return translated == target_amino_acids


def sequence_to_one_hot(sequence: str, device: str = 'cpu') -> torch.Tensor:
    """
    Convert RNA sequence string to one-hot encoding tensor.

    Args:
        sequence: RNA sequence string (e.g., "AUGC")
        device: Device for tensor (default: 'cpu')

    Returns:
        One-hot encoded tensor [seq_len, 4] where 4 corresponds to [A, C, G, U]
    """
    # Nucleotide to index mapping
    nucleotide_to_index = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}  # T is treated as U

    # Convert sequence to uppercase
    sequence = sequence.upper()

    # Create one-hot tensor
    seq_len = len(sequence)
    one_hot = torch.zeros(seq_len, 4, device=device)

    for i, nucleotide in enumerate(sequence):
        if nucleotide in nucleotide_to_index:
            idx = nucleotide_to_index[nucleotide]
            one_hot[i, idx] = 1.0
        else:
            # Handle unknown nucleotides (distribute equally)
            one_hot[i, :] = 0.25

    return one_hot