#!/usr/bin/env python3
"""
Sequence constants definition module

Defines all sequence-related constants including codons, UTR sequences, etc.
Dynamically retrieves UTR templates through UTR loader to avoid hardcoding
"""

from typing import Dict, Set, Optional
from id3.config.utr_loader import get_utr_loader


# =============================================================================
# Codon constants
# =============================================================================

# Start codon
START_CODON = 'ATG'
START_CODON_RNA = 'AUG'

# Stop codons (DNA)
STOP_CODONS_DNA = {'TAA', 'TAG', 'TGA'}

# Stop codons (RNA)
STOP_CODONS_RNA = {'UAA', 'UAG', 'UGA'}

# Non-synonymous codons (amino acids encoded by only one codon)
NON_SYNONYMOUS_CODONS = {'ATG', 'TGG'}  # Met, Trp

# Standard genetic code table (DNA)
GENETIC_CODE_DNA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# Standard genetic code table (RNA)
GENETIC_CODE_RNA = {
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


# =============================================================================
# UTR sequence retrieval functions
# =============================================================================

def get_default_utr5() -> str:
    """Get default 5'UTR sequence

    Dynamically loaded from UTR template file to avoid hardcoding

    Returns:
        str: 5'UTR sequence
    """
    loader = get_utr_loader()
    return loader.get_default_utr5()


def get_default_utr3() -> str:
    """Get default 3'UTR sequence

    Dynamically loaded from UTR template file to avoid hardcoding

    Returns:
        str: 3'UTR sequence
    """
    loader = get_utr_loader()
    return loader.get_default_utr3()


def get_utr_sequences() -> Dict[str, str]:
    """Get all default UTR sequences

    Returns:
        Dict[str, str]: Dictionary containing utr5 and utr3
    """
    loader = get_utr_loader()
    return {
        'utr5': loader.get_default_utr5(),
        'utr3': loader.get_default_utr3()
    }


def get_utr_lengths() -> Dict[str, int]:
    """Get UTR sequence length information

    Returns:
        Dict[str, int]: UTR length information
    """
    loader = get_utr_loader()
    info = loader.get_utr_info()
    return {
        'utr5_length': info['utr5_length'],
        'utr3_length': info['utr3_length'],
        'total_utr_length': info['total_utr_length']
    }


# =============================================================================
# Sequence utility functions
# =============================================================================

def dna_to_rna(sequence: str) -> str:
    """Convert DNA sequence to RNA sequence

    Args:
        sequence: DNA sequence

    Returns:
        str: RNA sequence
    """
    return sequence.replace('T', 'U')


def rna_to_dna(sequence: str) -> str:
    """Convert RNA sequence to DNA sequence

    Args:
        sequence: RNA sequence

    Returns:
        str: DNA sequence
    """
    return sequence.replace('U', 'T')


def translate_dna(sequence: str) -> str:
    """Translate DNA sequence to amino acid sequence

    Args:
        sequence: DNA sequence

    Returns:
        str: Amino acid sequence
    """
    if len(sequence) % 3 != 0:
        raise ValueError("DNA sequence length must be a multiple of 3")

    amino_acids = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE_DNA.get(codon, 'X')  # X represents unknown amino acid
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)


def translate_rna(sequence: str) -> str:
    """Translate RNA sequence to amino acid sequence

    Args:
        sequence: RNA sequence

    Returns:
        str: Amino acid sequence
    """
    if len(sequence) % 3 != 0:
        raise ValueError("RNA sequence length must be a multiple of 3")

    amino_acids = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE_RNA.get(codon, 'X')  # X represents unknown amino acid
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)


def is_start_codon(codon: str, rna: bool = False) -> bool:
    """Check if codon is a start codon

    Args:
        codon: Codon sequence
        rna: Whether it is an RNA sequence

    Returns:
        bool: Whether it is a start codon
    """
    if rna:
        return codon.upper() == START_CODON_RNA
    else:
        return codon.upper() == START_CODON


def is_stop_codon(codon: str, rna: bool = False) -> bool:
    """Check if codon is a stop codon

    Args:
        codon: Codon sequence
        rna: Whether it is an RNA sequence

    Returns:
        bool: Whether it is a stop codon
    """
    codon = codon.upper()
    if rna:
        return codon in STOP_CODONS_RNA
    else:
        return codon in STOP_CODONS_DNA


def validate_nucleotide_sequence(sequence: str, rna: bool = False) -> bool:
    """Validate the validity of nucleotide sequence

    Args:
        sequence: Nucleotide sequence
        rna: Whether it is an RNA sequence

    Returns:
        bool: Whether the sequence is valid
    """
    sequence = sequence.upper()
    if rna:
        valid_nucleotides = set('AUCG')
    else:
        valid_nucleotides = set('ATCG')
    
    return set(sequence).issubset(valid_nucleotides)


# =============================================================================
# Configuration validation functions
# =============================================================================

def validate_sequence_constants() -> Dict[str, bool]:
    """Validate the validity of all sequence constants

    Returns:
        Dict[str, bool]: Validation results
    """
    results = {}

    try:
        # Validate UTR sequences
        utr5 = get_default_utr5()
        utr3 = get_default_utr3()

        results['utr5_valid'] = validate_nucleotide_sequence(utr5) and len(utr5) > 0
        results['utr3_valid'] = validate_nucleotide_sequence(utr3) and len(utr3) > 0

        # Validate codon constants
        results['start_codon_valid'] = len(START_CODON) == 3 and validate_nucleotide_sequence(START_CODON)
        results['stop_codons_valid'] = all(len(codon) == 3 and validate_nucleotide_sequence(codon)
                                         for codon in STOP_CODONS_DNA)

        # Validate genetic code table
        results['genetic_code_complete'] = len(GENETIC_CODE_DNA) == 64
        
        results['overall_valid'] = all(results.values())
        
    except Exception as e:
        results['error'] = str(e)
        results['overall_valid'] = False
    
    return results


# =============================================================================
# Compatibility functions (retained for backward compatibility)
# =============================================================================

def get_default_5utr() -> str:
    """Get default 5'UTR sequence (compatibility function)"""
    return get_default_utr5()


def get_default_3utr() -> str:
    """Get default 3'UTR sequence (compatibility function)"""
    return get_default_utr3()


if __name__ == "__main__":
    # Test code
    print("Sequence constants module test:")
    print(f"Start codon: {START_CODON}")
    print(f"Stop codons: {STOP_CODONS_DNA}")

    try:
        utr5 = get_default_utr5()
        utr3 = get_default_utr3()

        print(f"5'UTR length: {len(utr5)}")
        print(f"3'UTR length: {len(utr3)}")

        lengths = get_utr_lengths()
        print(f"UTR length information: {lengths}")

        validation = validate_sequence_constants()
        print(f"Sequence validation results: {validation}")

    except Exception as e:
        print(f"Test failed: {e}")