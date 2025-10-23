"""


"""

import json
from pathlib import Path
from typing import Dict, Any
import numpy as np


def verify_amino_acid_constraint(dna_sequence: str, 
                                amino_sequence: str) -> bool:
    """

    
    Args:


        
    Returns:

    """
    codon_table = get_codon_table()
    translated = translate_dna(dna_sequence, codon_table)
    return translated == amino_sequence


def translate_dna(dna_sequence: str, 
                  codon_table: Dict[str, str] = None) -> str:
    """

    
    Args:


        
    Returns:

    """
    if codon_table is None:
        codon_table = get_codon_table()
    
    protein = []
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        if codon in codon_table:
            aa = codon_table[codon]
            if aa == '*':
                break
            protein.append(aa)
    
    return ''.join(protein)


def get_codon_table() -> Dict[str, str]:
    """

    
    Returns:

    """
    return {
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


def check_sequence_quality(dna_sequence: str) -> Dict[str, Any]:
    """

    
    Args:

        
    Returns:

    """
    gc_content = calculate_gc_content(dna_sequence)
    has_start = dna_sequence.startswith('ATG')
    stop_codons = ['TAA', 'TAG', 'TGA']
    has_stop = any(dna_sequence[-3:] == stop for stop in stop_codons)
    
    return {
        'length': len(dna_sequence),
        'gc_content': gc_content,
        'has_start_codon': has_start,
        'has_stop_codon': has_stop,
        'is_multiple_of_3': len(dna_sequence) % 3 == 0
    }


def calculate_gc_content(dna_sequence: str) -> float:
    """

    
    Args:

        
    Returns:

    """
    if not dna_sequence:
        return 0.0
    
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    return (gc_count / len(dna_sequence)) * 100


def compute_cai_from_sequence(rna_sequence: str) -> float:
    """

    
    Args:

        
    Returns:

    """

    weights_file = Path(__file__).parent.parent.parent / 'data' / 'codon_references' / 'ecoli_bl21de3_wi_weights_comparison.json'
    
    if not weights_file.exists():
        return 0.0
        
    with open(weights_file, 'r') as f:
        data = json.load(f)
        wi_table = data['wi_table']
    

    codons = []
    for i in range(0, len(rna_sequence) - 2, 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            codons.append(codon)
    
    if not codons:
        return 0.0
    

    cai_product = 1.0
    valid_codons = 0
    
    for codon in codons:

        dna_codon = codon.replace('U', 'T')
        
        if dna_codon in wi_table:
            weight = wi_table[dna_codon]
            if weight > 0:
                cai_product *= weight
                valid_codons += 1
    
    if valid_codons == 0:
        return 0.0
    

    cai = cai_product ** (1.0 / valid_codons)
    
    return cai