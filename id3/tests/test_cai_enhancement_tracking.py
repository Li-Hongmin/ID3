#!/usr/bin/env python
"""
CAI Enhancement Operator Tracking Test

"""

import torch
import logging
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root))

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons, codons
import numpy as np


logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_cai_enhancement_tracking():

    
    print("\n" + "="*80)
    print("CAI Enhancement Operator Tracking Test")
    print("="*80 + "\n")
    


    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    


    print(f"Amino acid sequence: {amino_acid_sequence}")
    

    enhancer = CAIEnhancementOperator(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_acid_sequence
    )


    


    seq_length = len(amino_acid_sequence)
    max_codons = 6
    

    valid_codon_mask = torch.zeros(seq_length, max_codons, dtype=torch.bool, device=device)
    codon_indices = torch.zeros(seq_length, max_codons, dtype=torch.long, device=device)
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            valid_codons = amino_acids_to_codons[aa]
            for i, codon in enumerate(valid_codons):
                if i < max_codons:
                    valid_codon_mask[pos, i] = True

                    codon_indices[pos, i] = codons.index(codon)
    
    print(f"✓ Valid codon mask created: {valid_codon_mask.shape}")
    print(f"✓ Codon indices created: {codon_indices.shape}")
    

    pi_accessibility = torch.rand(seq_length, max_codons, device=device)

    for pos in range(seq_length):
        valid_mask = valid_codon_mask[pos]
        if valid_mask.any():
            pi_accessibility[pos][~valid_mask] = 0
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    
    print(f"✓ Initial probability distribution created")
    


    initial_discrete = enhancer.discretize_distribution(pi_accessibility, valid_codon_mask)
    initial_cai = enhancer.compute_discrete_cai(
        initial_discrete, amino_acid_sequence, valid_codon_mask, 
        enhancer._build_standard_codon_indices(amino_acid_sequence, valid_codon_mask)
    )

    


    target_cai = 0.8

    

    original_binary_search = enhancer.binary_search_gamma
    
    def tracked_binary_search(*args, **kwargs):

        result = original_binary_search(*args, **kwargs)

        return result
    
    enhancer.binary_search_gamma = tracked_binary_search
    

    enhanced_dist, metadata = enhancer.apply_cai_enhancement(
        pi_accessibility,
        amino_acid_sequence,
        valid_codon_mask,
        codon_indices,
        target_cai
    )
    









    


    if metadata['final_cai'] < target_cai * 0.8:


        

        if enhancer.max_achievable_cai < target_cai:

        else:

    else:

    


    test_targets = [0.3, 0.5, 0.7, 0.8, 0.9]
    
    for target in test_targets:
        _, meta = enhancer.apply_cai_enhancement(
            pi_accessibility,
            amino_acid_sequence,
            valid_codon_mask,
            codon_indices,
            target
        )

    
    print("\n" + "="*80)

    print("="*80 + "\n")
    
    return metadata

def test_with_real_protein():
    """使用真实蛋白质测试"""
    print("\n" + "="*80)
    print("Real Protein CAI Enhancement Test")
    print("="*80 + "\n")
    

    amino_acid_sequence = "MIVPVIHFNSQSGTTVDIGLGPQRASLGGAIGSPGLRGRGRGGTGGRGR"
    print(f"Testing with O15263 (first 50 aa): {amino_acid_sequence}")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    enhancer = CAIEnhancementOperator(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_acid_sequence
    )
    
    print(f"Maximum achievable CAI: {enhancer.max_achievable_cai:.4f}")
    

    seq_length = len(amino_acid_sequence)
    max_codons = 6
    
    valid_codon_mask = torch.zeros(seq_length, max_codons, dtype=torch.bool, device=device)
    codon_indices = torch.zeros(seq_length, max_codons, dtype=torch.long, device=device)
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            valid_codons = amino_acids_to_codons[aa]
            for i, codon in enumerate(valid_codons):
                if i < max_codons:
                    valid_codon_mask[pos, i] = True
                    codon_indices[pos, i] = codons.index(codon)
    

    pi_accessibility = torch.rand(seq_length, max_codons, device=device)
    for pos in range(seq_length):
        valid_mask = valid_codon_mask[pos]
        if valid_mask.any():
            pi_accessibility[pos][~valid_mask] = 0
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    

    print("\nTesting different target CAI values:")
    test_targets = [0.3, 0.5, 0.7, 0.8]
    
    for target in test_targets:
        enhanced_dist, metadata = enhancer.apply_cai_enhancement(
            pi_accessibility,
            amino_acid_sequence,
            valid_codon_mask,
            codon_indices,
            target
        )
        
        print(f"\nTarget CAI: {target:.1f}")
        print(f"  Original CAI: {metadata['original_cai']:.4f}")
        print(f"  Final CAI: {metadata['final_cai']:.4f}")
        print(f"  Improvement: {metadata['final_cai'] - metadata['original_cai']:.4f}")
        print(f"  Optimal gamma: {metadata['optimal_gamma']:.4f}")
        print(f"  Constraint satisfied: {metadata['constraint_satisfied']}")

if __name__ == "__main__":

    test_cai_enhancement_tracking()
    

    test_with_real_protein()