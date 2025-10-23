#!/usr/bin/env python3
"""

"""

import torch
import numpy as np
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.constraints.amino_matching import AminoMatchingSoftmax

def test_cai_enhancement():

    

    test_sequences = {



    }
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    

    cai_operator = CAIEnhancementOperator(species='ecoli_bl21de3', device=device)
    

    print("="*80)
    
    for name, protein_seq in test_sequences.items():

        print("-"*60)
        

        constraint = AminoMatchingSoftmax(protein_seq, verbose=False)
        


        result = constraint.forward(alpha=1.0, beta=0.0, tau=1.0)
        


        if 'codon_probs' in result:
            pi_accessibility = result['codon_probs']
        else:

            for key in ['prob', 'probs', 'pi', 'codon_distribution']:
                if key in result:
                    pi_accessibility = result[key]
                    break
            else:


        

        valid_codon_mask = constraint.valid_codon_mask
        

        if pi_accessibility.dim() == 3:

            pi_accessibility = pi_accessibility.squeeze(0)
        


        

        initial_discrete = cai_operator.discretize_distribution(pi_accessibility, valid_codon_mask)
        initial_cai = compute_cai_from_discrete(initial_discrete, protein_seq, cai_operator)

        

        target_cais = [0.3, 0.5, 0.7, 0.8]
        
        for target_cai in target_cais:
            try:

                enhanced_dist, metadata = cai_operator.apply_cai_enhancement(
                    pi_accessibility=pi_accessibility,
                    amino_acid_sequence=protein_seq,
                    valid_codon_mask=valid_codon_mask,
                    target_cai=target_cai
                )
                

                actual_cai = metadata['final_cai']
                gamma = metadata['optimal_gamma']
                


                

                if target_cai == 0.8:
                    analyze_enhancement_details(
                        cai_operator, protein_seq, pi_accessibility, 
                        enhanced_dist, valid_codon_mask, metadata
                    )
                    
            except Exception as e:

        

        if cai_operator.max_achievable_cai is not None:


def compute_cai_from_discrete(discrete_dist, protein_seq, cai_operator):
    """从离散分布计算CAI"""
    

    codon_indices = build_codon_indices(protein_seq, discrete_dist.shape[1])
    valid_codon_mask = build_valid_mask(protein_seq, discrete_dist.shape[1])
    

    cai = cai_operator.compute_discrete_cai(
        discrete_dist, protein_seq, valid_codon_mask, codon_indices
    )
    
    return cai

def build_codon_indices(protein_seq, max_codons):

    genetic_code = {
        'F': ['UUU', 'UUC'],
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
        'Y': ['UAU', 'UAC'],
        'C': ['UGU', 'UGC'],
        'W': ['UGG'],
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'H': ['CAU', 'CAC'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'I': ['AUU', 'AUC', 'AUA'],
        'M': ['AUG'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'N': ['AAU', 'AAC'],
        'K': ['AAA', 'AAG'],
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'D': ['GAU', 'GAC'],
        'E': ['GAA', 'GAG'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG']
    }
    
    seq_length = len(protein_seq)
    codon_indices = torch.zeros(seq_length, max_codons, dtype=torch.long)
    
    for pos, aa in enumerate(protein_seq):
        if aa in genetic_code:
            codons = genetic_code[aa]
            for i, codon in enumerate(codons):
                if i < max_codons:

                    codon_idx = codon_to_standard_index(codon)
                    codon_indices[pos, i] = codon_idx
    
    return codon_indices

def build_valid_mask(protein_seq, max_codons):
    """构建有效密码子掩码"""
    genetic_code = {
        'F': 2, 'L': 6, 'S': 6, 'Y': 2, 'C': 2, 'W': 1,
        'P': 4, 'H': 2, 'Q': 2, 'R': 6, 'I': 3, 'M': 1,
        'T': 4, 'N': 2, 'K': 2, 'V': 4, 'A': 4, 'D': 2,
        'E': 2, 'G': 4
    }
    
    seq_length = len(protein_seq)
    valid_mask = torch.zeros(seq_length, max_codons, dtype=torch.bool)
    
    for pos, aa in enumerate(protein_seq):
        if aa in genetic_code:
            num_codons = genetic_code[aa]
            valid_mask[pos, :num_codons] = True
    
    return valid_mask

def codon_to_standard_index(codon):

    nucleotide_map = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
    
    if len(codon) != 3:
        return 0
    
    try:
        index = (nucleotide_map[codon[0]] * 16 + 
                nucleotide_map[codon[1]] * 4 + 
                nucleotide_map[codon[2]])
        return index
    except KeyError:
        return 0

def analyze_enhancement_details(cai_operator, protein_seq, pi_accessibility, 
                               enhanced_dist, valid_codon_mask, metadata):
    """详细分析CAI增强的效果"""
    print("\n  详细分析(目标CAI=0.8):")
    

    optimal_count = 0
    total_positions = len(protein_seq)
    
    for pos in range(total_positions):
        selected_slot = torch.argmax(enhanced_dist[pos]).item()
        if valid_codon_mask[pos, selected_slot]:


            if selected_slot == 0:
                optimal_count += 1
    
    print(f"    选择最优密码子的位置: {optimal_count}/{total_positions} "
          f"({optimal_count/total_positions*100:.1f}%)")
    print(f"    插值参数gamma: {metadata['optimal_gamma']:.4f}")
    print(f"    CAI提升: {metadata['final_cai'] - metadata['original_cai']:.4f}")

if __name__ == "__main__":
    print("CAI增强操作器测试")
    print("="*80)
    test_cai_enhancement()