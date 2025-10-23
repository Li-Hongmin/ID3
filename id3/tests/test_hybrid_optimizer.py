#!/usr/bin/env python3
"""








"""

import torch
import numpy as np
import time
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons

def test_all_methods():

    


    

    print("="*80)
    

    torch.manual_seed(42)
    num_positions = len(amino_acid_sequence)
    max_codons = 6
    

    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long)
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            n_codons = len(amino_acids_to_codons[aa])
            valid_codon_mask[pos, :n_codons] = True
            codon_indices[pos, :n_codons] = torch.arange(n_codons)
    
    target_cai = 0.8
    

    methods = ['binary_search', 'sado', 'hybrid_bs_sado']
    results = {}
    
    for method in methods:
        print(f"\n{'='*60}")

        print('='*60)
        
        try:
            operator = CAIEnhancementOperator(
                method=method,
                species='ecoli_bl21de3',
                amino_acid_sequence=amino_acid_sequence
            )
            

            sequence_hashes = []
            switch_count = 0
            
            for iteration in range(10):

                pi_accessibility = torch.rand(num_positions, max_codons)
                pi_accessibility = pi_accessibility * valid_codon_mask.float()
                

                for pos in range(num_positions):
                    if valid_codon_mask[pos].any():
                        pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
                

                if iteration > 3:

                    pi_accessibility = pi_accessibility ** 2
                    for pos in range(num_positions):
                        if valid_codon_mask[pos].any():
                            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
                
                start_time = time.time()
                discrete_dist, metadata = operator.apply_cai_enhancement(
                    pi_accessibility=pi_accessibility,
                    amino_acid_sequence=amino_acid_sequence,
                    valid_codon_mask=valid_codon_mask,
                    codon_indices=codon_indices,
                    target_cai=target_cai
                )
                elapsed_time = time.time() - start_time
                

                if discrete_dist.dim() > 1:
                    indices = discrete_dist.argmax(dim=-1)
                else:
                    indices = discrete_dist
                    
                seq_hash = hash(tuple(indices.cpu().numpy().tolist()))
                sequence_hashes.append(seq_hash)
                

                if method == 'hybrid_bs_sado' and metadata.get('switched_to_sado', False):
                    switch_count += 1

                
                if iteration == 0:




                    if 'optimal_gamma' in metadata:

                    if 'bs_gamma' in metadata:

            

            unique_sequences = len(set(sequence_hashes))
            repetition_rate = 1 - (unique_sequences / len(sequence_hashes))
            




            
            if method == 'hybrid_bs_sado':


                

                if hasattr(operator.optimizer, 'get_statistics'):
                    stats = operator.optimizer.get_statistics()


            
            results[method] = {
                'success': True,
                'unique_sequences': unique_sequences,
                'repetition_rate': repetition_rate,
                'switch_count': switch_count if method == 'hybrid_bs_sado' else 0
            }
            
        except Exception as e:

            results[method] = {'success': False, 'error': str(e)}
    

    print("\n" + "="*80)

    print("="*80)
    
    for method, result in results.items():
        if result['success']:



            if method == 'hybrid_bs_sado' and result['switch_count'] > 0:

        else:

    

    if all(results[m]['success'] for m in ['binary_search', 'hybrid_bs_sado']):
        bs_unique = results['binary_search']['unique_sequences']
        hybrid_unique = results['hybrid_bs_sado']['unique_sequences']
        
        if hybrid_unique > bs_unique:

        elif hybrid_unique == bs_unique:

        else:

    

    
    return results


def test_repetition_handling():
    """专门测试重复处理机制"""
    print("\n" + "="*80)
    print("重复处理机制测试")
    print("="*80)
    

    amino_acid_sequence = "MMM"
    
    print(f"测试序列: {amino_acid_sequence}")
    print("(所有M都只有一个密码子AUG，必然产生重复)")
    

    num_positions = len(amino_acid_sequence)
    max_codons = 6
    
    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long)
    
    for pos in range(num_positions):
        valid_codon_mask[pos, 0] = True
        codon_indices[pos, 0] = 0
    

    pi_accessibility = torch.zeros(num_positions, max_codons)
    pi_accessibility[:, 0] = 1.0
    

    operator = CAIEnhancementOperator(
        method='hybrid_bs_sado',
        species='ecoli_bl21de3',
        amino_acid_sequence=amino_acid_sequence
    )
    
    print("\n运行5次迭代（应该在第3次触发SADO）...")
    
    for i in range(5):
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=amino_acid_sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        
        switched = metadata.get('switched_to_sado', False)
        print(f"迭代 {i+1}: {'切换到SADO' if switched else '使用Binary Search'}")
        if switched:
            print(f"  切换原因: {metadata.get('switch_reason', 'unknown')}")
            print(f"  BS gamma传递: {metadata.get('bs_gamma', 'N/A')}")


if __name__ == "__main__":

    test_all_methods()
    

    test_repetition_handling()