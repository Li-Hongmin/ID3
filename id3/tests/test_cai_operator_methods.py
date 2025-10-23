#!/usr/bin/env python3
"""



"""

import torch
import numpy as np
import time
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons

def test_cai_methods():

    


    

    print("="*60)
    

    torch.manual_seed(42)
    num_positions = len(amino_acid_sequence)
    max_codons = 6
    pi_accessibility = torch.rand(num_positions, max_codons)
    

    mask = torch.zeros_like(pi_accessibility, dtype=torch.bool)
    for i, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            n_codons = len(amino_acids_to_codons[aa])
            mask[i, :n_codons] = True

            if n_codons > 0:
                pi_accessibility[i, :n_codons] = pi_accessibility[i, :n_codons] / pi_accessibility[i, :n_codons].sum()
    

    pi_accessibility[~mask] = 0
    
    valid_codon_mask = mask
    codon_indices = torch.arange(max_codons).unsqueeze(0).expand(num_positions, -1)
    
    target_cai = 0.8
    
    results = {}
    


    print("-"*40)
    
    try:
        operator_bs = CAIEnhancementOperator(
            method='binary_search',
            species='ecoli_bl21de3',
            amino_acid_sequence=amino_acid_sequence
        )
        
        start_time = time.time()
        discrete_dist_bs, metadata_bs = operator_bs.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=amino_acid_sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=target_cai
        )
        bs_time = time.time() - start_time
        
        results['binary_search'] = {
            'time': bs_time,
            'metadata': metadata_bs,
            'success': True
        }
        




        
    except Exception as e:

        results['binary_search'] = {'success': False, 'error': str(e)}
    


    print("-"*40)
    
    try:
        operator_sado = CAIEnhancementOperator(
            method='sado',
            species='ecoli_bl21de3',
            amino_acid_sequence=amino_acid_sequence
        )
        
        start_time = time.time()
        discrete_dist_sado, metadata_sado = operator_sado.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=amino_acid_sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=target_cai
        )
        sado_time = time.time() - start_time
        
        results['sado'] = {
            'time': sado_time,
            'metadata': metadata_sado,
            'success': True
        }
        



        if 'unique' in metadata_sado:

        
    except Exception as e:

        results['sado'] = {'success': False, 'error': str(e)}
    


    print("-"*40)
    
    try:
        operator = CAIEnhancementOperator(
            method='binary_search',
            species='ecoli_bl21de3',
            amino_acid_sequence=amino_acid_sequence
        )
        

        

        operator.switch_method('sado')

        

        operator.switch_method('binary_search')

        

        
    except Exception as e:

    

    print("\n" + "="*60)

    print("="*60)
    
    if results.get('binary_search', {}).get('success') and results.get('sado', {}).get('success'):
        bs_time = results['binary_search']['time']
        sado_time = results['sado']['time']
        


        
        if sado_time > 0:
            speedup = bs_time / sado_time
            if speedup > 1:

            else:

    

    
    return results


if __name__ == "__main__":
    test_cai_methods()