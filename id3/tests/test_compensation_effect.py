#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import hashlib

from id3.optimizers.cai.sado_incremental import SADOIncrementalOptimizer


def test_compensation_effectiveness():

    print("\n" + "="*70)

    print("="*70)
    

    seq_length = 300
    num_iterations = 20
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    




    

    optimizer = SADOIncrementalOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    cai_before_comp = []
    cai_after_comp = []
    compensation_positions = []
    unique_sequences = set()
    

    print("-" * 50)
    
    for i in range(num_iterations):

        distribution, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8
        )
        

        selected = torch.argmax(distribution, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        unique_sequences.add(seq_hash)
        

        before_cai = metadata.get('cai_before_compensation', metadata['final_cai'])
        after_cai = metadata['final_cai']
        comp_pos = metadata.get('compensation_positions', [])
        
        cai_before_comp.append(before_cai)
        cai_after_comp.append(after_cai)
        compensation_positions.append(comp_pos)
        
        if (i + 1) <= 5 or (i + 1) % 5 == 0:
            improvement = after_cai - before_cai





    

    print("\n" + "="*70)

    print("="*70)
    

    avg_before = np.mean(cai_before_comp)
    avg_after = np.mean(cai_after_comp)
    avg_improvement = avg_after - avg_before
    






    

    all_positions = []
    for pos_list in compensation_positions:
        all_positions.extend(pos_list)
    
    if all_positions:
        from collections import Counter
        pos_counter = Counter(all_positions)




        

        most_common = pos_counter.most_common(10)


    




    


    if avg_after >= 0.75 and len(unique_sequences)/num_iterations >= 0.9:

    elif avg_after >= 0.70 and len(unique_sequences)/num_iterations >= 0.8:

    else:

    


    if avg_after < 0.75:


    if len(unique_sequences)/num_iterations < 0.9:


    if len(pos_counter) < seq_length * 0.3:



if __name__ == '__main__':
    test_compensation_effectiveness()