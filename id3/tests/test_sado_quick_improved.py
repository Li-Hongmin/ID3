#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import hashlib
import time

from id3.optimizers.cai.sado import SADOOptimizer


def test_quick():

    print("\n" + "="*60)

    print("="*60)
    

    seq_length = 300
    num_iterations = 10
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    hashes = []
    cai_values = []
    prob_scores = []
    
    for i in range(num_iterations):
        start = time.time()
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True,
            gamma=0.3
        )
        elapsed = time.time() - start
        

        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        hashes.append(seq_hash)
        

        cai_values.append(meta['final_cai'])
        

        log_prob = 0.0
        for pos in range(seq_length):
            codon_idx = selected[pos].item()
            if pos < len(optimizer.codon_choices) and optimizer.codon_choices[pos]:
                for choice in optimizer.codon_choices[pos]:
                    if choice['local_index'] == codon_idx:
                        orig_idx = choice['original_local_index']
                        if orig_idx < pi_accessibility.shape[1]:
                            prob = pi_accessibility[pos, orig_idx].item()
                            if prob > 0:
                                log_prob += np.log(prob)
                        break
        prob_score = np.exp(log_prob / seq_length)
        prob_scores.append(prob_score)
        

    

    unique_count = len(set(hashes))
    avg_prob = np.mean(prob_scores)
    





    

    hash_counts = {}
    for h in hashes:
        hash_counts[h] = hash_counts.get(h, 0) + 1
    
    max_repeat = max(hash_counts.values())

    


    if unique_count >= 8:

    elif unique_count >= 5:

    else:

    
    if avg_prob > 0.15:

    else:



if __name__ == '__main__':
    test_quick()