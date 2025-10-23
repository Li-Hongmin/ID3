#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import time
import random

from id3.optimizers.cai.sado import SADOOptimizer


def test_length(length: int, max_iter: int = 30):

    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    results = {}
    

    optimizer1 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    original_greedy = optimizer1._greedy_cai_optimization
    def limited_greedy(indices, target_cai, max_iterations=100, pi_accessibility=None):
        return original_greedy(indices, target_cai, min(max_iterations, max_iter), pi_accessibility)
    optimizer1._greedy_cai_optimization = limited_greedy
    
    start = time.time()
    dist1, meta1 = optimizer1.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=False
    )
    time1 = time.time() - start
    

    

    optimizer2 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    original_diff = optimizer2._difference_driven_greedy_optimization
    def limited_diff(initial_indices, pi_accessibility, target_cai, max_iterations=50):
        return original_diff(initial_indices, pi_accessibility, target_cai, min(max_iterations, max_iter))
    optimizer2._difference_driven_greedy_optimization = limited_diff
    
    start = time.time()
    dist2, meta2 = optimizer2.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True
    )
    time2 = time.time() - start
    

    
    return results


def main():
    print("\n" + "="*70)

    print("="*70 + "\n")
    
    lengths = [100, 300, 500]
    

    print("-"*70)
    
    for length in lengths:
        results = test_length(length, max_iter=30)
        


        
        ratio = std['time'] / enh['time'] if enh['time'] > 0 else 0
        
        print(f"{length:<10} "
              f"CAI={std['cai']:.3f}, {std['time']:.2f}s   "
              f"CAI={enh['cai']:.3f}, {enh['time']:.2f}s   "
              f"{ratio:.2f}x")
    
    print("\n" + "="*70)

    print("="*70)





if __name__ == '__main__':
    main()