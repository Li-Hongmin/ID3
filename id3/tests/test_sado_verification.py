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

def verify_optimization(seq_length: int):

    
    print(f"\n{'='*60}")

    print(f"{'='*60}")
    


    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(seq_length))


    


    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)


    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    


    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )


    


    start_time = time.time()
    
    distribution, metadata = optimizer.optimize(
        pi_accessibility=pi_accessibility,
        target_cai=0.8,
        use_binary_search=False,

        gamma=0.3
    )
    
    elapsed = time.time() - start_time
    








    




    




    selected_codons = torch.argmax(distribution, dim=1)


    
    return elapsed, metadata['final_cai'], metadata['constraint_satisfied']


def compare_methods(seq_length: int):
    """比较标准版和增强版"""
    
    print(f"\n{'='*60}")
    print(f"比较不同方法 ({seq_length} AA)")
    print(f"{'='*60}")
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(seq_length))
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    results = []
    

    print("\n标准SADO:")
    optimizer1 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    original = optimizer1._greedy_cai_optimization
    def limited(indices, target_cai, max_iterations=100, pi_accessibility=None):
        return original(indices, target_cai, min(max_iterations, 10), pi_accessibility)
    optimizer1._greedy_cai_optimization = limited
    
    start = time.time()
    dist1, meta1 = optimizer1.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=False
    )
    time1 = time.time() - start
    
    print(f"  时间: {time1:.4f}秒, CAI: {meta1['final_cai']:.4f}")
    results.append(('标准', time1, meta1['final_cai']))
    

    print("\n增强版SADO:")
    optimizer2 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    start = time.time()
    dist2, meta2 = optimizer2.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True
    )
    time2 = time.time() - start
    
    print(f"  时间: {time2:.4f}秒, CAI: {meta2['final_cai']:.4f}")
    results.append(('增强', time2, meta2['final_cai']))
    

    print("\n比较结果:")
    print(f"  速度提升: {time1/time2:.2f}x")
    print(f"  CAI差异: {meta2['final_cai'] - meta1['final_cai']:+.4f}")
    
    return results


def main():

    
    print("\n" + "="*80)

    print("="*80)
    

    time_300, cai_300, satisfied_300 = verify_optimization(300)
    

    time_3000, cai_3000, satisfied_3000 = verify_optimization(3000)
    

    compare_methods(200)
    

    print("\n" + "="*80)

    print("="*80)
    



    



    



if __name__ == '__main__':
    main()