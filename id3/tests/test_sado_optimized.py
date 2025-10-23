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


def test_optimized_sado():

    
    print("\n" + "="*70)

    print("="*70 + "\n")
    

    test_cases = [




    ]
    
    for length, desc in test_cases:
        print(f"\n{desc} ({length} AA)")
        print("-"*50)
        

        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        amino_sequence = ''.join(random.choice(amino_acids) for _ in range(length))
        

        torch.manual_seed(42)
        pi_accessibility = torch.rand(length, 6)
        pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
        
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        optimizer1 = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        
        start = time.time()
        dist1, meta1 = optimizer1.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False,
            gamma=0.3
        )
        time_standard = time.time() - start
        

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

            gamma=0.3
        )
        time_enhanced = time.time() - start
        



        
        if time_standard > 0:
            speedup = time_standard / time_enhanced
            if speedup > 1:

            else:



def test_300aa_detailed():
    """详细测试300 AA序列"""
    
    print("\n" + "="*70)
    print("300 AA序列详细测试")
    print("="*70 + "\n")
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(300))
    

    torch.manual_seed(123)
    pi_accessibility = torch.rand(300, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    num_runs = 3
    times_standard = []
    times_enhanced = []
    cai_standard = []
    cai_enhanced = []
    
    print(f"运行{num_runs}次测试...")
    
    for i in range(num_runs):

        optimizer1 = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        
        start = time.time()
        dist1, meta1 = optimizer1.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False
        )
        times_standard.append(time.time() - start)
        cai_standard.append(meta1['final_cai'])
        

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
        times_enhanced.append(time.time() - start)
        cai_enhanced.append(meta2['final_cai'])
        
        print(f"  运行{i+1}: 标准={times_standard[-1]:.3f}s, 增强={times_enhanced[-1]:.3f}s")
    

    print("\n统计结果:")
    print(f"  标准SADO:")
    print(f"    平均时间: {np.mean(times_standard):.3f} ± {np.std(times_standard):.3f}s")
    print(f"    平均CAI:  {np.mean(cai_standard):.4f} ± {np.std(cai_standard):.4f}")
    
    print(f"  增强版(优化):")
    print(f"    平均时间: {np.mean(times_enhanced):.3f} ± {np.std(times_enhanced):.3f}s")
    print(f"    平均CAI:  {np.mean(cai_enhanced):.4f} ± {np.std(cai_enhanced):.4f}")
    
    avg_speedup = np.mean(times_standard) / np.mean(times_enhanced)
    if avg_speedup > 1:
        print(f"\n⚡ 优化后增强版平均速度提升: {avg_speedup:.2f}x")
    else:
        print(f"\n标准版平均更快: {1/avg_speedup:.2f}x")


def main():

    
    try:

        test_optimized_sado()
        

        test_300aa_detailed()
        
        print("\n" + "="*70)

        print("="*70)
        





        
    except Exception as e:

        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()