#!/usr/bin/env python3
"""


"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import time

from id3.optimizers.cai.sado import SADOOptimizer


def test_basic_functionality():

    
    print("\n" + "="*60)

    print("="*60 + "\n")
    

    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLK"
    seq_len = len(amino_sequence)

    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_len, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    
    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    print("\n" + "-"*40)

    print("-"*40)
    
    optimizer.reset()
    start_time = time.time()
    dist1, meta1 = optimizer.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=False,
        gamma=0.3
    )
    time1 = time.time() - start_time
    




    
    print("\n" + "-"*40)

    print("-"*40)
    
    optimizer.reset()
    start_time = time.time()
    dist2, meta2 = optimizer.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True,
        gamma=0.3
    )
    time2 = time.time() - start_time
    




    
    print("\n" + "-"*40)

    print("-"*40)
    
    optimizer.reset()
    start_time = time.time()
    dist3, meta3 = optimizer.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=True,
        use_difference_driven=True,
        gamma=0.3
    )
    time3 = time.time() - start_time
    




    if 'binary_search_alpha' in meta3:

    
    print("\n" + "="*40)

    print("="*40)
    



    
    if time1 > 0:



    



def test_iteration_stability():
    """测试多次迭代的稳定性"""
    
    print("\n" + "="*60)
    print("测试迭代稳定性")
    print("="*60 + "\n")
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVS"
    seq_len = len(amino_sequence)
    

    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_len, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        amino_acid_sequence=amino_sequence
    )
    
    cai_values = []
    
    print("运行5次迭代...")
    for i in range(5):
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=True,
            use_difference_driven=True
        )
        cai_values.append(meta['final_cai'])
        print(f"  迭代 {i+1}: CAI={meta['final_cai']:.4f}, 初始化={meta.get('init_method', 'unknown')}")
    
    avg_cai = np.mean(cai_values)
    std_cai = np.std(cai_values)
    
    print(f"\n平均CAI: {avg_cai:.4f} ± {std_cai:.4f}")
    print(f"序列多样性: {optimizer.get_diversity_stats()['unique_ratio']*100:.1f}%")
    
    print("\n✅ 稳定性测试完成!")


def main():

    
    try:

        test_basic_functionality()
        

        test_iteration_stability()
        
        print("\n" + "="*60)

        print("="*60 + "\n")
        
    except Exception as e:

        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()