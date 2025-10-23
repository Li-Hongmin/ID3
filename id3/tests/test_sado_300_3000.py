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


def test_specific_length(length: int):

    
    print(f"\n{'='*60}")

    print(f"{'='*60}\n")
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    


    
    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    original_greedy = optimizer._greedy_cai_optimization
    def limited_greedy(indices, target_cai, max_iterations=20, pi_accessibility=None):
        return original_greedy(indices, target_cai, min(max_iterations, 20), pi_accessibility)
    optimizer._greedy_cai_optimization = limited_greedy
    
    start_time = time.time()
    try:
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False,
            gamma=0.3
        )
        elapsed = time.time() - start_time
        

        print(f"  CAI: {meta['final_cai']:.4f}")


        
        return elapsed, meta['final_cai']
        
    except Exception as e:
        elapsed = time.time() - start_time


        return elapsed, 0


def main():
    """主函数"""
    
    print("\n" + "="*80)
    print("SADO性能快速测试 - 300 vs 3000 AA")
    print("="*80)
    

    time_300, cai_300 = test_specific_length(300)
    

    time_3000, cai_3000 = test_specific_length(3000)
    

    print("\n" + "="*60)
    print("总结")
    print("="*60)
    
    print(f"\n300 AA:  时间={time_300:.3f}秒, CAI={cai_300:.4f}")
    print(f"3000 AA: 时间={time_3000:.3f}秒, CAI={cai_3000:.4f}")
    
    if time_3000 > 0 and time_300 > 0:
        print(f"\n时间比: 3000 AA是300 AA的 {time_3000/time_300:.1f}倍")
        print(f"平均每1000 AA需要: {time_3000/3:.2f}秒")
    
    print("\n💡 结论:")
    if time_3000 < 5:
        print("• 标准SADO在3000 AA序列上表现良好（<5秒）")
        print("• 可以继续使用标准版本")
    elif time_3000 < 30:
        print("• 标准SADO在3000 AA序列上速度可接受")
        print("• 考虑使用增强版以获得更好的CAI")
    else:
        print("• 3000 AA序列较慢，建议使用增强版")
        print("• 或考虑减少最大迭代次数")


if __name__ == '__main__':
    main()