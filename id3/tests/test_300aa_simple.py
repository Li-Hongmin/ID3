#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import time
import random

from id3.optimizers.cai.sado import SADOOptimizer



amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
amino_sequence = ''.join(random.choice(amino_acids) for _ in range(300))

print(f"\n测试300 AA序列")
print("="*50)


torch.manual_seed(42)
pi_accessibility = torch.rand(300, 6)
pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


optimizer = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)


original_greedy = optimizer._greedy_cai_optimization
def quick_greedy(indices, target_cai, max_iterations=100, pi_accessibility=None):
    return original_greedy(indices, target_cai, min(max_iterations, 10), pi_accessibility)
optimizer._greedy_cai_optimization = quick_greedy

print("\n标准SADO (最多10次迭代):")
start = time.time()
dist, meta = optimizer.optimize(
    pi_accessibility=pi_accessibility,
    target_cai=0.8,
    use_binary_search=False,
    use_difference_driven=False,
    gamma=0.3
)
elapsed = time.time() - start

print(f"  时间: {elapsed:.3f}秒")
print(f"  CAI: {meta['final_cai']:.4f}")
print(f"  满足约束: {meta['constraint_satisfied']}")

print(f"\n结论: 300 AA序列使用标准SADO需要 {elapsed:.3f}秒")
print(f"预估3000 AA序列需要: ~{elapsed*10:.1f}秒")