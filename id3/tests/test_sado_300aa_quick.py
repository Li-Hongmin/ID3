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

print(f"\n测试300 AA序列 - 优化后的SADO")
print("="*50)


torch.manual_seed(42)
pi_accessibility = torch.rand(300, 6)
pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


print("\n1. 标准SADO:")
optimizer1 = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)


original = optimizer1._greedy_cai_optimization
def limited(indices, target_cai, max_iterations=100, pi_accessibility=None):
    return original(indices, target_cai, min(max_iterations, 20), pi_accessibility)
optimizer1._greedy_cai_optimization = limited

start = time.time()
dist1, meta1 = optimizer1.optimize(
    pi_accessibility=pi_accessibility.clone(),
    target_cai=0.8,
    use_binary_search=False,
    use_difference_driven=False
)
time1 = time.time() - start

print(f"  时间: {time1:.3f}秒")
print(f"  CAI: {meta1['final_cai']:.4f}")


print("\n2. 增强版(优化后):")
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

print(f"  时间: {time2:.3f}秒")
print(f"  CAI: {meta2['final_cai']:.4f}")
print(f"  初始化方法: {meta2.get('init_method', 'unknown')}")

print("\n" + "="*50)
print("比较:")
print(f"  标准版: {time1:.3f}秒")
print(f"  增强版: {time2:.3f}秒")

if time1 > 0 and time2 > 0:
    if time1 < time2:
        print(f"  ✅ 标准版更快 {time2/time1:.2f}x")
    else:
        print(f"  ⚡ 增强版更快 {time1/time2:.2f}x")

print("\n结论:")
if time2 < 5:
    print("• 优化后的增强版速度可接受（<5秒）")
    if time2 < time1 * 1.5:
        print("• 性能已接近标准版")
else:
    print("• 增强版仍需进一步优化")