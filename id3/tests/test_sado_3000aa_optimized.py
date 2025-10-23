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

print("\n" + "="*60)
print("3000 AA序列测试 - 优化后的SADO")
print("="*60)


amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
amino_sequence = ''.join(random.choice(amino_acids) for _ in range(3000))

print(f"\n序列长度: 3000 AA")
print(f"设备: {'cuda' if torch.cuda.is_available() else 'cpu'}")


print("\n生成概率分布...")
torch.manual_seed(42)
pi_accessibility = torch.rand(3000, 6)
pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


print("\n" + "-"*50)
print("1. 标准SADO (限制20次迭代):")
print("-"*50)

optimizer1 = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)


original_greedy = optimizer1._greedy_cai_optimization
def limited_greedy(indices, target_cai, max_iterations=100, pi_accessibility=None):
    return original_greedy(indices, target_cai, min(max_iterations, 20), pi_accessibility)
optimizer1._greedy_cai_optimization = limited_greedy

start = time.time()
try:
    dist1, meta1 = optimizer1.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=False,
        gamma=0.3
    )
    time1 = time.time() - start
    
    print(f"✅ 成功完成")
    print(f"  时间: {time1:.3f}秒")
    print(f"  CAI: {meta1['final_cai']:.4f}")
    print(f"  满足约束: {meta1['constraint_satisfied']}")
except Exception as e:
    time1 = time.time() - start
    print(f"❌ 失败: {str(e)[:100]}")
    print(f"  时间: {time1:.3f}秒")
    meta1 = {'final_cai': 0}


print("\n" + "-"*50)
print("2. 增强版SADO (优化后的差异驱动):")
print("-"*50)

optimizer2 = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)

start = time.time()
try:
    dist2, meta2 = optimizer2.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True
    )
    time2 = time.time() - start
    
    print(f"✅ 成功完成")
    print(f"  时间: {time2:.3f}秒")
    print(f"  CAI: {meta2['final_cai']:.4f}")
    print(f"  满足约束: {meta2['constraint_satisfied']}")
    print(f"  初始化方法: {meta2.get('init_method', 'unknown')}")
except Exception as e:
    time2 = time.time() - start
    print(f"❌ 失败: {str(e)[:100]}")
    print(f"  时间: {time2:.3f}秒")
    meta2 = {'final_cai': 0}


print("\n" + "="*60)
print("性能总结 (3000 AA序列)")
print("="*60)

print(f"\n标准SADO:  {time1:>6.3f}秒, CAI={meta1['final_cai']:.4f}")
print(f"增强版:    {time2:>6.3f}秒, CAI={meta2['final_cai']:.4f}")

if time1 > 0 and time2 > 0:
    if time1 < time2:
        print(f"\n标准版更快: {time2/time1:.2f}x")
    else:
        print(f"\n⚡ 增强版更快: {time1/time2:.2f}x")


print("\n" + "="*60)
print("💡 建议")
print("="*60)

if time2 < 10:
    print("\n✅ 优化后的增强版在3000 AA序列上表现优秀！")
    print(f"  - 执行时间: {time2:.3f}秒")
    print(f"  - CAI达到: {meta2['final_cai']:.4f}")
    print("  - 推荐使用增强版")
elif time2 < 30:
    print("\n✅ 优化后的增强版速度可接受")
    print(f"  - 执行时间: {time2:.3f}秒")
    print("  - 可根据需求选择")
else:
    print("\n⚠️ 3000 AA序列仍然较慢")
    print("  - 考虑进一步优化或分段处理")


print(f"\n预估性能:")
print(f"  5000 AA: ~{time2 * 5/3:.1f}秒")
print(f"  10000 AA: ~{time2 * 10/3:.1f}秒")