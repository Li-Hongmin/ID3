#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import time
import random
import logging


logging.basicConfig(level=logging.DEBUG)

from id3.optimizers.cai.sado import SADOOptimizer

print("\n测试3000 AA序列（仅测试增强版）")
print("="*50)


amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
amino_sequence = ''.join(random.choice(amino_acids) for _ in range(3000))

print(f"序列长度: 3000 AA")


torch.manual_seed(42)
pi_accessibility = torch.rand(3000, 6)
pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


print("\n测试增强版SADO...")

optimizer = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)

print("开始优化...")
start = time.time()

try:
    dist, meta = optimizer.optimize(
        pi_accessibility=pi_accessibility,
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True,
        gamma=0.3
    )
    elapsed = time.time() - start
    
    print(f"\n✅ 成功完成!")
    print(f"  时间: {elapsed:.3f}秒")
    print(f"  CAI: {meta['final_cai']:.4f}")
    print(f"  初始化方法: {meta.get('init_method', 'unknown')}")
    
except Exception as e:
    elapsed = time.time() - start
    print(f"\n❌ 失败!")
    print(f"  错误: {str(e)}")
    print(f"  时间: {elapsed:.3f}秒")
    
    import traceback
    traceback.print_exc()

print(f"\n结论: 3000 AA序列运行时间={elapsed:.3f}秒")