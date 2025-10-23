#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import hashlib

from id3.optimizers.cai.sado import SADOOptimizer



seq_length = 300
num_iterations = 20


amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))


torch.manual_seed(42)
pi_accessibility = torch.rand(seq_length, 6)
pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

print(f"\n测试重复性（{seq_length} AA, {num_iterations}次迭代）")
print("="*50)


print("\n场景1: 同一个优化器（有历史记忆）")
optimizer1 = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)

hashes1 = []
for i in range(num_iterations):
    dist, meta = optimizer1.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True
    )
    

    selected = torch.argmax(dist, dim=1)
    seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
    hashes1.append(seq_hash)
    
    if i < 5 or (i+1) % 5 == 0:
        unique = len(set(hashes1))
        print(f"  迭代{i+1:2d}: CAI={meta['final_cai']:.4f}, 唯一序列={unique}/{i+1}")

print(f"\n结果: {len(set(hashes1))}/{num_iterations} 唯一序列")
if len(set(hashes1)) == 1:
    print("❌ 所有迭代产生相同序列")
elif len(set(hashes1)) == num_iterations:
    print("✅ 每次迭代都产生不同序列")
else:
    print(f"⚠️ 有重复，重复率={(1-len(set(hashes1))/num_iterations)*100:.1f}%")


print("\n场景2: 每次新建优化器（无历史记忆）")
hashes2 = []
for i in range(10):
    optimizer2 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    dist, meta = optimizer2.optimize(
        pi_accessibility=pi_accessibility.clone(),
        target_cai=0.8,
        use_binary_search=False,
        use_difference_driven=True
    )
    
    selected = torch.argmax(dist, dim=1)
    seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
    hashes2.append(seq_hash)

print(f"结果: {len(set(hashes2))}/10 唯一序列")
if len(set(hashes2)) == 1:
    print("✅ 确定性算法：相同输入产生相同输出")
else:
    print(f"⚠️ 有随机性：产生了{len(set(hashes2))}个不同序列")


print("\n概率优化测试:")
print("-"*30)


random_prob = 0.0
for pos in range(seq_length):

    random_prob += np.log(1.0/6)
random_prob = np.exp(random_prob / seq_length)
print(f"随机选择的平均概率: {random_prob:.6f}")


optimizer3 = SADOOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=amino_sequence
)

dist, meta = optimizer3.optimize(
    pi_accessibility=pi_accessibility.clone(),
    target_cai=0.8,
    use_binary_search=False,
    use_difference_driven=True
)


selected = torch.argmax(dist, dim=1)
opt_prob = 0.0
for pos in range(seq_length):
    codon_idx = selected[pos].item()
    if codon_idx < pi_accessibility.shape[1]:
        prob = pi_accessibility[pos, codon_idx].item()
        if prob > 0:
            opt_prob += np.log(prob)
opt_prob = np.exp(opt_prob / seq_length)

print(f"优化后的平均概率: {opt_prob:.6f}")
print(f"概率提升: {opt_prob/random_prob:.2f}倍")

print("\n" + "="*50)
print("结论:")
if len(set(hashes1)) < num_iterations * 0.5:
    print("• 重复率较高，需要改进多样性机制")
else:
    print("• 重复率可接受")

if opt_prob > random_prob * 1.5:
    print("• 概率优化效果良好")
else:
    print("• 概率优化效果一般")