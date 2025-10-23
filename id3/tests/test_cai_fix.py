#!/usr/bin/env python3

import torch
import sys
import os
sys.path.insert(0, '/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.lagrangian import LagrangianConstraint
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper

print("🧪 测试Lagrangian约束的CAI修复效果")
print("=" * 60)


print("📦 初始化DeepRaccess模型...")
deepraccess = DeepRaccessID3Wrapper()


test_sequence = "MRVLYLLFSFLFIFLMPLPG"
print(f"🧬 测试序列: {test_sequence}")


print("🔧 创建Lagrangian约束...")
constraint = LagrangianConstraint(
    amino_acid_sequence=test_sequence,
    enable_cai=True,
    cai_target=0.8,
    cai_weight=1.0,
    verbose=True
)

print("\n🔄 运行优化测试...")
print("=" * 40)


results = []
for iteration in range(5):
    print(f"\n📊 迭代 {iteration + 1}/5:")
    

    result = constraint.forward(beta=1.0)
    

    discrete_cai = result.get('discrete_cai_value', 'N/A')
    if isinstance(discrete_cai, torch.Tensor):
        discrete_cai = discrete_cai.item()
    
    accessibility = result.get('sequence', None)
    if accessibility is not None:

        with torch.no_grad():
            atg_pos = test_sequence.find('M')
            accessibility_loss = deepraccess.compute_atg_window_accessibility(
                accessibility, atg_position=atg_pos, discrete=False
            )
            if isinstance(accessibility_loss, torch.Tensor):
                accessibility_loss = accessibility_loss.item()
    else:
        accessibility_loss = float('inf')
    
    results.append({
        'iteration': iteration + 1,
        'discrete_cai': discrete_cai,
        'accessibility_loss': accessibility_loss
    })
    
    print(f"   🎯 离散CAI值: {discrete_cai}")
    print(f"   🔬 可达性损失: {accessibility_loss:.4f}")
    

    if isinstance(discrete_cai, (int, float)) and discrete_cai >= 0.8:
        print(f"   ✅ CAI目标达成! ({discrete_cai:.4f} >= 0.8)")
    else:
        print(f"   ⚠️  CAI未达目标 ({discrete_cai} < 0.8)")

print("\n" + "=" * 60)
print("📈 测试结果总结:")
print("=" * 60)

successful_cai = 0
total_tests = len(results)

for r in results:
    cai_val = r['discrete_cai']
    if isinstance(cai_val, (int, float)) and cai_val >= 0.8:
        successful_cai += 1

print(f"🎯 CAI目标达成率: {successful_cai}/{total_tests} ({100*successful_cai/total_tests:.1f}%)")

if successful_cai > 0:
    print("✅ 修复成功! Lagrangian约束现在能够产生高CAI值")
    max_cai = max([r['discrete_cai'] for r in results if isinstance(r['discrete_cai'], (int, float))])
    print(f"🏆 最高CAI值: {max_cai:.4f}")
else:
    print("❌ 修复可能不完整，需要进一步调试")

print("\n🎉 测试完成!")