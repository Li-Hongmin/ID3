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


def generate_random_amino_sequence(length: int) -> str:

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))


def test_3000aa_sequence():
    """测试3000 AA序列的性能"""
    
    print("\n" + "="*80)
    print("3000 AA长序列SADO性能测试")
    print("="*80 + "\n")
    

    seq_length = 3000
    amino_sequence = generate_random_amino_sequence(seq_length)
    print(f"序列长度: {seq_length} AA")
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"使用设备: {device}\n")
    

    print("-"*60)
    print("测试1: 标准SADO (Gamma初始化 + 效率比值优化)")
    print("-"*60)
    
    optimizer_standard = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    start_time = time.time()
    try:
        dist_standard, meta_standard = optimizer_standard.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False,
            gamma=0.3
        )
        time_standard = time.time() - start_time
        
        print(f"✅ 成功完成")
        print(f"  最终CAI: {meta_standard['final_cai']:.4f}")
        print(f"  目标满足: {meta_standard['constraint_satisfied']}")
        print(f"  执行时间: {time_standard:.2f}秒")
        print(f"  初始化方法: {meta_standard.get('init_method', 'unknown')}")
    except Exception as e:
        time_standard = time.time() - start_time
        print(f"❌ 失败: {str(e)[:100]}")
        print(f"  运行时间: {time_standard:.2f}秒")
        meta_standard = {'final_cai': 0, 'constraint_satisfied': False}
    

    print("\n" + "-"*60)
    print("测试2: 增强版SADO (平衡状态 + 差异驱动)")
    print("-"*60)
    
    optimizer_enhanced = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    start_time = time.time()
    try:
        dist_enhanced, meta_enhanced = optimizer_enhanced.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True,
            gamma=0.3
        )
        time_enhanced = time.time() - start_time
        
        print(f"✅ 成功完成")
        print(f"  最终CAI: {meta_enhanced['final_cai']:.4f}")
        print(f"  目标满足: {meta_enhanced['constraint_satisfied']}")
        print(f"  执行时间: {time_enhanced:.2f}秒")
        print(f"  初始化方法: {meta_enhanced.get('init_method', 'unknown')}")
    except Exception as e:
        time_enhanced = time.time() - start_time
        print(f"❌ 失败: {str(e)[:100]}")
        print(f"  运行时间: {time_enhanced:.2f}秒")
        meta_enhanced = {'final_cai': 0, 'constraint_satisfied': False}
    

    print("\n" + "-"*60)
    print("测试3: 完整增强版 (二分查找 + 差异驱动)")
    print("-"*60)
    
    optimizer_full = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    start_time = time.time()
    try:
        dist_full, meta_full = optimizer_full.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=True,
            use_difference_driven=True,
            gamma=0.3
        )
        time_full = time.time() - start_time
        
        print(f"✅ 成功完成")
        print(f"  最终CAI: {meta_full['final_cai']:.4f}")
        print(f"  目标满足: {meta_full['constraint_satisfied']}")
        print(f"  执行时间: {time_full:.2f}秒")
        print(f"  初始化方法: {meta_full.get('init_method', 'unknown')}")
        if 'binary_search_alpha' in meta_full:
            print(f"  二分查找α: {meta_full['binary_search_alpha']:.4f}")
    except Exception as e:
        time_full = time.time() - start_time
        print(f"❌ 失败: {str(e)[:100]}")
        print(f"  运行时间: {time_full:.2f}秒")
        meta_full = {'final_cai': 0, 'constraint_satisfied': False}
    

    print("\n" + "="*60)
    print("性能总结 (3000 AA序列)")
    print("="*60)
    
    print(f"\n{'方法':<30} {'CAI':<10} {'时间(秒)':<12} {'速度'}")
    print("-"*60)
    
    print(f"{'标准SADO':<30} {meta_standard['final_cai']:<10.4f} {time_standard:<12.2f} {'基准'}")
    
    if time_standard > 0:
        speed_enhanced = time_standard / time_enhanced if time_enhanced > 0 else 0
        speed_full = time_standard / time_full if time_full > 0 else 0
        print(f"{'增强版(平衡+差异)':<30} {meta_enhanced['final_cai']:<10.4f} {time_enhanced:<12.2f} {speed_enhanced:.2f}x")
        print(f"{'完整增强(二分+差异)':<30} {meta_full['final_cai']:<10.4f} {time_full:<12.2f} {speed_full:.2f}x")
    

    print("\n" + "="*60)
    print("💡 建议")
    print("="*60)
    
    if time_standard < time_enhanced and time_standard < time_full:
        print("\n✅ 对于3000 AA序列，标准SADO表现最佳")
        print(f"  - 速度最快: {time_standard:.2f}秒")
        print(f"  - CAI达标: {meta_standard['final_cai']:.4f}")
    elif time_enhanced < time_full:
        print("\n✅ 对于3000 AA序列，增强版SADO（平衡+差异）表现最佳")
        print(f"  - 速度较快: {time_enhanced:.2f}秒")
        print(f"  - CAI达标: {meta_enhanced['final_cai']:.4f}")
    else:
        print("\n✅ 对于3000 AA序列，完整增强版表现最佳")
        print(f"  - 最高CAI: {meta_full['final_cai']:.4f}")
        print(f"  - 时间合理: {time_full:.2f}秒")


def test_multiple_iterations():

    
    print("\n" + "="*60)

    print("="*60 + "\n")
    
    seq_length = 500
    amino_sequence = generate_random_amino_sequence(seq_length)
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        amino_acid_sequence=amino_sequence
    )
    

    times = []
    cai_values = []
    
    for i in range(5):
        start = time.time()
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False
        )
        elapsed = time.time() - start
        times.append(elapsed)
        cai_values.append(meta['final_cai'])

    




def main():
    """主测试函数"""
    
    try:

        test_3000aa_sequence()
        

        test_multiple_iterations()
        
        print("\n" + "="*60)
        print("✅ 所有测试完成!")
        print("="*60 + "\n")
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()