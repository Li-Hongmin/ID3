#!/usr/bin/env python3
"""


"""

import torch
import sys
import os
from pathlib import Path


sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from id3.constraints.lagrangian import LagrangianConstraint
from id3.utils.sequence_utils import rna_to_amino_acids

def test_lagrangian_discrete_without_cai():

    

    print("=" * 60)
    



    

    constraint = LagrangianConstraint(
        amino_acid_sequence=amino_acid_sequence,

        lambda_init=0.1
    )
    




    print()
    

    test_cases = [



    ]
    
    all_passed = True
    
    for case in test_cases:


        

        result = constraint.forward(
            alpha=case['alpha'],
            beta=case['beta'],
            tau=1.0,
            compute_penalty=True
        )
        

        discrete_sequence = result.get('discrete_sequence', '')
        
        if not discrete_sequence:

            all_passed = False
            continue
            

        if len(discrete_sequence) != seq_len:

            all_passed = False
            continue
            

        translated = rna_to_amino_acids(discrete_sequence)
        constraint_satisfied = (translated == amino_acid_sequence)
        
        if constraint_satisfied:



        else:



            all_passed = False
            

        if 'enhanced_sequence' in result:
            enhanced = result['enhanced_sequence']

        
        print()
    

    print("=" * 60)
    if all_passed:


    else:

    
    return all_passed

def test_with_cai_comparison():
    """对比测试：CAI启用 vs 未启用"""
    
    print("\n" + "=" * 60)
    print("📊 对比测试: CAI启用 vs 未启用")
    print("=" * 60)
    
    amino_acid_sequence = "MSTGAV"
    

    print("\n1️⃣ CAI未启用:")
    constraint_no_cai = LagrangianConstraint(
        amino_acid_sequence=amino_acid_sequence,
        enable_cai=False,
        lambda_init=0.1
    )
    
    result_no_cai = constraint_no_cai.forward(alpha=0.0, beta=0.0, tau=1.0)
    discrete_no_cai = result_no_cai.get('discrete_sequence', '')
    translated_no_cai = rna_to_amino_acids(discrete_no_cai)
    match_no_cai = (translated_no_cai == amino_acid_sequence)
    
    print(f"   离散序列: {discrete_no_cai}")
    print(f"   翻译结果: {translated_no_cai}")
    print(f"   约束满足: {'✅ 是' if match_no_cai else '❌ 否'}")
    

    print("\n2️⃣ CAI启用:")
    constraint_with_cai = LagrangianConstraint(
        amino_acid_sequence=amino_acid_sequence,
        enable_cai=True,
        lambda_init=0.1,
        cai_target=0.8,
        lambda_cai=0.1
    )
    
    result_with_cai = constraint_with_cai.forward(alpha=0.0, beta=0.0, tau=1.0)
    discrete_with_cai = result_with_cai.get('discrete_sequence', '')
    translated_with_cai = rna_to_amino_acids(discrete_with_cai)
    match_with_cai = (translated_with_cai == amino_acid_sequence)
    
    print(f"   离散序列: {discrete_with_cai}")
    print(f"   翻译结果: {translated_with_cai}")
    print(f"   约束满足: {'✅ 是' if match_with_cai else '❌ 否'}")
    

    print("\n📈 对比结果:")
    print(f"   CAI未启用约束满足: {'✅' if match_no_cai else '❌'}")
    print(f"   CAI启用约束满足: {'✅' if match_with_cai else '❌'}")
    
    if match_no_cai and match_with_cai:
        print("\n🎉 修复验证成功！两种模式下离散序列都满足约束")
        return True
    else:
        print("\n⚠️ 存在问题，需要进一步调试")
        return False

if __name__ == "__main__":
    print("🚀 开始测试Lagrangian离散路径修复")
    print("=" * 60)
    

    test1_passed = test_lagrangian_discrete_without_cai()
    

    test2_passed = test_with_cai_comparison()
    

    print("\n" + "=" * 60)
    if test1_passed and test2_passed:
        print("🎊 所有测试通过！修复完全成功")
        print("   Lagrangian现在在所有情况下都生成满足约束的离散序列")
    else:
        print("⚠️ 部分测试失败，请检查修复实现")