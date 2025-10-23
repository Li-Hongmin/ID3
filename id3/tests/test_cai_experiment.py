#!/usr/bin/env python3
"""


"""

import torch
import time
import json
from pathlib import Path
from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint

def test_single_constraint(constraint_type, sequence="MKTAYGGSVKP"):

    

    print("-" * 40)
    
    try:

        if constraint_type == "lagrangian":
            constraint = LagrangianConstraint(
                amino_acid_sequence=sequence,
                enable_cai=True,
                cai_target=0.8,
                device=torch.device('cpu'),
                verbose=False
            )
        elif constraint_type == "ams":
            constraint = AminoMatchingSoftmax(
                amino_acid_sequence=sequence,
                enable_cai=True,
                cai_target=0.8,
                device='cpu',
                verbose=False
            )
        elif constraint_type == "cpc":
            constraint = CodonProfileConstraint(
                amino_acid_sequence=sequence,
                enable_cai=True,
                cai_target=0.8,
                device='cpu',
                verbose=False
            )
        else:
            raise ValueError(f"Unknown constraint type: {constraint_type}")
        

        with torch.no_grad():
            result = constraint.forward(alpha=1.0 if constraint_type == "lagrangian" else 0.0, 
                                       beta=1.0 if constraint_type == "lagrangian" else 0.0, 
                                       tau=1.0)
        

        if 'discrete_cai_value' in result:
            cai_value = result['discrete_cai_value']
        elif 'final_cai' in result:
            cai_value = result['final_cai']
        else:
            cai_value = None
            
        if cai_value is not None:

            if cai_value > 0.7:

                return True, cai_value
            else:

                return False, cai_value
        else:

            return False, None
            
    except Exception as e:

        return False, None

def run_mini_experiment():
    """运行迷你12x12实验（1个蛋白质，3个约束，4个变体）"""
    
    print("\n" + "="*60)
    print("运行迷你CAI实验（1x3x4 = 12个组合）")
    print("="*60)
    

    protein = "O15263"
    constraints = ["lagrangian", "ams", "cpc"]
    variants = ["00", "01", "10", "11"]
    
    results = []
    success_count = 0
    
    for constraint in constraints:
        for variant in variants:
            config_name = f"{protein}-{constraint}-{variant}"
            print(f"\n实验 {config_name}:")
            

            success, cai_value = test_single_constraint(constraint, "MKTAYGGSVKP")
            
            if success:
                success_count += 1
                
            results.append({
                "config": config_name,
                "constraint": constraint,
                "variant": variant,
                "success": success,
                "cai_value": cai_value
            })
            
            time.sleep(0.1)
    

    print("\n" + "="*60)
    print("实验总结")
    print("="*60)
    print(f"总实验数: 12")
    print(f"成功数: {success_count}")
    print(f"成功率: {success_count/12*100:.1f}%")
    

    for constraint in constraints:
        constraint_results = [r for r in results if r["constraint"] == constraint]
        avg_cai = sum(r["cai_value"] for r in constraint_results if r["cai_value"]) / len(constraint_results)
        print(f"\n{constraint.upper()}:")
        print(f"  平均CAI: {avg_cai:.6f}")
        print(f"  成功率: {sum(r['success'] for r in constraint_results)/len(constraint_results)*100:.0f}%")
    

    output_file = f"mini_cai_experiment_{int(time.time())}.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n结果已保存到: {output_file}")
    
    return results

if __name__ == "__main__":
    print("CAI实验测试脚本")
    print("="*60)
    

    print("\n步骤1: 验证CAI功能")
    print("="*60)
    
    for constraint_type in ["lagrangian", "ams", "cpc"]:
        test_single_constraint(constraint_type)
    

    print("\n步骤2: 运行迷你12实验")
    results = run_mini_experiment()
    
    print("\n✅ 测试完成!")