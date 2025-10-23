#!/usr/bin/env python3
"""


"""

import torch
import logging
from typing import Dict, List


from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_constraint_ste_consistency(constraint_class, constraint_name: str, amino_acid_sequence: str) -> Dict:

    
    logger.info(f"\n{'='*60}")

    logger.info(f"{'='*60}")
    
    test_configs = [




    ]
    
    results = []
    constraint_results = {"constraint": constraint_name, "tests": []}
    
    for config in test_configs:

        logger.info(f"-" * 40)
        
        try:

            constraint = constraint_class(
                amino_acid_sequence=amino_acid_sequence,
                species='ecoli_bl21de3',
                device='cuda',
                enable_cai=config['enable_cai']
            )
            

            if constraint_name == "Lagrangian":
                result = constraint.forward(
                    alpha=0.0, 
                    beta=config['beta'], 
                    tau=1.0, 
                    compute_penalty=True
                )
            else:
                result = constraint.forward(
                    alpha=0.0, 
                    beta=config['beta'], 
                    tau=1.0
                )
            
            rna_probs = result.get('probabilities')
            enhanced_sequence = result.get('enhanced_sequence')
            
            if rna_probs is None or enhanced_sequence is None:

                test_result = {
                    'name': config['name'],
                    'consistent': False,
                    'expected': config['expected_consistent'],
                    'diff': float('inf'),
                    'result': 'FAIL',
                    'error': 'incomplete_data'
                }
            else:

                if enhanced_sequence.dim() == 2 and rna_probs.dim() == 3:
                    enhanced_with_batch = enhanced_sequence.unsqueeze(0)
                else:
                    enhanced_with_batch = enhanced_sequence
                

                diff = torch.sum(torch.abs(rna_probs - enhanced_with_batch)).item()
                is_consistent = torch.allclose(rna_probs, enhanced_with_batch, atol=1e-5)
                





                

                if is_consistent == config['expected_consistent']:

                    result_status = "PASS"
                else:



                    result_status = "FAIL"
                
                test_result = {
                    'name': config['name'],
                    'consistent': is_consistent,
                    'expected': config['expected_consistent'],
                    'diff': diff,
                    'result': result_status,
                    'error': None
                }
        
        except Exception as e:

            test_result = {
                'name': config['name'],
                'consistent': False,
                'expected': config['expected_consistent'],
                'diff': float('inf'),
                'result': 'ERROR',
                'error': str(e)
            }
        
        constraint_results["tests"].append(test_result)
    

    passed = sum(1 for t in constraint_results["tests"] if t['result'] == 'PASS')
    total = len(constraint_results["tests"])
    

    for test in constraint_results["tests"]:
        status = "✅" if test['result'] == 'PASS' else ("❌" if test['result'] == 'FAIL' else "💥")
        logger.info(f"   {status} {test['name']}: {test['result']}")
        if test['error']:

    
    constraint_results["passed"] = passed
    constraint_results["total"] = total
    constraint_results["success_rate"] = passed / total if total > 0 else 0.0
    

    
    return constraint_results

def main():
    """主测试函数"""
    
    logger.info("🚀 开始验证所有约束类型的STE+CAI一致性")
    

    amino_acid_sequence = "MKAI"
    

    constraint_types = [
        (LagrangianConstraint, "Lagrangian"),
        (AminoMatchingSoftmax, "AMS"),
        (CodonProfileConstraint, "CPC"),
    ]
    
    all_results = []
    
    for constraint_class, constraint_name in constraint_types:
        try:
            result = test_constraint_ste_consistency(constraint_class, constraint_name, amino_acid_sequence)
            all_results.append(result)
        except Exception as e:
            logger.error(f"💥 {constraint_name} 约束测试完全失败: {str(e)}")
            all_results.append({
                "constraint": constraint_name,
                "tests": [],
                "passed": 0,
                "total": 0,
                "success_rate": 0.0,
                "error": str(e)
            })
    

    logger.info(f"\n{'='*80}")
    logger.info(f"📊 全局测试汇总")
    logger.info(f"{'='*80}")
    
    total_passed = 0
    total_tests = 0
    
    for result in all_results:
        constraint_name = result["constraint"]
        passed = result["passed"]
        total = result["total"]
        success_rate = result["success_rate"]
        
        status = "✅" if success_rate >= 1.0 else ("⚠️" if success_rate >= 0.75 else "❌")
        logger.info(f"{status} {constraint_name}: {passed}/{total} ({success_rate:.1%})")
        
        if "error" in result:
            logger.info(f"   💥 错误: {result['error']}")
        

        for test in result.get("tests", []):
            if test["result"] != "PASS":
                logger.warning(f"   ⚠️ {test['name']}: {test['result']}")
                if test.get("error"):
                    logger.warning(f"      错误: {test['error']}")
                elif test["result"] == "FAIL":
                    logger.warning(f"      预期: {test['expected']}, 实际: {test['consistent']}, 差异: {test['diff']:.2e}")
        
        total_passed += passed
        total_tests += total
    
    overall_success_rate = total_passed / total_tests if total_tests > 0 else 0.0
    
    logger.info(f"\n🎯 总体结果: {total_passed}/{total_tests} ({overall_success_rate:.1%})")
    

    problematic_constraints = [r for r in all_results if r["success_rate"] < 1.0]
    
    if problematic_constraints:
        logger.warning(f"\n⚠️ 发现问题的约束类型:")
        for result in problematic_constraints:
            logger.warning(f"   • {result['constraint']}: STE+CAI一致性问题")
            

            ste_cai_issues = [t for t in result.get("tests", []) if "CAI启用+STE" in t["name"] and t["result"] != "PASS"]
            ste_no_cai_issues = [t for t in result.get("tests", []) if "CAI未启用+STE" in t["name"] and t["result"] != "PASS"]
            
            if ste_cai_issues:
                logger.warning(f"     - CAI启用+STE模式有问题")
            if ste_no_cai_issues:
                logger.warning(f"     - CAI未启用+STE模式有问题")
    else:
        logger.info(f"\n🎉 所有约束类型的STE一致性都正常!")
    
    return overall_success_rate >= 1.0

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)