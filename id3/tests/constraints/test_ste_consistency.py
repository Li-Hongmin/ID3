#!/usr/bin/env python3
"""


"""

import torch
import logging
from typing import Dict, List, Tuple

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_constraint_ste_consistency(
    constraint_type: str,
    enable_cai: bool,
    beta: float,
    amino_acid_sequence: str = "MKAI"
) -> Tuple[bool, float]:

    

    if constraint_type == "lagrangian":
        constraint = LagrangianConstraint(
            amino_acid_sequence=amino_acid_sequence,
            species='ecoli_bl21de3',
            device='cuda',
            enable_cai=enable_cai
        )
    elif constraint_type == "ams":
        constraint = AminoMatchingSoftmax(
            amino_acid_sequence=amino_acid_sequence,
            species='ecoli_bl21de3',
            device='cuda',
            enable_cai=enable_cai
        )
    elif constraint_type == "cpc":
        constraint = CodonProfileConstraint(
            amino_acid_sequence=amino_acid_sequence,
            species='ecoli_bl21de3',
            device='cuda',
            enable_cai=enable_cai
        )
    else:
        raise ValueError(f"Unknown constraint type: {constraint_type}")
    

    result = constraint.forward(alpha=0.0, beta=beta, tau=1.0)
    

    rna_probs = result.get('probabilities')
    enhanced_sequence = result.get('enhanced_sequence')
    
    if rna_probs is None or enhanced_sequence is None:

        return False, float('inf')
    

    if enhanced_sequence.dim() == 2 and rna_probs.dim() == 3:
        enhanced_with_batch = enhanced_sequence.unsqueeze(0)
    elif enhanced_sequence.dim() == 3 and rna_probs.dim() == 2:
        rna_probs = rna_probs.unsqueeze(0)
        enhanced_with_batch = enhanced_sequence
    else:
        enhanced_with_batch = enhanced_sequence
    

    diff = torch.sum(torch.abs(rna_probs - enhanced_with_batch)).item()
    

    if beta == 1.0:

        consistent = torch.allclose(rna_probs, enhanced_with_batch, atol=1e-5)
        expected_consistent = True
    else:

        consistent = torch.allclose(rna_probs, enhanced_with_batch, atol=1e-5)
        expected_consistent = False
    
    test_pass = consistent == expected_consistent
    
    return test_pass, diff

def run_all_tests():
    """运行所有12个测试用例"""
    
    logger.info("🧪 测试所有约束类型的STE一致性")
    logger.info("="*60)
    
    constraint_types = ["lagrangian", "ams", "cpc"]
    cai_settings = [True, False]
    beta_settings = [1.0, 0.0]
    
    results = []
    total_tests = 0
    passed_tests = 0
    
    for constraint_type in constraint_types:
        for enable_cai in cai_settings:
            for beta in beta_settings:
                total_tests += 1
                

                cai_str = "CAI启用" if enable_cai else "CAI未启用"
                beta_str = "STE" if beta == 1.0 else "软概率"
                test_name = f"{constraint_type.upper()}-{cai_str}-{beta_str}"
                
                logger.info(f"\n测试 {total_tests}/12: {test_name}")
                
                try:

                    test_pass, diff = test_constraint_ste_consistency(
                        constraint_type, enable_cai, beta
                    )
                    

                    results.append({
                        'name': test_name,
                        'passed': test_pass,
                        'diff': diff,
                        'constraint': constraint_type,
                        'cai': enable_cai,
                        'beta': beta
                    })
                    
                    if test_pass:
                        passed_tests += 1
                        logger.info(f"   ✅ 通过 (L1差异: {diff:.10f})")
                    else:
                        logger.error(f"   ❌ 失败 (L1差异: {diff:.10f})")
                        
                except Exception as e:
                    logger.error(f"   ❌ 测试出错: {str(e)}")
                    results.append({
                        'name': test_name,
                        'passed': False,
                        'diff': float('inf'),
                        'constraint': constraint_type,
                        'cai': enable_cai,
                        'beta': beta
                    })
    

    logger.info("\n" + "="*60)
    logger.info("📊 测试汇总")
    logger.info("="*60)
    

    for constraint_type in constraint_types:
        logger.info(f"\n{constraint_type.upper()}:")
        constraint_results = [r for r in results if r['constraint'] == constraint_type]
        for r in constraint_results:
            status = "✅" if r['passed'] else "❌"
            logger.info(f"   {status} {r['name']}: L1差异={r['diff']:.10f}")
    

    logger.info("\n" + "="*60)
    logger.info(f"🎯 总结: {passed_tests}/{total_tests} 测试通过")
    
    if passed_tests == total_tests:
        logger.info("🎉 所有测试通过！STE一致性修复成功！")
    else:
        logger.error(f"⚠️ 有 {total_tests - passed_tests} 个测试失败")
    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)