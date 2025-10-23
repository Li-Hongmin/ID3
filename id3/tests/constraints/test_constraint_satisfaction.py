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
from id3.utils.constants import amino_acids_to_codons


CODON_TO_AMINO_ACID = {}
for amino_acid, codon_list in amino_acids_to_codons.items():
    for codon in codon_list:
        CODON_TO_AMINO_ACID[codon] = amino_acid

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_discrete_satisfaction(
    constraint_type: str,
    enable_cai: bool,
    beta: float,
    amino_acid_sequence: str = "MKAI"
) -> Tuple[bool, str, str, Dict]:

    

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
    

    discrete_sequence_str = result.get('discrete_sequence', '')
    

    if not discrete_sequence_str and 'enhanced_sequence' in result:
        enhanced = result['enhanced_sequence']
        if enhanced is not None:

            if enhanced.dim() == 3:
                indices = torch.argmax(enhanced[0], dim=-1)
            elif enhanced.dim() == 2:
                indices = torch.argmax(enhanced, dim=-1)
            else:
                indices = enhanced
            
            nucleotides = ['A', 'C', 'G', 'U']
            discrete_sequence_str = ''.join([nucleotides[idx.item()] for idx in indices.cpu()])
    

    if discrete_sequence_str:
        codons = [discrete_sequence_str[i:i+3] for i in range(0, len(discrete_sequence_str), 3)]
        amino_acids = [CODON_TO_AMINO_ACID.get(codon, '?') for codon in codons if len(codon) == 3]
        generated_amino = ''.join(amino_acids)
    else:
        generated_amino = "NO_DISCRETE_SEQ"
    

    constraint_satisfied = (generated_amino == amino_acid_sequence)
    
    return constraint_satisfied, generated_amino, amino_acid_sequence, result

def run_discrete_tests():
    """运行离散序列约束满足测试"""
    
    logger.info("🧪 测试离散化序列的约束满足情况")
    logger.info("="*60)
    
    constraint_types = ["lagrangian", "ams", "cpc"]
    cai_settings = [True, False]
    beta_settings = [1.0, 0.0]
    
    total_tests = 0
    passed_tests = 0
    detailed_results = []
    
    for constraint_type in constraint_types:
        logger.info(f"\n📊 {constraint_type.upper()} 约束的离散序列")
        logger.info("-"*40)
        
        for enable_cai in cai_settings:
            for beta in beta_settings:
                total_tests += 1
                

                cai_str = "CAI启用" if enable_cai else "CAI未启用"
                beta_str = "STE" if beta == 1.0 else "软概率"
                test_name = f"{constraint_type}-{cai_str}-{beta_str}"
                
                try:

                    satisfied, generated, target, result = test_discrete_satisfaction(
                        constraint_type, enable_cai, beta
                    )
                    

                    detailed_results.append({
                        'test': test_name,
                        'satisfied': satisfied,
                        'generated': generated,
                        'target': target,
                        'has_discrete': 'discrete_sequence' in result,
                        'has_enhanced': 'enhanced_sequence' in result and result['enhanced_sequence'] is not None
                    })
                    
                    if satisfied:
                        passed_tests += 1
                        logger.info(f"   ✅ {test_name}: 离散序列 '{generated}' == '{target}'")
                    else:
                        logger.error(f"   ❌ {test_name}: 离散序列 '{generated}' != '{target}'")
                        
                except Exception as e:
                    logger.error(f"   ❌ {test_name}: 错误 - {str(e)}")
                    detailed_results.append({
                        'test': test_name,
                        'satisfied': False,
                        'error': str(e)
                    })
    

    logger.info("\n" + "="*60)
    logger.info("📊 离散序列约束满足汇总")
    logger.info("="*60)
    

    for constraint_type in constraint_types:
        constraint_results = [r for r in detailed_results if constraint_type in r['test']]
        
        logger.info(f"\n{constraint_type.upper()}:")
        total = len(constraint_results)
        passed = sum(1 for r in constraint_results if r.get('satisfied', False))
        logger.info(f"  通过率: {passed}/{total} ({100*passed/total:.0f}%)")
        
        for r in constraint_results:
            if 'error' not in r:
                status = "✅" if r['satisfied'] else "❌"
                logger.info(f"    {status} {r['test']}: {r['generated']}")
                if not r['satisfied']:
                    logger.info(f"       - has_discrete_sequence: {r['has_discrete']}")
                    logger.info(f"       - has_enhanced_sequence: {r['has_enhanced']}")
    
    logger.info("\n" + "="*60)
    logger.info(f"🎯 总结: {passed_tests}/{total_tests} ({100*passed_tests/total_tests:.1f}%) 离散序列满足约束")
    

    logger.info("\n📌 CAI增强对离散序列的影响:")
    cai_enabled = [r for r in detailed_results if 'CAI启用' in r['test'] and 'error' not in r]
    cai_disabled = [r for r in detailed_results if 'CAI未启用' in r['test'] and 'error' not in r]
    
    cai_enabled_pass = sum(1 for r in cai_enabled if r['satisfied'])
    cai_disabled_pass = sum(1 for r in cai_disabled if r['satisfied'])
    
    logger.info(f"  - CAI启用时: {cai_enabled_pass}/{len(cai_enabled)} 满足约束")
    logger.info(f"  - CAI未启用时: {cai_disabled_pass}/{len(cai_disabled)} 满足约束")
    
    if cai_enabled_pass > cai_disabled_pass:
        logger.info("  → CAI增强显著提高了约束满足率")
    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = run_discrete_tests()
    exit(0 if success else 1)