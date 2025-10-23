#!/usr/bin/env python3
"""


"""

import torch
import logging
from typing import Dict, List, Tuple

from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.utils.constants import amino_acids_to_codons


CODON_TO_AMINO_ACID = {}
for amino_acid, codon_list in amino_acids_to_codons.items():
    for codon in codon_list:
        CODON_TO_AMINO_ACID[codon] = amino_acid

CODON_TO_AMINO_ACID['UAA'] = '*'
CODON_TO_AMINO_ACID['UAG'] = '*'
CODON_TO_AMINO_ACID['UGA'] = '*'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def rna_to_amino_acid(rna_sequence: torch.Tensor) -> str:


    if rna_sequence.dim() == 3:
        # [batch, seq_len, 4] -> [batch, seq_len]
        rna_indices = torch.argmax(rna_sequence, dim=-1)

    elif rna_sequence.dim() == 2:
        # [seq_len, 4] -> [seq_len]
        rna_indices = torch.argmax(rna_sequence, dim=-1)
        rna_sequence = rna_indices
    

    nucleotides = ['A', 'C', 'G', 'U']
    

    rna_str = ''.join([nucleotides[idx.item()] for idx in rna_sequence])
    

    codons = [rna_str[i:i+3] for i in range(0, len(rna_str), 3)]
    

    amino_acids = []
    for codon in codons:
        if len(codon) == 3:
            amino_acids.append(CODON_TO_AMINO_ACID.get(codon, '?'))
    
    return ''.join(amino_acids)

def test_constraint_satisfaction(
    constraint_type: str,
    enable_cai: bool,
    beta: float,
    amino_acid_sequence: str = "MKAI"
) -> Tuple[bool, str, str]:
    """测试约束满足情况"""
    

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
    

    if beta == 1.0:

        sequence = result.get('enhanced_sequence') if enable_cai else result.get('probabilities')
    else:

        sequence = result.get('probabilities')
    

    generated_amino = rna_to_amino_acid(sequence)
    

    constraint_satisfied = (generated_amino == amino_acid_sequence)
    
    return constraint_satisfied, generated_amino, amino_acid_sequence

def run_constraint_tests():

    

    logger.info("="*60)
    
    constraint_types = ["lagrangian", "ams", "cpc"]
    cai_settings = [True, False]
    beta_settings = [1.0, 0.0]
    
    total_tests = 0
    passed_tests = 0
    failed_details = []
    
    for constraint_type in constraint_types:

        logger.info("-"*40)
        
        for enable_cai in cai_settings:
            for beta in beta_settings:
                total_tests += 1
                



                test_name = f"{constraint_type}-{cai_str}-{beta_str}"
                
                try:

                    satisfied, generated, target = test_constraint_satisfaction(
                        constraint_type, enable_cai, beta
                    )
                    
                    if satisfied:
                        passed_tests += 1
                        logger.info(f"   ✅ {test_name}: {generated} == {target}")
                    else:
                        logger.error(f"   ❌ {test_name}: {generated} != {target}")
                        failed_details.append({
                            'test': test_name,
                            'generated': generated,
                            'target': target,
                            'constraint': constraint_type,
                            'cai': enable_cai,
                            'beta': beta
                        })
                        
                except Exception as e:

                    failed_details.append({
                        'test': test_name,
                        'error': str(e),
                        'constraint': constraint_type,
                        'cai': enable_cai,
                        'beta': beta
                    })
    

    logger.info("\n" + "="*60)

    logger.info("="*60)

    
    if failed_details:

        for fail in failed_details:
            if 'error' in fail:
                logger.error(f"   {fail['test']}: {fail['error']}")
            else:

    
    if passed_tests == total_tests:

    else:

    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = run_constraint_tests()
    exit(0 if success else 1)