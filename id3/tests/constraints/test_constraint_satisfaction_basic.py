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
    """Convert RNA sequence tensor to amino acid string"""

    if rna_sequence.dim() == 3:
        # [batch, seq_len, 4] -> [batch, seq_len]
        rna_indices = torch.argmax(rna_sequence, dim=-1)

    elif rna_sequence.dim() == 2:
        # [seq_len, 4] -> [seq_len]
        rna_indices = torch.argmax(rna_sequence, dim=-1)
        rna_sequence = rna_indices

    # Convert indices to nucleotides
    nucleotides = ['A', 'C', 'G', 'U']

    # Build RNA string
    rna_str = ''.join([nucleotides[idx.item()] for idx in rna_sequence])

    # Split into codons
    codons = [rna_str[i:i+3] for i in range(0, len(rna_str), 3)]

    # Convert to amino acids
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
    """Test constraint satisfaction"""

    # Create constraint
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

    # Forward pass
    result = constraint.forward(alpha=0.0, beta=beta, tau=1.0)

    # Get sequence based on beta
    if beta == 1.0:
        # Discrete mode: use enhanced_sequence if available
        sequence = result.get('enhanced_sequence') if enable_cai else result.get('probabilities')
    else:
        # Continuous mode: use probabilities
        sequence = result.get('probabilities')

    # Convert to amino acids
    generated_amino = rna_to_amino_acid(sequence)

    # Check constraint satisfaction
    constraint_satisfied = (generated_amino == amino_acid_sequence)
    
    return constraint_satisfied, generated_amino, amino_acid_sequence

def run_constraint_tests():
    """Run constraint satisfaction tests"""

    # Test setup
    logger.info("🧪 Testing Constraint Satisfaction")
    logger.info("="*60)
    
    constraint_types = ["lagrangian", "ams", "cpc"]
    cai_settings = [True, False]
    beta_settings = [1.0, 0.0]
    
    total_tests = 0
    passed_tests = 0
    failed_details = []
    
    for constraint_type in constraint_types:
        logger.info(f"\n📊 {constraint_type.upper()} Constraint")
        logger.info("-"*40)
        
        for enable_cai in cai_settings:
            for beta in beta_settings:
                total_tests += 1

                # Test name
                cai_str = "CAI-enabled" if enable_cai else "CAI-disabled"
                beta_str = "STE" if beta == 1.0 else "soft-prob"
                test_name = f"{constraint_type}-{cai_str}-{beta_str}"
                
                try:
                    # Run test
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
                    logger.error(f"   ❌ {test_name}: Error - {str(e)}")
                    failed_details.append({
                        'test': test_name,
                        'error': str(e),
                        'constraint': constraint_type,
                        'cai': enable_cai,
                        'beta': beta
                    })

    # Summary
    logger.info("\n" + "="*60)
    logger.info(f"📊 Summary: {passed_tests}/{total_tests} tests passed")
    logger.info("="*60)

    if failed_details:
        logger.info(f"\n❌ Failed tests:")
        for fail in failed_details:
            if 'error' in fail:
                logger.error(f"   {fail['test']}: {fail['error']}")
            else:
                logger.error(f"   {fail['test']}: {fail['generated']} != {fail['target']}")

    if passed_tests == total_tests:
        logger.info(f"\n🎉 All tests passed!")
    else:
        logger.warning(f"\n⚠️ {total_tests - passed_tests} tests failed")

    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = run_constraint_tests()
    exit(0 if success else 1)