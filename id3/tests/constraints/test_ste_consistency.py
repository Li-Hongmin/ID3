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
    """Test STE consistency between probabilities and enhanced_sequence"""

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

    # Get both sequences
    rna_probs = result.get('probabilities')
    enhanced_sequence = result.get('enhanced_sequence')

    if rna_probs is None or enhanced_sequence is None:
        # Missing required outputs
        return False, float('inf')

    # Align dimensions
    if enhanced_sequence.dim() == 2 and rna_probs.dim() == 3:
        enhanced_with_batch = enhanced_sequence.unsqueeze(0)
    elif enhanced_sequence.dim() == 3 and rna_probs.dim() == 2:
        rna_probs = rna_probs.unsqueeze(0)
        enhanced_with_batch = enhanced_sequence
    else:
        enhanced_with_batch = enhanced_sequence

    # Compute L1 difference
    diff = torch.sum(torch.abs(rna_probs - enhanced_with_batch)).item()

    # Check consistency based on beta
    if beta == 1.0:
        # STE mode: should be identical
        consistent = torch.allclose(rna_probs, enhanced_with_batch, atol=1e-5)
        expected_consistent = True
    else:
        # Soft probability mode: may differ
        consistent = torch.allclose(rna_probs, enhanced_with_batch, atol=1e-5)
        expected_consistent = False
    
    test_pass = consistent == expected_consistent
    
    return test_pass, diff

def run_all_tests():
    """Run all 12 test cases"""

    logger.info("🧪 Testing STE Consistency for All Constraint Types")
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

                # Test name
                cai_str = "CAI-enabled" if enable_cai else "CAI-disabled"
                beta_str = "STE" if beta == 1.0 else "soft-prob"
                test_name = f"{constraint_type.upper()}-{cai_str}-{beta_str}"

                logger.info(f"\nTest {total_tests}/12: {test_name}")
                
                try:
                    # Run test
                    test_pass, diff = test_constraint_ste_consistency(
                        constraint_type, enable_cai, beta
                    )

                    # Store results
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
                        logger.info(f"   ✅ Passed (L1 diff: {diff:.10f})")
                    else:
                        logger.error(f"   ❌ Failed (L1 diff: {diff:.10f})")

                except Exception as e:
                    logger.error(f"   ❌ Test error: {str(e)}")
                    results.append({
                        'name': test_name,
                        'passed': False,
                        'diff': float('inf'),
                        'constraint': constraint_type,
                        'cai': enable_cai,
                        'beta': beta
                    })

    # Summary
    logger.info("\n" + "="*60)
    logger.info("📊 Test Summary")
    logger.info("="*60)

    # Per-constraint summary
    for constraint_type in constraint_types:
        logger.info(f"\n{constraint_type.upper()}:")
        constraint_results = [r for r in results if r['constraint'] == constraint_type]
        for r in constraint_results:
            status = "✅" if r['passed'] else "❌"
            logger.info(f"   {status} {r['name']}: L1 diff={r['diff']:.10f}")

    # Overall summary
    logger.info("\n" + "="*60)
    logger.info(f"🎯 Summary: {passed_tests}/{total_tests} tests passed")

    if passed_tests == total_tests:
        logger.info("🎉 All tests passed! STE consistency fix successful!")
    else:
        logger.error(f"⚠️ {total_tests - passed_tests} tests failed")
    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = run_all_tests()
    exit(0 if success else 1)