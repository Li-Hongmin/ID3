"""
Comprehensive test to demonstrate and fix the CAI enhancement operator index mismatch bug.

The bug: 
- codon_indices uses standard 64-codon indexing (0-63)
- cached_codon_indices uses constants.py indexing (0-60, no stop codons)
- This causes matching failures in compute_cai_optimal_distribution
"""

import torch
import numpy as np
import logging
import json
from pathlib import Path
from typing import Dict, List, Tuple

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import codons, amino_acids_to_codons


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def create_standard_to_constants_mapping() -> Dict[int, int]:
    """
    Create mapping from standard 64-codon indices to constants.py indices.
    
    Standard indexing: AAA=0, AAC=1, AAG=2, AAU=3, ACA=4, ..., UUU=60, UUC=61, UUG=62, UUU=63
    Constants.py: Only 61 codons (no stop codons)
    """
    nucleotides = ['A', 'C', 'G', 'U']
    standard_to_constants = {}
    
    # Build standard codon list
    standard_codons = []
    for n1 in nucleotides:
        for n2 in nucleotides:
            for n3 in nucleotides:
                codon = n1 + n2 + n3
                standard_codons.append(codon)
    
    # Map to constants.py indices
    for std_idx, codon in enumerate(standard_codons):
        if codon in codons:
            const_idx = codons.index(codon)
            standard_to_constants[std_idx] = const_idx
        else:
            # Stop codons: UAA, UAG, UGA
            standard_to_constants[std_idx] = -1  # Mark as invalid
    
    return standard_to_constants, standard_codons


def test_index_systems():
    """Test the difference between index systems."""
    print("\n" + "="*80)
    print("Testing Index Systems")
    print("="*80)
    
    standard_to_constants, standard_codons = create_standard_to_constants_mapping()
    
    # Show some examples
    test_codons = ['AUG', 'UUU', 'UUC', 'CUG', 'UAA', 'UAG', 'UGA']
    
    print("\nIndex Comparison:")
    print(f"{'Codon':<8} {'Standard':<12} {'Constants.py':<15} {'In constants?':<15}")
    print("-"*50)
    
    for codon in test_codons:
        std_idx = standard_codons.index(codon) if codon in standard_codons else -1
        const_idx = codons.index(codon) if codon in codons else -1
        in_constants = "Yes" if codon in codons else "No (stop codon)"
        print(f"{codon:<8} {std_idx:<12} {const_idx:<15} {in_constants:<15}")
    
    print(f"\nTotal codons in standard system: {len(standard_codons)}")
    print(f"Total codons in constants.py: {len(codons)}")
    print(f"Missing codons (stop codons): {[c for c in standard_codons if c not in codons]}")
    
    return standard_to_constants


def test_cai_enhancement_bug():
    """Test the CAI enhancement operator bug."""
    print("\n" + "="*80)
    print("Testing CAI Enhancement Operator Bug")
    print("="*80)
    
    # Create operator
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    operator = CAIEnhancementOperator(species='ecoli_bl21de3', device=device)
    
    # Test with simple sequences
    test_cases = [
        ("M", "Single M (AUG only)"),
        ("MM", "Double M (AUG only)"),
        ("L", "Single L (6 codons)"),
        ("LL", "Double L (6 codons)"),
        ("MLM", "Mixed M-L-M"),
        ("MLFPY", "Mixed sequence"),
    ]
    
    for amino_acid_sequence, description in test_cases:
        print(f"\n--- Testing: {description} ---")
        print(f"Amino acid sequence: {amino_acid_sequence}")
        
        # Build valid codon mask and indices
        num_positions = len(amino_acid_sequence)
        max_codons = 6
        valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=device)
        codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=device)
        
        # Fill using standard indexing
        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                aa_codons = amino_acids_to_codons[aa]
                for i, codon in enumerate(aa_codons):
                    if i < max_codons:
                        valid_codon_mask[pos, i] = True
                        codon_indices[pos, i] = operator._codon_to_standard_index(codon)
        
        # Test the buggy function
        print("\nDEBUG: Testing compute_cai_optimal_distribution...")
        cai_optimal = operator.compute_cai_optimal_distribution(
            amino_acid_sequence, valid_codon_mask, codon_indices
        )
        
        # Check if it's all zeros (bug indicator)
        is_all_zeros = torch.allclose(cai_optimal, torch.zeros_like(cai_optimal))
        print(f"Result is all zeros: {is_all_zeros}")
        
        if not is_all_zeros:
            # Compute CAI
            discrete_dist = operator.discretize_distribution(cai_optimal, valid_codon_mask)
            cai_value = operator.compute_discrete_cai(
                discrete_dist, amino_acid_sequence, valid_codon_mask, codon_indices
            )
            print(f"CAI value: {cai_value:.4f}")
            
            # Show the distribution
            for pos in range(num_positions):
                aa = amino_acid_sequence[pos]
                aa_codons = amino_acids_to_codons.get(aa, [])
                print(f"  Position {pos} ({aa}):")
                for i, codon in enumerate(aa_codons):
                    if i < max_codons:
                        weight = cai_optimal[pos, i].item()
                        print(f"    {codon}: {weight:.4f}")
        else:
            print("BUG CONFIRMED: Distribution is all zeros!")
            
            # Debug: Show what's happening in the cache
            aa = amino_acid_sequence[0]
            if aa in CAIEnhancementOperator._amino_acid_weights_cache:
                cached = CAIEnhancementOperator._amino_acid_weights_cache[aa]
                print(f"\nDEBUG: Cache for '{aa}':")
                print(f"  Cached codon indices: {cached['codon_indices'].tolist()}")
                print(f"  Cached weights: {cached['weights'].tolist()}")
                
                # Compare with actual indices
                print(f"\nDEBUG: Actual codon indices in position 0:")
                for i in range(valid_codon_mask[0].sum()):
                    print(f"  Slot {i}: index={codon_indices[0, i].item()}")
                
                print("\nDEBUG: Index mismatch detected!")
                print("The cached indices use constants.py indexing,")
                print("but codon_indices uses standard 64-codon indexing!")


def test_fixed_implementation():
    """Test the fixed implementation."""
    print("\n" + "="*80)
    print("Testing Fixed Implementation")
    print("="*80)
    
    # Create a fixed version of the operator
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    operator = CAIEnhancementOperator(species='ecoli_bl21de3', device=device)
    
    # Clear cache to ensure fresh computation
    CAIEnhancementOperator._amino_acid_weights_cache.clear()
    CAIEnhancementOperator._amino_acid_sequence_cache.clear()
    
    # Re-precompute with fixed implementation
    operator._precompute_amino_acid_weights()
    
    # Test sequences
    test_sequences = [
        "M",      # Single M - should achieve CAI ~1.0
        "MMM",    # All M - should achieve CAI ~1.0
        "L",      # Single L - should achieve good CAI
        "LLL",    # All L - should achieve good CAI
        "MLFPY",  # Mixed - should achieve reasonable CAI
        "ACDEFGHIKLMNPQRSTVWY",  # All amino acids
    ]
    
    print("\nFixed Implementation Results:")
    print(f"{'Sequence':<25} {'Max CAI':<12} {'Status':<20}")
    print("-"*60)
    
    for amino_acid_sequence in test_sequences:
        # Precompute for this sequence
        operator._load_or_compute_amino_acid_cache(amino_acid_sequence)
        
        max_cai = operator.max_achievable_cai
        status = "✓ Good" if max_cai > 0.3 else "✗ Too low"
        
        display_seq = amino_acid_sequence if len(amino_acid_sequence) <= 20 else amino_acid_sequence[:17] + "..."
        print(f"{display_seq:<25} {max_cai:.4f}        {status}")
    
    # Test CAI enhancement
    print("\n" + "="*80)
    print("Testing CAI Enhancement")
    print("="*80)
    
    test_sequence = "MLFPY"
    target_cai = 0.8
    
    print(f"\nSequence: {test_sequence}")
    print(f"Target CAI: {target_cai}")
    
    # Build inputs
    num_positions = len(test_sequence)
    max_codons = 6
    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=device)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=device)
    
    for pos, aa in enumerate(test_sequence):
        if aa in amino_acids_to_codons:
            aa_codons = amino_acids_to_codons[aa]
            for i, codon in enumerate(aa_codons):
                if i < max_codons:
                    valid_codon_mask[pos, i] = True
                    codon_indices[pos, i] = operator._codon_to_standard_index(codon)
    
    # Create a random initial distribution
    pi_accessibility = torch.rand(num_positions, max_codons, device=device)
    pi_accessibility = pi_accessibility * valid_codon_mask.float()
    # Normalize
    for pos in range(num_positions):
        if valid_codon_mask[pos].any():
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    
    # Apply CAI enhancement
    enhanced_dist, metadata = operator.apply_cai_enhancement(
        pi_accessibility, test_sequence, valid_codon_mask, codon_indices, target_cai
    )
    
    print(f"\nResults:")
    print(f"  Original CAI: {metadata['original_cai']:.4f}")
    print(f"  Target CAI: {metadata['target_cai']:.4f}")
    print(f"  Final CAI: {metadata['final_cai']:.4f}")
    print(f"  Improvement: {metadata['final_cai'] - metadata['original_cai']:.4f}")
    print(f"  Optimal gamma: {metadata['optimal_gamma']:.4f}")
    print(f"  Constraint satisfied: {metadata['constraint_satisfied']}")


def main():
    """Main test runner."""
    print("\n" + "="*80)
    print("CAI Enhancement Operator Index Mismatch Bug Test")
    print("="*80)
    
    # Test index systems
    standard_to_constants = test_index_systems()
    
    # Test the bug
    print("\n\n🔍 DEMONSTRATING THE BUG:")
    test_cai_enhancement_bug()
    
    # Test the fixed implementation
    print("\n\n✅ TESTING FIXED IMPLEMENTATION:")
    test_fixed_implementation()
    
    print("\n" + "="*80)
    print("Test Complete")
    print("="*80)


if __name__ == "__main__":
    main()