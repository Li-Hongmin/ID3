"""
Quick test to verify CAI enhancement operator type fixes work correctly.
"""

import torch
import sys
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent.parent))

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator, create_cai_enhancement_operator


def test_cai_enhancement_type_fixes():
    """Test that type fixes in CAI enhancement operator work correctly."""
    
    # Test 1: Create with None device (should use Optional[torch.device])
    operator1 = CAIEnhancementOperator(species='ecoli_bl21de3', device=None)
    assert operator1.device == torch.device('cpu')
    print("✓ Test 1 passed: None device handled correctly")
    
    # Test 2: Create with explicit device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    operator2 = CAIEnhancementOperator(species='ecoli_bl21de3', device=device)
    assert operator2.device == device
    print("✓ Test 2 passed: Explicit device handled correctly")
    
    # Test 3: Factory function with None device
    operator3 = create_cai_enhancement_operator(species='ecoli_bl21de3', device=None)
    assert operator3.device == torch.device('cpu')
    print("✓ Test 3 passed: Factory function with None device works")
    
    # Test 4: Test integer indexing (the main type fix)
    # Create dummy data
    num_positions = 10
    max_codons = 6
    
    codon_probs = torch.randn(num_positions, max_codons)
    codon_probs = torch.softmax(codon_probs, dim=-1)
    valid_mask = torch.ones(num_positions, max_codons, dtype=torch.bool)
    
    # Test discretization (which uses int() conversion)
    discrete_dist = operator1.discretize_distribution(codon_probs, valid_mask)
    
    # Check that it's a valid one-hot distribution
    assert discrete_dist.shape == (num_positions, max_codons)
    assert torch.allclose(discrete_dist.sum(dim=1), torch.ones(num_positions))
    assert torch.all((discrete_dist == 0) | (discrete_dist == 1))
    print("✓ Test 4 passed: Discretization with integer indexing works")
    
    # Test 5: Test CAI computation with integer indexing
    amino_acid_sequence = "MSKGEELFTG"  # Sample sequence
    codon_indices = torch.randint(0, 61, (num_positions, max_codons))
    
    cai = operator1.compute_discrete_cai(
        discrete_dist, 
        amino_acid_sequence,
        valid_mask,
        codon_indices
    )
    
    assert isinstance(cai, float)
    assert 0 <= cai <= 1
    print(f"✓ Test 5 passed: CAI computation works (CAI = {cai:.4f})")
    
    print("\n✅ All type fix tests passed successfully!")
    return True


if __name__ == "__main__":
    try:
        test_cai_enhancement_type_fixes()
    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)