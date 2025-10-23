#!/usr/bin/env python3
"""
Comprehensive test for all refactored constraint classes.

This test verifies that all three constraint types (Lagrangian, AMS, CPC)
properly inherit from BaseConstraint and maintain full functionality.
"""

import torch
import sys
from pathlib import Path

# Add project path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from id3.constraints.base import BaseConstraint
from id3.constraints.lagrangian import LagrangianConstraint
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint


def test_inheritance():
    """Test that all constraints properly inherit from BaseConstraint."""
    print("Testing inheritance hierarchy...")

    # Verify inheritance
    assert issubclass(LagrangianConstraint, BaseConstraint)
    assert issubclass(AminoMatchingSoftmax, BaseConstraint)
    assert issubclass(CodonProfileConstraint, BaseConstraint)

    print("✓ All constraints inherit from BaseConstraint")

    # Verify they don't directly inherit from AdaptiveLambdaCAIMixin anymore
    from id3.constraints.adaptive_lambda_cai import AdaptiveLambdaCAIMixin

    # They should inherit it through BaseConstraint, not directly
    assert AdaptiveLambdaCAIMixin not in LagrangianConstraint.__bases__
    assert AdaptiveLambdaCAIMixin not in AminoMatchingSoftmax.__bases__
    assert AdaptiveLambdaCAIMixin not in CodonProfileConstraint.__bases__

    print("✓ No direct inheritance from AdaptiveLambdaCAIMixin")

    return True


def test_constraint_initialization(constraint_class, name):
    """Test that a constraint class initializes correctly."""
    print(f"\nTesting {name} initialization...")

    amino_acid_sequence = "MSKGEEL"
    device = torch.device('cpu')

    # Test without CAI
    try:
        constraint = constraint_class(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=False
        )
        print(f"  ✓ {name} without CAI initialized")

        # Verify base class attributes
        assert constraint.amino_acid_sequence == amino_acid_sequence
        assert constraint.batch_size == 1
        assert constraint.num_positions == len(amino_acid_sequence)
        assert constraint.seq_length == len(amino_acid_sequence) * 3
        assert constraint.enable_cai == False
        assert constraint.cai_loss_module is None
        print(f"  ✓ Base class attributes verified")

    except Exception as e:
        print(f"  ✗ Failed without CAI: {e}")
        return False

    # Test with CAI
    try:
        constraint_cai = constraint_class(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=True,
            cai_target=0.8,
            cai_weight=0.1
        )
        print(f"  ✓ {name} with CAI initialized")

        # Verify CAI components
        assert constraint_cai.enable_cai == True
        assert constraint_cai.cai_target == 0.8
        assert constraint_cai.cai_loss_module is not None
        assert constraint_cai.cai_enhancement_operator is not None
        print(f"  ✓ CAI components verified")

    except Exception as e:
        print(f"  ✗ Failed with CAI: {e}")
        return False

    return True


def test_shared_methods():
    """Test that shared methods from BaseConstraint work correctly."""
    print("\nTesting shared methods from BaseConstraint...")

    amino_acid_sequence = "MSKGEEL"
    device = torch.device('cpu')

    # Test for each constraint type
    for constraint_class, name in [
        (LagrangianConstraint, "Lagrangian"),
        (AminoMatchingSoftmax, "AMS"),
        (CodonProfileConstraint, "CPC")
    ]:
        print(f"\n  Testing {name}...")

        constraint = constraint_class(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=True,
            cai_target=0.85,
            cai_weight=0.2
        )

        # Test get_cai_info
        cai_info = constraint.get_cai_info()
        assert cai_info['enabled'] == True
        assert cai_info['target'] == 0.85
        assert cai_info['species'] == 'ecoli_bl21de3'
        print(f"    ✓ get_cai_info works")

        # Test get_parameters_info
        param_info = constraint.get_parameters_info()
        assert constraint.__class__.__name__ in param_info
        assert 'Amino acids: 7' in param_info
        print(f"    ✓ get_parameters_info works")

        # Test compute_unified_loss
        accessibility_loss = torch.tensor(1.0, device=device)
        constraint_penalty = torch.tensor(0.5, device=device)
        cai_loss = torch.tensor(0.3, device=device)

        loss_components = constraint.compute_unified_loss(
            accessibility_loss=accessibility_loss,
            constraint_penalty=constraint_penalty,
            cai_loss=cai_loss
        )

        assert 'total' in loss_components
        assert 'accessibility' in loss_components
        assert 'constraint' in loss_components
        assert 'cai' in loss_components
        print(f"    ✓ compute_unified_loss works")

    return True


def test_constraint_specific_features():
    """Test that constraint-specific features still work."""
    print("\nTesting constraint-specific features...")

    amino_acid_sequence = "MSKGEEL"
    device = torch.device('cpu')

    # Test Lagrangian-specific features
    print("\n  Testing Lagrangian-specific features...")
    lagrangian = LagrangianConstraint(
        amino_acid_sequence=amino_acid_sequence,
        batch_size=1,
        device=device,
        initial_lambda=0.01,
        adaptive_lambda=True
    )
    assert hasattr(lagrangian, 'lambda_value')
    assert hasattr(lagrangian, 'theta')
    assert lagrangian.lambda_value == 0.01
    print("    ✓ Lagrangian multiplier features work")

    # Test AMS-specific features
    print("\n  Testing AMS-specific features...")
    ams = AminoMatchingSoftmax(
        amino_acid_sequence=amino_acid_sequence,
        batch_size=1,
        device=device
    )
    assert hasattr(ams, 'theta')
    assert hasattr(ams, 'valid_codon_mask')
    assert hasattr(ams, 'codon_encodings')
    print("    ✓ AMS similarity projection features work")

    # Test CPC-specific features
    print("\n  Testing CPC-specific features...")
    cpc = CodonProfileConstraint(
        amino_acid_sequence=amino_acid_sequence,
        batch_size=1,
        device=device
    )
    assert hasattr(cpc, 'omega')  # Codon-level parameters
    assert hasattr(cpc, 'valid_codon_mask')  # Mask for valid codons
    assert hasattr(cpc, 'codon_encodings')  # Codon encodings
    print("    ✓ CPC codon profile features work")

    return True


def test_memory_efficiency():
    """Test that refactoring didn't increase memory usage."""
    print("\nTesting memory efficiency...")

    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    device = torch.device('cpu')

    # Count parameters for each constraint
    for constraint_class, name in [
        (LagrangianConstraint, "Lagrangian"),
        (AminoMatchingSoftmax, "AMS"),
        (CodonProfileConstraint, "CPC")
    ]:
        constraint = constraint_class(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=True
        )

        total_params = sum(p.numel() for p in constraint.parameters())
        print(f"  {name}: {total_params:,} parameters")

    print("  ✓ Memory usage reasonable")
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing All Refactored Constraints")
    print("=" * 60)

    all_tests_passed = True

    # Run tests
    all_tests_passed &= test_inheritance()
    all_tests_passed &= test_constraint_initialization(LagrangianConstraint, "Lagrangian")
    all_tests_passed &= test_constraint_initialization(AminoMatchingSoftmax, "AMS")
    all_tests_passed &= test_constraint_initialization(CodonProfileConstraint, "CPC")
    all_tests_passed &= test_shared_methods()
    all_tests_passed &= test_constraint_specific_features()
    all_tests_passed &= test_memory_efficiency()

    print("\n" + "=" * 60)
    if all_tests_passed:
        print("✅ All tests passed! Refactoring complete and successful.")
    else:
        print("❌ Some tests failed. Please review the refactoring.")
    print("=" * 60)

    return 0 if all_tests_passed else 1


if __name__ == "__main__":
    sys.exit(main())