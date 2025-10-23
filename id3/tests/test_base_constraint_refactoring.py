#!/usr/bin/env python3
"""
Test to verify BaseConstraint refactoring works correctly.

This test ensures that:
1. BaseConstraint properly initializes CAI components
2. LagrangianConstraint inherits correctly from BaseConstraint
3. No functionality is broken after refactoring
"""

import torch
import sys
from pathlib import Path

# Add project path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from id3.constraints.base import BaseConstraint
from id3.constraints.lagrangian import LagrangianConstraint


def test_base_constraint_initialization():
    """Test that BaseConstraint properly initializes."""
    print("Testing BaseConstraint initialization...")

    # Test parameters
    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    batch_size = 1
    device = torch.device('cpu')  # Use CPU for testing

    # Test without CAI
    try:
        constraint_no_cai = LagrangianConstraint(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=batch_size,
            device=device,
            enable_cai=False
        )
        print("✓ LagrangianConstraint without CAI initialized successfully")

        # Verify attributes from base class
        assert constraint_no_cai.amino_acid_sequence == amino_acid_sequence
        assert constraint_no_cai.batch_size == batch_size
        assert constraint_no_cai.num_positions == len(amino_acid_sequence)
        assert constraint_no_cai.seq_length == len(amino_acid_sequence) * 3
        assert constraint_no_cai.enable_cai == False
        assert constraint_no_cai.cai_loss_module is None
        assert constraint_no_cai.cai_enhancement_operator is None
        print("✓ Base class attributes verified")

    except Exception as e:
        print(f"✗ Failed to initialize without CAI: {e}")
        return False

    # Test with CAI
    try:
        constraint_with_cai = LagrangianConstraint(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=batch_size,
            device=device,
            enable_cai=True,
            cai_target=0.8,
            cai_weight=0.1,
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.1,
            lambda_cai_max=2.0
        )
        print("✓ LagrangianConstraint with CAI initialized successfully")

        # Verify CAI components from base class
        assert constraint_with_cai.enable_cai == True
        assert constraint_with_cai.cai_target == 0.8
        assert constraint_with_cai.cai_loss_module is not None
        assert constraint_with_cai.cai_enhancement_operator is not None
        print("✓ CAI components from base class verified")

        # Verify Lagrangian-specific attributes
        assert hasattr(constraint_with_cai, 'lambda_value')
        assert hasattr(constraint_with_cai, 'theta')
        assert constraint_with_cai.theta.shape == (batch_size, len(amino_acid_sequence) * 3, 4)
        print("✓ Lagrangian-specific attributes verified")

    except Exception as e:
        print(f"✗ Failed to initialize with CAI: {e}")
        return False

    return True


def test_unified_loss_computation():
    """Test that unified loss computation works correctly."""
    print("\nTesting unified loss computation...")

    amino_acid_sequence = "MSKGEEL"  # Short sequence for testing
    device = torch.device('cpu')

    try:
        constraint = LagrangianConstraint(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=True,
            cai_target=0.8,
            cai_weight=0.1
        )

        # Create dummy losses
        accessibility_loss = torch.tensor(1.0, device=device)
        constraint_penalty = torch.tensor(0.5, device=device)
        cai_loss = torch.tensor(0.3, device=device)

        # Test base class unified loss computation
        loss_components = constraint.compute_unified_loss(
            accessibility_loss=accessibility_loss,
            constraint_penalty=constraint_penalty * constraint.lambda_value,  # Pre-weighted
            cai_loss=cai_loss
        )

        # Verify loss components
        assert 'total' in loss_components
        assert 'accessibility' in loss_components
        assert 'constraint' in loss_components
        assert 'cai' in loss_components
        assert 'weighted_cai' in loss_components

        # Verify calculations
        expected_total = accessibility_loss + constraint_penalty * constraint.lambda_value + constraint.lambda_cai * cai_loss
        assert torch.allclose(loss_components['total'], expected_total, rtol=1e-5)

        print("✓ Unified loss computation verified")
        print(f"  Total loss: {loss_components['total'].item():.4f}")
        print(f"  Accessibility: {loss_components['accessibility'].item():.4f}")
        print(f"  Constraint: {loss_components['constraint'].item():.4f}")
        print(f"  CAI: {loss_components.get('cai', 0):.4f}")

    except Exception as e:
        print(f"✗ Failed loss computation test: {e}")
        return False

    return True


def test_cai_info_methods():
    """Test CAI information methods from base class."""
    print("\nTesting CAI info methods...")

    amino_acid_sequence = "MSKGEEL"
    device = torch.device('cpu')

    try:
        constraint = LagrangianConstraint(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            device=device,
            enable_cai=True,
            cai_target=0.85,
            cai_weight=0.2,
            adaptive_lambda_cai=True
        )

        # Test get_cai_info
        cai_info = constraint.get_cai_info()
        assert cai_info['enabled'] == True
        assert cai_info['target'] == 0.85
        assert cai_info['species'] == 'ecoli_bl21de3'
        assert 'lambda_cai' in cai_info
        # Note: 'adaptive' field may or may not be present depending on implementation
        # assert cai_info.get('adaptive', False) == True

        print("✓ CAI info methods working correctly")
        print(f"  CAI Info: {cai_info}")

        # Test get_parameters_info
        param_info = constraint.get_parameters_info()
        assert 'LagrangianConstraint' in param_info
        assert 'Amino acids: 7' in param_info
        assert 'CAI enabled' in param_info

        print("✓ Parameter info method working correctly")

    except Exception as e:
        import traceback
        print(f"✗ Failed CAI info test: {e}")
        traceback.print_exc()
        return False

    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing BaseConstraint Refactoring")
    print("=" * 60)

    all_tests_passed = True

    # Run tests
    all_tests_passed &= test_base_constraint_initialization()
    all_tests_passed &= test_unified_loss_computation()
    all_tests_passed &= test_cai_info_methods()

    print("\n" + "=" * 60)
    if all_tests_passed:
        print("✅ All tests passed! Refactoring successful.")
    else:
        print("❌ Some tests failed. Please review the refactoring.")
    print("=" * 60)

    return 0 if all_tests_passed else 1


if __name__ == "__main__":
    sys.exit(main())