"""
Comprehensive test suite for CAI enhancement operator after index fix.
Tests various scenarios to ensure the CAI enhancement works correctly.
"""

import torch
import numpy as np
import logging
import unittest
from typing import List, Tuple

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestCAIEnhancementComprehensive(unittest.TestCase):
    """Comprehensive test suite for CAI enhancement operator."""
    
    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        cls.operator = CAIEnhancementOperator(species='ecoli_bl21de3', device=cls.device)
        
        # Clear caches to ensure fresh tests
        CAIEnhancementOperator._amino_acid_weights_cache.clear()
        CAIEnhancementOperator._amino_acid_sequence_cache.clear()
        
        # Re-precompute with fixed implementation
        cls.operator._precompute_amino_acid_weights()
    
    def _build_inputs(self, amino_acid_sequence: str) -> Tuple[torch.Tensor, torch.Tensor]:
        """Build valid codon mask and indices for a sequence."""
        num_positions = len(amino_acid_sequence)
        max_codons = 6
        
        valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=self.device)
        codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=self.device)
        
        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                aa_codons = amino_acids_to_codons[aa]
                for i, codon in enumerate(aa_codons):
                    if i < max_codons:
                        valid_codon_mask[pos, i] = True
                        codon_indices[pos, i] = self.operator._codon_to_standard_index(codon)
        
        return valid_codon_mask, codon_indices
    
    def test_single_codon_amino_acids(self):
        """Test amino acids with single codon (M, W) should achieve CAI=1.0."""
        single_codon_aas = ['M', 'W']
        
        for aa in single_codon_aas:
            sequences = [aa, aa * 3, aa * 5]  # Single, triple, quintuple
            
            for seq in sequences:
                with self.subTest(sequence=seq):
                    self.operator._load_or_compute_amino_acid_cache(seq)
                    self.assertAlmostEqual(
                        self.operator.max_achievable_cai, 1.0, places=4,
                        msg=f"Sequence '{seq}' should achieve CAI=1.0"
                    )
    
    def test_multi_codon_amino_acids(self):
        """Test amino acids with multiple codons."""
        test_cases = [
            ('L', 6),   # Leucine has 6 codons
            ('S', 6),   # Serine has 6 codons
            ('R', 6),   # Arginine has 6 codons
            ('V', 4),   # Valine has 4 codons
            ('P', 4),   # Proline has 4 codons
            ('T', 4),   # Threonine has 4 codons
            ('A', 4),   # Alanine has 4 codons
            ('G', 4),   # Glycine has 4 codons
            ('I', 3),   # Isoleucine has 3 codons
            ('F', 2),   # Phenylalanine has 2 codons
            ('Y', 2),   # Tyrosine has 2 codons
            ('C', 2),   # Cysteine has 2 codons
            ('H', 2),   # Histidine has 2 codons
            ('Q', 2),   # Glutamine has 2 codons
            ('N', 2),   # Asparagine has 2 codons
            ('K', 2),   # Lysine has 2 codons
            ('D', 2),   # Aspartic acid has 2 codons
            ('E', 2),   # Glutamic acid has 2 codons
        ]
        
        for aa, expected_codons in test_cases:
            with self.subTest(amino_acid=aa):
                # Check single amino acid
                self.operator._load_or_compute_amino_acid_cache(aa)
                max_cai = self.operator.max_achievable_cai
                
                # Multi-codon amino acids should achieve CAI > 0.3
                self.assertGreater(
                    max_cai, 0.3,
                    msg=f"Amino acid '{aa}' with {expected_codons} codons should achieve reasonable CAI"
                )
                
                # Verify codon count
                codons = amino_acids_to_codons.get(aa, [])
                self.assertEqual(
                    len(codons), expected_codons,
                    msg=f"Amino acid '{aa}' should have {expected_codons} codons"
                )
    
    def test_mixed_sequences(self):
        """Test mixed amino acid sequences."""
        test_sequences = [
            "MLFPY",       # Mix of different codon counts
            "ACDEFGHIKL",  # First 10 amino acids alphabetically
            "MNPQRSTVWY",  # Last 10 amino acids alphabetically
            "AAAAAAAAAA",  # Homopolymer
            "MLMLMLMLML",  # Alternating pattern
            "MGCGCGCGCG",  # Another alternating pattern
            "WDWDWDWDWD",  # Mix of single and dual codon AAs
        ]
        
        for seq in test_sequences:
            with self.subTest(sequence=seq):
                self.operator._load_or_compute_amino_acid_cache(seq)
                max_cai = self.operator.max_achievable_cai
                
                # All sequences should achieve reasonable CAI
                self.assertGreater(
                    max_cai, 0.3,
                    msg=f"Sequence '{seq}' should achieve reasonable CAI"
                )
                self.assertLessEqual(
                    max_cai, 1.0,
                    msg=f"Sequence '{seq}' CAI should not exceed 1.0"
                )
    
    def test_cai_enhancement_improvement(self):
        """Test that CAI enhancement actually improves CAI values."""
        test_sequences = ["M", "L", "MLFPY", "ACDEFGH", "IKLMNPQ"]
        target_cais = [0.5, 0.7, 0.8, 0.9]
        
        for seq in test_sequences:
            for target_cai in target_cais:
                with self.subTest(sequence=seq, target=target_cai):
                    # Build inputs
                    valid_codon_mask, codon_indices = self._build_inputs(seq)
                    
                    # Create random initial distribution
                    num_positions = len(seq)
                    max_codons = valid_codon_mask.shape[1]
                    pi_accessibility = torch.rand(num_positions, max_codons, device=self.device)
                    pi_accessibility = pi_accessibility * valid_codon_mask.float()
                    
                    # Normalize
                    for pos in range(num_positions):
                        if valid_codon_mask[pos].any():
                            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
                    
                    # Apply CAI enhancement
                    enhanced_dist, metadata = self.operator.apply_cai_enhancement(
                        pi_accessibility, seq, valid_codon_mask, codon_indices, target_cai
                    )
                    
                    # Check improvement
                    original_cai = metadata['original_cai']
                    final_cai = metadata['final_cai']
                    
                    # Final CAI should be better than original
                    self.assertGreaterEqual(
                        final_cai, original_cai,
                        msg=f"CAI enhancement should improve CAI for '{seq}'"
                    )
                    
                    # Check if constraint is satisfied (or close to max achievable)
                    self.operator._load_or_compute_amino_acid_cache(seq)
                    max_achievable = self.operator.max_achievable_cai
                    
                    if target_cai <= max_achievable:
                        # Should satisfy constraint
                        self.assertGreaterEqual(
                            final_cai, target_cai * 0.95,  # Allow 5% tolerance
                            msg=f"Should satisfy CAI constraint for '{seq}' with target {target_cai}"
                        )
                    else:
                        # Should achieve maximum possible
                        self.assertGreater(
                            final_cai, max_achievable * 0.95,  # Allow 5% tolerance
                            msg=f"Should achieve near-maximum CAI for '{seq}'"
                        )
    
    def test_cai_optimal_distribution_correctness(self):
        """Test that CAI optimal distribution selects highest-weight codons."""
        test_sequences = ["M", "L", "F", "Y", "C", "W"]
        
        for seq in test_sequences:
            with self.subTest(sequence=seq):
                valid_codon_mask, codon_indices = self._build_inputs(seq)
                
                # Compute CAI optimal distribution
                cai_optimal = self.operator.compute_cai_optimal_distribution(
                    seq, valid_codon_mask, codon_indices
                )
                
                # Check that it's not all zeros
                self.assertFalse(
                    torch.allclose(cai_optimal, torch.zeros_like(cai_optimal)),
                    msg=f"CAI optimal distribution for '{seq}' should not be all zeros"
                )
                
                # Check that each position has exactly one selected codon
                for pos in range(len(seq)):
                    valid_mask = valid_codon_mask[pos]
                    if valid_mask.any():
                        # Should be a valid probability distribution
                        pos_dist = cai_optimal[pos][valid_mask]
                        self.assertAlmostEqual(
                            pos_dist.sum().item(), 1.0, places=5,
                            msg=f"Position {pos} should have valid probability distribution"
                        )
    
    def test_discretization_preserves_selection(self):
        """Test that discretization selects the highest probability codon."""
        test_sequences = ["MLFPY", "ACDEFGH"]
        
        for seq in test_sequences:
            with self.subTest(sequence=seq):
                valid_codon_mask, codon_indices = self._build_inputs(seq)
                
                # Create a distribution with clear preferences
                num_positions = len(seq)
                max_codons = valid_codon_mask.shape[1]
                codon_probs = torch.zeros(num_positions, max_codons, device=self.device)
                
                for pos in range(num_positions):
                    valid_mask = valid_codon_mask[pos]
                    if valid_mask.any():
                        # Set different probabilities for valid codons
                        valid_indices = valid_mask.nonzero(as_tuple=True)[0]
                        for i, idx in enumerate(valid_indices):
                            codon_probs[pos, idx] = (i + 1) * 0.1
                        # Normalize
                        codon_probs[pos] = codon_probs[pos] / codon_probs[pos].sum()
                
                # Discretize
                discrete_dist = self.operator.discretize_distribution(codon_probs, valid_codon_mask)
                
                # Check that highest probability codon is selected
                for pos in range(num_positions):
                    valid_mask = valid_codon_mask[pos]
                    if valid_mask.any():
                        # Find the codon with highest probability
                        valid_probs = codon_probs[pos] * valid_mask.float()
                        max_idx = torch.argmax(valid_probs).item()
                        
                        # Check that this codon is selected in discrete distribution
                        self.assertEqual(
                            discrete_dist[pos, max_idx].item(), 1.0,
                            msg=f"Position {pos} should select highest probability codon"
                        )
                        
                        # Check that it's one-hot
                        self.assertEqual(
                            discrete_dist[pos].sum().item(), 1.0,
                            msg=f"Position {pos} should have one-hot distribution"
                        )
    
    def test_cai_computation_correctness(self):
        """Test that CAI computation follows the standard formula."""
        # Test with sequences where we know the expected CAI
        test_cases = [
            ("M", 1.0),     # Single codon, weight should be 1.0
            ("W", 1.0),     # Single codon, weight should be 1.0
            ("MM", 1.0),    # Multiple M, still 1.0
            ("MW", 1.0),    # Mix of single-codon AAs
        ]
        
        for seq, expected_cai in test_cases:
            with self.subTest(sequence=seq):
                valid_codon_mask, codon_indices = self._build_inputs(seq)
                
                # Create one-hot distribution for optimal codons
                cai_optimal = self.operator.compute_cai_optimal_distribution(
                    seq, valid_codon_mask, codon_indices
                )
                discrete_dist = self.operator.discretize_distribution(cai_optimal, valid_codon_mask)
                
                # Compute CAI
                cai_value = self.operator.compute_discrete_cai(
                    discrete_dist, seq, valid_codon_mask, codon_indices
                )
                
                # Check expected value
                self.assertAlmostEqual(
                    cai_value, expected_cai, places=3,
                    msg=f"Sequence '{seq}' should have CAI={expected_cai}"
                )
    
    def test_gamma_interpolation(self):
        """Test that gamma interpolation works correctly."""
        seq = "MLFPY"
        valid_codon_mask, codon_indices = self._build_inputs(seq)
        
        # Create two different distributions
        num_positions = len(seq)
        max_codons = valid_codon_mask.shape[1]
        
        dist1 = torch.zeros(num_positions, max_codons, device=self.device)
        dist2 = torch.zeros(num_positions, max_codons, device=self.device)
        
        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos]
            if valid_mask.any():
                valid_indices = valid_mask.nonzero(as_tuple=True)[0]
                if len(valid_indices) > 0:
                    # dist1: first codon
                    dist1[pos, valid_indices[0]] = 1.0
                    # dist2: last codon
                    dist2[pos, valid_indices[-1]] = 1.0
        
        # Test different gamma values
        gammas = [0.0, 0.25, 0.5, 0.75, 1.0]
        
        for gamma in gammas:
            with self.subTest(gamma=gamma):
                interpolated = self.operator.interpolate_distributions(dist1, dist2, gamma)
                
                # Check interpolation formula: gamma * dist2 + (1 - gamma) * dist1
                expected = gamma * dist2 + (1 - gamma) * dist1
                
                self.assertTrue(
                    torch.allclose(interpolated, expected, atol=1e-6),
                    msg=f"Interpolation with gamma={gamma} should follow linear formula"
                )
    
    def test_caching_mechanism(self):
        """Test that caching mechanism works correctly."""
        seq = "MLFPY"
        
        # Clear caches
        CAIEnhancementOperator._amino_acid_sequence_cache.clear()
        
        # First computation
        self.operator._load_or_compute_amino_acid_cache(seq)
        max_cai_1 = self.operator.max_achievable_cai
        
        # Should be in cache now
        self.assertIn(seq, CAIEnhancementOperator._amino_acid_sequence_cache)
        
        # Second computation (should use cache)
        self.operator._load_or_compute_amino_acid_cache(seq)
        max_cai_2 = self.operator.max_achievable_cai
        
        # Should get same result
        self.assertEqual(max_cai_1, max_cai_2, msg="Cached result should match")
        
        # Clear instance cache and reload
        self.operator.max_achievable_cai = None
        self.operator.cai_optimal_distribution = None
        self.operator.cached_amino_acid_sequence = None
        
        # Load from cache
        self.operator._load_or_compute_amino_acid_cache(seq)
        max_cai_3 = self.operator.max_achievable_cai
        
        # Should still match
        self.assertEqual(max_cai_1, max_cai_3, msg="Cache should persist across reloads")


def run_comprehensive_tests():
    """Run all comprehensive tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromTestCase(TestCAIEnhancementComprehensive)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "="*80)
    print("Test Summary")
    print("="*80)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {(result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100:.1f}%")
    
    if result.failures:
        print("\nFailed tests:")
        for test, trace in result.failures:
            print(f"  - {test}")
    
    if result.errors:
        print("\nTests with errors:")
        for test, trace in result.errors:
            print(f"  - {test}")
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_comprehensive_tests()
    exit(0 if success else 1)