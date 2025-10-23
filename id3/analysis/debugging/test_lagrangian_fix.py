#!/usr/bin/env python3
"""
Test script to verify Lagrangian constraint Gumbel noise fix
"""

import torch
import sys
from pathlib import Path

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.constraints.lagrangian import SimplePiFunction

def test_gumbel_noise():
    """Test that Gumbel noise is properly added when alpha > 0"""
    
    print("=" * 80)
    print("Testing SimplePiFunction with Gumbel Noise Fix")
    print("=" * 80)
    
    # Create some test logits
    batch_size = 2
    seq_length = 9
    vocab_size = 4
    
    theta = torch.randn(batch_size, seq_length, vocab_size)
    
    # Test 1: Without noise (alpha=0)
    print("\n1. Testing without noise (alpha=0):")
    prob1 = SimplePiFunction.forward(theta, alpha=0.0, tau=1.0)
    prob2 = SimplePiFunction.forward(theta, alpha=0.0, tau=1.0)
    
    # Should be identical
    diff_no_noise = (prob1 - prob2).abs().max().item()
    print(f"   Max difference between two runs: {diff_no_noise:.10f}")
    print(f"   ✓ PASS: Deterministic when alpha=0" if diff_no_noise < 1e-6 else "   ✗ FAIL: Should be deterministic")
    
    # Test 2: With noise (alpha=1)
    print("\n2. Testing with noise (alpha=1):")
    prob3 = SimplePiFunction.forward(theta, alpha=1.0, tau=1.0)
    prob4 = SimplePiFunction.forward(theta, alpha=1.0, tau=1.0)
    
    # Should be different
    diff_with_noise = (prob3 - prob4).abs().max().item()
    print(f"   Max difference between two runs: {diff_with_noise:.10f}")
    print(f"   ✓ PASS: Stochastic when alpha=1" if diff_with_noise > 1e-4 else "   ✗ FAIL: Should be stochastic")
    
    # Test 3: Multiple runs to check diversity
    print("\n3. Testing sequence diversity with multiple runs:")
    num_runs = 100
    sequences = []
    
    for _ in range(num_runs):
        prob = SimplePiFunction.forward(theta, alpha=1.0, tau=1.0)
        # Convert to discrete sequence
        seq = prob.argmax(dim=-1)
        sequences.append(seq[0].tolist())  # Take first batch element
    
    # Count unique sequences
    unique_sequences = len(set(tuple(s) for s in sequences))
    uniqueness_rate = unique_sequences / num_runs
    
    print(f"   Unique sequences: {unique_sequences}/{num_runs}")
    print(f"   Uniqueness rate: {uniqueness_rate:.1%}")
    print(f"   ✓ PASS: Good diversity" if uniqueness_rate > 0.5 else "   ✗ FAIL: Poor diversity")
    
    # Test 4: Verify temperature effect
    print("\n4. Testing temperature effect:")
    prob_high_temp = SimplePiFunction.forward(theta, alpha=0.0, tau=10.0)
    prob_low_temp = SimplePiFunction.forward(theta, alpha=0.0, tau=0.1)
    
    entropy_high = -(prob_high_temp * (prob_high_temp + 1e-10).log()).sum(-1).mean().item()
    entropy_low = -(prob_low_temp * (prob_low_temp + 1e-10).log()).sum(-1).mean().item()
    
    print(f"   Entropy at tau=10.0: {entropy_high:.4f}")
    print(f"   Entropy at tau=0.1: {entropy_low:.4f}")
    print(f"   ✓ PASS: Temperature works correctly" if entropy_high > entropy_low else "   ✗ FAIL: Temperature effect incorrect")
    
    print("\n" + "=" * 80)
    print("Summary:")
    
    all_pass = (
        diff_no_noise < 1e-6 and 
        diff_with_noise > 1e-4 and 
        uniqueness_rate > 0.5 and 
        entropy_high > entropy_low
    )
    
    if all_pass:
        print("✅ ALL TESTS PASSED - Gumbel noise is working correctly!")
    else:
        print("❌ SOME TESTS FAILED - Please check the implementation")
    
    print("=" * 80)
    
    return all_pass

if __name__ == "__main__":
    success = test_gumbel_noise()
    sys.exit(0 if success else 1)