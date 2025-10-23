#!/usr/bin/env python
"""
Single iteration performance test for SADO v4 on 3000-length sequence.
Quick test to measure execution time.
"""

import sys
import time
import numpy as np
import torch
from pathlib import Path

# Add paths
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'experiments' / 'cai_enhancement_theory' / '03_SADO_algorithm'))

from algorithms.sado_v4_incremental import IncrementalSADO

def generate_protein_sequence(length: int) -> str:
    """Generate random protein sequence."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(np.random.choice(list(amino_acids), size=length))

def generate_target_distribution(seq_len: int) -> torch.Tensor:
    """Generate target distribution."""
    logits = torch.randn(seq_len, 6) * 2.0
    return torch.softmax(logits, dim=-1)

def test_single_iteration():
    """Test single iteration performance."""
    
    print("="*60)
    print("SADO v4 Single Iteration Performance Test")
    print("="*60)
    
    # Test different sequence lengths
    test_lengths = [100, 500, 1000]  # amino acids (x3 for nucleotides)
    
    for aa_length in test_lengths:
        nt_length = aa_length * 3
        print(f"\n🧬 Testing {nt_length} nt sequence ({aa_length} amino acids)...")
        
        # Generate sequence
        protein_seq = generate_protein_sequence(aa_length)
        
        # Initialize optimizer
        print("  Initializing optimizer...")
        init_start = time.time()
        optimizer = IncrementalSADO(protein_seq, target_cai=0.8)
        init_time = time.time() - init_start
        print(f"  ✓ Initialization: {init_time:.3f}s")
        
        # Generate target distribution
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        target_dist = generate_target_distribution(aa_length).to(device)
        
        # Run single optimization
        print("  Running optimization...")
        opt_start = time.time()
        
        try:
            indices, final_cai, log_prob = optimizer.optimize(target_dist, alpha=0.6)
            opt_time = time.time() - opt_start
            
            print(f"  ✓ Optimization: {opt_time:.3f}s")
            print(f"  ✓ Final CAI: {final_cai:.4f}")
            print(f"  ✓ Total time: {init_time + opt_time:.3f}s")
            print(f"  ✓ Speed: {nt_length/(init_time + opt_time):.0f} nt/s")
            
            # Estimate for 3000 nt
            if nt_length != 3000:
                estimated_3000 = (init_time + opt_time) * (3000 / nt_length)
                print(f"  📊 Estimated for 3000 nt: {estimated_3000:.1f}s")
                
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    print("\n" + "="*60)
    print("Test completed!")
    print("="*60)

if __name__ == '__main__':
    test_single_iteration()