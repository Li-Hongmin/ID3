#!/usr/bin/env python
"""
Simple performance test for SADO v4 incremental algorithm on 3000-length sequences.
"""

import os
import sys
import time
import numpy as np
import torch
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / 'experiments' / 'cai_enhancement_theory' / '03_SADO_algorithm'))

from algorithms.sado_v4_incremental import IncrementalSADO

def generate_protein_sequence(length: int = 1000) -> str:
    """Generate a random protein sequence (amino acids)."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(np.random.choice(list(amino_acids), size=length, replace=True))

def test_performance():
    """Test SADO v4 performance on 1000 amino acid (3000 nt) sequence."""
    
    print("=" * 80)
    print("SADO v4 Incremental - Performance Test")
    print("=" * 80)
    
    # Generate test sequence (protein sequence, not RNA)
    protein_length = 1000  # This will be 3000 nt after codon encoding
    test_sequence = generate_protein_sequence(protein_length)
    
    print(f"\n📊 Test Configuration:")
    print(f"  - Protein length: {protein_length} amino acids")
    print(f"  - RNA length: ~{protein_length * 3} nucleotides")
    print(f"  - Target CAI: 0.8")
    print(f"  - Device: {'CUDA' if torch.cuda.is_available() else 'CPU'}")
    
    print(f"\n⏱️  Starting optimization...")
    print("-" * 40)
    
    # Initialize and run optimizer
    start_time = time.time()
    
    try:
        # Create optimizer instance
        optimizer = IncrementalSADO(
            sequence=test_sequence,
            target_cai=0.8
        )
        
        # Run optimization
        print("Running SADO optimization...")
        optimized_sequence = optimizer.optimize()
        
        elapsed_time = time.time() - start_time
        
        # Results
        print(f"\n✅ Optimization completed!")
        print(f"\n📊 Performance Metrics:")
        print(f"  - Total time: {elapsed_time:.2f} seconds")
        print(f"  - Time per amino acid: {elapsed_time/protein_length*1000:.2f} ms")
        print(f"  - Processing speed: {protein_length*3/elapsed_time:.0f} nt/s")
        
        # Estimate for different sequence lengths
        print(f"\n📈 Time Estimates (linear scaling):")
        time_per_aa = elapsed_time / protein_length
        for aa_length in [500, 1000, 2000, 5000]:
            nt_length = aa_length * 3
            estimated_time = time_per_aa * aa_length
            print(f"  - {aa_length:4d} aa ({nt_length:5d} nt): {estimated_time:6.1f}s")
            if estimated_time > 60:
                print(f"     ({estimated_time/60:.1f} minutes)")
        
        print("\n" + "=" * 80)
        print("Test completed successfully!")
        print("=" * 80)
        
        return {
            'protein_length': protein_length,
            'rna_length': protein_length * 3,
            'optimization_time': elapsed_time,
            'success': True
        }
        
    except Exception as e:
        print(f"\n❌ Error during optimization: {e}")
        import traceback
        traceback.print_exc()
        return {
            'protein_length': protein_length,
            'rna_length': protein_length * 3,
            'optimization_time': time.time() - start_time,
            'success': False,
            'error': str(e)
        }

if __name__ == '__main__':
    # Run the test
    results = test_performance()
    
    if results['success']:
        print(f"\n🎯 Summary:")
        print(f"  - {results['rna_length']} nt sequence optimized in {results['optimization_time']:.1f} seconds")
        print(f"  - Performance: {results['rna_length']/results['optimization_time']:.0f} nt/s")