#!/usr/bin/env python
"""
Complete performance test for SADO v4 Incremental algorithm on 3000-length sequences.
Tests both IncrementalSADO and ID3_SADO_Interface with proper target distributions.
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
from id3_sado_interface import ID3_SADO_Interface

def generate_protein_sequence(length: int = 1000) -> str:
    """Generate a random protein sequence."""
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(np.random.choice(list(amino_acids), size=length, replace=True))

def generate_target_distribution(seq_len: int, max_codons: int = 6) -> torch.Tensor:
    """
    Generate a realistic target distribution for testing.
    Simulates softmax-like probabilities.
    """
    # Generate random logits
    logits = torch.randn(seq_len, max_codons) * 2.0
    
    # Apply softmax to get probabilities
    probs = torch.softmax(logits, dim=-1)
    
    return probs

def test_incremental_sado(protein_sequence: str):
    """Test the IncrementalSADO algorithm."""
    print("\n" + "="*60)
    print("Testing IncrementalSADO Algorithm")
    print("="*60)
    
    seq_len = len(protein_sequence)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Initialize algorithm
    optimizer = IncrementalSADO(
        sequence=protein_sequence,
        target_cai=0.8
    )
    
    # Generate target distribution
    target_dist = generate_target_distribution(seq_len).to(device)
    
    print(f"Sequence length: {seq_len} amino acids ({seq_len*3} nucleotides)")
    print(f"Device: {device}")
    print(f"Target CAI: 0.8")
    
    # Run multiple iterations to test incremental nature
    num_iterations = 10
    times = []
    cai_values = []
    
    print(f"\nRunning {num_iterations} incremental iterations...")
    print("-" * 40)
    
    for i in range(num_iterations):
        start_time = time.time()
        
        # Run optimization (alpha = 0.6 for balanced optimization)
        indices, final_cai, log_prob = optimizer.optimize(target_dist, alpha=0.6)
        
        elapsed = time.time() - start_time
        times.append(elapsed)
        cai_values.append(final_cai)
        
        print(f"Iteration {i+1:2d}: CAI={final_cai:.4f}, Time={elapsed:.3f}s")
    
    # Statistics
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    print("\n" + "-" * 40)
    print(f"Average time per iteration: {avg_time:.3f}s ± {std_time:.3f}s")
    print(f"Average CAI achieved: {np.mean(cai_values):.4f}")
    print(f"Processing speed: {seq_len*3/avg_time:.0f} nt/s")
    
    return {
        'algorithm': 'IncrementalSADO',
        'avg_time': avg_time,
        'std_time': std_time,
        'avg_cai': np.mean(cai_values),
        'times': times,
        'cai_values': cai_values
    }

def test_id3_interface(protein_sequence: str):
    """Test the ID3_SADO_Interface."""
    print("\n" + "="*60)
    print("Testing ID3_SADO_Interface")
    print("="*60)
    
    seq_len = len(protein_sequence)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Initialize interface
    interface = ID3_SADO_Interface(
        amino_acid_sequence=protein_sequence,
        target_cai=0.8,
        species='ecoli_bl21de3',
        device=device
    )
    
    # Generate target distribution
    target_dist = generate_target_distribution(seq_len).to(device)
    
    print(f"Sequence length: {seq_len} amino acids ({seq_len*3} nucleotides)")
    print(f"Device: {device}")
    print(f"Target CAI: 0.8")
    
    # Run multiple iterations
    num_iterations = 10
    times = []
    cai_values = []
    
    print(f"\nRunning {num_iterations} iterations with ID3 interface...")
    print("-" * 40)
    
    for i in range(num_iterations):
        start_time = time.time()
        
        # Run optimization (gamma = 0.6 for balanced optimization)
        result = interface.optimize_discrete(target_dist, gamma=0.6)
        
        elapsed = time.time() - start_time
        times.append(elapsed)
        cai_values.append(result['final_cai'])
        
        print(f"Iteration {i+1:2d}: CAI={result['final_cai']:.4f}, Time={elapsed:.3f}s, Unique={result['is_unique']}")
    
    # Statistics
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    print("\n" + "-" * 40)
    print(f"Average time per iteration: {avg_time:.3f}s ± {std_time:.3f}s")
    print(f"Average CAI achieved: {np.mean(cai_values):.4f}")
    print(f"Processing speed: {seq_len*3/avg_time:.0f} nt/s")
    
    return {
        'algorithm': 'ID3_SADO_Interface',
        'avg_time': avg_time,
        'std_time': std_time,
        'avg_cai': np.mean(cai_values),
        'times': times,
        'cai_values': cai_values
    }

def main():
    """Main test function."""
    print("\n" + "="*80)
    print("SADO v4 Performance Test - 3000 Nucleotide Sequences")
    print("="*80)
    
    # Generate test sequence (1000 amino acids = 3000 nucleotides)
    protein_length = 1000
    test_sequence = generate_protein_sequence(protein_length)
    
    print(f"\n📊 Test Configuration:")
    print(f"  - Protein length: {protein_length} amino acids")
    print(f"  - RNA length: {protein_length * 3} nucleotides")
    print(f"  - Algorithm: SADO v4 Incremental")
    print(f"  - Target CAI: 0.8")
    
    # Test both implementations
    results = []
    
    try:
        # Test IncrementalSADO
        result1 = test_incremental_sado(test_sequence)
        results.append(result1)
    except Exception as e:
        print(f"\n❌ Error testing IncrementalSADO: {e}")
        import traceback
        traceback.print_exc()
    
    try:
        # Test ID3_SADO_Interface
        result2 = test_id3_interface(test_sequence)
        results.append(result2)
    except Exception as e:
        print(f"\n❌ Error testing ID3_SADO_Interface: {e}")
        import traceback
        traceback.print_exc()
    
    # Summary
    if results:
        print("\n" + "="*80)
        print("PERFORMANCE SUMMARY")
        print("="*80)
        
        for result in results:
            print(f"\n{result['algorithm']}:")
            print(f"  - Average time: {result['avg_time']:.3f}s ± {result['std_time']:.3f}s")
            print(f"  - Processing speed: {3000/result['avg_time']:.0f} nt/s")
            print(f"  - Average CAI: {result['avg_cai']:.4f}")
        
        # Time estimates for different lengths
        print("\n📈 Time Estimates for Different Sequence Lengths:")
        print("-" * 60)
        print("Length (nt) | IncrementalSADO | ID3_Interface")
        print("-" * 60)
        
        for nt_length in [1000, 3000, 5000, 10000, 20000]:
            line = f"{nt_length:11d} |"
            for result in results:
                time_per_nt = result['avg_time'] / 3000
                estimated = time_per_nt * nt_length
                if estimated < 60:
                    line += f" {estimated:14.1f}s |"
                else:
                    line += f" {estimated/60:12.1f}min |"
            print(line)
    
    print("\n" + "="*80)
    print("Test completed!")
    print("="*80)

if __name__ == '__main__':
    main()