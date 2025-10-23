#!/usr/bin/env python
"""
Performance test for SADO v4 incremental algorithm on 3000-length sequences.
Tests runtime, memory usage, and optimization effectiveness.
"""

import os
import sys
import time
import tracemalloc
import numpy as np
import torch
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

from id3.cai.unified_calculator import UnifiedCAICalculator
sys.path.insert(0, str(project_root / 'experiments' / 'cai_enhancement_theory' / '03_SADO_algorithm'))
from algorithms.sado_v4_incremental import IncrementalSADO

def generate_test_sequence(length: int = 3000) -> str:
    """Generate a random test sequence of specified length."""
    codons = ['AAA', 'AAC', 'AAG', 'AAU', 'ACA', 'ACC', 'ACG', 'ACU',
              'AGA', 'AGC', 'AGG', 'AGU', 'AUA', 'AUC', 'AUG', 'AUU',
              'CAA', 'CAC', 'CAG', 'CAU', 'CCA', 'CCC', 'CCG', 'CCU',
              'CGA', 'CGC', 'CGG', 'CGU', 'CUA', 'CUC', 'CUG', 'CUU',
              'GAA', 'GAC', 'GAG', 'GAU', 'GCA', 'GCC', 'GCG', 'GCU',
              'GGA', 'GGC', 'GGG', 'GGU', 'GUA', 'GUC', 'GUG', 'GUU',
              'UAC', 'UAU', 'UCA', 'UCC', 'UCG', 'UCU',
              'UGC', 'UGG', 'UGU', 'UUA', 'UUC', 'UUG', 'UUU']
    
    # Ensure length is divisible by 3
    num_codons = length // 3
    selected_codons = np.random.choice(codons, size=num_codons, replace=True)
    return ''.join(selected_codons)

def test_performance():
    """Test SADO v4 performance on 3000-length sequence."""
    
    print("=" * 80)
    print("SADO v4 Performance Test - 3000 Length Sequence")
    print("=" * 80)
    
    # Generate test sequence
    sequence_length = 3000
    test_sequence = generate_test_sequence(sequence_length)
    print(f"\n📊 Test Configuration:")
    print(f"  - Sequence length: {sequence_length} nt ({sequence_length//3} codons)")
    print(f"  - Algorithm: SADO v4 Incremental")
    print(f"  - Device: {'CUDA' if torch.cuda.is_available() else 'CPU'}")
    
    # Initialize CAI calculator
    cai_calculator = UnifiedCAICalculator(
        species='ecoli_bl21de3',
        device='cuda' if torch.cuda.is_available() else 'cpu'
    )
    
    # Algorithm parameters
    params = {
        'cai_calculator': cai_calculator,
        'target_cai': 0.8,
        'patience': 50,
        'min_improvement': 0.001,
        'window_size': 30,
        'learning_rate': 0.1,
        'momentum': 0.9,
        'batch_size': 100,
        'adaptive_lr': True,
        'lr_decay': 0.95,
        'max_iterations': 1000,
        'verbose': False
    }
    
    # Initialize algorithm
    optimizer = IncrementalSADO(**params)
    
    print("\n⏱️  Starting performance measurement...")
    print("-" * 40)
    
    # Start memory tracking
    tracemalloc.start()
    
    # Measure initialization time
    init_start = time.time()
    initial_cai = cai_calculator.compute_cai(test_sequence)
    init_time = time.time() - init_start
    
    print(f"\n1️⃣  Initialization:")
    print(f"   - Initial CAI: {initial_cai:.4f}")
    print(f"   - Init time: {init_time:.3f}s")
    
    # Measure optimization time
    opt_start = time.time()
    
    # Run optimization with progress tracking
    print(f"\n2️⃣  Running optimization...")
    optimized_sequence, metrics = optimizer.optimize(test_sequence)
    
    opt_time = time.time() - opt_start
    
    # Get memory usage
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Calculate final metrics
    final_cai = cai_calculator.compute_cai(optimized_sequence)
    
    print(f"\n3️⃣  Optimization Results:")
    print(f"   - Iterations: {metrics['iterations']}")
    print(f"   - Final CAI: {final_cai:.4f}")
    print(f"   - CAI improvement: {final_cai - initial_cai:.4f}")
    print(f"   - Target achieved: {'✅ Yes' if final_cai >= 0.8 else '❌ No'}")
    
    print(f"\n⏱️  Performance Metrics:")
    print(f"   - Total time: {opt_time:.2f}s")
    print(f"   - Time per iteration: {opt_time/metrics['iterations']*1000:.1f}ms")
    print(f"   - Time per codon: {opt_time/(sequence_length/3)*1000:.2f}ms")
    print(f"   - Processing speed: {sequence_length/opt_time:.0f} nt/s")
    
    print(f"\n💾 Memory Usage:")
    print(f"   - Current: {current / 1024 / 1024:.2f} MB")
    print(f"   - Peak: {peak / 1024 / 1024:.2f} MB")
    
    # Estimate for different sequence lengths
    print(f"\n📈 Time Estimates (based on measured performance):")
    time_per_nt = opt_time / sequence_length
    for length in [1000, 3000, 5000, 10000]:
        estimated_time = time_per_nt * length
        print(f"   - {length:5d} nt: {estimated_time:6.1f}s ({estimated_time/60:.1f} min)")
    
    # Detailed iteration analysis
    if 'history' in metrics and metrics['history']:
        print(f"\n📊 Convergence Analysis:")
        history = metrics['history']
        
        # Check early iterations
        if len(history) >= 10:
            early_improvement = history[9] - history[0]
            print(f"   - First 10 iterations improvement: {early_improvement:.4f}")
        
        # Check convergence speed
        target_reached_iter = None
        for i, cai in enumerate(history):
            if cai >= 0.8:
                target_reached_iter = i + 1
                break
        
        if target_reached_iter:
            print(f"   - Target reached at iteration: {target_reached_iter}")
            print(f"   - Time to target: {opt_time * target_reached_iter / metrics['iterations']:.2f}s")
        
        # Final convergence
        if len(history) >= 20:
            final_20_std = np.std(history[-20:])
            print(f"   - Final 20 iterations std: {final_20_std:.6f}")
            print(f"   - Converged: {'✅ Yes' if final_20_std < 0.001 else '❌ No'}")
    
    print("\n" + "=" * 80)
    print("Test completed successfully!")
    print("=" * 80)
    
    return {
        'sequence_length': sequence_length,
        'optimization_time': opt_time,
        'iterations': metrics['iterations'],
        'initial_cai': initial_cai,
        'final_cai': final_cai,
        'memory_peak_mb': peak / 1024 / 1024
    }

if __name__ == '__main__':
    # Run the test
    results = test_performance()
    
    # Summary
    print(f"\n🎯 Summary:")
    print(f"   3000 nt sequence optimized in {results['optimization_time']:.1f} seconds")
    print(f"   Performance: {3000/results['optimization_time']:.0f} nt/s")