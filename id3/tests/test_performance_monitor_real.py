#!/usr/bin/env python3
"""
Real-world test for performance monitoring system.

Tests performance monitoring with actual ID3 operations.
"""

import torch
import time
import json
from pathlib import Path
import sys

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from id3.utils.performance_monitor import (
    get_performance_monitor, monitor_performance, ContextTimer, benchmark
)
from id3.optimizers.cai.incremental import IncrementalCAIOptimizer
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper


def test_real_world_monitoring():
    """Test performance monitoring with real ID3 components."""

    print("=" * 60)
    print("Real-World Performance Monitoring Test")
    print("=" * 60)

    monitor = get_performance_monitor()
    monitor.reset()  # Start fresh

    # Initialize real components
    print("\n📦 Initializing ID3 components...")
    deepraccess = DeepRaccessID3Wrapper()
    cai_optimizer = IncrementalCAIOptimizer()

    # Test 1: Monitor CAI optimization
    print("\n--- Test 1: CAI Optimization Monitoring ---")

    @monitor_performance
    def optimize_cai(sequence_length=500):
        """Simulate CAI optimization."""
        # Create random sequence and required inputs
        sequence = torch.randint(0, 64, (sequence_length,))
        valid_mask = torch.ones(sequence_length, dtype=torch.bool)
        pi_accessibility = torch.randn(sequence_length)  # Mock accessibility values

        # Create mock amino acid sequence (one AA per codon)
        amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                      'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        amino_acid_sequence = ''.join([amino_acids[i % len(amino_acids)]
                                      for i in range(sequence_length)])

        # Run optimization
        optimized = cai_optimizer.optimize(
            sequence=sequence,
            valid_codon_mask=valid_mask,
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=amino_acid_sequence,
            target_cai=0.8,
            max_iterations=10
        )
        return optimized

    # Run multiple optimizations
    for i in range(3):
        print(f"  Run {i+1}: ", end="")
        with ContextTimer(f"CAI optimization {i+1}", log_result=False) as timer:
            result = optimize_cai(500 + i * 100)
        print(f"{timer.elapsed_time:.3f}s")

    # Test 2: Monitor DeepRaccess computation
    print("\n--- Test 2: DeepRaccess Monitoring ---")

    @monitor_performance
    def compute_accessibility(sequence_length=1000):
        """Compute RNA accessibility."""
        # Create one-hot encoded sequence
        sequence = torch.zeros(1, sequence_length, 64, device=deepraccess.device)
        for i in range(sequence_length):
            codon_idx = torch.randint(0, 64, (1,)).item()
            sequence[0, i, codon_idx] = 1.0

        # Compute accessibility
        result = deepraccess.forward(sequence)
        return result

    # Run multiple computations
    for i in range(3):
        print(f"  Run {i+1}: ", end="")
        with ContextTimer(f"DeepRaccess {i+1}", log_result=False) as timer:
            result = compute_accessibility(1000 + i * 500)
        print(f"{timer.elapsed_time:.3f}s")

    # Test 3: Benchmark comparison
    print("\n--- Test 3: Benchmarking ---")

    def simple_operation():
        """Simple operation for benchmarking."""
        return torch.sum(torch.randn(100, 100))

    bench_results = benchmark(simple_operation, iterations=100)
    print(f"  Function: {bench_results['function']}")
    print(f"  Iterations: {bench_results['iterations']}")
    print(f"  Average time: {bench_results['average_time']*1000:.3f}ms")
    print(f"  Std deviation: {bench_results['std_time']*1000:.3f}ms")

    # Test 4: Profile sequence processing
    print("\n--- Test 4: Sequence Profiling ---")
    for seq_len in [500, 1000, 2000]:
        profile = monitor.profile_sequence_processing(seq_len)
        print(f"\n  Sequence length: {seq_len}bp")
        print(f"  Estimated time: {profile['estimated_time']*1000:.2f}ms")
        print(f"  Estimated memory: {profile['estimated_memory']:.1f}MB")
        if profile['optimization_suggestions']:
            print("  Suggestions:")
            for suggestion in profile['optimization_suggestions']:
                print(f"    • {suggestion}")

    # Get comprehensive statistics
    print("\n" + "=" * 60)
    print("Performance Statistics")
    print("=" * 60)

    stats = monitor.get_statistics()

    print(f"\n📊 Overall Metrics:")
    print(f"  Total executions: {stats['total_executions']}")
    print(f"  Total time: {stats['total_time']:.3f}s")
    print(f"  Average time: {stats['average_time']:.3f}s")
    print(f"  Max time: {stats['max_time']:.3f}s")

    print(f"\n📈 Function Breakdown:")
    for func_name, func_stats in stats['function_breakdown'].items():
        print(f"  {func_name}:")
        print(f"    Calls: {func_stats['calls']}")
        print(f"    Total: {func_stats['total_time']:.3f}s")
        print(f"    Average: {func_stats['average_time']:.3f}s")

    if stats['recommendations']:
        print(f"\n💡 Recommendations:")
        for rec in stats['recommendations']:
            print(f"  • {rec}")

    # Save report
    report_path = Path("/tmp/performance_report.json")
    monitor.save_report(str(report_path))
    print(f"\n📝 Report saved to: {report_path}")

    # Verify report was saved
    if report_path.exists():
        with open(report_path) as f:
            report = json.load(f)
        print(f"✅ Report contains {report['metrics_count']} metrics")
    else:
        print("❌ Failed to save report")

    return stats


def test_memory_tracking():
    """Test memory usage tracking."""

    print("\n" + "=" * 60)
    print("Memory Tracking Test")
    print("=" * 60)

    monitor = get_performance_monitor()

    print("\n🧠 Memory Usage:")

    # Baseline memory
    baseline = monitor._get_memory_usage()
    print(f"  Baseline: {baseline:.2f} MB")

    # Allocate some memory
    large_tensor = torch.randn(10000, 10000)
    current = monitor._get_memory_usage()
    print(f"  After allocation: {current:.2f} MB")
    print(f"  Increase: {current - baseline:.2f} MB")

    # GPU memory if available
    if torch.cuda.is_available():
        gpu_memory = monitor._get_gpu_memory()
        print(f"\n🎮 GPU Memory:")
        print(f"  Allocated: {gpu_memory:.2f} MB")

        # Allocate GPU tensor
        gpu_tensor = torch.randn(5000, 5000, device='cuda')
        gpu_memory_after = monitor._get_gpu_memory()
        print(f"  After GPU allocation: {gpu_memory_after:.2f} MB")
        print(f"  Increase: {gpu_memory_after - gpu_memory:.2f} MB")

    # Clean up
    del large_tensor
    if torch.cuda.is_available():
        del gpu_tensor
        torch.cuda.empty_cache()


if __name__ == "__main__":
    print("🚀 Starting Performance Monitoring System Tests\n")

    # Run real-world monitoring tests
    stats = test_real_world_monitoring()

    # Run memory tracking tests
    test_memory_tracking()

    print("\n✅ All performance monitoring tests completed!")

    # Final summary
    if stats['total_executions'] > 0:
        print(f"\n📊 Final Summary:")
        print(f"   Monitored {stats['total_executions']} operations")
        print(f"   Total runtime: {stats['total_time']:.2f}s")
        print(f"   System is working correctly!")
    else:
        print("\n⚠️ No operations were monitored")