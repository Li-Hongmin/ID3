#!/usr/bin/env python3
"""
Test script for ID3 framework improvements.

Verifies that all improvements are working correctly:
1. Refactored modules
2. Global caching
3. Error handling
4. Performance monitoring
5. Long sequence optimization
"""

import sys
import torch
import numpy as np
import time
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))


def test_refactored_modules():
    """Test refactored experiment modules."""
    print("\n=== Testing Refactored Modules ===")

    try:
        # Test base class
        from id3.experiments.core.experiment_base import ExperimentBase
        config = {'device': 'cpu'}
        base = ExperimentBase(config)
        assert base.device == torch.device('cpu')
        print("✓ ExperimentBase initialized successfully")

        # Test constraint factory
        from id3.experiments.core.constraint_factory import ConstraintFactory
        variant_info = ConstraintFactory.get_variant_info("01")
        assert not variant_info['use_gumbel']
        assert variant_info['use_ste']
        print("✓ ConstraintFactory working correctly")

        # Test progress tracker
        from id3.experiments.core.progress_tracker import ExperimentProgressTracker
        tracker = ExperimentProgressTracker(Path("/tmp/test_tracker"))
        tracker.initialize(total_experiments=10)
        assert tracker.get_stats()['total'] == 10
        print("✓ Progress tracker functional")

        return True

    except Exception as e:
        print(f"✗ Refactored modules test failed: {e}")
        return False


def test_global_caching():
    """Test global caching system."""
    print("\n=== Testing Global Caching ===")

    try:
        from id3.utils.global_cache import get_global_cache, cache_result

        cache = get_global_cache()

        # Test basic caching
        @cache_result('general')
        def expensive_function(x):
            time.sleep(0.1)  # Simulate expensive operation
            return x * 2

        # First call - should be slow
        start = time.time()
        result1 = expensive_function(5)
        time1 = time.time() - start

        # Second call - should be fast (cached)
        start = time.time()
        result2 = expensive_function(5)
        time2 = time.time() - start

        assert result1 == result2 == 10
        assert time2 < time1 / 2  # Cached call should be much faster
        print(f"✓ Caching working: {time1:.3f}s -> {time2:.3f}s")

        # Test cache stats
        stats = cache.get_stats()
        assert stats['general']['hits'] > 0
        print(f"✓ Cache stats: {stats['general']['hit_rate']:.1%} hit rate")

        return True

    except Exception as e:
        print(f"✗ Global caching test failed: {e}")
        return False


def test_error_handling():
    """Test improved error handling."""
    print("\n=== Testing Error Handling ===")

    try:
        from id3.utils.improved_error_handling import (
            safe_execute, validate_tensor, ErrorContext,
            DataError, handle_numerical_errors
        )

        # Test safe execution
        @safe_execute(default_return=0)
        def risky_function(x):
            if x < 0:
                raise ValueError("Negative value")
            return x * 2

        assert risky_function(5) == 10
        assert risky_function(-5) == 0  # Should return default on error
        print("✓ Safe execution working")

        # Test tensor validation
        valid_tensor = torch.randn(10, 10)
        validate_tensor(valid_tensor, "test_tensor")
        print("✓ Tensor validation passed")

        # Test numerical error handling
        bad_tensor = torch.tensor([1.0, float('nan'), float('inf')])
        fixed_tensor = handle_numerical_errors(bad_tensor)
        assert not torch.isnan(fixed_tensor).any()
        assert not torch.isinf(fixed_tensor).any()
        print("✓ Numerical error handling working")

        # Test error context
        with ErrorContext("test_operation", suppress_errors=True):
            # This error will be suppressed
            raise RuntimeError("Test error")
        print("✓ Error context manager working")

        return True

    except Exception as e:
        print(f"✗ Error handling test failed: {e}")
        return False


def test_performance_monitoring():
    """Test performance monitoring system."""
    print("\n=== Testing Performance Monitoring ===")

    try:
        from id3.utils.performance_monitor import (
            get_performance_monitor, monitor_performance,
            ContextTimer, benchmark
        )

        monitor = get_performance_monitor()
        monitor.reset()  # Clear any previous metrics

        # Test function monitoring
        @monitor_performance
        def test_function(n):
            time.sleep(0.01)
            return sum(range(n))

        result = test_function(1000)
        assert result == 499500
        print("✓ Function monitoring working")

        # Test context timer
        with ContextTimer("Test operation", log_result=False) as timer:
            time.sleep(0.05)

        assert 0.04 < timer.elapsed_time < 0.06
        print(f"✓ Context timer: {timer.elapsed_time:.3f}s")

        # Test performance stats
        stats = monitor.get_statistics()
        assert stats['total_executions'] > 0
        print(f"✓ Performance stats collected: {stats['total_executions']} executions")

        # Test benchmarking
        bench_results = benchmark(lambda x: x**2, 10, iterations=50)
        assert bench_results['iterations'] == 50
        print(f"✓ Benchmarking: {bench_results['average_time']*1000:.3f}ms average")

        return True

    except Exception as e:
        print(f"✗ Performance monitoring test failed: {e}")
        return False


def test_deepraccess_sliding_window():
    """Test DeepRaccess built-in sliding window processing."""
    print("\n=== Testing DeepRaccess Sliding Window ===")

    try:
        from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper

        # Initialize DeepRaccess
        wrapper = DeepRaccessID3Wrapper()

        # Test with long sequence that requires sliding window
        seq_len = 1000  # Longer than 440 max window
        test_sequence = torch.zeros(1, seq_len, 64, device=wrapper.device)

        # Create valid one-hot encoding
        for i in range(seq_len):
            codon_idx = torch.randint(0, 64, (1,)).item()
            test_sequence[0, i, codon_idx] = 1.0

        # Process through DeepRaccess (will use sliding window internally)
        result = wrapper.forward(test_sequence)

        assert result is not None
        assert result.shape[0] == 1  # Batch dimension
        assert result.shape[1] == seq_len  # Sequence length preserved
        print(f"✓ Sliding window processing for {seq_len}bp sequence")

        # Test window_size and overlap parameters
        assert wrapper.window_size == 440
        assert wrapper.step_size == 330
        assert wrapper.overlap == 110
        print(f"✓ Window parameters: size={wrapper.window_size}, step={wrapper.step_size}, overlap={wrapper.overlap}")

        # Verify the method exists
        assert hasattr(wrapper, 'process_sliding_window')
        print("✓ process_sliding_window method exists")

        return True

    except Exception as e:
        print(f"✗ DeepRaccess sliding window test failed: {e}")
        return False


def run_all_tests():
    """Run all improvement tests."""
    print("="*60)
    print("ID3 Framework Improvements Test Suite")
    print("="*60)

    tests = [
        ("Refactored Modules", test_refactored_modules),
        ("Global Caching", test_global_caching),
        ("Error Handling", test_error_handling),
        ("Performance Monitoring", test_performance_monitoring),
        ("DeepRaccess Sliding Window", test_deepraccess_sliding_window)
    ]

    results = []
    for test_name, test_func in tests:
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"\n✗ {test_name} failed with exception: {e}")
            results.append((test_name, False))

    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)

    passed = sum(1 for _, success in results if success)
    total = len(results)

    for test_name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{test_name:30} {status}")

    print("-"*60)
    print(f"Total: {passed}/{total} tests passed ({passed/total*100:.0f}%)")

    if passed == total:
        print("\n🎉 All improvements verified successfully!")
        return 0
    else:
        print(f"\n⚠️ {total - passed} tests failed. Please review.")
        return 1


if __name__ == "__main__":
    sys.exit(run_all_tests())