"""
Performance monitoring system for ID3 framework.

Provides comprehensive performance tracking, profiling, and optimization
recommendations for the framework.
"""

import time
import psutil
import torch
import functools
import logging
from typing import Dict, List, Optional, Callable, Any
from dataclasses import dataclass, field
from collections import defaultdict
from datetime import datetime
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class PerformanceMetrics:
    """Container for performance metrics."""
    execution_time: float = 0.0
    memory_used: float = 0.0  # MB
    gpu_memory_used: float = 0.0  # MB
    cpu_percent: float = 0.0
    gpu_utilization: float = 0.0
    function_name: str = ""
    timestamp: datetime = field(default_factory=datetime.now)
    additional_metrics: Dict[str, Any] = field(default_factory=dict)


class PerformanceMonitor:
    """
    Singleton performance monitor for tracking system performance.
    """

    _instance = None
    _initialized = False

    def __new__(cls):
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize performance monitor."""
        if not self._initialized:
            self.metrics_history: List[PerformanceMetrics] = []
            self.function_stats: Dict[str, List[float]] = defaultdict(list)
            self.memory_baseline = self._get_memory_usage()
            self.gpu_available = torch.cuda.is_available()
            self._initialized = True
            self.monitoring_enabled = True
            logger.info("Performance monitor initialized")

    def _get_memory_usage(self) -> float:
        """Get current memory usage in MB."""
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024

    def _get_gpu_memory(self) -> float:
        """Get current GPU memory usage in MB."""
        if self.gpu_available:
            return torch.cuda.memory_allocated() / 1024 / 1024
        return 0.0

    def _get_gpu_utilization(self) -> float:
        """Get GPU utilization percentage."""
        if self.gpu_available:
            try:
                import pynvml
                pynvml.nvmlInit()
                handle = pynvml.nvmlDeviceGetHandleByIndex(0)
                util = pynvml.nvmlDeviceGetUtilizationRates(handle)
                return util.gpu
            except:
                return 0.0
        return 0.0

    def measure(self, func: Callable) -> Callable:
        """
        Decorator to measure function performance.

        Args:
            func: Function to measure

        Returns:
            Wrapped function with performance measurement
        """
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not self.monitoring_enabled:
                return func(*args, **kwargs)

            # Pre-execution measurements
            start_time = time.perf_counter()
            start_memory = self._get_memory_usage()
            start_gpu_memory = self._get_gpu_memory()
            cpu_percent_before = psutil.cpu_percent(interval=None)

            try:
                # Execute function
                result = func(*args, **kwargs)

                # Post-execution measurements
                end_time = time.perf_counter()
                end_memory = self._get_memory_usage()
                end_gpu_memory = self._get_gpu_memory()
                cpu_percent_after = psutil.cpu_percent(interval=None)

                # Calculate metrics
                metrics = PerformanceMetrics(
                    execution_time=end_time - start_time,
                    memory_used=end_memory - start_memory,
                    gpu_memory_used=end_gpu_memory - start_gpu_memory,
                    cpu_percent=(cpu_percent_before + cpu_percent_after) / 2,
                    gpu_utilization=self._get_gpu_utilization(),
                    function_name=func.__name__
                )

                # Store metrics
                self.metrics_history.append(metrics)
                self.function_stats[func.__name__].append(metrics.execution_time)

                # Log if execution is slow
                if metrics.execution_time > 1.0:
                    logger.warning(
                        f"Slow execution: {func.__name__} took {metrics.execution_time:.2f}s"
                    )

                return result

            except Exception as e:
                # Still track failed executions
                end_time = time.perf_counter()
                metrics = PerformanceMetrics(
                    execution_time=end_time - start_time,
                    function_name=func.__name__,
                    additional_metrics={'error': str(e)}
                )
                self.metrics_history.append(metrics)
                raise

        return wrapper

    def profile_sequence_processing(self, sequence_length: int) -> Dict[str, float]:
        """
        Profile performance for different sequence lengths.

        Args:
            sequence_length: Length of sequence being processed

        Returns:
            Performance metrics for the sequence
        """
        metrics = {
            'sequence_length': sequence_length,
            'estimated_time': self._estimate_processing_time(sequence_length),
            'estimated_memory': self._estimate_memory_usage(sequence_length),
            'optimization_suggestions': []
        }

        # Add optimization suggestions based on sequence length
        if sequence_length > 1000:
            metrics['optimization_suggestions'].append("Consider chunking for long sequences")
        if sequence_length > 500:
            metrics['optimization_suggestions'].append("Enable caching for better performance")

        return metrics

    def _estimate_processing_time(self, sequence_length: int) -> float:
        """Estimate processing time based on sequence length."""
        # Based on empirical observations (can be refined)
        base_time = 0.001  # 1ms base
        scaling_factor = 0.0001  # Linear scaling
        quadratic_factor = 0.0000001  # Quadratic scaling for very long sequences

        return base_time + scaling_factor * sequence_length + quadratic_factor * sequence_length ** 2

    def _estimate_memory_usage(self, sequence_length: int) -> float:
        """Estimate memory usage based on sequence length."""
        # Rough estimation in MB
        base_memory = 10  # Base memory usage
        per_position_memory = 0.01  # Memory per sequence position

        return base_memory + per_position_memory * sequence_length

    def get_statistics(self) -> Dict:
        """
        Get comprehensive performance statistics.

        Returns:
            Dictionary with performance statistics
        """
        if not self.metrics_history:
            return {'message': 'No metrics collected yet'}

        # Calculate statistics
        execution_times = [m.execution_time for m in self.metrics_history]
        memory_usage = [m.memory_used for m in self.metrics_history]

        stats = {
            'total_executions': len(self.metrics_history),
            'total_time': sum(execution_times),
            'average_time': np.mean(execution_times),
            'median_time': np.median(execution_times),
            'max_time': max(execution_times),
            'min_time': min(execution_times),
            'average_memory': np.mean(memory_usage),
            'max_memory': max(memory_usage),
            'function_breakdown': self._get_function_breakdown(),
            'recommendations': self._generate_recommendations()
        }

        if self.gpu_available:
            gpu_memory = [m.gpu_memory_used for m in self.metrics_history]
            stats['average_gpu_memory'] = np.mean(gpu_memory)
            stats['max_gpu_memory'] = max(gpu_memory)

        return stats

    def _get_function_breakdown(self) -> Dict[str, Dict]:
        """Get performance breakdown by function."""
        breakdown = {}
        for func_name, times in self.function_stats.items():
            if times:
                breakdown[func_name] = {
                    'calls': len(times),
                    'total_time': sum(times),
                    'average_time': np.mean(times),
                    'max_time': max(times),
                    'min_time': min(times)
                }
        return breakdown

    def _generate_recommendations(self) -> List[str]:
        """Generate performance optimization recommendations."""
        recommendations = []

        # Analyze metrics for recommendations
        if self.metrics_history:
            avg_time = np.mean([m.execution_time for m in self.metrics_history])
            max_memory = max([m.memory_used for m in self.metrics_history])

            if avg_time > 1.0:
                recommendations.append("Consider optimizing slow functions or enabling caching")

            if max_memory > 1000:  # > 1GB
                recommendations.append("High memory usage detected, consider batch processing")

            # Check for specific slow functions
            for func_name, times in self.function_stats.items():
                if times and np.mean(times) > 2.0:
                    recommendations.append(f"Function '{func_name}' is slow, consider optimization")

        return recommendations

    def reset(self):
        """Reset all collected metrics."""
        self.metrics_history.clear()
        self.function_stats.clear()
        self.memory_baseline = self._get_memory_usage()
        logger.info("Performance monitor reset")

    def enable(self):
        """Enable performance monitoring."""
        self.monitoring_enabled = True

    def disable(self):
        """Disable performance monitoring."""
        self.monitoring_enabled = False

    def save_report(self, filepath: str):
        """
        Save performance report to file.

        Args:
            filepath: Path to save the report
        """
        import json
        from pathlib import Path

        stats = self.get_statistics()
        report = {
            'timestamp': datetime.now().isoformat(),
            'statistics': stats,
            'metrics_count': len(self.metrics_history)
        }

        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2, default=str)

        logger.info(f"Performance report saved to {filepath}")


class ContextTimer:
    """Context manager for timing code blocks."""

    def __init__(self, name: str = "Operation", log_result: bool = True):
        """
        Initialize timer context.

        Args:
            name: Name of the operation being timed
            log_result: Whether to log the result
        """
        self.name = name
        self.log_result = log_result
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        """Enter timing context."""
        self.start_time = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit timing context."""
        self.end_time = time.perf_counter()
        self.elapsed = self.end_time - self.start_time

        if self.log_result:
            logger.info(f"{self.name} took {self.elapsed:.4f} seconds")

    @property
    def elapsed_time(self) -> float:
        """Get elapsed time."""
        if self.end_time is None:
            return time.perf_counter() - self.start_time
        return self.elapsed


def benchmark(func: Callable, *args, iterations: int = 100, **kwargs) -> Dict:
    """
    Benchmark a function with multiple iterations.

    Args:
        func: Function to benchmark
        *args: Function arguments
        iterations: Number of iterations
        **kwargs: Function keyword arguments

    Returns:
        Benchmark results
    """
    times = []

    # Warmup
    for _ in range(min(10, iterations // 10)):
        func(*args, **kwargs)

    # Actual benchmark
    for _ in range(iterations):
        start = time.perf_counter()
        func(*args, **kwargs)
        end = time.perf_counter()
        times.append(end - start)

    return {
        'function': func.__name__,
        'iterations': iterations,
        'total_time': sum(times),
        'average_time': np.mean(times),
        'median_time': np.median(times),
        'std_time': np.std(times),
        'min_time': min(times),
        'max_time': max(times)
    }


# Global instance
_performance_monitor = PerformanceMonitor()


def get_performance_monitor() -> PerformanceMonitor:
    """Get global performance monitor instance."""
    return _performance_monitor


# Convenience decorator
def monitor_performance(func: Callable) -> Callable:
    """Convenience decorator for performance monitoring."""
    return get_performance_monitor().measure(func)