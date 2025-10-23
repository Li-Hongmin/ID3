#!/usr/bin/env python3
"""
Optimization Tracker Module

Handles recording and managing metrics during optimization process
"""

import time
from typing import Dict, Optional, List, Any
from collections import defaultdict


class OptimizationTracker:
    """Optimization process tracker"""

    def __init__(
        self,
        record_interval: int = 10,
        verbose_interval: int = 100
    ):
        """
        Initialize tracker

        Args:
            record_interval: Data recording interval
            verbose_interval: Print output interval
        """
        self.record_interval = record_interval
        self.verbose_interval = verbose_interval
        self.data = defaultdict(list)
        self.start_time = None
        self.iteration_count = 0

    def start(self):
        """Start tracking"""
        self.start_time = time.time()
        self.data.clear()
        self.iteration_count = 0

    def record(self, iteration: int, **metrics):
        """
        Record data for one iteration

        Args:
            iteration: Iteration number
            **metrics: Arbitrary metric key-value pairs
        """
        self.iteration_count = iteration

        if iteration % self.record_interval == 0:
            self.data['iterations'].append(iteration)
            self.data['timestamps'].append(time.time() - self.start_time)

            for key, value in metrics.items():
                if value is not None:
                    # Handle tensor types
                    if hasattr(value, 'item'):
                        value = value.item()
                    self.data[key].append(value)

    def print_progress(self, iteration: int, **metrics):
        """
        Print progress information

        Args:
            iteration: Iteration number
            **metrics: Metrics to print
        """
        if iteration % self.verbose_interval == 0:
            # Build message
            parts = [f"Iter {iteration:4d}"]

            for key, value in metrics.items():
                if value is not None:
                    # Handle different value types
                    if hasattr(value, 'item'):
                        value = value.item()

                    if isinstance(value, float):
                        parts.append(f"{key}={value:.4f}")
                    elif isinstance(value, bool):
                        parts.append(f"{key}={value}")
                    else:
                        parts.append(f"{key}={value}")

            print(" | ".join(parts))

    def get_results(self) -> Dict[str, List]:
        """Get tracking results"""
        return dict(self.data)

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics"""
        elapsed_time = time.time() - self.start_time if self.start_time else 0

        summary = {
            'total_iterations': self.iteration_count,
            'elapsed_time': elapsed_time,
            'iterations_per_second': self.iteration_count / elapsed_time if elapsed_time > 0 else 0
        }

        # Add final values for each metric
        for key, values in self.data.items():
            if key not in ['iterations', 'timestamps'] and values:
                summary[f'final_{key}'] = values[-1]
                summary[f'initial_{key}'] = values[0]

                # For numeric types, calculate improvement
                if isinstance(values[0], (int, float)):
                    summary[f'{key}_improvement'] = values[0] - values[-1]

        return summary


class ProgressTracker:
    """Simple progress tracker for experiment batches"""
    
    def __init__(self, total: int):
        """
        Initialize progress tracker
        
        Args:
            total: Total number of items to track
        """
        self.total = total
        self.current = 0
        self.start_time = time.time()
    
    def update(self, description: str = ""):
        """Update progress"""
        self.current += 1
        elapsed = time.time() - self.start_time
        rate = self.current / elapsed if elapsed > 0 else 0
        
        print(f"   Progress: {self.current}/{self.total} ({self.current/self.total*100:.1f}%) - "
              f"{rate:.2f} items/s - {description}")
    
    def reset(self):
        """Reset tracker"""
        self.current = 0
        self.start_time = time.time()