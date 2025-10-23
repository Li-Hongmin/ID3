"""
Memory Utilities

Provides memory monitoring, optimization, and management tools for
efficient tensor operations and memory usage tracking.
"""

import torch
import gc
import psutil
import os
from typing import Optional, Dict, List, Any, Callable, Union
from contextlib import contextmanager
from functools import wraps
import warnings


class MemoryMonitor:
    """Monitor memory usage for CPU and GPU"""

    def __init__(self):
        self.gpu_available = torch.cuda.is_available()
        self.reset_stats()

    def reset_stats(self):
        """Reset statistics"""
        self.peak_memory = 0
        self.initial_memory = self.get_current_memory()
    
    def get_current_memory(self) -> Dict[str, float]:

        memory_info = {}
        

        process = psutil.Process(os.getpid())
        cpu_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_info['cpu'] = cpu_memory
        

        if self.gpu_available:
            gpu_memory = torch.cuda.memory_allocated() / 1024 / 1024  # MB
            gpu_reserved = torch.cuda.memory_reserved() / 1024 / 1024  # MB
            memory_info['gpu_allocated'] = gpu_memory
            memory_info['gpu_reserved'] = gpu_reserved
        else:
            memory_info['gpu_allocated'] = 0
            memory_info['gpu_reserved'] = 0
        
        return memory_info
    
    def get_peak_memory(self) -> Dict[str, float]:
        """Get peak memory usage"""
        current = self.get_current_memory()
        
        if self.gpu_available:
            gpu_peak = torch.cuda.max_memory_allocated() / 1024 / 1024
            current['gpu_peak'] = gpu_peak
        
        return current
    
    def print_memory_stats(self, prefix: str = ""):

        current = self.get_current_memory()
        peak = self.get_peak_memory()
        
        print(f"{prefix}Memory Usage:")
        print(f"  CPU: {current['cpu']:.1f} MB")
        
        if self.gpu_available:
            print(f"  GPU Allocated: {current['gpu_allocated']:.1f} MB")
            print(f"  GPU Reserved: {current['gpu_reserved']:.1f} MB")
            print(f"  GPU Peak: {peak.get('gpu_peak', 0):.1f} MB")


@contextmanager
def memory_efficient_context():
    """Memory-efficient context manager"""

    torch.cuda.empty_cache() if torch.cuda.is_available() else None
    gc.collect()
    
    try:
        yield
    finally:

        torch.cuda.empty_cache() if torch.cuda.is_available() else None
        gc.collect()


def memory_efficient_decorator(cleanup_after: bool = True):
    """

    
    Args:

    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):

            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            gc.collect()
            
            result = func(*args, **kwargs)
            

            if cleanup_after:
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                gc.collect()
            
            return result
        return wrapper
    return decorator


class BatchProcessor:
    """Process data in batches to manage memory usage"""

    def __init__(self, batch_size: int, max_memory_mb: float = 1000):
        """
        Initialize batch processor

        Args:
            batch_size: Batch size
            max_memory_mb: Maximum memory usage limit (MB)
        """
        self.batch_size = batch_size
        self.max_memory_mb = max_memory_mb
        self.monitor = MemoryMonitor()
    
    def process_in_batches(self,
                          data: Union[List, torch.Tensor],
                          process_func: Callable,
                          combine_func: Optional[Callable] = None) -> Any:
        """
        Process data in batches

        Args:
            data: Data to process
            process_func: Processing function
            combine_func: Result combination function (default: torch.cat)

        Returns:
            Processed results
        """
        if isinstance(data, torch.Tensor):
            total_size = data.shape[0]
            data_batches = [data[i:i+self.batch_size] for i in range(0, total_size, self.batch_size)]
        else:
            total_size = len(data)
            data_batches = [data[i:i+self.batch_size] for i in range(0, total_size, self.batch_size)]
        
        results = []
        
        for i, batch in enumerate(data_batches):

            current_memory = self.monitor.get_current_memory()
            if current_memory['cpu'] > self.max_memory_mb:
                warnings.warn(f"Memory usage ({current_memory['cpu']:.1f}MB) exceeds limit ({self.max_memory_mb}MB)")
            

            with memory_efficient_context():
                batch_result = process_func(batch)
                results.append(batch_result)
            

            if (i + 1) % 10 == 0:
                print(f"Processed {i+1}/{len(data_batches)} batches")
        

        if combine_func is None:
            if isinstance(results[0], torch.Tensor):
                return torch.cat(results, dim=0)
            else:
                return results
        else:
            return combine_func(results)


class TensorCache:
    """Tensor cache manager"""
    
    def __init__(self, max_cache_size_mb: float = 500):
        """

        
        Args:

        """
        self.max_cache_size_mb = max_cache_size_mb
        self.cache = {}
        self.cache_sizes = {}
        self.access_count = {}
    
    def get(self, key: str) -> Optional[torch.Tensor]:

        if key in self.cache:
            self.access_count[key] += 1
            return self.cache[key]
        return None
    
    def set(self, key: str, tensor: torch.Tensor) -> bool:
        """Set cached tensor"""

        tensor_size = tensor.numel() * tensor.element_size() / 1024 / 1024
        

        current_size = sum(self.cache_sizes.values())
        if current_size + tensor_size > self.max_cache_size_mb:
            self._cleanup_cache(tensor_size)
        

        self.cache[key] = tensor
        self.cache_sizes[key] = tensor_size
        self.access_count[key] = 1
        
        return True
    
    def _cleanup_cache(self, required_space: float):


        sorted_keys = sorted(self.cache.keys(), key=lambda k: self.access_count[k])
        
        freed_space = 0
        for key in sorted_keys:
            freed_space += self.cache_sizes[key]
            del self.cache[key]
            del self.cache_sizes[key]
            del self.access_count[key]
            
            if freed_space >= required_space:
                break
    
    def clear(self):
        """Clear cache"""
        self.cache.clear()
        self.cache_sizes.clear()
        self.access_count.clear()


class MemoryOptimizer:

    
    def __init__(self):
        self.monitor = MemoryMonitor()
        self.tensor_cache = TensorCache()
    
    def optimize_tensor_operations(self,
                                 tensors: List[torch.Tensor],
                                 operations: List[str] = None) -> List[torch.Tensor]:
        """
        Optimize memory usage for tensor operations

        Args:
            tensors: List of tensors
            operations: List of operations (optional)

        Returns:
            List of optimized tensors
        """
        optimized = []
        
        for tensor in tensors:

            if tensor.requires_grad:

                optimized_tensor = tensor.clone()
            else:

                optimized_tensor = tensor
            

            if tensor.dtype == torch.float64 and not tensor.requires_grad:

                optimized_tensor = optimized_tensor.to(torch.float32)
            
            optimized.append(optimized_tensor)
        
        return optimized
    
    @contextmanager
    def no_grad_inference(self):
        """Inference mode context manager"""
        with torch.no_grad():

            torch.backends.cudnn.benchmark = True
            
            try:
                yield
            finally:

                if torch.cuda.is_available():
                    torch.cuda.empty_cache()


def reduce_memory_usage(tensor: torch.Tensor, 
                       strategy: str = "auto") -> torch.Tensor:
    """

    
    Args:


        
    Returns:

    """
    if strategy == "auto" or strategy == "precision":

        if tensor.dtype == torch.float64:
            tensor = tensor.to(torch.float32)
        elif tensor.dtype == torch.int64 and tensor.max() < 2**31:
            tensor = tensor.to(torch.int32)
    
    if strategy == "auto" or strategy == "device":

        if tensor.device.type == 'cpu' and torch.cuda.is_available():
            current_memory = torch.cuda.memory_allocated()
            tensor_size = tensor.numel() * tensor.element_size()
            

            if current_memory + tensor_size < torch.cuda.get_device_properties(0).total_memory * 0.8:
                tensor = tensor.cuda()
    
    return tensor



_global_optimizer = None

def get_memory_optimizer() -> MemoryOptimizer:

    global _global_optimizer
    if _global_optimizer is None:
        _global_optimizer = MemoryOptimizer()
    return _global_optimizer


def print_memory_summary(prefix: str = ""):
    """Print memory summary"""
    monitor = MemoryMonitor()
    monitor.print_memory_stats(prefix)


def clear_all_caches():

    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    gc.collect()
    

    global _global_optimizer
    if _global_optimizer is not None:
        _global_optimizer.tensor_cache.clear()