"""
Global caching system for ID3 framework.

This module provides a unified caching mechanism to reduce redundant
computations across the framework, with memory management and persistence options.
"""

import hashlib
import json
import pickle
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Union, Callable
from functools import wraps
from collections import OrderedDict
import numpy as np
import torch

logger = logging.getLogger(__name__)


class LRUCache:
    """Least Recently Used (LRU) cache implementation."""

    def __init__(self, max_size: int = 1000):
        """
        Initialize LRU cache.

        Args:
            max_size: Maximum number of items to store
        """
        self.max_size = max_size
        self.cache = OrderedDict()
        self.hits = 0
        self.misses = 0

    def get(self, key: str) -> Optional[Any]:
        """
        Get item from cache.

        Args:
            key: Cache key

        Returns:
            Cached value or None if not found
        """
        if key in self.cache:
            # Move to end (most recently used)
            self.cache.move_to_end(key)
            self.hits += 1
            return self.cache[key]
        self.misses += 1
        return None

    def put(self, key: str, value: Any):
        """
        Put item in cache.

        Args:
            key: Cache key
            value: Value to cache
        """
        if key in self.cache:
            # Update existing and move to end
            self.cache.move_to_end(key)
        else:
            # Add new item
            if len(self.cache) >= self.max_size:
                # Remove least recently used
                self.cache.popitem(last=False)
        self.cache[key] = value

    def clear(self):
        """Clear all cached items."""
        self.cache.clear()
        self.hits = 0
        self.misses = 0

    def get_stats(self) -> Dict:
        """Get cache statistics."""
        total = self.hits + self.misses
        hit_rate = self.hits / total if total > 0 else 0
        return {
            'size': len(self.cache),
            'max_size': self.max_size,
            'hits': self.hits,
            'misses': self.misses,
            'hit_rate': hit_rate
        }


class GlobalCache:
    """
    Global caching system with multiple cache types and persistence.

    Provides different caches for different types of data:
    - CAI calculations
    - DeepRaccess predictions
    - Sequence conversions
    - Constraint computations
    """

    _instance = None
    _initialized = False

    def __new__(cls):
        """Singleton pattern to ensure single global cache."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        """Initialize global cache system."""
        if not self._initialized:
            self.caches = {
                'cai': LRUCache(max_size=10000),
                'deepraccess': LRUCache(max_size=1000),
                'sequence': LRUCache(max_size=5000),
                'constraint': LRUCache(max_size=2000),
                'general': LRUCache(max_size=5000)
            }
            self.cache_dir = Path('cache')
            self.cache_dir.mkdir(exist_ok=True)
            self._initialized = True
            logger.info("Global cache system initialized")

    def _generate_key(self, *args, **kwargs) -> str:
        """
        Generate cache key from arguments.

        Args:
            *args: Positional arguments
            **kwargs: Keyword arguments

        Returns:
            MD5 hash as cache key
        """
        # Convert arguments to string representation
        key_parts = []

        for arg in args:
            if isinstance(arg, (torch.Tensor, np.ndarray)):
                # For tensors/arrays, use shape and sample of values
                key_parts.append(f"{arg.shape}_{arg.dtype}_{hash(arg.tobytes())}")
            elif isinstance(arg, (list, tuple)):
                key_parts.append(str(arg))
            else:
                key_parts.append(str(arg))

        for k, v in sorted(kwargs.items()):
            key_parts.append(f"{k}={v}")

        key_str = "_".join(key_parts)
        return hashlib.md5(key_str.encode()).hexdigest()

    def get(self, cache_type: str, key: str) -> Optional[Any]:
        """
        Get item from specified cache.

        Args:
            cache_type: Type of cache ('cai', 'deepraccess', etc.)
            key: Cache key

        Returns:
            Cached value or None
        """
        if cache_type in self.caches:
            return self.caches[cache_type].get(key)
        return None

    def put(self, cache_type: str, key: str, value: Any):
        """
        Put item in specified cache.

        Args:
            cache_type: Type of cache
            key: Cache key
            value: Value to cache
        """
        if cache_type in self.caches:
            self.caches[cache_type].put(key, value)

    def cache_cai(self, func: Callable) -> Callable:
        """
        Decorator for caching CAI calculations.

        Args:
            func: Function to cache

        Returns:
            Wrapped function with caching
        """
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Generate cache key
            cache_key = self._generate_key(*args, **kwargs)

            # Check cache
            cached = self.get('cai', cache_key)
            if cached is not None:
                return cached

            # Compute and cache
            result = func(*args, **kwargs)
            self.put('cai', cache_key, result)
            return result

        return wrapper

    def cache_deepraccess(self, func: Callable) -> Callable:
        """
        Decorator for caching DeepRaccess predictions.

        Args:
            func: Function to cache

        Returns:
            Wrapped function with caching
        """
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Special handling for DeepRaccess inputs
            if len(args) > 1 and isinstance(args[1], torch.Tensor):
                # Create key from tensor shape and sample
                tensor = args[1]
                cache_key = f"deepraccess_{tensor.shape}_{tensor.sum().item():.6f}"
            else:
                cache_key = self._generate_key(*args, **kwargs)

            # Check cache
            cached = self.get('deepraccess', cache_key)
            if cached is not None:
                logger.debug(f"DeepRaccess cache hit: {cache_key[:8]}")
                return cached

            # Compute and cache
            result = func(*args, **kwargs)
            self.put('deepraccess', cache_key, result)
            return result

        return wrapper

    def clear_cache(self, cache_type: Optional[str] = None):
        """
        Clear specified cache or all caches.

        Args:
            cache_type: Specific cache to clear, or None for all
        """
        if cache_type:
            if cache_type in self.caches:
                self.caches[cache_type].clear()
                logger.info(f"Cleared {cache_type} cache")
        else:
            for cache in self.caches.values():
                cache.clear()
            logger.info("Cleared all caches")

    def get_stats(self) -> Dict[str, Dict]:
        """
        Get statistics for all caches.

        Returns:
            Dictionary with stats for each cache type
        """
        return {
            cache_type: cache.get_stats()
            for cache_type, cache in self.caches.items()
        }

    def save_to_disk(self, cache_type: str):
        """
        Save cache to disk for persistence.

        Args:
            cache_type: Cache type to save
        """
        if cache_type not in self.caches:
            return

        cache_file = self.cache_dir / f"{cache_type}_cache.pkl"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(self.caches[cache_type].cache, f)
            logger.info(f"Saved {cache_type} cache to disk")
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")

    def load_from_disk(self, cache_type: str):
        """
        Load cache from disk.

        Args:
            cache_type: Cache type to load
        """
        if cache_type not in self.caches:
            return

        cache_file = self.cache_dir / f"{cache_type}_cache.pkl"
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    data = pickle.load(f)
                    if isinstance(data, dict):
                        self.caches[cache_type].cache = OrderedDict(data)
                        logger.info(f"Loaded {cache_type} cache from disk")
            except Exception as e:
                logger.error(f"Failed to load cache: {e}")


# Global instance
_global_cache = GlobalCache()


def get_global_cache() -> GlobalCache:
    """Get the global cache instance."""
    return _global_cache


def cache_result(cache_type: str = 'general'):
    """
    Decorator factory for caching function results.

    Args:
        cache_type: Type of cache to use

    Returns:
        Decorator function
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            cache = get_global_cache()
            cache_key = cache._generate_key(func.__name__, *args, **kwargs)

            # Check cache
            cached = cache.get(cache_type, cache_key)
            if cached is not None:
                return cached

            # Compute and cache
            result = func(*args, **kwargs)
            cache.put(cache_type, cache_key, result)
            return result

        return wrapper
    return decorator