"""
Incremental CAI Calculator for Efficient Binary Search

This module implements an incremental CAI calculation approach that leverages
the logarithmic additivity of CAI to avoid recomputing the entire sequence
when only a few positions change.

Mathematical Foundation:
    CAI = (∏ w_i)^(1/L) = exp(1/L * ∑ log(w_i))
    
When positions change, we can update incrementally:
    1. Subtract log(w_old) for changed positions
    2. Add log(w_new) for new codons
    3. Recompute CAI from updated log sum

Performance: O(k) per update instead of O(L), where k << L
"""

import torch
import numpy as np
from typing import List, Tuple, Dict, Optional, Set
from dataclasses import dataclass
import time


@dataclass
class PositionChange:
    """Represents a codon change at a specific position"""
    position: int
    old_codon_idx: int
    new_codon_idx: int


class IncrementalCAICalculator:
    """
    Efficient incremental CAI calculator using logarithmic additivity.
    
    This calculator maintains a running log sum and updates it incrementally
    when codons change, avoiding full sequence recomputation.
    """
    
    def __init__(self, cai_weights: List[float], device: str = 'cuda'):
        """
        Initialize incremental CAI calculator.
        
        Args:
            cai_weights: CAI weights for each codon
            device: Computation device
        """
        self.cai_weights = cai_weights
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')
        
        # Pre-compute log weights to avoid repeated log calculations
        self.log_weights = {}
        for i, w in enumerate(cai_weights):
            if w > 0:
                self.log_weights[i] = np.log(w)
            else:
                # Use a very negative value for zero weights
                self.log_weights[i] = -20.0  # exp(-20) ≈ 2e-9
        
        # Convert to tensor for GPU acceleration if needed
        self.log_weights_tensor = torch.tensor(
            [self.log_weights[i] for i in range(len(cai_weights))],
            dtype=torch.float32,
            device=self.device
        )
        
        # State tracking
        self.current_log_sum = 0.0
        self.current_sequence = None
        self.sequence_length = 0
        self.position_to_codon = {}  # Track which codon is at each position
        
        # Performance tracking
        self.stats = {
            'full_calculations': 0,
            'incremental_updates': 0,
            'total_positions_updated': 0,
            'time_saved': 0.0
        }
    
    def initialize(self, codon_indices: List[int]) -> float:
        """
        Initialize calculator with a sequence.
        
        Args:
            codon_indices: List of codon indices for the sequence
            
        Returns:
            Initial CAI value
        """
        self.stats['full_calculations'] += 1
        
        self.sequence_length = len(codon_indices)
        self.current_sequence = list(codon_indices)
        self.position_to_codon = {i: idx for i, idx in enumerate(codon_indices)}
        
        # Calculate initial log sum
        self.current_log_sum = 0.0
        for codon_idx in codon_indices:
            if codon_idx in self.log_weights:
                self.current_log_sum += self.log_weights[codon_idx]
        
        if self.sequence_length > 0:
            cai = np.exp(self.current_log_sum / self.sequence_length)
        else:
            cai = 0.0
        
        return cai
    
    def update_single(self, position: int, new_codon_idx: int) -> float:
        """
        Update a single position incrementally.
        
        Args:
            position: Position to update
            new_codon_idx: New codon index
            
        Returns:
            Updated CAI value
        """
        if position >= self.sequence_length or position < 0:
            raise ValueError(f"Position {position} out of range [0, {self.sequence_length})")
        
        old_codon_idx = self.position_to_codon.get(position, -1)
        
        # Skip if no change
        if old_codon_idx == new_codon_idx:
            return self.get_current_cai()
        
        self.stats['incremental_updates'] += 1
        self.stats['total_positions_updated'] += 1
        
        # Remove old codon's contribution
        if old_codon_idx in self.log_weights:
            self.current_log_sum -= self.log_weights[old_codon_idx]
        
        # Add new codon's contribution
        if new_codon_idx in self.log_weights:
            self.current_log_sum += self.log_weights[new_codon_idx]
        
        # Update tracking
        self.position_to_codon[position] = new_codon_idx
        self.current_sequence[position] = new_codon_idx
        
        return self.get_current_cai()
    
    def batch_update(self, changes: List[PositionChange]) -> float:
        """
        Update multiple positions in batch.
        
        Args:
            changes: List of position changes
            
        Returns:
            Updated CAI value
        """
        if not changes:
            return self.get_current_cai()
        
        self.stats['incremental_updates'] += 1
        self.stats['total_positions_updated'] += len(changes)
        
        # Process all changes
        for change in changes:
            # Remove old contribution
            if change.old_codon_idx in self.log_weights:
                self.current_log_sum -= self.log_weights[change.old_codon_idx]
            
            # Add new contribution
            if change.new_codon_idx in self.log_weights:
                self.current_log_sum += self.log_weights[change.new_codon_idx]
            
            # Update tracking
            self.position_to_codon[change.position] = change.new_codon_idx
            self.current_sequence[change.position] = change.new_codon_idx
        
        return self.get_current_cai()
    
    def update_from_switching_events(self,
                                    cai_factor: float,
                                    switching_events: List,
                                    pi_best_codons: Dict[int, int],
                                    w_best_codons: Dict[int, int]) -> float:
        """
        Update CAI based on which switching events are active at given CAI interpolation factor.

        Args:
            cai_factor: Current CAI interpolation parameter (alpha value)
            switching_events: List of switching events with alpha* values
            pi_best_codons: Best codon in π distribution for each position
            w_best_codons: Best codon in w distribution for each position

        Returns:
            CAI value at this CAI interpolation factor
        """
        # Determine which positions have switched
        changes = []
        for event in switching_events:
            if event.alpha_star <= cai_factor:
                # This position has switched from π to w
                pos = event.position
                old_codon = pi_best_codons.get(pos, 0)
                new_codon = w_best_codons.get(pos, 0)
                
                # Only add if actually different from current
                if self.position_to_codon.get(pos) != new_codon:
                    changes.append(PositionChange(pos, old_codon, new_codon))
            else:
                # This position stays with π
                pos = event.position
                expected_codon = pi_best_codons.get(pos, 0)
                if self.position_to_codon.get(pos) != expected_codon:
                    current = self.position_to_codon.get(pos, 0)
                    changes.append(PositionChange(pos, current, expected_codon))
        
        if changes:
            return self.batch_update(changes)
        else:
            return self.get_current_cai()
    
    def get_current_cai(self) -> float:
        """
        Get current CAI value without recomputation.
        
        Returns:
            Current CAI value
        """
        if self.sequence_length > 0:
            return np.exp(self.current_log_sum / self.sequence_length)
        else:
            return 0.0
    
    def get_current_sequence(self) -> List[int]:
        """
        Get current sequence of codon indices.
        
        Returns:
            List of codon indices
        """
        return self.current_sequence.copy()
    
    def recompute_full(self, codon_indices: List[int]) -> float:
        """
        Recompute CAI from scratch (for validation).
        
        Args:
            codon_indices: Codon indices
            
        Returns:
            CAI value
        """
        log_sum = 0.0
        for idx in codon_indices:
            if idx in self.log_weights:
                log_sum += self.log_weights[idx]
        
        if len(codon_indices) > 0:
            return np.exp(log_sum / len(codon_indices))
        else:
            return 0.0
    
    def validate_consistency(self) -> bool:
        """
        Validate that incremental calculation matches full calculation.
        
        Returns:
            True if consistent, False otherwise
        """
        incremental_cai = self.get_current_cai()
        full_cai = self.recompute_full(self.current_sequence)
        
        diff = abs(incremental_cai - full_cai)
        if diff > 1e-6:
            print(f"⚠️ Inconsistency detected: incremental={incremental_cai:.6f}, full={full_cai:.6f}")
            return False
        return True
    
    def get_statistics(self) -> Dict:
        """
        Get performance statistics.
        
        Returns:
            Dictionary of statistics
        """
        avg_positions_per_update = 0
        if self.stats['incremental_updates'] > 0:
            avg_positions_per_update = (
                self.stats['total_positions_updated'] / 
                self.stats['incremental_updates']
            )
        
        return {
            'full_calculations': self.stats['full_calculations'],
            'incremental_updates': self.stats['incremental_updates'],
            'total_positions_updated': self.stats['total_positions_updated'],
            'avg_positions_per_update': avg_positions_per_update,
            'time_saved_estimate': self.stats['time_saved'],
            'efficiency_ratio': f"1:{self.sequence_length / max(avg_positions_per_update, 1):.1f}"
        }
    
    def reset_statistics(self):
        """Reset performance statistics."""
        self.stats = {
            'full_calculations': 0,
            'incremental_updates': 0,
            'total_positions_updated': 0,
            'time_saved': 0.0
        }


class CachedIncrementalCalculator(IncrementalCAICalculator):
    """
    Extended calculator with caching for binary search optimization.
    
    Caches CAI values at different α values to avoid redundant calculations
    during binary search iterations.
    """
    
    def __init__(self, cai_weights: List[float], device: str = 'cuda', cache_size: int = 100):
        """
        Initialize cached calculator.
        
        Args:
            cai_weights: CAI weights
            device: Computation device
            cache_size: Maximum cache size
        """
        super().__init__(cai_weights, device)

        self.cache_hits = 0
        self.cache_misses = 0
    
    def get_cai_at_cai_factor(self,
                         cai_factor: float,
                         switching_events: List,
                         pi_best_codons: Dict[int, int],
                         w_best_codons: Dict[int, int]) -> Tuple[float, bool]:
        """
        Get CAI value at given CAI interpolation factor with incremental update.

        This method uses incremental calculation to efficiently compute CAI
        by only updating positions that have changed.

        Args:
            cai_factor: CAI interpolation parameter (alpha value)
            switching_events: Switching events
            pi_best_codons: π-optimal codons
            w_best_codons: w-optimal codons

        Returns:
            (CAI value, always_false_for_cache_hit)
        """
        # Track cache statistics
        self.cache_misses += 1

        # Compute CAI using incremental update
        cai = self.update_from_switching_events(
            cai_factor, switching_events, pi_best_codons, w_best_codons
        )

        return cai, False
    
    def get_cache_statistics(self) -> Dict:
        """Get cache performance statistics."""
        total_accesses = self.cache_hits + self.cache_misses
        hit_rate = self.cache_hits / max(total_accesses, 1)

        return {
            'total_accesses': total_accesses,
            'cache_hits': self.cache_hits,
            'cache_misses': self.cache_misses,
            'hit_rate': f"{hit_rate:.1%}",
            'cache_efficiency': hit_rate
        }

    def clear_cache(self):
        """Clear the cache statistics."""
        # Reset cache hit/miss counters
        self.cache_hits = 0
        self.cache_misses = 0