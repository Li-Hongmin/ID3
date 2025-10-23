"""
Discrete Binary Search for CAI Optimization

This module implements a discretized version of binary search that operates
on switching events rather than continuous α values, improving efficiency
from O(log(1/ε)) to O(log n) where n is the sequence length.

Key Innovation:
Instead of searching in continuous α ∈ [0,1], we pre-compute all switching
points where positions change from π to w distributions, then search in this
discrete space.
"""

import torch
import numpy as np
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from .incremental_cai_calculator import IncrementalCAICalculator, CachedIncrementalCalculator


@dataclass
class SwitchingEvent:
    """Represents a switching event where a position changes from π to w"""
    position: int  # Position in the sequence
    alpha_star: float  # α value where switching occurs
    pi_best: int  # Best codon index in π distribution
    w_best: int  # Best codon index in w distribution
    cai_gain: float  # CAI improvement from this switch


class DiscreteCAISearcher:
    """
    Discrete binary search for CAI optimization.
    
    This approach pre-computes all possible switching points and searches
    in discrete space, significantly improving efficiency over continuous search.
    """
    
    def __init__(self, cai_weights: List[float], device: str = 'cuda', use_incremental: bool = True):
        """
        Initialize discrete CAI searcher.
        
        Args:
            cai_weights: CAI weights for each codon
            device: Computation device
            use_incremental: Whether to use incremental CAI calculation
        """
        self.cai_weights = cai_weights
        self.device = torch.device(device if torch.cuda.is_available() else 'cpu')
        self.use_incremental = use_incremental
        
        # Initialize incremental calculator if enabled
        if self.use_incremental:
            self.incremental_calc = CachedIncrementalCalculator(cai_weights, device)
    
    def compute_switching_events(self,
                                 pi_probs: torch.Tensor,
                                 valid_codon_mask: torch.Tensor,
                                 codon_indices: torch.Tensor,
                                 amino_sequence: str) -> List[SwitchingEvent]:
        """
        Compute all switching events for the sequence.
        
        For each position, calculate the α value where it switches from
        choosing the π-optimal codon to the w-optimal codon.
        
        Args:
            pi_probs: π probability distribution [positions, max_codons]
            valid_codon_mask: Valid codon mask [positions, max_codons]
            codon_indices: Global codon indices [positions, max_codons]
            amino_sequence: Target amino acid sequence
            
        Returns:
            List of switching events sorted by α value
        """
        events = []
        num_positions = pi_probs.shape[0]
        
        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos]
            if not valid_mask.any():
                continue
            
            # Get π-optimal codon (highest probability in accessibility distribution)
            pi_pos_probs = pi_probs[pos, valid_mask]
            pi_best_local = torch.argmax(pi_pos_probs).item()
            
            # Get w-optimal codon (highest CAI weight)
            w_best_local = -1
            max_cai_weight = -1
            valid_slots = torch.where(valid_mask)[0]
            
            for local_idx, slot in enumerate(valid_slots.cpu().numpy()):
                global_codon_idx = codon_indices[pos, slot].item()
                if global_codon_idx < len(self.cai_weights):
                    if self.cai_weights[global_codon_idx] > max_cai_weight:
                        max_cai_weight = self.cai_weights[global_codon_idx]
                        w_best_local = local_idx
            
            # Skip if no valid w codon or if π and w choose the same codon
            if w_best_local == -1 or pi_best_local == w_best_local:
                continue
            
            # Calculate switching point α*
            # At α*, the interpolated distribution switches its argmax
            # α*w_best + (1-α*)π_best = α*w_pi + (1-α*)π_pi
            # Solving for when w_best becomes larger than π_best
            
            p_pi = pi_pos_probs[pi_best_local].item()  # π prob of π-best codon
            p_w = pi_pos_probs[w_best_local].item() if w_best_local < len(pi_pos_probs) else 0  # π prob of w-best codon
            
            # Get CAI weights for both codons
            pi_best_global = codon_indices[pos, valid_slots[pi_best_local]].item()
            w_best_global = codon_indices[pos, valid_slots[w_best_local]].item()
            
            cai_pi = self.cai_weights[pi_best_global] if pi_best_global < len(self.cai_weights) else 0.1
            cai_w = self.cai_weights[w_best_global] if w_best_global < len(self.cai_weights) else 0.1
            
            # Calculate switching point
            # The interpolated probability of codon c is: α*w_c + (1-α)*π_c
            # For π-best: α*cai_pi_norm + (1-α)*p_pi
            # For w-best: α*cai_w_norm + (1-α)*p_w
            # Switch occurs when these are equal
            
            # Simplified calculation for switching point
            if abs(p_pi - p_w) > 1e-10:
                # Calculate normalized CAI probabilities (after softmax)
                # This is approximate since we need full w distribution
                cai_diff = cai_w - cai_pi
                if cai_diff > 0:
                    # Switching point where w-best overtakes π-best
                    alpha_star = (p_pi - p_w) / (p_pi - p_w + cai_diff)
                    alpha_star = max(0.0, min(1.0, alpha_star))
                else:
                    # w never overtakes π
                    continue
            else:
                # Special case: equal π probabilities
                alpha_star = 0.5
            
            event = SwitchingEvent(
                position=pos,
                alpha_star=alpha_star,
                pi_best=pi_best_local,
                w_best=w_best_local,
                cai_gain=cai_w - cai_pi
            )
            events.append(event)
        
        # Sort events by switching point
        events.sort(key=lambda e: e.alpha_star)
        
        # Handle simultaneous switches by adding small perturbations
        events = self._break_ties(events)
        
        return events
    
    def _break_ties(self, events: List[SwitchingEvent]) -> List[SwitchingEvent]:
        """
        Break ties for events with identical switching points.
        
        When multiple positions have the same α*, we assign them slightly
        different values to maintain monotonicity.
        
        Args:
            events: List of switching events sorted by α*
            
        Returns:
            Events with unique switching points
        """
        if len(events) <= 1:
            return events
        
        epsilon = 1e-10
        modified_events = []
        
        i = 0
        while i < len(events):
            # Find group with same α*
            current_alpha = events[i].alpha_star
            group_start = i
            
            while i < len(events) and abs(events[i].alpha_star - current_alpha) < 1e-8:
                i += 1
            
            group_size = i - group_start
            
            if group_size > 1:
                # Randomly shuffle the group to break ties
                group = events[group_start:i]
                np.random.shuffle(group)
                
                # Assign incrementally different α values
                for j, event in enumerate(group):
                    event.alpha_star = current_alpha + j * epsilon
                    modified_events.append(event)
            else:
                modified_events.append(events[group_start])
        
        return modified_events
    
    def discrete_binary_search(self,
                               pi_probs: torch.Tensor,
                               valid_codon_mask: torch.Tensor,
                               codon_indices: torch.Tensor,
                               amino_sequence: str,
                               target_cai: float = 0.8,
                               verbose: bool = False) -> Tuple[float, Dict]:
        """
        Perform binary search in discrete switching event space.
        
        Instead of searching α ∈ [0,1], we search among pre-computed
        switching events, making the search O(log n) instead of O(log(1/ε)).
        
        Args:
            pi_probs: π probability distribution
            valid_codon_mask: Valid codon mask
            codon_indices: Global codon indices
            amino_sequence: Target amino acid sequence
            target_cai: Target CAI value
            verbose: Whether to print debug info
            
        Returns:
            optimal_alpha: Optimal α value
            metadata: Search metadata and statistics
        """
        # Compute all switching events
        events = self.compute_switching_events(
            pi_probs, valid_codon_mask, codon_indices, amino_sequence
        )
        
        if verbose:
            print(f"🎯 Discrete Binary Search {'with Incremental CAI' if self.use_incremental else ''}")
            print(f"   Found {len(events)} switching events")
            if len(events) > 0:
                print(f"   α range: [{events[0].alpha_star:.6f}, {events[-1].alpha_star:.6f}]")
        
        # Special case: no switching events
        if len(events) == 0:
            if verbose:
                print("   No switching events - all positions have same optimal codon")
            return 0.0, {'num_events': 0, 'iterations': 0}
        
        # Pre-compute best codons for all positions if using incremental
        pi_best_codons = {}
        w_best_codons = {}
        initial_sequence = []
        
        num_positions = pi_probs.shape[0]
        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos]
            if not valid_mask.any():
                continue
            
            # Get π-optimal codon
            valid_probs = pi_probs[pos, valid_mask]
            pi_best_local = torch.argmax(valid_probs).item()

            pi_best_global = codon_indices[pos, valid_slots[pi_best_local]].item()
            pi_best_codons[pos] = pi_best_global
            initial_sequence.append(pi_best_global)  # Start with π-optimal
            
            # Get w-optimal codon
            w_best_local = -1
            max_weight = -1
            for local_idx, slot in enumerate(valid_slots.cpu().numpy()):
                global_idx = codon_indices[pos, slot].item()
                if global_idx < len(self.cai_weights):
                    if self.cai_weights[global_idx] > max_weight:
                        max_weight = self.cai_weights[global_idx]
                        w_best_local = local_idx
            
            if w_best_local >= 0:
                w_best_global = codon_indices[pos, valid_slots[w_best_local]].item()
                w_best_codons[pos] = w_best_global
            else:
                w_best_codons[pos] = pi_best_global  # Fallback
        
        # Initialize incremental calculator if enabled
        if self.use_incremental:
            self.incremental_calc.initialize(initial_sequence)
        
        # Add boundary events for α=0 and α=1
        boundary_events = [
            SwitchingEvent(position=-1, alpha_star=0.0, pi_best=-1, w_best=-1, cai_gain=0),
            *events,
            SwitchingEvent(position=-1, alpha_star=1.0, pi_best=-1, w_best=-1, cai_gain=0)
        ]
        
        # Binary search in discrete space
        left = 0
        right = len(boundary_events) - 1
        best_idx = right  # Default to maximum CAI
        iterations = 0
        
        while left <= right:
            mid = (left + right) // 2
            cai_factor_mid = boundary_events[mid].alpha_star
            

            if self.use_incremental:
                # Use incremental calculation with caching
                cai_mid, cache_hit = self.incremental_calc.get_cai_at_cai_factor(
                    cai_factor_mid, events, pi_best_codons, w_best_codons
                )
            else:
                # Use original full computation
                cai_mid = self._compute_cai_at_cai_factor(
                    cai_factor_mid, pi_probs, valid_codon_mask, 

                )
            
            iterations += 1
            
            if verbose:
                status = "✅" if cai_mid >= target_cai else "❌"

            
            if cai_mid >= target_cai:
                # Requirement satisfied, try smaller α
                best_idx = mid
                right = mid - 1
            else:
                # Need higher CAI
                left = mid + 1
        
        optimal_alpha = boundary_events[best_idx].alpha_star
        
        if verbose:
            print(f"✅ Found optimal: α={optimal_alpha:.6f} at event {best_idx}")
        
        metadata = {
            'num_events': len(events),
            'iterations': iterations,
            'event_index': best_idx,
            'search_space_reduction': f"{iterations}/{len(boundary_events)} = {iterations/len(boundary_events):.1%}"
        }
        
        # Add incremental calculator statistics if used
        if self.use_incremental:
            metadata['incremental_stats'] = self.incremental_calc.get_statistics()
            metadata['cache_stats'] = self.incremental_calc.get_cache_statistics()
        
        return optimal_alpha, metadata
    
    def _compute_cai_at_cai_factor(self,
                             cai_factor: float,
                             pi_probs: torch.Tensor,
                             valid_codon_mask: torch.Tensor,
                             codon_indices: torch.Tensor,
                             active_events: List[SwitchingEvent]) -> float:
        """
        Compute CAI value at given CAI interpolation factor by applying active switching events.

        Args:
            cai_factor: CAI interpolation parameter (alpha value)
            pi_probs: π probability distribution
            valid_codon_mask: Valid codon mask
            codon_indices: Global codon indices
            active_events: Events that have switched (alpha* < cai_factor)

        Returns:
            CAI value at this α
        """
        # Build switched positions set for efficiency
        switched_positions = {event.position for event in active_events}
        
        # Compute discrete sequence
        from id3.cai.validator import compute_cai_from_sequence
        discrete_seq = self._get_discrete_sequence(
            cai_factor, pi_probs, valid_codon_mask, codon_indices, switched_positions
        )
        
        return compute_cai_from_sequence(discrete_seq)
    
    def _get_discrete_sequence(self,
                               alpha: float,
                               pi_probs: torch.Tensor,
                               valid_codon_mask: torch.Tensor,
                               codon_indices: torch.Tensor,
                               switched_positions: set) -> str:
        """
        Get discrete sequence at given α with specified switched positions.
        
        Args:
            alpha: Interpolation parameter
            pi_probs: π probability distribution
            valid_codon_mask: Valid codon mask
            codon_indices: Global codon indices
            switched_positions: Set of positions that have switched to w
            
        Returns:
            Discrete RNA sequence
        """
        from id3.utils.constants import codon_to_rna_matrix
        
        num_positions = pi_probs.shape[0]
        selected_codons = []
        
        for pos in range(num_positions):
            valid_mask = valid_codon_mask[pos]
            if not valid_mask.any():
                continue
            
            if pos in switched_positions:
                # Use w-optimal codon (highest CAI)
                best_idx = -1
                max_weight = -1
                valid_slots = torch.where(valid_mask)[0]
                
                for local_idx, slot in enumerate(valid_slots.cpu().numpy()):
                    global_idx = codon_indices[pos, slot].item()
                    if global_idx < len(self.cai_weights):
                        if self.cai_weights[global_idx] > max_weight:
                            max_weight = self.cai_weights[global_idx]
                            best_idx = local_idx
                
                if best_idx >= 0:
                    actual_slot = valid_slots[best_idx].item()
                    global_codon_idx = codon_indices[pos, actual_slot].item()
                else:
                    # Fallback to first valid
                    actual_slot = valid_slots[0].item()
                    global_codon_idx = codon_indices[pos, actual_slot].item()
            else:
                # Use π-optimal codon (highest probability)
                valid_probs = pi_probs[pos, valid_mask]
                best_idx = torch.argmax(valid_probs).item()
                valid_slots = torch.where(valid_mask)[0]
                actual_slot = valid_slots[best_idx].item()
                global_codon_idx = codon_indices[pos, actual_slot].item()
            
            selected_codons.append(global_codon_idx)
        
        # Convert to RNA sequence
        rna_parts = []
        for codon_idx in selected_codons:
            codon_rna = codon_to_rna_matrix[codon_idx]  # [3, 4] one-hot
            for pos in range(3):
                nt_idx = torch.argmax(codon_rna[pos]).item()
                nucleotides = ['A', 'C', 'G', 'U']
                rna_parts.append(nucleotides[nt_idx])
        
        return ''.join(rna_parts)
