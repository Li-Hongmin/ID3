"""
Incremental CAI Optimizer

Combines the stability of Binary Search with the incremental optimization efficiency of SADO:
- First iteration: Use Binary Search to obtain a high-quality initial solution
- Subsequent iterations: Use incremental optimization to perform local improvements based on previous results

This approach significantly improves performance for long sequences and multiple iterations.
"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, Tuple, Optional, Any, Set, List

from ..base import BaseCAIOptimizer
from .binary_search import BinarySearchCAIOptimizer
from .utils import load_cai_weights, compute_cai
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


class IncrementalCAIOptimizer(BaseCAIOptimizer):
    """
    Incremental CAI Optimizer

    Core innovations:
    1. First use Binary Search to establish a high-quality baseline
    2. Subsequent iterations use incremental optimization, adjusting only necessary positions
    3. Maintain historical best solutions to avoid performance degradation
    4. Intelligently determine when reinitialization is needed
    """

    # Optimization parameter constants
    MIN_K_POSITIONS = 5                    # Minimum number of positions to optimize
    MAX_K_POSITIONS = 20                   # Maximum number of positions to optimize
    K_POSITION_DIVISOR = 50                # Position count divisor
    MIN_K_PROBABILITY = 3                  # Minimum positions for probability optimization
    MAX_K_PROBABILITY = 10                 # Maximum positions for probability optimization
    K_PROBABILITY_DIVISOR = 100            # Probability position count divisor

    # Reinitialization thresholds
    CONSECUTIVE_FAILURE_THRESHOLD = 5      # Consecutive failure count threshold
    RANDOM_REINIT_PROBABILITY = 0.05       # Random reinitialization probability

    # Perturbation parameters
    MIN_PERTURB_RATIO = 0.05               # Minimum perturbation ratio
    MAX_PERTURB_RATIO = 0.10               # Maximum perturbation ratio
    MAX_PERTURB_POSITIONS = 20             # Maximum number of positions to perturb

    # Numerical stability
    NUMERICAL_EPSILON = np.finfo(float).eps  # Machine precision
    CAI_EPSILON = 1e-8                       # CAI computation epsilon
    PROBABILITY_EPSILON = 1e-10              # Probability computation epsilon

    # Scoring weights
    CAI_SCORE_WEIGHT = 0.3                  # CAI score weight
    PROBABILITY_SCORE_WEIGHT = 0.7          # Probability score weight
    
    def __init__(self,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None,
                 random_seed: Optional[int] = None,
                 enable_caching: bool = True):
        """
        Initialize incremental optimizer

        Args:
            species: Species name
            device: Computing device
            amino_acid_sequence: Amino acid sequence
            random_seed: Random seed for reproducible results
            enable_caching: Whether to enable caching optimization
        """
        super().__init__(species, device, amino_acid_sequence)

        # Setup random number generator
        self.random_seed = random_seed
        self.rng = np.random.RandomState(random_seed) if random_seed is not None else np.random

        # Initialize Binary Search optimizer
        self.bs_optimizer = BinarySearchCAIOptimizer(species, device, amino_acid_sequence)

        # Load CAI weights
        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)

        # State management
        self.last_result = None  # Last optimization result
        self.last_indices = None  # Last codon indices
        self.last_gamma = None  # Last optimal gamma
        self.iteration_count = 0  # Iteration counter
        self.consecutive_failures = 0  # Consecutive failure count

        # History tracking
        self.history_hashes: Set[str] = set()  # Avoid duplicates
        self.performance_history: List[Dict] = []  # Performance history

        # Caching optimization
        self.enable_caching = enable_caching
        if enable_caching:
            self._cai_cache: Dict[str, float] = {}  # CAI computation cache
            self._prob_cache: Dict[Tuple, float] = {}  # Probability computation cache
            self._cache_hits = 0
            self._cache_misses = 0

        # Build codon information
        if amino_acid_sequence:
            self._build_codon_info(amino_acid_sequence)
    
    def _build_codon_info(self, amino_acid_sequence: str):
        """Build codon information for incremental optimization"""
        self.codon_choices = []

        # Build standard 64 codon index mapping
        bases = ['U', 'C', 'A', 'G']
        codon_to_global_index = {}
        for i, base1 in enumerate(bases):
            for j, base2 in enumerate(bases):
                for k, base3 in enumerate(bases):
                    codon = base1 + base2 + base3
                    global_index = i * 16 + j * 4 + k
                    codon_to_global_index[codon] = global_index

        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                choices = []

                for i, codon in enumerate(codons):
                    weight = self.wi_table.get(codon, 0.0)
                    global_index = codon_to_global_index.get(codon, 0)
                    choices.append({
                        'local_index': i,
                        'global_index': global_index,
                        'codon': codon,
                        'weight': weight,
                        'aa': aa
                    })

                # Sort by weight
                choices.sort(key=lambda x: x['weight'], reverse=True)
                self.codon_choices.append(choices)
            else:
                self.codon_choices.append([])
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float,
                 amino_acid_sequence: str,
                 valid_codon_mask: Optional[torch.Tensor] = None,
                 **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Execute incremental optimization

        First call uses Binary Search, subsequent calls use incremental optimization
        """
        # Input validation
        try:
            self._validate_inputs(pi_accessibility, amino_acid_sequence, target_cai)
        except ValueError as e:
            logger.error(f"Input validation failed: {e}")
            raise

        self.iteration_count += 1

        # Check if reinitialization is needed
        should_reinit = self._should_reinitialize(pi_accessibility)

        if self.last_result is None or should_reinit:
            # First time or reinitialization needed: use Binary Search
            logger.debug(f"Iteration {self.iteration_count}: Using Binary Search initialization")

            discrete_dist, metadata = self.bs_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )

            # Save results
            self.last_result = discrete_dist
            self.last_indices = discrete_dist.argmax(dim=-1).cpu().numpy()
            self.last_gamma = metadata.get('optimal_gamma', 0.5)
            self.consecutive_failures = 0

            # Record to history
            seq_hash = self._hash_sequence(self.last_indices)
            self.history_hashes.add(seq_hash)

            metadata['method'] = 'incremental_bs'
            metadata['iteration'] = self.iteration_count

            return discrete_dist, metadata

        else:
            # Subsequent iterations: use incremental optimization
            logger.debug(f"Iteration {self.iteration_count}: Using incremental optimization")
            
            # Perform incremental optimization
            optimized_indices = self._incremental_optimize(
                self.last_indices.copy(),
                pi_accessibility,
                target_cai,
                amino_acid_sequence,
                valid_codon_mask
            )

            # Check if a new sequence was generated
            seq_hash = self._hash_sequence(optimized_indices)
            is_duplicate = seq_hash in self.history_hashes

            if is_duplicate:
                # If duplicate, try perturbing
                optimized_indices = self._perturb_sequence(
                    optimized_indices,
                    pi_accessibility,
                    target_cai,
                    valid_codon_mask
                )
                seq_hash = self._hash_sequence(optimized_indices)
                is_duplicate = seq_hash in self.history_hashes

            # Get correct codon count (handle batch dimension)
            if pi_accessibility.dim() == 3:
                num_codons = pi_accessibility.shape[2]
            else:
                num_codons = pi_accessibility.shape[1]

            # Convert to discrete distribution
            discrete_dist = self._indices_to_distribution(
                optimized_indices,
                num_codons
            )

            # Compute CAI
            final_cai = self._compute_cai_from_indices(optimized_indices)

            # Update state
            if final_cai >= target_cai:
                self.last_result = discrete_dist
                self.last_indices = optimized_indices
                self.consecutive_failures = 0
            else:
                self.consecutive_failures += 1

            # Record to history
            if not is_duplicate:
                self.history_hashes.add(seq_hash)

            # Build metadata
            metadata = {
                'method': 'incremental_inc',
                'iteration': self.iteration_count,
                'final_cai': final_cai,
                'target_cai': target_cai,
                'is_duplicate': is_duplicate,
                'incremental_changes': self._count_changes(self.last_indices, optimized_indices),
                'constraint_satisfied': final_cai >= target_cai
            }

            return discrete_dist, metadata
    
    def _incremental_optimize(self,
                              indices: np.ndarray,
                              pi_accessibility: torch.Tensor,
                              target_cai: float,
                              amino_acid_sequence: str,
                              valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        Perform incremental optimization: only optimize positions that need improvement

        Args:
            indices: Current codon indices
            pi_accessibility: RNA accessibility probability
            target_cai: Target CAI value
            amino_acid_sequence: Amino acid sequence
            valid_codon_mask: Valid codon mask

        Returns:
            Optimized codon indices
        """
        current_cai = self._compute_cai_from_indices(indices)

        # If CAI requirement is already met, try optimizing probability
        if current_cai >= target_cai:
            return self._optimize_probability(
                indices, pi_accessibility, target_cai, valid_codon_mask
            )

        # Select positions to optimize
        positions_to_optimize = self._select_positions_for_optimization(
            indices, pi_accessibility, target_cai
        )

        # Optimize selected positions
        return self._optimize_selected_positions(
            indices, positions_to_optimize, pi_accessibility, target_cai
        )

    def _select_positions_for_optimization(self,
                                          indices: np.ndarray,
                                          pi_accessibility: torch.Tensor,
                                          target_cai: float) -> np.ndarray:
        """
        Select positions to optimize

        Returns:
            Array of position indices to optimize
        """
        # Calculate improvement potential for each position
        improvements = self._calculate_improvement_potential(
            indices, pi_accessibility, target_cai
        )

        # Select top-k positions for optimization
        k = self._calculate_k_positions(len(indices))
        return np.argsort(improvements)[-k:]

    def _optimize_selected_positions(self,
                                    indices: np.ndarray,
                                    positions: np.ndarray,
                                    pi_accessibility: torch.Tensor,
                                    target_cai: float) -> np.ndarray:
        """
        Optimize selected positions

        Args:
            indices: Current codon indices
            positions: Positions to optimize
            pi_accessibility: RNA accessibility probability
            target_cai: Target CAI value

        Returns:
            Optimized codon indices
        """
        for pos in positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            # Find best codon for this position
            best_idx = self._find_best_codon_for_position(
                indices, pos, pi_accessibility, target_cai
            )
            indices[pos] = best_idx

        return indices

    def _find_best_codon_for_position(self,
                                     indices: np.ndarray,
                                     pos: int,
                                     pi_accessibility: torch.Tensor,
                                     target_cai: float) -> int:
        """
        Find best codon for a specific position

        Returns:
            Local index of the best codon
        """
        best_idx = indices[pos]
        best_score = self._evaluate_position(
            indices, pos, best_idx, pi_accessibility, target_cai
        )

        # Try all possible codons
        for choice in self.codon_choices[pos]:
            new_idx = choice['local_index']
            if new_idx == best_idx:
                continue

            score = self._evaluate_position(
                indices, pos, new_idx, pi_accessibility, target_cai
            )

            if score > best_score:
                best_idx = new_idx
                best_score = score

        return best_idx
    
    def _optimize_probability(self,
                              indices: np.ndarray,
                              pi_accessibility: torch.Tensor,
                              target_cai: float,
                              valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        Optimize probability while satisfying CAI constraints
        """
        # Handle batch dimension
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)

        # Calculate current probability score
        current_prob = self._compute_probability_score(indices, pi_accessibility)

        # Try to improve positions with lowest probability
        prob_scores = []
        for pos in range(len(indices)):
            if pos < pi_accessibility.shape[0] and indices[pos] < pi_accessibility.shape[1]:
                prob_scores.append(pi_accessibility[pos, indices[pos]].item())
            else:
                prob_scores.append(1.0)

        # Select positions with lowest probability
        k = self._calculate_k_probability_positions(len(indices))
        worst_positions = np.argsort(prob_scores)[:k]

        for pos in worst_positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            current_idx = indices[pos]
            best_idx = current_idx
            best_prob = prob_scores[pos]

            # Try codons with higher probability
            for choice in self.codon_choices[pos]:
                new_idx = choice['local_index']
                if new_idx == current_idx:
                    continue

                # Check if CAI constraint is still satisfied
                test_indices = indices.copy()
                test_indices[pos] = new_idx
                test_cai = self._compute_cai_from_indices(test_indices)

                if test_cai >= target_cai:  # Strict constraint, no tolerance
                    if new_idx < pi_accessibility.shape[1]:
                        new_prob = pi_accessibility[pos, new_idx].item()
                        if new_prob > best_prob:
                            best_idx = new_idx
                            best_prob = new_prob

            indices[pos] = best_idx

        return indices
    
    def _calculate_improvement_potential(self,
                                        indices: np.ndarray,
                                        pi_accessibility: torch.Tensor,
                                        target_cai: float) -> np.ndarray:
        """
        Calculate improvement potential for each position
        """
        potentials = np.zeros(len(indices))

        for pos in range(len(indices)):
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            current_idx = indices[pos]
            current_weight = 0.0
            max_weight = 0.0

            # Find current and maximum weights
            for choice in self.codon_choices[pos]:
                if choice['local_index'] == current_idx:
                    current_weight = choice['weight']
                max_weight = max(max_weight, choice['weight'])

            # Calculate CAI improvement potential
            cai_potential = (max_weight - current_weight) / (max_weight + self.CAI_EPSILON)

            # Calculate probability score
            prob_score = 0.0
            if pos < pi_accessibility.shape[0] and current_idx < pi_accessibility.shape[1]:
                prob_score = pi_accessibility[pos, current_idx].item()

            # Combined potential: CAI improvement potential * (1 - current probability)
            potentials[pos] = cai_potential * (1 - prob_score)

        return potentials
    
    def _evaluate_position(self,
                          indices: np.ndarray,
                          pos: int,
                          new_idx: int,
                          pi_accessibility: torch.Tensor,
                          target_cai: float) -> float:
        """
        Evaluate score for using a specific codon at a position
        """
        # Handle batch dimension
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)

        # Create test sequence
        test_indices = indices.copy()
        test_indices[pos] = new_idx

        # Compute CAI
        test_cai = self._compute_cai_from_indices(test_indices)

        # Return negative score if CAI constraint is not satisfied
        if test_cai < target_cai:
            return -1.0

        # Calculate probability score
        prob_score = 0.0
        if pos < pi_accessibility.shape[0] and new_idx < pi_accessibility.shape[1]:
            prob_score = pi_accessibility[pos, new_idx].item()

        # Combined score: CAI achievement + probability score
        cai_score = min(1.0, test_cai / target_cai)
        return cai_score * self.CAI_SCORE_WEIGHT + prob_score * self.PROBABILITY_SCORE_WEIGHT
    
    def _perturb_sequence(self,
                         indices: np.ndarray,
                         pi_accessibility: torch.Tensor,
                         target_cai: float,
                         valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        Perturb sequence to generate a new unique sequence
        """
        # Handle batch dimension
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)

        perturbed = indices.copy()
        num_positions = len(indices)

        # Randomly select positions to perturb
        num_perturb = self._calculate_perturb_positions(num_positions)
        positions = self.rng.choice(num_positions, num_perturb, replace=False)

        for pos in positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            choices = self.codon_choices[pos]
            if len(choices) <= 1:
                continue

            # Select a different codon
            current_idx = perturbed[pos]
            candidates = [c['local_index'] for c in choices if c['local_index'] != current_idx]

            if candidates:
                # Select based on probability
                probs = []
                for idx in candidates:
                    if idx < pi_accessibility.shape[1]:
                        probs.append(pi_accessibility[pos, idx].item())
                    else:
                        probs.append(0.1)

                probs = np.array(probs)

                # Improved numerical stability handling
                probs = self._normalize_probabilities(probs)

                new_idx = self.rng.choice(candidates, p=probs)
                perturbed[pos] = new_idx

        # Ensure CAI constraint is still satisfied
        perturbed_cai = self._compute_cai_from_indices(perturbed)
        if perturbed_cai < target_cai:
            # If not satisfied, return original sequence
            return indices

        return perturbed
    
    def _should_reinitialize(self, pi_accessibility: torch.Tensor) -> bool:
        """
        Determine if reinitialization is needed
        """
        # Too many consecutive failures
        if self.consecutive_failures > self.CONSECUTIVE_FAILURE_THRESHOLD:
            logger.debug("Reinitializing due to consecutive failures")
            return True

        # Input distribution changed significantly
        if self.last_result is not None:
            # Simple change detection (can be more sophisticated)
            if self.rng.random() < self.RANDOM_REINIT_PROBABILITY:  # Probabilistic reinitialization
                logger.debug("Random reinitialization")
                return True

        return False
    
    def _compute_cai_from_indices(self, indices: np.ndarray) -> float:
        """Compute CAI value with caching optimization"""
        # Generate cache key
        if self.enable_caching:
            cache_key = hashlib.md5(indices.tobytes()).hexdigest()
            if cache_key in self._cai_cache:
                self._cache_hits += 1
                return self._cai_cache[cache_key]
            self._cache_misses += 1

        cai_product = 1.0
        valid_count = 0

        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                choices = self.codon_choices[pos]
                for choice in choices:
                    if choice['local_index'] == idx:
                        weight = choice['weight']
                        if weight > 0:
                            cai_product *= weight
                            valid_count += 1
                        break

        result = cai_product ** (1.0 / valid_count) if valid_count > 0 else 0.0

        # Store in cache
        if self.enable_caching:
            # Limit cache size
            if len(self._cai_cache) > 10000:
                # Clear half of cache
                keys_to_remove = list(self._cai_cache.keys())[:5000]
                for key in keys_to_remove:
                    del self._cai_cache[key]
            self._cai_cache[cache_key] = result

        return result
    
    def _compute_probability_score(self,
                                  indices: np.ndarray,
                                  pi_accessibility: torch.Tensor) -> float:
        """Compute probability score"""
        # Handle batch dimension
        if pi_accessibility.dim() == 3:
            # If batch dimension exists, squeeze it (assuming batch_size=1)
            pi_accessibility = pi_accessibility.squeeze(0)

        log_prob = 0.0
        count = 0

        for pos, idx in enumerate(indices):
            if pos < pi_accessibility.shape[0] and idx < pi_accessibility.shape[1]:
                prob = pi_accessibility[pos, idx].item()
                if prob > 0:
                    log_prob += np.log(prob)
                    count += 1

        if count > 0:
            return np.exp(log_prob / count)
        return 0.0
    
    def _count_changes(self, old_indices: np.ndarray, new_indices: np.ndarray) -> int:
        """Count the number of changed positions"""
        if old_indices is None or new_indices is None:
            return 0
        return np.sum(old_indices != new_indices)

    def _hash_sequence(self, indices: np.ndarray) -> str:
        """Compute sequence hash"""
        return hashlib.md5(indices.tobytes()).hexdigest()

    def _indices_to_distribution(self, indices: np.ndarray, num_codons: int) -> torch.Tensor:
        """Convert indices to one-hot distribution"""
        seq_len = len(indices)
        distribution = torch.zeros(seq_len, num_codons, device=self.device)

        for pos, idx in enumerate(indices):
            if idx < num_codons:
                distribution[pos, idx] = 1.0

        return distribution

    def get_statistics(self) -> Dict[str, Any]:
        """Get statistics"""
        stats = {
            'iteration_count': self.iteration_count,
            'unique_sequences': len(self.history_hashes),
            'consecutive_failures': self.consecutive_failures,
            'last_gamma': self.last_gamma
        }

        if self.enable_caching:
            stats['cache_hits'] = self._cache_hits
            stats['cache_misses'] = self._cache_misses
            stats['cache_hit_rate'] = self._cache_hits / (self._cache_hits + self._cache_misses + 1e-10)

        return stats

    # ========== Helper Methods ==========

    def _calculate_k_positions(self, seq_length: int) -> int:
        """Calculate number of positions to optimize"""
        return min(self.MAX_K_POSITIONS,
                  max(self.MIN_K_POSITIONS, seq_length // self.K_POSITION_DIVISOR))

    def _calculate_k_probability_positions(self, seq_length: int) -> int:
        """Calculate number of positions for probability optimization"""
        return min(self.MAX_K_PROBABILITY,
                  max(self.MIN_K_PROBABILITY, seq_length // self.K_PROBABILITY_DIVISOR))

    def _calculate_perturb_positions(self, num_positions: int) -> int:
        """Calculate number of positions to perturb"""
        min_perturb = max(1, int(num_positions * self.MIN_PERTURB_RATIO))
        max_perturb = min(int(num_positions * self.MAX_PERTURB_RATIO), self.MAX_PERTURB_POSITIONS)
        return max(min_perturb, min(max_perturb, num_positions // 10))

    def _normalize_probabilities(self, probs: np.ndarray) -> np.ndarray:
        """
        Safely normalize probability distribution

        Args:
            probs: Unnormalized probability array

        Returns:
            Normalized probability distribution
        """
        # Handle empty array
        if len(probs) == 0:
            return probs

        # Handle NaN and Inf
        probs = np.nan_to_num(probs, nan=0.0, posinf=1.0, neginf=0.0)

        # Handle negative values
        probs = np.maximum(probs, 0)

        # Calculate sum
        prob_sum = probs.sum()

        # Handle zero sum or very small sum
        if prob_sum < self.NUMERICAL_EPSILON:
            # Use uniform distribution
            return np.ones(len(probs)) / len(probs)

        # Normalize
        normalized = probs / prob_sum

        # Ensure sum equals 1 (handle floating point errors)
        normalized = normalized / normalized.sum()

        return normalized

    def _validate_inputs(self,
                        pi_accessibility: torch.Tensor,
                        amino_acid_sequence: str,
                        target_cai: float) -> None:
        """
        Validate input parameters

        Raises:
            ValueError: If input is invalid
        """
        # Validate CAI target value
        if not 0 < target_cai <= 1:
            raise ValueError(f"target_cai must be in (0, 1], got {target_cai}")

        # Validate amino acid sequence
        if not amino_acid_sequence:
            raise ValueError("amino_acid_sequence cannot be empty")

        # Validate pi_accessibility tensor
        if pi_accessibility.dim() not in [2, 3]:
            raise ValueError(f"pi_accessibility must be 2D or 3D tensor, got {pi_accessibility.dim()}D")

        # Check for NaN or Inf
        if torch.isnan(pi_accessibility).any() or torch.isinf(pi_accessibility).any():
            logger.warning("pi_accessibility contains NaN or Inf values, will be handled")

    def clear_cache(self) -> None:
        """Clear cache"""
        if self.enable_caching:
            self._cai_cache.clear()
            self._prob_cache.clear()
            self._cache_hits = 0
            self._cache_misses = 0
            logger.debug("Cache cleared")