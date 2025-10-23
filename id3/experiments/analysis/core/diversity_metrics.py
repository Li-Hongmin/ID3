"""
Diversity Metrics Module

Provides various diversity metrics for sequence analysis.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from collections import Counter
import logging

logger = logging.getLogger(__name__)


class DiversityMetrics:
    """
    Calculate various diversity metrics for sequence collections.
    """
    
    @staticmethod
    def uniqueness_rate(sequences: List[str]) -> float:
        """
        Calculate the uniqueness rate (proportion of unique sequences).
        
        Args:
            sequences: List of sequences
            
        Returns:
            Uniqueness rate (0 to 1)
        """
        if not sequences:
            return 0.0
        
        unique_sequences = set(sequences)
        return len(unique_sequences) / len(sequences)
    
    @staticmethod
    def effective_exploration_rate(sequences: List[str], theoretical_max: Optional[int] = None) -> float:
        """
        Calculate effective exploration rate.
        
        Args:
            sequences: List of sequences
            theoretical_max: Theoretical maximum number of unique sequences
            
        Returns:
            Effective exploration rate
        """
        if not sequences:
            return 0.0
        
        unique_count = len(set(sequences))
        
        if theoretical_max:
            return unique_count / theoretical_max
        else:
            # Use sequence length to estimate theoretical maximum
            if sequences[0]:
                seq_length = len(sequences[0])
                # For RNA: 4^length, but use a more realistic estimate
                theoretical_max = min(4 ** min(seq_length, 10), len(sequences) * 10)
                return unique_count / theoretical_max
            return unique_count / len(sequences)
    
    @staticmethod
    def shannon_entropy(sequences: List[str]) -> float:
        """
        Calculate Shannon entropy of sequence distribution.
        
        Higher entropy indicates more diversity.
        
        Args:
            sequences: List of sequences
            
        Returns:
            Shannon entropy
        """
        if not sequences:
            return 0.0
        
        # Count sequence frequencies
        counter = Counter(sequences)
        total = len(sequences)
        
        # Calculate entropy
        entropy = 0.0
        for count in counter.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        return entropy
    
    @staticmethod
    def simpson_diversity_index(sequences: List[str]) -> float:
        """
        Calculate Simpson's diversity index.
        
        Values range from 0 (no diversity) to 1 (maximum diversity).
        
        Args:
            sequences: List of sequences
            
        Returns:
            Simpson's diversity index
        """
        if not sequences or len(sequences) == 1:
            return 0.0
        
        counter = Counter(sequences)
        n = len(sequences)
        
        # Calculate sum of n_i * (n_i - 1)
        sum_n = sum(count * (count - 1) for count in counter.values())
        
        # Simpson's diversity index
        if n * (n - 1) == 0:
            return 0.0
        
        return 1 - (sum_n / (n * (n - 1)))
    
    @staticmethod
    def repetition_factor(sequences: List[str]) -> float:
        """
        Calculate repetition factor (inverse of uniqueness rate).
        
        Args:
            sequences: List of sequences
            
        Returns:
            Repetition factor (>=1, where 1 means all unique)
        """
        if not sequences:
            return 1.0
        
        unique_count = len(set(sequences))
        if unique_count == 0:
            return float('inf')
        
        return len(sequences) / unique_count
    
    @staticmethod
    def most_repeated_sequences(sequences: List[str], top_k: int = 10) -> List[Tuple[str, int, float]]:
        """
        Find the most repeated sequences.
        
        Args:
            sequences: List of sequences
            top_k: Number of top sequences to return
            
        Returns:
            List of (sequence, count, percentage) tuples
        """
        if not sequences:
            return []
        
        counter = Counter(sequences)
        total = len(sequences)
        
        most_common = counter.most_common(top_k)
        
        return [
            (seq, count, count / total * 100)
            for seq, count in most_common
        ]
    
    @staticmethod
    def longest_repetition_streak(sequences: List[str]) -> Tuple[str, int, int]:
        """
        Find the longest consecutive repetition of the same sequence.
        
        Args:
            sequences: List of sequences
            
        Returns:
            Tuple of (sequence, streak_length, start_position)
        """
        if not sequences:
            return ("", 0, 0)
        
        max_streak = 1
        max_seq = sequences[0] if sequences else ""
        max_pos = 0
        
        current_streak = 1
        current_seq = sequences[0] if sequences else ""
        current_pos = 0
        
        for i in range(1, len(sequences)):
            if sequences[i] == sequences[i-1]:
                current_streak += 1
            else:
                if current_streak > max_streak:
                    max_streak = current_streak
                    max_seq = current_seq
                    max_pos = current_pos
                
                current_streak = 1
                current_seq = sequences[i]
                current_pos = i
        
        # Check final streak
        if current_streak > max_streak:
            max_streak = current_streak
            max_seq = current_seq
            max_pos = current_pos
        
        return (max_seq, max_streak, max_pos)
    
    @staticmethod
    def temporal_diversity_decay(sequences: List[str], window_size: int = 100) -> List[float]:
        """
        Calculate how diversity changes over time using sliding windows.
        
        Args:
            sequences: List of sequences in temporal order
            window_size: Size of sliding window
            
        Returns:
            List of uniqueness rates for each window
        """
        if not sequences or len(sequences) < window_size:
            return []
        
        diversity_curve = []
        
        for i in range(0, len(sequences) - window_size + 1, window_size // 2):
            window = sequences[i:i + window_size]
            uniqueness = len(set(window)) / len(window)
            diversity_curve.append(uniqueness)
        
        return diversity_curve
    
    @staticmethod
    def hamming_distance_matrix(sequences: List[str], sample_size: Optional[int] = None) -> np.ndarray:
        """
        Calculate pairwise Hamming distances between sequences.
        
        Args:
            sequences: List of sequences
            sample_size: If provided, sample this many sequences for analysis
            
        Returns:
            Matrix of Hamming distances
        """
        if not sequences:
            return np.array([])
        
        # Sample if needed (for computational efficiency)
        if sample_size and len(sequences) > sample_size:
            import random
            sequences = random.sample(sequences, sample_size)
        
        n = len(sequences)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                # Calculate Hamming distance
                if len(sequences[i]) == len(sequences[j]):
                    dist = sum(c1 != c2 for c1, c2 in zip(sequences[i], sequences[j]))
                    distances[i, j] = dist
                    distances[j, i] = dist
        
        return distances
    
    @staticmethod
    def mean_pairwise_distance(sequences: List[str], sample_size: Optional[int] = 100) -> float:
        """
        Calculate mean pairwise Hamming distance.
        
        Args:
            sequences: List of sequences
            sample_size: Sample size for large sequence sets
            
        Returns:
            Mean pairwise distance
        """
        distances = DiversityMetrics.hamming_distance_matrix(sequences, sample_size)
        
        if distances.size == 0:
            return 0.0
        
        # Get upper triangle (excluding diagonal)
        n = distances.shape[0]
        if n <= 1:
            return 0.0
        
        upper_triangle = distances[np.triu_indices(n, k=1)]
        return np.mean(upper_triangle) if upper_triangle.size > 0 else 0.0
    
    @staticmethod
    def convergence_point(sequences: List[str], threshold: float = 0.1) -> Optional[int]:
        """
        Find the iteration where diversity drops below threshold.
        
        Args:
            sequences: List of sequences in temporal order
            threshold: Uniqueness rate threshold
            
        Returns:
            Iteration index where convergence occurs, or None
        """
        if not sequences:
            return None
        
        window_size = min(100, len(sequences) // 10)
        
        for i in range(window_size, len(sequences)):
            window = sequences[max(0, i - window_size):i]
            uniqueness = len(set(window)) / len(window)
            
            if uniqueness < threshold:
                return i
        
        return None
    
    @staticmethod
    def calculate_all_metrics(sequences: List[str]) -> Dict[str, any]:
        """
        Calculate all diversity metrics.
        
        Args:
            sequences: List of sequences
            
        Returns:
            Dictionary of all metrics
        """
        if not sequences:
            return {
                'uniqueness_rate': 0.0,
                'effective_exploration': 0.0,
                'shannon_entropy': 0.0,
                'simpson_diversity': 0.0,
                'repetition_factor': 1.0,
                'total_sequences': 0,
                'unique_sequences': 0,
            }
        
        unique_seqs = set(sequences)
        
        metrics = {
            'total_sequences': len(sequences),
            'unique_sequences': len(unique_seqs),
            'uniqueness_rate': DiversityMetrics.uniqueness_rate(sequences),
            'effective_exploration': DiversityMetrics.effective_exploration_rate(sequences),
            'shannon_entropy': DiversityMetrics.shannon_entropy(sequences),
            'simpson_diversity': DiversityMetrics.simpson_diversity_index(sequences),
            'repetition_factor': DiversityMetrics.repetition_factor(sequences),
        }
        
        # Add temporal analysis
        convergence = DiversityMetrics.convergence_point(sequences)
        if convergence is not None:
            metrics['convergence_iteration'] = convergence
        
        # Add streak analysis
        seq, streak, pos = DiversityMetrics.longest_repetition_streak(sequences)
        metrics['max_repetition_streak'] = streak
        metrics['max_streak_position'] = pos
        
        return metrics