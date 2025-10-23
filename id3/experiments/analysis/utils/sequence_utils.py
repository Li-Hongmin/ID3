"""
Sequence Processing Utilities

Helper functions for sequence analysis.
"""

from typing import List, Dict, Tuple, Optional
import numpy as np


def validate_sequences(sequences: List[str]) -> bool:
    """
    Validate that all sequences are valid RNA sequences.
    
    Args:
        sequences: List of sequences
        
    Returns:
        True if all sequences are valid
    """
    if not sequences:
        return False
    
    valid_bases = set('ACGU')
    
    for seq in sequences:
        if not seq:
            return False
        if not all(base in valid_bases for base in seq.upper()):
            return False
    
    return True


def get_sequence_statistics(sequences: List[str]) -> Dict[str, any]:
    """
    Get basic statistics about sequences.
    
    Args:
        sequences: List of sequences
        
    Returns:
        Dictionary of statistics
    """
    if not sequences:
        return {
            'count': 0,
            'unique_count': 0,
            'min_length': 0,
            'max_length': 0,
            'mean_length': 0,
        }
    
    lengths = [len(seq) for seq in sequences]
    
    return {
        'count': len(sequences),
        'unique_count': len(set(sequences)),
        'min_length': min(lengths),
        'max_length': max(lengths),
        'mean_length': np.mean(lengths),
        'std_length': np.std(lengths),
    }


def find_sequence_changes(sequences: List[str]) -> List[int]:
    """
    Find iterations where sequence changes.
    
    Args:
        sequences: List of sequences in temporal order
        
    Returns:
        List of iteration indices where changes occur
    """
    if not sequences or len(sequences) < 2:
        return []
    
    changes = []
    for i in range(1, len(sequences)):
        if sequences[i] != sequences[i-1]:
            changes.append(i)
    
    return changes


def calculate_base_composition(sequence: str) -> Dict[str, float]:
    """
    Calculate base composition of a sequence.
    
    Args:
        sequence: RNA sequence
        
    Returns:
        Dictionary with base frequencies
    """
    if not sequence:
        return {'A': 0, 'C': 0, 'G': 0, 'U': 0}
    
    sequence = sequence.upper()
    total = len(sequence)
    
    composition = {
        'A': sequence.count('A') / total,
        'C': sequence.count('C') / total,
        'G': sequence.count('G') / total,
        'U': sequence.count('U') / total,
    }
    
    # GC content
    composition['GC'] = composition['G'] + composition['C']
    
    return composition


def find_repetitive_regions(sequence: str, min_repeat_length: int = 3) -> List[Tuple[int, int, str]]:
    """
    Find repetitive regions in a sequence.
    
    Args:
        sequence: RNA sequence
        min_repeat_length: Minimum length of repeat to detect
        
    Returns:
        List of (start, end, pattern) tuples
    """
    if not sequence or len(sequence) < min_repeat_length * 2:
        return []
    
    repeats = []
    
    for length in range(min_repeat_length, len(sequence) // 2 + 1):
        for start in range(len(sequence) - length * 2 + 1):
            pattern = sequence[start:start + length]
            
            # Check if pattern repeats immediately after
            if sequence[start + length:start + length * 2] == pattern:
                # Extend to find full repeat
                end = start + length * 2
                while end + length <= len(sequence) and sequence[end:end + length] == pattern:
                    end += length
                
                repeats.append((start, end, pattern))
    
    # Remove overlapping repeats (keep longest)
    repeats.sort(key=lambda x: x[1] - x[0], reverse=True)
    
    non_overlapping = []
    for repeat in repeats:
        overlap = False
        for existing in non_overlapping:
            if not (repeat[1] <= existing[0] or repeat[0] >= existing[1]):
                overlap = True
                break
        if not overlap:
            non_overlapping.append(repeat)
    
    return sorted(non_overlapping, key=lambda x: x[0])


def compress_sequence_list(sequences: List[str]) -> List[Tuple[str, int, int]]:
    """
    Compress a list of sequences by grouping consecutive repeats.
    
    Args:
        sequences: List of sequences
        
    Returns:
        List of (sequence, start_index, count) tuples
    """
    if not sequences:
        return []
    
    compressed = []
    current_seq = sequences[0]
    start_idx = 0
    count = 1
    
    for i in range(1, len(sequences)):
        if sequences[i] == current_seq:
            count += 1
        else:
            compressed.append((current_seq, start_idx, count))
            current_seq = sequences[i]
            start_idx = i
            count = 1
    
    # Add the last group
    compressed.append((current_seq, start_idx, count))
    
    return compressed


def sample_diverse_sequences(sequences: List[str], sample_size: int, 
                            strategy: str = 'uniform') -> List[str]:
    """
    Sample sequences to maintain diversity.
    
    Args:
        sequences: List of sequences
        sample_size: Number of sequences to sample
        strategy: Sampling strategy ('uniform', 'unique', 'weighted')
        
    Returns:
        Sampled sequences
    """
    if not sequences or sample_size <= 0:
        return []
    
    if len(sequences) <= sample_size:
        return sequences.copy()
    
    if strategy == 'uniform':
        # Uniform sampling across timeline
        indices = np.linspace(0, len(sequences) - 1, sample_size, dtype=int)
        return [sequences[i] for i in indices]
    
    elif strategy == 'unique':
        # Prefer unique sequences
        unique_seqs = list(set(sequences))
        if len(unique_seqs) <= sample_size:
            # Include all unique sequences and fill with repeats
            remaining = sample_size - len(unique_seqs)
            indices = np.random.choice(len(sequences), remaining, replace=False)
            return unique_seqs + [sequences[i] for i in indices]
        else:
            # Sample from unique sequences
            import random
            return random.sample(unique_seqs, sample_size)
    
    elif strategy == 'weighted':
        # Weight by inverse frequency (prefer rare sequences)
        from collections import Counter
        counts = Counter(sequences)
        weights = [1.0 / counts[seq] for seq in sequences]
        weights = np.array(weights) / sum(weights)
        
        indices = np.random.choice(len(sequences), sample_size, 
                                  replace=False, p=weights)
        return [sequences[i] for i in indices]
    
    else:
        raise ValueError(f"Unknown sampling strategy: {strategy}")