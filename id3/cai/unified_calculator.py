"""
Unified CAI Calculator

Provides a unified interface for computing CAI (Codon Adaptation Index) scores
using different calculation methods (standard, differentiable, discrete, batch).
"""

import torch
import numpy as np
from typing import Union, Optional, Dict, Any, Tuple, List
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class UnifiedCAICalculator:
    """
    Unified CAI (Codon Adaptation Index) Calculator.

    Provides unified interface for computing CAI scores using different methods:
    - standard: Traditional discrete CAI from RNA sequences
    - differentiable: Soft CAI from codon probability distributions
    - discrete: CAI from codon indices
    - batch: Batch CAI computation for multiple sequences
    """
    
    def __init__(self, species: str = 'ecoli_bl21de3', device: Optional[torch.device] = None):
        """
        Initialize unified CAI calculator.

        Args:
            species: Reference species for CAI weights
            device: Computation device (CPU or CUDA)
        """
        self.species = species
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Load CAI weights
        self.wi_table = self._load_wi_table()
        self.weights_tensor = self._create_weights_tensor()

        # Initialize cache
        self._cache = {}

    def _load_wi_table(self) -> Dict[str, float]:
        """Load codon weight table (wi) from reference data."""

        from id3.optimizers.cai.utils import load_cai_weights
        wi_table, _ = load_cai_weights(self.species)
        return wi_table
    
    def _create_weights_tensor(self) -> torch.Tensor:
        """Create codon weight tensor (64-dimensional)."""
        from id3.optimizers.cai.utils import load_cai_weights
        _, weights_tensor = load_cai_weights(self.species)
        return weights_tensor.to(self.device)
    
    def compute_cai(self,
                    sequence: Union[str, torch.Tensor, np.ndarray],
                    method: str = 'standard',
                    return_log: bool = False,
                    **kwargs) -> Union[float, torch.Tensor, Dict[str, Any]]:
        """
        Compute CAI score using specified method.

        Args:
            sequence: Input sequence (RNA string, codon probabilities, or codon indices)
            method: Calculation method
                - 'standard': Traditional CAI from RNA sequence
                - 'differentiable': Soft CAI from codon probabilities
                - 'discrete': CAI from codon indices
                - 'batch': Batch computation
            return_log: Whether to return detailed log information
            **kwargs: Additional method-specific parameters

        Returns:
            CAI score (float or tensor) or detailed information (dict)
        """
        if method == 'standard':
            return self._compute_standard_cai(sequence, return_log)
        elif method == 'differentiable':
            return self._compute_differentiable_cai(sequence, **kwargs)
        elif method == 'discrete':
            return self._compute_discrete_cai(sequence, **kwargs)
        elif method == 'batch':
            return self._compute_batch_cai(sequence, **kwargs)
        else:
            raise ValueError(f"Unknown method: {method}. Available: standard, differentiable, discrete, batch")
    
    def _compute_standard_cai(self, rna_sequence: str, return_log: bool = False) -> Union[float, Dict]:
        """
        Compute standard CAI from RNA sequence.

        Args:
            rna_sequence: RNA sequence string
            return_log: Whether to return detailed information

        Returns:
            CAI score (float) or detailed log (dict)
        """
        # Normalize sequence (DNA -> RNA)
        rna_sequence = rna_sequence.replace('T', 'U')

        # Validate sequence length
        if len(rna_sequence) % 3 != 0:
            raise ValueError(f"RNA sequence length must be multiple of 3, got {len(rna_sequence)}")

        # Split into codons
        codons = [rna_sequence[i:i+3] for i in range(0, len(rna_sequence), 3)]

        # Compute CAI using geometric mean of weights
        cai_product = 1.0
        weights = []
        
        for codon in codons:
            weight = self.wi_table.get(codon, 0.0)
            if weight == 0:
                logger.warning(f"Codon {codon} has zero weight or not found in weight table")
            cai_product *= weight
            weights.append(weight)

        # Compute CAI = (product of weights) ^ (1/N)
        num_codons = len(codons)
        cai = cai_product ** (1.0 / num_codons) if num_codons > 0 else 0.0
        
        if return_log:
            return {
                'cai': cai,
                'num_codons': num_codons,
                'weights': weights,
                'codons': codons,
                'geometric_mean': cai
            }
        
        return cai
    
    def _compute_differentiable_cai(self,
                                    codon_probs: torch.Tensor,
                                    amino_acid_sequence: Optional[str] = None,
                                    temperature: float = 1.0,
                                    epsilon: float = 1e-8) -> torch.Tensor:
        """
        Compute differentiable CAI from codon probability distributions.

        Args:
            codon_probs: Codon probabilities [seq_len, num_codons] or [batch, seq_len, num_codons]
            amino_acid_sequence: Amino acid sequence for position-specific weights
            temperature: Softmax temperature (unused, kept for compatibility)
            epsilon: Small constant for numerical stability

        Returns:
            Differentiable CAI score tensor
        """
        # Ensure tensor is on correct device
        if not isinstance(codon_probs, torch.Tensor):
            codon_probs = torch.tensor(codon_probs, dtype=torch.float32, device=self.device)
        else:
            codon_probs = codon_probs.to(self.device)

        # Add batch dimension if needed
        if codon_probs.dim() == 2:
            codon_probs = codon_probs.unsqueeze(0)  # [1, seq_len, num_codons]

        batch_size, seq_len, num_codons = codon_probs.shape

        # Build position-specific weights or use global weights

        if amino_acid_sequence:
            weights = self._build_position_weights(amino_acid_sequence, num_codons)
        else:
            # Use global codon weights
            if self.weights_tensor is not None and len(self.weights_tensor) >= num_codons:
                weights = self.weights_tensor[:num_codons].unsqueeze(0).expand(seq_len, -1)
            else:
                # Fallback: uniform weights
                weights = torch.ones(seq_len, num_codons, device=self.device) * 0.5

        weights = weights.to(self.device)

        # Compute soft CAI using weighted probabilities
        # CAI = exp(mean(log(sum(p_i * w_i))))
        weighted_probs = codon_probs * weights.unsqueeze(0)  # [batch, seq_len, num_codons]
        position_cai = weighted_probs.sum(dim=-1)  # [batch, seq_len]

        # Clamp for numerical stability
        position_cai = position_cai.clamp(min=epsilon)

        # Compute geometric mean via log-space
        log_cai = position_cai.log().mean(dim=-1)  # [batch]
        cai = log_cai.exp()

        # Remove batch dimension if single sample
        if batch_size == 1:
            return cai.squeeze(0)
        
        return cai
    
    def _compute_discrete_cai(self,
                              codon_indices: torch.Tensor,
                              amino_acid_sequence: Optional[str] = None) -> float:
        """
        Compute CAI from discrete codon indices.

        Args:
            codon_indices: Codon indices tensor [seq_len] or [batch, seq_len]
            amino_acid_sequence: Amino acid sequence (unused, kept for compatibility)

        Returns:
            CAI score (float or array)
        """
        # Convert to numpy
        if isinstance(codon_indices, torch.Tensor):
            indices = codon_indices.cpu().numpy()
        else:
            indices = np.array(codon_indices)

        # Add batch dimension if needed
        if indices.ndim == 1:
            indices = indices[np.newaxis, :]

        batch_size, seq_len = indices.shape
        cai_values = []

        for b in range(batch_size):
            cai_product = 1.0
            for pos in range(seq_len):
                # Look up weight for each codon index
                codon_idx = indices[b, pos]
                if codon_idx < len(self.weights_tensor):
                    weight = self.weights_tensor[codon_idx].item()
                else:
                    weight = 0.0
                    logger.warning(f"Invalid codon index: {codon_idx}")

                cai_product *= max(weight, 1e-8)

            # Geometric mean
            cai = cai_product ** (1.0 / seq_len) if seq_len > 0 else 0.0
            cai_values.append(cai)

        # Return scalar if single sample
        if batch_size == 1:
            return cai_values[0]
        
        return np.array(cai_values)
    
    def _compute_batch_cai(self,
                          sequences: Union[List[str], torch.Tensor],
                          parallel: bool = True,
                          **kwargs) -> np.ndarray:
        """
        Compute CAI for batch of sequences.

        Args:
            sequences: List of RNA strings or tensor of probabilities/indices
            parallel: Whether to use parallel computation (unused)
            **kwargs: Additional parameters

        Returns:
            Array of CAI scores
        """
        if isinstance(sequences, list) and isinstance(sequences[0], str):
            # Batch of RNA strings
            cai_values = []
            for seq in sequences:
                cai = self._compute_standard_cai(seq)
                cai_values.append(cai)
            return np.array(cai_values)

        elif isinstance(sequences, torch.Tensor):
            # Batch of probability distributions or indices
            if sequences.dim() == 3:
                # Codon probabilities [batch, seq_len, num_codons]
                return self._compute_differentiable_cai(sequences, **kwargs).cpu().numpy()
            elif sequences.dim() == 2:
                # Codon indices [batch, seq_len]
                return self._compute_discrete_cai(sequences, **kwargs)
        
        raise ValueError(f"Unsupported input type for batch computation")
    
    def _build_position_weights(self, amino_acid_sequence: str, num_codons: int) -> torch.Tensor:
        """
        Build position-specific codon weights based on amino acid sequence.

        Args:
            amino_acid_sequence: Target amino acid sequence
            num_codons: Number of codons per position

        Returns:
            Position-specific weight tensor [seq_len, num_codons]
        """
        from id3.utils.constants import amino_acids_to_codons
        
        seq_len = len(amino_acid_sequence)
        weights = torch.zeros(seq_len, num_codons, device=self.device)
        
        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                for i, codon in enumerate(codons[:num_codons]):
                    weight = self.wi_table.get(codon, 0.0)
                    weights[pos, i] = weight
        
        return weights
    
    def get_max_achievable_cai(self, amino_acid_sequence: str) -> float:
        """
        Compute maximum achievable CAI for given amino acid sequence.

        Args:
            amino_acid_sequence: Target amino acid sequence

        Returns:
            Maximum possible CAI score
        """
        from id3.utils.constants import amino_acids_to_codons

        max_cai_product = 1.0

        for aa in amino_acid_sequence:
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                # Find best codon for this amino acid
                max_weight = 0.0
                for codon in codons:
                    weight = self.wi_table.get(codon, 0.0)
                    max_weight = max(max_weight, weight)
                max_cai_product *= max_weight
            else:
                logger.warning(f"Unknown amino acid: {aa}")

        # Geometric mean
        seq_len = len(amino_acid_sequence)
        max_cai = max_cai_product ** (1.0 / seq_len) if seq_len > 0 else 0.0
        
        return max_cai
    
    def validate_cai_constraint(self,
                                rna_sequence: str,
                                target_cai: float,
                                tolerance: float = 0.01) -> Tuple[bool, float]:
        """
        Validate whether sequence meets CAI constraint.

        Args:
            rna_sequence: RNA sequence to validate
            target_cai: Target CAI threshold
            tolerance: Allowed tolerance below target

        Returns:
            (constraint_satisfied, actual_cai)
        """
        actual_cai = self._compute_standard_cai(rna_sequence)
        satisfied = actual_cai >= target_cai - tolerance
        return satisfied, actual_cai
    
    def clear_cache(self):
        """Clear computation cache."""
        self._cache.clear()
        logger.debug("CAI calculator cache cleared")


# Convenience function
def compute_cai(sequence: Union[str, torch.Tensor],
                species: str = 'ecoli_bl21de3',
                method: str = 'auto',
                **kwargs) -> Union[float, torch.Tensor]:
    """
    Convenience function for CAI computation.

    Automatically detects input type and selects appropriate calculation method.

    Args:
        sequence: Input sequence (RNA string, probabilities, or indices)
        species: Reference species
        method: Computation method ('auto' for automatic detection)
        **kwargs: Additional method-specific parameters

    Returns:
        CAI score
    """
    # Auto-detect method based on input type
    if method == 'auto':
        if isinstance(sequence, str):
            method = 'standard'
        elif isinstance(sequence, torch.Tensor):
            if sequence.dtype in [torch.float32, torch.float64]:
                method = 'differentiable'
            else:
                method = 'discrete'
        else:
            method = 'standard'

    # Use global calculator instance
    global _global_calculator
    if _global_calculator is None or _global_calculator.species != species:
        _global_calculator = UnifiedCAICalculator(species)
    
    return _global_calculator.compute_cai(sequence, method=method, **kwargs)



_global_calculator = None