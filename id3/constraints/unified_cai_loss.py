"""
Unified CAI Loss Implementation

This module implements the ECAI (Expected CAI) formula from the paper for
unified optimization with accessibility loss. This enables gradient-based
optimization of both accessibility and CAI simultaneously.

Mathematical Framework from Paper:
    ECAI(π) = (∏_{j=1}^M ∑_{c∈C(y_j)} π_{j,c} · w_c)^{1/M}
    L_CAI(Θ) = (ECAI(T_{ID3-Xαβ}(Θ)) - CAI_target)²
    L_unified = L_Access + λ_CAI * L_CAI

Where:
    - π: codon probability distributions from constraint mechanisms
    - w_c: CAI relative adaptiveness weights for codon c
    - M: number of amino acid positions
    - T_{ID3-Xαβ}: ID3 transformation (varies by constraint type)
"""

# Copyright (c) 2025 University of Tokyo
# Licensed under CC BY-NC-SA 4.0
# For commercial use, contact: lihongmin@edu.k.u-tokyo.ac.jp

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, Tuple, Optional, Union
import json
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class ECAIComputation:
    """
    Expected CAI (ECAI) computation from paper formula.
    
    Implements: ECAI(π) = (∏_{j=1}^M ∑_{c∈C(y_j)} π_{j,c} · w_c)^{1/M}
    """
    
    def __init__(self, species: str = 'ecoli_bl21de3', device: torch.device = None):
        """
        Initialize ECAI computation with CAI weights.
        
        Args:
            species: Species for CAI weights (default: ecoli_bl21de3)
            device: Computation device
        """
        self.species = species
        self.device = device if device is not None else torch.device('cpu')
        
        # Load CAI weights
        self.cai_weights = self._load_cai_weights()
        
        # Create codon to weight mapping
        self._create_codon_weight_mapping()
    
    def _load_cai_weights(self) -> torch.Tensor:
        """Load CAI weights from reference data."""
        # __file__ is id3/constraints/unified_cai_loss.py
        # .parent.parent.parent gets us to project root
        weights_file = Path(__file__).parent.parent.parent / 'data' / 'codon_references' / f'{self.species}_wi_weights_comparison.json'
        
        if not weights_file.exists():
            raise FileNotFoundError(f"CAI weights file not found: {weights_file}")
        
        with open(weights_file, 'r') as f:
            data = json.load(f)
        
        # Convert wi_values_array to tensor
        wi_values = torch.tensor(data['wi_values_array'], dtype=torch.float32, device=self.device)
        codon_order = data['codon_order']
        
        return wi_values, codon_order
    
    def _create_codon_weight_mapping(self):
        """Create mapping from codons to CAI weights."""
        wi_values, codon_order = self.cai_weights
        
        # Convert DNA codons to RNA codons (T -> U)
        rna_codon_order = [codon.replace('T', 'U') for codon in codon_order]
        
        # Create codon index mapping
        self.codon_to_idx = {codon: idx for idx, codon in enumerate(rna_codon_order)}
        
        # Store weights tensor
        self.weights_tensor = wi_values
    
    def compute_ecai_from_codon_probs(self, 
                                      codon_probs: torch.Tensor,
                                      amino_acid_sequence: str,
                                      valid_codon_mask: torch.Tensor,
                                      codon_indices: torch.Tensor) -> torch.Tensor:
        """
        Compute ECAI from codon probability distributions.
        
        Args:
            codon_probs: Codon probabilities [batch_size, num_positions, max_codons]
            amino_acid_sequence: Target amino acid sequence
            valid_codon_mask: Mask for valid codons [num_positions, max_codons]
            codon_indices: Indices of codons for each position [num_positions, max_codons]
            
        Returns:
            ECAI value [batch_size] or scalar if batch_size=1
        """
        # 输入验证：检查codon_probs是否为空或形状错误
        if codon_probs is None or codon_probs.numel() == 0:
            logger.warning(f"codon_probs为空或None: {codon_probs}")
            return torch.tensor(0.0, device=codon_probs.device if codon_probs is not None else 'cpu')
        
        if codon_probs.dim() == 0:
            logger.warning(f"codon_probs维度为0: {codon_probs.shape}")
            return torch.tensor(0.0, device=codon_probs.device)
        
        if codon_probs.dim() == 2:
            codon_probs = codon_probs.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False
        
        if codon_probs.dim() != 3:
            logger.error(f"codon_probs维度错误，期望3维，实际{codon_probs.dim()}维，形状：{codon_probs.shape}")
            return torch.tensor(0.0, device=codon_probs.device)
        
        batch_size, num_positions, max_codons = codon_probs.shape
        
        # Compute expected CAI for each position
        position_expected_cai = []
        
        for pos in range(num_positions):
            if pos >= len(amino_acid_sequence):
                # Handle case where sequence is shorter than expected
                position_expected_cai.append(torch.ones(batch_size, device=self.device))
                continue
            
            # Get valid codons and their probabilities for this position
            valid_mask = valid_codon_mask[pos]  # [max_codons]
            pos_codon_probs = codon_probs[:, pos, :]  # [batch_size, max_codons]
            pos_codon_indices = codon_indices[pos]  # [max_codons]
            
            # Compute weighted sum: Σ_{c∈C(y_j)} π_{j,c} · w_c
            expected_cai_pos = torch.zeros(batch_size, device=self.device)
            
            for codon_slot in range(max_codons):
                if not valid_mask[codon_slot]:
                    continue
                
                codon_idx = pos_codon_indices[codon_slot].item()
                
                # Get CAI weight for this codon
                if 0 <= codon_idx < len(self.weights_tensor):
                    cai_weight = self.weights_tensor[codon_idx]
                    
                    # Add weighted contribution: π_{j,c} · w_c
                    prob_contribution = pos_codon_probs[:, codon_slot] * cai_weight
                    expected_cai_pos += prob_contribution
            
            position_expected_cai.append(expected_cai_pos)
        
        # Stack all positions: [batch_size, num_positions]
        if position_expected_cai:
            position_cai_tensor = torch.stack(position_expected_cai, dim=1)
            
            # Compute geometric mean: (∏_{j=1}^M ∑_{c∈C(y_j)} π_{j,c} · w_c)^{1/M}
            # Using log-exp trick for numerical stability
            log_position_cai = torch.log(torch.clamp(position_cai_tensor, min=1e-10))
            mean_log_cai = torch.mean(log_position_cai, dim=1)
            ecai = torch.exp(mean_log_cai)
        else:
            # Fallback for empty sequence
            ecai = torch.ones(batch_size, device=self.device)
        
        if squeeze_output:
            ecai = ecai.squeeze(0)
        
        return ecai
    
    def compute_ecai_from_rna_probs(self,
                                    rna_probs: torch.Tensor,
                                    amino_acid_sequence: str) -> torch.Tensor:
        """
        Compute ECAI from RNA probabilities by reconstructing codon probabilities.
        
        This is needed for constraints that don't directly provide codon probabilities.
        
        Args:
            rna_probs: RNA probabilities [batch_size, seq_len, 4] or [seq_len, 4]
            amino_acid_sequence: Target amino acid sequence
            
        Returns:
            ECAI value [batch_size] or scalar
        """
        if rna_probs.dim() == 2:
            rna_probs = rna_probs.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False
        
        batch_size, seq_len, _ = rna_probs.shape
        num_positions = seq_len // 3
        
        # Ensure we don't exceed amino acid sequence length
        num_positions = min(num_positions, len(amino_acid_sequence))
        
        # Compute ECAI by directly evaluating codon probabilities
        position_expected_cai = []
        
        # Import amino acid to codons mapping
        from id3.utils.constants import amino_acids_to_codons, codons
        
        for pos in range(num_positions):
            aa = amino_acid_sequence[pos]
            valid_codons = amino_acids_to_codons.get(aa, [])
            
            if not valid_codons:
                position_expected_cai.append(torch.ones(batch_size, device=self.device))
                continue
            
            # Get RNA probabilities for this codon position
            rna_start = pos * 3
            pos_rna_probs = rna_probs[:, rna_start:rna_start+3, :]  # [batch_size, 3, 4]
            
            # Compute probability of each valid codon
            expected_cai_pos = torch.zeros(batch_size, device=self.device)
            
            for codon in valid_codons:
                if codon not in codons:
                    continue
                
                codon_idx = codons.index(codon)
                
                # Get CAI weight
                if codon_idx < len(self.weights_tensor):
                    cai_weight = self.weights_tensor[codon_idx]
                    
                    # Compute codon probability as product of nucleotide probabilities
                    codon_prob = torch.ones(batch_size, device=self.device)
                    
                    for i, nuc in enumerate(codon):
                        nuc_idx = {'A': 0, 'C': 1, 'G': 2, 'U': 3}.get(nuc, 0)
                        codon_prob *= pos_rna_probs[:, i, nuc_idx]
                    
                    # Add weighted contribution
                    expected_cai_pos += codon_prob * cai_weight
            
            position_expected_cai.append(expected_cai_pos)
        
        # Compute geometric mean
        if position_expected_cai:
            position_cai_tensor = torch.stack(position_expected_cai, dim=1)
            log_position_cai = torch.log(torch.clamp(position_cai_tensor, min=1e-10))
            mean_log_cai = torch.mean(log_position_cai, dim=1)
            ecai = torch.exp(mean_log_cai)
        else:
            ecai = torch.ones(batch_size, device=self.device)
        
        if squeeze_output:
            ecai = ecai.squeeze(0)
        
        return ecai


class UnifiedCAILoss(nn.Module):
    """
    Unified CAI Loss implementation using ECAI formula.
    
    Implements: L_CAI(Θ) = (ECAI(T_{ID3-Xαβ}(Θ)) - CAI_target)²
    """
    
    def __init__(self,
                 cai_target: float = 0.8,
                 lambda_cai: float = 1.0,
                 species: str = 'ecoli_bl21de3',
                 device: torch.device = None):
        """
        Initialize unified CAI loss.
        
        Args:
            cai_target: Target CAI value (default: 0.8)
            lambda_cai: Weight for CAI loss vs accessibility loss
            species: Species for CAI weights
            device: Computation device
        """
        super().__init__()
        
        self.cai_target = cai_target
        self.lambda_cai = lambda_cai
        self.device = device if device is not None else torch.device('cpu')
        
        # Initialize ECAI computation
        self.ecai_computer = ECAIComputation(species=species, device=self.device)
    
    def compute_cai_loss_from_codon_probs(self,
                                          codon_probs: torch.Tensor,
                                          amino_acid_sequence: str,
                                          valid_codon_mask: torch.Tensor,
                                          codon_indices: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute CAI loss from codon probability distributions.
        
        Returns:
            cai_loss: L_CAI(Θ) = (ECAI - target)²
            ecai_value: Current ECAI value for monitoring
        """
        # Compute ECAI
        ecai = self.ecai_computer.compute_ecai_from_codon_probs(
            codon_probs, amino_acid_sequence, valid_codon_mask, codon_indices
        )
        
        # Compute squared error loss (without lambda - applied externally)
        cai_loss = (ecai - self.cai_target) ** 2
        
        return cai_loss, ecai
    
    def compute_cai_loss_from_rna_probs(self,
                                        rna_probs: torch.Tensor,
                                        amino_acid_sequence: str) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute CAI loss from RNA probability distributions.
        
        Returns:
            cai_loss: L_CAI(Θ) = (ECAI - target)²
            ecai_value: Current ECAI value for monitoring
        """
        # Compute ECAI
        ecai = self.ecai_computer.compute_ecai_from_rna_probs(rna_probs, amino_acid_sequence)
        
        # Compute squared error loss (without lambda - applied externally)
        cai_loss = (ecai - self.cai_target) ** 2
        
        return cai_loss, ecai
    
    def forward(self,
                rna_probs: torch.Tensor,
                amino_acid_indices: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward method for nn.Module compatibility.
        
        Args:
            rna_probs: RNA probability distributions [batch_size, seq_len, 4]
            amino_acid_indices: Amino acid indices [num_positions]
            
        Returns:
            Tuple of (cai_loss, ecai_value)
        """
        # Convert amino acid indices to sequence string
        from id3.utils.constants import amino_acids_to_codons
        amino_acid_list = list(amino_acids_to_codons.keys())
        
        # Handle batch dimension
        if amino_acid_indices.dim() == 2:
            amino_acid_indices = amino_acid_indices[0]
        
        # Convert indices to amino acid sequence
        amino_acid_sequence = ''.join([
            amino_acid_list[int(idx.item())] 
            for idx in amino_acid_indices
        ])
        
        # Compute CAI loss
        return self.compute_cai_loss_from_rna_probs(rna_probs, amino_acid_sequence)
    
    def compute_unified_loss(self,
                            accessibility_loss: torch.Tensor,
                            cai_loss: torch.Tensor) -> Dict[str, torch.Tensor]:
        """
        Compute unified loss: L_unified = L_Access + λ_CAI * L_CAI
        
        Args:
            accessibility_loss: Accessibility loss from DeepRaccess
            cai_loss: CAI loss from ECAI computation
            
        Returns:
            Dictionary with loss components and total loss
        """
        weighted_cai_loss = self.lambda_cai * cai_loss
        unified_loss = accessibility_loss + weighted_cai_loss
        
        return {
            'unified_loss': unified_loss,
            'accessibility_loss': accessibility_loss,
            'cai_loss': cai_loss,
            'weighted_cai_loss': weighted_cai_loss,
            'lambda_cai': torch.tensor(self.lambda_cai, device=self.device)
        }
    
    def set_lambda_cai(self, lambda_cai: float):
        """Update λ_CAI coefficient."""
        self.lambda_cai = lambda_cai
    
    def get_cai_target(self) -> float:
        """Get current CAI target."""
        return self.cai_target
    
    def set_cai_target(self, cai_target: float):
        """Update CAI target."""
        self.cai_target = cai_target


def create_unified_cai_loss(cai_target: float = 0.8,
                           lambda_cai: float = 1.0,
                           species: str = 'ecoli_bl21de3',
                           device: Union[str, torch.device] = 'cuda') -> UnifiedCAILoss:
    """
    Factory function to create unified CAI loss.
    
    Args:
        cai_target: Target CAI value
        lambda_cai: Weight for CAI loss
        species: Species for CAI weights
        device: Computation device
        
    Returns:
        Configured UnifiedCAILoss instance
    """
    if isinstance(device, str):
        device = torch.device(device)
    
    return UnifiedCAILoss(
        cai_target=cai_target,
        lambda_cai=lambda_cai,
        species=species,
        device=device
    )