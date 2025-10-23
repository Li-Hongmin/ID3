"""
CAI-Accessibility Joint Constraint for ID3 Framework

This constraint combines CAI optimization using binary search with
accessibility optimization through gradient descent. It extends the
Lagrangian constraint with CAI-aware discretization.

Mathematical Framework:
    L_total = f_accessibility(S_CAI) + λ·C_amino + γ·(CAI_target - CAI(S))
    
where:
    - S_CAI: CAI-optimized sequence from binary search
    - f_accessibility: DeepRaccess accessibility prediction
    - C_amino: Amino acid constraint penalty
    - γ: CAI constraint weight
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, Optional, Tuple, Any
import numpy as np

# Import base constraint
from .lagrangian import LagrangianConstraint

# Import CAI components
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from id3.cai.cai_psi_function import CAIBinarySearchPsiFunction
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper


class CAIAccessibilityConstraint(LagrangianConstraint):
    """
    Joint CAI-Accessibility Constraint
    
    This constraint optimizes both CAI and accessibility by:
    1. Using binary search for optimal CAI discretization
    2. Computing accessibility on CAI-optimized sequences
    3. Gradient descent on continuous parameters for accessibility
    4. Maintaining amino acid constraints through Lagrangian
    """
    
    def __init__(self,
                 amino_acid_sequence: str,
                 utr5: str,
                 utr3: str,
                 lambda_value: float = 0.01,
                 cai_target: float = 0.8,
                 cai_weight: float = 0.1,
                 device: str = 'cuda',
                 deepraccess_model_path: Optional[str] = None):
        """
        Initialize CAI-Accessibility constraint
        
        Args:
            amino_acid_sequence: Target amino acid sequence
            utr5: 5' UTR sequence
            utr3: 3' UTR sequence
            lambda_value: Lagrangian penalty weight
            cai_target: Target CAI value
            cai_weight: Weight for CAI constraint
            device: Computation device
            deepraccess_model_path: Path to DeepRaccess model
        """
        # Initialize base Lagrangian constraint
        super().__init__(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=1,
            initial_lambda=lambda_value,
            adaptive_lambda=True,
            device=device
        )
        
        # Store UTRs
        self.utr5 = utr5
        self.utr3 = utr3
        
        self.cai_target = cai_target
        self.cai_weight = cai_weight
        
        # Initialize CAI Psi function with binary search
        self.cai_psi = CAIBinarySearchPsiFunction(
            cai_target=cai_target,
            device=device,
            deepraccess_model_path=deepraccess_model_path,
            enable_accessibility=True
        )
        
        # Initialize DeepRaccess wrapper
        self.deepraccess = DeepRaccessID3Wrapper(
            deepraccess_model_path=deepraccess_model_path,
            device=device
        )
        
        # Tracking metrics
        self.current_cai = 0.0
        self.current_accessibility = 0.0
        self.optimization_history = []
        
    def forward(self,
                alpha: float = 0.0,
                beta: float = 1.0,
                tau: float = 1.0,
                compute_penalty: bool = True) -> Dict[str, torch.Tensor]:
        """
        Forward pass with CAI-optimized discretization
        
        Args:
            alpha: Gumbel noise weight (0=deterministic, 1=stochastic)
            beta: Hard/soft discretization (0=soft, 1=hard)
            tau: Temperature for Gumbel softmax
            compute_penalty: Whether to compute constraint penalties
            
        Returns:
            Dictionary containing sequences and metrics
        """
        # Step 1: Get continuous probabilities from Lagrangian
        lagrangian_result = super().forward(
            alpha=alpha,
            beta=0.0,  # Always soft for probability computation
            tau=tau,
            compute_penalty=False  # Compute penalty separately
        )
        
        continuous_probs = lagrangian_result['rna_sequence']
        
        # Remove batch dimension if present
        if continuous_probs.dim() == 3:
            batch_size = continuous_probs.shape[0]
            continuous_probs = continuous_probs.squeeze(0) if batch_size == 1 else continuous_probs
        else:
            batch_size = 1
            continuous_probs = continuous_probs.unsqueeze(0)
        
        # Step 2: Apply CAI optimization with binary search
        cai_sequences, cai_metadata = self.cai_psi.forward(
            rna_probs=continuous_probs,
            amino_sequence=self.amino_acid_sequence,
            hard=(beta > 0.5),  # Hard if beta > 0.5
            compute_accessibility=False  # We compute it separately with UTRs
        )
        
        # Step 3: Apply straight-through estimator for gradient flow
        if beta > 0.5:  # Hard mode
            # Use CAI-optimized discrete sequences but maintain gradients
            sequences = cai_sequences.detach() + continuous_probs - continuous_probs.detach()
        else:  # Soft mode
            # Blend CAI-optimized with continuous
            blend_factor = 0.9  # How much to trust CAI optimization
            sequences = blend_factor * cai_sequences + (1 - blend_factor) * continuous_probs
        
        # Step 4: Compute accessibility
        # Add UTRs to the sequence
        if sequences.dim() == 3:
            cds_seq = sequences[0]  # Take first batch element
        else:
            cds_seq = sequences
        full_sequence = self._add_utrs(cds_seq)
        
        # Add batch dimension for DeepRaccess
        full_sequence_batch = full_sequence.unsqueeze(0)  # [1, seq_len, 4]
        
        # Compute accessibility at ATG position
        atg_position = len(self.utr5)  # ATG starts after UTR5
        accessibility = self.deepraccess.compute_atg_window_accessibility(
            full_sequence_batch,
            atg_position=atg_position,
            window_size=31,  # Default window
            discrete=(beta > 0.5)
        )
        
        # Step 5: Compute penalties if requested
        constraint_penalty = torch.tensor(0.0, device=self.device)
        cai_penalty = torch.tensor(0.0, device=self.device)
        
        if compute_penalty:
            # Amino acid constraint penalty
            with torch.no_grad():
                # Get discrete sequence for checking
                discrete_indices = torch.argmax(cai_sequences[0], dim=-1)
                nucleotides = ['A', 'C', 'G', 'U']
                discrete_seq = ''.join([nucleotides[idx] for idx in discrete_indices.cpu().numpy()])
                
                # Check amino acid sequence
                amino_acids = self._rna_string_to_amino_acids(discrete_seq)
                
                # Count mismatches
                mismatches = sum(1 for a, b in zip(amino_acids, self.amino_acid_sequence) if a != b)
                constraint_penalty = torch.tensor(mismatches, dtype=torch.float32, device=self.device)
            
            # CAI constraint penalty
            current_cai = cai_metadata.get('avg_cai', 0.0)
            if current_cai < self.cai_target:
                cai_penalty = torch.tensor(
                    self.cai_target - current_cai,
                    dtype=torch.float32,
                    device=self.device
                )
        
        # Update tracking
        self.current_cai = cai_metadata.get('avg_cai', 0.0)
        self.current_accessibility = accessibility.item() if isinstance(accessibility, torch.Tensor) else accessibility
        self.lambda_value = self.lagrangian_lambda.item()  # Get current lambda value
        
        # Build result dictionary
        result = {
            'rna_sequence': sequences,
            'constraint_penalty': constraint_penalty,
            'cai_penalty': cai_penalty,
            'accessibility': accessibility,
            'cai_score': torch.tensor(self.current_cai, device=self.device),
            'discrete_indices': torch.argmax(cai_sequences[0], dim=-1) if batch_size == 1 else torch.argmax(cai_sequences, dim=-1)
        }
        
        # Add metadata
        result.update({
            'cai_metadata': cai_metadata,
            'amino_acid_match': amino_acids == self.amino_acid_sequence if compute_penalty else None,
            'optimization_step': len(self.optimization_history)
        })
        
        # Record history
        self.optimization_history.append({
            'step': len(self.optimization_history),
            'cai': self.current_cai,
            'accessibility': self.current_accessibility,
            'constraint_penalty': constraint_penalty.item(),
            'cai_penalty': cai_penalty.item()
        })
        
        return result
    
    def compute_total_loss(self,
                          accessibility: torch.Tensor,
                          constraint_penalty: torch.Tensor,
                          cai_penalty: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Compute total loss combining accessibility, amino acid, and CAI constraints
        
        Args:
            accessibility: Accessibility score
            constraint_penalty: Amino acid constraint penalty
            cai_penalty: CAI constraint penalty
            
        Returns:
            Total loss
        """
        # Base loss from accessibility and amino acid constraint
        loss = super().compute_total_loss(accessibility, constraint_penalty)
        
        # Add CAI penalty if provided
        if cai_penalty is not None and self.cai_weight > 0:
            loss = loss + self.cai_weight * cai_penalty
        
        return loss
    
    def get_optimization_summary(self) -> Dict[str, Any]:
        """
        Get summary of optimization progress
        
        Returns:
            Dictionary with optimization metrics
        """
        if not self.optimization_history:
            return {}
        
        history = self.optimization_history
        
        return {
            'initial_cai': history[0]['cai'] if history else 0.0,
            'final_cai': history[-1]['cai'] if history else 0.0,
            'initial_accessibility': history[0]['accessibility'] if history else 0.0,
            'final_accessibility': history[-1]['accessibility'] if history else 0.0,
            'cai_improvement': history[-1]['cai'] - history[0]['cai'] if len(history) > 1 else 0.0,
            'accessibility_improvement': history[0]['accessibility'] - history[-1]['accessibility'] if len(history) > 1 else 0.0,
            'total_steps': len(history),
            'cai_target': self.cai_target,
            'cai_achieved': history[-1]['cai'] >= self.cai_target if history else False,
            'history': history
        }
    
    def reset_history(self):
        """Reset optimization history for new experiment"""
        self.optimization_history = []
        self.current_cai = 0.0
        self.current_accessibility = 0.0
    
    def _rna_string_to_amino_acids(self, rna_seq: str) -> str:
        """
        Convert RNA sequence string to amino acid sequence
        
        Args:
            rna_seq: RNA sequence string
            
        Returns:
            Amino acid sequence string
        """
        codon_table = {
            'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
            'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
            'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        # Convert T to U if needed
        rna_seq = rna_seq.replace('T', 'U')
        
        # Translate to amino acids
        amino_acids = []
        for i in range(0, len(rna_seq) - 2, 3):
            codon = rna_seq[i:i+3]
            amino_acids.append(codon_table.get(codon, 'X'))
        
        return ''.join(amino_acids)
    
    def _add_utrs(self, cds_sequence: torch.Tensor) -> torch.Tensor:
        """
        Add UTRs to CDS sequence
        
        Args:
            cds_sequence: CDS one-hot tensor [seq_len, 4]
            
        Returns:
            Full sequence with UTRs [full_len, 4]
        """
        # Convert UTRs to one-hot
        rna_to_idx = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
        
        # Create UTR5 one-hot
        utr5_len = len(self.utr5)
        utr5_tensor = torch.zeros(utr5_len, 4, device=self.device)
        for i, nt in enumerate(self.utr5):
            if nt in rna_to_idx:
                utr5_tensor[i, rna_to_idx[nt]] = 1.0
        
        # Create UTR3 one-hot
        utr3_len = len(self.utr3)
        utr3_tensor = torch.zeros(utr3_len, 4, device=self.device)
        for i, nt in enumerate(self.utr3):
            if nt in rna_to_idx:
                utr3_tensor[i, rna_to_idx[nt]] = 1.0
        
        # Concatenate UTR5 + CDS + UTR3
        full_sequence = torch.cat([utr5_tensor, cds_sequence, utr3_tensor], dim=0)
        
        return full_sequence


def create_cai_accessibility_constraint(amino_sequence: str,
                                       utr5: str,
                                       utr3: str,
                                       cai_target: float = 0.8,
                                       lambda_value: float = 0.01,
                                       cai_weight: float = 0.1,
                                       device: str = 'cuda',
                                       deepraccess_path: Optional[str] = None) -> CAIAccessibilityConstraint:
    """
    Create a CAI-Accessibility constraint
    
    Args:
        amino_sequence: Target amino acid sequence
        utr5: 5' UTR sequence
        utr3: 3' UTR sequence
        cai_target: Target CAI value
        lambda_value: Lagrangian penalty weight
        cai_weight: CAI constraint weight
        device: Computation device
        deepraccess_path: Path to DeepRaccess model
        
    Returns:
        Configured constraint
    """
    return CAIAccessibilityConstraint(
        amino_acid_sequence=amino_sequence,
        utr5=utr5,
        utr3=utr3,
        lambda_value=lambda_value,
        cai_target=cai_target,
        cai_weight=cai_weight,
        device=device,
        deepraccess_model_path=deepraccess_path
    )