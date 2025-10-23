"""
Codon Profile Constraint (CPC) Implementation

This module implements the Codon Profile Constraint mechanism as described
in Section 2.3.3 of the paper, which operates at the codon level to ensure
amino acid constraint satisfaction by construction.

Mathematical Framework:
    ID3-C: Ω → Π_codon → π → Ψ_π→S → S → f_model → L → ∇ → Ω^(t+1)

where:
    - Ω: Codon-level parameters (instead of nucleotide-level Θ)
    - π: Codon probability distributions
    - Constraint satisfaction guaranteed by construction
"""

from math import e
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, List, Tuple, Optional
from id3.utils.functions import (
    amino_acid_token_map,
    codon_to_rna_matrix,
    amino_acid_to_codon_matrix
)
from id3.utils.sequence_utils import rna_to_amino_acids
from id3.utils.constants import NUCLEOTIDES
from id3.constraints.base import BaseConstraint


class CodonPiFunction:
    """
    Codon-level Pi Function: Ω → π

    Adapts the standard Pi function to operate on codon-level parameters,
    computing probabilities only over valid codons for each amino acid.
    """

    @staticmethod
    def forward(omega: torch.Tensor,
                amino_indices: torch.Tensor,
                valid_codon_mask: torch.Tensor,
                alpha: float = 0.0,
                tau: float = 1.0) -> torch.Tensor:
        """
        Transform codon parameters to codon probability distributions.

        Args:
            omega: Codon-level parameters [batch_size, num_positions, max_codons]
            amino_indices: Amino acid indices for each position [num_positions]
            valid_codon_mask: Mask for valid codons [num_positions, max_codons]
            alpha: Stochasticity control
            tau: Temperature parameter

        Returns:
            π: Codon probabilities [batch_size, num_positions, max_codons]
        """
        # Add Gumbel noise if stochastic
        if alpha > 0:
            uniform = torch.rand_like(omega)
            gumbel = -torch.log(-torch.log(uniform + 1e-20) + 1e-20)
            omega_noisy = omega + alpha * gumbel
        else:
            omega_noisy = omega

        # Apply temperature scaling
        omega_scaled = omega_noisy / tau

        # Mask invalid codons with -inf before softmax
        omega_masked = omega_scaled.clone()
        if omega.dim() == 2:
            omega_masked = torch.where(valid_codon_mask, omega_masked, torch.tensor(-float('inf')))
        else:  # batch dimension
            mask_expanded = valid_codon_mask.unsqueeze(0).expand(omega.shape[0], -1, -1)
            omega_masked = torch.where(mask_expanded, omega_masked, torch.tensor(-float('inf')))

        # Compute softmax only over valid codons
        codon_probs = F.softmax(omega_masked, dim=-1)

        # Zero out any numerical errors for invalid codons (using out-of-place operation)
        if omega.dim() == 2:
            codon_probs = torch.where(valid_codon_mask, codon_probs, torch.zeros_like(codon_probs))
        else:
            mask_expanded = valid_codon_mask.unsqueeze(0).expand(omega.shape[0], -1, -1)
            codon_probs = torch.where(mask_expanded, codon_probs, torch.zeros_like(codon_probs))

        return codon_probs


class CodonReconstruction:
    """
    Reconstruction operator R: π → P

    Transforms codon probability distributions to RNA nucleotide probabilities
    using the reconstruction operator from Eq. (8) in the paper.
    """

    @staticmethod
    def forward(codon_probs: torch.Tensor,
                codon_encodings: torch.Tensor) -> torch.Tensor:
        """
        Reconstruct RNA probabilities from codon probabilities.

        R(π_j) = Σ_{c ∈ C(y_j)} π_{j,c} · E[c]

        Args:
            codon_probs: Codon probabilities [batch_size, num_positions, max_codons]
            codon_encodings: Binary codon encodings [num_positions, max_codons, 3, 4]

        Returns:
            P: RNA nucleotide probabilities [batch_size, seq_len, 4]
        """
        batch_size = codon_probs.shape[0] if codon_probs.dim() == 3 else 1
        num_positions = codon_probs.shape[-2]

        if codon_probs.dim() == 2:
            codon_probs = codon_probs.unsqueeze(0)

        # Reshape for batch multiplication
        # codon_probs: [batch, positions, codons] -> [batch, positions, codons, 1, 1]
        codon_probs_expanded = codon_probs.unsqueeze(-1).unsqueeze(-1)

        # codon_encodings: [positions, codons, 3, 4] -> [1, positions, codons, 3, 4]
        encodings_expanded = codon_encodings.unsqueeze(0)

        # Weighted sum: [batch, positions, 3, 4]
        rna_probs_3d = (codon_probs_expanded * encodings_expanded).sum(dim=2)

        # Reshape to [batch, seq_len, 4]
        rna_probs = rna_probs_3d.reshape(batch_size, num_positions * 3, 4)

        if batch_size == 1:
            rna_probs = rna_probs.squeeze(0)

        return rna_probs


class CodonPsiFunction:
    """
    Codon-level Psi Function: π → S

    Transforms codon probabilities to RNA sequence representation,
    adapted from Eq. (11) in the paper.
    """

    @staticmethod
    def forward(codon_probs: torch.Tensor,
                codon_encodings: torch.Tensor,
                beta: float = 0.0) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Transform codon probabilities to RNA sequences.

        Args:
            codon_probs: Codon probabilities [batch_size, num_positions, max_codons]
            codon_encodings: Binary codon encodings [num_positions, max_codons, 3, 4]
            beta: Output type (0=soft, 1=hard)

        Returns:
            S: RNA sequence representation [batch_size, seq_len, 4]
            selected_codons: Indices of selected codons [batch_size, num_positions]
        """
        # Get selected codon indices
        selected_codons = torch.argmax(codon_probs, dim=-1)

        if beta == 0:
            # Soft mode: use reconstruction operator
            rna_sequence = CodonReconstruction.forward(codon_probs, codon_encodings)
        else:
            # Hard mode: select discrete codons
            batch_size = codon_probs.shape[0] if codon_probs.dim() == 3 else 1
            num_positions = codon_probs.shape[-2]

            if codon_probs.dim() == 2:
                codon_probs = codon_probs.unsqueeze(0)
                selected_codons = selected_codons.unsqueeze(0)

            # Create one-hot for selected codons
            one_hot = F.one_hot(selected_codons, num_classes=codon_probs.shape[-1]).float()

            # Use STE for gradient flow
            codon_probs_ste = one_hot + codon_probs - codon_probs.detach()

            # Reconstruct with STE probabilities
            rna_sequence = CodonReconstruction.forward(codon_probs_ste, codon_encodings)

            if batch_size == 1:
                rna_sequence = rna_sequence.squeeze(0)
                selected_codons = selected_codons.squeeze(0)

        return rna_sequence, selected_codons


class CodonProfileConstraint(BaseConstraint):
    """
    Complete Codon Profile Constraint (CPC) mechanism.

    This implements the complete ID3-C variant from the paper,
    ensuring amino acid constraints are satisfied by construction.
    """

    def __init__(self, amino_acid_sequence: str, batch_size: int = 1,
                 enable_cai: bool = False, cai_target: float = 0.8,
                 cai_weight: float = 0.1, device: str = 'cuda',
                 verbose: bool = False,
                 adaptive_lambda_cai: bool = False,
                 lambda_cai_lr: float = 0.1,
                 lambda_cai_max: float = 2.0,
                 cai_tolerance: float = 0.05,
                 **kwargs):
        """
        Initialize CPC with target amino acid sequence.

        Args:
            amino_acid_sequence: Target amino acid sequence
            batch_size: Batch size for optimization
            enable_cai: Whether to enable CAI optimization
            cai_target: Target CAI value (0.0-1.0)
            cai_weight: Initial weight for CAI vs probability trade-off
            device: Computation device
            adaptive_lambda_cai: Whether to enable dynamic lambda_cai adjustment
            lambda_cai_lr: Learning rate for lambda_cai adaptation
            lambda_cai_max: Maximum allowed lambda_cai value
            cai_tolerance: Tolerance for CAI target satisfaction
        """
        # Call base class constructor with CAI parameters
        super().__init__(
            amino_acid_sequence=amino_acid_sequence,
            batch_size=batch_size,
            device=device,
            enable_cai=enable_cai,
            cai_target=cai_target,
            cai_weight=cai_weight,
            species='ecoli_bl21de3',
            adaptive_lambda_cai=adaptive_lambda_cai,
            lambda_cai_lr=lambda_cai_lr,
            lambda_cai_max=lambda_cai_max,
            cai_tolerance=cai_tolerance,
            verbose=verbose,
            **kwargs
        )

        # Precompute amino acid indices and valid codons
        self._precompute_constraints()

        # Initialize codon-level parameters Ω
        self._initialize_parameters()

        # Cache for last forward pass result
        self._last_result = None

    def _precompute_constraints(self):
        """Precompute constraint information for efficiency."""
        # Get amino acid indices
        self.amino_indices = torch.tensor([
            amino_acid_token_map[aa] for aa in self.amino_acid_sequence
        ])

        # Determine maximum number of codons for any amino acid
        max_codons = 6  # Maximum synonymous codons (e.g., Leucine, Serine)

        # Create valid codon mask and encodings
        self.valid_codon_mask = torch.zeros(self.num_positions, max_codons, dtype=torch.bool)
        self.codon_encodings = torch.zeros(self.num_positions, max_codons, 3, 4)
        self.codon_indices = torch.zeros(self.num_positions, max_codons, dtype=torch.long)

        for pos, aa in enumerate(self.amino_acid_sequence):
            aa_idx = amino_acid_token_map[aa]
            valid_codons = amino_acid_to_codon_matrix[aa_idx].nonzero(as_tuple=True)[0]

            num_valid = len(valid_codons)
            self.valid_codon_mask[pos, :num_valid] = True
            self.codon_indices[pos, :num_valid] = valid_codons

            for i, codon_idx in enumerate(valid_codons):
                self.codon_encodings[pos, i] = codon_to_rna_matrix[codon_idx]

    def _initialize_parameters(self):
        """Initialize codon-level parameters Ω."""
        max_codons = self.valid_codon_mask.shape[1]

        # Initialize parameters for valid codons only
        self.omega = nn.Parameter(
            torch.randn(self.batch_size, self.num_positions, max_codons) * 0.1
        )

    def _apply(self, fn):
        """Apply function to all tensors (for device transfer)."""
        super()._apply(fn)
        self.amino_indices = fn(self.amino_indices)
        self.valid_codon_mask = fn(self.valid_codon_mask)
        self.codon_encodings = fn(self.codon_encodings)
        self.codon_indices = fn(self.codon_indices)
        return self

    def forward(self,
                alpha: float = 0.0,
                beta: float = 0.0,
                tau: float = 1.0) -> Dict[str, torch.Tensor]:
        """
        Execute CPC transformation: Ω → π → S

        Args:
            alpha: Stochasticity control
            beta: Output type
            tau: Temperature

        Returns:
            Dictionary containing RNA sequence and intermediate values
        """
        # Step 1: Codon Pi function (Ω → π)
        codon_probs = CodonPiFunction.forward(
            self.omega,
            self.amino_indices,
            self.valid_codon_mask,
            alpha,
            tau
        )

        # Step 1.5: Prepare for unified CAI loss computation
        cai_metadata = {}
        if self.enable_cai and self.cai_loss_module is not None:
            cai_metadata = {
                'cai_enabled': True,
                'lambda_cai': self.lambda_cai,
                'cai_target': self.cai_target,
                'method': 'unified_loss_ecai'
            }



        cai_metadata = None


        rna_sequence, selected_codons = CodonPsiFunction.forward(
            codon_probs,
            self.codon_encodings,

        )


        if self.enable_cai and self.cai_enhancement_operator is not None:

            # Convert 3D codon_probs to 2D for CAI enhancement
            if codon_probs.dim() == 3:
                pi_accessibility = codon_probs[0]  # Take first batch: [positions, max_codons]
            else:
                pi_accessibility = codon_probs


            discrete_dist, cai_metadata = self.cai_enhancement_operator.apply_cai_enhancement(
                pi_accessibility, self.amino_acid_sequence,
                self.valid_codon_mask, self.codon_indices, self.cai_target
            )

            if self.verbose:
                print(f"\n🔍 CAI Enhancement Debug Info (beta={beta:.3f}):")
                print(f"   Target CAI: {self.cai_target}")
                print(f"   Achieved CAI: {cai_metadata.get('final_cai', 'N/A')}")
                print(f"   Optimal gamma: {cai_metadata.get('optimal_gamma', 'N/A')}")

            enhanced_codon_dist = discrete_dist.clone()
        else:

            if codon_probs.dim() == 3:
                discrete_codon_probs = codon_probs[0]  # Take first batch: [positions, max_codons]
            else:
                discrete_codon_probs = codon_probs
            

            discrete_codon_dist = torch.zeros_like(discrete_codon_probs)
            max_indices = torch.argmax(discrete_codon_probs, dim=-1)
            discrete_codon_dist.scatter_(1, max_indices.unsqueeze(-1), 1.0)
            enhanced_codon_dist = discrete_codon_dist

        result = {
            'rna_sequence': rna_sequence,
            'codon_probs': codon_probs,
            'selected_codons': selected_codons,
            'omega': self.omega
        }




        if enhanced_codon_dist is not None:
            result['enhanced_codon_dist'] = enhanced_codon_dist


            with torch.no_grad():

                if enhanced_codon_dist.dim() == 2:
                    enhanced_dist_batch = enhanced_codon_dist.unsqueeze(0)
                else:
                    enhanced_dist_batch = enhanced_codon_dist


                enhanced_rna_sequence = CodonReconstruction.forward(
                    enhanced_dist_batch,
                    self.codon_encodings
                )


                if enhanced_rna_sequence.dim() == 3 and enhanced_rna_sequence.shape[0] == 1:
                    enhanced_rna_sequence = enhanced_rna_sequence.squeeze(0)


                indices = torch.argmax(enhanced_rna_sequence, dim=-1)
                discrete_sequence_str = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])


            result['discrete_sequence'] = discrete_sequence_str
        else:

            with torch.no_grad():
                if rna_sequence.dim() == 3:
                    rna_seq_for_string = rna_sequence[0]
                else:
                    rna_seq_for_string = rna_sequence

                indices = torch.argmax(rna_seq_for_string, dim=-1)
                discrete_sequence_str = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
                result['discrete_sequence'] = discrete_sequence_str

        # Add CAI metadata if available
        if cai_metadata:
            result['cai_metadata'] = cai_metadata

            if isinstance(cai_metadata, dict) and 'final_cai' in cai_metadata:
                result['discrete_cai_value'] = cai_metadata['final_cai']

        # Cache the result for get_discrete_sequence
        self._last_result = result


        return result

    def get_discrete_sequence(self) -> str:
        """
        Get the current discrete RNA sequence.

        简化版本：直接从缓存的结果中获取离散序列字符串
        如果没有缓存，返回空字符串或抛出异常

        注意：forward()方法现在直接在结果中包含'discrete_sequence'字段
        """
        if hasattr(self, '_last_result') and self._last_result is not None:

            return self._last_result.get('discrete_sequence', '')
        else:







            return ''






    def compute_total_loss(self,
                          accessibility_loss: torch.Tensor,
                          codon_probs: torch.Tensor = None,
                          enhanced_codon_dist: torch.Tensor = None) -> Dict[str, torch.Tensor]:
        """
        Compute total loss with optional CAI loss for CPC constraint.

        Unified Loss: L_total = L_access + λ_CAI * L_CAI

        Args:
            accessibility_loss: Accessibility loss from DeepRaccess
            codon_probs: Codon probabilities for CAI loss computation

        Returns:
            Dictionary with loss components and total loss
        """
        result = {
            'total_loss': accessibility_loss,
            'accessibility_loss': accessibility_loss
        }

        # Add CAI loss if enabled
        if self.enable_cai and self.cai_loss_module is not None and codon_probs is not None:
            cai_loss, ecai_value = self.cai_loss_module.compute_cai_loss_from_codon_probs(
                codon_probs, self.amino_acid_sequence,
                self.valid_codon_mask, self.codon_indices
            )

            # Compute discrete CAI value for evaluation

            discrete_cai_value = None


            if hasattr(self, '_last_result') and self._last_result:
                if 'discrete_cai_value' in self._last_result:
                    discrete_cai_value = self._last_result['discrete_cai_value']
                elif 'cai_metadata' in self._last_result:
                    metadata = self._last_result['cai_metadata']
                    if isinstance(metadata, dict) and 'final_cai' in metadata:
                        discrete_cai_value = metadata['final_cai']


            if discrete_cai_value is None:
                if enhanced_codon_dist is not None:

                    discrete_codon_dist = enhanced_codon_dist
                    if discrete_codon_dist.dim() == 2:
                        discrete_codon_dist = discrete_codon_dist.unsqueeze(0)
                else:

                    discrete_codon_dist = torch.zeros_like(codon_probs)
                    max_indices = torch.argmax(codon_probs, dim=-1)
                    discrete_codon_dist.scatter_(-1, max_indices.unsqueeze(-1), 1.0)


                with torch.no_grad():
                    _, discrete_ecai = self.cai_loss_module.compute_cai_loss_from_codon_probs(
                        discrete_codon_dist, self.amino_acid_sequence,
                        self.valid_codon_mask, self.codon_indices
                    )
                    discrete_cai_value = discrete_ecai.item()
                
                # Update lambda_cai adaptively if enabled
                if hasattr(self, 'adaptive_lambda_cai_controller'):
                    update_stats = self.update_lambda_cai(discrete_cai_value, self.cai_target)
                    
                    # Update the CAI loss module with new lambda_cai
                    if 'status' in update_stats and update_stats['status'] == 'updated':
                        self.cai_loss_module.lambda_cai = self.lambda_cai
                        
                        if self.verbose:
                            print(f"🎯 CPC Lambda CAI adapted: {update_stats['lambda_cai']:.4f} "
                                  f"(CAI gap: {update_stats['cai_gap']:.4f})")
                    
                    # Store adaptation statistics
                    self._last_lambda_cai_stats = update_stats

            # Compute unified loss: L_access + λ_CAI * L_CAI
            unified_total_loss = accessibility_loss + self.lambda_cai * cai_loss

            result.update({
                'total_loss': unified_total_loss,
                'cai_loss': cai_loss,
                'weighted_cai_loss': self.lambda_cai * cai_loss,
                'ecai_value': ecai_value,
                'eval_cai': discrete_cai_value,  # Actual CAI of discrete sequence
                'lambda_cai': torch.tensor(self.lambda_cai, device=accessibility_loss.device),
                'cai_target': torch.tensor(self.cai_target, device=accessibility_loss.device)
            })
            
            # Add lambda_cai adaptation statistics if available
            if hasattr(self, '_last_lambda_cai_stats'):
                result['adaptive_lambda_stats'] = self._last_lambda_cai_stats

        return result
