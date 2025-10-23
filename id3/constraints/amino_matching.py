"""
Amino Matching Softmax (AMS) Constraint Implementation

This module implements the Amino Matching Softmax mechanism as described
in Section 2.3.4 of the paper, which uses similarity-based projection from
RNA parameter space to valid codons.

Mathematical Framework:
    ID3-A: Θ → Π_amino → π → Ψ_π→S → S → f_model → L → ∇ → Θ^(t+1)

where:
    - Θ: Standard RNA parameters (maintains full parameter space)
    - Π_amino: Similarity-based projection to valid codons
    - Inner product similarity: ⟨Θ, E[c]⟩ for codon selection
"""

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
from id3.constraints.base import BaseConstraint


class AminoPiFunction:
    """
    Amino Matching Pi Function: Θ → π

    Implements similarity-based projection from RNA parameters to codon probabilities
    using inner product similarity as described in Eq. (13) of the paper.
    """

    @staticmethod
    def forward(theta: torch.Tensor,
                codon_encodings: torch.Tensor,
                valid_codon_mask: torch.Tensor,
                alpha: float = 0.0,
                tau: float = 1.0) -> torch.Tensor:
        """
        Transform RNA parameters to codon probabilities via similarity matching.

        π_{j,c} = exp(⟨Θ^(j) + α·g^(j), E[c]⟩/τ) / Σ_{c'∈C(y_j)} exp(⟨Θ^(j) + α·g^(j), E[c']⟩/τ)

        Args:
            theta: RNA parameters [batch_size, seq_len, 4] or [seq_len, 4]
            codon_encodings: Valid codon encodings [num_positions, max_codons, 3, 4]
            valid_codon_mask: Mask for valid codons [num_positions, max_codons]
            alpha: Stochasticity control
            tau: Temperature parameter

        Returns:
            π: Codon probabilities [batch_size, num_positions, max_codons]
        """
        # Handle batch dimension
        if theta.dim() == 2:
            theta = theta.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False

        batch_size, seq_len, vocab_size = theta.shape
        num_positions = seq_len // 3
        max_codons = codon_encodings.shape[1]

        # Reshape theta to position-wise: [batch, num_positions, 3, 4]
        theta_reshaped = theta.reshape(batch_size, num_positions, 3, 4)

        # Add Gumbel noise if stochastic
        if alpha > 0:
            uniform = torch.rand_like(theta_reshaped)
            gumbel = -torch.log(-torch.log(uniform + 1e-20) + 1e-20)
            theta_noisy = theta_reshaped + alpha * gumbel
        else:
            theta_noisy = theta_reshaped

        # Compute inner products with all codons
        # theta_noisy: [batch, positions, 3, 4]
        # codon_encodings: [positions, max_codons, 3, 4]

        # Expand for batch computation
        theta_expanded = theta_noisy.unsqueeze(2)  # [batch, positions, 1, 3, 4]
        encodings_expanded = codon_encodings.unsqueeze(0)  # [1, positions, max_codons, 3, 4]

        # Compute inner product: sum over nucleotide positions (3) and bases (4)
        similarities = (theta_expanded * encodings_expanded).sum(dim=(3, 4))  # [batch, positions, max_codons]

        # Apply temperature scaling
        similarities_scaled = similarities / tau

        # Mask invalid codons
        mask_expanded = valid_codon_mask.unsqueeze(0).expand(batch_size, -1, -1)
        similarities_masked = torch.where(
            mask_expanded,
            similarities_scaled,
            torch.tensor(-float('inf'), dtype=similarities_scaled.dtype, device=similarities_scaled.device)
        )

        # Compute softmax over valid codons
        codon_probs = F.softmax(similarities_masked, dim=-1)

        # Zero out numerical errors for invalid codons
        codon_probs = torch.where(mask_expanded, codon_probs, torch.zeros_like(codon_probs))

        if squeeze_output:
            codon_probs = codon_probs.squeeze(0)

        return codon_probs

    @staticmethod
    def forward_with_similarity_scores(theta: torch.Tensor,
                                       codon_encodings: torch.Tensor,
                                       valid_codon_mask: torch.Tensor,
                                       alpha: float = 0.0,
                                       tau: float = 1.0) -> Dict[str, torch.Tensor]:
        """
        Forward pass with additional similarity information for analysis.
        """
        codon_probs = AminoPiFunction.forward(
            theta, codon_encodings, valid_codon_mask, alpha, tau
        )

        # Also return raw similarities for debugging
        if theta.dim() == 2:
            theta = theta.unsqueeze(0)

        batch_size, seq_len, _ = theta.shape
        num_positions = seq_len // 3
        theta_reshaped = theta.reshape(batch_size, num_positions, 3, 4)

        theta_expanded = theta_reshaped.unsqueeze(2)
        encodings_expanded = codon_encodings.unsqueeze(0)
        similarities = (theta_expanded * encodings_expanded).sum(dim=(3, 4))

        return {
            'codon_probs': codon_probs,
            'similarities': similarities,
            'max_similarity': similarities.max(dim=-1)[0],
            'selected_codons': codon_probs.argmax(dim=-1)
        }


class AminoMatchingSoftmax(BaseConstraint):
    """
    Complete Amino Matching Softmax (AMS) constraint mechanism.

    This implements the ID3-A variant from the paper, maintaining standard
    RNA parameters while ensuring amino acid constraints through similarity-based
    projection to valid codons.
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
        Initialize AMS with target amino acid sequence.

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

        # Precompute constraint information
        self._precompute_constraints()

        # Initialize standard RNA parameters Θ
        self.theta = nn.Parameter(
            torch.randn(batch_size, self.seq_length, 4) * 0.1
        )

        # Cache for last forward pass result
        self._last_result = None

    def _precompute_constraints(self):
        """Precompute codon encodings and masks for efficiency."""
        # Maximum number of codons
        max_codons = 6

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

    def _apply(self, fn):
        """Apply function to all tensors (for device transfer)."""
        super()._apply(fn)
        self.valid_codon_mask = fn(self.valid_codon_mask)
        self.codon_encodings = fn(self.codon_encodings)
        self.codon_indices = fn(self.codon_indices)
        return self

    def forward(self,
                alpha: float = 0.0,
                beta: float = 0.0,
                tau: float = 1.0) -> Dict[str, torch.Tensor]:
        """
        Execute AMS transformation: Θ → π → S

        Args:
            alpha: Stochasticity control
            beta: Output type (0=soft, 1=hard)
            tau: Temperature parameter

        Returns:
            Dictionary containing RNA sequence and intermediate values
        """
        # Step 1: Amino Pi function (Θ → π via similarity)
        codon_probs = AminoPiFunction.forward(
            self.theta,
            self.codon_encodings,
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


        enhanced_codon_dist = None


        from id3.constraints.codon_profile import CodonPsiFunction
        rna_sequence, selected_codons = CodonPsiFunction.forward(
            codon_probs,
            self.codon_encodings,
            beta
        )


        if self.enable_cai and self.cai_enhancement_operator is not None:

            # Convert 3D codon_probs to 2D for CAI enhancement
            if codon_probs.dim() == 3:
                pi_accessibility = codon_probs[0]  # Take first batch: [positions, max_codons]
            else:
                pi_accessibility = codon_probs


            discrete_dist, cai_metadata_inner = self.cai_enhancement_operator.apply_cai_enhancement(
                pi_accessibility, self.amino_acid_sequence,
                self.valid_codon_mask, self.codon_indices, self.cai_target
            )


            if cai_metadata_inner:
                cai_metadata.update(cai_metadata_inner)


            if hasattr(self, 'verbose') and self.verbose and beta > 0.5:
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
            'theta': self.theta
        }




        if enhanced_codon_dist is not None:
            result['enhanced_codon_dist'] = enhanced_codon_dist


            with torch.no_grad():

                if enhanced_codon_dist.dim() == 2:
                    enhanced_dist_batch = enhanced_codon_dist.unsqueeze(0)
                else:
                    enhanced_dist_batch = enhanced_codon_dist


                from id3.constraints.codon_profile import CodonReconstruction
                enhanced_rna_sequence = CodonReconstruction.forward(
                    enhanced_dist_batch,
                    self.codon_encodings
                )


                if enhanced_rna_sequence.dim() == 3 and enhanced_rna_sequence.shape[0] == 1:
                    enhanced_rna_sequence = enhanced_rna_sequence.squeeze(0)


                from id3.utils.constants import NUCLEOTIDES
                indices = torch.argmax(enhanced_rna_sequence, dim=-1)
                discrete_sequence_str = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])


            result['discrete_sequence'] = discrete_sequence_str
        else:

            from id3.utils.constants import NUCLEOTIDES
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
        self._last_beta = beta

        return result

    def compute_similarity_matrix(self) -> torch.Tensor:
        """
        Compute the similarity matrix between RNA parameters and valid codons.
        Useful for visualization and analysis.
        """
        with torch.no_grad():
            result = AminoPiFunction.forward_with_similarity_scores(
                self.theta,
                self.codon_encodings,
                self.valid_codon_mask,
                alpha=0.0,
                tau=1.0
            )
        return result['similarities']

    def get_discrete_sequence(self) -> str:
        """
        Get the current discrete RNA sequence.





        """
        if hasattr(self, '_last_result') and self._last_result is not None:

            return self._last_result.get('discrete_sequence', '')
        else:


            raise RuntimeError("No cached result available. Please call forward() first.")


    # Removed _apply_cai_to_sequence method - replaced with unified ECAI loss

    def project_to_valid_codons(self) -> torch.Tensor:
        """
        Project current RNA parameters to nearest valid codons.
        This shows the discrete sequence that would be selected.
        """
        with torch.no_grad():
            result = self.forward(alpha=0.0, beta=1.0, tau=0.01)  # Low temperature for sharp selection
        return result['selected_codons']

    def compute_total_loss(self,
                          accessibility_loss: torch.Tensor,
                          codon_probs: torch.Tensor = None,
                          enhanced_codon_dist: torch.Tensor = None) -> Dict[str, torch.Tensor]:
        """
        Compute total loss with optional CAI loss for AMS constraint.

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
                            print(f"🎯 AMS Lambda CAI adapted: {update_stats['lambda_cai']:.4f} "
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
