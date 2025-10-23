"""
Lagrangian Multiplier Constraint Implementation

ID3-L: Θ → Π → P → Ψ → S → f_model + λC → L → ∇ → Θ^(t+1)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, Tuple, Optional
from id3.utils.functions import (
    amino_acid_token_map,
    codon_to_rna_matrix,
    amino_acid_to_codon_matrix
)
from id3.utils.sequence_utils import rna_to_amino_acids
from id3.utils.constraint_satisfied_argmax import get_constraint_satisfied_argmax
from id3.constraints.base import BaseConstraint


class SimplePiFunction:
    """Pi function: Θ → P with Gumbel noise for exploration"""
    
    @staticmethod
    def forward(theta: torch.Tensor, alpha: float = 0.0, tau: float = 1.0) -> torch.Tensor:
        """
        Apply softmax with optional Gumbel noise.
        
        Args:
            theta: Input logits
            alpha: Gumbel noise weight (0=deterministic, 1=full noise)
            tau: Temperature parameter
            
        Returns:
            Probabilities after softmax (with optional Gumbel noise)
        """
        # Scale by temperature
        theta_scaled = theta / tau
        
        # Add Gumbel noise if alpha > 0 (matching AMS/CPC implementation)
        if alpha > 0:
            gumbel_noise = -torch.log(-torch.log(torch.rand_like(theta_scaled) + 1e-10) + 1e-10)
            theta_scaled = theta_scaled + alpha * gumbel_noise
        
        return F.softmax(theta_scaled, dim=-1)


class ConstraintPenalty:
    """Constraint penalty C(P): (1/M) Σ_j min_{c ∈ C(y_j)} ||P^(j) - E[c]||^2"""
    
    @staticmethod
    def compute_vectorized(rna_probs: torch.Tensor,
                          codon_encodings: torch.Tensor,
                          valid_codon_mask: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        if rna_probs.dim() == 2:
            rna_probs = rna_probs.unsqueeze(0)
            squeeze_output = True
        else:
            squeeze_output = False
        
        batch_size, seq_len, _ = rna_probs.shape
        num_positions = seq_len // 3
        max_codons = codon_encodings.shape[1]
        
        rna_reshaped = rna_probs.reshape(batch_size, num_positions, 3, 4)
        rna_expanded = rna_reshaped.unsqueeze(2)
        encodings_expanded = codon_encodings.unsqueeze(0)
        
        squared_diff = (rna_expanded - encodings_expanded) ** 2
        distances = squared_diff.sum(dim=(3, 4))
        
        mask_expanded = valid_codon_mask.unsqueeze(0).expand(batch_size, -1, -1)
        distances = torch.where(mask_expanded, distances, torch.tensor(1e10, dtype=distances.dtype, device=distances.device))
        
        # Find minimum distance for each position
        min_distances, best_codon_indices = distances.min(dim=-1)  # [batch, positions]
        
        # Average penalty across positions
        penalty = min_distances.mean(dim=1)  # [batch]
        
        if squeeze_output:
            penalty = penalty.squeeze(0)
            best_codon_indices = best_codon_indices.squeeze(0)
        
        return penalty, best_codon_indices
    
    @staticmethod
    def get_best_valid_codons(rna_probs: torch.Tensor,
                              codon_encodings: torch.Tensor,
                              valid_codon_mask: torch.Tensor) -> torch.Tensor:
        """
        Find the nearest valid codon for each position.
        
        Returns:
            best_rna: RNA sequence of nearest valid codons [batch_size, seq_len, 4]
        """
        _, best_indices = ConstraintPenalty.compute_vectorized(
            rna_probs, codon_encodings, valid_codon_mask
        )
        
        batch_size = rna_probs.shape[0] if rna_probs.dim() == 3 else 1
        num_positions = best_indices.shape[-1]
        
        if best_indices.dim() == 1:
            best_indices = best_indices.unsqueeze(0)
        
        # Gather the best codon encodings
        # Create indices for gathering
        batch_idx = torch.arange(batch_size).unsqueeze(1).expand(-1, num_positions)
        pos_idx = torch.arange(num_positions).unsqueeze(0).expand(batch_size, -1)
        
        # Select best codons: [batch, positions, 3, 4]
        best_encodings = codon_encodings[pos_idx, best_indices]
        
        # Reshape to RNA sequence: [batch, seq_len, 4]
        best_rna = best_encodings.reshape(batch_size, num_positions * 3, 4)
        
        if batch_size == 1:
            best_rna = best_rna.squeeze(0)
        
        return best_rna


class LagrangianPsiFunction:
    """
    Enhanced Lagrangian Psi function with dual-path CAI enhancement support.
    
    Implements dual-path architecture:
    - Continuous path: For gradient optimization (controlled by beta)
    - Discrete path: Always CAI-enhanced when CAI is enabled (for evaluation)
    """
    
    @staticmethod
    def forward(probabilities: torch.Tensor,
                amino_acid_sequence: str,
                beta: float = 0.0,
                tau: float = 1.0,
                cai_enhancement_operator=None,
                cai_target: float = 0.8,
                valid_codon_mask: torch.Tensor = None,
                codon_indices: torch.Tensor = None) -> Tuple[torch.Tensor, torch.Tensor, Optional[torch.Tensor], Optional[Dict]]:
        """
        Enhanced Psi function with dual-path CAI support: P → S
        
        Args:
            probabilities: RNA probability distributions [batch, seq_len, 4]
            amino_acid_sequence: Target amino acid sequence
            beta: Output type for continuous path (0=soft, 1=hard)
            cai_enhancement_operator: CAI enhancement operator
            cai_target: Target CAI value
            valid_codon_mask: Valid codon mask
            codon_indices: Codon indices
            
        Returns:
            sequence: RNA sequence representation (for gradient optimization)
            discrete_indices: Discrete sequence indices
            enhanced_sequence: CAI-enhanced discrete sequence (when CAI enabled)
            cai_metadata: CAI enhancement metadata
        """
        # Import at the beginning of the function to ensure availability
        from id3.utils.constraint_satisfied_argmax import get_constraint_satisfied_argmax
        
        # Continuous path: For gradient optimization (controlled by beta)
        if beta == 0:
            # Soft mode: Return continuous probabilities
            sequence = probabilities
            discrete_indices = torch.argmax(probabilities, dim=-1)
        else:
            # Hard mode: Use constraint-satisfied argmax with STE
            hard_onehot = get_constraint_satisfied_argmax(probabilities, amino_acid_sequence)
            # Use STE for gradient flow
            sequence = hard_onehot + probabilities - probabilities.detach()
            discrete_indices = torch.argmax(hard_onehot, dim=-1)
        

        enhanced_sequence = None
        cai_metadata = None
        
        if (cai_enhancement_operator is not None and 
            valid_codon_mask is not None and 
            codon_indices is not None):

            # Convert RNA probabilities to codon probabilities
            codon_probs = LagrangianPsiFunction._convert_rna_to_codon_probs(
                probabilities, amino_acid_sequence, valid_codon_mask, codon_indices, tau
            )
            
            # Apply Ψ_{π→τ} operator for CAI enhancement
            discrete_codon_dist, cai_metadata = cai_enhancement_operator.apply_cai_enhancement(
                codon_probs, amino_acid_sequence, valid_codon_mask, 
                codon_indices, cai_target
            )
            
            # Convert discrete codon distribution back to RNA representation
            enhanced_sequence = LagrangianPsiFunction._convert_codon_to_rna(
                discrete_codon_dist, valid_codon_mask, codon_indices
            )
        else:


            enhanced_sequence = get_constraint_satisfied_argmax(probabilities, amino_acid_sequence)
        
        return sequence, discrete_indices, enhanced_sequence, cai_metadata
    
    @staticmethod
    def _convert_rna_to_codon_probs(rna_probs: torch.Tensor,
                                  amino_acid_sequence: str,
                                  valid_codon_mask: torch.Tensor,
                                  codon_indices: torch.Tensor,
                                  tau: float = 1.0) -> torch.Tensor:
        """
        Convert RNA probabilities to codon probabilities using similarity matching.
        
        This implements Π_amino from the paper for Lagrangian constraint.
        🚀 FIXED: Now using proper inner product computation and temperature scaling like AMS.
        """
        # Handle batch dimension
        squeeze_output = False
        if rna_probs.dim() == 2:
            rna_probs = rna_probs.unsqueeze(0)  # Add batch dimension
            squeeze_output = True
        elif rna_probs.dim() == 3:
            pass  # Already has batch dimension
        else:
            raise ValueError(f"Invalid rna_probs dimensions: {rna_probs.shape}")
        
        batch_size, seq_len, vocab_size = rna_probs.shape
        num_positions = seq_len // 3
        max_codons = valid_codon_mask.shape[1]
        
        # Import codon encodings and build codon encoding matrix
        from id3.utils.functions import codon_to_rna_matrix
        codon_to_rna_device = codon_to_rna_matrix.to(rna_probs.device)
        
        # Build codon encodings tensor: [num_positions, max_codons, 3, 4]
        codon_encodings = torch.zeros(num_positions, max_codons, 3, 4, 
                                    dtype=rna_probs.dtype, device=rna_probs.device)
        
        for pos in range(num_positions):
            if pos >= len(amino_acid_sequence):
                break
            for i in range(max_codons):
                if valid_codon_mask[pos, i]:
                    codon_idx = codon_indices[pos, i]
                    if codon_idx < len(codon_to_rna_device):
                        codon_encodings[pos, i] = codon_to_rna_device[codon_idx]
        
        # Reshape rna_probs to position-wise: [batch, num_positions, 3, 4]
        rna_probs_reshaped = rna_probs.reshape(batch_size, num_positions, 3, 4)
        
        # 🚀 Compute inner products using proper batch operations like AMS
        # rna_probs_reshaped: [batch, positions, 3, 4]
        # codon_encodings: [positions, max_codons, 3, 4]
        # Expand for batch computation
        rna_expanded = rna_probs_reshaped.unsqueeze(2)  # [batch, positions, 1, 3, 4]
        encodings_expanded = codon_encodings.unsqueeze(0)  # [1, positions, max_codons, 3, 4]
        
        # Compute inner product: sum over nucleotide positions (3) and bases (4)
        similarities = (rna_expanded * encodings_expanded).sum(dim=(3, 4))  # [batch, positions, max_codons]
        
        # 🚀 Apply temperature scaling like AMS
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
    def _convert_codon_to_rna(discrete_codon_dist: torch.Tensor,
                            valid_codon_mask: torch.Tensor,
                            codon_indices: torch.Tensor) -> torch.Tensor:
        """
        Convert discrete codon distribution to RNA one-hot representation.
        """
        # Handle batch dimension if present (from incremental optimizer)
        if discrete_codon_dist.dim() == 3:
            # Remove batch dimension (assuming batch_size=1)
            discrete_codon_dist = discrete_codon_dist.squeeze(0)
        
        num_positions, max_codons = discrete_codon_dist.shape
        seq_len = num_positions * 3
        
        # Import codon encodings
        from id3.utils.functions import codon_to_rna_matrix
        codon_to_rna_device = codon_to_rna_matrix.to(discrete_codon_dist.device)
        
        rna_sequence = torch.zeros(seq_len, 4, device=discrete_codon_dist.device)
        
        for pos in range(num_positions):
            # Find selected codon
            codon_dist_at_pos = discrete_codon_dist[pos]
            # Ensure it's 1D (in case squeeze didn't work)
            if codon_dist_at_pos.dim() > 1:
                codon_dist_at_pos = codon_dist_at_pos.squeeze()
            selected_slot = torch.argmax(codon_dist_at_pos).item()
            
            if valid_codon_mask[pos, selected_slot]:
                codon_idx = codon_indices[pos, selected_slot]
                if codon_idx < len(codon_to_rna_device):
                    codon_encoding = codon_to_rna_device[codon_idx]  # [3, 4]
                    
                    # Copy codon encoding to RNA sequence
                    rna_start = pos * 3
                    rna_sequence[rna_start:rna_start+3] = codon_encoding
        
        return rna_sequence


class LagrangianConstraint(BaseConstraint):
    """
    Complete Lagrangian Multiplier constraint mechanism.
    
    This implements the ID3-L variant from the paper, using soft penalties
    to enforce amino acid constraints while maintaining full parameter flexibility.
    """
    
    def __init__(self,
                 amino_acid_sequence: str,
                 batch_size: int = 1,
                 initial_lambda: float = 0.01,
                 adaptive_lambda: bool = True,
                 lambda_lr: float = 0.01,
                 lambda_max: float = 10.0,
                 device: torch.device = None,
                 enable_cai: bool = False,
                 cai_target: float = 0.8,
                 cai_weight: float = 0.1,
                 verbose: bool = False,
                 adaptive_lambda_cai: bool = False,
                 lambda_cai_lr: float = 0.1,
                 lambda_cai_max: float = 2.0,
                 cai_tolerance: float = 0.05,
                 **kwargs):
        """
        Initialize Lagrangian constraint.
        
        Args:
            amino_acid_sequence: Target amino acid sequence
            batch_size: Batch size for optimization
            initial_lambda: Initial Lagrangian multiplier value
            adaptive_lambda: Whether to adapt lambda during optimization
            lambda_lr: Learning rate for lambda adaptation
            lambda_max: Maximum value for lambda
            device: Device for computation
            enable_cai: Whether to enable CAI optimization
            cai_target: Target CAI value (0.0-1.0)
            cai_weight: Initial weight for CAI vs probability trade-off
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

        # Lagrangian-specific parameters
        self.lambda_value = initial_lambda
        self.adaptive_lambda = adaptive_lambda
        self.lambda_lr = lambda_lr
        self.lambda_max = lambda_max
        
        # Precompute constraints
        self._precompute_constraints()
        
        # Initialize standard RNA parameters Θ
        self.theta = nn.Parameter(
            torch.randn(batch_size, self.seq_length, 4, device=self.device) * 0.1
        )
        
        # Initialize transformation functions
        self.pi_function = SimplePiFunction()
        
        # Cache for last forward pass result
        self._last_result = None
        # Cache for loss components
        self._last_loss_components = None

    def _precompute_constraints(self):
        """Precompute codon encodings and masks."""
        max_codons = 6
        
        self.valid_codon_mask = torch.zeros(self.num_positions, max_codons, dtype=torch.bool, device=self.device)
        self.codon_encodings = torch.zeros(self.num_positions, max_codons, 3, 4, device=self.device)
        

        codon_to_rna_device = codon_to_rna_matrix.to(self.device)
        amino_to_codon_device = amino_acid_to_codon_matrix.to(self.device)
        
        for pos, aa in enumerate(self.amino_acid_sequence):
            aa_idx = amino_acid_token_map[aa]
            valid_codons = amino_to_codon_device[aa_idx].nonzero(as_tuple=True)[0]
            
            num_valid = len(valid_codons)
            self.valid_codon_mask[pos, :num_valid] = True
            
            for i, codon_idx in enumerate(valid_codons):
                self.codon_encodings[pos, i] = codon_to_rna_device[codon_idx]
        

        self.codon_encodings = self.codon_encodings.to(self.device)
        self.valid_codon_mask = self.valid_codon_mask.to(self.device)
        
        # Generate codon indices for CAI enhancement
        self.codon_indices = self._generate_codon_indices()
    
    def _generate_codon_indices(self) -> torch.Tensor:
        """Generate codon indices tensor for CAI enhancement - FIXED to match AMS/CPC method."""
        max_codons = 6
        codon_indices = torch.zeros(self.num_positions, max_codons, dtype=torch.long, device=self.device)
        
        # 🚀 FIXED: Use same method as AMS/CPC constraints for consistency
        for pos, aa in enumerate(self.amino_acid_sequence):
            aa_idx = amino_acid_token_map[aa]
            valid_codons = amino_acid_to_codon_matrix[aa_idx].nonzero(as_tuple=True)[0]
            
            num_valid = len(valid_codons)
            if num_valid > max_codons:
                num_valid = max_codons
            
            # Use the same indexing system as AMS/CPC
            codon_indices[pos, :num_valid] = valid_codons[:num_valid]
        
        return codon_indices
    
    def forward(self,
                alpha: float = 0.0,
                beta: float = 0.0,
                tau: float = 1.0,
                compute_penalty: bool = True) -> Dict[str, torch.Tensor]:
        """
        Execute Lagrangian transformation: Θ → P → S with penalty.
        
        Args:
            alpha: ID3 stochasticity control (original ID3 parameter)
            beta: Output type
            tau: Temperature
            compute_penalty: Whether to compute constraint penalty
            
        Returns:
            Dictionary containing sequence, penalty, and intermediate values
        """
        # Step 1: Standard Pi function (Θ → P)
        probabilities = self.pi_function.forward(self.theta, alpha, tau)
        
        # Step 1.5: Prepare for unified CAI loss computation
        cai_metadata = None
        if self.enable_cai and self.cai_loss_module is not None:
            cai_metadata = {
                'cai_enabled': True,
                'lambda_cai': self.lambda_cai,
                'cai_target': self.cai_target,
                'method': 'unified_loss_ecai'
            }
        
        # Step 2: Dual-path architecture with CAI support (P → S)
        sequence, discrete_indices, enhanced_sequence, cai_metadata = LagrangianPsiFunction.forward(
            probabilities,
            self.amino_acid_sequence,
            beta,
            tau,
            cai_enhancement_operator=self.cai_enhancement_operator if self.enable_cai else None,
            cai_target=self.cai_target if self.enable_cai else 0.8,
            valid_codon_mask=self.valid_codon_mask if self.enable_cai else None,
            codon_indices=self.codon_indices if self.enable_cai else None
        )
        
        

        if enhanced_sequence is not None:

            discrete_sequence_for_eval = enhanced_sequence

            # cai_metadata is already the correct metadata from LagrangianPsiFunction
        else:

            discrete_sequence_for_eval = sequence
        

        from id3.utils.constants import NUCLEOTIDES
        with torch.no_grad():
            if discrete_sequence_for_eval.dim() == 3:
                seq_for_string = discrete_sequence_for_eval[0]
            elif discrete_sequence_for_eval.dim() == 2:
                seq_for_string = discrete_sequence_for_eval
            else:
                seq_for_string = discrete_sequence_for_eval
            
            if seq_for_string.dim() == 2:  # [seq_len, 4]
                indices = torch.argmax(seq_for_string, dim=-1)
                discrete_sequence_str = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
            else:
                discrete_sequence_str = ''
        
        result = {
            'rna_sequence': sequence,
            'probabilities': probabilities,
            'discrete_indices': discrete_indices,
            'discrete_sequence': discrete_sequence_str,
            'theta': self.theta,
            'lambda': torch.tensor(self.lambda_value)
        }
        
        # Add enhanced distribution if available
        if enhanced_sequence is not None:
            result['enhanced_sequence'] = enhanced_sequence
        
        # Add CAI metadata if available
        if cai_metadata:
            result['cai_metadata'] = cai_metadata
            if cai_metadata and 'final_cai' in cai_metadata:
                result['discrete_cai_value'] = cai_metadata['final_cai']
        
        # Step 3: Compute constraint penalty if needed
        if compute_penalty:

            penalty, best_codons = ConstraintPenalty.compute_vectorized(
                probabilities,
                self.codon_encodings.to(probabilities.device),
                self.valid_codon_mask.to(probabilities.device)
            )
            result['constraint_penalty'] = penalty
            result['best_valid_codons'] = best_codons
        
        # Cache the result for get_discrete_sequence
        self._last_result = result
        
        return result
    
    def compute_total_loss(self,
                          model_loss: torch.Tensor,
                          constraint_penalty: torch.Tensor,
                          probabilities: torch.Tensor = None,
                          enhanced_sequence: torch.Tensor = None,
                          cai_metadata: Dict = None) -> Dict[str, torch.Tensor]:
        """
        Compute total loss with Lagrangian penalty and optional CAI loss.
        
        Unified Loss: L_total = L_access + λ * C(P) + λ_CAI * L_CAI
        
        Args:
            model_loss: Accessibility loss from DeepRaccess
            constraint_penalty: Constraint violation penalty C(P)
            probabilities: RNA probabilities for CAI loss computation
            
        Returns:
            Dictionary with loss components and total loss
        """
        # Basic Lagrangian loss
        lagrangian_loss = model_loss + self.lambda_value * constraint_penalty
        
        result = {
            'total_loss': lagrangian_loss,
            'accessibility_loss': model_loss,
            'constraint_penalty': constraint_penalty,
            'lagrangian_lambda': torch.tensor(self.lambda_value, device=model_loss.device)
        }
        
        # Add CAI loss if enabled
        if self.enable_cai and self.cai_loss_module is not None and probabilities is not None:

            cai_loss, ecai_value = self.cai_loss_module.compute_cai_loss_from_rna_probs(
                probabilities, self.amino_acid_sequence
            )
            

            if cai_metadata is not None and 'final_cai' in cai_metadata:
                discrete_cai_value = cai_metadata['final_cai']
            else:

                discrete_rna = torch.zeros_like(probabilities)
                max_indices = torch.argmax(probabilities, dim=-1)
                discrete_rna.scatter_(-1, max_indices.unsqueeze(-1), 1.0)
                
                with torch.no_grad():
                    _, discrete_ecai = self.cai_loss_module.compute_cai_loss_from_rna_probs(
                        discrete_rna, self.amino_acid_sequence
                    )
                    discrete_cai_value = discrete_ecai.item()
            
            
            # Compute unified loss: L_access + λ * C(P) + λ_CAI * L_CAI
            unified_total_loss = lagrangian_loss + self.lambda_cai * cai_loss
            
            result.update({
                'total_loss': unified_total_loss,
                'cai_loss': cai_loss,
                'weighted_cai_loss': self.lambda_cai * cai_loss,
                'ecai_value': ecai_value,
                'eval_cai': discrete_cai_value,  # Actual CAI of discrete sequence
                'lambda_cai': torch.tensor(self.lambda_cai, device=model_loss.device),
                'cai_target': torch.tensor(self.cai_target, device=model_loss.device)
            })
            
            # Add lambda_cai adaptation statistics if available
            if hasattr(self, '_last_lambda_cai_stats'):
                result['adaptive_lambda_stats'] = self._last_lambda_cai_stats
        
        return result
    
    def update_lambda(self, constraint_penalty: float, tolerance: float = 0.01):
        """
        Update Lagrangian multiplier using standard subgradient method.
        


        

        λ^(t+1) = max(0, λ^(t) + α_t * g_t)


        
        Args:
            constraint_penalty: Current constraint violation C(P)
            tolerance: Target constraint satisfaction level
        """
        if not self.adaptive_lambda:
            return
        

        if not hasattr(self, 'lambda_iteration'):
            self.lambda_iteration = 0
        

        subgradient = constraint_penalty - tolerance
        


        import math
        step_size = self.lambda_lr / math.sqrt(self.lambda_iteration + 1)
        

        new_lambda = self.lambda_value + step_size * subgradient
        self.lambda_value = max(0, min(new_lambda, self.lambda_max))
        
        self.lambda_iteration += 1
        

        if not hasattr(self, 'lambda_convergence_stats'):
            self.lambda_convergence_stats = {
                'subgradients': [],
                'step_sizes': [],
                'lambda_values': []
            }
        
        self.lambda_convergence_stats['subgradients'].append(float(subgradient))
        self.lambda_convergence_stats['step_sizes'].append(float(step_size))
        self.lambda_convergence_stats['lambda_values'].append(float(self.lambda_value))
    
    def get_lambda_convergence_analysis(self):
        """Get lambda convergence statistics."""
        if not hasattr(self, 'lambda_convergence_stats') or len(self.lambda_convergence_stats['lambda_values']) < 2:
            return {'status': 'insufficient_data'}
        
        import numpy as np
        stats = self.lambda_convergence_stats
        lambdas = np.array(stats['lambda_values'])
        
        return {
            'iterations': len(lambdas),
            'current_lambda': float(lambdas[-1]),
            'lambda_stability': float(np.std(lambdas[-min(20, len(lambdas)):]))
        }
    
    def get_discrete_sequence(self) -> str:
        """Get the current discrete RNA sequence."""
        with torch.no_grad():
            # Use hard mode to get discrete sequence
            result = self.forward(alpha=0.0, beta=1.0, tau=1.0)
                    
            rna_probs = result['rna_sequence']
            
            # Convert to nucleotide string
            if rna_probs.dim() == 3:
                rna_probs = rna_probs[0]
            
            indices = torch.argmax(rna_probs, dim=-1)
            nucleotides = ['A', 'C', 'G', 'U']
            sequence = ''.join([nucleotides[idx] for idx in indices.cpu().numpy()])
            
        return sequence
    
    def check_constraint_satisfaction(self) -> Dict[str, float]:
        """Check current level of constraint satisfaction."""
        with torch.no_grad():
            result = self.forward(alpha=0.0, beta=0.0, tau=1.0, compute_penalty=True)
            penalty = result['constraint_penalty'].item()
            
            # Also check discrete sequence constraint satisfaction
            discrete_result = self.forward(alpha=0.0, beta=1.0, tau=1.0)
            discrete_sequence = self.get_discrete_sequence()
            
            # Verify amino acid sequence
            amino_acids = rna_to_amino_acids(discrete_sequence)
            correct = amino_acids == self.amino_acid_sequence
            
        return {
            'constraint_penalty': penalty,
            'lambda_value': self.lambda_value,
            'amino_acid_correct': correct,
            'discrete_sequence': discrete_sequence,
            'amino_acids': amino_acids
        }