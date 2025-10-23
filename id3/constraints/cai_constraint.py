"""
CAI Constraint Implementation for ID3 Framework

This module extends the ID3 constraint system with Codon Adaptation Index (CAI) optimization.
It integrates seamlessly with existing CPC/RAMS/Lagrangian constraints by enhancing the 
Psi function (Ψ) discretization step with CAI-aware codon selection.

Mathematical Framework:
    Enhanced ID3: Θ → Π → P → Ψ_CAI → S_CAI → f_model + C_CAI → L → ∇ → Θ^(t+1)
    
where:
    - Ψ_CAI: CAI-enhanced discretization with dual optimization
    - C_CAI: Combined constraint penalty (original + CAI)
    - S_CAI: CAI-optimized sequence in codon space
"""

import torch
import torch.nn as nn
from typing import Dict, Tuple, Optional, Any
from id3.core.transformations import PsiFunction
from id3.utils.functions import amino_acid_token_map, codon_to_rna_matrix
from id3.utils.sequence_utils import rna_to_amino_acids

# Import CAI optimizers from seq_design_integration
import sys
import os
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

try:
    from seq_design_integration.optimization.cai_optimizer import CAIOptimizer
    from seq_design_integration.optimization.cai_advance_optimizer import CAIAdvanceOptimizer
    from seq_design_integration.optimization.cai_dp_optimal import CAIDynamicProgrammingOptimal
    from id3.cai.module import CAIModule
    CAI_AVAILABLE = True
except ImportError as e:
    print(f"Warning: CAI modules not available: {e}")
    CAI_AVAILABLE = False


class CAIEnhancedPsiFunction:
    """
    CAI-Enhanced Psi Function: Probability-to-Sequence with CAI Optimization
    
    This extends the standard Psi function by adding CAI awareness to the discretization step.
    It maintains compatibility with all existing ID3 variants while providing CAI optimization.
    
    Design Principles:
    1. **Non-intrusive**: Only activates when CAI parameters are provided
    2. **Dual-path**: Works with original constraints + CAI objective
    3. **Codon-aware**: Operates in codon space for precise CAI control
    """
    
    def __init__(self,
                 cai_target: float = 0.8,
                 cai_weight: float = 0.1,
                 cai_method: str = 'retreat',  # 'retreat', 'advance', 'dp_optimal'
                 species: str = 'ecoli_bl21de3',
                 device: str = 'cuda'):
        """
        Initialize CAI-enhanced Psi function.
        
        Args:
            cai_target: Target CAI value (0.0-1.0)
            cai_weight: Weight for CAI vs. probability loss trade-off
            cai_method: CAI optimization method
            species: Reference species for CAI calculation
            device: Computation device
        """
        self.cai_target = cai_target
        self.cai_weight = cai_weight
        self.cai_method = cai_method
        self.device = device
        self.species = species
        
        # Initialize CAI optimizers if available
        if CAI_AVAILABLE:
            self.cai_optimizer = self._get_cai_optimizer(cai_method)
            self.cai_module = CAIModule(reference_species=species, device=torch.device(device))
        else:
            self.cai_optimizer = None
            self.cai_module = None
            
        # Standard Psi function for fallback
        self.standard_psi = PsiFunction()
    
    def _get_cai_optimizer(self, method: str):
        """Get the appropriate CAI optimizer based on method."""
        optimizer_map = {
            'retreat': CAIOptimizer,
            'advance': CAIAdvanceOptimizer,
            'dp_optimal': CAIDynamicProgrammingOptimal
        }
        
        if method in optimizer_map:
            return optimizer_map[method](
                target_cai=self.cai_target,
                species=self.species,
                device=self.device
            )
        else:
            return CAIOptimizer(  # Default fallback
                target_cai=self.cai_target,
                species=self.species,
                device=self.device
            )
    
    def forward(self, 
                probs: torch.Tensor,
                hard: bool = True,
                amino_acid_sequence: Optional[str] = None,
                enable_cai: bool = True,
                **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Enhanced forward pass with CAI optimization.
        
        Args:
            probs: Probability distributions [batch_size, seq_len, vocab_size]
            hard: Whether to use hard (discrete) or soft sampling
            amino_acid_sequence: Target amino acid sequence for CAI optimization
            enable_cai: Whether to enable CAI optimization
            **kwargs: Additional parameters
            
        Returns:
            sequences: Optimized sequences [batch_size, seq_len, vocab_size]
            metadata: CAI optimization details
        """
        metadata = {
            'cai_enabled': enable_cai and CAI_AVAILABLE and amino_acid_sequence is not None,
            'cai_target': self.cai_target,
            'cai_method': self.cai_method,
            'optimization_details': []
        }
        
        # If CAI is not enabled or not available, use standard Psi
        if not metadata['cai_enabled']:
            sequences = self.standard_psi.forward(probs, hard=hard)
            metadata['fallback_reason'] = 'CAI disabled or unavailable'
            return sequences, metadata
        
        # Enhanced CAI-aware discretization
        batch_size = probs.shape[0]
        optimized_sequences = []
        
        for batch_idx in range(batch_size):
            try:
                # Extract single sequence probabilities
                seq_probs = probs[batch_idx]  # [seq_len, vocab_size]
                
                # Convert to codon probability format expected by CAI optimizers
                codon_probabilities = self._convert_rna_probs_to_codon_probs(
                    seq_probs, amino_acid_sequence
                )
                
                # Apply CAI optimization
                optimized_rna_seq, final_cai, cai_details = self.cai_optimizer.optimize_iteratively(
                    amino_acid_sequence=amino_acid_sequence,
                    codon_probabilities=codon_probabilities,
                    max_iterations=100,
                    verbose=False
                )
                
                # Convert back to tensor format
                optimized_tensor = self._convert_rna_sequence_to_tensor(optimized_rna_seq, seq_probs.shape)
                optimized_sequences.append(optimized_tensor)
                
                # Store optimization details
                metadata['optimization_details'].append({
                    'batch_idx': batch_idx,
                    'final_cai': final_cai,
                    'success': cai_details.get('success', False),
                    'method_details': cai_details
                })
                
            except Exception as e:
                # Fallback to standard discretization on error
                print(f"CAI optimization failed for batch {batch_idx}: {e}")
                standard_seq = self.standard_psi.forward(seq_probs.unsqueeze(0), hard=hard)[0]
                optimized_sequences.append(standard_seq)
                
                metadata['optimization_details'].append({
                    'batch_idx': batch_idx,
                    'error': str(e),
                    'used_fallback': True
                })
        
        # Stack optimized sequences
        if optimized_sequences:
            final_sequences = torch.stack(optimized_sequences, dim=0)
        else:
            final_sequences = self.standard_psi.forward(probs, hard=hard)
            metadata['fallback_reason'] = 'All CAI optimizations failed'
        
        return final_sequences, metadata
    
    def _convert_rna_probs_to_codon_probs(self, 
                                         rna_probs: torch.Tensor, 
                                         amino_acid_sequence: str) -> Dict[int, Dict[str, float]]:
        """
        Convert RNA probability distributions to codon probability format.
        
        This is a crucial bridge function that converts ID3's RNA probability
        distributions to the codon-level probabilities expected by CAI optimizers.
        
        Args:
            rna_probs: RNA probabilities [seq_len, 4] (A,C,G,U)
            amino_acid_sequence: Target amino acid sequence
            
        Returns:
            codon_probabilities: {position: {codon: probability}} format
        """
        from id3.utils.constants import amino_acids_to_codons
        
        codon_probabilities = {}
        rna_length = len(amino_acid_sequence) * 3  # Each amino acid = 3 RNA positions
        
        for aa_pos, amino_acid in enumerate(amino_acid_sequence):
            if amino_acid not in amino_acids_to_codons:
                continue
                
            # RNA positions for this amino acid
            start_pos = aa_pos * 3
            end_pos = min(start_pos + 3, rna_probs.shape[0])
            
            if end_pos - start_pos < 3:
                continue  # Skip incomplete codons
                
            # Get RNA probabilities for this codon position
            codon_rna_probs = rna_probs[start_pos:end_pos]  # [3, 4]
            
            # Calculate probability for each possible codon
            available_codons = amino_acids_to_codons[amino_acid]
            codon_probs = {}
            
            # RNA alphabet: A=0, C=1, G=2, U=3
            rna_to_idx = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
            
            for codon_dna in available_codons:
                codon_rna = codon_dna.replace('T', 'U')  # Convert to RNA
                
                # Calculate joint probability P(codon) = P(pos1) × P(pos2) × P(pos3)
                codon_prob = 1.0
                for pos_in_codon, nucleotide in enumerate(codon_rna):
                    if pos_in_codon < codon_rna_probs.shape[0]:
                        nuc_idx = rna_to_idx.get(nucleotide, 3)  # Default to U
                        codon_prob *= codon_rna_probs[pos_in_codon, nuc_idx].item()
                    else:
                        codon_prob *= 0.25  # Uniform fallback
                
                codon_probs[codon_rna] = codon_prob
            
            # Normalize probabilities
            total_prob = sum(codon_probs.values())
            if total_prob > 0:
                codon_probs = {codon: prob / total_prob for codon, prob in codon_probs.items()}
                
            codon_probabilities[aa_pos] = codon_probs
        
        return codon_probabilities
    
    def _convert_rna_sequence_to_tensor(self, 
                                       rna_sequence: str, 
                                       target_shape: torch.Size) -> torch.Tensor:
        """
        Convert optimized RNA sequence back to tensor format.
        
        Args:
            rna_sequence: Optimized RNA sequence string
            target_shape: Target tensor shape [seq_len, vocab_size]
            
        Returns:
            sequence_tensor: One-hot encoded tensor
        """
        # RNA alphabet mapping
        rna_to_idx = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
        seq_len, vocab_size = target_shape
        
        # Create one-hot tensor
        sequence_tensor = torch.zeros(seq_len, vocab_size, device=self.device)
        
        # Fill in the optimized sequence
        for pos, nucleotide in enumerate(rna_sequence[:seq_len]):
            nuc_idx = rna_to_idx.get(nucleotide, 3)  # Default to U
            if nuc_idx < vocab_size:
                sequence_tensor[pos, nuc_idx] = 1.0
        
        # Handle any remaining positions with uniform distribution
        for pos in range(len(rna_sequence), seq_len):
            sequence_tensor[pos] = torch.ones(vocab_size, device=self.device) / vocab_size
        
        return sequence_tensor


class CAIMixedConstraint:
    """
    Mixed Constraint System: Original ID3 Constraints + CAI Optimization
    
    This combines the original constraint penalty (CPC/RAMS/Lagrangian) with CAI objectives
    to create a dual-optimization system that respects both amino acid constraints and CAI targets.
    """
    
    def __init__(self, 
                 original_constraint,
                 cai_target: float = 0.8,
                 cai_weight: float = 0.1,
                 species: str = 'ecoli_bl21de3',
                 device: str = 'cuda'):
        """
        Initialize mixed constraint system.
        
        Args:
            original_constraint: Existing ID3 constraint (CPC/RAMS/Lagrangian)
            cai_target: Target CAI value
            cai_weight: Relative weight for CAI vs. original constraint
            species: Reference species for CAI calculation
            device: Computation device
        """
        self.original_constraint = original_constraint
        self.cai_target = cai_target
        self.cai_weight = cai_weight
        self.device = device
        
        if CAI_AVAILABLE:
            self.cai_module = CAIModule(reference_species=species, device=torch.device(device))
        else:
            self.cai_module = None
    
    def compute_mixed_penalty(self, 
                            rna_probs: torch.Tensor,
                            codon_encodings: torch.Tensor,
                            valid_codon_mask: torch.Tensor,
                            amino_acid_sequence: Optional[str] = None) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        Compute combined constraint penalty.
        
        Args:
            rna_probs: RNA probability distributions
            codon_encodings: Valid codon encodings
            valid_codon_mask: Mask for valid codons
            amino_acid_sequence: Amino acid sequence for CAI calculation
            
        Returns:
            total_penalty: Combined penalty value
            penalty_details: Breakdown of penalty components
        """
        # Compute original constraint penalty
        if hasattr(self.original_constraint, 'compute_vectorized'):
            original_penalty, original_details = self.original_constraint.compute_vectorized(
                rna_probs, codon_encodings, valid_codon_mask
            )
        else:
            # Fallback for constraints without vectorized implementation
            original_penalty = torch.tensor(0.0, device=self.device)
            original_details = {}
        
        penalty_details = {
            'original_penalty': original_penalty.item(),
            'original_details': original_details,
            'cai_penalty': 0.0,
            'cai_details': {}
        }
        
        # Compute CAI penalty if enabled
        cai_penalty = torch.tensor(0.0, device=self.device)
        if self.cai_module and amino_acid_sequence:
            try:
                # Convert probabilities to discrete sequence for CAI calculation
                discrete_sequence = self._probabilities_to_rna_sequence(rna_probs)
                current_cai = self.cai_module.compute_cai_score(discrete_sequence)
                
                # CAI penalty: squared difference from target
                cai_penalty = (current_cai - self.cai_target) ** 2
                
                penalty_details['cai_penalty'] = cai_penalty.item()
                penalty_details['cai_details'] = {
                    'current_cai': current_cai,
                    'target_cai': self.cai_target,
                    'cai_difference': current_cai - self.cai_target
                }
                
            except Exception as e:
                penalty_details['cai_error'] = str(e)
        
        # Combine penalties with weights
        total_penalty = original_penalty + self.cai_weight * cai_penalty
        penalty_details['total_penalty'] = total_penalty.item()
        penalty_details['weights'] = {
            'original_weight': 1.0,
            'cai_weight': self.cai_weight
        }
        
        return total_penalty, penalty_details
    
    def _probabilities_to_rna_sequence(self, rna_probs: torch.Tensor) -> str:
        """
        Convert RNA probabilities to discrete sequence for CAI calculation.
        
        Args:
            rna_probs: RNA probability distributions [batch_size, seq_len, 4] or [seq_len, 4]
            
        Returns:
            rna_sequence: RNA sequence string
        """
        if rna_probs.dim() == 3:
            rna_probs = rna_probs[0]  # Take first batch
        
        # Convert to discrete sequence by taking argmax
        discrete_indices = torch.argmax(rna_probs, dim=-1)  # [seq_len]
        
        # Map indices to nucleotides
        idx_to_rna = ['A', 'C', 'G', 'U']
        rna_sequence = ''.join([idx_to_rna[idx.item()] for idx in discrete_indices])
        
        return rna_sequence


# Convenience functions for integrating with existing ID3 variants
def create_cai_enhanced_constraint(base_constraint_class, 
                                 cai_target: float = 0.8,
                                 cai_weight: float = 0.1,
                                 cai_method: str = 'retreat',
                                 species: str = 'ecoli_bl21de3',
                                 **base_kwargs):
    """
    Factory function to create CAI-enhanced versions of existing constraints.
    
    Usage:
        # Enhance Lagrangian constraint with CAI
        CAILagrangianConstraint = create_cai_enhanced_constraint(
            LagrangianConstraint, cai_target=0.8, cai_weight=0.1
        )
        
        # Use in experiment
        constraint = CAILagrangianConstraint(lambda_penalty=0.01)
    """
    
    class CAIEnhancedConstraint(base_constraint_class):
        def __init__(self, *args, **kwargs):
            # Initialize base constraint
            super().__init__(*args, **kwargs)
            
            # Add CAI enhancement
            self.cai_mixed_constraint = CAIMixedConstraint(
                original_constraint=self,
                cai_target=cai_target,
                cai_weight=cai_weight,
                species=species,
                device=getattr(self, 'device', 'cuda')
            )
            
            # Enhanced Psi function
            self.cai_psi = CAIEnhancedPsiFunction(
                cai_target=cai_target,
                cai_weight=cai_weight,
                cai_method=cai_method,
                species=species,
                device=getattr(self, 'device', 'cuda')
            )
        
        def compute_penalty_with_cai(self, 
                                   rna_probs: torch.Tensor,
                                   codon_encodings: torch.Tensor,
                                   valid_codon_mask: torch.Tensor,
                                   amino_acid_sequence: Optional[str] = None):
            """Compute constraint penalty with CAI enhancement."""
            return self.cai_mixed_constraint.compute_mixed_penalty(
                rna_probs, codon_encodings, valid_codon_mask, amino_acid_sequence
            )
        
        def discretize_with_cai(self, 
                              probs: torch.Tensor,
                              amino_acid_sequence: Optional[str] = None,
                              hard: bool = True):
            """Discretize probabilities with CAI optimization."""
            return self.cai_psi.forward(
                probs, hard=hard, amino_acid_sequence=amino_acid_sequence
            )
    
    CAIEnhancedConstraint.__name__ = f"CAI{base_constraint_class.__name__}"
    CAIEnhancedConstraint.__doc__ = f"""
    CAI-enhanced version of {base_constraint_class.__name__}.
    
    This constraint combines the original {base_constraint_class.__name__} with CAI optimization,
    allowing dual optimization of amino acid constraints and Codon Adaptation Index.
    
    CAI Parameters:
        - cai_target: {cai_target}
        - cai_weight: {cai_weight}
        - cai_method: {cai_method}
        - species: {species}
    """
    
    return CAIEnhancedConstraint