"""
Base Constraint Class for ID3 Framework

This module provides the base class for all ID3 constraint mechanisms,
centralizing CAI integration and common functionality to eliminate code duplication.

Author: ID3 Framework Team
Date: 2025-09-15
"""

import torch
import torch.nn as nn
from typing import Optional, Dict, Any, Tuple
from id3.constraints.adaptive_lambda_cai import AdaptiveLambdaCAIMixin


class BaseConstraint(nn.Module, AdaptiveLambdaCAIMixin):
    """
    Base class for all ID3 constraint mechanisms.

    This class provides:
    1. Unified CAI component initialization
    2. Common loss computation interface
    3. Shared precomputation and caching mechanisms

    All constraint types (Lagrangian, AMS, CPC) should inherit from this class
    to ensure consistent CAI integration and reduce code duplication.
    """

    def __init__(
        self,
        amino_acid_sequence: str,
        batch_size: int = 1,
        device: Optional[torch.device] = None,
        # CAI-related parameters
        enable_cai: bool = False,
        cai_target: float = 0.8,
        cai_weight: float = 0.1,
        species: str = 'ecoli_bl21de3',
        # Adaptive lambda_cai parameters
        adaptive_lambda_cai: bool = False,
        lambda_cai_lr: float = 0.1,
        lambda_cai_max: float = 2.0,
        lambda_cai_min: float = 0.01,
        cai_tolerance: float = 0.05,
        smoothing_factor: float = 0.9,
        # Other parameters
        verbose: bool = False,
        **kwargs
    ):
        """
        Initialize base constraint with CAI components.

        Args:
            amino_acid_sequence: Target amino acid sequence
            batch_size: Batch size for optimization
            device: Computation device (cuda/cpu)
            enable_cai: Whether to enable CAI optimization
            cai_target: Target CAI value (0.0-1.0)
            cai_weight: Initial weight for CAI loss
            species: Species for codon usage table
            adaptive_lambda_cai: Whether to enable dynamic lambda_cai adjustment
            lambda_cai_lr: Learning rate for lambda_cai adaptation
            lambda_cai_max: Maximum allowed lambda_cai value
            lambda_cai_min: Minimum allowed lambda_cai value
            cai_tolerance: Tolerance for CAI target satisfaction
            smoothing_factor: Exponential smoothing factor for CAI measurements
            verbose: Whether to print debug information
            **kwargs: Additional arguments for specific constraints
        """
        super().__init__()

        # Basic configuration
        self.amino_acid_sequence = amino_acid_sequence
        self.batch_size = batch_size
        self.num_positions = len(amino_acid_sequence)
        self.seq_length = self.num_positions * 3

        # Device configuration
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        elif isinstance(device, str):
            self.device = torch.device(device)
        else:
            self.device = device

        self.verbose = verbose

        # CAI configuration
        self.enable_cai = enable_cai
        self.cai_target = cai_target
        self.species = species

        # Initialize CAI components
        self._init_cai_components(
            cai_weight=cai_weight,
            adaptive_lambda_cai=adaptive_lambda_cai,
            lambda_cai_lr=lambda_cai_lr,
            lambda_cai_max=lambda_cai_max,
            lambda_cai_min=lambda_cai_min,
            cai_tolerance=cai_tolerance,
            smoothing_factor=smoothing_factor,
            **kwargs
        )

        # Cache for last forward pass result
        self._last_result = None

    def _init_cai_components(
        self,
        cai_weight: float,
        adaptive_lambda_cai: bool,
        lambda_cai_lr: float,
        lambda_cai_max: float,
        lambda_cai_min: float,
        cai_tolerance: float,
        smoothing_factor: float,
        **kwargs
    ):
        """
        Initialize all CAI-related components in a unified manner.

        This method centralizes the initialization of:
        1. Adaptive lambda_cai (if enabled)
        2. Unified CAI loss module
        3. CAI enhancement operator for discretization
        """
        # Initialize adaptive lambda_cai if CAI is enabled
        if self.enable_cai:
            self._init_adaptive_lambda_cai(
                initial_lambda_cai=cai_weight,
                adaptive_lambda_cai=adaptive_lambda_cai,
                lambda_cai_lr=lambda_cai_lr,
                lambda_cai_max=lambda_cai_max,
                cai_tolerance=cai_tolerance,
                verbose=self.verbose,
                smoothing_factor=smoothing_factor,
                lambda_cai_min=lambda_cai_min,
                **kwargs
            )
        else:
            self.lambda_cai = cai_weight  # Fixed weight when CAI disabled

        # Initialize unified CAI loss module
        self.cai_loss_module = None
        if self.enable_cai:
            from id3.constraints.unified_cai_loss import UnifiedCAILoss
            self.cai_loss_module = UnifiedCAILoss(
                cai_target=self.cai_target,
                lambda_cai=self.lambda_cai,  # Use potentially adaptive lambda_cai
                species=self.species,
                device=self.device
            )

        # Initialize CAI enhancement operator for β=1 discretization
        self.cai_enhancement_operator = None
        if self.enable_cai:
            from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
            self.cai_enhancement_operator = CAIEnhancementOperator(
                species=self.species,
                device=self.device,
                amino_acid_sequence=self.amino_acid_sequence  # Precomputation optimization
            )

    def compute_unified_loss(
        self,
        accessibility_loss: torch.Tensor,
        constraint_penalty: torch.Tensor,
        cai_loss: Optional[torch.Tensor] = None
    ) -> Dict[str, torch.Tensor]:
        """
        Compute the unified loss function: L_total = L_Access + λ_CAI * L_CAI + constraint_penalty

        This method implements the paper's unified loss formulation,
        combining accessibility, CAI, and constraint penalties.

        Args:
            accessibility_loss: DeepRaccess predicted accessibility loss
            constraint_penalty: Constraint violation penalty (may be weighted)
            cai_loss: CAI loss if CAI is enabled

        Returns:
            Dictionary containing:
                - 'total': Total combined loss
                - 'accessibility': Accessibility component
                - 'constraint': Constraint penalty component
                - 'cai': CAI loss (if enabled)
                - 'weighted_cai': λ_CAI * CAI loss (if enabled)
        """
        loss_components = {
            'accessibility': accessibility_loss,
            'constraint': constraint_penalty
        }

        # Base loss without CAI
        total_loss = accessibility_loss + constraint_penalty

        # Add CAI component if enabled
        if self.enable_cai and cai_loss is not None:
            # Update lambda_cai if adaptive (handled by mixin)
            if hasattr(self, 'adaptive_lambda_cai') and self.adaptive_lambda_cai:
                current_lambda = self.get_current_lambda_cai()
            else:
                current_lambda = self.lambda_cai

            weighted_cai_loss = current_lambda * cai_loss
            loss_components['cai'] = cai_loss
            loss_components['weighted_cai'] = weighted_cai_loss
            loss_components['lambda_cai'] = current_lambda
            total_loss = total_loss + weighted_cai_loss

        loss_components['total'] = total_loss
        return loss_components

    def get_cai_info(self) -> Dict[str, Any]:
        """
        Get current CAI-related information.

        Returns:
            Dictionary with CAI configuration and current state
        """
        info = {
            'enabled': self.enable_cai,
            'target': self.cai_target if self.enable_cai else None,
            'species': self.species if self.enable_cai else None
        }

        if self.enable_cai:
            if hasattr(self, 'adaptive_lambda_cai') and self.adaptive_lambda_cai:
                info['lambda_cai'] = self.get_current_lambda_cai()
                info['adaptive'] = True
            else:
                info['lambda_cai'] = self.lambda_cai
                info['adaptive'] = False

        return info

    def forward(self, *args, **kwargs):
        """
        Forward pass - must be implemented by subclasses.

        Each constraint type (Lagrangian, AMS, CPC) must implement
        its specific forward logic while utilizing the base class
        CAI components and loss computation.
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement forward() method"
        )

    def get_parameters_info(self) -> str:
        """
        Get a string representation of current parameters.

        Returns:
            String with parameter information for debugging
        """
        info = [
            f"Constraint: {self.__class__.__name__}",
            f"Amino acids: {len(self.amino_acid_sequence)}",
            f"Batch size: {self.batch_size}",
            f"Device: {self.device}",
        ]

        if self.enable_cai:
            cai_info = self.get_cai_info()
            info.append(f"CAI enabled: target={cai_info['target']:.2f}, λ={cai_info['lambda_cai']:.4f}")

        return "\n".join(info)

    def reset_cache(self):
        """Reset cached results."""
        self._last_result = None

    def _apply(self, fn):
        """
        Apply function to all tensors (for device transfer).

        This ensures all components move together when calling .to(device)
        """
        super()._apply(fn)

        # Note: CAI modules handle their own device transfer
        # They are nn.Module instances and will be handled automatically

        return self