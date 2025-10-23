"""
Factory for creating constraint mechanisms.

This module handles the creation and configuration of different
constraint types (Lagrangian, AMS, CPC) with CAI integration.
"""

import torch
import logging
from typing import Dict, Optional

logger = logging.getLogger(__name__)


class ConstraintFactory:
    """Factory class for creating constraint mechanisms."""

    @staticmethod
    def create_constraint(
        constraint_type: str,
        amino_acid_sequence: str,
        variant: str,
        config: Dict,
        device: torch.device
    ) -> torch.nn.Module:
        """
        Create constraint mechanism with unified CAI integration.

        Args:
            constraint_type: Type of constraint ('lagrangian', 'ams', 'cpc')
            amino_acid_sequence: Target amino acid sequence
            variant: Variant specification (00, 01, 10, 11)
            config: Configuration dictionary with CAI parameters
            device: Torch device for computation

        Returns:
            Constraint instance with unified CAI support

        Raises:
            ValueError: If variant format is invalid or constraint type unknown
        """
        # Validate variant format
        if len(variant) != 2 or not variant.isdigit():
            raise ValueError(f"Variant must be 2 digits, got: {variant}")

        # Parse variant: first digit = deterministic/gumbel, second digit = soft/ste
        use_gumbel = (variant[0] == '1')
        use_ste = (variant[1] == '1')

        # Extract CAI parameters from config
        enable_cai = config.get('enable_cai', False)
        cai_target = config.get('cai_target', 0.8)
        cai_weight = config.get('cai_weight', 0.1)
        species = config.get('species', 'ecoli_bl21de3')

        # Common parameters for all constraint types
        common_params = {
            'amino_acid_sequence': amino_acid_sequence,
            'batch_size': 1,
            'device': device,
            'enable_cai': enable_cai,
            'cai_target': cai_target,
            'cai_weight': cai_weight,
            'species': species,
            'use_gumbel_softmax': use_gumbel,
            'use_straight_through': use_ste,
            'verbose': config.get('verbose', False)
        }

        # Lazy import to avoid circular dependencies
        if constraint_type == 'lagrangian':
            from id3.constraints.lagrangian import LagrangianConstraint
            return LagrangianConstraint(**common_params)
        elif constraint_type == 'ams':
            from id3.constraints.amino_matching import AminoMatchingSoftmax
            return AminoMatchingSoftmax(**common_params)
        elif constraint_type == 'cpc':
            from id3.constraints.codon_profile import CodonProfileConstraint
            return CodonProfileConstraint(**common_params)
        else:
            raise ValueError(f"Unknown constraint type: {constraint_type}")

    @staticmethod
    def get_variant_info(variant: str) -> Dict[str, bool]:
        """
        Parse variant string into configuration.

        Args:
            variant: Two-digit variant string

        Returns:
            Dictionary with variant configuration
        """
        if len(variant) != 2 or not variant.isdigit():
            raise ValueError(f"Invalid variant format: {variant}")

        return {
            'use_gumbel': variant[0] == '1',
            'use_ste': variant[1] == '1',
            'variant_str': variant
        }