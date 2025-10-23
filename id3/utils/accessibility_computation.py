#!/usr/bin/env python3
"""
Accessibility Computation Module

Simple interface for computing RNA accessibility using the DeepRaccess model.
"""

import torch
from typing import Tuple, Optional


def compute_accessibility_simple(
    sequence_tensor: torch.Tensor,
    deepraccess_model,
    atg_position: int = 0
) -> torch.Tensor:
    """
    Compute RNA accessibility at ATG translation initiation site.

    Args:
        sequence_tensor: Sequence probability tensor from ID3 model
        deepraccess_model: DeepRaccess model wrapper instance
        atg_position: Position of ATG start codon (default: 0)

    Returns:
        accessibility: Accessibility score at ATG window
    """
    return deepraccess_model.compute_atg_window_accessibility(
        sequence_tensor, 
        atg_position=atg_position, 
        discrete=False
    )