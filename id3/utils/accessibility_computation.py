#!/usr/bin/env python3
"""




"""

import torch
from typing import Tuple, Optional


def compute_accessibility_simple(
    sequence_tensor: torch.Tensor,
    deepraccess_model,
    atg_position: int = 0
) -> torch.Tensor:
    """

    
    Args:



        
    Returns:

    """
    return deepraccess_model.compute_atg_window_accessibility(
        sequence_tensor, 
        atg_position=atg_position, 
        discrete=False
    )