#!/usr/bin/env python3
"""
Sequence Processing Utility Module

Handles sequence conversion, encoding, and processing
"""

import torch
from typing import List, Tuple, Optional


class SequenceProcessor:
    """Sequence processor"""

    NUCLEOTIDES = ['A', 'C', 'G', 'U']

    def __init__(self, device: Optional[torch.device] = None):
        """
        Initialize processor

        Args:
            device: Computation device
        """
        if device is None:
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = device

    def sequence_to_probs(self, sequence: str) -> torch.Tensor:
        """
        Convert nucleotide sequence to probability tensor

        Args:
            sequence: Nucleotide sequence string

        Returns:
            Probability tensor [L, 4]
        """
        probs = []
        for nt in sequence.upper().replace('T', 'U'):
            prob = torch.zeros(4, dtype=torch.float32, device=self.device)
            if nt in self.NUCLEOTIDES:
                prob[self.NUCLEOTIDES.index(nt)] = 1.0
            else:
                # Unknown nucleotides use uniform distribution
                prob[:] = 0.25
            probs.append(prob)
        return torch.stack(probs)

    def probs_to_sequence(self, probs: torch.Tensor) -> str:
        """
        Convert probability tensor to nucleotide sequence

        Args:
            probs: Probability tensor [L, 4]

        Returns:
            Nucleotide sequence string
        """
        # Get nucleotide with maximum probability
        indices = probs.argmax(dim=-1)
        sequence = ''.join([self.NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
        return sequence

    def prepare_utr_probs(
        self,
        utr5: str,
        utr3: str
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Prepare UTR probability tensors

        Args:
            utr5: 5'UTR sequence
            utr3: 3'UTR sequence

        Returns:
            (utr5_probs, utr3_probs)
        """
        utr5_probs = self.sequence_to_probs(utr5)
        utr3_probs = self.sequence_to_probs(utr3)
        return utr5_probs, utr3_probs

    def concat_mrna_probs(
        self,
        utr5_probs: torch.Tensor,
        cds_probs: torch.Tensor,
        utr3_probs: torch.Tensor
    ) -> torch.Tensor:
        """
        Concatenate complete mRNA probabilities

        Args:
            utr5_probs: 5'UTR probabilities [L5, 4]
            cds_probs: CDS probabilities [L_cds, 4]
            utr3_probs: 3'UTR probabilities [L3, 4]

        Returns:
            Complete mRNA probabilities [1, L_total, 4]
        """
        full_probs = torch.cat([utr5_probs, cds_probs, utr3_probs], dim=0)
        return full_probs.unsqueeze(0)  # Add batch dimension
