#!/usr/bin/env python3
"""

"""

import torch
import sys
from pathlib import Path

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.utils.constants import amino_acid_token_map, amino_acid_to_codon_matrix, codon_to_rna_matrix, NUCLEOTIDES

def debug_constraint_argmax():

    


    batch_size = 2
    seq_len = 3
    


    

    aa_idx = amino_acid_token_map[amino_acid_sequence[0]]
    valid_codon_indices = amino_acid_to_codon_matrix[aa_idx].nonzero(as_tuple=True)[0]


    

    for idx in valid_codon_indices:
        codon_encoding = codon_to_rna_matrix[idx]


        indices = torch.argmax(codon_encoding, dim=-1)
        seq = ''.join([NUCLEOTIDES[i] for i in indices])

    

    torch.manual_seed(42)
    rna_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)
    

    from id3.utils.constraint_satisfied_argmax import get_constraint_satisfied_argmax
    result = get_constraint_satisfied_argmax(rna_probs, amino_acid_sequence)
    

    for b in range(batch_size):
        indices = torch.argmax(result[b], dim=-1)
        dna_seq = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])

        

        if dna_seq == "ATG":

        else:


if __name__ == "__main__":
    debug_constraint_argmax()