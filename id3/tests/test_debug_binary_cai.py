#!/usr/bin/env python3
"""

"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.utils.constants import amino_acids_to_codons
from id3.optimizers.cai import BinarySearchCAIOptimizer


seq_length = 100
amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
test_sequence = ''.join(np.random.choice(amino_acids, seq_length))

print(f"测试序列长度: {seq_length}")


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
optimizer = BinarySearchCAIOptimizer(
    species='ecoli_bl21de3',
    device=device,
    amino_acid_sequence=test_sequence
)


pi_accessibility = torch.rand(seq_length, 6, device=device)
for pos in range(seq_length):
    if pi_accessibility[pos].sum() > 0:
        pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()


valid_mask = torch.zeros(seq_length, 6, dtype=torch.bool, device=device)
for pos, aa in enumerate(test_sequence):
    if aa in amino_acids_to_codons:
        num_codons = len(amino_acids_to_codons[aa])
        for i in range(min(num_codons, 6)):
            valid_mask[pos, i] = True

print(f"最大可达CAI: {optimizer.max_achievable_cai:.4f}")


for target_cai in [0.9, 0.8, 0.7, 0.6, 0.5]:
    try:
        result, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            amino_acid_sequence=test_sequence,
            valid_codon_mask=valid_mask
        )
        
        print(f"\nTarget CAI: {target_cai}")
        print(f"  Achieved CAI: {metadata['final_cai']:.4f}")
        print(f"  Optimal gamma: {metadata['optimal_gamma']:.4f}")
        print(f"  Constraint satisfied: {metadata['constraint_satisfied']}")
    except Exception as e:
        print(f"\nTarget CAI: {target_cai}")
        print(f"  Error: {e}")