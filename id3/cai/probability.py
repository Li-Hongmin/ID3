"""
概率转换模块
只负责概率分布转换
"""

import torch
import numpy as np
from typing import Dict, List


def rna_to_codon_probabilities(rna_probs: torch.Tensor,
                               amino_sequence: str,
                               codon_table: Dict) -> Dict[int, Dict[str, float]]:
    """
    将RNA概率转换为密码子概率
    
    Args:
        rna_probs: RNA概率张量 [seq_len, 4]
        amino_sequence: 氨基酸序列
        codon_table: 氨基酸到密码子映射
        
    Returns:
        密码子概率分布字典
    """
    rna_indices = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    codon_probabilities = {}
    
    for aa_pos, amino_acid in enumerate(amino_sequence):
        if amino_acid not in codon_table:
            continue
        
        rna_start = aa_pos * 3
        rna_end = min(rna_start + 3, rna_probs.shape[0])
        
        if rna_end - rna_start < 3:
            continue
        
        codon_rna_probs = rna_probs[rna_start:rna_end]
        
        # 转换为numpy便于计算
        if isinstance(codon_rna_probs, torch.Tensor):
            codon_rna_probs = codon_rna_probs.detach().cpu().numpy()
        
        available_codons = codon_table[amino_acid]
        codon_probs = {}
        
        # 计算每个密码子的联合概率
        for codon in available_codons:
            prob = 1.0
            for pos, nuc in enumerate(codon):
                nuc_idx = rna_indices.get(nuc, 3)
                prob *= codon_rna_probs[pos, nuc_idx]
            codon_probs[codon] = float(prob)
        
        # 归一化
        total = sum(codon_probs.values())
        if total > 1e-10:
            codon_probs = {c: p/total for c, p in codon_probs.items()}
        
        codon_probabilities[aa_pos] = codon_probs
    
    return codon_probabilities


def generate_uniform_codon_probs(amino_sequence: str,
                                 codon_table: Dict) -> Dict[int, Dict[str, float]]:
    """
    生成均匀的密码子概率分布
    
    Args:
        amino_sequence: 氨基酸序列
        codon_table: 氨基酸到密码子映射
        
    Returns:
        均匀密码子概率分布
    """
    codon_probs = {}
    
    for pos, aa in enumerate(amino_sequence):
        if aa in codon_table:
            codons = codon_table[aa]
            uniform_prob = 1.0 / len(codons)
            codon_probs[pos] = {c: uniform_prob for c in codons}
    
    return codon_probs


def discretize_codon_probs(codon_probs: Dict[int, Dict[str, float]],
                          amino_sequence: str,
                          codon_table: Dict) -> str:
    """
    离散化密码子概率（选择最高概率）
    
    Args:
        codon_probs: 密码子概率分布
        amino_sequence: 氨基酸序列
        codon_table: 氨基酸到密码子映射
        
    Returns:
        DNA序列
    """
    sequence = []
    
    for pos, aa in enumerate(amino_sequence):
        if pos in codon_probs and codon_probs[pos]:
            # 选择概率最高的密码子
            best_codon = max(codon_probs[pos].items(), key=lambda x: x[1])[0]
            sequence.append(best_codon.replace('U', 'T'))
        elif aa in codon_table:
            # 备用：第一个密码子
            sequence.append(codon_table[aa][0].replace('U', 'T'))
    
    return ''.join(sequence)