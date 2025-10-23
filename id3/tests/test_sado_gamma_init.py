#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import hashlib
import time

from id3.optimizers.cai.sado import SADOOptimizer


def test_gamma_initialization():

    print("\n" + "="*60)

    print("="*60)
    

    seq_length = 300
    num_iterations = 5
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    gamma_values = [0.0, 0.1, 0.2, 0.3, 0.5]
    
    for gamma in gamma_values:

        

        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        
        hashes = []
        cai_values = []
        prob_scores = []
        
        for i in range(num_iterations):
            dist, meta = optimizer.optimize(
                pi_accessibility=pi_accessibility.clone(),
                target_cai=0.8,
                use_binary_search=False,
                use_difference_driven=True,

            )
            

            selected = torch.argmax(dist, dim=1)
            seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
            hashes.append(seq_hash)
            

            cai_values.append(meta['final_cai'])
            

            log_prob = 0.0
            for pos in range(seq_length):
                codon_idx = selected[pos].item()
                if pos < len(optimizer.codon_choices) and optimizer.codon_choices[pos]:
                    for choice in optimizer.codon_choices[pos]:
                        if choice['local_index'] == codon_idx:
                            orig_idx = choice['original_local_index']
                            if orig_idx < pi_accessibility.shape[1]:
                                prob = pi_accessibility[pos, orig_idx].item()
                                if prob > 0:
                                    log_prob += np.log(prob)
                            break
            prob_score = np.exp(log_prob / seq_length)
            prob_scores.append(prob_score)
        

        unique_count = len(set(hashes))
        avg_cai = np.mean(cai_values)
        avg_prob = np.mean(prob_scores)
        






def test_without_difference_driven():
    """测试不使用差异驱动优化"""
    print("\n" + "="*60)
    print("测试标准贪心优化（不使用差异驱动）")
    print("="*60)
    

    seq_length = 300
    num_iterations = 10
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    hashes = []
    cai_values = []
    prob_scores = []
    
    for i in range(num_iterations):
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=False,
            gamma=0.1
        )
        

        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        hashes.append(seq_hash)
        

        cai_values.append(meta['final_cai'])
        

        log_prob = 0.0
        for pos in range(seq_length):
            codon_idx = selected[pos].item()
            if pos < len(optimizer.codon_choices) and optimizer.codon_choices[pos]:
                for choice in optimizer.codon_choices[pos]:
                    if choice['local_index'] == codon_idx:
                        orig_idx = choice['original_local_index']
                        if orig_idx < pi_accessibility.shape[1]:
                            prob = pi_accessibility[pos, orig_idx].item()
                            if prob > 0:
                                log_prob += np.log(prob)
                        break
        prob_score = np.exp(log_prob / seq_length)
        prob_scores.append(prob_score)
        
        print(f"  迭代{i+1:2d}: CAI={meta['final_cai']:.4f}, 概率={prob_score:.6f}")
    

    unique_count = len(set(hashes))
    avg_prob = np.mean(prob_scores)
    
    print(f"\n结果:")
    print(f"  唯一序列: {unique_count}/{num_iterations} ({unique_count*10}%)")
    print(f"  平均CAI: {np.mean(cai_values):.4f}")
    print(f"  平均概率: {avg_prob:.6f}")
    print(f"  概率提升: {avg_prob/0.167:.2f}x")


if __name__ == '__main__':
    test_gamma_initialization()
    test_without_difference_driven()