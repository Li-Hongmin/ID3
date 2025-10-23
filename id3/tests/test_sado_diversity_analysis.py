#!/usr/bin/env python3
"""


"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import hashlib
from collections import Counter

from id3.optimizers.cai.sado import SADOOptimizer


def compute_sequence_hash(distribution: torch.Tensor) -> str:


    selected = torch.argmax(distribution, dim=1)

    return hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()


def compute_probability_score(distribution: torch.Tensor, pi_accessibility: torch.Tensor) -> float:
    """计算序列的概率得分"""

    selected = torch.argmax(distribution, dim=1)
    

    log_prob_sum = 0.0
    for pos in range(len(selected)):
        codon_idx = selected[pos].item()
        if codon_idx < pi_accessibility.shape[1]:
            prob = pi_accessibility[pos, codon_idx].item()
            if prob > 0:
                log_prob_sum += np.log(prob)
    

    return np.exp(log_prob_sum / len(selected))


def test_diversity_and_probability():

    
    print("\n" + "="*70)

    print("="*70)
    

    seq_length = 3000
    num_iterations = 100
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    




    

    all_hashes = []
    all_cai_values = []
    all_prob_scores = []
    

    optimizer_with_history = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    print("-" * 50)
    
    for i in range(num_iterations):

        distribution, metadata = optimizer_with_history.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True,
            gamma=0.3
        )
        

        seq_hash = compute_sequence_hash(distribution)
        all_hashes.append(seq_hash)
        

        all_cai_values.append(metadata['final_cai'])
        

        prob_score = compute_probability_score(distribution, pi_accessibility)
        all_prob_scores.append(prob_score)
        

        if (i + 1) % 20 == 0:
            unique_count = len(set(all_hashes))

                  f"CAI={metadata['final_cai']:.4f}, "

    

    print("\n" + "="*70)

    print("="*70)
    

    unique_hashes = set(all_hashes)
    hash_counts = Counter(all_hashes)
    




    

    most_common = hash_counts.most_common(5)

    for hash_val, count in most_common:
        if count > 1:

    







    






    


    first_10_avg = np.mean(all_prob_scores[:10])
    last_10_avg = np.mean(all_prob_scores[-10:])



    

    diversity_stats = optimizer_with_history.get_diversity_stats()




    

    print("\n" + "="*70)

    print("="*70)
    
    if len(unique_hashes) == num_iterations:

    elif len(unique_hashes) >= num_iterations * 0.95:

    elif len(unique_hashes) >= num_iterations * 0.8:

    else:

    
    if np.mean(all_prob_scores) > 0.1:

    else:

    
    return len(unique_hashes), all_prob_scores


def test_without_history_reset():
    """测试每次重置优化器的情况"""
    
    print("\n" + "="*70)
    print("测试每次重置优化器（无历史记忆）")
    print("="*70)
    
    seq_length = 1000
    num_iterations = 20
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    all_hashes = []
    
    print(f"\n运行{num_iterations}次（每次新建优化器）...")
    
    for i in range(num_iterations):

        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        
        distribution, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True
        )
        
        seq_hash = compute_sequence_hash(distribution)
        all_hashes.append(seq_hash)
        
        if (i + 1) % 5 == 0:
            unique_count = len(set(all_hashes))
            print(f"  迭代 {i+1}: 唯一序列={unique_count}/{i+1}")
    
    unique_hashes = set(all_hashes)
    print(f"\n结果: {len(unique_hashes)}/{num_iterations} 唯一序列")
    
    if len(unique_hashes) == 1:
        print("❌ 每次都产生相同序列（确定性算法）")
    else:
        print(f"✅ 产生了{len(unique_hashes)}个不同序列")


if __name__ == '__main__':

    unique_count, prob_scores = test_diversity_and_probability()
    

    test_without_history_reset()
    
    print("\n" + "="*70)
    print("测试完成！")
    print("="*70)