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
from collections import Counter

from id3.optimizers.cai.sado import SADOOptimizer


def compute_sequence_hash(distribution: torch.Tensor) -> str:

    selected = torch.argmax(distribution, dim=1)
    return hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()


def compute_probability_score(distribution: torch.Tensor, pi_accessibility: torch.Tensor, 
                             codon_choices) -> float:
    """计算序列的概率得分"""
    selected = torch.argmax(distribution, dim=1)
    
    log_prob_sum = 0.0
    for pos in range(len(selected)):
        codon_idx = selected[pos].item()
        if pos < len(codon_choices) and codon_choices[pos]:

            for choice in codon_choices[pos]:
                if choice['local_index'] == codon_idx:
                    orig_idx = choice['original_local_index']
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        if prob > 0:
                            log_prob_sum += np.log(prob)
                    break
    

    return np.exp(log_prob_sum / len(selected))


def test_quick_diversity():

    print("\n" + "="*70)

    print("="*70)
    

    seq_length = 300
    num_iterations = 20
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
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
            use_difference_driven=True,
            gamma=0.3
        )
        

        seq_hash = compute_sequence_hash(dist)
        hashes.append(seq_hash)
        

        cai_values.append(meta['final_cai'])
        

        prob_score = compute_probability_score(dist, pi_accessibility, optimizer.codon_choices)
        prob_scores.append(prob_score)
        
        if i < 5 or (i+1) % 5 == 0:
            unique = len(set(hashes))


    

    unique_hashes = len(set(hashes))



    





    


    avg_prob = np.mean(prob_scores)




    


    if unique_hashes == num_iterations:

    elif unique_hashes >= num_iterations * 0.8:

    elif unique_hashes >= num_iterations * 0.5:

    else:

    
    if avg_prob > random_prob * 0.9:

    else:

    
    return unique_hashes, avg_prob


def test_performance():
    """测试性能（3000 AA序列）"""
    print("\n" + "="*70)
    print("性能测试（3000 AA）")
    print("="*70)
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), 3000))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(3000, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    times = []
    for i in range(5):
        start = time.time()
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True,
            gamma=0.3
        )
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"  运行{i+1}: {elapsed:.3f}秒, CAI={meta['final_cai']:.4f}")
    
    avg_time = np.mean(times)
    print(f"\n平均时间: {avg_time:.3f}秒")
    
    if avg_time < 0.5:
        print("✅ 性能优秀！")
    elif avg_time < 2.0:
        print("✅ 性能良好")
    else:
        print("⚠️ 性能一般")
    
    return avg_time


def test_100_iterations():

    print("\n" + "="*70)

    print("="*70)
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), 3000))
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(3000, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    
    hashes = []
    cai_values = []
    times = []
    
    start_total = time.time()
    
    for i in range(100):
        start = time.time()
        dist, meta = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,
            use_difference_driven=True,
            gamma=0.3
        )
        elapsed = time.time() - start
        
        seq_hash = compute_sequence_hash(dist)
        hashes.append(seq_hash)
        cai_values.append(meta['final_cai'])
        times.append(elapsed)
        
        if (i+1) % 20 == 0:
            unique = len(set(hashes))
            avg_time = np.mean(times)


    
    total_time = time.time() - start_total
    

    unique_count = len(set(hashes))
    





    

    if unique_count >= 95:

    elif unique_count >= 80:

    elif unique_count >= 50:

    else:

    
    if total_time < 30:

    elif total_time < 60:

    else:

    
    return unique_count, total_time


if __name__ == '__main__':
    print("\n" + "="*80)

    print("="*80)
    

    unique_20, prob_score = test_quick_diversity()
    

    avg_time_3000 = test_performance()
    

    unique_100, total_time_100 = test_100_iterations()
    

    print("\n" + "="*80)

    print("="*80)
    






    
    if unique_20 >= 18 and unique_100 >= 80 and avg_time_3000 < 2.0:

    else:

