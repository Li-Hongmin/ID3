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

from id3.optimizers.cai.sado_incremental import SADOIncrementalOptimizer


def test_incremental_optimization():

    print("\n" + "="*70)

    print("="*70)
    

    seq_length = 300
    num_iterations = 30
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    




    

    optimizer = SADOIncrementalOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    hashes = []
    cai_values = []
    prob_scores = []
    gamma_values = []
    times = []
    

    print("-" * 50)
    
    for i in range(num_iterations):
        start = time.time()
        

        distribution, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8
        )
        
        elapsed = time.time() - start
        times.append(elapsed)
        

        selected = torch.argmax(distribution, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        hashes.append(seq_hash)
        

        cai_values.append(metadata['final_cai'])
        prob_scores.append(metadata.get('prob_score', 0))
        gamma_values.append(metadata.get('gamma', 0))
        

        if i < 5 or (i + 1) % 5 == 0:
            unique_count = len(set(hashes))
            gamma_range = metadata.get('gamma_range', [])




    

    stats = optimizer.get_statistics()
    

    print("\n" + "="*70)

    print("="*70)
    

    unique_hashes = set(hashes)
    hash_counts = Counter(hashes)
    


    

    most_common = hash_counts.most_common(3)
    max_repeat = most_common[0][1] if most_common else 1

    







    

    early_cai = np.mean(cai_values[:10])
    late_cai = np.mean(cai_values[-10:])



    





    






    

    gamma_std_early = np.std(gamma_values[:10]) if len(gamma_values) >= 10 else 0
    gamma_std_late = np.std(gamma_values[-10:]) if len(gamma_values) >= 10 else 0


    if gamma_std_late < gamma_std_early * 0.5:

    else:

    





    

    early_time = np.mean(times[:10])
    late_time = np.mean(times[-10:])


    if late_time < early_time:

    

    print("\n" + "="*70)

    print("="*70)
    
    diversity_score = len(unique_hashes) / num_iterations
    cai_score = sum(1 for c in cai_values if c >= 0.8) / len(cai_values)
    prob_score = np.mean(prob_scores) / 0.167
    




    
    success_count = sum([
        diversity_score >= 0.9,
        cai_score >= 0.7,
        prob_score >= 0.75
    ])
    
    if success_count == 3:

    elif success_count >= 2:

    else:

    
    return {
        'diversity': diversity_score,
        'cai': cai_score,
        'prob': prob_score,
        'time': np.mean(times)
    }


def test_perturbation_effectiveness():
    """测试扰动机制的有效性"""
    print("\n" + "="*70)
    print("扰动机制测试")
    print("="*70)
    

    seq_length = 300
    num_iterations = 20
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOIncrementalOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    print("\n第一轮：正常运行...")
    hashes_round1 = []
    
    for i in range(num_iterations):
        dist, meta = optimizer.optimize(pi_accessibility.clone(), target_cai=0.8)
        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        hashes_round1.append(seq_hash)
    
    unique_round1 = len(set(hashes_round1))
    print(f"  唯一序列: {unique_round1}/{num_iterations}")
    

    print("\n第二轮：继续优化（测试扰动）...")
    hashes_round2 = []
    
    for i in range(num_iterations):
        dist, meta = optimizer.optimize(pi_accessibility.clone(), target_cai=0.8)
        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        hashes_round2.append(seq_hash)
    
    unique_round2 = len(set(hashes_round2))
    print(f"  唯一序列: {unique_round2}/{num_iterations}")
    

    all_hashes = set(hashes_round1) | set(hashes_round2)
    total_unique = len(all_hashes)
    
    print(f"\n扰动效果:")
    print(f"  总唯一序列: {total_unique}/{num_iterations*2}")
    print(f"  跨轮重复: {num_iterations*2 - total_unique}")
    
    if total_unique >= num_iterations * 1.8:
        print("  ✅ 扰动机制非常有效")
    elif total_unique >= num_iterations * 1.5:
        print("  ⚠️ 扰动机制部分有效")
    else:
        print("  ❌ 扰动机制需要改进")


def test_gamma_convergence():

    print("\n" + "="*70)

    print("="*70)
    

    num_experiments = 3
    all_gamma_histories = []
    
    for exp in range(num_experiments):

        

        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        amino_sequence = ''.join(np.random.choice(list(amino_acids), 300))
        
        torch.manual_seed(42 + exp)
        pi_accessibility = torch.rand(300, 6)
        pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
        
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        optimizer = SADOIncrementalOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        

        gamma_history = []
        for i in range(30):
            dist, meta = optimizer.optimize(pi_accessibility.clone(), target_cai=0.8)
            gamma_history.append(meta.get('gamma', 0))
        
        all_gamma_histories.append(gamma_history)
        

        final_gamma = gamma_history[-1]
        convergence_window = gamma_history[-10:]
        convergence_std = np.std(convergence_window)
        


        
        if convergence_std < 0.05:

        elif convergence_std < 0.1:

        else:

    


    all_final_gammas = [hist[-1] for hist in all_gamma_histories]
    consensus_gamma = np.mean(all_final_gammas)
    consensus_std = np.std(all_final_gammas)
    

    
    if consensus_std < 0.1:

    else:



if __name__ == '__main__':
    print("\n" + "="*80)

    print("="*80)
    

    results = test_incremental_optimization()
    

    test_perturbation_effectiveness()
    

    test_gamma_convergence()
    
    print("\n" + "="*80)

    print("="*80)
    






