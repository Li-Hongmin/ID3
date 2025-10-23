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

from id3.optimizers.cai.sado_v2 import SADOv2Optimizer


def test_sado_v2():

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
    




    

    optimizer = SADOv2Optimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    hashes = []
    cai_values = []
    prob_scores = []
    diversity_scores = []
    pareto_sizes = []
    weights_history = []
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
        diversity_scores.append(metadata.get('diversity_score', 0))
        pareto_sizes.append(metadata.get('pareto_size', 0))
        weights_history.append(metadata.get('weights', (0, 0, 0)))
        

        if (i + 1) <= 5 or (i + 1) % 5 == 0:
            unique_count = len(set(hashes))
            weights = metadata.get('weights', (0, 0, 0))





    

    print("\n" + "="*70)

    print("="*70)
    

    unique_hashes = set(hashes)
    hash_counts = Counter(hashes)
    


    

    most_common = hash_counts.most_common(3)
    if most_common[0][1] > 1:

    







    

    valid_probs = [p for p in prob_scores if p > 0]
    if valid_probs:






    




    


    final_weights = weights_history[-1] if weights_history else (0, 0, 0)


    

    prob_weights = [w[0] for w in weights_history]
    cai_weights = [w[1] for w in weights_history]


    




    

    print("\n" + "="*70)

    print("="*70)
    
    diversity_score = len(unique_hashes) / num_iterations
    cai_score = sum(1 for c in cai_values if c >= 0.8) / len(cai_values)
    prob_score = np.mean(valid_probs) / 0.167 if valid_probs else 0
    




    
    overall_success = diversity_score >= 0.8 and cai_score >= 0.5 and prob_score >= 0.7
    
    if overall_success:

    else:

    
    return {
        'diversity': diversity_score,
        'cai': cai_score,
        'prob': prob_score,
        'time': np.mean(times)
    }


def compare_with_v1():
    """与SADO v1比较"""
    print("\n" + "="*70)
    print("SADO版本对比测试")
    print("="*70)
    

    seq_length = 300
    num_iterations = 10
    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(np.random.choice(list(amino_acids), seq_length))
    
    torch.manual_seed(123)
    pi_accessibility = torch.rand(seq_length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    print("\n测试SADO 2.0...")
    optimizer_v2 = SADOv2Optimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    v2_hashes = []
    v2_cais = []
    v2_times = []
    
    for i in range(num_iterations):
        start = time.time()
        dist, meta = optimizer_v2.optimize(pi_accessibility.clone(), target_cai=0.8)
        elapsed = time.time() - start
        
        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        
        v2_hashes.append(seq_hash)
        v2_cais.append(meta['final_cai'])
        v2_times.append(elapsed)
    

    print("\n测试SADO 1.0...")
    from id3.optimizers.cai.sado import SADOOptimizer
    
    optimizer_v1 = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    v1_hashes = []
    v1_cais = []
    v1_times = []
    
    for i in range(num_iterations):
        start = time.time()
        dist, meta = optimizer_v1.optimize(
            pi_accessibility.clone(), 
            target_cai=0.8,
            use_difference_driven=True
        )
        elapsed = time.time() - start
        
        selected = torch.argmax(dist, dim=1)
        seq_hash = hashlib.md5(selected.cpu().numpy().tobytes()).hexdigest()
        
        v1_hashes.append(seq_hash)
        v1_cais.append(meta['final_cai'])
        v1_times.append(elapsed)
    

    print("\n" + "="*70)
    print("对比结果")
    print("="*70)
    
    print(f"\n指标           SADO 1.0        SADO 2.0        改进")
    print("-" * 60)
    
    v1_diversity = len(set(v1_hashes)) / num_iterations
    v2_diversity = len(set(v2_hashes)) / num_iterations
    print(f"多样性         {v1_diversity:.2%}           {v2_diversity:.2%}           "
          f"{'+' if v2_diversity > v1_diversity else ''}{(v2_diversity - v1_diversity):.2%}")
    
    v1_avg_cai = np.mean(v1_cais)
    v2_avg_cai = np.mean(v2_cais)
    print(f"平均CAI        {v1_avg_cai:.4f}          {v2_avg_cai:.4f}          "
          f"{'+' if v2_avg_cai > v1_avg_cai else ''}{(v2_avg_cai - v1_avg_cai):.4f}")
    
    v1_avg_time = np.mean(v1_times)
    v2_avg_time = np.mean(v2_times)
    print(f"平均时间       {v1_avg_time:.3f}s          {v2_avg_time:.3f}s          "
          f"{'+' if v2_avg_time < v1_avg_time else ''}{(v2_avg_time - v1_avg_time):.3f}s")
    
    print("\n结论:")
    improvements = 0
    if v2_diversity > v1_diversity:
        print("  ✅ SADO 2.0在多样性方面有改进")
        improvements += 1
    if v2_avg_cai > v1_avg_cai:
        print("  ✅ SADO 2.0在CAI优化方面有改进")
        improvements += 1
    if v2_avg_time < v1_avg_time * 1.5:
        print("  ✅ SADO 2.0的性能可接受")
        improvements += 1
    
    if improvements >= 2:
        print("\n🎉 SADO 2.0总体优于SADO 1.0")
    else:
        print("\n⚠️ SADO 2.0需要进一步优化")


if __name__ == '__main__':

    results = test_sado_v2()
    

    compare_with_v1()
    
    print("\n" + "="*70)
    print("测试完成")
    print("="*70)