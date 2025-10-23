#!/usr/bin/env python3
"""

"""
import torch
import time
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons

def benchmark_cai_enhancement(sequence_length=1000, num_runs=10):

    

    amino_acid_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK" * (sequence_length // 238 + 1)
    amino_acid_sequence = amino_acid_sequence[:sequence_length]
    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    
    operator = CAIEnhancementOperator(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_acid_sequence
    )
    

    max_codons = 6
    num_positions = len(amino_acid_sequence)
    

    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool, device=device)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long, device=device)
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            codons = amino_acids_to_codons[aa]
            for i, codon in enumerate(codons):
                if i < max_codons:
                    valid_codon_mask[pos, i] = True

                    codon_indices[pos, i] = i
    

    pi_accessibility = torch.rand(num_positions, max_codons, device=device)

    pi_accessibility = pi_accessibility * valid_codon_mask.float()
    pi_accessibility = pi_accessibility / (pi_accessibility.sum(dim=-1, keepdim=True) + 1e-10)
    

    target_cai = 0.7
    




    


    for _ in range(2):
        _, _ = operator.apply_cai_enhancement(
            pi_accessibility, amino_acid_sequence, 
            valid_codon_mask, codon_indices, target_cai
        )
    


    times = []
    
    for run in range(num_runs):
        start_time = time.time()
        
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility, amino_acid_sequence, 
            valid_codon_mask, codon_indices, target_cai
        )
        
        elapsed = time.time() - start_time
        times.append(elapsed)
        
        if run == 0:





    

    avg_time = sum(times) / len(times)
    min_time = min(times)
    max_time = max(times)
    




    



    speedup = estimated_original_time / avg_time
    

    
    return avg_time

def compare_sequence_lengths():
    """比较不同序列长度的性能"""
    lengths = [100, 500, 1000, 2000]
    
    print("\n" + "="*50)
    print("不同序列长度性能比较")
    print("="*50)
    
    results = []
    for length in lengths:
        print(f"\n序列长度: {length}")
        avg_time = benchmark_cai_enhancement(sequence_length=length, num_runs=5)
        results.append((length, avg_time))
    
    print("\n" + "="*50)
    print("总结:")
    print("-"*50)
    print("长度\t时间(ms)\t每个位置(μs)")
    print("-"*50)
    for length, time_sec in results:
        time_ms = time_sec * 1000
        time_per_pos = (time_sec * 1e6) / length
        print(f"{length}\t{time_ms:.2f}\t\t{time_per_pos:.2f}")

if __name__ == "__main__":
    print("CAI Enhancement性能优化测试")
    print("="*50)
    

    benchmark_cai_enhancement()
    

    compare_sequence_lengths()