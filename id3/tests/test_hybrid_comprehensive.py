#!/usr/bin/env python3
"""
Comprehensive test for hybrid_bs_sado optimizer







"""

import torch
import numpy as np
import time
import hashlib
from collections import defaultdict
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons
from id3.experiments.utils.data_loader import ProteinDataLoader


def generate_test_inputs(amino_acid_sequence):

    num_positions = len(amino_acid_sequence)
    max_codons = 6
    
    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long)
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            n_codons = len(amino_acids_to_codons[aa])
            valid_codon_mask[pos, :n_codons] = True
            codon_indices[pos, :n_codons] = torch.arange(n_codons)
    
    return valid_codon_mask, codon_indices


def test_sequence_diversity(method, sequence, iterations=20):
    """测试序列多样性"""
    print(f"\n测试 {method} - 序列长度: {len(sequence)}")
    print("-" * 60)
    
    valid_codon_mask, codon_indices = generate_test_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method=method,
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    
    sequences = []
    switch_events = []
    gamma_values = []
    times = []
    
    torch.manual_seed(42)
    
    for i in range(iterations):

        num_positions = len(sequence)
        pi_accessibility = torch.rand(num_positions, 6)
        pi_accessibility = pi_accessibility * valid_codon_mask.float()
        

        for pos in range(num_positions):
            if valid_codon_mask[pos].any():
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        if i > iterations // 2:
            pi_accessibility = pi_accessibility ** 3
            for pos in range(num_positions):
                if valid_codon_mask[pos].any():
                    pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        
        start_time = time.time()
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        elapsed_time = time.time() - start_time
        times.append(elapsed_time)
        

        if discrete_dist.dim() > 1:
            indices = discrete_dist.argmax(dim=-1)
        else:
            indices = discrete_dist
        
        seq_hash = hashlib.md5(indices.cpu().numpy().tobytes()).hexdigest()
        sequences.append(seq_hash)
        

        if method == 'hybrid_bs_sado':
            if metadata.get('switched_to_sado', False):
                switch_events.append(i)
                if 'bs_gamma' in metadata:
                    gamma_values.append(metadata['bs_gamma'])
        

        if 'optimal_gamma' in metadata:
            gamma_values.append(metadata['optimal_gamma'])
    

    unique_sequences = len(set(sequences))
    repetition_rate = 1 - (unique_sequences / len(sequences))
    avg_time = np.mean(times) * 1000
    
    print(f"  生成序列: {len(sequences)}")
    print(f"  唯一序列: {unique_sequences}")
    print(f"  重复率: {repetition_rate:.1%}")
    print(f"  平均时间: {avg_time:.2f}ms")
    
    if method == 'hybrid_bs_sado' and switch_events:
        print(f"  切换到SADO: {len(switch_events)}次 (在迭代: {switch_events})")
        if gamma_values:
            print(f"  传递的gamma值: {[f'{g:.3f}' for g in gamma_values[:5]]}...")
    
    if hasattr(operator.optimizer, 'get_statistics'):
        stats = operator.optimizer.get_statistics()
        if method == 'hybrid_bs_sado':
            print(f"  BS调用: {stats.get('bs_calls', 0)}")
            print(f"  SADO调用: {stats.get('sado_calls', 0)}")
            print(f"  切换率: {stats.get('switch_rate', 0):.1%}")
    
    return {
        'unique_sequences': unique_sequences,
        'repetition_rate': repetition_rate,
        'avg_time': avg_time,
        'switch_events': len(switch_events) if method == 'hybrid_bs_sado' else 0
    }


def test_gamma_transmission():

    print("\n" + "="*80)

    print("="*80)
    

    valid_codon_mask, codon_indices = generate_test_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method='hybrid_bs_sado',
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    

    num_positions = len(sequence)
    pi_accessibility = torch.ones(num_positions, 6) * 0.1

    

    pi_accessibility = pi_accessibility * valid_codon_mask.float()
    for pos in range(num_positions):
        if valid_codon_mask[pos].any():
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    

    
    for i in range(5):
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        
        switched = metadata.get('switched_to_sado', False)
        method_used = "SADO" if switched else "Binary Search"
        

        
        if 'optimal_gamma' in metadata:
            print(f"  BS optimal_gamma: {metadata['optimal_gamma']:.4f}")
        
        if switched and 'bs_gamma' in metadata:




def test_real_proteins():
    """测试真实蛋白质序列"""
    print("\n" + "="*80)
    print("真实蛋白质序列测试")
    print("="*80)
    
    loader = ProteinDataLoader()
    

    test_proteins = {
        'O15263': None,
        'P0CG48': None,
    }
    
    for protein_id in test_proteins:
        try:
            sequence = loader.load_protein_sequence(protein_id)
            test_proteins[protein_id] = sequence[:100]
        except:

            if protein_id == 'O15263':
                test_proteins[protein_id] = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAI"[:50]
            else:
                test_proteins[protein_id] = "MKWVTFISLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS"[:50]
    
    results = {}
    methods = ['binary_search', 'sado', 'hybrid_bs_sado']
    
    for protein_id, sequence in test_proteins.items():
        if sequence:
            print(f"\n\n蛋白质 {protein_id} (长度={len(sequence)})")
            print("="*60)
            
            for method in methods:
                result = test_sequence_diversity(method, sequence, iterations=10)
                results[f"{protein_id}_{method}"] = result
    

    print("\n" + "="*80)
    print("性能比较总结")
    print("="*80)
    
    for protein_id in test_proteins:
        if test_proteins[protein_id]:
            print(f"\n{protein_id}:")
            print(f"  {'方法':<15} {'唯一序列':<10} {'重复率':<10} {'平均时间':<10} {'切换次数':<10}")
            print("  " + "-"*60)
            
            for method in methods:
                key = f"{protein_id}_{method}"
                if key in results:
                    r = results[key]
                    switch_str = f"{r['switch_events']}" if method == 'hybrid_bs_sado' else "N/A"
                    print(f"  {method:<15} {r['unique_sequences']:<10} {r['repetition_rate']:<10.1%} {r['avg_time']:<10.2f}ms {switch_str:<10}")


def test_extreme_cases():

    print("\n" + "="*80)

    print("="*80)
    


    sequence = "MWMWMW"
    valid_codon_mask, codon_indices = generate_test_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method='hybrid_bs_sado',
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    

    num_positions = len(sequence)
    pi_accessibility = torch.zeros(num_positions, 6)

    
    for i in range(3):
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        
        switched = metadata.get('switched_to_sado', False)

    


    sequence = "LLLLLLLLLL"
    result = test_sequence_diversity('hybrid_bs_sado', sequence, iterations=15)
    




if __name__ == "__main__":
    print("="*80)
    print("Hybrid BS-SADO Comprehensive Test")
    print("="*80)
    

    test_gamma_transmission()
    

    test_real_proteins()
    

    test_extreme_cases()
    
    print("\n" + "="*80)

    print("="*80)