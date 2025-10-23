#!/usr/bin/env python3
"""







"""

import torch
import numpy as np
import time
import hashlib
from collections import defaultdict
from typing import Dict, List
import pandas as pd
from tabulate import tabulate

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons
from id3.experiments.utils.data_loader import ProteinDataLoader


def prepare_inputs(sequence: str):

    num_positions = len(sequence)
    max_codons = 6
    
    valid_codon_mask = torch.zeros(num_positions, max_codons, dtype=torch.bool)
    codon_indices = torch.zeros(num_positions, max_codons, dtype=torch.long)
    
    for pos, aa in enumerate(sequence):
        if aa in amino_acids_to_codons:
            n_codons = len(amino_acids_to_codons[aa])
            valid_codon_mask[pos, :n_codons] = True
            codon_indices[pos, :n_codons] = torch.arange(n_codons)
    
    return valid_codon_mask, codon_indices


def test_multiple_iterations(method: str, sequence: str, iterations: int = 50):
    """测试多次迭代的性能"""
    
    valid_codon_mask, codon_indices = prepare_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method=method,
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    
    results = {
        'times': [],
        'cai_values': [],
        'unique_sequences': set(),
        'iteration_types': []
    }
    
    torch.manual_seed(42)
    
    for i in range(iterations):

        num_positions = len(sequence)
        pi_accessibility = torch.rand(num_positions, 6)
        

        pi_accessibility = pi_accessibility * valid_codon_mask.float()
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
        

        results['times'].append(elapsed_time * 1000)
        results['cai_values'].append(metadata.get('final_cai', 0))
        

        if discrete_dist.dim() > 1:
            indices = discrete_dist.argmax(dim=-1)
        else:
            indices = discrete_dist
        seq_hash = hashlib.md5(indices.cpu().numpy().tobytes()).hexdigest()
        results['unique_sequences'].add(seq_hash)
        

        if method == 'incremental':
            iteration_type = metadata.get('method', 'unknown')
            results['iteration_types'].append(iteration_type)
    
    return results


def compare_methods():

    
    print("="*80)

    print("="*80)
    

    test_sequences = {



    }
    
    methods = ['binary_search', 'sado', 'hybrid_bs_sado', 'incremental']
    
    all_results = []
    
    for seq_name, sequence in test_sequences.items():

        print("-"*60)
        
        for method in methods:

            
            try:
                results = test_multiple_iterations(method, sequence, iterations=30)
                

                avg_time = np.mean(results['times'])
                std_time = np.std(results['times'])
                unique_count = len(results['unique_sequences'])
                repetition_rate = 1 - (unique_count / len(results['times']))
                avg_cai = np.mean(results['cai_values'])
                

                if method == 'incremental' and results['iteration_types']:
                    bs_count = sum(1 for t in results['iteration_types'] if 'bs' in t)
                    inc_count = sum(1 for t in results['iteration_types'] if 'inc' in t)

                

                if len(results['times']) > 10:
                    first_time = results['times'][0]
                    avg_subsequent = np.mean(results['times'][1:])
                    speedup = first_time / avg_subsequent if avg_subsequent > 0 else 1.0
                    


                




                

                all_results.append({
                    'sequence': seq_name,
                    'method': method,
                    'avg_time_ms': avg_time,
                    'std_time_ms': std_time,
                    'unique_sequences': unique_count,
                    'repetition_rate': repetition_rate * 100,
                    'avg_cai': avg_cai,
                    'first_time_ms': results['times'][0] if results['times'] else 0,
                    'subsequent_avg_ms': np.mean(results['times'][1:]) if len(results['times']) > 1 else 0
                })
                
            except Exception as e:

    

    print("\n" + "="*80)

    print("="*80)
    
    df = pd.DataFrame(all_results)
    
    for seq_name in test_sequences.keys():
        print(f"\n{seq_name}:")
        seq_df = df[df['sequence'] == seq_name]
        

        table_data = []
        for _, row in seq_df.iterrows():
            table_data.append([
                row['method'],
                f"{row['avg_time_ms']:.2f}",
                f"{row['first_time_ms']:.2f}",
                f"{row['subsequent_avg_ms']:.2f}",
                f"{row['unique_sequences']}",
                f"{row['repetition_rate']:.1f}%",
                f"{row['avg_cai']:.4f}"
            ])
        

        print(tabulate(table_data, headers=headers, tablefmt='grid'))
    

    print("\n" + "="*80)

    print("="*80)
    
    incremental_df = df[df['method'] == 'incremental']
    if not incremental_df.empty:



        

        bs_df = df[df['method'] == 'binary_search']
        for seq_name in test_sequences.keys():
            inc_row = incremental_df[incremental_df['sequence'] == seq_name]
            bs_row = bs_df[bs_df['sequence'] == seq_name]
            
            if not inc_row.empty and not bs_row.empty:
                speedup = bs_row.iloc[0]['avg_time_ms'] / inc_row.iloc[0]['subsequent_avg_ms']
                print(f"\n{seq_name}:")

    
    return df


def test_convergence():
    """测试收敛性和稳定性"""
    print("\n" + "="*80)
    print("收敛性测试")
    print("="*80)
    
    sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTT"
    valid_codon_mask, codon_indices = prepare_inputs(sequence)
    

    torch.manual_seed(42)
    num_positions = len(sequence)
    pi_accessibility = torch.rand(num_positions, 6)
    pi_accessibility = pi_accessibility * valid_codon_mask.float()
    for pos in range(num_positions):
        if valid_codon_mask[pos].any():
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    
    operator = CAIEnhancementOperator(
        method='incremental',
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    
    print("\n使用相同输入连续迭代10次：")
    for i in range(10):
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        
        print(f"  迭代 {i+1}: CAI={metadata.get('final_cai', 0):.4f}, "
              f"方法={metadata.get('method', 'unknown')}, "
              f"变化数={metadata.get('incremental_changes', 0)}")


if __name__ == "__main__":
    print("开始测试增量式CAI优化器...")
    

    results_df = compare_methods()
    

    test_convergence()
    

    results_df.to_csv('incremental_optimizer_results.csv', index=False)
    print("\n结果已保存到 incremental_optimizer_results.csv")
    
    print("\n" + "="*80)
    print("测试完成！")
    print("="*80)