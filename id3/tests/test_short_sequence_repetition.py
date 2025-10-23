#!/usr/bin/env python3
"""








"""

import torch
import numpy as np
import hashlib
import time
from collections import defaultdict, Counter
from typing import Dict, List
import pandas as pd
from tabulate import tabulate

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


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


def test_repetition_with_convergence(method: str, sequence: str, iterations: int = 100):
    """测试不同收敛程度下的重复率"""
    
    valid_codon_mask, codon_indices = prepare_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method=method,
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    
    results_by_convergence = {}
    

    convergence_levels = [
        ('uniform', 1),
        ('mild', 2),
        ('moderate', 4),
        ('high', 8),
        ('extreme', 16),
    ]
    
    for conv_name, conv_power in convergence_levels:
        sequences_generated = []
        switch_events = 0
        
        torch.manual_seed(42)
        
        for i in range(iterations):

            num_positions = len(sequence)
            pi_accessibility = torch.rand(num_positions, 6)
            

            pi_accessibility = pi_accessibility ** conv_power
            

            pi_accessibility = pi_accessibility * valid_codon_mask.float()
            for pos in range(num_positions):
                if valid_codon_mask[pos].any():
                    pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
            
            try:
                discrete_dist, metadata = operator.apply_cai_enhancement(
                    pi_accessibility=pi_accessibility,
                    amino_acid_sequence=sequence,
                    valid_codon_mask=valid_codon_mask,
                    codon_indices=codon_indices,
                    target_cai=0.8
                )
                

                if discrete_dist.dim() > 1:
                    indices = discrete_dist.argmax(dim=-1)
                else:
                    indices = discrete_dist
                
                seq_hash = hashlib.md5(indices.cpu().numpy().tobytes()).hexdigest()
                sequences_generated.append(seq_hash)
                

                if method == 'hybrid_bs_sado' and metadata.get('switched_to_sado', False):
                    switch_events += 1
                    
            except Exception as e:
                continue
        

        unique_sequences = len(set(sequences_generated))
        repetition_rate = 1 - (unique_sequences / len(sequences_generated)) if sequences_generated else 0
        

        seq_counter = Counter(sequences_generated)
        most_common = seq_counter.most_common(1)[0] if seq_counter else (None, 0)
        
        results_by_convergence[conv_name] = {
            'unique_sequences': unique_sequences,
            'total_sequences': len(sequences_generated),
            'repetition_rate': repetition_rate * 100,
            'most_repeated_count': most_common[1],
            'switch_events': switch_events
        }
    
    return results_by_convergence


def test_different_sequence_types():

    
    print("="*80)

    print("="*80)
    

    test_sequences = {








    }
    
    methods = ['binary_search', 'sado', 'hybrid_bs_sado']
    all_results = []
    
    for seq_name, sequence in test_sequences.items():


        print("-"*60)
        
        for method in methods:

            convergence_results = test_repetition_with_convergence(method, sequence, iterations=50)
            

            for conv_level, result in convergence_results.items():



                      end="")
                if method == 'hybrid_bs_sado' and result['switch_events'] > 0:

                else:
                    print()
            

            for conv_level, result in convergence_results.items():
                all_results.append({
                    'sequence_type': seq_name,
                    'sequence_length': len(sequence),
                    'method': method,
                    'convergence': conv_level,
                    'unique_sequences': result['unique_sequences'],
                    'repetition_rate': result['repetition_rate'],
                    'switch_events': result['switch_events'] if method == 'hybrid_bs_sado' else 0
                })
    
    return all_results


def create_summary_table(results: List[Dict]):
    """创建汇总表格"""
    
    df = pd.DataFrame(results)
    
    print("\n" + "="*80)
    print("重复率汇总（按收敛程度）")
    print("="*80)
    

    for convergence in ['uniform', 'mild', 'moderate', 'high', 'extreme']:
        print(f"\n收敛程度: {convergence}")
        print("-"*70)
        
        conv_df = df[df['convergence'] == convergence]
        

        summary_data = []
        for method in ['binary_search', 'sado', 'hybrid_bs_sado']:
            method_df = conv_df[conv_df['method'] == method]
            
            avg_repetition = method_df['repetition_rate'].mean()
            max_repetition = method_df['repetition_rate'].max()
            min_unique = method_df['unique_sequences'].min()
            

            if method == 'hybrid_bs_sado':
                total_switches = method_df['switch_events'].sum()
                switch_rate = (method_df['switch_events'] > 0).mean() * 100
                summary_data.append([
                    method,
                    f"{avg_repetition:.1f}%",
                    f"{max_repetition:.1f}%",
                    f"{min_unique}",
                    f"{total_switches} ({switch_rate:.0f}%)"
                ])
            else:
                summary_data.append([
                    method,
                    f"{avg_repetition:.1f}%",
                    f"{max_repetition:.1f}%",
                    f"{min_unique}",
                    "N/A"
                ])
        
        headers = ['方法', '平均重复率', '最大重复率', '最少唯一序列', '切换事件']
        print(tabulate(summary_data, headers=headers, tablefmt='grid'))
    

    print("\n" + "="*80)
    print("特殊序列重复率分析")
    print("="*80)
    
    special_sequences = ['单密码子(M和W)', '低多样性(单一AA)', '重复模式']
    
    for seq_type in special_sequences:
        print(f"\n{seq_type}:")
        seq_df = df[df['sequence_type'] == seq_type]
        

        extreme_df = seq_df[seq_df['convergence'] == 'extreme']
        
        for _, row in extreme_df.iterrows():
            print(f"  {row['method']:15}: 重复率={row['repetition_rate']:5.1f}%, "
                  f"唯一序列={row['unique_sequences']}/50")
            
            if row['method'] == 'hybrid_bs_sado' and row['switch_events'] > 0:
                print(f"    → 触发SADO切换{row['switch_events']}次，有效降低重复")


def test_switch_effectiveness():

    
    print("\n" + "="*80)

    print("="*80)
    


    valid_codon_mask, codon_indices = prepare_inputs(sequence)
    
    operator = CAIEnhancementOperator(
        method='hybrid_bs_sado',
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    

    num_positions = len(sequence)
    pi_accessibility = torch.zeros(num_positions, 6)

    



    
    for i in range(10):
        discrete_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=pi_accessibility,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_codon_mask,
            codon_indices=codon_indices,
            target_cai=0.8
        )
        
        switched = metadata.get('switched_to_sado', False)

        
        if switched:

            if 'bs_gamma' in metadata:

        else:

    

    if hasattr(operator.optimizer, 'get_statistics'):
        stats = operator.optimizer.get_statistics()







if __name__ == "__main__":


    

    results = test_different_sequence_types()
    

    create_summary_table(results)
    

    test_switch_effectiveness()
    

    df = pd.DataFrame(results)
    df.to_csv('short_sequence_repetition_results.csv', index=False)

    
    print("\n" + "="*80)

    print("="*80)