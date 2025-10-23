#!/usr/bin/env python3
"""












"""

import torch
import numpy as np
import time
import hashlib
import random
import string
from collections import defaultdict
from typing import Dict, List, Tuple
import pandas as pd
from tabulate import tabulate
import matplotlib.pyplot as plt
import seaborn as sns

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


torch.manual_seed(42)
np.random.seed(42)
random.seed(42)


def generate_random_protein(length: int) -> str:


    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    weights = [0.074, 0.042, 0.044, 0.059, 0.037, 0.074, 0.022, 0.058, 0.051, 
               0.090, 0.057, 0.036, 0.051, 0.043, 0.052, 0.073, 0.056, 0.063, 
               0.032, 0.013]
    
    sequence = ''.join(random.choices(amino_acids, weights=weights, k=length))
    return sequence


def calculate_probability_score(pi_accessibility: torch.Tensor, 
                               discrete_dist: torch.Tensor,
                               valid_mask: torch.Tensor) -> float:
    """计算概率保持得分（选中密码子的概率之和）"""
    if discrete_dist.dim() > 1:
        indices = discrete_dist.argmax(dim=-1)
    else:
        indices = discrete_dist
    
    score = 0.0
    for pos in range(len(indices)):
        if valid_mask[pos].any():
            selected_codon = indices[pos].item()
            if selected_codon < pi_accessibility.shape[1]:
                prob = pi_accessibility[pos, selected_codon].item()
                score += np.log(prob + 1e-10)
    
    return score / len(indices)


def prepare_inputs(sequence: str) -> Tuple[torch.Tensor, torch.Tensor]:

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


def benchmark_method(method: str, 
                     sequence: str, 
                     iterations: int = 100,
                     target_cai: float = 0.8) -> Dict:
    """对单个方法进行基准测试"""
    
    print(f"  测试 {method}...")
    
    valid_codon_mask, codon_indices = prepare_inputs(sequence)
    

    operator = CAIEnhancementOperator(
        method=method,
        species='ecoli_bl21de3',
        amino_acid_sequence=sequence
    )
    

    times = []
    sequences_generated = []
    cai_values = []
    prob_scores = []
    switch_events = 0
    gamma_values = []
    

    num_positions = len(sequence)
    
    for i in range(iterations):

        if i < iterations // 3:

            pi_accessibility = torch.rand(num_positions, 6)
        elif i < 2 * iterations // 3:

            pi_accessibility = torch.rand(num_positions, 6) ** 2
        else:

            pi_accessibility = torch.rand(num_positions, 6) ** 4
        

        pi_accessibility = pi_accessibility * valid_codon_mask.float()
        for pos in range(num_positions):
            if valid_codon_mask[pos].any():
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        start_time = time.time()
        try:
            discrete_dist, metadata = operator.apply_cai_enhancement(
                pi_accessibility=pi_accessibility,
                amino_acid_sequence=sequence,
                valid_codon_mask=valid_codon_mask,
                codon_indices=codon_indices,
                target_cai=target_cai
            )
            elapsed_time = time.time() - start_time
            times.append(elapsed_time)
            

            if discrete_dist.dim() > 1:
                indices = discrete_dist.argmax(dim=-1)
            else:
                indices = discrete_dist
            
            seq_hash = hashlib.md5(indices.cpu().numpy().tobytes()).hexdigest()
            sequences_generated.append(seq_hash)
            

            if 'achieved_cai' in metadata:
                cai_values.append(metadata['achieved_cai'])
            elif 'final_cai' in metadata:
                cai_values.append(metadata['final_cai'])
            

            prob_score = calculate_probability_score(pi_accessibility, discrete_dist, valid_codon_mask)
            prob_scores.append(prob_score)
            

            if method == 'hybrid_bs_sado' and metadata.get('switched_to_sado', False):
                switch_events += 1
                if 'bs_gamma' in metadata:
                    gamma_values.append(metadata['bs_gamma'])
            
        except Exception as e:
            print(f"    警告：迭代 {i} 失败: {e}")
            continue
    

    unique_sequences = len(set(sequences_generated))
    repetition_rate = 1 - (unique_sequences / len(sequences_generated)) if sequences_generated else 0
    
    result = {
        'method': method,
        'sequence_length': len(sequence),
        'iterations': len(times),
        'avg_time_ms': np.mean(times) * 1000 if times else 0,
        'std_time_ms': np.std(times) * 1000 if times else 0,
        'total_time_s': sum(times) if times else 0,
        'unique_sequences': unique_sequences,
        'repetition_rate': repetition_rate * 100,
        'avg_cai': np.mean(cai_values) if cai_values else 0,
        'cai_success_rate': sum(1 for c in cai_values if c >= target_cai) / len(cai_values) * 100 if cai_values else 0,
        'avg_prob_score': np.mean(prob_scores) if prob_scores else 0,
        'switch_events': switch_events,
        'switch_rate': switch_events / len(times) * 100 if times else 0,
    }
    

    if method == 'hybrid_bs_sado' and hasattr(operator.optimizer, 'get_statistics'):
        stats = operator.optimizer.get_statistics()
        result['bs_calls'] = stats.get('bs_calls', 0)
        result['sado_calls'] = stats.get('sado_calls', 0)
    
    return result


def run_comprehensive_benchmark():

    
    print("="*80)

    print("="*80)
    

    sequence_lengths = [200, 500, 1000, 2000, 3000]
    methods = ['binary_search', 'sado', 'hybrid_bs_sado']

    
    all_results = []
    
    for length in sequence_lengths:

        print("-"*60)
        

        sequence = generate_random_protein(length)
        
        for method in methods:
            result = benchmark_method(method, sequence, iterations)
            all_results.append(result)
            

            print(f"    {method}: {result['avg_time_ms']:.2f}ms, "


    
    return all_results


def create_comparison_table(results: List[Dict]):
    """创建比较表格"""
    
    print("\n" + "="*100)
    print("性能比较表格")
    print("="*100)
    

    df = pd.DataFrame(results)
    

    for length in df['sequence_length'].unique():
        print(f"\n序列长度: {length} aa")
        print("-"*90)
        
        length_df = df[df['sequence_length'] == length].copy()
        

        display_cols = [
            'method',
            'avg_time_ms',
            'unique_sequences', 
            'repetition_rate',
            'avg_cai',
            'cai_success_rate',
            'avg_prob_score',
            'switch_rate'
        ]
        

        table_data = []
        for _, row in length_df.iterrows():
            formatted_row = [
                row['method'],
                f"{row['avg_time_ms']:.2f}",
                f"{row['unique_sequences']}/{row['iterations']}",
                f"{row['repetition_rate']:.1f}%",
                f"{row['avg_cai']:.3f}",
                f"{row['cai_success_rate']:.1f}%",
                f"{row['avg_prob_score']:.3f}",
                f"{row['switch_rate']:.1f}%" if row['method'] == 'hybrid_bs_sado' else "N/A"
            ]
            table_data.append(formatted_row)
        
        headers = ['方法', '平均时间(ms)', '唯一序列', '重复率', '平均CAI', 'CAI成功率', '概率得分', '切换率']
        print(tabulate(table_data, headers=headers, tablefmt='grid'))
    

    print("\n" + "="*100)
    print("汇总统计")
    print("="*100)
    
    summary_data = []
    for method in df['method'].unique():
        method_df = df[df['method'] == method]
        summary_data.append([
            method,
            f"{method_df['avg_time_ms'].mean():.2f}",
            f"{method_df['repetition_rate'].mean():.1f}%",
            f"{method_df['cai_success_rate'].mean():.1f}%",
            f"{method_df['avg_prob_score'].mean():.3f}"
        ])
    
    summary_headers = ['方法', '平均时间(ms)', '平均重复率', '平均CAI成功率', '平均概率得分']
    print(tabulate(summary_data, headers=summary_headers, tablefmt='grid'))
    
    return df


def plot_performance_comparison(df: pd.DataFrame):

    
    try:
        import matplotlib

        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        

        ax = axes[0, 0]
        for method in df['method'].unique():
            method_df = df[df['method'] == method]
            ax.plot(method_df['sequence_length'], method_df['avg_time_ms'], 
                   marker='o', label=method, linewidth=2)



        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[0, 1]
        for method in df['method'].unique():
            method_df = df[df['method'] == method]
            ax.plot(method_df['sequence_length'], method_df['repetition_rate'],
                   marker='s', label=method, linewidth=2)



        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[0, 2]
        for method in df['method'].unique():
            method_df = df[df['method'] == method]
            ax.plot(method_df['sequence_length'], method_df['cai_success_rate'],
                   marker='^', label=method, linewidth=2)



        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[1, 0]
        width = 0.25
        x = np.arange(len(df['sequence_length'].unique()))
        for i, method in enumerate(df['method'].unique()):
            method_df = df[df['method'] == method]
            ax.bar(x + i*width, method_df['unique_sequences'], 
                  width, label=method)



        ax.set_xticks(x + width)
        ax.set_xticklabels(df['sequence_length'].unique())
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        

        ax = axes[1, 1]
        for method in df['method'].unique():
            method_df = df[df['method'] == method]
            ax.plot(method_df['sequence_length'], method_df['avg_prob_score'],
                   marker='d', label=method, linewidth=2)



        ax.legend()
        ax.grid(True, alpha=0.3)
        

        ax = axes[1, 2]
        hybrid_df = df[df['method'] == 'hybrid_bs_sado']
        if not hybrid_df.empty:
            ax.bar(hybrid_df['sequence_length'], hybrid_df['switch_rate'],
                  color='orange', alpha=0.7)



            ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig('cai_optimizer_benchmark.png', dpi=150, bbox_inches='tight')

        
    except ImportError:

    except Exception as e:



def analyze_switch_behavior(results: List[Dict]):
    """分析混合方法的切换行为"""
    
    print("\n" + "="*100)
    print("Hybrid方法切换行为分析")
    print("="*100)
    
    hybrid_results = [r for r in results if r['method'] == 'hybrid_bs_sado']
    
    if not hybrid_results:
        print("没有hybrid方法的结果")
        return
    
    for result in hybrid_results:
        print(f"\n序列长度: {result['sequence_length']} aa")
        print(f"  总迭代次数: {result['iterations']}")
        print(f"  切换到SADO: {result['switch_events']}次")
        print(f"  切换率: {result['switch_rate']:.1f}%")
        
        if 'bs_calls' in result:
            print(f"  Binary Search调用: {result['bs_calls']}次")
        if 'sado_calls' in result:
            print(f"  SADO调用: {result['sado_calls']}次")
        

        if result['switch_events'] > 0:
            print(f"  切换效果：")
            print(f"    - 生成{result['unique_sequences']}个唯一序列")
            print(f"    - 重复率降至{result['repetition_rate']:.1f}%")


if __name__ == "__main__":
    print("开始大规模CAI优化器基准测试...")
    print("配置：100次迭代，序列长度200-3000aa")
    print("这可能需要几分钟时间...\n")
    

    results = run_comprehensive_benchmark()
    

    df = create_comparison_table(results)
    

    analyze_switch_behavior(results)
    

    plot_performance_comparison(df)
    

    df.to_csv('cai_optimizer_benchmark_results.csv', index=False)
    print("\n详细结果已保存到 cai_optimizer_benchmark_results.csv")
    
    print("\n" + "="*100)
    print("测试完成！")
    print("="*100)