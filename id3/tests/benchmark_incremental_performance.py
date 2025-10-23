#!/usr/bin/env python3
"""






"""

import sys
import time
import torch
import numpy as np
import logging
from pathlib import Path
from typing import Dict, List, Tuple
import gc


project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from id3.optimizers.cai.incremental import IncrementalCAIOptimizer


logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)


def generate_test_sequence(length: int) -> str:

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    np.random.seed(42)
    return ''.join(np.random.choice(amino_acids, length))


def benchmark_optimizer(
    seq_length: int,
    num_iterations: int = 10,
    enable_caching: bool = True
) -> Dict:
    """
    对优化器进行基准测试

    Args:
        seq_length: 序列长度
        num_iterations: 迭代次数
        enable_caching: 是否启用缓存

    Returns:
        性能指标字典
    """

    amino_acid_sequence = generate_test_sequence(seq_length)


    optimizer = IncrementalCAIOptimizer(
        amino_acid_sequence=amino_acid_sequence,
        random_seed=42,
        enable_caching=enable_caching
    )


    pi_accessibility = torch.rand(seq_length, 64)
    valid_codon_mask = torch.ones(seq_length, 64)
    target_cai = 0.8


    optimizer.optimize(
        pi_accessibility,
        target_cai,
        amino_acid_sequence,
        valid_codon_mask=valid_codon_mask
    )


    times = []
    for i in range(num_iterations):

        pi_accessibility_iter = pi_accessibility + torch.randn_like(pi_accessibility) * 0.01
        pi_accessibility_iter = torch.clamp(pi_accessibility_iter, 0, 1)

        start_time = time.perf_counter()

        result, metadata = optimizer.optimize(
            pi_accessibility_iter,
            target_cai,
            amino_acid_sequence,
            valid_codon_mask=valid_codon_mask
        )

        end_time = time.perf_counter()
        times.append(end_time - start_time)


    stats = optimizer.get_statistics()


    return {
        'seq_length': seq_length,
        'caching_enabled': enable_caching,
        'avg_time': np.mean(times),
        'std_time': np.std(times),
        'min_time': np.min(times),
        'max_time': np.max(times),
        'median_time': np.median(times),
        'cache_hit_rate': stats.get('cache_hit_rate', 0) if enable_caching else 0,
        'total_iterations': stats['iteration_count'],
        'unique_sequences': stats['unique_sequences']
    }


def compare_caching_performance(seq_lengths: List[int], num_iterations: int = 10):
    """比较启用和禁用缓存的性能差异"""
    results = []

    for seq_length in seq_lengths:
        print(f"\n测试序列长度: {seq_length}")


        gc.collect()
        torch.cuda.empty_cache() if torch.cuda.is_available() else None


        print("  - 测试无缓存版本...")
        result_no_cache = benchmark_optimizer(seq_length, num_iterations, enable_caching=False)


        gc.collect()
        torch.cuda.empty_cache() if torch.cuda.is_available() else None


        print("  - 测试有缓存版本...")
        result_with_cache = benchmark_optimizer(seq_length, num_iterations, enable_caching=True)


        speedup = result_no_cache['avg_time'] / result_with_cache['avg_time']

        results.append({
            'seq_length': seq_length,
            'no_cache_time': result_no_cache['avg_time'],
            'with_cache_time': result_with_cache['avg_time'],
            'speedup': speedup,
            'cache_hit_rate': result_with_cache['cache_hit_rate']
        })

        print(f"  ✓ 完成 - 加速比: {speedup:.2f}x")

    return results


def test_incremental_vs_binary_search(seq_length: int = 100, num_iterations: int = 20):



    amino_acid_sequence = generate_test_sequence(seq_length)


    optimizer = IncrementalCAIOptimizer(
        amino_acid_sequence=amino_acid_sequence,
        random_seed=42,
        enable_caching=True
    )


    pi_accessibility = torch.rand(seq_length, 64)
    valid_codon_mask = torch.ones(seq_length, 64)
    target_cai = 0.8

    times = []
    methods_used = []

    for i in range(num_iterations):

        pi_iter = pi_accessibility + torch.randn_like(pi_accessibility) * 0.01
        pi_iter = torch.clamp(pi_iter, 0, 1)

        start_time = time.perf_counter()
        result, metadata = optimizer.optimize(
            pi_iter,
            target_cai,
            amino_acid_sequence,
            valid_codon_mask=valid_codon_mask
        )
        end_time = time.perf_counter()

        times.append(end_time - start_time)
        methods_used.append(metadata['method'])


    bs_times = [t for t, m in zip(times, methods_used) if 'bs' in m]
    inc_times = [t for t, m in zip(times, methods_used) if 'inc' in m]

    if len(inc_times) > 0 and len(bs_times) > 0:
        avg_bs_time = np.mean(bs_times)
        avg_inc_time = np.mean(inc_times)
        speedup = avg_bs_time / avg_inc_time if avg_inc_time > 0 else 1.0




    else:


    return optimizer.get_statistics()


def print_performance_table(results: List[Dict]):
    """打印性能对比表格"""
    print("\n" + "="*80)
    print("性能基准测试结果")
    print("="*80)
    print(f"{'序列长度':<10} {'无缓存(ms)':<12} {'有缓存(ms)':<12} {'加速比':<10} {'缓存命中率':<12}")
    print("-"*80)

    for r in results:
        print(f"{r['seq_length']:<10} "
              f"{r['no_cache_time']*1000:<12.2f} "
              f"{r['with_cache_time']*1000:<12.2f} "
              f"{r['speedup']:<10.2f}x "
              f"{r['cache_hit_rate']:<12.1%}")


    avg_speedup = np.mean([r['speedup'] for r in results])
    print("-"*80)
    print(f"平均加速比: {avg_speedup:.2f}x")
    print("="*80)
    return avg_speedup


def main():

    print("="*80)

    print("="*80)


    seq_lengths = [50, 100, 200, 500, 1000]


    print("-"*40)
    cache_results = compare_caching_performance(seq_lengths, num_iterations=10)
    avg_speedup = print_performance_table(cache_results)



    print("-"*40)
    inc_stats = test_incremental_vs_binary_search(seq_length=200, num_iterations=20)



    print("-"*40)

    optimizer = IncrementalCAIOptimizer(
        amino_acid_sequence=generate_test_sequence(100),
        enable_caching=True
    )


    import timeit


    k_calc_time = timeit.timeit(
        lambda: optimizer._calculate_k_positions(200),
        number=10000


    norm_prob_time = timeit.timeit(
        lambda: optimizer._normalize_probabilities(np.random.rand(100)),
        number=10000
    ) / 10000 * 1000







    print("-"*40)

    cache_size = len(optimizer._cai_cache) if hasattr(optimizer, '_cai_cache') else 0








    print("\n" + "="*80)

    print("="*80)

    if avg_speedup > 1.0:




    else:







    print("="*80)


if __name__ == "__main__":
    main()