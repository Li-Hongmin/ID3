#!/usr/bin/env python3
"""


"""

import torch
import time
import logging
from typing import Dict, Any
import numpy as np
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


TEST_SEQUENCE = "MQYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE"
NUM_ITERATIONS = 20
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'


def test_cpc_version(version: str, enable_cai: bool = False) -> Dict[str, Any]:
    """

    
    Args:


        
    Returns:

    """
    logger.info(f"\n{'='*60}")
    logger.info(f"测试CPC {version} (CAI={'启用' if enable_cai else '禁用'})")
    logger.info(f"{'='*60}")
    

    deepraccess_model = DeepRaccessID3Wrapper(device=DEVICE)
    

    if version == 'v2':
        from id3.constraints.cpc_v2_efficient import CodonProfileConstraintV2Efficient
        constraint = CodonProfileConstraintV2Efficient(
            amino_acid_sequence=TEST_SEQUENCE,
            deepraccess_model=deepraccess_model,
            enable_cai=enable_cai,
            device=DEVICE
        )
    elif version == 'v3':
        from id3.constraints.cpc_v3_stable import CodonProfileConstraintV3Stable
        constraint = CodonProfileConstraintV3Stable(
            amino_acid_sequence=TEST_SEQUENCE,
            deepraccess_model=deepraccess_model,
            enable_cai=enable_cai,
            device=DEVICE
        )
    elif version == 'v4':
        from id3.constraints.cpc_v4_hybrid import CodonProfileConstraintV4Hybrid

        strategies = ['smart', 'fast', 'safe']
        strategy = strategies[0]
        
        constraint = CodonProfileConstraintV4Hybrid(
            amino_acid_sequence=TEST_SEQUENCE,
            deepraccess_model=deepraccess_model,
            enable_cai=enable_cai,
            strategy=strategy,
            enable_pool=True,
            device=DEVICE
        )
        logger.info(f"V4策略: {strategy}")
    else:
        raise ValueError(f"未知版本: {version}")
    

    forward_times = []
    loss_values = []
    accessibility_values = []
    cai_values = [] if enable_cai else None
    

    logger.info("预热中...")
    for _ in range(3):
        _ = constraint.forward_with_loss(alpha=0.5, tau=1.0, beta=0.0)
    

    logger.info(f"开始测试 {NUM_ITERATIONS} 次迭代...")
    
    for i in range(NUM_ITERATIONS):
        start_time = time.time()
        

        result = constraint.forward_with_loss(
            alpha=0.5,
            tau=1.0,
            beta=0.0
        )
        
        forward_time = time.time() - start_time
        forward_times.append(forward_time)
        

        loss_values.append(result.get('loss_access', 0))
        accessibility_values.append(result.get('eval_access', 0))
        if enable_cai and 'eval_cai' in result:
            cai_values.append(result['eval_cai'])
        

        if (i + 1) % 5 == 0:
            avg_time = np.mean(forward_times) * 1000
            logger.info(f"  进度: {i+1}/{NUM_ITERATIONS}, 平均时间: {avg_time:.2f}ms")
    

    stats = {
        'version': version,
        'enable_cai': enable_cai,
        'num_iterations': NUM_ITERATIONS,
        'forward_times_ms': {
            'mean': np.mean(forward_times) * 1000,
            'std': np.std(forward_times) * 1000,
            'min': np.min(forward_times) * 1000,
            'max': np.max(forward_times) * 1000,
            'p50': np.percentile(forward_times, 50) * 1000,
            'p95': np.percentile(forward_times, 95) * 1000,
            'p99': np.percentile(forward_times, 99) * 1000
        },
        'loss_values': {
            'mean': np.mean(loss_values),
            'std': np.std(loss_values),
            'min': np.min(loss_values),
            'max': np.max(loss_values)
        },
        'accessibility_values': {
            'mean': np.mean(accessibility_values),
            'std': np.std(accessibility_values),
            'min': np.min(accessibility_values),
            'max': np.max(accessibility_values)
        }
    }
    
    if cai_values:
        stats['cai_values'] = {
            'mean': np.mean(cai_values),
            'std': np.std(cai_values),
            'min': np.min(cai_values),
            'max': np.max(cai_values)
        }
    

    if hasattr(constraint, 'get_performance_stats'):
        version_stats = constraint.get_performance_stats()
        stats['version_specific'] = version_stats
    
    return stats


def print_comparison_table(results: list):

    print("\n" + "="*80)

    print("="*80)
    


    print("-"*60)

    print("-"*60)
    
    for r in results:
        t = r['forward_times_ms']
        cai_str = " (CAI)" if r['enable_cai'] else ""
        print(f"{r['version']+cai_str:<10} {t['mean']:<12.2f} {t['std']:<12.2f} "
              f"{t['p50']:<12.2f} {t['p95']:<12.2f} {t['p99']:<12.2f}")
    


    print("-"*60)

    print("-"*60)
    
    for r in results:
        l = r['loss_values']
        cai_str = " (CAI)" if r['enable_cai'] else ""
        print(f"{r['version']+cai_str:<10} {l['mean']:<15.6f} {l['std']:<15.6f} "
              f"{l['min']:<15.6f} {l['max']:<15.6f}")
    


    print("-"*60)

    print("-"*60)
    
    for r in results:
        a = r['accessibility_values']
        cai_str = " (CAI)" if r['enable_cai'] else ""
        print(f"{r['version']+cai_str:<10} {a['mean']:<15.6f} {a['std']:<15.6f} "
              f"{a['min']:<15.6f} {a['max']:<15.6f}")
    

    cai_results = [r for r in results if 'cai_values' in r]
    if cai_results:

        print("-"*60)

        print("-"*60)
        
        for r in cai_results:
            c = r['cai_values']
            print(f"{r['version']:<10} {c['mean']:<15.6f} {c['std']:<15.6f} "
                  f"{c['min']:<15.6f} {c['max']:<15.6f}")
    


    print("-"*60)
    
    for r in results:
        if 'version_specific' in r:
            cai_str = " (CAI)" if r['enable_cai'] else ""
            print(f"\n{r['version']}{cai_str}:")
            vs = r['version_specific']
            

            if 'success_rate' in vs:

            if 'reset_count' in vs:

            if 'strategies' in vs:

                for strategy, stats in vs['strategies'].items():
                    if stats['calls'] > 0:
                        success_rate = stats['success'] / stats['calls'] * 100
                        avg_time = stats['time'] / stats['calls'] * 1000




def main():
    """主测试函数"""
    print("\n" + "="*80)
    print("CPC版本性能对比测试")
    print("="*80)
    print(f"测试序列长度: {len(TEST_SEQUENCE)} 氨基酸")
    print(f"测试迭代次数: {NUM_ITERATIONS}")
    print(f"设备: {DEVICE}")
    
    results = []
    

    for version in ['v2', 'v3', 'v4']:
        try:
            stats = test_cpc_version(version, enable_cai=False)
            results.append(stats)
        except Exception as e:
            logger.error(f"测试{version}失败: {e}")
    

    print("\n" + "="*80)
    print("测试CAI模式")
    print("="*80)
    
    for version in ['v2', 'v3', 'v4']:
        try:
            stats = test_cpc_version(version, enable_cai=True)
            results.append(stats)
        except Exception as e:
            logger.error(f"测试{version} (CAI)失败: {e}")
    

    print_comparison_table(results)
    

    print("\n" + "="*80)
    print("💡 性能分析与建议")
    print("="*80)
    

    no_cai_results = [r for r in results if not r['enable_cai']]
    if no_cai_results:
        fastest = min(no_cai_results, key=lambda x: x['forward_times_ms']['mean'])
        slowest = max(no_cai_results, key=lambda x: x['forward_times_ms']['mean'])
        
        speed_ratio = slowest['forward_times_ms']['mean'] / fastest['forward_times_ms']['mean']
        
        print(f"\n无CAI模式:")
        print(f"  最快: {fastest['version']} ({fastest['forward_times_ms']['mean']:.2f}ms)")
        print(f"  最慢: {slowest['version']} ({slowest['forward_times_ms']['mean']:.2f}ms)")
        print(f"  速度差异: {speed_ratio:.1f}倍")
    
    cai_results = [r for r in results if r['enable_cai']]
    if cai_results:
        fastest = min(cai_results, key=lambda x: x['forward_times_ms']['mean'])
        slowest = max(cai_results, key=lambda x: x['forward_times_ms']['mean'])
        
        speed_ratio = slowest['forward_times_ms']['mean'] / fastest['forward_times_ms']['mean']
        
        print(f"\nCAI模式:")
        print(f"  最快: {fastest['version']} ({fastest['forward_times_ms']['mean']:.2f}ms)")
        print(f"  最慢: {slowest['version']} ({slowest['forward_times_ms']['mean']:.2f}ms)")
        print(f"  速度差异: {speed_ratio:.1f}倍")
    
    print("\n建议:")
    print("  • 生产环境: 使用v4 (smart策略) 或 v3，平衡性能和可靠性")
    print("  • 开发调试: 使用v2，代码简单易调试")
    print("  • 长时间训练: 使用v4 (enable_pool=True)，最佳性能")
    print("  • 高可靠性要求: 使用v4 (safe策略) 或 v2")


if __name__ == "__main__":
    main()