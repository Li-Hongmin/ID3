#!/usr/bin/env python3
"""


"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import time
import random

from id3.optimizers.cai.sado import SADOOptimizer


def generate_random_amino_sequence(length: int) -> str:

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    return ''.join(random.choice(amino_acids) for _ in range(length))


def test_sequence_length(length: int, timeout: float = 10.0):
    """测试特定长度序列的性能"""
    
    print(f"\n{'='*60}")
    print(f"测试序列长度: {length} AA")
    print(f"{'='*60}")
    

    amino_sequence = generate_random_amino_sequence(length)
    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(length, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    results = {}
    

    configs = [
        ('标准SADO', False, False),
        ('增强SADO(差异驱动)', False, True),
        ('完整增强(二分+差异)', True, True)
    ]
    
    for name, use_binary, use_diff in configs:
        print(f"\n{name}:")
        
        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        
        start_time = time.time()
        try:

            import signal
            
            def timeout_handler(signum, frame):
                raise TimeoutError("超时")
            

            if hasattr(signal, 'SIGALRM'):
                signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(int(timeout))
            
            dist, meta = optimizer.optimize(
                pi_accessibility=pi_accessibility.clone(),
                target_cai=0.8,
                use_binary_search=use_binary,
                use_difference_driven=use_diff,
                gamma=0.3
            )
            
            if hasattr(signal, 'SIGALRM'):
                signal.alarm(0)
            
            elapsed = time.time() - start_time
            
            results[name] = {
                'success': True,
                'time': elapsed,
                'cai': meta['final_cai'],
                'satisfied': meta['constraint_satisfied']
            }
            
            print(f"  ✅ 成功: CAI={meta['final_cai']:.4f}, 时间={elapsed:.3f}秒")
            
        except (TimeoutError, Exception) as e:
            elapsed = time.time() - start_time
            results[name] = {
                'success': False,
                'time': elapsed,
                'cai': 0,
                'satisfied': False
            }
            print(f"  ❌ 失败: {str(e)[:50]}, 时间={elapsed:.3f}秒")
            
            if hasattr(signal, 'SIGALRM'):
                signal.alarm(0)
    
    return results


def main():

    
    print("\n" + "="*80)

    print("="*80)
    

    test_lengths = [50, 100, 200, 300, 500, 1000, 1500, 2000, 3000]
    
    all_results = {}
    
    for length in test_lengths:

        timeout = 30.0 if length > 1000 else 10.0
        
        try:
            results = test_sequence_length(length, timeout=timeout)
            all_results[length] = results
        except Exception as e:

            continue
    

    print("\n" + "="*80)

    print("="*80)
    

    print("-"*70)
    
    for length in sorted(all_results.keys()):
        results = all_results[length]
        
        line = f"{length:<10}"
        

            if method in results:
                r = results[method]
                if r['success']:
                    line += f"{r['time']:.2f}s/{r['cai']:.3f}  "
                else:

            else:
                line += "-"*18 + "  "
        
        print(line)
    

    print("\n" + "="*80)

    print("="*80)
    







if __name__ == '__main__':
    main()