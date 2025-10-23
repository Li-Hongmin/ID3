#!/usr/bin/env python3
"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np
from typing import List, Dict, Tuple

from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging

logger = setup_logging(level='INFO', name='binary_vs_greedy')

class BinarySearchCAIOptimizer:
    """

    


    """
    
    def __init__(self, sado_optimizer: SADOOptimizer):
        self.sado = sado_optimizer
    
    def binary_search_optimize(self, pi_accessibility: torch.Tensor, 
                              target_cai: float = 0.8, 
                              epsilon: float = 1e-6,
                              max_iterations: int = 100) -> Tuple[np.ndarray, Dict]:
        """

        

        """
        seq_len = pi_accessibility.shape[0]
        indices = np.zeros(seq_len, dtype=np.int32)
        
        start_time = time.time()
        

        for pos in range(seq_len):
            if pos >= len(self.sado.codon_choices) or not self.sado.codon_choices[pos]:
                continue
                
            choices = self.sado.codon_choices[pos]
            probs = pi_accessibility[pos].detach().cpu().numpy()
            

            best_score = -1
            best_idx = 0
            

            for choice_idx, choice in enumerate(choices):
                orig_idx = choice['original_local_index']
                prob = probs[orig_idx] if orig_idx < len(probs) else 1e-8
                cai_weight = choice['weight']
                
                if prob > 0 and cai_weight > 0:

                    gamma = self._binary_search_gamma(prob, cai_weight, epsilon)
                    score = (prob ** (1 - gamma)) * (cai_weight ** gamma)
                    
                    if score > best_score:
                        best_score = score
                        best_idx = choice_idx
            
            indices[pos] = best_idx
        
        execution_time = (time.time() - start_time) * 1000
        

        final_cai = self.sado._compute_cai_from_indices(indices)
        
        metadata = {
            'execution_time_ms': execution_time,
            'final_cai': final_cai,
            'target_cai': target_cai,
            'constraint_satisfied': final_cai >= target_cai,
            'method': 'binary_search',
            'theoretical_complexity': f'O({seq_len} × {4} × {int(np.log2(1/epsilon))}) ≈ {seq_len * 4 * int(np.log2(1/epsilon))}'
        }
        
        return indices, metadata
    
    def _binary_search_gamma(self, prob: float, cai_weight: float, epsilon: float) -> float:
        """

        """
        left, right = 0.0, 1.0
        best_gamma = 0.5
        best_score = 0.0
        

        while right - left > epsilon:
            gamma1 = left + (right - left) / 3
            gamma2 = right - (right - left) / 3
            
            score1 = (prob ** (1 - gamma1)) * (cai_weight ** gamma1)
            score2 = (prob ** (1 - gamma2)) * (cai_weight ** gamma2)
            
            if score1 > score2:
                right = gamma2
                if score1 > best_score:
                    best_score = score1
                    best_gamma = gamma1
            else:
                left = gamma1
                if score2 > best_score:
                    best_score = score2
                    best_gamma = gamma2
        
        return best_gamma


class PerformanceComparator:

    
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    
    def generate_test_sequence(self, length: int) -> str:
        """生成测试序列"""
        amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
        return ''.join(np.random.choice(amino_acids, length))
    
    def create_test_distribution(self, seq_length: int, amino_sequence: str):

        num_codons = 6
        valid_mask = torch.zeros(seq_length, num_codons, dtype=torch.bool, device=self.device)
        
        for pos, aa in enumerate(amino_sequence):
            if aa in amino_acids_to_codons:
                num_valid = len(amino_acids_to_codons[aa])
                valid_mask[pos, :min(num_valid, num_codons)] = True
        
        pi = torch.rand(seq_length, num_codons, device=self.device)
        pi = pi * valid_mask.float()
        
        for pos in range(seq_length):
            if valid_mask[pos].any():
                pi[pos] = pi[pos] / pi[pos].sum()
        
        return pi
    
    def compare_methods(self, test_lengths: List[int]):
        """对比不同方法的性能"""
        logger.info("=" * 80)
        logger.info("二分查找 vs 贪婪搜索性能对比")
        logger.info("=" * 80)
        
        results = []
        
        for length in test_lengths:
            logger.info(f"\n测试序列长度: {length} 氨基酸")
            logger.info("-" * 60)
            

            test_sequence = self.generate_test_sequence(length)
            pi_distribution = self.create_test_distribution(length, test_sequence)
            

            sado_optimizer = SADOOptimizer(
                species='ecoli_bl21de3',
                device=self.device,
                amino_acid_sequence=test_sequence
            )
            
            binary_optimizer = BinarySearchCAIOptimizer(sado_optimizer)
            

            logger.info("测试二分查找方法...")
            binary_indices, binary_meta = binary_optimizer.binary_search_optimize(
                pi_distribution, target_cai=0.8
            )
            
            logger.info(f"  二分查找: {binary_meta['execution_time_ms']:.1f}ms, "
                       f"CAI={binary_meta['final_cai']:.4f}, "
                       f"满足约束={'✅' if binary_meta['constraint_satisfied'] else '❌'}")
            

            logger.info("测试SADO贪婪方法...")
            start_time = time.time()
            timeout_seconds = min(30, length / 100)
            
            try:
                sado_optimizer.reset()
                _, sado_meta = sado_optimizer.optimize(
                    pi_accessibility=pi_distribution,
                    target_cai=0.8,
                    amino_acid_sequence=test_sequence,
                    gamma=0.0
                )
                sado_time = (time.time() - start_time) * 1000
                
                logger.info(f"  SADO贪婪: {sado_time:.1f}ms, "
                           f"CAI={sado_meta['final_cai']:.4f}, "
                           f"满足约束={'✅' if sado_meta['constraint_satisfied'] else '❌'}")
                

                if sado_time > 0:
                    speedup = sado_time / binary_meta['execution_time_ms']
                    logger.info(f"  二分查找加速比: {speedup:.1f}x")
                
            except Exception as e:
                logger.error(f"  SADO贪婪方法失败: {e}")
                sado_time = float('inf')
                speedup = float('inf')
            

            results.append({
                'length': length,
                'binary_time_ms': binary_meta['execution_time_ms'],
                'binary_cai': binary_meta['final_cai'],
                'binary_satisfied': binary_meta['constraint_satisfied'],
                'sado_time_ms': sado_time,
                'sado_cai': sado_meta.get('final_cai', 0.0) if 'sado_meta' in locals() else 0.0,
                'sado_satisfied': sado_meta.get('constraint_satisfied', False) if 'sado_meta' in locals() else False,
                'speedup': speedup if speedup != float('inf') else 'timeout'
            })
        

        self.analyze_results(results)
        
        return results
    
    def analyze_results(self, results: List[Dict]):

        logger.info("=" * 80)

        logger.info("=" * 80)
        


        logger.info("-" * 60)
        
        for result in results:
            binary_time = f"{result['binary_time_ms']:.1f}ms"
            sado_time = f"{result['sado_time_ms']:.1f}ms" if result['sado_time_ms'] != float('inf') else "timeout"
            speedup = f"{result['speedup']:.1f}x" if isinstance(result['speedup'], float) and result['speedup'] != float('inf') else "∞"
            cai_diff = result['binary_cai'] - result['sado_cai']
            
            logger.info(f"{result['length']:<8} {binary_time:<12} {sado_time:<12} {speedup:<10} {cai_diff:+.3f}")
        


        

        if len(results) >= 2:
            for i in range(1, len(results)):
                prev_result = results[i-1]
                curr_result = results[i]
                
                length_ratio = curr_result['length'] / prev_result['length']
                time_ratio = curr_result['binary_time_ms'] / prev_result['binary_time_ms']
                

                


        


        working_results = [r for r in results if isinstance(r['speedup'], float) and r['speedup'] != float('inf')]
        if working_results:
            avg_speedup = np.mean([r['speedup'] for r in working_results])

            
            if avg_speedup > 10:

            elif avg_speedup > 3:

            else:

        



def main():
    """主函数"""
    logger.info("🚀 开始二分查找vs贪婪搜索性能对比")
    
    comparator = PerformanceComparator()
    

    test_lengths = [100, 300, 500, 1000, 2000, 3000]
    
    results = comparator.compare_methods(test_lengths)
    
    logger.info("✅ 对比测试完成!")
    
    return results


if __name__ == "__main__":
    main()