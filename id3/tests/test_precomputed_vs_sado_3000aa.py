#!/usr/bin/env python3
"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np
from typing import List, Dict


from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


import sys
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/experiments/cai_enhancement_theory/04_method_comparison')
from precomputed_switching_search import PrecomputedSwitchingSearch

logger = setup_logging(level='INFO', name='precomputed_vs_sado')


class PerformanceComparator3000AA:

    
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    
    def generate_test_sequence(self, length: int) -> str:
        """生成测试序列"""

        amino_acids = list('ACDEFGHIKLMNPQRSTVWY')

        weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
        weights = weights / weights.sum()
        
        return ''.join(np.random.choice(amino_acids, length, p=weights))
    
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

                if valid_mask[pos, 0]:
                    pi[pos, 0] *= 1.5
                pi[pos] = pi[pos] / pi[pos].sum()
        
        return pi
    
    def test_methods_comparison(self):
        """对比测试不同方法"""
        logger.info("=" * 80)
        logger.info("预计算切换搜索 vs SADO (3000氨基酸)")
        logger.info("=" * 80)
        

        logger.info("生成3000氨基酸测试序列...")
        test_sequence = self.generate_test_sequence(3000)
        pi_distribution = self.create_test_distribution(3000, test_sequence)
        
        logger.info(f"测试序列长度: {len(test_sequence)}")
        logger.info(f"概率分布形状: {pi_distribution.shape}")
        
        results = {}
        

        logger.info("\n🔍 测试预计算切换搜索...")
        try:
            precomputed = PrecomputedSwitchingSearch(test_sequence, target_cai=0.8)
            
            start_time = time.time()
            precomputed_result = precomputed.precomputed_search(pi_distribution)
            precomputed_time = time.time() - start_time
            
            logger.info(f"✅ 预计算搜索完成: {precomputed_time*1000:.1f}ms")
            logger.info(f"   CAI: {precomputed_result['cai']:.4f}")
            logger.info(f"   log(P): {precomputed_result['log_prob']:.2f}")
            logger.info(f"   切换事件: {precomputed_result['operations']['switching_events_computed']}")
            logger.info(f"   搜索的γ值: {precomputed_result['searched_gammas']}")
            
            results['precomputed'] = {
                'time_ms': precomputed_time * 1000,
                'cai': precomputed_result['cai'],
                'log_prob': precomputed_result['log_prob'],
                'constraint_satisfied': precomputed_result['cai'] >= 0.8,
                'operations': precomputed_result['operations'],
                'method': 'PrecomputedSwitchingSearch'
            }
            
        except Exception as e:
            logger.error(f"❌ 预计算搜索失败: {e}")
            results['precomputed'] = {
                'time_ms': float('inf'),
                'error': str(e)
            }
        

        logger.info("\n🚀 测试SADO算法...")
        try:
            sado_optimizer = SADOOptimizer(
                species='ecoli_bl21de3',
                device=self.device,
                amino_acid_sequence=test_sequence
            )
            
            start_time = time.time()
            _, sado_metadata = sado_optimizer.optimize(
                pi_accessibility=pi_distribution,
                target_cai=0.8,
                amino_acid_sequence=test_sequence,
                gamma=0.0
            )
            sado_time = time.time() - start_time
            

            log_prob = -float('inf')
            if sado_optimizer.last_indices is not None:
                log_prob = 0.0
                for pos, idx in enumerate(sado_optimizer.last_indices):
                    if pos < len(sado_optimizer.codon_choices) and idx < len(sado_optimizer.codon_choices[pos]):
                        choice = sado_optimizer.codon_choices[pos][idx]
                        orig_idx = choice.get('original_local_index', idx)
                        if orig_idx < pi_distribution.shape[1]:
                            p = pi_distribution[pos, orig_idx].item()
                            if p > 0:
                                log_prob += np.log(p)
                            else:
                                log_prob = -float('inf')
                                break
            
            logger.info(f"✅ SADO完成: {sado_time*1000:.1f}ms")
            logger.info(f"   CAI: {sado_metadata['final_cai']:.4f}")
            logger.info(f"   log(P): {log_prob:.2f}")
            logger.info(f"   唯一序列: {sado_metadata['unique_sequences']}")
            
            results['sado'] = {
                'time_ms': sado_time * 1000,
                'cai': sado_metadata['final_cai'],
                'log_prob': log_prob,
                'constraint_satisfied': sado_metadata['constraint_satisfied'],
                'unique_sequences': sado_metadata['unique_sequences'],
                'method': 'SADO'
            }
            
        except Exception as e:
            logger.error(f"❌ SADO失败: {e}")
            results['sado'] = {
                'time_ms': float('inf'),
                'error': str(e)
            }
        

        logger.info("\n🎯 测试SADO + Gamma初始化...")
        try:
            sado_optimizer_gamma = SADOOptimizer(
                species='ecoli_bl21de3',
                device=self.device,
                amino_acid_sequence=test_sequence
            )
            
            start_time = time.time()
            _, sado_gamma_metadata = sado_optimizer_gamma.optimize(
                pi_accessibility=pi_distribution,
                target_cai=0.8,
                amino_acid_sequence=test_sequence,
                gamma=0.3
            )
            sado_gamma_time = time.time() - start_time
            

            log_prob_gamma = -float('inf')
            if sado_optimizer_gamma.last_indices is not None:
                log_prob_gamma = 0.0
                for pos, idx in enumerate(sado_optimizer_gamma.last_indices):
                    if pos < len(sado_optimizer_gamma.codon_choices) and idx < len(sado_optimizer_gamma.codon_choices[pos]):
                        choice = sado_optimizer_gamma.codon_choices[pos][idx]
                        orig_idx = choice.get('original_local_index', idx)
                        if orig_idx < pi_distribution.shape[1]:
                            p = pi_distribution[pos, orig_idx].item()
                            if p > 0:
                                log_prob_gamma += np.log(p)
                            else:
                                log_prob_gamma = -float('inf')
                                break
            
            logger.info(f"✅ SADO+Gamma完成: {sado_gamma_time*1000:.1f}ms")
            logger.info(f"   CAI: {sado_gamma_metadata['final_cai']:.4f}")
            logger.info(f"   log(P): {log_prob_gamma:.2f}")
            logger.info(f"   Gamma: {sado_gamma_metadata['gamma']}")
            
            results['sado_gamma'] = {
                'time_ms': sado_gamma_time * 1000,
                'cai': sado_gamma_metadata['final_cai'],
                'log_prob': log_prob_gamma,
                'constraint_satisfied': sado_gamma_metadata['constraint_satisfied'],
                'gamma': sado_gamma_metadata['gamma'],
                'method': 'SADO+Gamma'
            }
            
        except Exception as e:
            logger.error(f"❌ SADO+Gamma失败: {e}")
            results['sado_gamma'] = {
                'time_ms': float('inf'),
                'error': str(e)
            }
        

        self.analyze_results(results)
        
        return results
    
    def analyze_results(self, results: Dict):

        logger.info("=" * 80)

        logger.info("=" * 80)
        


        logger.info("-" * 70)
        
        valid_results = {}
        
        for method_name, result in results.items():
            if 'error' not in result:
                time_str = f"{result['time_ms']:.1f}" if result['time_ms'] != float('inf') else "timeout"
                cai_str = f"{result['cai']:.4f}"
                log_prob_str = f"{result['log_prob']:.2f}" if result['log_prob'] != -float('inf') else "-∞"
                satisfied = "✅" if result['constraint_satisfied'] else "❌"

                
                valid_results[method_name] = result
            else:

                cai_str = "N/A"
                log_prob_str = "N/A"
                satisfied = "❌"

            
            logger.info(f"{method_name:<20} {time_str:<12} {cai_str:<8} {log_prob_str:<10} {satisfied:<8} {status}")
        

        if len(valid_results) >= 2:

            

            fastest_method = min(valid_results.keys(), key=lambda k: valid_results[k]['time_ms'])
            fastest_time = valid_results[fastest_method]['time_ms']
            

            

            for method_name, result in valid_results.items():
                if method_name != fastest_method:
                    speedup = result['time_ms'] / fastest_time

            

            best_cai_method = max(valid_results.keys(), key=lambda k: valid_results[k]['cai'])
            best_prob_method = max(valid_results.keys(), key=lambda k: valid_results[k]['log_prob'] if valid_results[k]['log_prob'] != -float('inf') else -1e10)
            



            





            
            if 'precomputed' in valid_results:
                precomputed_ops = valid_results['precomputed']['operations']['switching_events_computed']


            





        
        if len([r for r in valid_results.values() if r['constraint_satisfied']]) > 0:

        
        if 'sado_gamma' in valid_results and 'sado' in valid_results:
            gamma_improvement = valid_results['sado']['time_ms'] / valid_results['sado_gamma']['time_ms']



def main():
    """主函数"""
    logger.info("🚀 开始3000氨基酸性能对比测试")
    
    comparator = PerformanceComparator3000AA()
    results = comparator.test_methods_comparison()
    
    logger.info("✅ 对比测试完成!")
    
    return results


if __name__ == "__main__":
    main()