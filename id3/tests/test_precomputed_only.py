#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from typing import List

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/experiments/cai_enhancement_theory/04_method_comparison')
from precomputed_switching_search import PrecomputedSwitchingSearch

logger = setup_logging(level='INFO', name='precomputed_test')


class PrecomputedTester:

    
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
                pi[pos] = pi[pos] / pi[pos].sum()
        
        return pi
    
    def test_scaling_performance(self):
        """测试不同序列长度的性能表现"""
        logger.info("=" * 80)
        logger.info("预计算切换搜索性能测试")
        logger.info("=" * 80)
        

        test_lengths = [100, 300, 500, 1000, 2000, 3000]
        results = []
        
        for length in test_lengths:
            logger.info(f"\n📊 测试序列长度: {length} 氨基酸")
            logger.info("-" * 50)
            

            test_sequence = self.generate_test_sequence(length)
            pi_distribution = self.create_test_distribution(length, test_sequence)
            
            try:

                precomputed = PrecomputedSwitchingSearch(test_sequence, target_cai=0.8)
                

                result = precomputed.precomputed_search(pi_distribution)
                

                result_summary = {
                    'length': length,
                    'time_ms': result['time'] * 1000,
                    'cai': result['cai'],
                    'log_prob': result['log_prob'],
                    'constraint_satisfied': result['cai'] >= 0.8,
                    'switching_events': result['operations']['switching_events_computed'],
                    'total_gammas': result['total_gammas'],
                    'searched_gammas': result['searched_gammas'],
                    'sequences_generated': result['operations']['sequences_generated']
                }
                
                results.append(result_summary)
                
                logger.info(f"✅ 完成: {result_summary['time_ms']:.1f}ms")
                logger.info(f"   CAI: {result_summary['cai']:.4f}")
                logger.info(f"   约束满足: {'✅' if result_summary['constraint_satisfied'] else '❌'}")
                logger.info(f"   切换事件数: {result_summary['switching_events']}")
                logger.info(f"   搜索序列数: {result_summary['searched_gammas']}")
                
            except Exception as e:
                logger.error(f"❌ 长度 {length} 测试失败: {e}")
                results.append({
                    'length': length,
                    'time_ms': float('inf'),
                    'error': str(e)
                })
        

        self.analyze_scaling_results(results)
        
        return results
    
    def analyze_scaling_results(self, results: List[dict]):

        logger.info("=" * 80)

        logger.info("=" * 80)
        

        valid_results = [r for r in results if 'error' not in r]
        
        if len(valid_results) < 2:

            return
        


        logger.info("-" * 70)
        
        for result in valid_results:
            length = result['length']
            time_ms = result['time_ms']
            cai = result['cai']
            satisfied = "✅" if result['constraint_satisfied'] else "❌"
            switching_events = result['switching_events']
            


            theoretical = length * k * np.log(length * k)
            ratio = switching_events / theoretical if theoretical > 0 else 0
            
            logger.info(f"{length:<6} {time_ms:<10.1f} {cai:<8.4f} {satisfied:<4} "
                       f"{switching_events:<8} {theoretical:<12.0f} {ratio:<.3f}")
        


        

        for i in range(1, len(valid_results)):
            prev_result = valid_results[i-1]
            curr_result = valid_results[i]
            
            length_ratio = curr_result['length'] / prev_result['length']
            time_ratio = curr_result['time_ms'] / prev_result['time_ms']
            

            theoretical_ratio = length_ratio * np.log(curr_result['length']) / np.log(prev_result['length'])
            
            logger.info(f"  {prev_result['length']} → {curr_result['length']}: "

        


        

        result_3000 = next((r for r in valid_results if r['length'] == 3000), None)
        if result_3000:
            time_3000 = result_3000['time_ms']
            




            else:

            




        

        satisfied_count = sum(1 for r in valid_results if r['constraint_satisfied'])
        satisfaction_rate = satisfied_count / len(valid_results)
        


        
        avg_cai = np.mean([r['cai'] for r in valid_results])

        
        if satisfaction_rate >= 0.9:

        elif satisfaction_rate >= 0.7:

        else:

        





        
        if result_3000 and result_3000['time_ms'] < 10000:



def main():
    """主函数"""
    logger.info("🚀 开始预计算切换搜索性能测试")
    
    tester = PrecomputedTester()
    results = tester.test_scaling_performance()
    
    logger.info("✅ 测试完成!")
    
    return results


if __name__ == "__main__":
    main()