#!/usr/bin/env python3
"""







"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np
import pandas as pd
from typing import List, Dict
import matplotlib.pyplot as plt

from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging

logger = setup_logging(level='INFO', name='sado_long_iter')

class SADOLongSequenceIterationTester:
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        logger.info(f"设备: {self.device}")
    
    def compute_sequence_probability(self, indices, pi_distribution, optimizer):

        log_prob = 0.0
        
        for pos, idx in enumerate(indices):
            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):
                choice = optimizer.codon_choices[pos][idx]
                orig_idx = choice.get('original_local_index', idx)
                
                if orig_idx < pi_distribution.shape[1]:
                    p = pi_distribution[pos, orig_idx].item()
                    if p > 0:
                        log_prob += np.log(p)
                    else:
                        return float('-inf')
        
        return log_prob
    
    def generate_random_amino_sequence(self, length: int) -> str:
        """生成随机氨基酸序列"""
        amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
        return ''.join(np.random.choice(amino_acids, length))
    
    def create_random_distribution(self, seq_length: int, num_codons: int, amino_sequence: str):


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
    
    def test_long_sequence_iteration(self):
        """长序列迭代测试"""
        logger.info("=" * 80)
        logger.info("SADO长序列迭代性能测试")
        logger.info("=" * 80)
        

        seq_length = 3000
        num_iterations = 1000
        num_codons = 6
        target_cai = 0.8
        
        logger.info(f"序列长度: {seq_length} 氨基酸")
        logger.info(f"迭代次数: {num_iterations}")
        logger.info(f"目标CAI: {target_cai}")
        

        logger.info("生成3000个氨基酸的随机序列...")
        test_sequence = self.generate_random_amino_sequence(seq_length)
        

        logger.info("初始化SADO优化器...")
        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=test_sequence
        )
        

        results = {
            'iteration': [],
            'execution_time_ms': [],
            'final_cai': [],
            'constraint_satisfied': [],
            'log_probability': [],
            'unique_sequences': [],
            'is_first_run': []
        }
        
        logger.info("开始1000次迭代测试...")
        logger.info("-" * 80)
        

        base_seed = 42
        
        for i in range(num_iterations):

            torch.manual_seed(base_seed + i)
            np.random.seed(base_seed + i)
            

            pi_distribution = self.create_random_distribution(seq_length, num_codons, test_sequence)
            

            start_time = time.time()
            

            _, metadata = optimizer.optimize(
                pi_accessibility=pi_distribution,
                target_cai=target_cai,
                amino_acid_sequence=test_sequence
            )
            
            end_time = time.time()
            execution_time = (end_time - start_time) * 1000  # ms
            

            log_prob = float('-inf')
            if optimizer.last_indices is not None:
                log_prob = self.compute_sequence_probability(
                    optimizer.last_indices, pi_distribution, optimizer
                )
            

            results['iteration'].append(i + 1)
            results['execution_time_ms'].append(execution_time)
            results['final_cai'].append(metadata['final_cai'])
            results['constraint_satisfied'].append(metadata['constraint_satisfied'])
            results['log_probability'].append(log_prob)
            results['unique_sequences'].append(metadata['unique_sequences'])
            results['is_first_run'].append(i == 0)
            

            if (i + 1) % 100 == 0:
                avg_time = np.mean(results['execution_time_ms'][-100:])
                avg_cai = np.mean(results['final_cai'][-100:])
                satisfaction_rate = sum(results['constraint_satisfied'][-100:])
                logger.info(f"迭代 {i+1:4d}: 平均用时={avg_time:6.1f}ms, "
                           f"平均CAI={avg_cai:.4f}, 满足率={satisfaction_rate}/100")
        

        df = pd.DataFrame(results)
        

        self.analyze_iteration_results(df)
        
        return df
    
    def analyze_iteration_results(self, df: pd.DataFrame):

        logger.info("=" * 80)

        logger.info("=" * 80)
        

        total_iterations = len(df)
        first_run_time = df[df['is_first_run']]['execution_time_ms'].iloc[0]
        subsequent_times = df[~df['is_first_run']]['execution_time_ms']
        



        

        speedup = first_run_time / subsequent_times.mean()

        

        satisfaction_rate = df['constraint_satisfied'].sum() / total_iterations
        avg_cai = df['final_cai'].mean()


        

        valid_probs = df[df['log_probability'] != float('-inf')]['log_probability']
        if len(valid_probs) > 0:


        

        final_unique = df['unique_sequences'].iloc[-1]


        


        

        stages = [



        ]
        
        for stage_name, stage_data in stages:
            avg_time = stage_data['execution_time_ms'].mean()
            avg_cai = stage_data['final_cai'].mean()
            satisfaction = stage_data['constraint_satisfied'].sum()

        


        




        else:

        

        
        if satisfaction_rate >= 0.95:

        elif satisfaction_rate >= 0.8:

        else:

        

        

        if speedup > 10:

        elif speedup > 5:

        elif speedup > 2:

        else:

        


def main():
    """主函数"""
    logger.info("🚀 开始SADO长序列迭代性能测试")
    
    tester = SADOLongSequenceIterationTester()
    results_df = tester.test_long_sequence_iteration()
    
    logger.info("✅ 测试完成!")
    logger.info(f"详细结果已保存，共 {len(results_df)} 次迭代数据")
    
    return results_df

if __name__ == "__main__":
    main()