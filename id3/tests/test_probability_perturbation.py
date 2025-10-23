#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
import hashlib
from typing import Dict, List, Tuple

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import BinarySearchCAIOptimizer

logger = setup_logging(level='INFO', name='prob_perturbation')


class ProbabilityPerturbationOptimizer:

    
    def __init__(self, sequence: str, target_cai: float = 0.8):
        self.sequence = sequence
        self.target_cai = target_cai
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self.base_optimizer = BinarySearchCAIOptimizer(
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=sequence
        )
        

        self.history_hashes = set()
        self.unique_sequences = 0
        
    def perturb_probability(self, 
                           pi_original: torch.Tensor, 
                           perturbation_strength: float = 0.1,
                           positions_ratio: float = 0.1) -> torch.Tensor:
        """
        对概率分布进行轻微扰动
        
        Args:
            pi_original: 原始概率分布
            perturbation_strength: 扰动强度 (0-1)
            positions_ratio: 扰动位置的比例
        
        Returns:
            扰动后的概率分布
        """
        pi_perturbed = pi_original.clone()
        seq_len = pi_perturbed.shape[0]
        

        num_positions = max(1, int(seq_len * positions_ratio))
        positions = np.random.choice(seq_len, num_positions, replace=False)
        
        for pos in positions:

            current_probs = pi_perturbed[pos]
            

            noise = torch.randn_like(current_probs) * perturbation_strength
            

            new_probs = current_probs + noise
            new_probs = torch.clamp(new_probs, min=0.0)
            

            if new_probs.sum() > 0:
                new_probs = new_probs / new_probs.sum()
            else:

            
            pi_perturbed[pos] = new_probs
        
        return pi_perturbed
    
    def optimize_with_perturbation(self, 
                                  pi_original: torch.Tensor,
                                  max_attempts: int = 10) -> Dict:
        """
        使用概率扰动策略优化
        
        1. 对原始概率进行轻微扰动
        2. 用二分查找找到满足CAI≥0.8的序列
        3. 如果重复，增加扰动强度重试
        """
        best_result = None

        
        for attempt in range(max_attempts):

            if attempt == 0:

                pi_perturbed = pi_original
            else:

                current_strength = perturbation_strength * (1 + attempt * 0.5)

                


                
                pi_perturbed = self.perturb_probability(
                    pi_original, 
                    perturbation_strength=current_strength,
                    positions_ratio=positions_ratio
                )
            

            result, metadata = self.base_optimizer.optimize(
                pi_accessibility=pi_perturbed,
                target_cai=self.target_cai,
                amino_acid_sequence=self.sequence
            )
            

            indices = result.argmax(dim=-1).cpu().numpy()
            seq_hash = hashlib.md5(indices.tobytes()).hexdigest()
            

            if seq_hash not in self.history_hashes:

                self.history_hashes.add(seq_hash)
                self.unique_sequences += 1
                
                return {
                    'sequence': result,
                    'indices': indices,
                    'cai': metadata['final_cai'],
                    'gamma': metadata['gamma'],
                    'constraint_satisfied': metadata['constraint_satisfied'],
                    'perturbation_strength': current_strength if attempt > 0 else 0,
                    'attempt': attempt + 1,
                    'unique_sequences': self.unique_sequences,
                    'is_duplicate': False
                }
            

            if best_result is None or metadata['final_cai'] > best_result['cai']:
                best_result = {
                    'sequence': result,
                    'indices': indices,
                    'cai': metadata['final_cai'],
                    'gamma': metadata['gamma'],
                    'constraint_satisfied': metadata['constraint_satisfied'],
                    'perturbation_strength': current_strength if attempt > 0 else 0,
                    'attempt': attempt + 1,
                    'unique_sequences': self.unique_sequences,
                    'is_duplicate': True
                }
        

        return best_result


def test_perturbation_strategy():
    """测试概率扰动策略的效果"""
    logger.info("="*80)
    logger.info("概率扰动策略测试")
    logger.info("="*80)
    

    seq_length = 500
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, seq_length, p=weights))
    
    logger.info(f"测试序列长度: {seq_length} aa")
    logger.info(f"目标CAI: 0.8")
    

    optimizer = ProbabilityPerturbationOptimizer(test_sequence, target_cai=0.8)
    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    pi_original = torch.rand(seq_length, 6, device=device)
    for pos in range(seq_length):
        if pi_original[pos].sum() > 0:
            pi_original[pos] = pi_original[pos] / pi_original[pos].sum()
    

    num_iterations = 20
    results = []
    
    logger.info(f"\n开始{num_iterations}次迭代测试...")
    logger.info("-"*60)
    
    for i in range(num_iterations):

        if i % 5 == 0 and i > 0:

            pi_original = torch.rand(seq_length, 6, device=device)
            for pos in range(seq_length):
                if pi_original[pos].sum() > 0:
                    pi_original[pos] = pi_original[pos] / pi_original[pos].sum()
            logger.info(f"  更换基础概率分布")
        

        result = optimizer.optimize_with_perturbation(pi_original)
        results.append(result)
        
        logger.info(f"迭代 {i+1}:")
        logger.info(f"  CAI: {result['cai']:.4f} {'✅' if result['constraint_satisfied'] else '❌'}")
        logger.info(f"  扰动强度: {result['perturbation_strength']*100:.1f}%")
        logger.info(f"  尝试次数: {result['attempt']}")
        logger.info(f"  是否重复: {'是' if result['is_duplicate'] else '否'}")
        logger.info(f"  唯一序列总数: {result['unique_sequences']}")
    

    logger.info("\n" + "="*80)
    logger.info("效果分析")
    logger.info("="*80)
    

    unique_count = optimizer.unique_sequences
    duplicate_count = sum(1 for r in results if r['is_duplicate'])
    
    logger.info(f"\n📊 多样性:")
    logger.info(f"  唯一序列: {unique_count}/{num_iterations}")
    logger.info(f"  重复次数: {duplicate_count}")
    logger.info(f"  多样性率: {unique_count/num_iterations*100:.1f}%")
    

    cais = [r['cai'] for r in results]
    satisfied = sum(1 for r in results if r['constraint_satisfied'])
    
    logger.info(f"\n📈 CAI性能:")
    logger.info(f"  约束满足率: {satisfied}/{num_iterations} ({satisfied/num_iterations*100:.1f}%)")
    logger.info(f"  平均CAI: {np.mean(cais):.4f}")
    logger.info(f"  CAI范围: [{min(cais):.4f}, {max(cais):.4f}]")
    

    perturbations = [r['perturbation_strength'] for r in results if r['perturbation_strength'] > 0]
    if perturbations:
        logger.info(f"\n🎯 扰动策略:")
        logger.info(f"  平均扰动强度: {np.mean(perturbations)*100:.1f}%")
        logger.info(f"  最大扰动强度: {max(perturbations)*100:.1f}%")
    

    attempts = [r['attempt'] for r in results]
    logger.info(f"\n⚡ 效率:")
    logger.info(f"  平均尝试次数: {np.mean(attempts):.1f}")
    logger.info(f"  最多尝试次数: {max(attempts)}")
    

    logger.info(f"\n💡 对比分析:")
    logger.info(f"  原始二分查找: 重复率100%")
    logger.info(f"  概率扰动策略: 重复率{duplicate_count/num_iterations*100:.1f}%")
    logger.info(f"  改进效果: {(1-duplicate_count/num_iterations)*100:.1f}%多样性提升")
    
    if satisfied == num_iterations:
        logger.info(f"  ✅ 完美！所有序列都满足CAI≥0.8约束")
    
    return results


def test_perturbation_strength_impact():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    seq_length = 300
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, seq_length, p=weights))
    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    pi_original = torch.rand(seq_length, 6, device=device)
    for pos in range(seq_length):
        if pi_original[pos].sum() > 0:
            pi_original[pos] = pi_original[pos] / pi_original[pos].sum()
    

    strengths = [0.01, 0.05, 0.1, 0.15, 0.2, 0.3]
    

    logger.info("-"*60)
    
    for strength in strengths:
        optimizer = ProbabilityPerturbationOptimizer(test_sequence, target_cai=0.8)
        

        pi_perturbed = optimizer.perturb_probability(
            pi_original, 
            perturbation_strength=strength,
            positions_ratio=0.2
        )
        
        result, metadata = optimizer.base_optimizer.optimize(
            pi_accessibility=pi_perturbed,
            target_cai=0.8,
            amino_acid_sequence=test_sequence
        )
        

        result_original, _ = optimizer.base_optimizer.optimize(
            pi_accessibility=pi_original,
            target_cai=0.8,
            amino_acid_sequence=test_sequence
        )
        
        indices_perturbed = result.argmax(dim=-1).cpu().numpy()
        indices_original = result_original.argmax(dim=-1).cpu().numpy()
        
        difference = np.sum(indices_perturbed != indices_original)
        

        logger.info(f"  CAI: {metadata['final_cai']:.4f}")


    






if __name__ == "__main__":

    

    results = test_perturbation_strategy()
    

    test_perturbation_strength_impact()
    

