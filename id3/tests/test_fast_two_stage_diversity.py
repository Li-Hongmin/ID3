#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
import hashlib
import time
from typing import Dict, Tuple, Optional, Set

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import BinarySearchCAIOptimizer

logger = setup_logging(level='INFO', name='fast_two_stage')


class FastTwoStageOptimizer:

    
    def __init__(self, sequence: str, target_cai: float = 0.8):
        self.sequence = sequence
        self.target_cai = target_cai
        self.seq_len = len(sequence)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self.binary_optimizer = BinarySearchCAIOptimizer(
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=sequence
        )
        

        self.wi_table = self.binary_optimizer.wi_table
        

        self._build_codon_info()
        

        self.history_hashes = set()
        self.unique_sequences = 0
        
    def _build_codon_info(self):
        """构建密码子选择信息"""
        self.codon_choices = []
        self.swappable_positions = []
        
        for pos, aa in enumerate(self.sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                choices = []
                for i, codon in enumerate(codons):
                    choices.append({
                        'codon': codon,
                        'index': i,
                        'weight': self.wi_table.get(codon, 0.0)
                    })

                choices.sort(key=lambda x: x['weight'], reverse=True)
                self.codon_choices.append(choices)
                

                if len(choices) > 1:
                    self.swappable_positions.append(pos)
            else:
                self.codon_choices.append([])
    
    def stage1_binary_search(self, pi_accessibility: torch.Tensor) -> Tuple[np.ndarray, float, float]:
        """

        """

        valid_mask = torch.zeros(self.seq_len, 6, dtype=torch.bool, device=self.device)
        for pos, aa in enumerate(self.sequence):
            if aa in amino_acids_to_codons:
                num_codons = len(amino_acids_to_codons[aa])
                for i in range(min(num_codons, 6)):
                    valid_mask[pos, i] = True
        
        result, metadata = self.binary_optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=self.target_cai,
            amino_acid_sequence=self.sequence,
            valid_codon_mask=valid_mask
        )
        

        indices = result.argmax(dim=-1).cpu().numpy()
        
        return indices, metadata['final_cai'], metadata['optimal_gamma']
    
    def stage2_simplest_diversity(self, initial_indices: np.ndarray, 
                                  max_attempts: int = 5) -> Tuple[np.ndarray, float, bool]:
        """


        """

        initial_hash = hashlib.md5(initial_indices.tobytes()).hexdigest()
        if initial_hash not in self.history_hashes:
            initial_cai = self._calculate_cai(initial_indices)
            return initial_indices, initial_cai, True
        

        if not self.swappable_positions:
            initial_cai = self._calculate_cai(initial_indices)
            return initial_indices, initial_cai, False
        

        for attempt in range(max_attempts):
            new_indices = initial_indices.copy()
            

            num_swaps = min(3, len(self.swappable_positions))
            if attempt > 2:

                num_swaps = min(5, len(self.swappable_positions))
            
            positions = np.random.choice(
                self.swappable_positions, 
                min(num_swaps, len(self.swappable_positions)), 
                replace=False
            )
            
            for pos in positions:
                choices = self.codon_choices[pos]
                current_idx = new_indices[pos]
                


                if len(choices) > 1:
                    if choices[0]['index'] == current_idx:

                        new_indices[pos] = choices[1]['index']
                    else:

                        if np.random.random() < 0.7:

                            new_indices[pos] = choices[0]['index']
                        else:

                            available = [c['index'] for c in choices if c['index'] != current_idx]
                            if available:
                                new_indices[pos] = np.random.choice(available)
            

            new_hash = hashlib.md5(new_indices.tobytes()).hexdigest()
            if new_hash not in self.history_hashes:

                new_cai = self._calculate_cai(new_indices)
                if new_cai >= self.target_cai:

                    return new_indices, new_cai, True
        

        initial_cai = self._calculate_cai(initial_indices)
        return initial_indices, initial_cai, False
    
    def _calculate_cai(self, indices: np.ndarray) -> float:

        log_sum = 0.0
        count = 0
        
        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                if idx < len(self.codon_choices[pos]):
                    weight = self.codon_choices[pos][idx]['weight']
                    if weight > 0:
                        log_sum += np.log(weight)
                        count += 1
        
        if count == 0:
            return 0.0
        
        return np.exp(log_sum / count)
    
    def optimize_with_diversity(self, pi_accessibility: torch.Tensor) -> Dict:
        """
        主优化函数：保证多样性
        """
        start_time = time.time()
        

        stage1_start = time.time()
        indices_bs, cai_bs, gamma = self.stage1_binary_search(pi_accessibility)
        stage1_time = time.time() - stage1_start
        

        if cai_bs < self.target_cai:
            logger.debug(f"Warning: Binary search CAI={cai_bs:.4f} < target={self.target_cai}, gamma={gamma:.4f}")
        

        stage2_start = time.time()
        final_indices, final_cai, is_new = self.stage2_simplest_diversity(indices_bs)
        stage2_time = time.time() - stage2_start
        

        if is_new:
            seq_hash = hashlib.md5(final_indices.tobytes()).hexdigest()
            self.history_hashes.add(seq_hash)
            self.unique_sequences += 1
        
        total_time = time.time() - start_time
        
        return {
            'indices': final_indices,
            'cai': final_cai,
            'gamma': gamma,
            'is_new': is_new,
            'unique_count': self.unique_sequences,
            'stage1_time': stage1_time * 1000,  # ms
            'stage2_time': stage2_time * 1000,  # ms
            'total_time': total_time * 1000,    # ms
        }


def test_diversity_performance():
    """测试多样性和性能"""
    logger.info("="*80)
    logger.info("快速两阶段优化测试")
    logger.info("="*80)
    

    seq_length = 500
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, seq_length, p=weights))
    
    logger.info(f"测试序列: {seq_length} aa")
    logger.info(f"目标CAI: 0.8")
    

    optimizer = FastTwoStageOptimizer(test_sequence, target_cai=0.8)
    logger.info(f"可交换位置数: {len(optimizer.swappable_positions)}")
    logger.info(f"最大可达CAI: {optimizer.binary_optimizer.max_achievable_cai:.4f}")
    

    if optimizer.binary_optimizer.max_achievable_cai < 0.8:
        new_target = optimizer.binary_optimizer.max_achievable_cai * 0.9
        logger.info(f"调整目标CAI: {0.8:.4f} -> {new_target:.4f}")
        optimizer.target_cai = new_target
    

    num_iterations = 100
    results = []
    
    logger.info(f"\n运行 {num_iterations} 次迭代...")
    logger.info("-"*60)
    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    for i in range(num_iterations):

        pi_accessibility = torch.rand(seq_length, 6, device=device)
        for pos in range(seq_length):
            if pi_accessibility[pos].sum() > 0:
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        result = optimizer.optimize_with_diversity(pi_accessibility)
        results.append(result)
        

        if (i + 1) % 20 == 0:
            unique_rate = optimizer.unique_sequences / (i + 1) * 100
            avg_time = np.mean([r['total_time'] for r in results[-20:]])
            logger.info(f"迭代 {i+1}: 唯一序列 {optimizer.unique_sequences}, "
                       f"唯一率 {unique_rate:.1f}%, "
                       f"平均时间 {avg_time:.1f}ms")
    

    logger.info("\n" + "="*80)
    logger.info("性能分析")
    logger.info("="*80)
    

    unique_count = optimizer.unique_sequences
    duplicate_count = sum(1 for r in results if not r['is_new'])
    
    logger.info(f"\n📊 多样性:")
    logger.info(f"  唯一序列: {unique_count}/{num_iterations}")
    logger.info(f"  重复次数: {duplicate_count}")
    logger.info(f"  多样性率: {unique_count/num_iterations*100:.1f}%")
    

    cais = [r['cai'] for r in results]
    constraint_satisfied = sum(1 for cai in cais if cai >= 0.8)
    
    logger.info(f"\n📈 CAI性能:")
    logger.info(f"  约束满足率: {constraint_satisfied}/{num_iterations} "
                f"({constraint_satisfied/num_iterations*100:.1f}%)")
    logger.info(f"  平均CAI: {np.mean(cais):.4f}")
    logger.info(f"  CAI范围: [{min(cais):.4f}, {max(cais):.4f}]")
    

    stage1_times = [r['stage1_time'] for r in results]
    stage2_times = [r['stage2_time'] for r in results]
    total_times = [r['total_time'] for r in results]
    
    logger.info(f"\n⚡ 速度性能:")
    logger.info(f"  Stage1平均: {np.mean(stage1_times):.1f}ms")
    logger.info(f"  Stage2平均: {np.mean(stage2_times):.1f}ms")
    logger.info(f"  总平均时间: {np.mean(total_times):.1f}ms")
    logger.info(f"  最快: {min(total_times):.1f}ms")
    logger.info(f"  最慢: {max(total_times):.1f}ms")
    

    if unique_count > 0:
        estimated_1000 = min(1000, unique_count * 1000 / num_iterations)
        estimated_time = np.mean(total_times) * 1000 / 1000
        logger.info(f"\n💡 1000次迭代预测:")
        logger.info(f"  预计唯一序列: ~{estimated_1000:.0f}")
        logger.info(f"  预计总时间: ~{estimated_time:.1f}秒")
    

    logger.info(f"\n🎯 结论:")
    if unique_count / num_iterations > 0.9:
        logger.info(f"  ✅ 多样性优秀 (>90%唯一)")
    elif unique_count / num_iterations > 0.7:
        logger.info(f"  ✅ 多样性良好 (>70%唯一)")
    else:
        logger.info(f"  ⚠️ 多样性不足，需要更多策略")
    
    if np.mean(total_times) < 20:
        logger.info(f"  ✅ 速度优秀 (<20ms/iteration)")
    elif np.mean(total_times) < 50:
        logger.info(f"  ✅ 速度良好 (<50ms/iteration)")
    else:
        logger.info(f"  ⚠️ 速度较慢")
    
    if constraint_satisfied == num_iterations:
        logger.info(f"  ✅ CAI约束完美满足")
    
    return results


def test_real_protein():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    data_dir = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/data')
    spike_file = data_dir / 'P0DTC2.fasta'
    
    with open(spike_file, 'r') as f:
        lines = f.readlines()
    
    sequence = ''
    for line in lines:
        if not line.startswith('>'):
            sequence += line.strip()
    
    logger.info(f"SARS-CoV-2 Spike: {len(sequence)} aa")
    

    optimizer = FastTwoStageOptimizer(sequence, target_cai=0.8)

    

    num_iterations = 20
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    pi_accessibility = torch.rand(len(sequence), 6, device=device)
    for pos in range(len(sequence)):
        if pi_accessibility[pos].sum() > 0:
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    
    times = []
    for i in range(num_iterations):
        start = time.time()
        result = optimizer.optimize_with_diversity(pi_accessibility)
        times.append((time.time() - start) * 1000)
        
        if i % 5 == 0:



    






if __name__ == "__main__":

    

    results = test_diversity_performance()
    

    test_real_protein()
    

