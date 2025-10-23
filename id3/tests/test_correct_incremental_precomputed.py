#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np
import hashlib
from typing import Dict, List, Tuple, Optional, Set

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import SADOOptimizer

logger = setup_logging(level='INFO', name='correct_incremental')


class CorrectIncrementalPrecomputed:

    
    def __init__(self, sequence: str, target_cai: float = 0.8):
        self.sequence = sequence
        self.target_cai = target_cai
        self.seq_len = len(sequence)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self._init_cai_weights()
        

        self.cached_switching_events = {}  # {position: [(gamma, from_idx, to_idx), ...]}
        self.cached_all_gammas = None
        self.cached_sequences = {}  # {gamma: sequence_indices}
        self.cached_cais = {}  # {gamma: cai_value}
        

        self.history_hashes: Set[str] = set()
        self.used_gammas: Set[float] = set()
        

        self.iteration_count = 0
        self.unique_sequences = 0
        
    def _init_cai_weights(self):
        """初始化CAI权重表和密码子信息"""
        sado = SADOOptimizer(
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=self.sequence
        )
        self.wi_table = sado.wi_table
        

        self.codon_choices = []
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
            else:
                self.codon_choices.append([])
        

        self.cai_optimal_indices = []
        for choices in self.codon_choices:
            if choices:

                self.cai_optimal_indices.append(choices[0]['index'])
            else:
                self.cai_optimal_indices.append(0)
    
    def _compute_switching_events_for_position(self, 
                                              pos: int, 
                                              target_dist: torch.Tensor) -> List[Tuple[float, int, int]]:
        """


        """
        events = []
        
        if not self.codon_choices[pos]:
            return events
        
        choices = self.codon_choices[pos]
        

        for i in range(len(choices)):
            for j in range(len(choices)):
                if i == j:
                    continue
                    
                codon_i = choices[i]
                codon_j = choices[j]
                

                p_i = target_dist[pos, codon_i['index']].item() if codon_i['index'] < target_dist.shape[1] else 0
                p_j = target_dist[pos, codon_j['index']].item() if codon_j['index'] < target_dist.shape[1] else 0
                w_i = codon_i['weight']
                w_j = codon_j['weight']
                

                # γ * w_i + (1 - γ) * p_i = γ * w_j + (1 - γ) * p_j
                # γ * (w_i - w_j) = (1 - γ) * (p_j - p_i)
                # γ = (p_j - p_i) / ((p_j - p_i) + (w_i - w_j))
                
                denom = (p_j - p_i) + (w_i - w_j)
                if abs(denom) > 1e-10:
                    gamma_switch = (p_j - p_i) / denom
                    if 0 <= gamma_switch <= 1:
                        events.append((gamma_switch, codon_i['index'], codon_j['index']))
        
        return events
    
    def _generate_sequence_at_gamma(self, gamma: float, target_dist: torch.Tensor) -> np.ndarray:
        """


        """
        sequence_indices = np.zeros(self.seq_len, dtype=int)
        
        for pos in range(self.seq_len):
            if not self.codon_choices[pos]:
                continue
            
            best_score = -float('inf')
            best_idx = 0
            
            for choice in self.codon_choices[pos]:
                idx = choice['index']
                

                cai_score = choice['weight']
                

                if idx < target_dist.shape[1]:
                    prob_score = target_dist[pos, idx].item()
                else:
                    prob_score = 0
                

                interpolated_score = gamma * cai_score + (1 - gamma) * prob_score
                
                if interpolated_score > best_score:
                    best_score = interpolated_score
                    best_idx = idx
            
            sequence_indices[pos] = best_idx
        
        return sequence_indices
    
    def _calculate_cai_from_indices(self, indices: np.ndarray) -> float:

        log_sum = 0.0
        count = 0
        
        for pos, aa in enumerate(self.sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                if indices[pos] < len(codons):
                    codon = codons[indices[pos]]
                    weight = self.wi_table.get(codon, 1e-10)
                    log_sum += np.log(weight)
                    count += 1
        
        if count == 0:
            return 0.0
        
        return np.exp(log_sum / count)
    
    def _compute_sequence_hash(self, indices: np.ndarray) -> str:
        """计算序列的MD5哈希"""
        return hashlib.md5(indices.tobytes()).hexdigest()
    
    def incremental_search_with_diversity(self, 
                                         target_dist: torch.Tensor,
                                         changed_positions: Optional[List[int]] = None) -> Dict:
        """

        """
        start_time = time.time()
        self.iteration_count += 1
        
        operations = {
            'changed_positions': 0,
            'new_switching_events': 0,
            'reused_switching_events': 0,
            'total_gammas': 0,
            'sequences_generated': 0,
            'cai_calculations': 0,
            'duplicates_avoided': 0
        }
        

        if changed_positions is None:

            changed_positions = list(range(self.seq_len))
            logger.info(f"  首次计算，处理所有{self.seq_len}个位置")
        else:
            logger.info(f"  增量更新，处理{len(changed_positions)}个变化位置")
        
        operations['changed_positions'] = len(changed_positions)
        

        for pos in changed_positions:
            events = self._compute_switching_events_for_position(pos, target_dist)
            self.cached_switching_events[pos] = events
            operations['new_switching_events'] += len(events)
        

        for pos in range(self.seq_len):
            if pos not in changed_positions and pos in self.cached_switching_events:
                operations['reused_switching_events'] += len(self.cached_switching_events[pos])
        

        all_gammas = set([0.0, 1.0])
        for events in self.cached_switching_events.values():
            for gamma, _, _ in events:
                all_gammas.add(gamma)
        
        all_gammas = sorted(all_gammas)
        operations['total_gammas'] = len(all_gammas)
        self.cached_all_gammas = all_gammas
        

        sequences_to_evaluate = []
        cais_to_evaluate = []
        
        for gamma in all_gammas:

            if gamma in self.cached_sequences and gamma not in self.used_gammas:

                sequence_indices = self.cached_sequences[gamma]
                cai = self.cached_cais[gamma]
            else:

                sequence_indices = self._generate_sequence_at_gamma(gamma, target_dist)
                operations['sequences_generated'] += 1
                

                seq_hash = self._compute_sequence_hash(sequence_indices)
                if seq_hash in self.history_hashes:
                    operations['duplicates_avoided'] += 1
                    continue
                

                cai = self._calculate_cai_from_indices(sequence_indices)
                operations['cai_calculations'] += 1
                

                self.cached_sequences[gamma] = sequence_indices
                self.cached_cais[gamma] = cai
            
            sequences_to_evaluate.append((gamma, sequence_indices, cai))
            cais_to_evaluate.append(cai)
        

        best_sequence = None
        best_gamma = None
        best_cai = 0
        

        sequences_to_evaluate.sort(key=lambda x: x[2])
        

        for gamma, sequence_indices, cai in sequences_to_evaluate:
            if cai >= self.target_cai:

                seq_hash = self._compute_sequence_hash(sequence_indices)
                if seq_hash not in self.history_hashes:
                    best_gamma = gamma
                    best_sequence = sequence_indices
                    best_cai = cai
                    

                    self.history_hashes.add(seq_hash)
                    self.used_gammas.add(gamma)
                    self.unique_sequences += 1
                    break
        

        if best_sequence is None:
            for gamma, sequence_indices, cai in reversed(sequences_to_evaluate):
                seq_hash = self._compute_sequence_hash(sequence_indices)
                if seq_hash not in self.history_hashes:
                    best_gamma = gamma
                    best_sequence = sequence_indices
                    best_cai = cai
                    
                    self.history_hashes.add(seq_hash)
                    self.used_gammas.add(gamma)
                    self.unique_sequences += 1
                    break
        
        elapsed_time = time.time() - start_time
        
        return {
            'time': elapsed_time,
            'gamma': best_gamma,
            'sequence': best_sequence,
            'cai': best_cai,
            'operations': operations,
            'unique_sequences': self.unique_sequences,
            'iteration': self.iteration_count,
            'constraint_satisfied': best_cai >= self.target_cai
        }


def test_correct_incremental():

    logger.info("="*80)

    logger.info("="*80)
    


    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, seq_length, p=weights))
    


    

    searcher = CorrectIncrementalPrecomputed(test_sequence, target_cai=0.8)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    num_iterations = 20
    results = []
    

    logger.info("-"*60)
    
    for i in range(num_iterations):

        if i == 0:

            current_dist = torch.rand(seq_length, 6, device=device)
            changed_positions = None
        else:

            current_dist = torch.rand(seq_length, 6, device=device)
            num_changes = seq_length // 10
            changed_positions = np.random.choice(seq_length, num_changes, replace=False).tolist()
        

        for pos in range(seq_length):
            if current_dist[pos].sum() > 0:
                current_dist[pos] = current_dist[pos] / current_dist[pos].sum()
        

        result = searcher.incremental_search_with_diversity(current_dist, changed_positions)
        results.append(result)
        


        logger.info(f"  CAI: {result['cai']:.4f} {'✅' if result['constraint_satisfied'] else '❌'}")
        logger.info(f"  Gamma: {result['gamma']:.4f}")




    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    first_time = results[0]['time']
    incremental_times = [r['time'] for r in results[1:]]
    avg_incremental = np.mean(incremental_times) if incremental_times else 0
    



    if avg_incremental > 0:

    

    final_unique = results[-1]['unique_sequences']



    

    satisfied_count = sum(1 for r in results if r['constraint_satisfied'])



    

    total_new = sum(r['operations']['new_switching_events'] for r in results)
    total_reused = sum(r['operations']['reused_switching_events'] for r in results)
    if total_new + total_reused > 0:
        reuse_rate = total_reused / (total_new + total_reused)


    
    return results


if __name__ == "__main__":

    results = test_correct_incremental()

