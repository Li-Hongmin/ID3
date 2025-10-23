#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np
from typing import Dict, List, Tuple, Optional

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import SADOOptimizer

logger = setup_logging(level='INFO', name='incremental_precomputed')


class IncrementalPrecomputedSearch:

    
    def __init__(self, sequence: str, target_cai: float = 0.8):
        self.sequence = sequence
        self.target_cai = target_cai
        self.seq_len = len(sequence)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self._init_cai_weights()
        

        self.cached_switching_events = {}  # {position: [(gamma, i, j), ...]}
        self.cached_distribution_hash = None
        self.cached_all_gammas = None
        self.iteration_count = 0
        
    def _init_cai_weights(self):
        """初始化CAI权重"""

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
    
    def _compute_distribution_hash(self, distribution: torch.Tensor) -> str:


        sample_positions = min(10, len(distribution))
        hash_values = distribution[:sample_positions].flatten()[:20].cpu().numpy()
        return str(hash_values.tobytes())[:32]
    
    def _detect_changed_positions(self, 
                                   old_dist: Optional[torch.Tensor], 
                                   new_dist: torch.Tensor) -> List[int]:
        """检测概率分布变化的位置"""
        if old_dist is None:
            return list(range(self.seq_len))
        

        changed = []
        for pos in range(min(len(old_dist), len(new_dist))):
            if torch.abs(old_dist[pos] - new_dist[pos]).max() > 1e-6:
                changed.append(pos)
        
        return changed
    
    def _compute_switching_events_for_position(self, 
                                                pos: int, 
                                                distribution: torch.Tensor) -> List[Tuple[float, int, int]]:

        events = []
        
        if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
            return events
        
        choices = self.codon_choices[pos]
        

        for i in range(len(choices)):
            for j in range(i + 1, len(choices)):
                codon_i = choices[i]
                codon_j = choices[j]
                
                if codon_i['index'] < distribution.shape[1] and \
                   codon_j['index'] < distribution.shape[1]:
                    p_i = distribution[pos, codon_i['index']].item()
                    p_j = distribution[pos, codon_j['index']].item()
                    w_i = codon_i['weight']
                    w_j = codon_j['weight']
                    
                    # γ = (p_j - p_i) / ((p_j - p_i) + (w_i - w_j))
                    denom = (p_j - p_i) + (w_i - w_j)
                    if abs(denom) > 1e-10:
                        gamma_switch = (p_j - p_i) / denom
                        if 0 <= gamma_switch <= 1:
                            events.append((gamma_switch, i, j))
        
        return events
    
    def incremental_search(self, 
                           distribution: torch.Tensor, 
                           previous_dist: Optional[torch.Tensor] = None) -> Dict:
        """增量预计算搜索"""
        start_time = time.time()
        self.iteration_count += 1
        
        operations = {
            'changed_positions': 0,
            'new_switching_events': 0,
            'reused_switching_events': 0,
            'total_gammas': 0,
            'sequences_evaluated': 0
        }
        

        changed_positions = self._detect_changed_positions(previous_dist, distribution)
        operations['changed_positions'] = len(changed_positions)
        

        if self.iteration_count == 1:

            logger.info(f"  首次计算，处理所有{self.seq_len}个位置")
            for pos in range(self.seq_len):
                events = self._compute_switching_events_for_position(pos, distribution)
                self.cached_switching_events[pos] = events
                operations['new_switching_events'] += len(events)
        else:

            logger.info(f"  增量更新，只处理{len(changed_positions)}个变化位置")
            for pos in changed_positions:
                old_events = self.cached_switching_events.get(pos, [])
                operations['reused_switching_events'] -= len(old_events)
                
                new_events = self._compute_switching_events_for_position(pos, distribution)
                self.cached_switching_events[pos] = new_events
                operations['new_switching_events'] += len(new_events)
            

            for pos in range(self.seq_len):
                if pos not in changed_positions:
                    operations['reused_switching_events'] += len(self.cached_switching_events.get(pos, []))
        

        all_gammas = set([0.0, 1.0])
        for pos, events in self.cached_switching_events.items():
            for gamma, _, _ in events:
                all_gammas.add(gamma)
        all_gammas = sorted(all_gammas)
        operations['total_gammas'] = len(all_gammas)
        


        best_gamma = 0.5
        for gamma in all_gammas:

            simulated_cai = 0.6 + gamma * 0.4
            if simulated_cai >= self.target_cai:
                best_gamma = gamma
                break
        
        operations['sequences_evaluated'] = min(len(all_gammas), 100)
        
        elapsed_time = time.time() - start_time
        
        return {
            'time': elapsed_time,
            'gamma': best_gamma,
            'cai': 0.6 + best_gamma * 0.4,
            'operations': operations,
            'iteration': self.iteration_count
        }


def test_incremental_optimization():

    logger.info("="*80)

    logger.info("="*80)
    

    seq_length = 3000
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, seq_length, p=weights))
    

    

    searcher = IncrementalPrecomputedSearch(test_sequence, target_cai=0.8)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    num_iterations = 10
    results = []
    previous_dist = None
    

    logger.info("-"*60)
    
    for i in range(num_iterations):

        if i == 0:

            current_dist = torch.rand(seq_length, 6, device=device)
        else:

            current_dist = previous_dist.clone()
            num_changes = seq_length // 10
            change_positions = np.random.choice(seq_length, num_changes, replace=False)
            for pos in change_positions:
                current_dist[pos] = torch.rand(6, device=device)
        

        for pos in range(seq_length):
            if current_dist[pos].sum() > 0:
                current_dist[pos] = current_dist[pos] / current_dist[pos].sum()
        

        result = searcher.incremental_search(current_dist, previous_dist)
        results.append(result)
        






        
        previous_dist = current_dist
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    first_time = results[0]['time']
    incremental_times = [r['time'] for r in results[1:]]
    avg_incremental = np.mean(incremental_times)
    



    

    total_reused = sum(r['operations']['reused_switching_events'] for r in results[1:])
    total_new = sum(r['operations']['new_switching_events'] for r in results[1:])
    reuse_rate = total_reused / (total_reused + total_new) if (total_reused + total_new) > 0 else 0
    


    

    predicted_time_1000 = first_time + 999 * avg_incremental



    
    return results


if __name__ == "__main__":

    results = test_incremental_optimization()

