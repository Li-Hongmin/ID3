#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from typing import List, Dict

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import SADOOptimizer

logger = setup_logging(level='INFO', name='switching_count')


class SwitchingEventCounter:

    
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.seq_len = len(sequence)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        sado = SADOOptimizer(
            species='ecoli_bl21de3',
            device=self.device,
            amino_acid_sequence=sequence
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
    
    def count_switching_events(self, target_dist: torch.Tensor) -> Dict:
        """统计切换事件（不计算CAI）"""
        all_gammas = []
        
        for pos in range(self.seq_len):
            if not self.codon_choices[pos]:
                continue
            
            choices = self.codon_choices[pos]
            

            for i in range(len(choices)):
                for j in range(i+1, len(choices)):
                    codon_i = choices[i]
                    codon_j = choices[j]
                    

                    p_i = target_dist[pos, codon_i['index']].item() if codon_i['index'] < target_dist.shape[1] else 0
                    p_j = target_dist[pos, codon_j['index']].item() if codon_j['index'] < target_dist.shape[1] else 0
                    

                    w_i = codon_i['weight']
                    w_j = codon_j['weight']
                    

                    denom = (p_j - p_i) + (w_i - w_j)
                    if abs(denom) > 1e-10:
                        gamma_switch = (p_j - p_i) / denom
                        if 0 <= gamma_switch <= 1:
                            all_gammas.append(gamma_switch)
        

        all_gammas.extend([0.0, 1.0])
        

        all_gammas = sorted(set(all_gammas))
        
        return {
            'total_events': len(all_gammas),
            'gammas': all_gammas
        }
    
    def estimate_feasible_region(self, all_gammas: List[float], target_dist: torch.Tensor, 
                                target_cai: float = 0.8, sample_size: int = 100) -> Dict:

        

        if len(all_gammas) <= sample_size:
            sample_gammas = all_gammas
        else:

            indices = np.linspace(0, len(all_gammas)-1, sample_size, dtype=int)
            sample_gammas = [all_gammas[i] for i in indices]
        

        feasible_count = 0
        min_cai = 1.0
        max_cai = 0.0
        boundary_gamma = None
        
        for gamma in sample_gammas:

            cai = self._quick_cai_estimate(gamma, target_dist, sample_positions=min(100, self.seq_len))
            
            if cai < min_cai:
                min_cai = cai
            if cai > max_cai:
                max_cai = cai
            
            if cai >= target_cai:
                feasible_count += 1
                if boundary_gamma is None or gamma < boundary_gamma:
                    boundary_gamma = gamma
        

        estimated_feasible = int(feasible_count * len(all_gammas) / len(sample_gammas))
        
        return {
            'sample_size': len(sample_gammas),
            'feasible_in_sample': feasible_count,
            'estimated_feasible_total': estimated_feasible,
            'estimated_ratio': estimated_feasible / len(all_gammas) if len(all_gammas) > 0 else 0,
            'min_cai': min_cai,
            'max_cai': max_cai,
            'boundary_gamma': boundary_gamma
        }
    
    def _quick_cai_estimate(self, gamma: float, target_dist: torch.Tensor, sample_positions: int = 100) -> float:
        """快速估算CAI（只用部分位置）"""
        log_sum = 0.0
        count = 0
        

        positions = np.random.choice(min(self.seq_len, sample_positions), 
                                   min(sample_positions, self.seq_len), 
                                   replace=False)
        
        for pos in positions:
            if not self.codon_choices[pos]:
                continue
            

            best_score = -float('inf')
            best_codon = None
            
            for choice in self.codon_choices[pos]:
                idx = choice['index']
                cai_score = choice['weight']
                prob_score = target_dist[pos, idx].item() if idx < target_dist.shape[1] else 0
                
                interpolated_score = gamma * cai_score + (1 - gamma) * prob_score
                
                if interpolated_score > best_score:
                    best_score = interpolated_score
                    best_codon = choice['codon']
            
            if best_codon:
                weight = self.wi_table.get(best_codon, 1e-10)
                log_sum += np.log(weight)
                count += 1
        
        if count == 0:
            return 0.0
        

        return np.exp(log_sum / count)


def test_real_sequences():

    logger.info("="*80)

    logger.info("="*80)
    

    data_dir = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/data')
    test_files = [





    ]
    
    results = []
    
    for file_name, description in test_files:
        file_path = data_dir / file_name
        if not file_path.exists():

            continue
        

        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        sequence = ''
        for line in lines:
            if not line.startswith('>'):
                sequence += line.strip()
        
        if not sequence:
            continue
        


        logger.info("-"*60)
        

        counter = SwitchingEventCounter(sequence)
        

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        target_dist = torch.rand(len(sequence), 6, device=device)
        

        for pos in range(len(sequence)):
            if target_dist[pos].sum() > 0:
                target_dist[pos] = target_dist[pos] / target_dist[pos].sum()
        


        event_stats = counter.count_switching_events(target_dist)

        



            sample_size = event_stats['total_events']


            sample_size = min(200, event_stats['total_events'])
        
        feasible_stats = counter.estimate_feasible_region(
            event_stats['gammas'], target_dist, target_cai=0.8, sample_size=sample_size
        )
        






        
        results.append({
            'name': description,
            'length': len(sequence),
            'total_events': event_stats['total_events'],
            'estimated_feasible': feasible_stats['estimated_feasible_total'],
            'feasible_ratio': feasible_stats['estimated_ratio'],
            'max_cai': feasible_stats['max_cai']
        })
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    


    logger.info("-"*80)
    
    for r in results:
        logger.info(f"{r['name']:<25} {r['length']:<8} {r['total_events']:<10} "
                   f"{r['estimated_feasible']:<10} {r['feasible_ratio']*100:<10.1f}% "
                   f"{r['max_cai']:<10.4f}")
    


    
    if results:

        short = [r for r in results if r['length'] < 100]
        medium = [r for r in results if 100 <= r['length'] < 500]
        long = [r for r in results if r['length'] >= 500]
        
        if short:
            avg_events = np.mean([r['total_events'] for r in short])
            avg_feasible = np.mean([r['estimated_feasible'] for r in short])

        
        if medium:
            avg_events = np.mean([r['total_events'] for r in medium])
            avg_feasible = np.mean([r['estimated_feasible'] for r in medium])

        
        if long:
            avg_events = np.mean([r['total_events'] for r in long])
            avg_feasible = np.mean([r['estimated_feasible'] for r in long])

        

        avg_ratio = np.mean([r['feasible_ratio'] for r in results])




    
    return results


if __name__ == "__main__":

    results = test_real_sequences()

