#!/usr/bin/env python3
"""


"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
import json
from typing import List, Tuple, Dict

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import SADOOptimizer

logger = setup_logging(level='INFO', name='feasible_region')


class FeasibleRegionAnalyzer:

    
    def __init__(self, sequence: str, target_cai: float = 0.8):
        self.sequence = sequence
        self.target_cai = target_cai
        self.seq_len = len(sequence)
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self._init_cai_weights()
        
    def _init_cai_weights(self):
        """初始化CAI权重表"""
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
    
    def compute_switching_events(self, target_dist: torch.Tensor) -> List[float]:
        """

        """
        all_gammas = []
        
        for pos in range(self.seq_len):
            if not self.codon_choices[pos]:
                continue
            
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
                    

                    denom = (p_j - p_i) + (w_i - w_j)
                    if abs(denom) > 1e-10:
                        gamma_switch = (p_j - p_i) / denom
                        if 0 <= gamma_switch <= 1:
                            all_gammas.append(gamma_switch)
        

        all_gammas.extend([0.0, 1.0])
        

        all_gammas = sorted(set(all_gammas))
        
        return all_gammas
    
    def generate_sequence_at_gamma(self, gamma: float, target_dist: torch.Tensor) -> np.ndarray:

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
    
    def calculate_cai_from_indices(self, indices: np.ndarray) -> float:
        """计算CAI值"""
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
    
    def analyze_feasible_region(self, target_dist: torch.Tensor) -> Dict:
        """


        """

        logger.info(f"  计算切换事件...")
        all_gammas = self.compute_switching_events(target_dist)
        logger.info(f"  总共{len(all_gammas)}个唯一gamma值")
        

        logger.info(f"  计算所有序列的CAI...")
        gamma_cai_pairs = []
        
        for i, gamma in enumerate(all_gammas):
            if i % 500 == 0:
                logger.info(f"    进度: {i}/{len(all_gammas)}")
            

            indices = self.generate_sequence_at_gamma(gamma, target_dist)
            

            cai = self.calculate_cai_from_indices(indices)
            
            gamma_cai_pairs.append((gamma, cai))
        

        gamma_cai_pairs.sort(key=lambda x: x[1])
        

        boundary_idx = -1
        for i, (gamma, cai) in enumerate(gamma_cai_pairs):
            if cai >= self.target_cai:
                boundary_idx = i
                break
        

        if boundary_idx >= 0:
            feasible_count = len(gamma_cai_pairs) - boundary_idx
            feasible_gammas = [g for g, c in gamma_cai_pairs[boundary_idx:]]
            feasible_cais = [c for g, c in gamma_cai_pairs[boundary_idx:]]
            

            best_idx = boundary_idx
            best_gamma = gamma_cai_pairs[best_idx][0]
            best_cai = gamma_cai_pairs[best_idx][1]
        else:
            feasible_count = 0
            feasible_gammas = []
            feasible_cais = []
            best_gamma = gamma_cai_pairs[-1][0]
            best_cai = gamma_cai_pairs[-1][1]
        

        results = {
            'sequence_length': self.seq_len,
            'total_gammas': len(all_gammas),
            'feasible_count': feasible_count,
            'feasible_ratio': feasible_count / len(all_gammas) if len(all_gammas) > 0 else 0,
            'best_gamma': best_gamma,
            'best_cai': best_cai,
            'min_cai': gamma_cai_pairs[0][1] if gamma_cai_pairs else 0,
            'max_cai': gamma_cai_pairs[-1][1] if gamma_cai_pairs else 0,
            'boundary_cai': gamma_cai_pairs[boundary_idx][1] if boundary_idx >= 0 else None,
            'feasible_gamma_range': (min(feasible_gammas), max(feasible_gammas)) if feasible_gammas else (None, None),
            'feasible_cai_range': (min(feasible_cais), max(feasible_cais)) if feasible_cais else (None, None)
        }
        
        return results


def test_real_sequences():

    logger.info("="*80)

    logger.info("="*80)
    

    data_dir = Path('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper/data')
    test_files = [





    ]
    
    all_results = {}
    
    for file_name in test_files:
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
        
        protein_name = file_name.replace('.fasta.txt', '')


        logger.info("-"*60)
        

        analyzer = FeasibleRegionAnalyzer(sequence, target_cai=0.8)
        

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        target_dist = torch.rand(len(sequence), 6, device=device)
        

        for pos in range(len(sequence)):
            if target_dist[pos].sum() > 0:
                target_dist[pos] = target_dist[pos] / target_dist[pos].sum()
        

        results = analyzer.analyze_feasible_region(target_dist)
        






        
        if results['feasible_count'] > 0:



        else:


        
        all_results[protein_name] = results
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    


    logger.info("-"*60)
    
    for protein_name, results in all_results.items():
        logger.info(f"{protein_name:<20} {results['sequence_length']:<8} "
                   f"{results['total_gammas']:<10} {results['feasible_count']:<10} "
                   f"{results['feasible_ratio']*100:<10.1f}%")
    

    if all_results:
        avg_ratio = np.mean([r['feasible_ratio'] for r in all_results.values()])
        avg_feasible = np.mean([r['feasible_count'] for r in all_results.values()])
        



        





    
    return all_results


if __name__ == "__main__":

    results = test_real_sequences()

