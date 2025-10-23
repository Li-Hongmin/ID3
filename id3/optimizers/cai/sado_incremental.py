"""







"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from collections import deque
from scipy.interpolate import interp1d

from ..base import BaseCAIOptimizer
from .utils import load_cai_weights, compute_cai
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


@dataclass
class OptimizationRecord:

    gamma: float
    cai: float
    prob_score: float
    sequence_hash: str
    indices: np.ndarray


class GammaRangePredictor:
    """Gamma区间预测器"""
    
    def __init__(self, target_cai: float = 0.8):
        self.target_cai = target_cai
        self.history: List[OptimizationRecord] = []
        self.gamma_cai_mapping = {}  # {gamma: avg_cai}
        
    def add_record(self, record: OptimizationRecord):

        self.history.append(record)
        

        if record.gamma not in self.gamma_cai_mapping:
            self.gamma_cai_mapping[record.gamma] = []
        self.gamma_cai_mapping[record.gamma].append(record.cai)
    
    def predict_range(self, iteration: int) -> np.ndarray:
        """预测最优gamma范围"""
        if len(self.history) < 3:

            if iteration == 0:
                return np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
            else:

                return np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        

        successful_gammas = []
        for record in self.history[-20:]:
            if abs(record.cai - self.target_cai) < 0.05:
                successful_gammas.append(record.gamma)
        
        if successful_gammas:

            center = np.median(successful_gammas)
            std = np.std(successful_gammas) if len(successful_gammas) > 1 else 0.1
            

            low = max(0, center - 1.5 * std)
            high = min(1, center + 1.5 * std)
            

            return np.linspace(low, high, 7)
        else:

            recent_cais = [r.cai for r in self.history[-5:]]
            avg_cai = np.mean(recent_cais)
            
            if avg_cai < self.target_cai * 0.85:

                last_gamma = self.history[-1].gamma
                return np.linspace(last_gamma, min(1.0, last_gamma + 0.4), 5)
            elif avg_cai < self.target_cai * 0.95:

                last_gamma = self.history[-1].gamma
                return np.linspace(last_gamma, min(1.0, last_gamma + 0.3), 5)
            elif avg_cai < self.target_cai:

                last_gamma = self.history[-1].gamma
                return np.linspace(last_gamma, min(1.0, last_gamma + 0.2), 5)
            else:

                last_gamma = self.history[-1].gamma

                return np.linspace(max(0, last_gamma - 0.1), min(1.0, last_gamma + 0.1), 5)
    
    def get_best_gamma(self) -> float:

        if not self.history:
            return 0.3
        

        best_record = min(self.history, 
                         key=lambda r: abs(r.cai - self.target_cai) - 0.1 * r.prob_score)
        return best_record.gamma


class SmartPerturbationEngine:
    """智能扰动引擎"""
    
    def __init__(self):
        self.perturbation_history = []
        
    def minimal_perturbation(self, 
                            indices: np.ndarray,
                            pi_accessibility: torch.Tensor,
                            codon_choices: List[List[Dict]],
                            num_changes: int = 2) -> np.ndarray:
        """

        """
        seq_len = len(indices)
        

        perturbation_candidates = []
        
        for pos in range(seq_len):
            if not codon_choices[pos]:
                continue
                
            current_idx = indices[pos]
            current_choice = codon_choices[pos][current_idx]
            current_orig_idx = current_choice['original_local_index']
            
            if current_orig_idx >= pi_accessibility.shape[1]:
                continue
                
            current_prob = pi_accessibility[pos, current_orig_idx].item()
            

            alternatives = []
            for choice in codon_choices[pos]:
                if choice['local_index'] == current_idx:
                    continue
                    
                orig_idx = choice['original_local_index']
                if orig_idx < pi_accessibility.shape[1]:
                    alt_prob = pi_accessibility[pos, orig_idx].item()
                    prob_ratio = alt_prob / (current_prob + 1e-8)
                    
                    if prob_ratio > 0.85:
                        alternatives.append({
                            'idx': choice['local_index'],
                            'prob_loss': current_prob - alt_prob,
                            'cai_weight': choice['weight']
                        })
            
            if alternatives:

                alternatives.sort(key=lambda x: x['cai_weight'], reverse=True)
                best_alt = alternatives[0]
                
                perturbation_candidates.append({
                    'pos': pos,
                    'new_idx': best_alt['idx'],
                    'prob_loss': best_alt['prob_loss'],
                    'score': best_alt['prob_loss'] / (best_alt['cai_weight'] + 0.1)
                })
        

        perturbation_candidates.sort(key=lambda x: x['score'])
        

        new_indices = indices.copy()
        actual_changes = min(num_changes, len(perturbation_candidates))
        
        for i in range(actual_changes):
            candidate = perturbation_candidates[i]
            new_indices[candidate['pos']] = candidate['new_idx']
            logger.debug(f"扰动位置{candidate['pos']}: "
                        f"概率损失={candidate['prob_loss']:.4f}")
        
        return new_indices
    
    def diversity_perturbation(self,
                             indices: np.ndarray,
                             history_hashes: set,
                             codon_choices: List[List[Dict]]) -> np.ndarray:
        """

        """
        seq_hash = hashlib.md5(indices.tobytes()).hexdigest()
        
        if seq_hash not in history_hashes:
            return indices
        

        seq_len = len(indices)
        changeable_positions = []
        
        for pos in range(seq_len):
            if codon_choices[pos] and len(codon_choices[pos]) > 1:
                changeable_positions.append(pos)
        
        if not changeable_positions:
            return indices
        

        num_changes = min(3, len(changeable_positions))
        positions = np.random.choice(changeable_positions, num_changes, replace=False)
        
        new_indices = indices.copy()
        for pos in positions:
            current_idx = indices[pos]
            choices = list(range(len(codon_choices[pos])))
            choices.remove(current_idx)
            
            if choices:

                choices.sort(key=lambda idx: codon_choices[pos][idx]['weight'], reverse=True)
                new_indices[pos] = choices[0]
        
        return new_indices
    
    def cai_compensation(self,
                        indices: np.ndarray,
                        target_cai: float,
                        codon_choices: List[List[Dict]],
                        perturbed_positions: set,
                        max_adjustments: int = 5) -> np.ndarray:
        """

        

        """

        current_cai = self._compute_cai_simple(indices, codon_choices)
        
        if current_cai >= target_cai * 0.95:
            return indices
        

        improvement_candidates = []
        
        for pos in range(len(indices)):
            if pos in perturbed_positions:
                continue
                
            if not codon_choices[pos]:
                continue
                
            current_idx = indices[pos]
            current_weight = codon_choices[pos][current_idx]['weight']
            

            for choice in codon_choices[pos]:
                if choice['local_index'] < current_idx:
                    potential_gain = choice['weight'] - current_weight
                    if potential_gain > 0:
                        improvement_candidates.append({
                            'pos': pos,
                            'new_idx': choice['local_index'],
                            'cai_gain': potential_gain,
                            'new_weight': choice['weight']
                        })
                        break
        
        if not improvement_candidates:
            return indices
        

        compensated_indices = indices.copy()
        adjustments_made = 0
        

        gains = np.array([c['cai_gain'] for c in improvement_candidates])

        temperature = 0.2
        exp_gains = np.exp(gains / temperature)
        probabilities = exp_gains / exp_gains.sum()
        

        num_to_adjust = min(max_adjustments, len(improvement_candidates))

        cai_gap = target_cai - current_cai
        if cai_gap < 0.05:
            num_to_adjust = min(3, num_to_adjust)
        elif cai_gap < 0.1:
            num_to_adjust = min(5, num_to_adjust)
        elif cai_gap < 0.15:
            num_to_adjust = min(7, num_to_adjust)
        

        selected_indices = np.random.choice(
            len(improvement_candidates),
            size=min(num_to_adjust, len(improvement_candidates)),
            replace=False,
            p=probabilities
        )
        

        for idx in selected_indices:
            candidate = improvement_candidates[idx]
            

            compensated_indices[candidate['pos']] = candidate['new_idx']
            adjustments_made += 1
            

            new_cai = self._compute_cai_simple(compensated_indices, codon_choices)
            
            logger.debug(f"CAI补偿 位置{candidate['pos']}: "
                        f"CAI {current_cai:.4f} → {new_cai:.4f}")
            
            current_cai = new_cai
            

            if new_cai >= target_cai * 0.95:
                break
        
        if adjustments_made > 0:
            logger.debug(f"CAI补偿完成: {adjustments_made}个位置调整, "
                        f"最终CAI={current_cai:.4f}")
        
        return compensated_indices
    
    def _compute_cai_simple(self, indices: np.ndarray, codon_choices: List[List[Dict]]) -> float:

        log_sum = 0.0
        valid_count = 0
        
        for pos, idx in enumerate(indices):
            if pos < len(codon_choices) and codon_choices[pos]:
                weight = codon_choices[pos][idx]['weight']
                if weight > 0:
                    log_sum += np.log(weight)
                    valid_count += 1
        
        return np.exp(log_sum / valid_count) if valid_count > 0 else 0.0


class SADOIncrementalOptimizer(BaseCAIOptimizer):
    """
    SADO 3.0 - 增量优化版本
    """
    
    def __init__(self,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        super().__init__(species, device, amino_acid_sequence)
        

        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)
        


        self.perturbation_engine = SmartPerturbationEngine()
        self.history_hashes = set()
        self.optimization_history = []
        

        self.cai_cache = {}  # {seq_hash: cai}
        self.prob_cache = {}  # {seq_hash: prob_score}
        
        if amino_acid_sequence:
            self._build_codon_info(amino_acid_sequence)
    
    def _build_codon_info(self, amino_acid_sequence: str):
        """构建密码子信息"""
        self.codon_choices = []
        

        bases = ['U', 'C', 'A', 'G']
        codon_to_global_index = {}
        for i, base1 in enumerate(bases):
            for j, base2 in enumerate(bases):
                for k, base3 in enumerate(bases):
                    codon = base1 + base2 + base3
                    global_index = i * 16 + j * 4 + k
                    codon_to_global_index[codon] = global_index
        
        for pos, aa in enumerate(amino_acid_sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                choices = []
                
                for i, codon in enumerate(codons):
                    weight = self.wi_table.get(codon, 0.0)
                    global_index = codon_to_global_index.get(codon, 0)
                    choices.append({
                        'original_local_index': i,
                        'global_index': global_index,
                        'codon': codon,
                        'weight': weight,
                        'aa': aa
                    })
                

                choices.sort(key=lambda x: x['weight'], reverse=True)
                

                for i, choice in enumerate(choices):
                    choice['local_index'] = i
                
                self.codon_choices.append(choices)
            else:
                self.codon_choices.append([])
    
    def _batch_compute_cai(self, indices_batch: np.ndarray) -> np.ndarray:

        batch_size = indices_batch.shape[0]
        cai_values = np.zeros(batch_size)
        
        for b in range(batch_size):
            indices = indices_batch[b]
            seq_hash = hashlib.md5(indices.tobytes()).hexdigest()
            

            if seq_hash in self.cai_cache:
                cai_values[b] = self.cai_cache[seq_hash]
            else:

                log_sum = 0.0
                valid_count = 0
                
                for pos, idx in enumerate(indices):
                    if pos < len(self.codon_choices) and self.codon_choices[pos]:
                        weight = self.codon_choices[pos][idx]['weight']
                        if weight > 0:
                            log_sum += np.log(weight)
                            valid_count += 1
                
                cai = np.exp(log_sum / valid_count) if valid_count > 0 else 0.0
                cai_values[b] = cai
                self.cai_cache[seq_hash] = cai
        
        return cai_values
    
    def _batch_evaluate_gammas(self, 
                              gamma_range: np.ndarray,
                              pi_accessibility: torch.Tensor) -> List[Dict]:
        """批量评估多个gamma值"""
        batch_size = len(gamma_range)
        seq_len = len(self.codon_choices)
        

        indices_batch = np.zeros((batch_size, seq_len), dtype=np.int32)
        
        for b, gamma in enumerate(gamma_range):

            for pos in range(seq_len):
                if not self.codon_choices[pos]:
                    continue
                
                best_score = -1
                best_idx = 0
                
                for choice in self.codon_choices[pos]:
                    orig_idx = choice['original_local_index']
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        weight = choice['weight']
                        
                        if prob > 0 and weight > 0:

                            score = (prob ** (1 - gamma)) * (weight ** gamma)
                            if score > best_score:
                                best_score = score
                                best_idx = choice['local_index']
                
                indices_batch[b, pos] = best_idx
        

        cai_values = self._batch_compute_cai(indices_batch)
        

        prob_scores = np.zeros(batch_size)
        for b in range(batch_size):
            indices = indices_batch[b]
            log_prob = 0.0
            
            for pos, idx in enumerate(indices):
                if pos < len(self.codon_choices) and self.codon_choices[pos]:
                    orig_idx = self.codon_choices[pos][idx]['original_local_index']
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        if prob > 0:
                            log_prob += np.log(prob)
            
            prob_scores[b] = np.exp(log_prob / seq_len)
        

        results = []
        for b in range(batch_size):
            results.append({
                'gamma': gamma_range[b],
                'indices': indices_batch[b],
                'cai': cai_values[b],
                'prob_score': prob_scores[b],
                'hash': hashlib.md5(indices_batch[b].tobytes()).hexdigest()
            })
        
        return results
    
    def _interpolate_optimal_gamma(self, 
                                  gamma_range: np.ndarray,
                                  cai_values: np.ndarray,
                                  target_cai: float) -> float:

        if len(gamma_range) < 3:

            idx = np.argmin(np.abs(cai_values - target_cai))
            return gamma_range[idx]
        
        try:

            f = interp1d(gamma_range, cai_values, kind='cubic', 
                        bounds_error=False, fill_value='extrapolate')
            

            gamma_fine = np.linspace(gamma_range.min(), gamma_range.max(), 100)
            cai_fine = f(gamma_fine)
            

            idx = np.argmin(np.abs(cai_fine - target_cai))
            return gamma_fine[idx]
        except:

            idx = np.argmin(np.abs(cai_values - target_cai))
            return gamma_range[idx]
    
    def optimize(self, 
                pi_accessibility: torch.Tensor, 
                target_cai: float = 0.8,
                **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """主优化函数"""
        

        if self.gamma_predictor is None:
            self.gamma_predictor = GammaRangePredictor(target_cai)
        
        iteration = len(self.optimization_history)
        logger.debug(f"增量优化第{iteration+1}次迭代，目标CAI={target_cai}")
        

        gamma_range = self.gamma_predictor.predict_range(iteration)
        logger.debug(f"Gamma搜索范围: {gamma_range}")
        

        candidates = self._batch_evaluate_gammas(gamma_range, pi_accessibility)
        


        best_candidate = None
        best_score = -float('inf')
        
        for candidate in candidates:


            if candidate['cai'] >= target_cai:

                cai_score = 1.0 + (candidate['cai'] - target_cai) * 0.5
            else:

                cai_score = candidate['cai'] / target_cai
            
            prob_score = candidate['prob_score']
            

            if candidate['cai'] < target_cai * 0.9:

                combined_score = cai_score * 0.8 + prob_score * 0.2
            elif candidate['cai'] < target_cai:

                combined_score = cai_score * 0.6 + prob_score * 0.4
            else:

                combined_score = cai_score * 0.3 + prob_score * 0.7
            
            if combined_score > best_score:
                best_score = combined_score
                best_candidate = candidate
        

        original_hash = best_candidate['hash']
        attempts = 0
        max_attempts = 10
        perturbed_positions = set()
        
        while best_candidate['hash'] in self.history_hashes and attempts < max_attempts:
            logger.debug(f"检测到重复（尝试{attempts+1}），执行扰动")
            

            num_changes = min(2 + attempts // 2, 5)
            

            original_indices = best_candidate['indices'].copy()
            

            if attempts % 2 == 0:

                perturbed_indices = self.perturbation_engine.minimal_perturbation(
                    best_candidate['indices'],
                    pi_accessibility,
                    self.codon_choices,
                    num_changes=num_changes
                )
            else:

                perturbed_indices = self.perturbation_engine.diversity_perturbation(
                    best_candidate['indices'],
                    self.history_hashes,
                    self.codon_choices
                )
            

            for pos in range(len(original_indices)):
                if original_indices[pos] != perturbed_indices[pos]:
                    perturbed_positions.add(pos)
            

            cai_values = self._batch_compute_cai(perturbed_indices.reshape(1, -1))
            new_hash = hashlib.md5(perturbed_indices.tobytes()).hexdigest()
            

            best_candidate['indices'] = perturbed_indices
            best_candidate['cai'] = cai_values[0]
            best_candidate['hash'] = new_hash
            

            log_prob = 0.0
            for pos, idx in enumerate(perturbed_indices):
                if pos < len(self.codon_choices) and self.codon_choices[pos]:
                    orig_idx = self.codon_choices[pos][idx]['original_local_index']
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        if prob > 0:
                            log_prob += np.log(prob)
            best_candidate['prob_score'] = np.exp(log_prob / len(perturbed_indices))
            
            attempts += 1
        
        if attempts > 0:
            logger.debug(f"扰动成功，CAI: {best_candidate['cai']:.4f}, "
                        f"概率: {best_candidate['prob_score']:.6f}")
        


        cai_before_compensation = best_candidate['cai']
        compensation_positions = []
        
        if best_candidate['cai'] < target_cai * 0.95:
            logger.debug(f"CAI不足({best_candidate['cai']:.4f})，执行补偿")
            

            if len(perturbed_positions) == 0:


                num_protected = min(20, len(best_candidate['indices']) // 10)
                protected_positions = set(np.random.choice(
                    len(best_candidate['indices']), 
                    size=num_protected, 
                    replace=False
                ))
            else:
                protected_positions = perturbed_positions
            

            compensated_indices = self.perturbation_engine.cai_compensation(
                best_candidate['indices'],
                target_cai,
                self.codon_choices,
                protected_positions,
                max_adjustments=15
            )
            

            for pos in range(len(best_candidate['indices'])):
                if best_candidate['indices'][pos] != compensated_indices[pos]:
                    compensation_positions.append(pos)
            

            cai_values = self._batch_compute_cai(compensated_indices.reshape(1, -1))
            
            if cai_values[0] > best_candidate['cai']:

                best_candidate['indices'] = compensated_indices
                best_candidate['cai'] = cai_values[0]
                best_candidate['hash'] = hashlib.md5(compensated_indices.tobytes()).hexdigest()
                

                log_prob = 0.0
                for pos, idx in enumerate(compensated_indices):
                    if pos < len(self.codon_choices) and self.codon_choices[pos]:
                        orig_idx = self.codon_choices[pos][idx]['original_local_index']
                        if orig_idx < pi_accessibility.shape[1]:
                            prob = pi_accessibility[pos, orig_idx].item()
                            if prob > 0:
                                log_prob += np.log(prob)
                best_candidate['prob_score'] = np.exp(log_prob / len(compensated_indices))
                
                logger.debug(f"CAI补偿成功: {cai_before_compensation:.4f} → {cai_values[0]:.4f}")
        

        record = OptimizationRecord(
            gamma=best_candidate['gamma'],
            cai=best_candidate['cai'],
            prob_score=best_candidate['prob_score'],
            sequence_hash=best_candidate['hash'],
            indices=best_candidate['indices']
        )
        
        self.optimization_history.append(record)
        self.gamma_predictor.add_record(record)
        self.history_hashes.add(best_candidate['hash'])
        

        distribution = self._indices_to_distribution(
            best_candidate['indices'],
            pi_accessibility.shape[1]
        )
        

        metadata = {
            'final_cai': best_candidate['cai'],
            'target_cai': target_cai,
            'constraint_satisfied': best_candidate['cai'] >= target_cai,
            'prob_score': best_candidate['prob_score'],
            'gamma': best_candidate['gamma'],
            'iteration': iteration,
            'unique_sequences': len(self.history_hashes),
            'gamma_range': gamma_range.tolist(),
            'method': 'sado_incremental',
            'cai_before_compensation': cai_before_compensation,
            'compensation_positions': compensation_positions
        }
        
        logger.debug(f"优化完成: CAI={best_candidate['cai']:.4f}, "
                    f"概率={best_candidate['prob_score']:.6f}, "
                    f"Gamma={best_candidate['gamma']:.3f}")
        
        return distribution, metadata
    
    def _indices_to_distribution(self, indices: np.ndarray, num_codons: int) -> torch.Tensor:

        seq_len = len(indices)
        distribution = torch.zeros(seq_len, num_codons, device=self.device)
        
        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                orig_idx = self.codon_choices[pos][idx]['original_local_index']
                if orig_idx < num_codons:
                    distribution[pos, orig_idx] = 1.0
        
        return distribution
    
    def reset(self):
        """重置优化器状态"""
        self.gamma_predictor = None
        self.perturbation_engine = SmartPerturbationEngine()
        self.history_hashes.clear()
        self.optimization_history.clear()
        self.cai_cache.clear()
        self.prob_cache.clear()
        logger.debug("SADO增量优化器已重置")
    
    def get_statistics(self) -> Dict[str, Any]:

        if not self.optimization_history:
            return {}
        
        cai_values = [r.cai for r in self.optimization_history]
        prob_scores = [r.prob_score for r in self.optimization_history]
        gamma_values = [r.gamma for r in self.optimization_history]
        
        return {
            'num_iterations': len(self.optimization_history),
            'unique_sequences': len(self.history_hashes),
            'avg_cai': np.mean(cai_values),
            'std_cai': np.std(cai_values),
            'avg_prob': np.mean(prob_scores),
            'std_prob': np.std(prob_scores),
            'avg_gamma': np.mean(gamma_values),
            'convergence_gamma': gamma_values[-1] if gamma_values else None,
            'cache_size': len(self.cai_cache)
        }