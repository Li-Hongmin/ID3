"""
SADO (Sequence Adaptive Discretization Optimizer)


"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, Tuple, Optional, Any, Set, List

from ..base import BaseCAIOptimizer
from .utils import load_cai_weights, compute_cai
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


class SADOOptimizer(BaseCAIOptimizer):
    """
    SADO - Sequence Adaptive Discretization Optimizer










    """

    def __init__(self,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        """


        Args:



        """
        super().__init__(species, device, amino_acid_sequence)


        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)


        self.last_indices = None
        self.history_hashes: Set[str] = set()
        self.all_sequences: List[np.ndarray] = []
        

        self.discrete_searcher = None
        self.last_binary_search_result = None


        if amino_acid_sequence:
            self._build_codon_info(amino_acid_sequence)

    def _build_codon_info(self, amino_acid_sequence: str):

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

    def _hash_sequence(self, indices: np.ndarray) -> str:
        """计算序列的哈希值"""
        return hashlib.md5(indices.tobytes()).hexdigest()

    def _indices_to_distribution(self, local_indices: np.ndarray, num_codons: int = 6) -> torch.Tensor:

        seq_len = len(local_indices)
        distribution = torch.zeros(seq_len, num_codons, device=self.device)

        for pos, local_idx in enumerate(local_indices):
            if pos < len(self.codon_choices) and local_idx < len(self.codon_choices[pos]):

                original_idx = self.codon_choices[pos][local_idx]['original_local_index']
                if original_idx < num_codons:
                    distribution[pos, original_idx] = 1.0

        return distribution

    def _compute_cai_from_indices(self, local_indices: np.ndarray) -> float:
        """从本地密码子索引计算CAI值（使用对数提高数值稳定性）"""
        if len(local_indices) == 0:
            return 0.0

        log_weights_sum = 0.0
        valid_positions = 0

        for pos, local_idx in enumerate(local_indices):
            if (pos < len(self.codon_choices) and
                local_idx < len(self.codon_choices[pos])):
                weight = self.codon_choices[pos][local_idx]['weight']
                if weight > 0:
                    log_weights_sum += np.log(weight)
                    valid_positions += 1


        if valid_positions > 0:
            return np.exp(log_weights_sum / valid_positions)
        return 0.0

    def _equilibrium_state_initialization(self, target_cai: float = 0.8) -> np.ndarray:
        """

        




        
        Args:

            
        Returns:

        """
        seq_len = len(self.codon_choices)
        indices = np.zeros(seq_len, dtype=np.int32)
        

        position_info = []
        total_log_sum = 0.0
        valid_count = 0
        
        for pos in range(seq_len):
            if self.codon_choices[pos]:
                choices = self.codon_choices[pos]

                sorted_choices = sorted(choices, key=lambda x: x['weight'], reverse=True)
                

                if target_cai > 0.75:

                    base_idx = len(sorted_choices) // 3
                else:

                    base_idx = len(sorted_choices) // 2
                

                offset_range = max(1, int(len(sorted_choices) * 0.4))
                offset = np.random.randint(-offset_range, offset_range + 1)
                selected_idx = max(0, min(len(sorted_choices) - 1, base_idx + offset))
                indices[pos] = sorted_choices[selected_idx]['local_index']
                mid_idx = selected_idx
                

                weight = sorted_choices[mid_idx]['weight']
                if weight > 0:
                    total_log_sum += np.log(weight)
                    valid_count += 1
                
                position_info.append({
                    'pos': pos,
                    'choices': sorted_choices,
                    'current_level': mid_idx
                })
        

        current_cai = np.exp(total_log_sum / valid_count) if valid_count > 0 else 0.0
        

        max_iterations = 10
        
        for _ in range(max_iterations):
            if abs(current_cai - target_cai) < 0.01:
                break
                
            if current_cai < target_cai:

                positions_to_upgrade = []
                for info in position_info:
                    if info['current_level'] > 0:
                        positions_to_upgrade.append(info)
                

                upgrade_ratio = 0.2 + np.random.random() * 0.2
                upgrade_count = max(1, int(len(positions_to_upgrade) * upgrade_ratio))
                np.random.shuffle(positions_to_upgrade)
                for info in positions_to_upgrade[:upgrade_count]:
                    pos = info['pos']
                    old_level = info['current_level']
                    new_level = old_level - 1
                    

                    old_weight = info['choices'][old_level]['weight']
                    new_weight = info['choices'][new_level]['weight']
                    
                    indices[pos] = info['choices'][new_level]['local_index']
                    info['current_level'] = new_level
                    

                    if old_weight > 0:
                        total_log_sum -= np.log(old_weight)
                    if new_weight > 0:
                        total_log_sum += np.log(new_weight)
            else:

                positions_to_downgrade = []
                for info in position_info:
                    if info['current_level'] < len(info['choices']) - 1:
                        positions_to_downgrade.append(info)
                

                downgrade_ratio = 0.2 + np.random.random() * 0.2
                downgrade_count = max(1, int(len(positions_to_downgrade) * downgrade_ratio))
                np.random.shuffle(positions_to_downgrade)
                for info in positions_to_downgrade[:downgrade_count]:
                    pos = info['pos']
                    old_level = info['current_level']
                    new_level = old_level + 1
                    

                    old_weight = info['choices'][old_level]['weight']
                    new_weight = info['choices'][new_level]['weight']
                    
                    indices[pos] = info['choices'][new_level]['local_index']
                    info['current_level'] = new_level
                    

                    if old_weight > 0:
                        total_log_sum -= np.log(old_weight)
                    if new_weight > 0:
                        total_log_sum += np.log(new_weight)
            

            current_cai = np.exp(total_log_sum / valid_count) if valid_count > 0 else 0.0
        
        logger.debug(f"快速平衡状态初始化: CAI={current_cai:.4f}, 目标={target_cai:.4f}")
        return indices
    
    def _gamma_initialization(self, pi_accessibility: torch.Tensor, gamma: float = 0.3) -> np.ndarray:
        """

        

        
        Args:






        
        Returns:

        """
        seq_len = pi_accessibility.shape[0]
        indices = np.zeros(seq_len, dtype=np.int32)

        for pos in range(seq_len):
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            choices = self.codon_choices[pos]
            probs = pi_accessibility[pos].detach().cpu().numpy()
            
            best_score = -1
            best_idx = 0

            for choice in choices:
                orig_idx = choice['original_local_index']
                prob = probs[orig_idx] if orig_idx < len(probs) else 1e-8
                cai_weight = choice['weight']
                

                if prob > 0 and cai_weight > 0:
                    score = (prob ** (1 - gamma)) * (cai_weight ** gamma)
                    
                    if score > best_score:
                        best_score = score
                        best_idx = choice['local_index']

            indices[pos] = best_idx

        return indices
    
    def _initial_optimization(self, pi_accessibility: torch.Tensor) -> np.ndarray:
        """


        """
        return self._gamma_initialization(pi_accessibility, gamma=0.0)

    def _calculate_efficiency_ratio(self, pos: int, current_idx: int, new_idx: int,
                                   delta_cai: float, choices: List[Dict],
                                   pi_accessibility: Optional[torch.Tensor]) -> float:

        if pi_accessibility is None:


        current_choice = choices[current_idx]
        new_choice = choices[new_idx]

        current_orig_idx = current_choice['original_local_index']
        new_orig_idx = new_choice['original_local_index']


        if (current_orig_idx >= pi_accessibility.shape[1] or
            new_orig_idx >= pi_accessibility.shape[1]):
            return delta_cai

        current_prob = pi_accessibility[pos, current_orig_idx].item()
        new_prob = pi_accessibility[pos, new_orig_idx].item()


        if current_prob > 0 and new_prob > 0:
            delta_log_prob = np.log(new_prob) - np.log(current_prob)


                return delta_cai / abs(delta_log_prob)

                return float('inf')



    def _compute_probability_differences(self, current_indices: np.ndarray, 
                                        pi_accessibility: torch.Tensor,
                                        equilibrium_indices: np.ndarray) -> np.ndarray:
        """
        计算当前序列与平衡状态的概率差异
        
        Args:
            current_indices: 当前密码子索引
            pi_accessibility: 概率分布
            equilibrium_indices: 平衡状态索引
            
        Returns:
            每个位置的概率差异值
        """
        seq_len = len(current_indices)
        differences = np.zeros(seq_len)
        
        for pos in range(seq_len):
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue
                
            current_idx = current_indices[pos]
            equilibrium_idx = equilibrium_indices[pos]
            

            current_orig = self.codon_choices[pos][current_idx]['original_local_index']
            equilibrium_orig = self.codon_choices[pos][equilibrium_idx]['original_local_index']
            

            if current_orig < pi_accessibility.shape[1] and equilibrium_orig < pi_accessibility.shape[1]:
                current_prob = pi_accessibility[pos, current_orig].item()
                equilibrium_prob = pi_accessibility[pos, equilibrium_orig].item()
                differences[pos] = abs(current_prob - equilibrium_prob)
        
        return differences
    
    def _difference_driven_greedy_optimization(self, 
                                              initial_indices: np.ndarray,
                                              pi_accessibility: torch.Tensor,
                                              target_cai: float,
                                              max_iterations: int = 20) -> np.ndarray:
        """
        基于差异的贪心优化（增强版：包含概率采样和历史避让）
        
        核心思想：
        1. 使用平衡状态作为参考
        2. 基于概率分布进行采样（而非确定性选择）
        3. 历史避让机制确保序列多样性
        4. 平衡CAI约束和概率最大化
        
        Args:
            initial_indices: 初始密码子索引（可以是二分查找结果或平衡状态）
            pi_accessibility: 概率分布
            target_cai: 目标CAI值
            max_iterations: 最大迭代次数
            
        Returns:
            优化后的密码子索引
        """
        current_indices = initial_indices.copy()
        current_cai = self._compute_cai_from_indices(current_indices)
        

        equilibrium_indices = self._equilibrium_state_initialization(target_cai)
        

        

        temperature = 1.0
        temperature_decay = 0.95
        
        for iteration in range(max_iterations):

            differences = self._compute_probability_differences(
                current_indices, pi_accessibility, equilibrium_indices
            )
            

            diff_probs = differences / (differences.sum() + 1e-8)
            

            sampling_probs = np.exp(diff_probs / temperature)
            sampling_probs = sampling_probs / sampling_probs.sum()
            

            num_positions = min(15, len(differences))
            selected_positions = np.random.choice(
                len(differences), 
                size=num_positions, 
                replace=False,
                p=sampling_probs
            )
            

            improved = False
            for pos in selected_positions:
                if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                    continue
                    
                current_idx = current_indices[pos]
                choices = self.codon_choices[pos]
                

                candidates = []
                
                for choice in choices:
                    new_idx = choice['local_index']
                    if new_idx == current_idx:
                        continue
                        

                    test_indices = current_indices.copy()
                    test_indices[pos] = new_idx
                    test_cai = self._compute_cai_from_indices(test_indices)
                    

                    test_hash = self._hash_sequence(test_indices)
                    is_duplicate = test_hash in self.history_hashes
                    



                        orig_idx = choice['original_local_index']
                        if orig_idx < pi_accessibility.shape[1]:
                            prob_score = pi_accessibility[pos, orig_idx].item()
                            


                            log_prob_score = np.log(prob_score + 1e-8)

                            norm_cai = (test_cai - 0.7) / 0.2
                            norm_cai = max(0, min(1, norm_cai))
                            

                            combined_score = np.exp(log_prob_score * 0.85 + np.log(norm_cai + 0.1) * 0.15)
                            
                            if is_duplicate:

                            
                            candidates.append({
                                'idx': new_idx,
                                'score': combined_score,
                                'cai': test_cai,
                                'prob': prob_score,
                                'is_duplicate': is_duplicate
                            })
                

                if candidates:

                    candidates.sort(key=lambda x: x['score'], reverse=True)
                    

                    scores = np.array([c['score'] for c in candidates])
                    sample_probs = np.exp(scores / temperature)
                    sample_probs = sample_probs / sample_probs.sum()
                    

                    selected_idx = np.random.choice(len(candidates), p=sample_probs)
                    selected = candidates[selected_idx]
                    

                    current_indices[pos] = selected['idx']
                    old_cai = current_cai
                    current_cai = selected['cai']
                    improved = True
                    
                    if iteration % 5 == 0:


                    break
            

            temperature *= temperature_decay

            
            if not improved or current_cai >= target_cai:
                break
        

        return current_indices
    
    def _greedy_cai_optimization(self, indices: np.ndarray, target_cai: float, max_iterations: int = 100, pi_accessibility: Optional[torch.Tensor] = None) -> np.ndarray:
        """
        效率比值优化 - SADO的核心机制

        使用效率比值 E = ΔCA I/ |Δlog P| 选择最优替换，真正实现概率最大化
        """
        current_cai = self._compute_cai_from_indices(indices)
        iterations = 0



        while current_cai < target_cai and iterations < max_iterations:
            efficiencies = []


            for pos in range(len(indices)):
                if pos >= len(self.codon_choices):
                    continue

                current_idx = indices[pos]
                choices = self.codon_choices[pos]


                for new_idx in range(len(choices)):
                    if new_idx == current_idx:
                        continue


                    test_indices = indices.copy()
                    test_indices[pos] = new_idx
                    test_cai = self._compute_cai_from_indices(test_indices)
                    delta_cai = test_cai - current_cai



                        efficiency = self._calculate_efficiency_ratio(
                            pos, current_idx, new_idx, delta_cai, choices, pi_accessibility
                        )

                        efficiencies.append({
                            'pos': pos,
                            'new_idx': new_idx,
                            'new_cai': test_cai,
                            'cai_gain': delta_cai,
                            'efficiency': efficiency
                        })

            if not efficiencies:

                break


            efficiencies.sort(key=lambda x: x['efficiency'], reverse=True)



            selected = None

            for candidate in efficiencies:

                test_indices = indices.copy()
                test_indices[candidate['pos']] = candidate['new_idx']
                test_hash = self._hash_sequence(test_indices)


                if test_hash not in self.history_hashes:
                    selected = candidate

                    break
                else:



            if selected is None:
                selected = efficiencies[0]



            if selected['new_cai'] >= target_cai:

            else:



            indices[selected['pos']] = selected['new_idx']
            old_cai = current_cai
            current_cai = selected['new_cai']
            iterations += 1

            if iterations % 10 == 0:



        return indices

    def _initialize_binary_searcher(self):
        """初始化二分查找器"""
        if self.discrete_searcher is None:
            try:
                from id3.cai.discrete_binary_search import DiscreteCAISearcher

                cai_weights = []
                from id3.utils.constants import codons
                for codon in codons:
                    dna_codon = codon.replace('U', 'T')
                    weight = self.wi_table.get(dna_codon, 0.1)
                    cai_weights.append(weight)
                self.discrete_searcher = DiscreteCAISearcher(
                    cai_weights, 
                    str(self.device),
                    use_incremental=True
                )
            except ImportError:
                logger.warning("无法导入DiscreteCAISearcher，将使用默认初始化")
                self.discrete_searcher = None
    
    def _get_sequence_from_binary_search(self, pi_accessibility: torch.Tensor, 
                                        target_cai: float) -> Optional[np.ndarray]:
        """

        
        Args:


            
        Returns:

        """
        if self.discrete_searcher is None:
            self._initialize_binary_searcher()
            
        if self.discrete_searcher is None:
            return None
            

        seq_len = pi_accessibility.shape[0]
        num_codons = pi_accessibility.shape[1]
        
        valid_codon_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=self.device)
        codon_indices = torch.zeros(seq_len, num_codons, dtype=torch.long, device=self.device)
        
        for pos in range(seq_len):
            if pos < len(self.codon_choices):
                for choice in self.codon_choices[pos]:
                    orig_idx = choice['original_local_index']
                    global_idx = choice['global_index']
                    if orig_idx < num_codons:
                        valid_codon_mask[pos, orig_idx] = True
                        codon_indices[pos, orig_idx] = global_idx
        

        try:
            optimal_alpha, metadata = self.discrete_searcher.discrete_binary_search(
                pi_probs=pi_accessibility,
                valid_codon_mask=valid_codon_mask,
                codon_indices=codon_indices,
                amino_sequence=self.amino_acid_sequence,
                target_cai=target_cai,
                verbose=False
            )
            

            indices = np.zeros(seq_len, dtype=np.int32)
            for pos in range(seq_len):
                if pos < len(self.codon_choices):
                    choices = self.codon_choices[pos]
                    if optimal_alpha > 0.5:

                        indices[pos] = 0
                    else:

                        probs = pi_accessibility[pos].detach().cpu().numpy()
                        best_prob = -1
                        best_idx = 0
                        for choice in choices:
                            orig_idx = choice['original_local_index']
                            if orig_idx < len(probs) and probs[orig_idx] > best_prob:
                                best_prob = probs[orig_idx]
                                best_idx = choice['local_index']
                        indices[pos] = best_idx
            
            self.last_binary_search_result = {
                'alpha': optimal_alpha,
                'metadata': metadata,
                'indices': indices.copy()
            }
            
            logger.debug(f"二分查找成功: α={optimal_alpha:.4f}, CAI={self._compute_cai_from_indices(indices):.4f}")
            return indices
            
        except Exception as e:
            logger.warning(f"二分查找失败: {str(e)}")
            return None
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float = 0.8,
                 amino_acid_sequence: Optional[str] = None,
                 gamma: float = 0.3,
                 use_binary_search: bool = True,
                 use_difference_driven: bool = True) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """







        Args:







        Returns:

        """

        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)


        if amino_acid_sequence and amino_acid_sequence != self.amino_acid_sequence:
            self.amino_acid_sequence = amino_acid_sequence
            self._build_codon_info(amino_acid_sequence)


        indices = None
        init_method = "unknown"
        

        if use_binary_search and self.last_indices is None:
            indices = self._get_sequence_from_binary_search(pi_accessibility, target_cai)
            if indices is not None:
                init_method = "binary_search"
                logger.debug(f"使用二分查找初始化，CAI={self._compute_cai_from_indices(indices):.4f}")
        

        if indices is None and self.last_indices is not None:
            indices = self.last_indices.copy()
            init_method = "history"
            logger.debug(f"使用历史序列，CAI={self._compute_cai_from_indices(indices):.4f}")
        

        if indices is None:
            if use_difference_driven:

                indices = self._equilibrium_state_initialization(target_cai)
                init_method = "equilibrium"
                logger.debug(f"使用平衡状态初始化，CAI={self._compute_cai_from_indices(indices):.4f}")
            else:

                indices = self._gamma_initialization(pi_accessibility, gamma=gamma)
                init_method = "gamma"
                logger.debug(f"使用Gamma={gamma}初始化，CAI={self._compute_cai_from_indices(indices):.4f}")


        current_cai = self._compute_cai_from_indices(indices)
        
        if use_difference_driven:

            indices = self._difference_driven_greedy_optimization(
                indices, pi_accessibility, target_cai, max_iterations=50
            )
        elif current_cai < target_cai:

            indices = self._greedy_cai_optimization(
                indices, target_cai, max_iterations=50, pi_accessibility=pi_accessibility
            )


        seq_hash = self._hash_sequence(indices)
        collision_detected = seq_hash in self.history_hashes

        if collision_detected:
            logger.warning("检测到序列重复 - 历史避让机制未完全生效")


        self.history_hashes.add(seq_hash)
        self.last_indices = indices
        self.all_sequences.append(indices.copy())


        distribution = self._indices_to_distribution(indices, pi_accessibility.shape[1])


        final_cai = self._compute_cai_from_indices(indices)


        metadata = {
            'final_cai': final_cai,
            'target_cai': target_cai,
            'constraint_satisfied': final_cai >= target_cai,
            'init_method': init_method,
            'use_binary_search': use_binary_search,
            'use_difference_driven': use_difference_driven,
            'unique_sequences': len(self.history_hashes),
            'collision_detected': collision_detected,
            'method': 'sado_enhanced'
        }
        

        if self.last_binary_search_result:
            metadata['binary_search_alpha'] = self.last_binary_search_result.get('alpha')

        return distribution, metadata

    def reset(self):

        self.last_indices = None
        self.history_hashes.clear()
        self.all_sequences.clear()
        self.reset_statistics()
        logger.debug("SADO optimizer state reset")

    def get_diversity_stats(self) -> Dict[str, Any]:
        """获取序列多样性统计"""
        if not self.all_sequences:
            return {'num_sequences': 0, 'unique_ratio': 0.0}

        unique_hashes = set()
        for seq in self.all_sequences:
            unique_hashes.add(self._hash_sequence(seq))

        return {
            'num_sequences': len(self.all_sequences),
            'num_unique': len(unique_hashes),
            'unique_ratio': len(unique_hashes) / len(self.all_sequences),
            'repetition_rate': 1.0 - (len(unique_hashes) / len(self.all_sequences))
        }
