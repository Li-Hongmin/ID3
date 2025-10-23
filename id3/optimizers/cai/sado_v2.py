"""
SADO 2.0 - 基于帕累托前沿的多目标优化器

核心创新：
1. 多解探索：同时维护多个候选解
2. 帕累托前沿：找到CAI-概率-多样性的最优权衡
3. 自适应权重：根据历史动态调整优化策略
4. 智能扰动：基于贡献度的局部优化
"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from collections import deque

from ..base import BaseCAIOptimizer
from .utils import load_cai_weights, compute_cai
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


@dataclass
class Solution:
    """解的数据结构"""
    indices: np.ndarray
    cai: float
    prob_score: float
    diversity_score: float
    hash: str
    
    @property
    def fitness(self) -> float:
        """综合适应度"""
        return self.prob_score * 0.5 + self.cai * 0.3 + self.diversity_score * 0.2


class ParetoFrontier:
    """帕累托前沿管理器"""
    
    def __init__(self):
        self.solutions: List[Solution] = []
    
    def add(self, solution: Solution) -> bool:
        """添加解到前沿（如果不被支配）"""
        # 移除被新解支配的解
        self.solutions = [s for s in self.solutions if not self._dominates(solution, s)]
        
        # 检查新解是否被支配
        for s in self.solutions:
            if self._dominates(s, solution):
                return False
        
        self.solutions.append(solution)
        return True
    
    def _dominates(self, s1: Solution, s2: Solution) -> bool:
        """检查s1是否支配s2"""
        better_in_at_least_one = False
        
        # CAI比较
        if s1.cai < s2.cai:
            return False
        elif s1.cai > s2.cai:
            better_in_at_least_one = True
        
        # 概率比较
        if s1.prob_score < s2.prob_score:
            return False
        elif s1.prob_score > s2.prob_score:
            better_in_at_least_one = True
        
        # 多样性比较
        if s1.diversity_score < s2.diversity_score:
            return False
        elif s1.diversity_score > s2.diversity_score:
            better_in_at_least_one = True
        
        return better_in_at_least_one
    
    def get_best(self, weights: Tuple[float, float, float]) -> Solution:
        """根据权重获取最佳解"""
        if not self.solutions:
            return None
        
        best_solution = None
        best_score = -float('inf')
        
        for solution in self.solutions:
            score = (weights[0] * solution.prob_score + 
                    weights[1] * solution.cai + 
                    weights[2] * solution.diversity_score)
            if score > best_score:
                best_score = score
                best_solution = solution
        
        return best_solution


class AdaptiveWeightManager:
    """自适应权重管理器"""
    
    def __init__(self, target_cai: float = 0.8):
        self.target_cai = target_cai
        self.history = deque(maxlen=20)  # 最近20个解的历史
        
    def add_solution(self, solution: Solution):
        """添加解到历史"""
        self.history.append(solution)
    
    def get_weights(self) -> Tuple[float, float, float]:
        """计算当前权重"""
        if len(self.history) < 3:
            # 初始权重：平衡策略
            return (0.4, 0.4, 0.2)  # (prob, cai, diversity)
        
        # 分析历史
        avg_cai = np.mean([s.cai for s in self.history])
        avg_prob = np.mean([s.prob_score for s in self.history])
        unique_hashes = len(set(s.hash for s in self.history))
        diversity_rate = unique_hashes / len(self.history)
        
        # 动态调整权重
        w_cai = 0.3
        w_prob = 0.5
        w_div = 0.2
        
        # 如果CAI不达标，增加CAI权重
        if avg_cai < self.target_cai * 0.95:
            w_cai = 0.5
            w_prob = 0.3
        
        # 如果多样性不足，增加多样性权重
        if diversity_rate < 0.8:
            w_div = 0.3
            w_cai = 0.3
            w_prob = 0.4
        
        # 如果概率太低，增加概率权重
        if avg_prob < 0.15:  # 低于随机基准
            w_prob = 0.6
            w_cai = 0.25
            w_div = 0.15
        
        # 归一化
        total = w_prob + w_cai + w_div
        return (w_prob/total, w_cai/total, w_div/total)


class SADOv2Optimizer(BaseCAIOptimizer):
    """
    SADO 2.0 - 多目标优化版本
    """
    
    def __init__(self,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        super().__init__(species, device, amino_acid_sequence)
        
        # 加载CAI权重
        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)
        
        # 优化组件
        self.pareto_frontier = ParetoFrontier()
        self.weight_manager = AdaptiveWeightManager()
        self.history_hashes = set()
        self.solution_cache = {}  # 缓存已评估的解
        
        if amino_acid_sequence:
            self._build_codon_info(amino_acid_sequence)
    
    def _build_codon_info(self, amino_acid_sequence: str):
        """构建密码子信息"""
        self.codon_choices = []
        
        # 构建标准64密码子索引映射
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
                
                # 按权重排序
                choices.sort(key=lambda x: x['weight'], reverse=True)
                
                # 重新分配local_index
                for i, choice in enumerate(choices):
                    choice['local_index'] = i
                
                self.codon_choices.append(choices)
            else:
                self.codon_choices.append([])
    
    def _compute_cai_from_indices(self, indices: np.ndarray) -> float:
        """计算CAI值"""
        log_sum = 0.0
        valid_count = 0
        
        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                weight = self.codon_choices[pos][idx]['weight']
                if weight > 0:
                    log_sum += np.log(weight)
                    valid_count += 1
        
        return np.exp(log_sum / valid_count) if valid_count > 0 else 0.0
    
    def _compute_prob_score(self, indices: np.ndarray, pi_accessibility: torch.Tensor) -> float:
        """计算概率得分"""
        log_prob = 0.0
        
        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                orig_idx = self.codon_choices[pos][idx]['original_local_index']
                if orig_idx < pi_accessibility.shape[1]:
                    prob = pi_accessibility[pos, orig_idx].item()
                    if prob > 0:
                        log_prob += np.log(prob)
        
        return np.exp(log_prob / len(indices))
    
    def _compute_diversity_score(self, indices: np.ndarray) -> float:
        """计算多样性得分"""
        # 基于哈希唯一性
        seq_hash = hashlib.md5(indices.tobytes()).hexdigest()
        if seq_hash in self.history_hashes:
            return 0.0
        else:
            # 计算与历史解的最小汉明距离
            if not self.history_hashes:
                return 1.0
            
            # 简化：返回0.5表示新解
            return 0.5
    
    def _generate_initial_solutions(self, pi_accessibility: torch.Tensor, 
                                   target_cai: float) -> List[Solution]:
        """生成初始解集"""
        solutions = []
        
        # 1. 纯概率优化解
        prob_indices = self._probability_greedy(pi_accessibility)
        solutions.append(self._create_solution(prob_indices, pi_accessibility))
        
        # 2. 纯CAI优化解
        cai_indices = self._cai_greedy(target_cai)
        solutions.append(self._create_solution(cai_indices, pi_accessibility))
        
        # 3. 多个gamma平衡解
        for gamma in [0.1, 0.3, 0.5, 0.7]:
            balanced = self._gamma_initialization(pi_accessibility, gamma)
            solutions.append(self._create_solution(balanced, pi_accessibility))
        
        # 4. 随机扰动解
        for _ in range(5):
            base_idx = np.random.randint(0, len(solutions))
            perturbed = self._random_perturbation(solutions[base_idx].indices.copy())
            solutions.append(self._create_solution(perturbed, pi_accessibility))
        
        return solutions
    
    def _probability_greedy(self, pi_accessibility: torch.Tensor) -> np.ndarray:
        """纯概率贪心"""
        indices = np.zeros(len(self.codon_choices), dtype=np.int32)
        
        for pos in range(len(self.codon_choices)):
            if self.codon_choices[pos]:
                best_prob = -1
                best_idx = 0
                
                for choice in self.codon_choices[pos]:
                    orig_idx = choice['original_local_index']
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        if prob > best_prob:
                            best_prob = prob
                            best_idx = choice['local_index']
                
                indices[pos] = best_idx
        
        return indices
    
    def _cai_greedy(self, target_cai: float) -> np.ndarray:
        """纯CAI贪心"""
        indices = np.zeros(len(self.codon_choices), dtype=np.int32)
        
        # 直接选择最高权重
        for pos in range(len(self.codon_choices)):
            if self.codon_choices[pos]:
                indices[pos] = 0  # 已按权重排序，第一个最高
        
        return indices
    
    def _gamma_initialization(self, pi_accessibility: torch.Tensor, gamma: float) -> np.ndarray:
        """Gamma加权初始化"""
        indices = np.zeros(len(self.codon_choices), dtype=np.int32)
        
        for pos in range(len(self.codon_choices)):
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
            
            indices[pos] = best_idx
        
        return indices
    
    def _random_perturbation(self, indices: np.ndarray, rate: float = 0.1) -> np.ndarray:
        """随机扰动"""
        num_changes = max(1, int(len(indices) * rate))
        positions = np.random.choice(len(indices), num_changes, replace=False)
        
        for pos in positions:
            if self.codon_choices[pos] and len(self.codon_choices[pos]) > 1:
                current = indices[pos]
                choices = list(range(len(self.codon_choices[pos])))
                choices.remove(current)
                indices[pos] = np.random.choice(choices)
        
        return indices
    
    def _create_solution(self, indices: np.ndarray, pi_accessibility: torch.Tensor) -> Solution:
        """创建解对象"""
        seq_hash = hashlib.md5(indices.tobytes()).hexdigest()
        
        # 检查缓存
        if seq_hash in self.solution_cache:
            return self.solution_cache[seq_hash]
        
        solution = Solution(
            indices=indices.copy(),
            cai=self._compute_cai_from_indices(indices),
            prob_score=self._compute_prob_score(indices, pi_accessibility),
            diversity_score=self._compute_diversity_score(indices),
            hash=seq_hash
        )
        
        self.solution_cache[seq_hash] = solution
        return solution
    
    def _local_optimization(self, solution: Solution, pi_accessibility: torch.Tensor,
                           target_cai: float, mode: str = 'hybrid') -> Solution:
        """局部优化"""
        indices = solution.indices.copy()
        
        if mode == 'prob':
            # 概率引导的局部优化
            for _ in range(10):
                pos = np.random.randint(0, len(indices))
                if not self.codon_choices[pos]:
                    continue
                
                best_prob = -1
                best_idx = indices[pos]
                
                for choice in self.codon_choices[pos]:
                    test_indices = indices.copy()
                    test_indices[pos] = choice['local_index']
                    test_cai = self._compute_cai_from_indices(test_indices)
                    
                    if test_cai >= target_cai * 0.8:  # 放松约束
                        orig_idx = choice['original_local_index']
                        if orig_idx < pi_accessibility.shape[1]:
                            prob = pi_accessibility[pos, orig_idx].item()
                            if prob > best_prob:
                                best_prob = prob
                                best_idx = choice['local_index']
                
                indices[pos] = best_idx
        
        elif mode == 'cai':
            # CAI引导的局部优化
            for _ in range(10):
                pos = np.random.randint(0, len(indices))
                if not self.codon_choices[pos]:
                    continue
                
                current_cai = self._compute_cai_from_indices(indices)
                if current_cai >= target_cai:
                    break
                
                # 尝试提升CAI
                for choice in self.codon_choices[pos]:
                    if choice['local_index'] < indices[pos]:  # 更高权重
                        indices[pos] = choice['local_index']
                        break
        
        else:  # hybrid
            # 混合优化
            for _ in range(10):
                pos = np.random.randint(0, len(indices))
                if not self.codon_choices[pos]:
                    continue
                
                best_score = -1
                best_idx = indices[pos]
                
                for choice in self.codon_choices[pos]:
                    test_indices = indices.copy()
                    test_indices[pos] = choice['local_index']
                    
                    # 计算综合得分
                    test_cai = self._compute_cai_from_indices(test_indices)
                    orig_idx = choice['original_local_index']
                    
                    if orig_idx < pi_accessibility.shape[1]:
                        prob = pi_accessibility[pos, orig_idx].item()
                        
                        # 自适应权重
                        if test_cai < target_cai * 0.9:
                            score = test_cai * 0.7 + prob * 0.3
                        else:
                            score = test_cai * 0.3 + prob * 0.7
                        
                        if score > best_score:
                            best_score = score
                            best_idx = choice['local_index']
                
                indices[pos] = best_idx
        
        return self._create_solution(indices, pi_accessibility)
    
    def optimize(self, pi_accessibility: torch.Tensor, target_cai: float = 0.8,
                 **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """主优化函数"""
        logger.debug(f"SADO 2.0优化开始: 目标CAI={target_cai}")
        
        # 阶段1：生成初始解集
        solutions = self._generate_initial_solutions(pi_accessibility, target_cai)
        logger.debug(f"生成{len(solutions)}个初始解")
        
        # 阶段2：局部优化
        optimized_solutions = []
        for solution in solutions:
            # 尝试不同的优化策略
            opt_prob = self._local_optimization(solution, pi_accessibility, target_cai, 'prob')
            opt_cai = self._local_optimization(solution, pi_accessibility, target_cai, 'cai')
            opt_hybrid = self._local_optimization(solution, pi_accessibility, target_cai, 'hybrid')
            
            optimized_solutions.extend([opt_prob, opt_cai, opt_hybrid])
        
        logger.debug(f"局部优化后有{len(optimized_solutions)}个解")
        
        # 构建帕累托前沿
        for solution in optimized_solutions:
            self.pareto_frontier.add(solution)
        
        logger.debug(f"帕累托前沿包含{len(self.pareto_frontier.solutions)}个解")
        
        # 阶段3：智能选择
        weights = self.weight_manager.get_weights()
        best_solution = self.pareto_frontier.get_best(weights)
        
        if best_solution is None:
            # 降级到简单解
            best_solution = solutions[0]
        
        # 更新历史
        self.history_hashes.add(best_solution.hash)
        self.weight_manager.add_solution(best_solution)
        
        # 转换为分布
        distribution = self._indices_to_distribution(best_solution.indices, 
                                                    pi_accessibility.shape[1])
        
        # 构建元数据
        metadata = {
            'final_cai': best_solution.cai,
            'target_cai': target_cai,
            'constraint_satisfied': best_solution.cai >= target_cai,
            'prob_score': best_solution.prob_score,
            'diversity_score': best_solution.diversity_score,
            'pareto_size': len(self.pareto_frontier.solutions),
            'weights': weights,
            'method': 'sado_v2'
        }
        
        logger.debug(f"优化完成: CAI={best_solution.cai:.4f}, "
                    f"概率={best_solution.prob_score:.6f}, "
                    f"权重={weights}")
        
        return distribution, metadata
    
    def _indices_to_distribution(self, indices: np.ndarray, num_codons: int) -> torch.Tensor:
        """将索引转换为one-hot分布"""
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
        self.pareto_frontier = ParetoFrontier()
        self.weight_manager = AdaptiveWeightManager()
        self.history_hashes.clear()
        self.solution_cache.clear()
        logger.debug("SADO 2.0 optimizer reset")