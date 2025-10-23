"""
增量式CAI优化器 (Incremental CAI Optimizer)

结合Binary Search的稳定性和SADO的增量优化效率：
- 第一次迭代：使用Binary Search获得高质量初始解
- 后续迭代：使用增量优化，基于上次结果进行局部改进

这种方法在长序列和多次迭代时显著提升性能。
"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, Tuple, Optional, Any, Set, List

from ..base import BaseCAIOptimizer
from .binary_search import BinarySearchCAIOptimizer
from .utils import load_cai_weights, compute_cai
from id3.utils.constants import amino_acids_to_codons

logger = logging.getLogger(__name__)


class IncrementalCAIOptimizer(BaseCAIOptimizer):
    """
    增量式CAI优化器

    核心创新：
    1. 首次使用Binary Search建立高质量基线
    2. 后续迭代使用增量优化，只调整必要的位置
    3. 保持历史最优解，避免性能退化
    4. 智能判断何时需要重新初始化
    """

    # 优化参数常量
    MIN_K_POSITIONS = 5                    # 最少优化位置数
    MAX_K_POSITIONS = 20                   # 最多优化位置数
    K_POSITION_DIVISOR = 50                # 位置数除数
    MIN_K_PROBABILITY = 3                  # 概率优化最少位置数
    MAX_K_PROBABILITY = 10                 # 概率优化最多位置数
    K_PROBABILITY_DIVISOR = 100            # 概率位置数除数

    # 重新初始化阈值
    CONSECUTIVE_FAILURE_THRESHOLD = 5      # 连续失败次数阈值
    RANDOM_REINIT_PROBABILITY = 0.05       # 随机重新初始化概率

    # 扰动参数
    MIN_PERTURB_RATIO = 0.05               # 最小扰动比例
    MAX_PERTURB_RATIO = 0.10               # 最大扰动比例
    MAX_PERTURB_POSITIONS = 20             # 最大扰动位置数

    # 数值稳定性
    NUMERICAL_EPSILON = np.finfo(float).eps  # 机器精度
    CAI_EPSILON = 1e-8                       # CAI计算epsilon
    PROBABILITY_EPSILON = 1e-10              # 概率计算epsilon

    # 评分权重
    CAI_SCORE_WEIGHT = 0.3                  # CAI得分权重
    PROBABILITY_SCORE_WEIGHT = 0.7          # 概率得分权重
    
    def __init__(self,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None,
                 random_seed: Optional[int] = None,
                 enable_caching: bool = True):
        """
        初始化增量优化器

        Args:
            species: 物种名称
            device: 计算设备
            amino_acid_sequence: 氨基酸序列
            random_seed: 随机种子，用于结果可重现
            enable_caching: 是否启用缓存优化
        """
        super().__init__(species, device, amino_acid_sequence)

        # 设置随机数生成器
        self.random_seed = random_seed
        self.rng = np.random.RandomState(random_seed) if random_seed is not None else np.random

        # 初始化Binary Search优化器
        self.bs_optimizer = BinarySearchCAIOptimizer(species, device, amino_acid_sequence)
        
        # 加载CAI权重
        self.wi_table, self.weights_tensor = load_cai_weights(species)
        self.weights_tensor = self.weights_tensor.to(self.device)
        
        # 状态管理
        self.last_result = None  # 上一次的优化结果
        self.last_indices = None  # 上一次的密码子索引
        self.last_gamma = None  # 上一次的最优gamma
        self.iteration_count = 0  # 迭代计数
        self.consecutive_failures = 0  # 连续失败次数
        
        # 历史记录
        self.history_hashes: Set[str] = set()  # 避免重复
        self.performance_history: List[Dict] = []  # 性能历史

        # 缓存优化
        self.enable_caching = enable_caching
        if enable_caching:
            self._cai_cache: Dict[str, float] = {}  # CAI计算缓存
            self._prob_cache: Dict[Tuple, float] = {}  # 概率计算缓存
            self._cache_hits = 0
            self._cache_misses = 0

        # 构建密码子信息
        if amino_acid_sequence:
            self._build_codon_info(amino_acid_sequence)
    
    def _build_codon_info(self, amino_acid_sequence: str):
        """构建密码子信息用于增量优化"""
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
                        'local_index': i,
                        'global_index': global_index,
                        'codon': codon,
                        'weight': weight,
                        'aa': aa
                    })
                
                # 按权重排序
                choices.sort(key=lambda x: x['weight'], reverse=True)
                self.codon_choices.append(choices)
            else:
                self.codon_choices.append([])
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float,
                 amino_acid_sequence: str,
                 valid_codon_mask: Optional[torch.Tensor] = None,
                 **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """
        执行增量优化

        第一次调用使用Binary Search，后续使用增量优化
        """
        # 输入验证
        try:
            self._validate_inputs(pi_accessibility, amino_acid_sequence, target_cai)
        except ValueError as e:
            logger.error(f"Input validation failed: {e}")
            raise

        self.iteration_count += 1
        
        # 检查是否需要重新初始化
        should_reinit = self._should_reinitialize(pi_accessibility)
        
        if self.last_result is None or should_reinit:
            # 首次或需要重新初始化：使用Binary Search
            logger.debug(f"Iteration {self.iteration_count}: Using Binary Search initialization")
            
            discrete_dist, metadata = self.bs_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )
            
            # 保存结果
            self.last_result = discrete_dist
            self.last_indices = discrete_dist.argmax(dim=-1).cpu().numpy()
            self.last_gamma = metadata.get('optimal_gamma', 0.5)
            self.consecutive_failures = 0
            
            # 记录到历史
            seq_hash = self._hash_sequence(self.last_indices)
            self.history_hashes.add(seq_hash)
            
            metadata['method'] = 'incremental_bs'
            metadata['iteration'] = self.iteration_count
            
            return discrete_dist, metadata
        
        else:
            # 后续迭代：使用增量优化
            logger.debug(f"Iteration {self.iteration_count}: Using incremental optimization")
            
            # 执行增量优化
            optimized_indices = self._incremental_optimize(
                self.last_indices.copy(),
                pi_accessibility,
                target_cai,
                amino_acid_sequence,
                valid_codon_mask
            )
            
            # 检查是否产生了新序列
            seq_hash = self._hash_sequence(optimized_indices)
            is_duplicate = seq_hash in self.history_hashes
            
            if is_duplicate:
                # 如果重复，尝试扰动
                optimized_indices = self._perturb_sequence(
                    optimized_indices,
                    pi_accessibility,
                    target_cai,
                    valid_codon_mask
                )
                seq_hash = self._hash_sequence(optimized_indices)
                is_duplicate = seq_hash in self.history_hashes
            
            # 获取正确的codon数量（处理batch维度）
            if pi_accessibility.dim() == 3:
                num_codons = pi_accessibility.shape[2]
            else:
                num_codons = pi_accessibility.shape[1]
            
            # 转换为离散分布
            discrete_dist = self._indices_to_distribution(
                optimized_indices,
                num_codons
            )
            
            # 计算CAI
            final_cai = self._compute_cai_from_indices(optimized_indices)
            
            # 更新状态
            if final_cai >= target_cai:
                self.last_result = discrete_dist
                self.last_indices = optimized_indices
                self.consecutive_failures = 0
            else:
                self.consecutive_failures += 1
            
            # 记录到历史
            if not is_duplicate:
                self.history_hashes.add(seq_hash)
            
            # 构建元数据
            metadata = {
                'method': 'incremental_inc',
                'iteration': self.iteration_count,
                'final_cai': final_cai,
                'target_cai': target_cai,
                'is_duplicate': is_duplicate,
                'incremental_changes': self._count_changes(self.last_indices, optimized_indices),
                'constraint_satisfied': final_cai >= target_cai
            }
            
            return discrete_dist, metadata
    
    def _incremental_optimize(self,
                              indices: np.ndarray,
                              pi_accessibility: torch.Tensor,
                              target_cai: float,
                              amino_acid_sequence: str,
                              valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        执行增量优化：只优化需要改进的位置

        Args:
            indices: 当前密码子索引
            pi_accessibility: RNA可及性概率
            target_cai: 目标CAI值
            amino_acid_sequence: 氨基酸序列
            valid_codon_mask: 有效密码子掩码

        Returns:
            优化后的密码子索引
        """
        current_cai = self._compute_cai_from_indices(indices)

        # 如果已经满足CAI要求，尝试优化概率
        if current_cai >= target_cai:
            return self._optimize_probability(
                indices, pi_accessibility, target_cai, valid_codon_mask
            )

        # 选择需要优化的位置
        positions_to_optimize = self._select_positions_for_optimization(
            indices, pi_accessibility, target_cai
        )

        # 优化选中的位置
        return self._optimize_selected_positions(
            indices, positions_to_optimize, pi_accessibility, target_cai
        )

    def _select_positions_for_optimization(self,
                                          indices: np.ndarray,
                                          pi_accessibility: torch.Tensor,
                                          target_cai: float) -> np.ndarray:
        """
        选择需要优化的位置

        Returns:
            需要优化的位置索引数组
        """
        # 计算每个位置的改进潜力
        improvements = self._calculate_improvement_potential(
            indices, pi_accessibility, target_cai
        )

        # 选择top-k个位置进行优化
        k = self._calculate_k_positions(len(indices))
        return np.argsort(improvements)[-k:]

    def _optimize_selected_positions(self,
                                    indices: np.ndarray,
                                    positions: np.ndarray,
                                    pi_accessibility: torch.Tensor,
                                    target_cai: float) -> np.ndarray:
        """
        优化选中的位置

        Args:
            indices: 当前密码子索引
            positions: 要优化的位置
            pi_accessibility: RNA可及性概率
            target_cai: 目标CAI值

        Returns:
            优化后的密码子索引
        """
        for pos in positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue

            # 找到该位置的最佳密码子
            best_idx = self._find_best_codon_for_position(
                indices, pos, pi_accessibility, target_cai
            )
            indices[pos] = best_idx

        return indices

    def _find_best_codon_for_position(self,
                                     indices: np.ndarray,
                                     pos: int,
                                     pi_accessibility: torch.Tensor,
                                     target_cai: float) -> int:
        """
        为特定位置找到最佳密码子

        Returns:
            最佳密码子的局部索引
        """
        best_idx = indices[pos]
        best_score = self._evaluate_position(
            indices, pos, best_idx, pi_accessibility, target_cai
        )

        # 尝试所有可能的密码子
        for choice in self.codon_choices[pos]:
            new_idx = choice['local_index']
            if new_idx == best_idx:
                continue

            score = self._evaluate_position(
                indices, pos, new_idx, pi_accessibility, target_cai
            )

            if score > best_score:
                best_idx = new_idx
                best_score = score

        return best_idx
    
    def _optimize_probability(self,
                              indices: np.ndarray,
                              pi_accessibility: torch.Tensor,
                              target_cai: float,
                              valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        在满足CAI约束的情况下优化概率
        """
        # 处理batch维度
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)
        
        # 计算当前概率得分
        current_prob = self._compute_probability_score(indices, pi_accessibility)
        
        # 尝试改进概率最低的位置
        prob_scores = []
        for pos in range(len(indices)):
            if pos < pi_accessibility.shape[0] and indices[pos] < pi_accessibility.shape[1]:
                prob_scores.append(pi_accessibility[pos, indices[pos]].item())
            else:
                prob_scores.append(1.0)
        
        # 选择概率最低的位置
        k = self._calculate_k_probability_positions(len(indices))
        worst_positions = np.argsort(prob_scores)[:k]
        
        for pos in worst_positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue
            
            current_idx = indices[pos]
            best_idx = current_idx
            best_prob = prob_scores[pos]
            
            # 尝试更高概率的密码子
            for choice in self.codon_choices[pos]:
                new_idx = choice['local_index']
                if new_idx == current_idx:
                    continue
                
                # 检查是否仍满足CAI约束
                test_indices = indices.copy()
                test_indices[pos] = new_idx
                test_cai = self._compute_cai_from_indices(test_indices)
                
                if test_cai >= target_cai:  # 严格约束，无容差
                    if new_idx < pi_accessibility.shape[1]:
                        new_prob = pi_accessibility[pos, new_idx].item()
                        if new_prob > best_prob:
                            best_idx = new_idx
                            best_prob = new_prob
            
            indices[pos] = best_idx
        
        return indices
    
    def _calculate_improvement_potential(self,
                                        indices: np.ndarray,
                                        pi_accessibility: torch.Tensor,
                                        target_cai: float) -> np.ndarray:
        """
        计算每个位置的改进潜力
        """
        potentials = np.zeros(len(indices))
        
        for pos in range(len(indices)):
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue
            
            current_idx = indices[pos]
            current_weight = 0.0
            max_weight = 0.0
            
            # 找到当前和最大权重
            for choice in self.codon_choices[pos]:
                if choice['local_index'] == current_idx:
                    current_weight = choice['weight']
                max_weight = max(max_weight, choice['weight'])
            
            # 计算CAI改进潜力
            cai_potential = (max_weight - current_weight) / (max_weight + self.CAI_EPSILON)
            
            # 计算概率得分
            prob_score = 0.0
            if pos < pi_accessibility.shape[0] and current_idx < pi_accessibility.shape[1]:
                prob_score = pi_accessibility[pos, current_idx].item()
            
            # 综合潜力：CAI改进潜力 * (1 - 当前概率)
            potentials[pos] = cai_potential * (1 - prob_score)
        
        return potentials
    
    def _evaluate_position(self,
                          indices: np.ndarray,
                          pos: int,
                          new_idx: int,
                          pi_accessibility: torch.Tensor,
                          target_cai: float) -> float:
        """
        评估某个位置使用特定密码子的得分
        """
        # 处理batch维度
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)
        
        # 创建测试序列
        test_indices = indices.copy()
        test_indices[pos] = new_idx
        
        # 计算CAI
        test_cai = self._compute_cai_from_indices(test_indices)
        
        # 如果不满足CAI约束，返回负分
        if test_cai < target_cai:
            return -1.0
        
        # 计算概率得分
        prob_score = 0.0
        if pos < pi_accessibility.shape[0] and new_idx < pi_accessibility.shape[1]:
            prob_score = pi_accessibility[pos, new_idx].item()
        
        # 综合得分：CAI达成度 + 概率得分
        cai_score = min(1.0, test_cai / target_cai)
        return cai_score * self.CAI_SCORE_WEIGHT + prob_score * self.PROBABILITY_SCORE_WEIGHT
    
    def _perturb_sequence(self,
                         indices: np.ndarray,
                         pi_accessibility: torch.Tensor,
                         target_cai: float,
                         valid_codon_mask: Optional[torch.Tensor]) -> np.ndarray:
        """
        扰动序列以生成新的唯一序列
        """
        # 处理batch维度
        if pi_accessibility.dim() == 3:
            pi_accessibility = pi_accessibility.squeeze(0)
        
        perturbed = indices.copy()
        num_positions = len(indices)
        
        # 随机选择位置进行扰动
        num_perturb = self._calculate_perturb_positions(num_positions)
        positions = self.rng.choice(num_positions, num_perturb, replace=False)
        
        for pos in positions:
            if pos >= len(self.codon_choices) or not self.codon_choices[pos]:
                continue
            
            choices = self.codon_choices[pos]
            if len(choices) <= 1:
                continue
            
            # 选择一个不同的密码子
            current_idx = perturbed[pos]
            candidates = [c['local_index'] for c in choices if c['local_index'] != current_idx]
            
            if candidates:
                # 基于概率选择
                probs = []
                for idx in candidates:
                    if idx < pi_accessibility.shape[1]:
                        probs.append(pi_accessibility[pos, idx].item())
                    else:
                        probs.append(0.1)
                
                probs = np.array(probs)
                
                # 改进的数值稳定性处理
                probs = self._normalize_probabilities(probs)

                new_idx = self.rng.choice(candidates, p=probs)
                perturbed[pos] = new_idx
        
        # 确保仍满足CAI约束
        perturbed_cai = self._compute_cai_from_indices(perturbed)
        if perturbed_cai < target_cai:
            # 如果不满足，返回原序列
            return indices
        
        return perturbed
    
    def _should_reinitialize(self, pi_accessibility: torch.Tensor) -> bool:
        """
        判断是否需要重新初始化
        """
        # 连续失败太多次
        if self.consecutive_failures > self.CONSECUTIVE_FAILURE_THRESHOLD:
            logger.debug("Reinitializing due to consecutive failures")
            return True
        
        # 输入分布变化太大
        if self.last_result is not None:
            # 简单的变化检测（可以更复杂）
            if self.rng.random() < self.RANDOM_REINIT_PROBABILITY:  # 概率重新初始化
                logger.debug("Random reinitialization")
                return True
        
        return False
    
    def _compute_cai_from_indices(self, indices: np.ndarray) -> float:
        """计算CAI值，带缓存优化"""
        # 生成缓存键
        if self.enable_caching:
            cache_key = hashlib.md5(indices.tobytes()).hexdigest()
            if cache_key in self._cai_cache:
                self._cache_hits += 1
                return self._cai_cache[cache_key]
            self._cache_misses += 1

        cai_product = 1.0
        valid_count = 0

        for pos, idx in enumerate(indices):
            if pos < len(self.codon_choices) and self.codon_choices[pos]:
                choices = self.codon_choices[pos]
                for choice in choices:
                    if choice['local_index'] == idx:
                        weight = choice['weight']
                        if weight > 0:
                            cai_product *= weight
                            valid_count += 1
                        break

        result = cai_product ** (1.0 / valid_count) if valid_count > 0 else 0.0

        # 存入缓存
        if self.enable_caching:
            # 限制缓存大小
            if len(self._cai_cache) > 10000:
                # 清理一半缓存
                keys_to_remove = list(self._cai_cache.keys())[:5000]
                for key in keys_to_remove:
                    del self._cai_cache[key]
            self._cai_cache[cache_key] = result

        return result
    
    def _compute_probability_score(self,
                                  indices: np.ndarray,
                                  pi_accessibility: torch.Tensor) -> float:
        """计算概率得分"""
        # 处理batch维度
        if pi_accessibility.dim() == 3:
            # 如果有batch维度，squeeze掉（假设batch_size=1）
            pi_accessibility = pi_accessibility.squeeze(0)
        
        log_prob = 0.0
        count = 0
        
        for pos, idx in enumerate(indices):
            if pos < pi_accessibility.shape[0] and idx < pi_accessibility.shape[1]:
                prob = pi_accessibility[pos, idx].item()
                if prob > 0:
                    log_prob += np.log(prob)
                    count += 1
        
        if count > 0:
            return np.exp(log_prob / count)
        return 0.0
    
    def _count_changes(self, old_indices: np.ndarray, new_indices: np.ndarray) -> int:
        """统计变化的位置数"""
        if old_indices is None or new_indices is None:
            return 0
        return np.sum(old_indices != new_indices)
    
    def _hash_sequence(self, indices: np.ndarray) -> str:
        """计算序列哈希"""
        return hashlib.md5(indices.tobytes()).hexdigest()
    
    def _indices_to_distribution(self, indices: np.ndarray, num_codons: int) -> torch.Tensor:
        """将索引转换为one-hot分布"""
        seq_len = len(indices)
        distribution = torch.zeros(seq_len, num_codons, device=self.device)
        
        for pos, idx in enumerate(indices):
            if idx < num_codons:
                distribution[pos, idx] = 1.0
        
        return distribution
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取统计信息"""
        stats = {
            'iteration_count': self.iteration_count,
            'unique_sequences': len(self.history_hashes),
            'consecutive_failures': self.consecutive_failures,
            'last_gamma': self.last_gamma
        }

        if self.enable_caching:
            stats['cache_hits'] = self._cache_hits
            stats['cache_misses'] = self._cache_misses
            stats['cache_hit_rate'] = self._cache_hits / (self._cache_hits + self._cache_misses + 1e-10)

        return stats

    # ========== 辅助方法 ==========

    def _calculate_k_positions(self, seq_length: int) -> int:
        """计算优化位置数"""
        return min(self.MAX_K_POSITIONS,
                  max(self.MIN_K_POSITIONS, seq_length // self.K_POSITION_DIVISOR))

    def _calculate_k_probability_positions(self, seq_length: int) -> int:
        """计算概率优化位置数"""
        return min(self.MAX_K_PROBABILITY,
                  max(self.MIN_K_PROBABILITY, seq_length // self.K_PROBABILITY_DIVISOR))

    def _calculate_perturb_positions(self, num_positions: int) -> int:
        """计算扰动位置数"""
        min_perturb = max(1, int(num_positions * self.MIN_PERTURB_RATIO))
        max_perturb = min(int(num_positions * self.MAX_PERTURB_RATIO), self.MAX_PERTURB_POSITIONS)
        return max(min_perturb, min(max_perturb, num_positions // 10))

    def _normalize_probabilities(self, probs: np.ndarray) -> np.ndarray:
        """
        安全地标准化概率分布

        Args:
            probs: 未标准化的概率数组

        Returns:
            标准化后的概率分布
        """
        # 处理空数组
        if len(probs) == 0:
            return probs

        # 处理NaN和Inf
        probs = np.nan_to_num(probs, nan=0.0, posinf=1.0, neginf=0.0)

        # 处理负数
        probs = np.maximum(probs, 0)

        # 计算总和
        prob_sum = probs.sum()

        # 处理零和或非常小的总和
        if prob_sum < self.NUMERICAL_EPSILON:
            # 使用均匀分布
            return np.ones(len(probs)) / len(probs)

        # 标准化
        normalized = probs / prob_sum

        # 确保和为1（处理浮点误差）
        normalized = normalized / normalized.sum()

        return normalized

    def _validate_inputs(self,
                        pi_accessibility: torch.Tensor,
                        amino_acid_sequence: str,
                        target_cai: float) -> None:
        """
        验证输入参数

        Raises:
            ValueError: 如果输入无效
        """
        # 验证CAI目标值
        if not 0 < target_cai <= 1:
            raise ValueError(f"target_cai must be in (0, 1], got {target_cai}")

        # 验证氨基酸序列
        if not amino_acid_sequence:
            raise ValueError("amino_acid_sequence cannot be empty")

        # 验证pi_accessibility张量
        if pi_accessibility.dim() not in [2, 3]:
            raise ValueError(f"pi_accessibility must be 2D or 3D tensor, got {pi_accessibility.dim()}D")

        # 检查是否包含NaN或Inf
        if torch.isnan(pi_accessibility).any() or torch.isinf(pi_accessibility).any():
            logger.warning("pi_accessibility contains NaN or Inf values, will be handled")

    def clear_cache(self) -> None:
        """清空缓存"""
        if self.enable_caching:
            self._cai_cache.clear()
            self._prob_cache.clear()
            self._cache_hits = 0
            self._cache_misses = 0
            logger.debug("Cache cleared")