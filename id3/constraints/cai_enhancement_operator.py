"""
CAI Enhancement Operator - Router Implementation







Mathematical Formulation from Paper:
    Ψ_{π→τ}(π, τ_CAI) = arg max_S P(S|π) s.t. CAI(S) ≥ τ_CAI
"""

import torch
import logging
from typing import Dict, Tuple, Optional, Any
from pathlib import Path


from id3.optimizers.cai import BinarySearchCAIOptimizer, SADOOptimizer
from id3.optimizers.cai.hybrid_bs_sado import HybridBSSADOOptimizer

logger = logging.getLogger(__name__)


class CAIEnhancementOperator:
    """

    





    


        operator = CAIEnhancementOperator()
        

        operator = CAIEnhancementOperator(method='binary_search')
        

        operator = CAIEnhancementOperator(method='sado')
    """
    

    OPTIMIZERS = {
        'binary_search': BinarySearchCAIOptimizer,
        'sado': SADOOptimizer,
        'hybrid_bs_sado': HybridBSSADOOptimizer,
        'incremental': None,
    }
    

    _amino_acid_sequence_cache = {}
    _amino_acid_weights_cache = {}
    
    def __init__(self, 
                 method: str = 'incremental',
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        """

        
        Args:








        """
        if method not in self.OPTIMIZERS:
            raise ValueError(
                f"Unknown optimization method: {method}. "
                f"Available methods: {list(self.OPTIMIZERS.keys())}"
            )
        
        self.method = method
        self.species = species
        self.device = device if device is not None else torch.device('cpu')
        self.amino_acid_sequence = amino_acid_sequence
        

        if method == 'incremental':

            from id3.optimizers.cai.incremental import IncrementalCAIOptimizer
            self.optimizer = IncrementalCAIOptimizer(
                species=species,
                device=self.device,
                amino_acid_sequence=amino_acid_sequence
            )
        else:
            optimizer_class = self.OPTIMIZERS[method]
            self.optimizer = optimizer_class(
                species=species,
                device=self.device,
                amino_acid_sequence=amino_acid_sequence
            )
        

        if hasattr(self.optimizer, 'wi_table'):
            self.wi_table = self.optimizer.wi_table
        if hasattr(self.optimizer, 'weights_tensor'):
            self.weights_tensor = self.optimizer.weights_tensor
        

        self._sync_cache_state()
        
        logger.info(f"CAIEnhancementOperator initialized with method: {method}")
    
    def _sync_cache_state(self):

        if hasattr(self.optimizer, 'max_achievable_cai'):
            self.max_achievable_cai = self.optimizer.max_achievable_cai
        else:
            self.max_achievable_cai = None
            
        if hasattr(self.optimizer, 'cai_optimal_distribution'):
            self.cai_optimal_distribution = self.optimizer.cai_optimal_distribution
        else:
            self.cai_optimal_distribution = None
            
        if hasattr(self.optimizer, 'cached_amino_acid_sequence'):
            self.cached_amino_acid_sequence = self.optimizer.cached_amino_acid_sequence
        else:
            self.cached_amino_acid_sequence = None
    
    def apply_cai_enhancement(self,
                             pi_accessibility: torch.Tensor,
                             amino_acid_sequence: str,
                             valid_codon_mask: torch.Tensor,
                             codon_indices: torch.Tensor,
                             target_cai: float,
                             **kwargs) -> Tuple[torch.Tensor, Dict]:
        """
        应用CAI增强操作符 Ψ_{π→τ}
        
        这是论文中Ψ_{π→τ}操作符的主要入口点。
        根据初始化时选择的方法，路由到相应的优化器。
        
        Args:
            pi_accessibility: 可及性优化的分布 π
            amino_acid_sequence: 目标氨基酸序列
            valid_codon_mask: 有效密码子掩码
            codon_indices: 密码子索引（某些优化器可能忽略）
            target_cai: 目标CAI值 τ_CAI
            **kwargs: 其他方法特定的参数
            
        Returns:
            Tuple of (离散分布, 元数据字典)
        """

        if self.method == 'sado':

            discrete_distribution, metadata = self.optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                **kwargs
            )
        elif self.method in ['binary_search', 'hybrid_bs_sado']:

            discrete_distribution, metadata = self.optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )
        else:

            discrete_distribution, metadata = self.optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )
        

        self._sync_cache_state()
        

        if 'method' not in metadata:
            metadata['method'] = self.method
        
        return discrete_distribution, metadata
    
    def enhance(self, 
                pi_accessibility: torch.Tensor,
                target_cai: float = 0.8,
                **kwargs) -> Tuple[torch.Tensor, Dict]:
        """
        简化的增强接口（向后兼容）
        
        Args:
            pi_accessibility: 可及性优化的分布
            target_cai: 目标CAI值
            **kwargs: 其他参数
            
        Returns:
            Tuple of (优化后的分布, 元数据)
        """
        return self.optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            **kwargs
        )
    
    def reset(self):
        """
        重置优化器状态
        
        某些优化器（如SADO）维护内部状态，需要在新的优化序列开始时重置。
        """
        if hasattr(self.optimizer, 'reset'):
            self.optimizer.reset()
            logger.debug(f"Optimizer {self.method} state reset")
    
    def get_statistics(self) -> Dict[str, Any]:
        """
        获取优化器的性能统计
        
        Returns:
            包含性能指标的字典
        """
        stats = {
            'method': self.method,
            'species': self.species,
        }
        

        if hasattr(self.optimizer, 'get_statistics'):
            optimizer_stats = self.optimizer.get_statistics()
            stats.update(optimizer_stats)
        

        if self.method == 'sado' and hasattr(self.optimizer, 'get_diversity_stats'):
            diversity_stats = self.optimizer.get_diversity_stats()
            stats['diversity'] = diversity_stats
        
        return stats
    
    def switch_method(self, new_method: str):
        """
        动态切换优化方法
        
        Args:
            new_method: 新的优化方法名称
        """
        if new_method not in self.OPTIMIZERS:
            raise ValueError(f"Unknown method: {new_method}")
        
        if new_method != self.method:
            logger.info(f"Switching optimization method from {self.method} to {new_method}")
            

            optimizer_class = self.OPTIMIZERS[new_method]
            self.optimizer = optimizer_class(
                species=self.species,
                device=self.device,
                amino_acid_sequence=self.amino_acid_sequence
            )
            
            self.method = new_method
            self._sync_cache_state()
    

    
    def _load_or_compute_amino_acid_cache(self, amino_acid_sequence: str):
        """向后兼容：加载或计算氨基酸序列缓存"""
        if hasattr(self.optimizer, '_load_or_compute_amino_acid_cache'):
            self.optimizer._load_or_compute_amino_acid_cache(amino_acid_sequence)
            self._sync_cache_state()
    
    def _precompute_amino_acid_weights(self):

        if hasattr(self.optimizer, '_precompute_amino_acid_weights'):
            self.optimizer._precompute_amino_acid_weights()
    
    def discretize_distribution(self, distribution: torch.Tensor, valid_mask: torch.Tensor) -> torch.Tensor:
        """向后兼容：离散化分布"""
        from id3.optimizers.cai.utils import discretize_distribution
        return discretize_distribution(distribution, valid_mask)
    
    def compute_discrete_cai(self, 
                            discrete_dist: torch.Tensor,
                            amino_acid_sequence: str,
                            valid_codon_mask: torch.Tensor,
                            codon_indices: torch.Tensor) -> float:

        if hasattr(self.optimizer, '_compute_cai_from_indices'):
            indices = discrete_dist.argmax(dim=-1)
            return self.optimizer._compute_cai_from_indices(indices, amino_acid_sequence)
        return 0.0
    
    def interpolate_distributions(self,
                                 dist1: torch.Tensor,
                                 dist2: torch.Tensor,
                                 gamma: float) -> torch.Tensor:
        """向后兼容：分布插值"""
        from id3.optimizers.cai.utils import interpolate_distributions
        return interpolate_distributions(dist1, dist2, gamma)
    
    @classmethod
    def available_methods(cls) -> list:

        return list(cls.OPTIMIZERS.keys())
    
    def __repr__(self) -> str:
        return (
            f"CAIEnhancementOperator(method='{self.method}', "
            f"species='{self.species}', device={self.device})"
        )