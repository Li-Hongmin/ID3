"""
CAI Constraint Implementation for ID3 Framework (Improved Version)



Mathematical Framework:
    Enhanced ID3: Θ → Π → P → Ψ_CAI → S_CAI → f_model + C_CAI → L → ∇ → Θ^(t+1)
    
where:
    - Ψ_CAI: CAI-enhanced discretization with dual optimization
    - C_CAI: Combined constraint penalty (original + CAI)
    - S_CAI: CAI-optimized sequence in codon space
"""

import torch
import torch.nn as nn
from typing import Dict, Tuple, Optional, Any
import logging


from id3.cai.unified_calculator import UnifiedCAICalculator, compute_cai
from id3.utils.logging_config import get_logger

logger = get_logger(__name__)


class CAIEnhancedPsiFunction:
    """
    CAI-Enhanced Psi Function: Probability-to-Sequence with CAI Optimization
    

    
    Design Principles:
    1. **Non-intrusive**: Only activates when CAI parameters are provided
    2. **Dual-path**: Works with original constraints + CAI objective
    3. **Codon-aware**: Operates in codon space for precise CAI control
    4. **Efficient**: Uses unified CAI calculator with caching
    """
    
    def __init__(self, 
                 amino_acid_sequence: str,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 optimizer_type: str = 'binary_search'):
        """

        
        Args:




        """
        self.amino_acid_sequence = amino_acid_sequence
        self.species = species
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.optimizer_type = optimizer_type
        

        self.cai_calculator = UnifiedCAICalculator(species=species, device=self.device)
        

        self.max_achievable_cai = self.cai_calculator.get_max_achievable_cai(amino_acid_sequence)
        logger.info(f"Max achievable CAI for sequence: {self.max_achievable_cai:.4f}")
        
    def apply(self, 
              codon_probabilities: torch.Tensor,
              target_cai: float = 0.8,
              beta: float = 1.0,
              enable_cai: bool = True,
              **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """

        
        Args:




            
        Returns:

        """
        metadata = {
            'target_cai': target_cai,
            'max_achievable_cai': self.max_achievable_cai,
            'enable_cai': enable_cai,
            'optimizer_type': self.optimizer_type
        }
        
        if not enable_cai:

            if beta == 1.0:

                discrete_indices = torch.argmax(codon_probabilities, dim=-1)
                discrete_sequence = self._indices_to_one_hot(discrete_indices, codon_probabilities.shape[-1])
            else:

                discrete_sequence = codon_probabilities
            
            metadata['actual_cai'] = 0.0
            metadata['constraint_satisfied'] = False
            
        else:

            discrete_sequence, cai_metadata = self._optimize_with_cai(
                codon_probabilities, 
                target_cai,
                beta
            )
            metadata.update(cai_metadata)
        
        return discrete_sequence, metadata
    
    def _optimize_with_cai(self,
                           codon_probabilities: torch.Tensor,
                           target_cai: float,
                           beta: float) -> Tuple[torch.Tensor, Dict]:
        """

        

        """


        
        if self.optimizer_type == 'binary_search':
            from id3.optimizers.cai.binary_search import BinarySearchCAIOptimizer
            optimizer = BinarySearchCAIOptimizer(
                species=self.species,
                device=self.device,
                amino_acid_sequence=self.amino_acid_sequence
            )
        elif self.optimizer_type == 'sado':
            from id3.optimizers.cai.sado import SADOOptimizer
            optimizer = SADOOptimizer(
                species=self.species,
                device=self.device,
                amino_acid_sequence=self.amino_acid_sequence
            )
        else:
            raise ValueError(f"Unknown optimizer type: {self.optimizer_type}")
        

        from id3.utils.constants import amino_acids_to_codons
        seq_len = len(self.amino_acid_sequence)
        num_codons = codon_probabilities.shape[-1]
        valid_codon_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=self.device)
        
        for pos, aa in enumerate(self.amino_acid_sequence):
            if aa in amino_acids_to_codons:
                num_valid = len(amino_acids_to_codons[aa])
                valid_codon_mask[pos, :min(num_valid, num_codons)] = True
        

        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=codon_probabilities,
            target_cai=target_cai,
            amino_acid_sequence=self.amino_acid_sequence,
            valid_codon_mask=valid_codon_mask
        )
        

        if beta == 1.0:
            # Straight-Through Estimator
            soft_dist = codon_probabilities
            hard_dist = optimized_dist

            discrete_sequence = hard_dist + soft_dist - soft_dist.detach()
        else:

            discrete_sequence = optimized_dist
        

        actual_cai = self.cai_calculator.compute_cai(
            optimized_dist,
            method='differentiable',
            amino_acid_sequence=self.amino_acid_sequence
        )
        
        metadata['actual_cai'] = actual_cai.item() if isinstance(actual_cai, torch.Tensor) else actual_cai
        metadata['constraint_satisfied'] = metadata['actual_cai'] >= target_cai
        
        return discrete_sequence, metadata
    
    def _indices_to_one_hot(self, indices: torch.Tensor, num_classes: int) -> torch.Tensor:

        shape = indices.shape + (num_classes,)
        one_hot = torch.zeros(shape, device=indices.device)
        one_hot.scatter_(-1, indices.unsqueeze(-1), 1)
        return one_hot


class CAIConstraint(nn.Module):
    """
    CAI约束模块
    
    作为ID3框架的插件，与现有约束（CPC/RAMS/Lagrangian）协同工作。
    """
    
    def __init__(self,
                 amino_acid_sequence: str,
                 target_cai: float = 0.8,
                 lambda_cai: float = 0.1,
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 optimizer_type: str = 'binary_search'):
        """
        初始化CAI约束
        
        Args:
            amino_acid_sequence: 目标氨基酸序列
            target_cai: 目标CAI值
            lambda_cai: CAI损失权重
            species: 物种
            device: 计算设备
            optimizer_type: 优化器类型
        """
        super().__init__()
        
        self.amino_acid_sequence = amino_acid_sequence
        self.target_cai = target_cai
        self.lambda_cai = lambda_cai
        self.species = species
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        

        self.psi_function = CAIEnhancedPsiFunction(
            amino_acid_sequence=amino_acid_sequence,
            species=species,
            device=self.device,
            optimizer_type=optimizer_type
        )
        

        self.cai_calculator = UnifiedCAICalculator(species=species, device=self.device)
        
        logger.info(f"CAIConstraint initialized: target={target_cai:.3f}, lambda={lambda_cai:.3f}")
    
    def forward(self, 
                codon_probabilities: torch.Tensor,
                beta: float = 1.0,
                compute_loss: bool = True) -> Dict[str, Any]:
        """
        前向传播
        
        Args:
            codon_probabilities: 密码子概率
            beta: 离散化强度
            compute_loss: 是否计算损失
            
        Returns:
            包含离散序列、损失和元数据的字典
        """

        discrete_sequence, metadata = self.psi_function.apply(
            codon_probabilities,
            target_cai=self.target_cai,
            beta=beta,
            enable_cai=True
        )
        
        result = {
            'discrete_sequence': discrete_sequence,
            'metadata': metadata
        }
        
        if compute_loss:

            cai_loss = self._compute_cai_loss(discrete_sequence)
            result['cai_loss'] = cai_loss
            result['total_loss'] = self.lambda_cai * cai_loss
        
        return result
    
    def _compute_cai_loss(self, sequence: torch.Tensor) -> torch.Tensor:
        """
        计算CAI损失
        
        使用统一的CAI计算器，保持原有的损失计算逻辑。
        """

        current_cai = self.cai_calculator.compute_cai(
            sequence,
            method='differentiable',
            amino_acid_sequence=self.amino_acid_sequence
        )
        


        cai_loss = torch.relu(self.target_cai - current_cai)
        
        return cai_loss
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取统计信息"""
        return {
            'target_cai': self.target_cai,
            'max_achievable_cai': self.psi_function.max_achievable_cai,
            'lambda_cai': self.lambda_cai,
            'species': self.species,
            'optimizer_type': self.psi_function.optimizer_type
        }


def create_cai_constraint(config: Dict[str, Any]) -> CAIConstraint:
    """

    
    Args:

        
    Returns:

    """
    return CAIConstraint(
        amino_acid_sequence=config['amino_acid_sequence'],
        target_cai=config.get('target_cai', 0.8),
        lambda_cai=config.get('lambda_cai', 0.1),
        species=config.get('species', 'ecoli_bl21de3'),
        device=config.get('device'),
        optimizer_type=config.get('cai_optimizer', 'binary_search')
    )