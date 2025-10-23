"""




"""

import torch
import numpy as np
import hashlib
import logging
from typing import Dict, Tuple, Optional, Set
from ..base import BaseCAIOptimizer
from .binary_search import BinarySearchCAIOptimizer
from .sado import SADOOptimizer

logger = logging.getLogger(__name__)


class HybridBSSADOOptimizer(BaseCAIOptimizer):
    """

    





    """
    
    def __init__(self, 
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None,
                 repetition_threshold: int = 3):
        """

        
        Args:




        """
        super().__init__(species, device, amino_acid_sequence)
        

        self.bs_optimizer = BinarySearchCAIOptimizer(
            species=species,
            device=device,
            amino_acid_sequence=amino_acid_sequence
        )
        
        self.sado_optimizer = SADOOptimizer(
            species=species,
            device=device,
            amino_acid_sequence=amino_acid_sequence
        )
        

        self.sequence_history: Set[str] = set()
        self.recent_sequences = []
        self.repetition_threshold = repetition_threshold
        self.repetition_count = 0
        self.last_sequence_hash = None
        

        self.bs_calls = 0
        self.sado_calls = 0
        self.switch_count = 0
        
        logger.info(f"HybridBSSADOOptimizer initialized with threshold={repetition_threshold}")
    
    def _compute_sequence_hash(self, indices: torch.Tensor) -> str:

        if indices.dim() > 1:
            indices = indices.argmax(dim=-1)
        indices_np = indices.cpu().numpy()
        return hashlib.md5(indices_np.tobytes()).hexdigest()
    
    def _check_repetition(self, sequence_hash: str) -> bool:
        """检查序列是否重复"""

        if sequence_hash == self.last_sequence_hash:
            self.repetition_count += 1
        else:
            self.repetition_count = 0
            self.last_sequence_hash = sequence_hash
        

        is_duplicate = sequence_hash in self.sequence_history
        

        self.sequence_history.add(sequence_hash)
        self.recent_sequences.append(sequence_hash)
        

        if len(self.recent_sequences) > 10:
            self.recent_sequences.pop(0)
        
        return is_duplicate or (self.repetition_count >= self.repetition_threshold)
    
    def optimize(self,
                 pi_accessibility: torch.Tensor,
                 target_cai: float = 0.8,
                 amino_acid_sequence: Optional[str] = None,
                 valid_codon_mask: Optional[torch.Tensor] = None,
                 **kwargs) -> Tuple[torch.Tensor, Dict]:
        """

        
        Args:





            
        Returns:

        """
        if amino_acid_sequence is None:
            amino_acid_sequence = self.amino_acid_sequence
            

        self.bs_calls += 1
        
        try:

            bs_result, bs_metadata = self.bs_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                valid_codon_mask=valid_codon_mask,
                **kwargs
            )
            

            bs_gamma = bs_metadata.get('optimal_gamma', 0.5)
            

            if bs_result.dim() > 1:
                indices = bs_result.argmax(dim=-1)
            else:
                indices = bs_result
                
            seq_hash = self._compute_sequence_hash(indices)
            

            if self._check_repetition(seq_hash):
                logger.info(f"Binary Search产生重复序列，切换到SADO (gamma={bs_gamma:.3f})")
                self.switch_count += 1
                self.sado_calls += 1
                


                sado_kwargs = {k: v for k, v in kwargs.items() if k != 'valid_codon_mask'}
                sado_kwargs['gamma'] = bs_gamma
                
                sado_result, sado_metadata = self.sado_optimizer.optimize(
                    pi_accessibility=pi_accessibility,
                    target_cai=target_cai,
                    amino_acid_sequence=amino_acid_sequence,
                    **sado_kwargs
                )
                

                sado_metadata['method'] = 'hybrid_bs_sado'
                sado_metadata['switched_to_sado'] = True
                sado_metadata['bs_gamma'] = bs_gamma
                sado_metadata['switch_reason'] = 'repetition_detected'
                sado_metadata['bs_calls'] = self.bs_calls
                sado_metadata['sado_calls'] = self.sado_calls
                sado_metadata['switch_count'] = self.switch_count
                
                return sado_result, sado_metadata
            
            else:

                bs_metadata['method'] = 'hybrid_bs_sado'
                bs_metadata['switched_to_sado'] = False
                bs_metadata['bs_calls'] = self.bs_calls
                bs_metadata['sado_calls'] = self.sado_calls
                bs_metadata['switch_count'] = self.switch_count
                
                return bs_result, bs_metadata
                
        except Exception as e:
            logger.error(f"Binary Search失败: {e}，切换到SADO")
            self.switch_count += 1
            self.sado_calls += 1
            

            sado_kwargs = {k: v for k, v in kwargs.items() if k != 'valid_codon_mask'}
            sado_kwargs['gamma'] = 0.5
            
            sado_result, sado_metadata = self.sado_optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_acid_sequence,
                **sado_kwargs
            )
            
            sado_metadata['method'] = 'hybrid_bs_sado'
            sado_metadata['switched_to_sado'] = True
            sado_metadata['switch_reason'] = 'bs_error'
            sado_metadata['error'] = str(e)
            sado_metadata['bs_calls'] = self.bs_calls
            sado_metadata['sado_calls'] = self.sado_calls
            sado_metadata['switch_count'] = self.switch_count
            
            return sado_result, sado_metadata
    
    def reset(self):

        self.sequence_history.clear()
        self.recent_sequences.clear()
        self.repetition_count = 0
        self.last_sequence_hash = None
        self.bs_calls = 0
        self.sado_calls = 0
        self.switch_count = 0
        

        if hasattr(self.bs_optimizer, 'reset'):
            self.bs_optimizer.reset()
        if hasattr(self.sado_optimizer, 'reset'):
            self.sado_optimizer.reset()
    
    def get_statistics(self) -> Dict:
        """获取统计信息"""
        return {
            'method': 'hybrid_bs_sado',
            'bs_calls': self.bs_calls,
            'sado_calls': self.sado_calls,
            'switch_count': self.switch_count,
            'switch_rate': self.switch_count / max(1, self.bs_calls),
            'unique_sequences': len(self.sequence_history),
            'repetition_threshold': self.repetition_threshold
        }