"""



"""

from abc import ABC, abstractmethod
from typing import Dict, Optional, Tuple, Any
import torch
import time
import logging

logger = logging.getLogger(__name__)


class BaseCAIOptimizer(ABC):
    """

    

    """
    
    def __init__(self, 
                 species: str = 'ecoli_bl21de3',
                 device: Optional[torch.device] = None,
                 amino_acid_sequence: Optional[str] = None):
        """

        
        Args:



        """
        self.species = species
        self.device = device if device is not None else torch.device('cpu')
        self.amino_acid_sequence = amino_acid_sequence
        

        self.optimization_count = 0
        self.total_time = 0.0
        
    @abstractmethod
    def optimize(self, 
                 pi_accessibility: torch.Tensor,
                 target_cai: float = 0.8,
                 **kwargs) -> Tuple[torch.Tensor, Dict[str, Any]]:
        """

        

        
        Args:



            
        Returns:

        """
        pass
    
    def validate_result(self, 
                       distribution: torch.Tensor,
                       target_cai: float,
                       computed_cai: float) -> bool:
        """

        
        Args:



            
        Returns:

        """
        if computed_cai < target_cai:
            logger.warning(f"CAI constraint not satisfied: {computed_cai:.4f} < {target_cai:.4f}")
            return False
        return True
    
    def reset_statistics(self):

        self.optimization_count = 0
        self.total_time = 0.0
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取性能统计信息"""
        avg_time = self.total_time / max(1, self.optimization_count)
        return {
            'optimization_count': self.optimization_count,
            'total_time': self.total_time,
            'average_time': avg_time,
            'optimizer_name': self.__class__.__name__
        }
    
    def _time_operation(self, operation_func, *args, **kwargs):

        start_time = time.time()
        result = operation_func(*args, **kwargs)
        elapsed_time = time.time() - start_time
        
        self.optimization_count += 1
        self.total_time += elapsed_time
        
        return result, elapsed_time