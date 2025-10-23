"""



"""

import torch
import numpy as np
from typing import Optional, List, Dict, Any, Union, Tuple, Callable
import warnings
import functools
import logging
from enum import Enum
from contextlib import contextmanager


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class ErrorSeverity(Enum):

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class ValidationError(Exception):
    """验证错误"""
    def __init__(self, message: str, severity: ErrorSeverity = ErrorSeverity.ERROR):
        self.severity = severity
        super().__init__(message)


class DeviceMismatchError(ValidationError):

    pass


class SequenceValidationError(ValidationError):
    """序列验证错误"""
    pass


class CAIOptimizationError(Exception):

    pass


class ConstraintError(Exception):
    """约束计算错误"""
    pass


class RecoveryStrategy:

    RETRY = "retry"
    FALLBACK = "fallback"
    REDUCE_COMPLEXITY = "reduce_complexity"
    SKIP = "skip"


class TensorValidator:
    """张量验证器"""
    
    @staticmethod
    def validate_tensor_shape(tensor: torch.Tensor, 
                            expected_shape: Optional[Tuple] = None,
                            min_dims: Optional[int] = None,
                            max_dims: Optional[int] = None) -> bool:
        """

        
        Args:




            
        Returns:

            
        Raises:

        """
        if not isinstance(tensor, torch.Tensor):
            raise ValidationError(f"Expected torch.Tensor, got {type(tensor)}")
        

        dims = tensor.dim()
        if min_dims is not None and dims < min_dims:
            raise ValidationError(f"Tensor has {dims} dimensions, expected at least {min_dims}")
        
        if max_dims is not None and dims > max_dims:
            raise ValidationError(f"Tensor has {dims} dimensions, expected at most {max_dims}")
        

        if expected_shape is not None:
            actual_shape = tensor.shape
            if len(expected_shape) != len(actual_shape):
                raise ValidationError(f"Shape mismatch: expected {len(expected_shape)} dims, got {len(actual_shape)}")
            
            for i, (expected, actual) in enumerate(zip(expected_shape, actual_shape)):
                if expected is not None and expected != actual:
                    raise ValidationError(f"Shape mismatch at dim {i}: expected {expected}, got {actual}")
        
        return True
    
    @staticmethod
    def validate_tensor_dtype(tensor: torch.Tensor, 
                            allowed_dtypes: List[torch.dtype]) -> bool:
        """

        
        Args:


            
        Returns:

        """
        if tensor.dtype not in allowed_dtypes:
            raise ValidationError(f"Tensor dtype {tensor.dtype} not in allowed types {allowed_dtypes}")
        return True
    
    @staticmethod
    def validate_tensor_device(tensor: torch.Tensor, 
                             expected_device: Optional[torch.device] = None) -> bool:
        """

        
        Args:


            
        Returns:

        """
        if expected_device is not None and tensor.device != expected_device:
            raise DeviceMismatchError(f"Tensor on device {tensor.device}, expected {expected_device}")
        return True


class SequenceValidator:

    
    VALID_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY*')
    VALID_NUCLEOTIDES = set('ACGTU')
    
    @staticmethod
    def validate_amino_acid_sequence(sequence: str) -> bool:
        """
        验证氨基酸序列
        
        Args:
            sequence: 氨基酸序列字符串
            
        Returns:
            是否验证通过
        """
        if not isinstance(sequence, str):
            raise SequenceValidationError("Amino acid sequence must be a string")
        
        if len(sequence) == 0:
            raise SequenceValidationError("Amino acid sequence cannot be empty")
        
        invalid_chars = set(sequence.upper()) - SequenceValidator.VALID_AMINO_ACIDS
        if invalid_chars:
            raise SequenceValidationError(f"Invalid amino acid characters: {invalid_chars}")
        
        return True
    
    @staticmethod
    def validate_nucleotide_sequence(sequence: str) -> bool:
        """
        验证核苷酸序列
        
        Args:
            sequence: 核苷酸序列字符串
            
        Returns:
            是否验证通过
        """
        if not isinstance(sequence, str):
            raise SequenceValidationError("Nucleotide sequence must be a string")
        
        if len(sequence) == 0:
            raise SequenceValidationError("Nucleotide sequence cannot be empty")
        
        if len(sequence) % 3 != 0:
            warnings.warn("Nucleotide sequence length not divisible by 3", UserWarning)
        
        invalid_chars = set(sequence.upper()) - SequenceValidator.VALID_NUCLEOTIDES
        if invalid_chars:
            raise SequenceValidationError(f"Invalid nucleotide characters: {invalid_chars}")
        
        return True


class DeviceManager:
    """设备管理器"""
    
    @staticmethod
    def ensure_same_device(*tensors: torch.Tensor) -> torch.device:
        """

        
        Args:

            
        Returns:

            
        Raises:

        """
        if not tensors:
            return torch.device('cpu')
        
        devices = [t.device for t in tensors]
        first_device = devices[0]
        
        for i, device in enumerate(devices[1:], 1):
            if device != first_device:
                raise DeviceMismatchError(f"Device mismatch: tensor 0 on {first_device}, tensor {i} on {device}")
        
        return first_device
    
    @staticmethod
    def move_to_device(tensors: Union[torch.Tensor, List[torch.Tensor]], 
                      target_device: torch.device) -> Union[torch.Tensor, List[torch.Tensor]]:
        """

        
        Args:


            
        Returns:

        """
        if isinstance(tensors, torch.Tensor):
            return tensors.to(target_device)
        else:
            return [t.to(target_device) for t in tensors]
    
    @staticmethod
    def get_optimal_device() -> torch.device:

        if torch.cuda.is_available():

            best_gpu = 0
            max_memory = 0
            for i in range(torch.cuda.device_count()):
                memory = torch.cuda.get_device_properties(i).total_memory
                if memory > max_memory:
                    max_memory = memory
                    best_gpu = i
            return torch.device(f'cuda:{best_gpu}')
        else:
            return torch.device('cpu')


class ErrorRecoveryManager:
    """错误恢复管理器"""
    
    def __init__(self, max_retries: int = 3, fallback_device: torch.device = None):
        self.max_retries = max_retries
        self.fallback_device = fallback_device or torch.device('cpu')
        self.error_counts = {}
    
    def with_retry(self, func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_error = None
            
            for attempt in range(self.max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except (RuntimeError, ValueError, ValidationError) as e:
                    last_error = e
                    error_type = type(e).__name__
                    

                    self.error_counts[error_type] = self.error_counts.get(error_type, 0) + 1
                    
                    if attempt < self.max_retries:
                        print(f"⚠️ Attempt {attempt + 1} failed: {e}. Retrying...")
                        

                        if "out of memory" in str(e).lower():
                            torch.cuda.empty_cache()
                        elif "device" in str(e).lower():

                            kwargs = self._move_kwargs_to_device(kwargs, self.fallback_device)
                    else:
                        print(f"❌ All {self.max_retries + 1} attempts failed")
                        break
            
            raise last_error
        
        return wrapper
    
    def _move_kwargs_to_device(self, kwargs: Dict, device: torch.device) -> Dict:
        """将kwargs中的张量移动到指定设备"""
        new_kwargs = {}
        for key, value in kwargs.items():
            if isinstance(value, torch.Tensor):
                new_kwargs[key] = value.to(device)
            else:
                new_kwargs[key] = value
        return new_kwargs


def validate_cai_inputs(codon_probs: torch.Tensor, 
                       amino_acid_sequence: Optional[str] = None) -> bool:
    """

    
    Args:


        
    Returns:

    """

    TensorValidator.validate_tensor_shape(codon_probs, min_dims=2, max_dims=3)
    TensorValidator.validate_tensor_dtype(codon_probs, [torch.float32, torch.float64])
    

    if codon_probs.shape[-1] != 64:
        raise ValidationError(f"Expected 64 codons, got {codon_probs.shape[-1]}")
    

    prob_sums = torch.sum(codon_probs, dim=-1)
    if not torch.allclose(prob_sums, torch.ones_like(prob_sums), atol=1e-6):
        warnings.warn("Codon probabilities do not sum to 1", UserWarning)
    

    if amino_acid_sequence is not None:
        SequenceValidator.validate_amino_acid_sequence(amino_acid_sequence)
        

        expected_codons = len(amino_acid_sequence)
        if codon_probs.dim() == 2:
            actual_codons = codon_probs.shape[0]
        else:
            actual_codons = codon_probs.shape[1]
        
        if actual_codons != expected_codons:
            raise ValidationError(f"Sequence length mismatch: {expected_codons} amino acids, {actual_codons} codon positions")
    
    return True


def validate_optimization_inputs(logits: torch.Tensor,
                               amino_acid_sequence: str,
                               utr_config: Dict) -> bool:
    """

    
    Args:



        
    Returns:

    """

    TensorValidator.validate_tensor_shape(logits, min_dims=2, max_dims=2)
    if logits.shape[-1] != 4:
        raise ValidationError(f"Expected 4 nucleotides, got {logits.shape[-1]}")
    

    SequenceValidator.validate_amino_acid_sequence(amino_acid_sequence)
    

    required_keys = ['utr5_template', 'utr3_template']
    for key in required_keys:
        if key not in utr_config:
            raise ValidationError(f"Missing required UTR config key: {key}")
    

    expected_length = len(utr_config['utr5_template']) + len(amino_acid_sequence) * 3 + len(utr_config['utr3_template'])
    actual_length = logits.shape[0]
    
    if actual_length != expected_length:
        raise ValidationError(f"Length mismatch: expected {expected_length}, got {actual_length}")
    
    return True



def with_error_handling(validator_func: Optional[callable] = None, 
                       recovery_manager: Optional[ErrorRecoveryManager] = None):
    """

    
    Args:


    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:

                if validator_func is not None:
                    validator_func(*args, **kwargs)
                

                if recovery_manager is not None:
                    return recovery_manager.with_retry(func)(*args, **kwargs)
                else:
                    return func(*args, **kwargs)
                    
            except ValidationError as e:
                print(f"🚫 Validation Error: {e}")
                raise
            except Exception as e:
                print(f"💥 Unexpected Error: {e}")
                raise
        
        return wrapper
    return decorator


# =============================================================================

# =============================================================================

def robust_cai_optimization(max_gpu_retries: int = 2):
    """

    
    Args:

    

    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            original_kwargs = kwargs.copy()
            
            for attempt in range(max_gpu_retries + 1):
                try:
                    return func(*args, **kwargs)
                    
                except torch.cuda.OutOfMemoryError as e:
                    if attempt < max_gpu_retries:
                        logger.warning(f"GPU内存不足，清理后重试 ({attempt + 1}/{max_gpu_retries}): {e}")
                        if torch.cuda.is_available():
                            torch.cuda.empty_cache()

                        kwargs = original_kwargs.copy()
                    else:
                        logger.error(f"GPU内存不足，已重试{max_gpu_retries}次仍失败")
                        raise CAIOptimizationError(f"GPU内存不足: {e}")
                    
                except Exception as e:

                    logger.error(f"CAI优化失败: {e}")
                    raise CAIOptimizationError(f"CAI优化失败: {e}")
            

            raise CAIOptimizationError("意外的执行路径")
        
        return wrapper
    return decorator


def constraint_error_handler(allow_cpu_fallback: bool = False):
    """

    
    Args:

    

    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
                
            except torch.cuda.OutOfMemoryError as e:
                if allow_cpu_fallback:
                    logger.warning(f"GPU内存不足，切换到CPU: {e}")

                    if 'device' in kwargs:
                        kwargs['device'] = 'cpu'
                    elif len(args) > 0 and hasattr(args[0], 'device_str'):
                        args[0].device_str = 'cpu'
                        if hasattr(args[0], 'to'):
                            args[0].to('cpu')
                    
                    torch.cuda.empty_cache()
                    return func(*args, **kwargs)
                else:
                    logger.error(f"约束计算GPU内存不足: {e}")
                    raise ConstraintError(f"GPU内存不足，请减少批大小或使用CPU: {e}")
                    
            except Exception as e:
                logger.error(f"约束计算失败: {e}")
                raise ConstraintError(f"约束计算失败: {e}")
        
        return wrapper
    return decorator


@contextmanager
def safe_tensor_operations(fallback_device: str = 'cpu'):
    """

    
    Args:

    """
    try:
        yield
    except torch.cuda.OutOfMemoryError:
        logger.warning(f"CUDA内存不足，切换到 {fallback_device}")
        torch.cuda.empty_cache()

        raise
    except RuntimeError as e:
        if "CUDA" in str(e):
            logger.warning(f"CUDA运行时错误: {e}，切换到 {fallback_device}")
            torch.cuda.empty_cache()
            raise
        else:
            raise


def safe_experiment_runner(experiment_func: Callable, 
                         config: Dict[str, Any],
                         fail_fast: bool = True) -> Dict[str, Any]:
    """

    
    Args:



    
    Returns:

    """
    results = {
        'success': [],
        'failed': [],
        'errors': []
    }
    

    readonly_config = config.copy()
    
    try:
        result = experiment_func(readonly_config)
        results['success'].append(result)
        
    except CAIOptimizationError as e:
        error_msg = f"CAI优化错误: {e}"
        logger.error(error_msg)
        results['errors'].append(error_msg)
        results['failed'].append(config.get('variant', 'unknown'))
        
        if fail_fast:
            raise
                
    except ConstraintError as e:
        error_msg = f"约束计算错误: {e}"
        logger.error(error_msg)
        results['errors'].append(error_msg)
        results['failed'].append(config.get('variant', 'unknown'))
        
        if fail_fast:
            raise
            
    except Exception as e:
        error_msg = f"实验执行错误: {e}"
        logger.error(error_msg)
        results['errors'].append(error_msg)
        results['failed'].append(config.get('variant', 'unknown'))
        
        if fail_fast:
            raise
    
    return results








robust_cai = robust_cai_optimization(max_gpu_retries=2)
strict_constraint = constraint_error_handler(allow_cpu_fallback=False)



def log_experiment_error(variant: str, error: Exception, context: Dict[str, Any] = None):


    if context:



def check_system_resources() -> Dict[str, Union[bool, float]]:
    """检查系统资源状态"""
    resources = {}
    

    if torch.cuda.is_available():
        gpu_memory = torch.cuda.get_device_properties(0).total_memory
        gpu_memory_used = torch.cuda.memory_allocated(0)
        resources['gpu_available'] = True
        resources['gpu_memory_free_pct'] = 1.0 - (gpu_memory_used / gpu_memory)
    else:
        resources['gpu_available'] = False
        resources['gpu_memory_free_pct'] = 0.0
    

    resources['cpu_available'] = True
    
    return resources


def get_error_details(error: Exception, context: Dict[str, Any]) -> Dict[str, str]:

    details = {
        'error_type': type(error).__name__,
        'error_message': str(error),
        'context': str(context) if context else 'None'
    }
    
    if isinstance(error, torch.cuda.OutOfMemoryError):


    elif isinstance(error, CAIOptimizationError):


    elif isinstance(error, ConstraintError):


    else:


    
    return details