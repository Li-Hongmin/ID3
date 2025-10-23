"""
Error Handling Module

Provides comprehensive error handling, validation, and recovery strategies for ID3.
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
    """Error severity levels"""
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class ValidationError(Exception):
    """Validation error"""
    def __init__(self, message: str, severity: ErrorSeverity = ErrorSeverity.ERROR):
        self.severity = severity
        super().__init__(message)


class DeviceMismatchError(ValidationError):
    """Device mismatch error"""
    pass


class SequenceValidationError(ValidationError):
    """Sequence validation error"""
    pass


class CAIOptimizationError(Exception):
    """CAI optimization error"""
    pass


class ConstraintError(Exception):
    """Constraint calculation error"""
    pass


class RecoveryStrategy:
    """Error recovery strategies"""
    RETRY = "retry"
    FALLBACK = "fallback"
    REDUCE_COMPLEXITY = "reduce_complexity"
    SKIP = "skip"


class TensorValidator:
    """Tensor validator"""

    @staticmethod
    def validate_tensor_shape(tensor: torch.Tensor,
                            expected_shape: Optional[Tuple] = None,
                            min_dims: Optional[int] = None,
                            max_dims: Optional[int] = None) -> bool:
        """
        Validate tensor shape

        Args:
            tensor: Tensor to validate
            expected_shape: Expected shape (None for any dimension)
            min_dims: Minimum number of dimensions
            max_dims: Maximum number of dimensions

        Returns:
            True if validation passes

        Raises:
            ValidationError: If validation fails
        """
        if not isinstance(tensor, torch.Tensor):
            raise ValidationError(f"Expected torch.Tensor, got {type(tensor)}")

        # Check number of dimensions
        dims = tensor.dim()
        if min_dims is not None and dims < min_dims:
            raise ValidationError(f"Tensor has {dims} dimensions, expected at least {min_dims}")

        if max_dims is not None and dims > max_dims:
            raise ValidationError(f"Tensor has {dims} dimensions, expected at most {max_dims}")

        # Check expected shape
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
        Validate tensor data type

        Args:
            tensor: Tensor to validate
            allowed_dtypes: List of allowed data types

        Returns:
            True if validation passes
        """
        if tensor.dtype not in allowed_dtypes:
            raise ValidationError(f"Tensor dtype {tensor.dtype} not in allowed types {allowed_dtypes}")
        return True
    
    @staticmethod
    def validate_tensor_device(tensor: torch.Tensor,
                             expected_device: Optional[torch.device] = None) -> bool:
        """
        Validate tensor device

        Args:
            tensor: Tensor to validate
            expected_device: Expected device

        Returns:
            True if validation passes
        """
        if expected_device is not None and tensor.device != expected_device:
            raise DeviceMismatchError(f"Tensor on device {tensor.device}, expected {expected_device}")
        return True


class SequenceValidator:
    """Sequence validator for amino acids and nucleotides"""

    VALID_AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWY*')
    VALID_NUCLEOTIDES = set('ACGTU')

    @staticmethod
    def validate_amino_acid_sequence(sequence: str) -> bool:
        """
        Validate amino acid sequence

        Args:
            sequence: Amino acid sequence string

        Returns:
            True if validation passes
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
        Validate nucleotide sequence

        Args:
            sequence: Nucleotide sequence string

        Returns:
            True if validation passes
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
    """Device manager"""

    @staticmethod
    def ensure_same_device(*tensors: torch.Tensor) -> torch.device:
        """
        Ensure all tensors are on the same device

        Args:
            *tensors: Tensors to check

        Returns:
            The common device

        Raises:
            DeviceMismatchError: If tensors are on different devices
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
        Move tensors to target device

        Args:
            tensors: Tensor or list of tensors to move
            target_device: Target device

        Returns:
            Moved tensor(s)
        """
        if isinstance(tensors, torch.Tensor):
            return tensors.to(target_device)
        else:
            return [t.to(target_device) for t in tensors]
    
    @staticmethod
    def get_optimal_device() -> torch.device:
        """Get the optimal computing device (GPU with most memory or CPU)"""
        if torch.cuda.is_available():
            # Select GPU with maximum memory
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
    """Error recovery manager"""

    def __init__(self, max_retries: int = 3, fallback_device: torch.device = None):
        self.max_retries = max_retries
        self.fallback_device = fallback_device or torch.device('cpu')
        self.error_counts = {}

    def with_retry(self, func):
        """Decorator for automatic retry on failure"""
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_error = None

            for attempt in range(self.max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except (RuntimeError, ValueError, ValidationError) as e:
                    last_error = e
                    error_type = type(e).__name__

                    # Record error count
                    self.error_counts[error_type] = self.error_counts.get(error_type, 0) + 1

                    if attempt < self.max_retries:
                        print(f"⚠️ Attempt {attempt + 1} failed: {e}. Retrying...")

                        # Apply recovery strategies
                        if "out of memory" in str(e).lower():
                            torch.cuda.empty_cache()
                        elif "device" in str(e).lower():
                            # Move tensors to fallback device
                            kwargs = self._move_kwargs_to_device(kwargs, self.fallback_device)
                    else:
                        print(f"❌ All {self.max_retries + 1} attempts failed")
                        break

            raise last_error

        return wrapper

    def _move_kwargs_to_device(self, kwargs: Dict, device: torch.device) -> Dict:
        """Move tensors in kwargs to specified device"""
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
    Validate CAI calculation inputs

    Args:
        codon_probs: Codon probability tensor
        amino_acid_sequence: Optional amino acid sequence

    Returns:
        True if validation passes
    """
    # Validate tensor shape and dtype
    TensorValidator.validate_tensor_shape(codon_probs, min_dims=2, max_dims=3)
    TensorValidator.validate_tensor_dtype(codon_probs, [torch.float32, torch.float64])

    # Validate codon dimension
    if codon_probs.shape[-1] != 64:
        raise ValidationError(f"Expected 64 codons, got {codon_probs.shape[-1]}")

    # Validate probability sums
    prob_sums = torch.sum(codon_probs, dim=-1)
    if not torch.allclose(prob_sums, torch.ones_like(prob_sums), atol=1e-6):
        warnings.warn("Codon probabilities do not sum to 1", UserWarning)

    # Validate amino acid sequence if provided
    if amino_acid_sequence is not None:
        SequenceValidator.validate_amino_acid_sequence(amino_acid_sequence)

        # Check sequence length matches codon positions
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
    Validate optimization inputs

    Args:
        logits: Nucleotide logits tensor
        amino_acid_sequence: Amino acid sequence
        utr_config: UTR configuration dictionary

    Returns:
        True if validation passes
    """
    # Validate logits shape
    TensorValidator.validate_tensor_shape(logits, min_dims=2, max_dims=2)
    if logits.shape[-1] != 4:
        raise ValidationError(f"Expected 4 nucleotides, got {logits.shape[-1]}")

    # Validate amino acid sequence
    SequenceValidator.validate_amino_acid_sequence(amino_acid_sequence)

    # Validate UTR config
    required_keys = ['utr5_template', 'utr3_template']
    for key in required_keys:
        if key not in utr_config:
            raise ValidationError(f"Missing required UTR config key: {key}")

    # Validate sequence length
    expected_length = len(utr_config['utr5_template']) + len(amino_acid_sequence) * 3 + len(utr_config['utr3_template'])
    actual_length = logits.shape[0]

    if actual_length != expected_length:
        raise ValidationError(f"Length mismatch: expected {expected_length}, got {actual_length}")

    return True



def with_error_handling(validator_func: Optional[callable] = None,
                       recovery_manager: Optional[ErrorRecoveryManager] = None):
    """
    Error handling decorator with validation and recovery

    Args:
        validator_func: Optional validation function
        recovery_manager: Optional error recovery manager
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                # Run validation if provided
                if validator_func is not None:
                    validator_func(*args, **kwargs)

                # Execute with recovery if manager provided
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
# Error Handling Decorators for ID3 Components
# =============================================================================

def robust_cai_optimization(max_gpu_retries: int = 2):
    """
    Decorator for robust CAI optimization with automatic GPU memory management

    Args:
        max_gpu_retries: Maximum number of GPU retries

    Usage:
        @robust_cai_optimization(max_gpu_retries=2)
        def optimize_cai(...): ...
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # Save original kwargs for retry
            original_kwargs = kwargs.copy()

            for attempt in range(max_gpu_retries + 1):
                try:
                    return func(*args, **kwargs)

                except torch.cuda.OutOfMemoryError as e:
                    if attempt < max_gpu_retries:
                        logger.warning(f"GPU out of memory, clearing cache and retrying ({attempt + 1}/{max_gpu_retries}): {e}")
                        if torch.cuda.is_available():
                            torch.cuda.empty_cache()
                        # Reset kwargs for retry
                        kwargs = original_kwargs.copy()
                    else:
                        logger.error(f"GPU out of memory after {max_gpu_retries} retries")
                        raise CAIOptimizationError(f"GPU out of memory: {e}")

                except Exception as e:
                    # Other errors are not retried
                    logger.error(f"CAI optimization failed: {e}")
                    raise CAIOptimizationError(f"CAI optimization failed: {e}")

            # This should never be reached
            raise CAIOptimizationError("Unexpected execution path")

        return wrapper
    return decorator


def constraint_error_handler(allow_cpu_fallback: bool = False):
    """
    Decorator for constraint calculation error handling

    Args:
        allow_cpu_fallback: Whether to allow fallback to CPU on GPU OOM

    Usage:
        @constraint_error_handler(allow_cpu_fallback=True)
        def calculate_constraint(...): ...
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)

            except torch.cuda.OutOfMemoryError as e:
                if allow_cpu_fallback:
                    logger.warning(f"GPU out of memory, switching to CPU: {e}")
                    # Try to switch to CPU
                    if 'device' in kwargs:
                        kwargs['device'] = 'cpu'
                    elif len(args) > 0 and hasattr(args[0], 'device_str'):
                        args[0].device_str = 'cpu'
                        if hasattr(args[0], 'to'):
                            args[0].to('cpu')

                    torch.cuda.empty_cache()
                    return func(*args, **kwargs)
                else:
                    logger.error(f"GPU out of memory in constraint calculation: {e}")
                    raise ConstraintError(f"GPU out of memory, please reduce batch size or use CPU: {e}")

            except Exception as e:
                logger.error(f"Constraint calculation failed: {e}")
                raise ConstraintError(f"Constraint calculation failed: {e}")

        return wrapper
    return decorator


@contextmanager
def safe_tensor_operations(fallback_device: str = 'cpu'):
    """
    Context manager for safe tensor operations with automatic fallback

    Args:
        fallback_device: Fallback device on error (default: 'cpu')
    """
    try:
        yield
    except torch.cuda.OutOfMemoryError:
        logger.warning(f"CUDA out of memory, switch to {fallback_device}")
        torch.cuda.empty_cache()
        # Re-raise to allow calling code to handle
        raise
    except RuntimeError as e:
        if "CUDA" in str(e):
            logger.warning(f"CUDA runtime error: {e}, switch to {fallback_device}")
            torch.cuda.empty_cache()
            raise
        else:
            raise


def safe_experiment_runner(experiment_func: Callable,
                         config: Dict[str, Any],
                         fail_fast: bool = True) -> Dict[str, Any]:
    """
    Safely run experiment with comprehensive error handling

    Args:
        experiment_func: Experiment function to run
        config: Experiment configuration
        fail_fast: Whether to raise exception on first error

    Returns:
        Dictionary with success, failed, and error results
    """
    results = {
        'success': [],
        'failed': [],
        'errors': []
    }

    # Use readonly config to prevent modification
    readonly_config = config.copy()

    try:
        result = experiment_func(readonly_config)
        results['success'].append(result)

    except CAIOptimizationError as e:
        error_msg = f"CAI optimization error: {e}"
        logger.error(error_msg)
        results['errors'].append(error_msg)
        results['failed'].append(config.get('variant', 'unknown'))

        if fail_fast:
            raise

    except ConstraintError as e:
        error_msg = f"Constraint calculation error: {e}"
        logger.error(error_msg)
        results['errors'].append(error_msg)
        results['failed'].append(config.get('variant', 'unknown'))

        if fail_fast:
            raise

    except Exception as e:
        error_msg = f"Experiment execution error: {e}"
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
    """Check system resource status"""
    resources = {}

    # Check GPU availability and memory
    if torch.cuda.is_available():
        gpu_memory = torch.cuda.get_device_properties(0).total_memory
        gpu_memory_used = torch.cuda.memory_allocated(0)
        resources['gpu_available'] = True
        resources['gpu_memory_free_pct'] = 1.0 - (gpu_memory_used / gpu_memory)
    else:
        resources['gpu_available'] = False
        resources['gpu_memory_free_pct'] = 0.0

    # CPU is always available
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