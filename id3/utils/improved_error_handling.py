"""
Improved error handling system for ID3 framework.

Provides consistent error handling, logging, and recovery mechanisms
across the entire framework.
"""

import logging
import traceback
import functools
from typing import Any, Callable, Optional, Type, Union, Dict
from enum import Enum
import torch

logger = logging.getLogger(__name__)


class ErrorSeverity(Enum):
    """Error severity levels."""
    DEBUG = "debug"
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


class ID3Error(Exception):
    """Base exception class for ID3 framework."""

    def __init__(self, message: str, severity: ErrorSeverity = ErrorSeverity.ERROR, **kwargs):
        """
        Initialize ID3 error.

        Args:
            message: Error message
            severity: Error severity level
            **kwargs: Additional context information
        """
        super().__init__(message)
        self.severity = severity
        self.context = kwargs

    def log(self):
        """Log the error with appropriate severity."""
        log_func = getattr(logger, self.severity.value)
        log_func(f"{self.__class__.__name__}: {str(self)}")
        if self.context:
            logger.debug(f"Context: {self.context}")


class ConfigurationError(ID3Error):
    """Configuration-related errors."""
    pass


class DataError(ID3Error):
    """Data processing errors."""
    pass


class ComputationError(ID3Error):
    """Computation and numerical errors."""
    pass


class ConstraintError(ID3Error):
    """Constraint satisfaction errors."""
    pass


class MemoryError(ID3Error):
    """Memory and resource errors."""
    pass


def safe_execute(
    default_return: Any = None,
    exceptions: tuple = (Exception,),
    log_errors: bool = True,
    raise_on_critical: bool = True
):
    """
    Decorator for safe function execution with error handling.

    Args:
        default_return: Default value to return on error
        exceptions: Tuple of exceptions to catch
        log_errors: Whether to log errors
        raise_on_critical: Whether to re-raise critical errors

    Returns:
        Decorator function
    """
    def decorator(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except exceptions as e:
                if log_errors:
                    logger.error(f"Error in {func.__name__}: {str(e)}")
                    logger.debug(traceback.format_exc())

                # Check if it's a critical error
                if isinstance(e, ID3Error) and e.severity == ErrorSeverity.CRITICAL:
                    if raise_on_critical:
                        raise
                    else:
                        logger.critical(f"Critical error suppressed: {str(e)}")

                return default_return

        return wrapper
    return decorator


def validate_tensor(
    tensor: torch.Tensor,
    name: str = "tensor",
    check_nan: bool = True,
    check_inf: bool = True,
    check_shape: Optional[tuple] = None,
    check_range: Optional[tuple] = None
) -> torch.Tensor:
    """
    Validate tensor properties and raise appropriate errors.

    Args:
        tensor: Tensor to validate
        name: Name for error messages
        check_nan: Check for NaN values
        check_inf: Check for Inf values
        check_shape: Expected shape (if provided)
        check_range: Expected value range (min, max)

    Returns:
        The validated tensor

    Raises:
        DataError: If validation fails
    """
    if not isinstance(tensor, torch.Tensor):
        raise DataError(f"{name} must be a torch.Tensor, got {type(tensor)}")

    if check_nan and torch.isnan(tensor).any():
        raise DataError(
            f"{name} contains NaN values",
            severity=ErrorSeverity.ERROR,
            shape=tensor.shape
        )

    if check_inf and torch.isinf(tensor).any():
        raise DataError(
            f"{name} contains Inf values",
            severity=ErrorSeverity.ERROR,
            shape=tensor.shape
        )

    if check_shape is not None and tensor.shape != check_shape:
        raise DataError(
            f"{name} shape mismatch: expected {check_shape}, got {tensor.shape}",
            severity=ErrorSeverity.ERROR
        )

    if check_range is not None:
        min_val, max_val = check_range
        if tensor.min() < min_val or tensor.max() > max_val:
            raise DataError(
                f"{name} values out of range [{min_val}, {max_val}]",
                severity=ErrorSeverity.WARNING,
                actual_range=(tensor.min().item(), tensor.max().item())
            )

    return tensor


def retry_on_error(
    max_retries: int = 3,
    delay: float = 1.0,
    backoff: float = 2.0,
    exceptions: tuple = (Exception,)
):
    """
    Decorator for retrying function execution on errors.

    Args:
        max_retries: Maximum number of retry attempts
        delay: Initial delay between retries (seconds)
        backoff: Delay multiplication factor for each retry
        exceptions: Tuple of exceptions to retry on

    Returns:
        Decorator function
    """
    def decorator(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            import time

            current_delay = delay
            last_exception = None

            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e
                    if attempt < max_retries:
                        logger.warning(
                            f"Attempt {attempt + 1}/{max_retries + 1} failed for {func.__name__}: {str(e)}"
                        )
                        time.sleep(current_delay)
                        current_delay *= backoff
                    else:
                        logger.error(
                            f"All {max_retries + 1} attempts failed for {func.__name__}"
                        )

            # Re-raise the last exception if all retries failed
            if last_exception:
                raise last_exception

        return wrapper
    return decorator


class ErrorContext:
    """Context manager for error handling with cleanup."""

    def __init__(
        self,
        operation_name: str,
        cleanup_func: Optional[Callable] = None,
        suppress_errors: bool = False
    ):
        """
        Initialize error context.

        Args:
            operation_name: Name of the operation for logging
            cleanup_func: Function to call for cleanup on error
            suppress_errors: Whether to suppress exceptions
        """
        self.operation_name = operation_name
        self.cleanup_func = cleanup_func
        self.suppress_errors = suppress_errors

    def __enter__(self):
        """Enter context."""
        logger.debug(f"Starting operation: {self.operation_name}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Exit context with error handling.

        Args:
            exc_type: Exception type
            exc_val: Exception value
            exc_tb: Exception traceback

        Returns:
            True to suppress exception, False to propagate
        """
        if exc_type is not None:
            logger.error(
                f"Error in {self.operation_name}: {exc_type.__name__}: {exc_val}"
            )

            # Perform cleanup if provided
            if self.cleanup_func:
                try:
                    self.cleanup_func()
                    logger.debug(f"Cleanup completed for {self.operation_name}")
                except Exception as cleanup_error:
                    logger.error(f"Cleanup failed: {cleanup_error}")

            # Return True to suppress the exception if configured
            return self.suppress_errors
        else:
            logger.debug(f"Completed operation: {self.operation_name}")
            return False


def handle_numerical_errors(tensor: torch.Tensor, name: str = "tensor") -> torch.Tensor:
    """
    Handle numerical errors in tensors by replacing NaN/Inf values.

    Args:
        tensor: Input tensor
        name: Tensor name for logging

    Returns:
        Cleaned tensor
    """
    if torch.isnan(tensor).any() or torch.isinf(tensor).any():
        logger.warning(f"Numerical issues detected in {name}, applying corrections")

        # Replace NaN with 0
        tensor = torch.nan_to_num(tensor, nan=0.0, posinf=1e10, neginf=-1e10)

    return tensor


class ErrorRecovery:
    """Strategies for recovering from errors."""

    @staticmethod
    def checkpoint_recovery(checkpoint_path: str, operation: Callable):
        """
        Execute operation with checkpoint-based recovery.

        Args:
            checkpoint_path: Path to save/load checkpoint
            operation: Operation to execute

        Returns:
            Operation result
        """
        import pickle
        from pathlib import Path

        checkpoint = Path(checkpoint_path)

        try:
            # Try to load from checkpoint first
            if checkpoint.exists():
                logger.info(f"Loading from checkpoint: {checkpoint}")
                with open(checkpoint, 'rb') as f:
                    return pickle.load(f)
        except Exception as e:
            logger.warning(f"Failed to load checkpoint: {e}")

        # Execute operation
        try:
            result = operation()

            # Save checkpoint
            try:
                with open(checkpoint, 'wb') as f:
                    pickle.dump(result, f)
                logger.info(f"Saved checkpoint: {checkpoint}")
            except Exception as e:
                logger.warning(f"Failed to save checkpoint: {e}")

            return result

        except Exception as e:
            logger.error(f"Operation failed: {e}")
            raise


# Global error statistics tracking
class ErrorStatistics:
    """Track error statistics for monitoring."""

    def __init__(self):
        self.errors: Dict[str, int] = {}
        self.warnings: Dict[str, int] = {}

    def record_error(self, error_type: str):
        """Record an error occurrence."""
        self.errors[error_type] = self.errors.get(error_type, 0) + 1

    def record_warning(self, warning_type: str):
        """Record a warning occurrence."""
        self.warnings[warning_type] = self.warnings.get(warning_type, 0) + 1

    def get_stats(self) -> Dict:
        """Get error statistics."""
        return {
            'errors': dict(self.errors),
            'warnings': dict(self.warnings),
            'total_errors': sum(self.errors.values()),
            'total_warnings': sum(self.warnings.values())
        }

    def reset(self):
        """Reset statistics."""
        self.errors.clear()
        self.warnings.clear()


# Global instance
_error_stats = ErrorStatistics()


def get_error_statistics() -> ErrorStatistics:
    """Get global error statistics instance."""
    return _error_stats