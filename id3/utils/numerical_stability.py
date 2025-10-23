"""
Numerical Stability Module

Provides stable mathematical operations to prevent numerical overflow/underflow.
"""

import torch
import torch.nn.functional as F
import math
from typing import Optional, Tuple


class NumericalConfig:
    """Numerical stability configuration parameters"""

    # Epsilon protection values
    EPS = 1e-8              # General epsilon
    LOG_EPS = 1e-10         # Log operation epsilon
    DIV_EPS = 1e-12         # Division operation epsilon

    # Value range limits
    MIN_LOGIT = -10.0       # Minimum logit value
    MAX_LOGIT = 10.0        # Maximum logit value
    MAX_EXP = 20.0          # Maximum exponent value

    # Probability range
    MIN_PROB = 1e-10        # Minimum probability
    MAX_PROB = 1.0 - 1e-10  # Maximum probability


def safe_log(x: torch.Tensor, eps: float = NumericalConfig.LOG_EPS) -> torch.Tensor:
    """
    Numerically stable logarithm function

    Args:
        x: Input tensor
        eps: Protection epsilon

    Returns:
        Stable log(x + eps)
    """
    return torch.log(torch.clamp(x, min=eps))


def safe_exp(x: torch.Tensor, max_val: float = NumericalConfig.MAX_EXP) -> torch.Tensor:
    """
    Numerically stable exponential function

    Args:
        x: Input tensor
        max_val: Maximum value limit

    Returns:
        Stable exp(clamp(x))
    """
    return torch.exp(torch.clamp(x, max=max_val))


def safe_softmax(logits: torch.Tensor, dim: int = -1, temperature: float = 1.0) -> torch.Tensor:
    """
    Numerically stable softmax function

    Args:
        logits: Logit tensor
        dim: Dimension to apply softmax
        temperature: Temperature parameter

    Returns:
        Stable softmax probabilities
    """
    # Apply temperature scaling
    if temperature != 1.0:
        logits = logits / temperature

    # Clip logits to prevent overflow
    logits = torch.clamp(logits, NumericalConfig.MIN_LOGIT, NumericalConfig.MAX_LOGIT)

    # Subtract maximum for numerical stability
    max_logits = torch.max(logits, dim=dim, keepdim=True)[0]
    stable_logits = logits - max_logits
    exp_logits = safe_exp(stable_logits)

    # Normalize
    sum_exp = torch.sum(exp_logits, dim=dim, keepdim=True)
    probs = exp_logits / torch.clamp(sum_exp, min=NumericalConfig.EPS)

    # Clip probabilities to valid range
    probs = torch.clamp(probs, NumericalConfig.MIN_PROB, NumericalConfig.MAX_PROB)

    return probs


def safe_gumbel_softmax(logits: torch.Tensor,
                       temperature: float = 1.0,
                       dim: int = -1,
                       hard: bool = False) -> torch.Tensor:
    """
    Numerically stable Gumbel-Softmax sampling

    Args:
        logits: Logit tensor
        temperature: Gumbel temperature
        dim: Softmax dimension
        hard: Whether to use hard version

    Returns:
        Gumbel-Softmax sampling result
    """
    # Sample Gumbel noise with numerical stability
    gumbel_noise = _sample_gumbel_stable(logits.shape, device=logits.device)

    # Add Gumbel noise to logits
    gumbel_logits = logits + gumbel_noise

    # Apply stable softmax
    soft_samples = safe_softmax(gumbel_logits, dim=dim, temperature=temperature)

    if hard:
        # Create hard one-hot samples
        hard_samples = _make_one_hot(soft_samples, dim=dim)
        # Use straight-through estimator: forward hard, backward soft
        return (hard_samples - soft_samples).detach() + soft_samples
    else:
        return soft_samples


def _sample_gumbel_stable(shape: torch.Size,
                         device: torch.device = None,
                         eps: float = NumericalConfig.EPS) -> torch.Tensor:
    """
    Numerically stable Gumbel noise sampling
    Uses double epsilon protection to avoid log(0)
    """
    U = torch.rand(shape, device=device)
    # Clip U to avoid log(0) and log(log(0))
    U = torch.clamp(U, eps, 1.0 - eps)
    return -safe_log(-safe_log(U, eps), eps)


def _make_one_hot(soft_samples: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """Convert soft samples to one-hot encoding"""
    _, max_indices = torch.max(soft_samples, dim=dim, keepdim=True)
    hard_samples = torch.zeros_like(soft_samples)
    hard_samples.scatter_(dim, max_indices, 1.0)
    return hard_samples


def safe_geometric_mean(weights: torch.Tensor, eps: float = NumericalConfig.EPS) -> torch.Tensor:
    """
    Numerically stable geometric mean calculation

    Args:
        weights: Weight tensor
        eps: Protection epsilon

    Returns:
        Geometric mean result
    """
    # Clip weights to positive values
    safe_weights = torch.clamp(weights, min=eps)

    # Compute mean in log space
    log_weights = safe_log(safe_weights, eps)
    mean_log = torch.mean(log_weights)

    # Convert back to original space
    geometric_mean = safe_exp(mean_log)

    return geometric_mean


def safe_division(numerator: torch.Tensor,
                 denominator: torch.Tensor,
                 eps: float = NumericalConfig.DIV_EPS) -> torch.Tensor:
    """
    Numerically stable division

    Args:
        numerator: Numerator tensor
        denominator: Denominator tensor
        eps: Protection epsilon

    Returns:
        Division result
    """
    safe_denominator = torch.clamp(torch.abs(denominator), min=eps)
    return numerator / safe_denominator


def clip_gradients_by_value(tensor: torch.Tensor,
                          clip_value: float = 10.0) -> torch.Tensor:
    """
    Clip gradients by value

    Args:
        tensor: Tensor to clip
        clip_value: Clipping value

    Returns:
        Clipped tensor
    """
    return torch.clamp(tensor, -clip_value, clip_value)


def check_tensor_health(tensor: torch.Tensor, name: str = "tensor") -> Tuple[bool, str]:
    """
    Check tensor for NaN, Inf, and extreme values

    Args:
        tensor: Tensor to check
        name: Tensor name for error messages

    Returns:
        Tuple of (is_healthy, message)
    """
    if tensor.numel() == 0:
        return False, f"{name} is empty"

    has_nan = torch.isnan(tensor).any().item()
    has_inf = torch.isinf(tensor).any().item()

    if has_nan:
        return False, f"{name} contains NaN values"

    if has_inf:
        return False, f"{name} contains Inf values"

    tensor_min = tensor.min().item()
    tensor_max = tensor.max().item()

    if abs(tensor_max) > 1e6:
        return False, f"{name} has extremely large values: max={tensor_max}"

    if abs(tensor_min) < 1e-10 and tensor_min != 0:
        return False, f"{name} has extremely small values: min={tensor_min}"

    return True, f"{name} is numerically healthy"


def apply_numerical_protection(tensor: torch.Tensor,
                             operation: str = "general") -> torch.Tensor:
    """
    Apply numerical protection based on operation type

    Args:
        tensor: Input tensor
        operation: Operation type (softmax, log, exp, prob, general)

    Returns:
        Protected tensor
    """
    if operation == "softmax":
        return torch.clamp(tensor, NumericalConfig.MIN_LOGIT, NumericalConfig.MAX_LOGIT)
    elif operation == "log":
        return torch.clamp(tensor, min=NumericalConfig.LOG_EPS)
    elif operation == "exp":
        return torch.clamp(tensor, max=NumericalConfig.MAX_EXP)
    elif operation == "prob":
        return torch.clamp(tensor, NumericalConfig.MIN_PROB, NumericalConfig.MAX_PROB)
    else:  # general
        # Replace NaN and Inf values
        tensor = torch.nan_to_num(tensor, nan=0.0, posinf=1e6, neginf=-1e6)
        return torch.clamp(tensor, -1e6, 1e6)


class NumericalStabilityMixin:
    """Mixin class providing numerical stability methods"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.numerical_config = NumericalConfig()

    def safe_log(self, x: torch.Tensor) -> torch.Tensor:
        return safe_log(x, self.numerical_config.LOG_EPS)

    def safe_exp(self, x: torch.Tensor) -> torch.Tensor:
        return safe_exp(x, self.numerical_config.MAX_EXP)

    def safe_softmax(self, logits: torch.Tensor, dim: int = -1) -> torch.Tensor:
        return safe_softmax(logits, dim)

    def check_tensor_health(self, tensor: torch.Tensor, name: str = "tensor") -> bool:
        is_healthy, message = check_tensor_health(tensor, name)
        if not is_healthy:
            print(f"⚠️ Numerical Warning: {message}")
        return is_healthy
