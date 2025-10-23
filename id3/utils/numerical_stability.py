"""


"""

import torch
import torch.nn.functional as F
import math
from typing import Optional, Tuple


class NumericalConfig:

    





    








def safe_log(x: torch.Tensor, eps: float = NumericalConfig.LOG_EPS) -> torch.Tensor:
    """
    数值稳定的对数函数
    
    Args:
        x: 输入张量
        eps: 保护epsilon
        
    Returns:
        稳定的log(x + eps)
    """
    return torch.log(torch.clamp(x, min=eps))


def safe_exp(x: torch.Tensor, max_val: float = NumericalConfig.MAX_EXP) -> torch.Tensor:
    """
    数值稳定的指数函数
    
    Args:
        x: 输入张量
        max_val: 最大值限制
        
    Returns:
        稳定的exp(clamp(x))
    """
    return torch.exp(torch.clamp(x, max=max_val))


def safe_softmax(logits: torch.Tensor, dim: int = -1, temperature: float = 1.0) -> torch.Tensor:
    """
    数值稳定的softmax函数
    
    Args:
        logits: logit张量
        dim: 应用softmax的维度
        temperature: 温度参数
        
    Returns:
        稳定的softmax概率
    """

    if temperature != 1.0:
        logits = logits / temperature
    

    logits = torch.clamp(logits, NumericalConfig.MIN_LOGIT, NumericalConfig.MAX_LOGIT)
    

    max_logits = torch.max(logits, dim=dim, keepdim=True)[0]
    stable_logits = logits - max_logits
    exp_logits = safe_exp(stable_logits)
    

    sum_exp = torch.sum(exp_logits, dim=dim, keepdim=True)
    probs = exp_logits / torch.clamp(sum_exp, min=NumericalConfig.EPS)
    

    probs = torch.clamp(probs, NumericalConfig.MIN_PROB, NumericalConfig.MAX_PROB)
    
    return probs


def safe_gumbel_softmax(logits: torch.Tensor, 
                       temperature: float = 1.0, 
                       dim: int = -1,
                       hard: bool = False) -> torch.Tensor:
    """
    数值稳定的Gumbel-Softmax采样
    
    Args:
        logits: logit张量
        temperature: Gumbel温度
        dim: softmax维度
        hard: 是否使用hard版本
        
    Returns:
        Gumbel-Softmax采样结果
    """

    gumbel_noise = _sample_gumbel_stable(logits.shape, device=logits.device)
    

    gumbel_logits = logits + gumbel_noise
    

    soft_samples = safe_softmax(gumbel_logits, dim=dim, temperature=temperature)
    
    if hard:

        hard_samples = _make_one_hot(soft_samples, dim=dim)

        return (hard_samples - soft_samples).detach() + soft_samples
    else:
        return soft_samples


def _sample_gumbel_stable(shape: torch.Size, 
                         device: torch.device = None,
                         eps: float = NumericalConfig.EPS) -> torch.Tensor:
    """
    数值稳定的Gumbel噪声采样
    使用双重epsilon保护避免log(0)
    """
    U = torch.rand(shape, device=device)

    U = torch.clamp(U, eps, 1.0 - eps)
    return -safe_log(-safe_log(U, eps), eps)


def _make_one_hot(soft_samples: torch.Tensor, dim: int = -1) -> torch.Tensor:
    """将soft采样转换为one-hot编码"""
    _, max_indices = torch.max(soft_samples, dim=dim, keepdim=True)
    hard_samples = torch.zeros_like(soft_samples)
    hard_samples.scatter_(dim, max_indices, 1.0)
    return hard_samples


def safe_geometric_mean(weights: torch.Tensor, eps: float = NumericalConfig.EPS) -> torch.Tensor:
    """

    
    Args:


        
    Returns:

    """

    safe_weights = torch.clamp(weights, min=eps)
    

    log_weights = safe_log(safe_weights, eps)
    mean_log = torch.mean(log_weights)
    

    geometric_mean = safe_exp(mean_log)
    
    return geometric_mean


def safe_division(numerator: torch.Tensor, 
                 denominator: torch.Tensor,
                 eps: float = NumericalConfig.DIV_EPS) -> torch.Tensor:
    """

    
    Args:



        
    Returns:

    """
    safe_denominator = torch.clamp(torch.abs(denominator), min=eps)
    return numerator / safe_denominator


def clip_gradients_by_value(tensor: torch.Tensor, 
                          clip_value: float = 10.0) -> torch.Tensor:
    """

    
    Args:


        
    Returns:

    """
    return torch.clamp(tensor, -clip_value, clip_value)


def check_tensor_health(tensor: torch.Tensor, name: str = "tensor") -> Tuple[bool, str]:
    """

    
    Args:


        
    Returns:

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

    
    Args:


        
    Returns:

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

        tensor = torch.nan_to_num(tensor, nan=0.0, posinf=1e6, neginf=-1e6)
        return torch.clamp(tensor, -1e6, 1e6)


class NumericalStabilityMixin:

    
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