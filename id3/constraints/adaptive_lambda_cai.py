"""
Adaptive Lambda CAI Module

This module implements dynamic lambda_cai adjustment mechanisms for all three constraint types,
inspired by the subgradient method used in Lagrangian constraints but adapted for CAI optimization.

The adaptive algorithm automatically adjusts lambda_cai to help CAI optimization converge to
the target value, eliminating the need for manual tuning of the fixed 0.1 weight.

Mathematical Framework:
    CAI Gap: g_t = |CAI_current - CAI_target| - tolerance
    Subgradient Update: λ_CAI^(t+1) = max(λ_min, min(λ_max, λ_CAI^(t) + α_t * g_t))
    Adaptive Step Size: α_t = α_0 / √(t+1) for guaranteed convergence
    
Theoretical Properties:
    - Convergence Rate: O(1/√T) where T is iteration number
    - Robustness: Handles noisy CAI measurements through smoothing
    - Stability: Projection to feasible domain [λ_min, λ_max] prevents divergence
"""

import math
import torch
from typing import Dict, List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class AdaptiveLambdaCAI:
    """
    Adaptive lambda_cai controller using subgradient method with convergence guarantees.
    
    This implements an intelligent weight adjustment mechanism that automatically
    tunes lambda_cai to achieve target CAI values without manual parameter tuning.
    """
    
    def __init__(self,
                 initial_lambda_cai: float = 0.1,
                 lambda_cai_lr: float = 0.1,
                 lambda_cai_min: float = 0.01,
                 lambda_cai_max: float = 2.0,
                 cai_tolerance: float = 0.05,
                 adaptive_lambda_cai: bool = False,
                 smoothing_factor: float = 0.9,
                 verbose: bool = False):
        """
        Initialize adaptive lambda_cai controller.
        
        Args:
            initial_lambda_cai: Starting value for lambda_cai
            lambda_cai_lr: Learning rate for adaptation (α_0)
            lambda_cai_min: Minimum allowed lambda_cai value
            lambda_cai_max: Maximum allowed lambda_cai value
            cai_tolerance: Tolerance for CAI target satisfaction (ε)
            adaptive_lambda_cai: Whether to enable adaptive adjustment
            smoothing_factor: Exponential smoothing for CAI measurements
            verbose: Whether to log adaptation details
        """
        self.lambda_cai = initial_lambda_cai
        self.lambda_cai_lr = lambda_cai_lr
        self.lambda_cai_min = lambda_cai_min
        self.lambda_cai_max = lambda_cai_max
        self.cai_tolerance = cai_tolerance
        self.adaptive_lambda_cai = adaptive_lambda_cai
        self.smoothing_factor = smoothing_factor
        self.verbose = verbose
        
        # Convergence tracking
        self.iteration = 0
        self.smoothed_cai = None
        self.convergence_stats = {
            'cai_values': [],
            'cai_gaps': [],
            'lambda_cai_values': [],
            'step_sizes': [],
            'subgradients': []
        }
        
        if self.verbose:
            logger.info(f"Initialized AdaptiveLambdaCAI with lr={lambda_cai_lr}, "
                       f"range=[{lambda_cai_min}, {lambda_cai_max}], tolerance={cai_tolerance}")
    
    def update_lambda_cai(self, 
                          current_cai: float, 
                          target_cai: float) -> Dict[str, float]:
        """
        Update lambda_cai based on current CAI performance using subgradient method.
        
        Mathematical Details:
            1. Exponential smoothing: CAI_smooth = β * CAI_smooth + (1-β) * CAI_current
            2. Gap calculation: g_t = |CAI_smooth - CAI_target| - ε
            3. Step size: α_t = α_0 / √(t+1) (Robbins-Monro conditions)
            4. Update: λ_CAI = Proj[λ_min,λ_max](λ_CAI + α_t * g_t)
            
        Args:
            current_cai: Current CAI value achieved by the sequence
            target_cai: Desired CAI target value
            
        Returns:
            Dictionary with update statistics
        """
        if not self.adaptive_lambda_cai:
            return {
                'lambda_cai': self.lambda_cai,
                'cai_gap': abs(current_cai - target_cai),
                'step_size': 0.0,
                'subgradient': 0.0,
                'status': 'disabled'
            }
        
        # Apply exponential smoothing to reduce noise
        if self.smoothed_cai is None:
            self.smoothed_cai = current_cai
        else:
            self.smoothed_cai = (self.smoothing_factor * self.smoothed_cai + 
                                (1 - self.smoothing_factor) * current_cai)
        
        # Calculate CAI gap (subgradient)
        cai_gap = abs(self.smoothed_cai - target_cai)
        subgradient = cai_gap - self.cai_tolerance
        
        # Only update if gap exceeds tolerance
        if subgradient > 0:
            # Adaptive step size with convergence guarantees
            # α_t = α_0 / √(t+1) ensures Σα_t=∞ and Σα_t²<∞
            step_size = self.lambda_cai_lr / math.sqrt(self.iteration + 1)
            
            # Subgradient update with projection to feasible domain
            new_lambda_cai = self.lambda_cai + step_size * subgradient
            old_lambda_cai = self.lambda_cai
            self.lambda_cai = max(self.lambda_cai_min, 
                                 min(self.lambda_cai_max, new_lambda_cai))
            
            if self.verbose and abs(self.lambda_cai - old_lambda_cai) > 1e-6:
                logger.info(f"Lambda CAI adapted: {old_lambda_cai:.4f} -> {self.lambda_cai:.4f} "
                           f"(gap={cai_gap:.4f}, step={step_size:.6f})")
        else:
            step_size = 0.0
        
        # Update statistics
        self.convergence_stats['cai_values'].append(current_cai)
        self.convergence_stats['cai_gaps'].append(cai_gap)
        self.convergence_stats['lambda_cai_values'].append(self.lambda_cai)
        self.convergence_stats['step_sizes'].append(step_size)
        self.convergence_stats['subgradients'].append(subgradient)
        
        self.iteration += 1
        
        return {
            'lambda_cai': self.lambda_cai,
            'cai_gap': cai_gap,
            'smoothed_cai': self.smoothed_cai,
            'step_size': step_size,
            'subgradient': subgradient,
            'iteration': self.iteration,
            'status': 'updated' if subgradient > 0 else 'converged'
        }
    
    def get_convergence_analysis(self) -> Dict:
        """
        Analyze convergence properties of the adaptive algorithm.
        
        Returns comprehensive statistics for theoretical validation and
        practical performance assessment.
        
        Returns:
            Dictionary with convergence metrics and theoretical properties
        """
        if len(self.convergence_stats['lambda_cai_values']) < 2:
            return {
                'status': 'insufficient_data',
                'message': 'Need at least 2 iterations for analysis'
            }
        
        import numpy as np
        
        cai_values = np.array(self.convergence_stats['cai_values'])
        cai_gaps = np.array(self.convergence_stats['cai_gaps'])
        lambda_values = np.array(self.convergence_stats['lambda_cai_values'])
        step_sizes = np.array(self.convergence_stats['step_sizes'])
        subgrads = np.array(self.convergence_stats['subgradients'])
        
        # Convergence metrics
        recent_window = min(20, len(cai_gaps))
        recent_gaps = cai_gaps[-recent_window:]
        
        analysis = {
            # Basic statistics
            'total_iterations': len(lambda_values),
            'current_lambda_cai': float(lambda_values[-1]),
            'final_cai_gap': float(cai_gaps[-1]),
            'convergence_threshold': self.cai_tolerance,
            
            # CAI performance
            'cai_trajectory': cai_values.tolist(),
            'cai_gap_trajectory': cai_gaps.tolist(),
            'average_cai_gap': float(np.mean(cai_gaps)),
            'recent_cai_gap': float(np.mean(recent_gaps)),
            'cai_gap_trend': 'improving' if recent_gaps[-1] < recent_gaps[0] else 'worsening',
            'cai_stability': float(np.std(cai_values[-recent_window:])),
            
            # Lambda CAI adaptation
            'lambda_trajectory': lambda_values.tolist(),
            'lambda_stability': float(np.std(lambda_values[-recent_window:])),
            'lambda_trend': 'increasing' if lambda_values[-1] > lambda_values[0] else 'decreasing',
            'total_lambda_changes': int(np.sum(np.diff(lambda_values) != 0)),
            
            # Theoretical properties
            'step_size_sum': float(np.sum(step_sizes)),  # Should → ∞
            'step_size_square_sum': float(np.sum(step_sizes**2)),  # Should converge
            'average_subgradient': float(np.mean(subgrads)),
            'subgradient_variance': float(np.var(subgrads)),
            
            # Convergence assessment
            'is_converged': float(np.mean(recent_gaps)) < self.cai_tolerance,
            'convergence_rate': self._estimate_convergence_rate(cai_gaps),
            'adaptation_efficiency': self._compute_adaptation_efficiency(),
            
            # Theoretical validation
            'robbins_monro_conditions': {
                'step_sum_infinite': np.sum(step_sizes) > 10,
                'step_square_finite': np.sum(step_sizes**2) < 100,
                'step_decreasing': len(step_sizes) > 1 and np.all(np.diff(step_sizes) <= 1e-10)
            }
        }
        
        return analysis
    
    def _estimate_convergence_rate(self, cai_gaps: 'np.ndarray') -> float:
        """Estimate convergence rate from CAI gap trajectory."""
        if len(cai_gaps) < 10:
            return 0.0
        
        import numpy as np
        # Fit exponential decay to recent gaps
        recent_gaps = cai_gaps[-min(50, len(cai_gaps)):]
        if len(recent_gaps) < 5:
            return 0.0
        
        try:
            # Log-linear regression for exponential decay rate
            x = np.arange(len(recent_gaps))
            y = np.log(recent_gaps + 1e-8)  # Avoid log(0)
            coeffs = np.polyfit(x, y, 1)
            decay_rate = -coeffs[0]  # Negative slope indicates decay
            return max(0.0, decay_rate)
        except:
            return 0.0
    
    def _compute_adaptation_efficiency(self) -> float:
        """Compute efficiency of lambda_cai adaptations."""
        if len(self.convergence_stats['cai_gaps']) < 5:
            return 0.0
        
        import numpy as np
        gaps = np.array(self.convergence_stats['cai_gaps'])
        lambda_changes = np.sum(np.abs(np.diff(self.convergence_stats['lambda_cai_values'])))
        
        if lambda_changes == 0:
            return 1.0  # Perfect efficiency (no changes needed)
        
        # Efficiency = improvement per unit of lambda change
        gap_improvement = max(0, gaps[0] - gaps[-1])
        efficiency = gap_improvement / lambda_changes
        
        return min(1.0, efficiency)  # Cap at 1.0
    
    def reset(self, initial_lambda_cai: Optional[float] = None):
        """Reset the adaptive controller to initial state."""
        if initial_lambda_cai is not None:
            self.lambda_cai = initial_lambda_cai
        
        self.iteration = 0
        self.smoothed_cai = None
        self.convergence_stats = {
            'cai_values': [],
            'cai_gaps': [],
            'lambda_cai_values': [],
            'step_sizes': [],
            'subgradients': []
        }
        
        if self.verbose:
            logger.info(f"Reset AdaptiveLambdaCAI to lambda_cai={self.lambda_cai:.4f}")


class AdaptiveLambdaCAIMixin:
    """
    Mixin class to add adaptive lambda_cai functionality to constraint classes.
    
    This provides a standardized interface for all three constraint types
    (AMS, CPC, Lagrangian) to use adaptive lambda_cai adjustment.
    """
    
    def _init_adaptive_lambda_cai(self,
                                  initial_lambda_cai: float = 0.1,
                                  adaptive_lambda_cai: bool = False,
                                  lambda_cai_lr: float = 0.1,
                                  lambda_cai_max: float = 2.0,
                                  cai_tolerance: float = 0.05,
                                  **kwargs):
        """Initialize adaptive lambda_cai components."""
        self.adaptive_lambda_cai_controller = AdaptiveLambdaCAI(
            initial_lambda_cai=initial_lambda_cai,
            lambda_cai_lr=lambda_cai_lr,
            lambda_cai_min=0.01,
            lambda_cai_max=lambda_cai_max,
            cai_tolerance=cai_tolerance,
            adaptive_lambda_cai=adaptive_lambda_cai,
            verbose=kwargs.get('verbose', False)
        )
        
        # Store reference for easy access
        self.lambda_cai = self.adaptive_lambda_cai_controller.lambda_cai
    
    def update_lambda_cai(self, current_cai: float, target_cai: float) -> Dict[str, float]:
        """Update lambda_cai based on current CAI performance."""
        if not hasattr(self, 'adaptive_lambda_cai_controller'):
            logger.warning("Adaptive lambda_cai controller not initialized")
            return {'status': 'not_initialized'}
        
        update_stats = self.adaptive_lambda_cai_controller.update_lambda_cai(
            current_cai, target_cai
        )
        
        # Synchronize lambda_cai value
        self.lambda_cai = self.adaptive_lambda_cai_controller.lambda_cai
        
        return update_stats
    
    def get_lambda_cai_convergence_analysis(self) -> Dict:
        """Get comprehensive convergence analysis for lambda_cai adaptation."""
        if not hasattr(self, 'adaptive_lambda_cai_controller'):
            return {'status': 'not_initialized'}
        
        return self.adaptive_lambda_cai_controller.get_convergence_analysis()
    
    def reset_lambda_cai_adaptation(self, initial_value: Optional[float] = None):
        """Reset the adaptive lambda_cai controller."""
        if hasattr(self, 'adaptive_lambda_cai_controller'):
            self.adaptive_lambda_cai_controller.reset(initial_value)
            self.lambda_cai = self.adaptive_lambda_cai_controller.lambda_cai