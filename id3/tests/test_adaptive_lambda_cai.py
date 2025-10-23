"""
Integration Tests for Adaptive Lambda CAI Functionality

This module provides comprehensive tests for the dynamic lambda_cai adjustment
mechanism across all three constraint types (AMS, CPC, Lagrangian).

Tests cover:
- Basic functionality and parameter validation
- Convergence behavior and theoretical properties
- Integration with existing constraint mechanisms  
- Backward compatibility and configuration management
- Performance and stability under various conditions
"""

import pytest
import torch
import numpy as np
from typing import Dict, Any, List, Tuple
import logging

# Import constraint classes
from id3.constraints.amino_matching import AminoMatchingSoftmax
from id3.constraints.codon_profile import CodonProfileConstraint
from id3.constraints.lagrangian import LagrangianConstraint

# Import adaptive lambda_cai components
from id3.constraints.adaptive_lambda_cai import AdaptiveLambdaCAI, AdaptiveLambdaCAIMixin
from id3.constraints.config import (
    ConstraintConfig, PresetConfigs, 
    get_backward_compatible_kwargs, load_config_from_dict
)

# Test configuration
TEST_AMINO_ACID_SEQUENCE = "MKTAYGGSVKPGVYH"  # Short test sequence
DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

logger = logging.getLogger(__name__)


class TestAdaptiveLambdaCAICore:
    """Test the core adaptive lambda_cai algorithm."""
    
    def test_initialization(self):
        """Test proper initialization of adaptive controller."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=True,
            lambda_cai_lr=0.1,
            lambda_cai_max=2.0,
            cai_tolerance=0.05
        )
        
        assert controller.lambda_cai == 0.1
        assert controller.adaptive_lambda_cai == True
        assert controller.iteration == 0
        assert controller.smoothed_cai is None
        assert len(controller.convergence_stats['lambda_cai_values']) == 0
    
    def test_disabled_adaptation(self):
        """Test that disabled adaptation doesn't change lambda_cai."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=False
        )
        
        # Multiple updates should not change lambda_cai
        for _ in range(10):
            stats = controller.update_lambda_cai(current_cai=0.6, target_cai=0.8)
            assert stats['status'] == 'disabled'
            assert controller.lambda_cai == 0.1
    
    def test_convergence_behavior(self):
        """Test convergence when CAI is close to target."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=True,
            cai_tolerance=0.05
        )
        
        # CAI close to target should not trigger updates
        stats = controller.update_lambda_cai(current_cai=0.78, target_cai=0.8)
        assert stats['status'] == 'converged'
        assert controller.lambda_cai == 0.1  # No change
    
    def test_adaptation_direction(self):
        """Test that adaptation moves in correct direction."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=True,
            cai_tolerance=0.01
        )
        
        initial_lambda = controller.lambda_cai
        
        # Low CAI should increase lambda_cai
        stats = controller.update_lambda_cai(current_cai=0.5, target_cai=0.8)
        assert stats['status'] == 'updated'
        assert controller.lambda_cai > initial_lambda
        
        # High CAI should decrease lambda_cai (if we had a way to simulate this)
        # Note: The current algorithm only increases lambda_cai when gap > tolerance
        # This is by design for stability
    
    def test_convergence_analysis(self):
        """Test convergence analysis functionality."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=True
        )
        
        # Need at least 2 iterations for analysis
        analysis = controller.get_convergence_analysis()
        assert analysis['status'] == 'insufficient_data'
        
        # Simulate several updates
        cai_values = [0.6, 0.65, 0.7, 0.75, 0.78, 0.8]
        for cai in cai_values:
            controller.update_lambda_cai(cai, target_cai=0.8)
        
        analysis = controller.get_convergence_analysis()
        assert analysis['status'] != 'insufficient_data'
        assert 'total_iterations' in analysis
        assert 'convergence_rate' in analysis
        assert 'robbins_monro_conditions' in analysis
        assert len(analysis['lambda_trajectory']) == len(cai_values)
    
    def test_parameter_bounds(self):
        """Test that lambda_cai stays within specified bounds."""
        controller = AdaptiveLambdaCAI(
            initial_lambda_cai=0.1,
            adaptive_lambda_cai=True,
            lambda_cai_lr=10.0,  # Very high learning rate
            lambda_cai_min=0.01,
            lambda_cai_max=2.0
        )
        
        # Many updates with large gaps
        for _ in range(100):
            controller.update_lambda_cai(current_cai=0.1, target_cai=0.9)
        
        assert controller.lambda_cai <= 2.0  # Should not exceed maximum
        assert controller.lambda_cai >= 0.01  # Should not go below minimum


class TestConstraintIntegration:
    """Test integration with all three constraint types."""
    
    @pytest.fixture(params=['ams', 'cpc', 'lagrangian'])
    def constraint_type(self, request):
        """Parameterized fixture for all constraint types."""
        return request.param
    
    @pytest.fixture
    def constraint_kwargs(self, constraint_type):
        """Generate constraint-specific kwargs."""
        base_kwargs = {
            'amino_acid_sequence': TEST_AMINO_ACID_SEQUENCE,
            'batch_size': 1,
            'enable_cai': True,
            'cai_target': 0.8,
            'cai_weight': 0.1,
            'device': DEVICE,
            'verbose': False,
            'adaptive_lambda_cai': True,
            'lambda_cai_lr': 0.1,
            'lambda_cai_max': 2.0,
            'cai_tolerance': 0.05
        }\n        \n        if constraint_type == 'lagrangian':\n            base_kwargs.update({\n                'device': torch.device(DEVICE),\n                'initial_lambda': 0.01,\n                'adaptive_lambda': True,\n                'lambda_lr': 0.01,\n                'lambda_max': 10.0\n            })\n        \n        return base_kwargs\n    \n    def get_constraint_class(self, constraint_type: str):\n        \"\"\"Get the constraint class for given type.\"\"\"\n        classes = {\n            'ams': AminoMatchingSoftmax,\n            'cpc': CodonProfileConstraint,\n            'lagrangian': LagrangianConstraint\n        }\n        return classes[constraint_type]\n    \n    def test_constraint_initialization(self, constraint_type, constraint_kwargs):\n        \"\"\"Test that all constraints initialize properly with adaptive lambda_cai.\"\"\"\n        ConstraintClass = self.get_constraint_class(constraint_type)\n        constraint = ConstraintClass(**constraint_kwargs)\n        \n        # Check basic properties\n        assert constraint.enable_cai == True\n        assert constraint.cai_target == 0.8\n        assert hasattr(constraint, 'lambda_cai')\n        assert constraint.lambda_cai == 0.1  # Initial value\n        \n        # Check adaptive controller\n        assert hasattr(constraint, 'adaptive_lambda_cai_controller')\n        assert constraint.adaptive_lambda_cai_controller.adaptive_lambda_cai == True\n        \n        # Check CAI loss module\n        assert hasattr(constraint, 'cai_loss_module')\n        assert constraint.cai_loss_module is not None\n    \n    def test_forward_pass_basic(self, constraint_type, constraint_kwargs):\n        \"\"\"Test basic forward pass with adaptive lambda_cai enabled.\"\"\"\n        ConstraintClass = self.get_constraint_class(constraint_type)\n        constraint = ConstraintClass(**constraint_kwargs)\n        \n        # Forward pass should work without errors\n        result = constraint.forward(alpha=0.0, beta=0.0, tau=1.0)\n        \n        # Check basic result structure\n        required_keys = ['rna_sequence']\n        if constraint_type == 'lagrangian':\n            required_keys.extend(['probabilities', 'discrete_indices'])\n        elif constraint_type in ['ams', 'cpc']:\n            required_keys.extend(['codon_probs', 'selected_codons'])\n        \n        for key in required_keys:\n            assert key in result, f\"Missing key {key} in {constraint_type} result\"\n    \n    def test_loss_computation_with_adaptation(self, constraint_type, constraint_kwargs):\n        \"\"\"Test loss computation with lambda_cai adaptation.\"\"\"\n        ConstraintClass = self.get_constraint_class(constraint_type)\n        constraint = ConstraintClass(**constraint_kwargs)\n        \n        # Get forward pass result\n        result = constraint.forward(alpha=0.0, beta=0.0, tau=1.0)\n        \n        # Simulate accessibility loss\n        accessibility_loss = torch.tensor(1.5, requires_grad=True)\n        \n        # Compute total loss\n        if constraint_type == 'lagrangian':\n            # Lagrangian needs constraint penalty\n            constraint_penalty = torch.tensor(0.1)\n            probabilities = result['probabilities']\n            enhanced_sequence = result.get('enhanced_sequence')\n            \n            loss_result = constraint.compute_total_loss(\n                accessibility_loss, constraint_penalty, probabilities, enhanced_sequence\n            )\n        else:\n            # AMS and CPC\n            codon_probs = result['codon_probs']\n            enhanced_codon_dist = result.get('enhanced_codon_dist')\n            \n            loss_result = constraint.compute_total_loss(\n                accessibility_loss, codon_probs, enhanced_codon_dist\n            )\n        \n        # Check loss result structure\n        assert 'total_loss' in loss_result\n        assert 'accessibility_loss' in loss_result\n        assert 'cai_loss' in loss_result\n        assert 'lambda_cai' in loss_result\n        assert 'eval_cai' in loss_result\n        \n        # Check adaptive statistics if available\n        if 'adaptive_lambda_stats' in loss_result:\n            stats = loss_result['adaptive_lambda_stats']\n            assert 'lambda_cai' in stats\n            assert 'status' in stats\n    \n    def test_backward_compatibility(self, constraint_type):\n        \"\"\"Test that constraints work with legacy parameters.\"\"\"\n        # Legacy kwargs without adaptive features\n        legacy_kwargs = {\n            'amino_acid_sequence': TEST_AMINO_ACID_SEQUENCE,\n            'batch_size': 1,\n            'enable_cai': True,\n            'cai_target': 0.8,\n            'cai_weight': 0.1,\n            'device': DEVICE if constraint_type != 'lagrangian' else torch.device(DEVICE),\n            'verbose': False\n            # No adaptive_lambda_cai parameter - should default to False\n        }\n        \n        if constraint_type == 'lagrangian':\n            legacy_kwargs.update({\n                'initial_lambda': 0.01,\n                'adaptive_lambda': True\n            })\n        \n        ConstraintClass = self.get_constraint_class(constraint_type)\n        constraint = ConstraintClass(**legacy_kwargs)\n        \n        # Should initialize without adaptive controller\n        if constraint.enable_cai:\n            # CAI enabled but lambda_cai should be fixed\n            assert constraint.lambda_cai == 0.1\n            # May or may not have adaptive controller depending on defaults\n        else:\n            assert constraint.lambda_cai == 0.1\n        \n        # Forward pass should still work\n        result = constraint.forward(alpha=0.0, beta=0.0, tau=1.0)\n        assert 'rna_sequence' in result or 'codon_probs' in result\n\n\nclass TestConfigurationSystem:\n    \"\"\"Test the configuration management system.\"\"\"\n    \n    def test_preset_configs(self):\n        \"\"\"Test predefined configuration presets.\"\"\"\n        # Test default fixed config\n        config = PresetConfigs.default_fixed()\n        assert config.adaptive_config.adaptive_lambda_cai == False\n        \n        # Test conservative adaptive config\n        config = PresetConfigs.adaptive_conservative()\n        assert config.adaptive_config.adaptive_lambda_cai == True\n        assert config.adaptive_config.lambda_cai_lr == 0.05\n        assert config.adaptive_config.smoothing_factor == 0.95\n        \n        # Test aggressive adaptive config\n        config = PresetConfigs.adaptive_aggressive()\n        assert config.adaptive_config.lambda_cai_lr == 0.2\n        assert config.adaptive_config.lambda_cai_max == 3.0\n        \n        # Test research mode\n        config = PresetConfigs.research_mode()\n        assert config.experimental_features['enhanced_convergence_tracking'] == True\n    \n    def test_config_to_kwargs(self):\n        \"\"\"Test conversion from config to constraint kwargs.\"\"\"\n        config = ConstraintConfig(\n            enable_cai=True,\n            cai_target=0.9,\n            cai_weight=0.15\n        )\n        config.adaptive_config.adaptive_lambda_cai = True\n        config.adaptive_config.lambda_cai_lr = 0.12\n        \n        # Test for each constraint type\n        for constraint_type in ['ams', 'cpc', 'lagrangian']:\n            kwargs = config.get_constraint_kwargs(constraint_type)\n            \n            assert kwargs['enable_cai'] == True\n            assert kwargs['cai_target'] == 0.9\n            assert kwargs['cai_weight'] == 0.15\n            assert kwargs['adaptive_lambda_cai'] == True\n            assert kwargs['lambda_cai_lr'] == 0.12\n    \n    def test_backward_compatibility_conversion(self):\n        \"\"\"Test conversion of legacy parameter names.\"\"\"\n        legacy_kwargs = {\n            'cai_lambda': 0.2,  # Old name\n            'adaptive_cai': True,  # Old name\n            'cai_lr': 0.15,  # Old name\n            'max_cai_weight': 3.0,  # Old name\n            'some_other_param': 'value'  # Should be preserved\n        }\n        \n        updated_kwargs = get_backward_compatible_kwargs(**legacy_kwargs)\n        \n        assert updated_kwargs['cai_weight'] == 0.2\n        assert updated_kwargs['adaptive_lambda_cai'] == True\n        assert updated_kwargs['lambda_cai_lr'] == 0.15\n        assert updated_kwargs['lambda_cai_max'] == 3.0\n        assert updated_kwargs['some_other_param'] == 'value'\n        \n        # Original keys should be removed\n        assert 'cai_lambda' not in updated_kwargs\n        assert 'adaptive_cai' not in updated_kwargs\n    \n    def test_config_from_dict(self):\n        \"\"\"Test loading configuration from dictionary.\"\"\"\n        config_dict = {\n            'enable_cai': True,\n            'cai_target': 0.85,\n            'cai_weight': 0.12,\n            'adaptive_config': {\n                'adaptive_lambda_cai': True,\n                'lambda_cai_lr': 0.08,\n                'lambda_cai_max': 2.5\n            },\n            'constraint_specific': {\n                'ams': {'device': 'cpu'},\n                'cpc': {'batch_size': 2}\n            }\n        }\n        \n        config = load_config_from_dict(config_dict)\n        \n        assert config.enable_cai == True\n        assert config.cai_target == 0.85\n        assert config.adaptive_config.adaptive_lambda_cai == True\n        assert config.adaptive_config.lambda_cai_lr == 0.08\n        assert config.constraint_specific['ams']['device'] == 'cpu'\n\n\nclass TestPerformanceAndStability:\n    \"\"\"Test performance characteristics and stability.\"\"\"\n    \n    def test_adaptation_stability(self):\n        \"\"\"Test that adaptation remains stable under various conditions.\"\"\"\n        controller = AdaptiveLambdaCAI(\n            initial_lambda_cai=0.1,\n            adaptive_lambda_cai=True,\n            lambda_cai_lr=0.1\n        )\n        \n        # Simulate noisy CAI measurements\n        np.random.seed(42)\n        target_cai = 0.8\n        \n        lambda_history = []\n        for i in range(100):\n            # Add noise to CAI measurements\n            noise = np.random.normal(0, 0.05)\n            noisy_cai = 0.7 + 0.1 * (i / 100) + noise  # Gradually improving CAI\n            \n            stats = controller.update_lambda_cai(noisy_cai, target_cai)\n            lambda_history.append(controller.lambda_cai)\n        \n        # Check stability: lambda_cai shouldn't oscillate wildly\n        lambda_array = np.array(lambda_history)\n        \n        # Standard deviation of later values should be reasonable\n        later_values = lambda_array[-20:]\n        stability = np.std(later_values)\n        \n        assert stability < 0.5, f\"Lambda CAI too unstable: std={stability}\"\n    \n    def test_convergence_properties(self):\n        \"\"\"Test theoretical convergence properties.\"\"\"\n        controller = AdaptiveLambdaCAI(\n            initial_lambda_cai=0.1,\n            adaptive_lambda_cai=True,\n            lambda_cai_lr=0.1\n        )\n        \n        # Simulate convergence scenario\n        for i in range(50):\n            # CAI gradually improves towards target\n            current_cai = 0.6 + 0.2 * (1 - np.exp(-i / 10))\n            controller.update_lambda_cai(current_cai, target_cai=0.8)\n        \n        analysis = controller.get_convergence_analysis()\n        \n        # Check Robbins-Monro conditions\n        conditions = analysis['robbins_monro_conditions']\n        \n        # Step sizes should be decreasing\n        assert conditions['step_decreasing'] == True\n        \n        # Should show convergence trend\n        if analysis['final_cai_gap'] < analysis['average_cai_gap']:\n            assert analysis['cai_gap_trend'] == 'improving'\n    \n    def test_computational_efficiency(self):\n        \"\"\"Test that adaptive mechanism doesn't significantly slow down training.\"\"\"\n        import time\n        \n        # Test with AMS constraint (similar for others)\n        constraint_fixed = AminoMatchingSoftmax(\n            amino_acid_sequence=TEST_AMINO_ACID_SEQUENCE,\n            enable_cai=True,\n            adaptive_lambda_cai=False,  # Fixed\n            device=DEVICE\n        )\n        \n        constraint_adaptive = AminoMatchingSoftmax(\n            amino_acid_sequence=TEST_AMINO_ACID_SEQUENCE,\n            enable_cai=True,\n            adaptive_lambda_cai=True,  # Adaptive\n            device=DEVICE\n        )\n        \n        # Time forward passes\n        def time_forward_passes(constraint, num_passes=50):\n            start_time = time.time()\n            for _ in range(num_passes):\n                result = constraint.forward(alpha=0.0, beta=0.0, tau=1.0)\n                \n                # Simulate loss computation\n                accessibility_loss = torch.tensor(1.0)\n                codon_probs = result['codon_probs']\n                constraint.compute_total_loss(accessibility_loss, codon_probs)\n            \n            return time.time() - start_time\n        \n        time_fixed = time_forward_passes(constraint_fixed)\n        time_adaptive = time_forward_passes(constraint_adaptive)\n        \n        # Adaptive should not be more than 50% slower\n        overhead = (time_adaptive - time_fixed) / time_fixed\n        assert overhead < 0.5, f\"Too much overhead: {overhead:.2%}\"\n\n\nclass TestEdgeCases:\n    \"\"\"Test edge cases and error handling.\"\"\"\n    \n    def test_extreme_cai_values(self):\n        \"\"\"Test behavior with extreme CAI values.\"\"\"\n        controller = AdaptiveLambdaCAI(\n            initial_lambda_cai=0.1,\n            adaptive_lambda_cai=True\n        )\n        \n        # Test with very low CAI\n        stats = controller.update_lambda_cai(current_cai=0.01, target_cai=0.8)\n        assert stats['status'] == 'updated'\n        assert controller.lambda_cai > 0.1\n        \n        # Test with CAI > 1.0 (shouldn't happen but should be handled)\n        stats = controller.update_lambda_cai(current_cai=1.1, target_cai=0.8)\n        # Should still work without crashing\n        assert 'status' in stats\n    \n    def test_invalid_parameters(self):\n        \"\"\"Test handling of invalid parameters.\"\"\"\n        # Test invalid constraint type in config\n        config = ConstraintConfig()\n        \n        with pytest.raises(ValueError):\n            config.get_constraint_kwargs('invalid_type')\n        \n        # Test invalid lambda_cai bounds\n        with pytest.raises((ValueError, AssertionError)):\n            AdaptiveLambdaCAI(\n                lambda_cai_min=2.0,  # Min > max\n                lambda_cai_max=1.0\n            )\n    \n    def test_zero_learning_rate(self):\n        \"\"\"Test behavior with zero learning rate.\"\"\"\n        controller = AdaptiveLambdaCAI(\n            initial_lambda_cai=0.1,\n            adaptive_lambda_cai=True,\n            lambda_cai_lr=0.0  # Zero learning rate\n        )\n        \n        initial_lambda = controller.lambda_cai\n        \n        # Updates should not change lambda_cai\n        for _ in range(10):\n            controller.update_lambda_cai(current_cai=0.5, target_cai=0.8)\n            assert controller.lambda_cai == initial_lambda\n\n\n# Helper functions for manual testing and validation\ndef run_convergence_experiment(constraint_type: str = 'ams', \n                               target_cai: float = 0.8,\n                               num_iterations: int = 100) -> Dict[str, Any]:\n    \"\"\"\n    Run a convergence experiment with adaptive lambda_cai.\n    \n    This function can be used for manual validation and performance analysis.\n    \"\"\"\n    # Create constraint with adaptive lambda_cai\n    kwargs = {\n        'amino_acid_sequence': TEST_AMINO_ACID_SEQUENCE,\n        'enable_cai': True,\n        'cai_target': target_cai,\n        'adaptive_lambda_cai': True,\n        'lambda_cai_lr': 0.1,\n        'verbose': True,\n        'device': DEVICE\n    }\n    \n    if constraint_type == 'lagrangian':\n        kwargs['device'] = torch.device(DEVICE)\n        constraint = LagrangianConstraint(**kwargs)\n    elif constraint_type == 'ams':\n        constraint = AminoMatchingSoftmax(**kwargs)\n    elif constraint_type == 'cpc':\n        constraint = CodonProfileConstraint(**kwargs)\n    else:\n        raise ValueError(f\"Unknown constraint type: {constraint_type}\")\n    \n    # Track metrics\n    lambda_history = []\n    cai_history = []\n    loss_history = []\n    \n    for i in range(num_iterations):\n        # Forward pass\n        result = constraint.forward(alpha=0.0, beta=0.0, tau=1.0)\n        \n        # Compute loss\n        accessibility_loss = torch.tensor(1.0 + 0.5 * np.random.random())\n        \n        if constraint_type == 'lagrangian':\n            constraint_penalty = torch.tensor(0.1)\n            loss_result = constraint.compute_total_loss(\n                accessibility_loss, constraint_penalty, \n                result['probabilities'], result.get('enhanced_sequence')\n            )\n        else:\n            loss_result = constraint.compute_total_loss(\n                accessibility_loss, result['codon_probs'], \n                result.get('enhanced_codon_dist')\n            )\n        \n        # Record metrics\n        lambda_history.append(constraint.lambda_cai)\n        cai_history.append(loss_result.get('eval_cai', 0.0))\n        loss_history.append(loss_result['total_loss'].item())\n    \n    # Get convergence analysis\n    if hasattr(constraint, 'adaptive_lambda_cai_controller'):\n        convergence_analysis = constraint.get_lambda_cai_convergence_analysis()\n    else:\n        convergence_analysis = {}\n    \n    return {\n        'constraint_type': constraint_type,\n        'target_cai': target_cai,\n        'lambda_history': lambda_history,\n        'cai_history': cai_history,\n        'loss_history': loss_history,\n        'convergence_analysis': convergence_analysis,\n        'final_lambda_cai': constraint.lambda_cai,\n        'final_cai': cai_history[-1] if cai_history else 0.0\n    }\n\n\nif __name__ == \"__main__\":\n    # Run basic functionality test when script is executed directly\n    print(\"Running basic adaptive lambda_cai functionality test...\")\n    \n    # Test each constraint type\n    for constraint_type in ['ams', 'cpc', 'lagrangian']:\n        print(f\"\\nTesting {constraint_type.upper()} constraint:\")\n        \n        try:\n            experiment_result = run_convergence_experiment(\n                constraint_type=constraint_type,\n                num_iterations=20\n            )\n            \n            print(f\"  Final lambda_cai: {experiment_result['final_lambda_cai']:.4f}\")\n            print(f\"  Final CAI: {experiment_result['final_cai']:.4f}\")\n            print(f\"  Target CAI: {experiment_result['target_cai']:.4f}\")\n            print(f\"  ✓ Success\")\n            \n        except Exception as e:\n            print(f\"  ✗ Error: {e}\")\n    \n    print(\"\\nBasic test completed.\")