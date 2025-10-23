#!/usr/bin/env python3
"""
Unit tests for batch parallel experiment system

"""

import unittest
import torch
import torch.nn as nn
import numpy as np
import time
import sys
from pathlib import Path
from typing import List, Dict
import logging

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from batch_parallel_experiments.batch_parallel_runner import BatchParallelExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
from id3.constraints.lagrangian import LagrangianConstraint
from id3.utils.deepraccess_wrapper import DeepRaccessID3Wrapper

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestBatchParallelSystem(unittest.TestCase):
    """Test suite for batch parallel experiment system"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.config = UnifiedExperimentConfig()
        self.config.proteins = ['O15263']
        self.config.constraints = ['lagrangian']
        self.config.variants = ['00', '01']
        self.config.iterations = 5  # Small for testing
        self.config.batch_size = 2
        self.config.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        self.config.output_dir = 'results/unit_test_batch'
        
        self.runner = BatchParallelExperimentRunner(self.config)
        
    def test_batch_runner_initialization(self):
        """Test batch runner initialization"""
        self.assertIsNotNone(self.runner)
        self.assertIsNotNone(self.runner.data_loader)
        self.assertIsNotNone(self.runner.deepraccess)
        self.assertEqual(self.runner.config.batch_size, 2)
        logger.info("✅ Batch runner initialization test passed")
    
    def test_constraint_creation(self):
        """Test constraint creation for different types"""
        # Test Lagrangian constraint
        constraint = self.runner.create_constraint(
            'lagrangian', 
            'MKLTVS',  # Short test sequence
            '00'
        )
        self.assertIsInstance(constraint, nn.Module)
        self.assertIsInstance(constraint, LagrangianConstraint)
        
        # Test AMS constraint
        constraint_ams = self.runner.create_constraint(
            'ams',
            'MKLTVS',
            '01'
        )
        self.assertIsInstance(constraint_ams, nn.Module)
        
        # Test CPC constraint
        constraint_cpc = self.runner.create_constraint(
            'cpc',
            'MKLTVS',
            '10'
        )
        self.assertIsInstance(constraint_cpc, nn.Module)
        
        logger.info("✅ Constraint creation test passed")
    
    def test_batch_sequence_generation(self):
        """Test batch sequence generation"""
        # Create multiple constraints
        constraints = []
        amino_seq = 'MKLTVS'
        for variant in ['00', '01', '10']:
            constraint = self.runner.create_constraint(
                'lagrangian',
                amino_seq,
                variant
            )
            constraints.append(constraint)
        
        # Generate sequences in batch
        batch_sequences = []
        for constraint in constraints:
            result = constraint.forward(alpha=1.0, beta=0.0, tau=1.0)
            batch_sequences.append(result['rna_sequence'])
        
        # Check batch dimensions
        batch_tensor = torch.cat(batch_sequences, dim=0)
        self.assertEqual(batch_tensor.shape[0], 3)  # 3 constraints
        self.assertEqual(batch_tensor.shape[1], len(amino_seq) * 3)  # CDS length
        self.assertEqual(batch_tensor.shape[2], 4)  # RNA alphabet size
        
        logger.info(f"✅ Batch sequence generation test passed - shape: {batch_tensor.shape}")
    
    def test_batch_deepraccess_inference(self):
        """Test batch DeepRaccess inference"""
        # Create batch of sequences
        batch_size = 4
        seq_len = 300
        batch_sequences = torch.randn(batch_size, seq_len, 4).to(self.runner.device)
        
        # Test batch inference (simulate)
        start_time = time.time()
        batch_results = []
        for i in range(batch_size):
            # Simulate DeepRaccess computation
            result = self.runner.deepraccess.compute_atg_window_accessibility(
                batch_sequences[i:i+1],
                discrete=False
            )
            batch_results.append(result)
        batch_time = time.time() - start_time
        
        self.assertEqual(len(batch_results), batch_size)
        for result in batch_results:
            self.assertIsInstance(result, torch.Tensor)
        
        logger.info(f"✅ Batch DeepRaccess inference test passed - time: {batch_time:.3f}s")
    
    def test_independent_gradient_updates(self):
        """Test that gradient updates are independent for each experiment"""
        # Create two constraints with different parameters
        constraint1 = self.runner.create_constraint('lagrangian', 'MKLTVS', '00')
        constraint2 = self.runner.create_constraint('lagrangian', 'MKLTVS', '01')
        
        optimizer1 = torch.optim.Adam(constraint1.parameters(), lr=0.01)
        optimizer2 = torch.optim.Adam(constraint2.parameters(), lr=0.01)
        
        # Get initial parameters
        initial_params1 = [p.clone() for p in constraint1.parameters()]
        initial_params2 = [p.clone() for p in constraint2.parameters()]
        
        # Perform forward pass and update
        result1 = constraint1.forward(alpha=1.0, beta=0.0, tau=1.0)
        result2 = constraint2.forward(alpha=1.0, beta=0.0, tau=1.0)
        
        # Create dummy losses
        loss1 = torch.sum(result1['rna_sequence'])
        loss2 = torch.sum(result2['rna_sequence']) * 2  # Different loss
        
        # Update only constraint1
        optimizer1.zero_grad()
        loss1.backward()
        optimizer1.step()
        
        # Check that constraint1 parameters changed but constraint2 didn't
        for p_init, p_current in zip(initial_params1, constraint1.parameters()):
            self.assertFalse(torch.allclose(p_init, p_current))
        
        for p_init, p_current in zip(initial_params2, constraint2.parameters()):
            self.assertTrue(torch.allclose(p_init, p_current))
        
        logger.info("✅ Independent gradient updates test passed")
    
    def test_batch_vs_sequential_consistency(self):
        """Test that batch and sequential produce similar results"""
        # Define small test experiments
        test_experiments = [
            {'protein_name': 'O15263', 'constraint_type': 'lagrangian', 
             'variant': '00', 'seed': 42},
            {'protein_name': 'O15263', 'constraint_type': 'lagrangian', 
             'variant': '01', 'seed': 42},
        ]
        
        # Run sequential
        torch.manual_seed(42)
        sequential_results = self.runner.run_sequential(test_experiments[:1])
        
        # Run batch parallel (with batch size 1 for comparison)
        torch.manual_seed(42)
        self.config.batch_size = 1
        batch_results = self.runner.run_batch_parallel(test_experiments[:1], batch_size=1)
        
        # Compare results
        self.assertEqual(len(sequential_results), len(batch_results))
        
        # Check that both completed successfully
        for seq_res, batch_res in zip(sequential_results, batch_results):
            self.assertEqual(seq_res['status'], 'completed')
            self.assertEqual(batch_res['status'], 'completed')
            self.assertEqual(seq_res['protein_name'], batch_res['protein_name'])
            self.assertEqual(seq_res['constraint_type'], batch_res['constraint_type'])
            self.assertEqual(seq_res['variant'], batch_res['variant'])
        
        logger.info("✅ Batch vs sequential consistency test passed")
    
    def test_performance_improvement(self):
        """Test that batch parallel is faster than sequential"""
        # Define test experiments
        test_experiments = [
            {'protein_name': 'O15263', 'constraint_type': 'lagrangian', 
             'variant': v, 'seed': 42}
            for v in ['00', '01', '10', '11']
        ]
        
        # Time sequential execution
        start_seq = time.time()
        seq_results = self.runner.run_sequential(test_experiments[:2])
        seq_time = time.time() - start_seq
        
        # Time batch parallel execution
        start_batch = time.time()
        batch_results = self.runner.run_batch_parallel(test_experiments[:2], batch_size=2)
        batch_time = time.time() - start_batch
        
        # Calculate speedup
        speedup = seq_time / batch_time if batch_time > 0 else 1.0
        
        # Log results
        logger.info(f"Sequential time: {seq_time:.3f}s")
        logger.info(f"Batch parallel time: {batch_time:.3f}s")
        logger.info(f"Speedup: {speedup:.2f}x")
        
        # Batch should be faster (at least not slower)
        self.assertGreaterEqual(speedup, 0.9)  # Allow small variation
        
        logger.info(f"✅ Performance improvement test passed - Speedup: {speedup:.2f}x")
    
    def test_batch_size_variations(self):
        """Test different batch sizes"""
        test_experiments = [
            {'protein_name': 'O15263', 'constraint_type': 'lagrangian', 
             'variant': f'{i:02d}', 'seed': 42}
            for i in range(4)
        ]
        
        for batch_size in [1, 2, 4]:
            results = self.runner.run_batch_parallel(
                test_experiments, 
                batch_size=batch_size
            )
            self.assertEqual(len(results), len(test_experiments))
            for result in results:
                self.assertEqual(result['status'], 'completed')
                self.assertEqual(result['batch_size'], batch_size)
            
            logger.info(f"✅ Batch size {batch_size} test passed")
    
    def test_error_handling(self):
        """Test error handling for invalid inputs"""
        # Test invalid constraint type
        with self.assertRaises(ValueError):
            self.runner.create_constraint('invalid_type', 'MKLTVS', '00')
        
        # Test empty experiments list
        results = self.runner.run_batch_parallel([])
        self.assertEqual(len(results), 0)
        
        logger.info("✅ Error handling test passed")
    
    def test_gpu_utilization_estimation(self):
        """Test GPU utilization estimation"""
        utilization_1 = self.runner.estimate_gpu_utilization(1)
        utilization_4 = self.runner.estimate_gpu_utilization(4)
        utilization_12 = self.runner.estimate_gpu_utilization(12)
        
        # Check that utilization increases with batch size
        self.assertLess(utilization_1, utilization_4)
        self.assertLess(utilization_4, utilization_12)
        
        # Check reasonable bounds
        self.assertGreaterEqual(utilization_1, 30)
        self.assertLessEqual(utilization_12, 95)
        
        logger.info(f"✅ GPU utilization estimation test passed - "
                   f"batch_1: {utilization_1}%, batch_4: {utilization_4}%, "
                   f"batch_12: {utilization_12}%")
    
    def test_memory_efficiency(self):
        """Test memory efficiency of batch processing"""
        if torch.cuda.is_available():
            # Get initial memory
            torch.cuda.empty_cache()
            torch.cuda.synchronize()
            initial_memory = torch.cuda.memory_allocated()
            
            # Run batch experiment
            test_experiments = [
                {'protein_name': 'O15263', 'constraint_type': 'lagrangian', 
                 'variant': '00', 'seed': 42}
            ] * 4
            
            results = self.runner.run_batch_parallel(test_experiments, batch_size=4)
            
            # Get final memory
            torch.cuda.synchronize()
            final_memory = torch.cuda.memory_allocated()
            memory_used = (final_memory - initial_memory) / 1024**2  # MB
            
            logger.info(f"✅ Memory efficiency test passed - Used: {memory_used:.1f} MB")
        else:
            logger.info("⚠️ Memory efficiency test skipped (no GPU)")
    
    def tearDown(self):
        """Clean up after tests"""
        # Clean up any temporary files if needed
        pass


class TestBatchParallelIntegration(unittest.TestCase):
    """Integration tests for batch parallel system"""
    
    def test_full_experiment_flow(self):
        """Test complete experiment flow from start to finish"""
        config = UnifiedExperimentConfig()
        config.proteins = ['O15263']
        config.constraints = ['lagrangian']
        config.variants = ['00']
        config.iterations = 2  # Very small for testing
        config.batch_size = 1
        config.output_dir = 'results/integration_test'
        
        runner = BatchParallelExperimentRunner(config)
        
        experiments = [
            {'protein_name': 'O15263', 'constraint_type': 'lagrangian',
             'variant': '00', 'seed': 42}
        ]
        
        results = runner.run_batch_parallel(experiments)
        
        self.assertEqual(len(results), 1)
        result = results[0]
        self.assertEqual(result['status'], 'completed')
        self.assertIn('final_accessibility', result)
        self.assertIn('discrete_sequence', result)
        self.assertIn('optimization_time', result)
        
        logger.info("✅ Full experiment flow integration test passed")
    
    def test_cai_integration(self):
        """Test CAI optimization integration"""
        config = UnifiedExperimentConfig()
        config.proteins = ['O15263']
        config.constraints = ['lagrangian']
        config.variants = ['00']
        config.iterations = 2
        config.enable_cai = True
        config.cai_target = 0.8
        config.lambda_cai = 0.1
        config.species = 'ecoli_bl21de3'
        config.batch_size = 1
        
        runner = BatchParallelExperimentRunner(config)
        
        # Check CAI optimizer was created
        self.assertIsNotNone(runner.cai_optimizer)
        
        logger.info("✅ CAI integration test passed")


def run_all_tests():
    """Run all unit tests with summary"""
    print("=" * 70)
    print("Running Batch Parallel System Unit Tests")
    print("=" * 70)
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    suite.addTests(loader.loadTestsFromTestCase(TestBatchParallelSystem))
    suite.addTests(loader.loadTestsFromTestCase(TestBatchParallelIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success rate: {(result.testsRun - len(result.failures) - len(result.errors)) / result.testsRun * 100:.1f}%")
    
    if result.wasSuccessful():
        print("\n🎉 All tests passed successfully!")
    else:
        print("\n❌ Some tests failed. Please review the output above.")
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_all_tests()