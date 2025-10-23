#!/usr/bin/env python3
"""







"""

import sys
import os
import time
import torch
import numpy as np
from pathlib import Path
import json
import logging
from typing import Dict, List, Any


sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.ray_runner.ray_experiment_runner import RayExperimentRunner
import ray

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CompatibilityTester:

    
    def __init__(self):
        self.original_runner = None
        self.ray_runner = None
        self.test_configs = []
        
    def setup(self):
        """初始化测试环境"""
        logger.info("Setting up test environment...")
        

        config = {
            'device': 'cuda' if torch.cuda.is_available() else 'cpu',
            'use_amp': True,
            'gradient_clip': 1.0,
            'enable_deferred_validation': False,
            'verbose': False
        }
        self.original_runner = UnifiedExperimentRunner(config)
        

        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
        self.ray_runner = RayExperimentRunner(num_workers=2, num_gpus_per_worker=0.5)
        

        self.test_configs = [

            {
                'protein': 'O15263',
                'constraint': 'lagrangian',
                'variant': '00',
                'iterations': 10,
                'seed': 42,
                'enable_cai': False
            },

            {
                'protein': 'O15263',
                'constraint': 'lagrangian',
                'variant': '11',
                'iterations': 10,
                'seed': 42,
                'enable_cai': True,
                'cai_target': 0.8,
                'lambda_cai': 0.1,
                'cai_species': 'ecoli_bl21de3'
            },

            {
                'protein': 'O15263',
                'constraint': 'ams',
                'variant': '01',
                'iterations': 10,
                'seed': 42,
                'enable_cai': False
            },
            {
                'protein': 'O15263',
                'constraint': 'cpc',
                'variant': '10',
                'iterations': 10,
                'seed': 42,
                'enable_cai': False
            }
        ]
        
        logger.info(f"Test environment ready with {len(self.test_configs)} test cases")
    
    def test_single_experiment_consistency(self, config: Dict[str, Any]) -> Dict[str, Any]:

        logger.info(f"Testing: {config['protein']}_{config['constraint']}_{config['variant']}")
        

        logger.info("Running original implementation...")
        

        self.original_runner.config['iterations'] = config['iterations']
        self.original_runner.config['enable_cai'] = config.get('enable_cai', False)
        self.original_runner.config['cai_target'] = config.get('cai_target', 0.8)
        self.original_runner.config['lambda_cai'] = config.get('lambda_cai', 0.1)
        self.original_runner.config['cai_species'] = config.get('cai_species', 'ecoli_bl21de3')
        
        start_time = time.time()
        original_result = self.original_runner.run_single_experiment(
            protein_name=config['protein'],
            constraint_type=config['constraint'],
            variant=config['variant'],
            seed=config['seed'],
            show_progress=False
        )
        original_time = time.time() - start_time
        

        logger.info("Running Ray implementation...")
        start_time = time.time()
        ray_results = self.ray_runner.run_experiments([config])
        ray_result = ray_results[0] if ray_results else None
        ray_time = time.time() - start_time
        

        comparison = {
            'config': config,
            'original_time': original_time,
            'ray_time': ray_time,
            'speedup': original_time / ray_time if ray_time > 0 else 0,
            'consistent': False,
            'differences': []
        }
        
        if not ray_result:
            comparison['differences'].append("Ray implementation failed to produce result")
            return comparison
        


        

        orig_acc = original_result.get('best_accessibility', float('inf'))
        ray_acc = ray_result.get('best_accessibility', float('inf'))
        
        if abs(orig_acc - ray_acc) > tolerance:
            comparison['differences'].append(
                f"Accessibility mismatch: original={orig_acc:.6f}, ray={ray_acc:.6f}"
            )
        

        orig_satisfied = original_result.get('amino_acids_match', False)
        ray_satisfied = ray_result.get('amino_acids_match', False)
        
        if orig_satisfied != ray_satisfied:
            comparison['differences'].append(
                f"Constraint satisfaction mismatch: original={orig_satisfied}, ray={ray_satisfied}"
            )
        

        if config.get('enable_cai', False):
            orig_cai = original_result.get('discrete_cai')
            ray_cai = ray_result.get('discrete_cai')
            
            if orig_cai and ray_cai and abs(orig_cai - ray_cai) > tolerance:
                comparison['differences'].append(
                    f"CAI mismatch: original={orig_cai:.4f}, ray={ray_cai:.4f}"
                )
        

        comparison['consistent'] = len(comparison['differences']) == 0
        

        if comparison['consistent']:
            logger.info(f"✅ PASSED - Results consistent (speedup: {comparison['speedup']:.2f}x)")
        else:
            logger.warning(f"❌ FAILED - Differences found:")
            for diff in comparison['differences']:
                logger.warning(f"  - {diff}")
        
        return comparison
    
    def test_parallel_execution(self) -> Dict[str, Any]:
        """测试并行执行效率"""
        logger.info("\n" + "="*50)
        logger.info("Testing parallel execution efficiency...")
        

        parallel_configs = []
        for i in range(8):
            config = self.test_configs[0].copy()
            config['seed'] = 42 + i
            parallel_configs.append(config)
        

        logger.info("Running 8 experiments serially with original implementation...")
        start_time = time.time()
        serial_results = []
        for config in parallel_configs:

            self.original_runner.config['iterations'] = config['iterations']
            self.original_runner.config['enable_cai'] = config.get('enable_cai', False)
            
            result = self.original_runner.run_single_experiment(
                protein_name=config['protein'],
                constraint_type=config['constraint'],
                variant=config['variant'],
                seed=config['seed'],
                show_progress=False
            )
            serial_results.append(result)
        serial_time = time.time() - start_time
        

        logger.info("Running 8 experiments in parallel with Ray...")
        start_time = time.time()
        ray_results = self.ray_runner.run_experiments(parallel_configs)
        parallel_time = time.time() - start_time
        
        speedup = serial_time / parallel_time if parallel_time > 0 else 0
        
        logger.info(f"Serial time: {serial_time:.2f}s")
        logger.info(f"Parallel time: {parallel_time:.2f}s")
        logger.info(f"Speedup: {speedup:.2f}x")
        
        return {
            'num_experiments': len(parallel_configs),
            'serial_time': serial_time,
            'parallel_time': parallel_time,
            'speedup': speedup,
            'efficiency': speedup / 2  # 2 workers
        }
    
    def test_gpu_memory_usage(self) -> Dict[str, Any]:

        if not torch.cuda.is_available():
            logger.info("GPU not available, skipping memory test")
            return {'gpu_available': False}
        
        logger.info("\n" + "="*50)
        logger.info("Testing GPU memory usage...")
        

        torch.cuda.empty_cache()
        initial_memory = torch.cuda.memory_allocated()
        

        config = self.test_configs[0]
        self.ray_runner.run_experiments([config])
        

        peak_memory = torch.cuda.max_memory_allocated()
        memory_used = peak_memory - initial_memory
        
        logger.info(f"Initial memory: {initial_memory / 1024**2:.2f} MB")
        logger.info(f"Peak memory: {peak_memory / 1024**2:.2f} MB")
        logger.info(f"Memory used: {memory_used / 1024**2:.2f} MB")
        
        return {
            'gpu_available': True,
            'initial_memory_mb': initial_memory / 1024**2,
            'peak_memory_mb': peak_memory / 1024**2,
            'memory_used_mb': memory_used / 1024**2
        }
    
    def run_all_tests(self) -> Dict[str, Any]:
        """运行所有测试"""
        logger.info("\n" + "="*60)
        logger.info("STARTING RAY COMPATIBILITY TESTS")
        logger.info("="*60)
        
        all_results = {
            'consistency_tests': [],
            'parallel_test': None,
            'gpu_memory_test': None,
            'summary': {
                'total_tests': 0,
                'passed': 0,
                'failed': 0
            }
        }
        

        self.setup()
        

        logger.info("\n" + "="*50)
        logger.info("Running consistency tests...")
        for config in self.test_configs:
            result = self.test_single_experiment_consistency(config)
            all_results['consistency_tests'].append(result)
            all_results['summary']['total_tests'] += 1
            if result['consistent']:
                all_results['summary']['passed'] += 1
            else:
                all_results['summary']['failed'] += 1
        

        all_results['parallel_test'] = self.test_parallel_execution()
        

        all_results['gpu_memory_test'] = self.test_gpu_memory_usage()
        

        logger.info("\n" + "="*60)
        logger.info("TEST SUMMARY")
        logger.info("="*60)
        logger.info(f"Total tests: {all_results['summary']['total_tests']}")
        logger.info(f"Passed: {all_results['summary']['passed']}")
        logger.info(f"Failed: {all_results['summary']['failed']}")
        
        if all_results['parallel_test']:
            logger.info(f"Parallel speedup: {all_results['parallel_test']['speedup']:.2f}x")
        
        if all_results['summary']['failed'] == 0:
            logger.info("\n✅ ALL TESTS PASSED - Ray implementation is compatible!")
        else:
            logger.warning(f"\n⚠️ {all_results['summary']['failed']} tests failed - needs investigation")
        

        output_file = Path("ray_compatibility_test_results.json")
        with open(output_file, 'w') as f:
            json.dump(all_results, f, indent=2, default=str)
        logger.info(f"\nDetailed results saved to: {output_file}")
        

        self.cleanup()
        
        return all_results
    
    def cleanup(self):

        logger.info("\nCleaning up...")
        if self.ray_runner:
            self.ray_runner.shutdown()
        logger.info("Cleanup complete")


def main():
    """主测试函数"""
    tester = CompatibilityTester()
    results = tester.run_all_tests()
    

    if results['summary']['failed'] == 0:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()