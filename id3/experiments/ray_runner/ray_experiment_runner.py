#!/usr/bin/env python3
"""
Ray-based Experiment Runner for ID3-DeepRaccess






"""

import ray
import torch
import logging
import json
from pathlib import Path
from typing import Dict, List, Any, Optional
from datetime import datetime
from tqdm import tqdm
import numpy as np
import time

logger = logging.getLogger(__name__)


class RayExperimentRunner:
    """

    





    """
    
    def __init__(self, num_workers: int = 4, num_gpus_per_worker: float = 0.25):
        """

        
        Args:


        """
        self.num_workers = num_workers
        self.num_gpus_per_worker = num_gpus_per_worker
        

        if not ray.is_initialized():
            ray.init(ignore_reinit_error=True)
            logger.info("Ray initialized successfully")
        

        resources = ray.available_resources()
        self.available_gpus = resources.get("GPU", 0)
        self.available_cpus = resources.get("CPU", 1)
        
        logger.info(f"Available resources - GPUs: {self.available_gpus}, CPUs: {self.available_cpus}")
        

        if self.available_gpus > 0 and self.num_gpus_per_worker > 0:
            max_workers_by_gpu = int(self.available_gpus / self.num_gpus_per_worker)
            if max_workers_by_gpu < self.num_workers:
                self.num_workers = max_workers_by_gpu
                logger.warning(f"Adjusted workers to {self.num_workers} based on available GPUs")
        elif self.available_gpus == 0 and self.num_gpus_per_worker > 0:
            logger.warning("No GPUs available, setting num_gpus_per_worker to 0")
            self.num_gpus_per_worker = 0
    
    def run_experiments(self, 
                        experiments: List[Dict[str, Any]], 
                        output_dir: Optional[Path] = None,
                        save_individual: bool = True) -> List[Dict[str, Any]]:
        """

        
        Args:



            
        Returns:

        """
        from .ray_worker import ExperimentWorker
        
        start_time = time.time()
        

        if output_dir is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_dir = Path(f"results/{timestamp}_ray_experiments")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Starting {len(experiments)} experiments with {self.num_workers} workers")
        logger.info(f"Output directory: {output_dir}")
        

        workers = [
            ExperimentWorker.remote(num_gpus=self.num_gpus_per_worker) 
            for _ in range(self.num_workers)
        ]
        

        futures = []
        experiment_info = []
        
        for i, exp_config in enumerate(experiments):
            worker = workers[i % self.num_workers]
            future = worker.run_experiment.remote(exp_config)
            futures.append(future)
            experiment_info.append(exp_config)
        

        results = []
        failed_experiments = []
        

        with tqdm(total=len(futures), desc="Running experiments") as pbar:
            pending = futures.copy()
            
            while pending:

                ready, pending = ray.wait(pending, num_returns=1, timeout=1.0)
                
                for future in ready:
                    try:
                        result = ray.get(future)
                        results.append(result)
                        

                        if save_individual and result:
                            self._save_single_result(result, output_dir)
                        
                        pbar.update(1)
                    except Exception as e:

                        idx = futures.index(future)
                        exp_config = experiment_info[idx]
                        
                        logger.error(f"Experiment failed: {exp_config}")
                        logger.error(f"Error: {str(e)}")
                        
                        failed_experiments.append({
                            'config': exp_config,
                            'error': str(e)
                        })
                        pbar.update(1)
        

        summary = self._generate_summary(results, failed_experiments, start_time)
        

        with open(output_dir / "summary.json", 'w') as f:
            json.dump(summary, f, indent=2)
        

        if failed_experiments:
            with open(output_dir / "failed_experiments.json", 'w') as f:
                json.dump(failed_experiments, f, indent=2)
        

        if results:
            self._generate_performance_table(results, output_dir)
        
        logger.info(f"Completed {len(results)}/{len(experiments)} experiments successfully")
        logger.info(f"Total time: {time.time() - start_time:.2f} seconds")
        
        return results
    
    def _save_single_result(self, result: Dict[str, Any], output_dir: Path):

        filename = f"{result['protein_name']}_{result['constraint_type']}_{result['variant']}_seed{result['seed']}.json"
        filepath = output_dir / filename
        
        with open(filepath, 'w') as f:
            json.dump(result, f, indent=2)
    
    def _generate_summary(self, results: List[Dict[str, Any]], 
                         failed_experiments: List[Dict[str, Any]],
                         start_time: float) -> Dict[str, Any]:
        """生成实验汇总"""
        if not results:
            return {
                'total_experiments': 0,
                'successful': 0,
                'failed': len(failed_experiments),
                'message': 'No successful experiments'
            }
        

        accessibilities = [r['best_accessibility'] for r in results if 'best_accessibility' in r and r['best_accessibility'] != float('inf')]
        
        summary = {
            'total_experiments': len(results) + len(failed_experiments),
            'successful': len(results),
            'failed': len(failed_experiments),
            'total_time': time.time() - start_time,
            'performance': {
                'average_accessibility': np.mean(accessibilities) if len(accessibilities) > 0 else None,
                'best_accessibility': min(accessibilities) if len(accessibilities) > 0 else None,
                'worst_accessibility': max(accessibilities) if len(accessibilities) > 0 else None,
                'std_accessibility': np.std(accessibilities) if len(accessibilities) > 0 else None
            }
        }
        

        constraint_stats = {}
        for result in results:
            constraint = result.get('constraint_type', 'unknown')
            if constraint not in constraint_stats:
                constraint_stats[constraint] = []
            if 'best_accessibility' in result:
                constraint_stats[constraint].append(result['best_accessibility'])
        
        summary['constraint_performance'] = {}
        for constraint, values in constraint_stats.items():
            if len(values) > 0:
                summary['constraint_performance'][constraint] = {
                    'count': len(values),
                    'average': np.mean(values),
                    'best': min(values),
                    'worst': max(values)
                }
            else:
                summary['constraint_performance'][constraint] = {
                    'count': 0,
                    'average': None,
                    'best': None,
                    'worst': None
                }
        
        return summary
    
    def _generate_performance_table(self, results: List[Dict[str, Any]], output_dir: Path):

        import pandas as pd
        

        data = []
        for result in results:
            data.append({
                'protein': result['protein_name'],
                'constraint': result['constraint_type'],
                'variant': result['variant'],
                'best_accessibility': result.get('best_accessibility', np.nan),
                'amino_acids_match': result.get('amino_acids_match', False),
                'seed': result.get('seed', 42)
            })
        
        df = pd.DataFrame(data)
        

        pivot_table = df.pivot_table(
            values='best_accessibility',
            index='protein',
            columns='constraint',

        )
        

        pivot_table.to_csv(output_dir / "performance_table.csv")
        

        stats_table = df.groupby(['constraint', 'variant']).agg({
            'best_accessibility': ['mean', 'min', 'max', 'std', 'count']
        }).round(4)
        
        stats_table.to_csv(output_dir / "statistics_table.csv")
        
        logger.info("Performance tables saved")
    
    def shutdown(self):
        """关闭Ray"""
        if ray.is_initialized():
            ray.shutdown()
            logger.info("Ray shutdown complete")


def run_batch_experiments_ray(preset: str = "quick-test", 
                              num_workers: int = 4,
                              num_gpus_per_worker: float = 0.25) -> List[Dict[str, Any]]:
    """

    
    Args:



        
    Returns:

    """
    import sys
    from pathlib import Path
    sys.path.append(str(Path(__file__).parent.parent.parent.parent))
    
    from id3.experiments.configs.unified_experiment_config import ExperimentPresets
    

    preset_map = {
        'quick-test': ExperimentPresets.quick_test,
        'quick-test-cai': ExperimentPresets.quick_test_cai,
        'full-12x12': ExperimentPresets.full_12x12,
        'full-12x12-cai': ExperimentPresets.full_12x12_cai,
    }
    

    if preset == 'quick-test-both':

        configs = ExperimentPresets.quick_test_both()
        all_results = []
        for config in configs:
            experiments = []
            for protein in config.proteins:
                for constraint in config.constraints:
                    for variant in config.variants:
                        experiments.append({
                            'protein': protein,
                            'constraint': constraint,
                            'variant': variant,
                            'iterations': config.iterations,
                            'seed': config.base_seed,
                            'enable_cai': config.enable_cai,
                            'cai_target': getattr(config, 'cai_target', 0.8),
                            'lambda_cai': getattr(config, 'lambda_cai', 0.1)
                        })
            
            runner = RayExperimentRunner(num_workers=num_workers, num_gpus_per_worker=num_gpus_per_worker)
            try:
                results = runner.run_experiments(experiments)
                all_results.extend(results)
            finally:
                runner.shutdown()
        return all_results
    
    if preset not in preset_map:
        raise ValueError(f"Unknown preset: {preset}. Available: {list(preset_map.keys())}")
    
    config = preset_map[preset]()
    

    experiments = []
    for protein in config.proteins:
        for constraint in config.constraints:
            for variant in config.variants:
                experiments.append({
                    'protein': protein,
                    'constraint': constraint,
                    'variant': variant,
                    'iterations': config.iterations,
                    'seed': config.base_seed,
                    'enable_cai': config.enable_cai,
                    'cai_target': getattr(config, 'cai_target', 0.8),
                    'lambda_cai': getattr(config, 'lambda_cai', 0.1)
                })
    
    logger.info(f"Generated {len(experiments)} experiments from preset '{preset}'")
    

    runner = RayExperimentRunner(num_workers=num_workers, num_gpus_per_worker=num_gpus_per_worker)
    
    try:
        results = runner.run_experiments(experiments)
        return results
    finally:
        runner.shutdown()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Ray-based ID3 experiment runner")
    parser.add_argument('--preset', default='quick-test', help='Experiment preset')
    parser.add_argument('--workers', type=int, default=4, help='Number of parallel workers')
    parser.add_argument('--gpu-fraction', type=float, default=0.25, help='GPU fraction per worker')
    
    args = parser.parse_args()
    
    results = run_batch_experiments_ray(
        preset=args.preset,
        num_workers=args.workers,
        num_gpus_per_worker=args.gpu_fraction
    )
    
    print(f"\n✅ Completed {len(results)} experiments")