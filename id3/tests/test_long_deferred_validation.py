#!/usr/bin/env python3
"""






"""

import time
import logging
import torch
import numpy as np
from pathlib import Path

# Import the unified experiment system
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_long_deferred_validation():

    

    test_configs = [
        {

            "enable_deferred_validation": False,
        },
        {

            "enable_deferred_validation": True,
        },
    ]
    

    protein_name = 'O15263'
    constraint_type = 'lagrangian'
    variant = '11'  # alpha=1, beta=1
    seed = 42
    test_iterations = 10000
    





    logger.info("=" * 70)
    
    results = {}
    
    for i, config in enumerate(test_configs):

        

        exp_config = UnifiedExperimentConfig()
        exp_config.iterations = test_iterations
        exp_config.learning_rate = 0.01


        exp_config.gradient_clip = 1.0
        

        runner = UnifiedExperimentRunner(exp_config)
        

        torch.cuda.empty_cache()
        start_memory = torch.cuda.memory_allocated()
        start_time = time.time()
        
        try:

            result = runner.run_single_experiment(
                protein_name=protein_name,
                constraint_type=constraint_type,
                variant=variant,
                seed=seed,

                enable_deferred_validation=config["enable_deferred_validation"]
            )
            

            elapsed_time = time.time() - start_time
            end_memory = torch.cuda.memory_allocated()
            memory_usage = (end_memory - start_memory) / 1024**2  # MB
            

            trajectory = result.get('trajectory', {})
            accessibility_values = trajectory.get('accessibility', [])
            
            if accessibility_values:
                initial_acc = accessibility_values[0]
                final_acc = accessibility_values[-1]
                min_acc = min(accessibility_values)
                improvement = initial_acc - final_acc
                

                acc_array = np.array(accessibility_values)
                std_dev = np.std(acc_array)
                has_nan = np.isnan(acc_array).any()
                has_inf = np.isinf(acc_array).any()
                

                extremely_low_values = (acc_array < -10).sum()
                negative_values = (acc_array < 0).sum()
                
                results[config["name"]] = {
                    'elapsed_time': elapsed_time,
                    'memory_usage_mb': memory_usage,
                    'iterations_per_second': test_iterations / elapsed_time,
                    'initial_accessibility': initial_acc,
                    'final_accessibility': final_acc,
                    'best_accessibility': min_acc,
                    'improvement': improvement,
                    'std_dev': std_dev,
                    'has_numerical_issues': has_nan or has_inf,
                    'extremely_low_count': int(extremely_low_values),
                    'negative_count': int(negative_values),
                    'total_iterations': len(accessibility_values),
                    'status': 'success'
                }
                






                

                if has_nan or has_inf:

                if extremely_low_values > 0:

                if negative_values > 0:

                    
            else:

                results[config["name"]] = {'status': 'no_trajectory'}
                
        except Exception as e:

            results[config["name"]] = {
                'elapsed_time': float('inf'),
                'status': f'failed: {str(e)}'
            }
    

    logger.info("\n" + "=" * 70)

    logger.info("=" * 70)
    
    traditional = results.get(test_configs[0]["name"], {})
    deferred = results.get(test_configs[1]["name"], {})
    
    if traditional.get('status') == 'success' and deferred.get('status') == 'success':
        traditional_time = traditional['elapsed_time']
        deferred_time = deferred['elapsed_time']
        
        speedup = traditional_time / deferred_time
        improvement_percent = (speedup - 1) * 100
        


        logger.info(f"")

        

        mem_traditional = traditional.get('memory_usage_mb', 0)
        mem_deferred = deferred.get('memory_usage_mb', 0)
        if mem_traditional > 0 and mem_deferred > 0:
            mem_saving = ((mem_traditional - mem_deferred) / mem_traditional) * 100

        

        traditional_final = traditional.get('final_accessibility', float('inf'))
        deferred_final = deferred.get('final_accessibility', float('inf'))
        
        if abs(traditional_final - deferred_final) < 1e-3:

        else:
            error = abs(traditional_final - deferred_final)



        


        
        for name, data in [(test_configs[0]["name"], traditional), (test_configs[1]["name"], deferred)]:
            logger.info(f"   {name}:")




        

        if improvement_percent > 30:

        elif improvement_percent > 10:

        else:

    
    else:

    
    return results

if __name__ == "__main__":
    results = test_long_deferred_validation()
    

    import json
    output_file = Path("long_deferred_validation_test_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    

