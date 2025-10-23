#!/usr/bin/env python3
"""



"""

import time
import logging
from pathlib import Path

# Import the unified experiment system
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_deferred_validation_performance():

    

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

    seed = 42

    





    logger.info("=" * 70)
    
    results = {}
    
    for i, config in enumerate(test_configs):

        logger.info(f"   {config['description']}")
        

        exp_config = UnifiedExperimentConfig()
        exp_config.iterations = test_iterations
        exp_config.learning_rate = 0.01

        

        runner = UnifiedExperimentRunner(exp_config)
        

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
            results[config["name"]] = {
                'elapsed_time': elapsed_time,
                'iterations_per_second': test_iterations / elapsed_time,
                'final_accessibility': result.get('final_accessibility', float('inf')),
                'status': 'success'
            }
            



            
        except Exception as e:
            results[config["name"]] = {
                'elapsed_time': float('inf'),
                'iterations_per_second': 0,
                'final_accessibility': float('inf'),
                'status': f'failed: {str(e)}'
            }

    

    logger.info("\n" + "=" * 70)

    logger.info("=" * 70)
    
    traditional_time = results[test_configs[0]["name"]]['elapsed_time']
    deferred_time = results[test_configs[1]["name"]]['elapsed_time']
    
    if traditional_time != float('inf') and deferred_time != float('inf'):
        speedup = traditional_time / deferred_time
        improvement_percent = (speedup - 1) * 100
        


        logger.info(f"")

        
        if improvement_percent > 30:

        elif improvement_percent > 10:

        elif improvement_percent > -5:

        else:

        

        traditional_acc = results[test_configs[0]["name"]]['final_accessibility']
        deferred_acc = results[test_configs[1]["name"]]['final_accessibility']
        
        if abs(traditional_acc - deferred_acc) < 1e-3:

        else:



    
    else:

    
    return results

if __name__ == "__main__":
    results = test_deferred_validation_performance()
    

    import json
    output_file = Path("deferred_validation_test_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    

