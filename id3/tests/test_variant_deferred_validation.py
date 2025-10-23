#!/usr/bin/env python3
"""






"""

import time
import logging
from pathlib import Path

from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_variant_deferred_validation():

    

    variants_to_test = [
        {
            'variant': '11', 


        },
        {
            'variant': '01', 


        },
        {
            'variant': '10', 


        }
    ]
    

    protein_name = 'O15263'
    constraint_type = 'lagrangian'
    seed = 42

    




    logger.info("=" * 80)
    
    all_results = {}
    
    for variant_config in variants_to_test:
        variant = variant_config['variant']
        variant_name = variant_config['name']
        

        logger.info(f"   {variant_config['description']}")
        logger.info("-" * 60)
        
        variant_results = {}
        

        for enable_deferred in [False, True]:

            logger.info(f"\n  📊 {mode_name}:")
            

            exp_config = UnifiedExperimentConfig()
            exp_config.iterations = test_iterations
            exp_config.learning_rate = 0.01
            exp_config.use_mixed_precision = True
            exp_config.enable_cai = False
            

            runner = UnifiedExperimentRunner(exp_config)
            

            start_time = time.time()
            
            try:
                result = runner.run_single_experiment(
                    protein_name=protein_name,
                    constraint_type=constraint_type,
                    variant=variant,
                    seed=seed,
                    show_progress=False,
                    enable_deferred_validation=enable_deferred
                )
                
                elapsed_time = time.time() - start_time
                
                variant_results[mode_name] = {
                    'elapsed_time': elapsed_time,
                    'iterations_per_second': test_iterations / elapsed_time,
                    'final_accessibility': result.get('final_accessibility', float('inf')),
                    'status': 'success'
                }
                



                
            except Exception as e:
                variant_results[mode_name] = {
                    'elapsed_time': float('inf'),
                    'iterations_per_second': 0,
                    'status': f'failed: {str(e)}'
                }

        



        
        if traditional.get('status') == 'success' and deferred.get('status') == 'success':
            traditional_time = traditional['elapsed_time']
            deferred_time = deferred['elapsed_time']
            speedup = traditional_time / deferred_time
            improvement_percent = (speedup - 1) * 100
            


            

            if variant == '11':
                if abs(improvement_percent) < 10:

                else:

            else:
                if improvement_percent > 30:

                elif improvement_percent > 10:

                else:

            

            traditional_acc = traditional.get('final_accessibility', float('inf'))
            deferred_acc = deferred.get('final_accessibility', float('inf'))
            error = abs(traditional_acc - deferred_acc)
            
            if error < 1e-3:

            else:

        
        all_results[variant_name] = variant_results
    

    logger.info("\n" + "=" * 80)

    logger.info("=" * 80)
    
    for variant_config in variants_to_test:
        variant_name = variant_config['name']
        variant_results = all_results.get(variant_name, {})
        


        
        if traditional.get('status') == 'success' and deferred.get('status') == 'success':
            speedup = traditional['elapsed_time'] / deferred['elapsed_time']
            improvement = (speedup - 1) * 100
            
            logger.info(f"\n{variant_name}:")



    




    
    return all_results

if __name__ == "__main__":
    results = test_variant_deferred_validation()
    

    import json
    output_file = Path("variant_deferred_validation_results.json")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    

