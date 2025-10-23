#!/usr/bin/env python3
"""

"""

import sys
import torch
from batch_parallel_experiments.massive_async_runner import MassiveAsyncExperiment
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig

def test_massive_async_cai():

    

    config = UnifiedExperimentConfig(




        device='cuda',

        cai_target=0.8,


    )
    






    

    runner = MassiveAsyncExperiment(config, max_batch_size=4, num_workers=2)
    

    num_experiments = runner.setup_all_experiments(
        config.proteins, 
        config.constraints, 
        config.variants
    )
    

    


    results = runner.run_massive_parallel()
    






    

    if results['experiment_overview']['successful'] > 0:

        



                  f"{exp_result['constraint_type']}-{exp_result['variant']} "

    else:

    
    return results

if __name__ == "__main__":
    test_massive_async_cai()