#!/usr/bin/env python3
"""

"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pathlib import Path
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner


def test_cai_output_directory():


    

    config = UnifiedExperimentConfig(

        constraints=['lagrangian'],
        variants=['00'],

        seeds=1,

        cai_target=0.8,
        lambda_cai=1.0,
        species='ecoli_bl21de3',
        verbose=True
    )
    

    output_dir = config.get_output_dir()

    



    

    runner_config = config.to_dict()
    runner_config['output_dir'] = str(output_dir)
    

    print(f"  - enable_cai: {runner_config.get('enable_cai')}")
    print(f"  - output_dir: {runner_config.get('output_dir')}")
    

    runner = UnifiedExperimentRunner(runner_config)
    

    



    

    


    config_no_cai = UnifiedExperimentConfig(
        proteins=['O15263'],
        constraints=['lagrangian'],
        variants=['00'],
        iterations=10,
        seeds=1,

        verbose=True
    )
    
    output_dir_no_cai = config_no_cai.get_output_dir()

    



    

    runner_config_no_cai = config_no_cai.to_dict()
    runner_config_no_cai['output_dir'] = str(output_dir_no_cai)
    

    runner_no_cai = UnifiedExperimentRunner(runner_config_no_cai)
    

    



    

    





    

if __name__ == "__main__":
    test_cai_output_directory()