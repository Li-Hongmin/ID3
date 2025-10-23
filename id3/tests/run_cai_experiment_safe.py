#!/usr/bin/env python3
"""

"""

import sys
import os
import logging
from pathlib import Path
from datetime import datetime


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


sys.path.insert(0, str(Path(__file__).parent))

def main():

    from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
    from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    logger.info("=" * 60)
    

    config = UnifiedExperimentConfig(




        learning_rate=0.05,

        device='cuda',
        output_dir=f'./results/cai_safe_{timestamp}',

        cai_target=0.8,
        lambda_cai=1.0,
        species='ecoli_bl21de3'
    )
    







    




    logger.info("=" * 60)
    

    runner = UnifiedExperimentRunner(config.to_dict())
    

    experiments = []
    for protein in config.proteins:
        for constraint in config.constraints:
            for variant in config.variants:
                for seed in config.seeds:
                    experiments.append({
                        'protein_name': protein,
                        'constraint_type': constraint,
                        'variant': variant,
                        'seed': seed
                    })
    
    total_experiments = len(experiments)


    



    
    results = runner.run_batch(experiments, num_workers=num_workers)
    

    success_count = sum(1 for r in results if r.get('status') == 'completed')
    failed_count = sum(1 for r in results if r.get('status') == 'failed')
    
    logger.info("=" * 60)



    

    if failed_count > 0:

        error_types = {}
        for r in results:
            if r.get('status') == 'failed':
                error = r.get('error', 'Unknown error')
                error_types[error] = error_types.get(error, 0) + 1
        
        for error, count in error_types.items():

    

    if success_count > 0:

        

        constraint_results = {}
        for r in results:
            if r.get('status') == 'completed':
                constraint = r.get('constraint_type', 'Unknown')
                if constraint not in constraint_results:
                    constraint_results[constraint] = {
                        'count': 0,
                        'total_accessibility': 0,
                        'total_cai': 0
                    }
                constraint_results[constraint]['count'] += 1
                constraint_results[constraint]['total_accessibility'] += r.get('final_accessibility', 0)

                cai_value = r.get('final_cai_discrete', 0)
                if cai_value != 'N/A' and cai_value is not None:
                    constraint_results[constraint]['total_cai'] += float(cai_value)
        
        for constraint, stats in constraint_results.items():
            avg_access = stats['total_accessibility'] / stats['count']
            avg_cai = stats['total_cai'] / stats['count'] if stats['total_cai'] > 0 else 0


            if avg_cai > 0:

    

    summary_file = Path(config.output_dir) / 'experiment_summary.txt'
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_file, 'w') as f:





    

    
    return success_count == total_experiments

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)