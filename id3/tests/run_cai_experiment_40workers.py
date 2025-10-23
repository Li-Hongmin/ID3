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
        proteins=['O15263', 'P00004', 'P01308', 'P01825',
                 'P04637', 'P0CG48', 'P0DTC2', 'P0DTC9',




        learning_rate=0.05,

        device='cuda',
        output_dir=f'./results/cai_formal_{timestamp}',

        cai_target=0.8,
        lambda_cai=1.0,
        species='ecoli_bl21de3'
    )
    



    logger.info(f"    {', '.join(config.proteins[:4])}...")
    logger.info(f"    {', '.join(config.proteins[4:8])}...")
    logger.info(f"    {', '.join(config.proteins[8:])}")




    

    num_workers = 40


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



    






    logger.info("")
    
    results = runner.run_batch(experiments, num_workers=num_workers)
    

    success_count = sum(1 for r in results if r.get('status') == 'completed')
    failed_count = sum(1 for r in results if r.get('status') == 'failed')
    
    logger.info("=" * 60)


    
    if failed_count > 0:

    

    if failed_count > 0:

        error_types = {}
        for r in results:
            if r.get('status') == 'failed':
                error = r.get('error', 'Unknown error')
                error_types[error] = error_types.get(error, 0) + 1
        
        for error, count in sorted(error_types.items(), key=lambda x: x[1], reverse=True)[:5]:

    

    if success_count > 0:

        

        constraint_results = {}
        for r in results:
            if r.get('status') == 'completed':
                constraint = r.get('constraint_type', 'Unknown')
                if constraint not in constraint_results:
                    constraint_results[constraint] = {
                        'count': 0,
                        'total_accessibility': 0,
                        'total_cai': 0,

                    }
                constraint_results[constraint]['count'] += 1
                constraint_results[constraint]['total_accessibility'] += r.get('final_accessibility', 0)
                

                cai_value = r.get('final_cai_discrete', 0)
                if cai_value != 'N/A' and cai_value is not None:
                    cai_float = float(cai_value)
                    constraint_results[constraint]['total_cai'] += cai_float

                        constraint_results[constraint]['cai_success'] += 1
        
        for constraint, stats in sorted(constraint_results.items()):
            avg_access = stats['total_accessibility'] / stats['count']
            avg_cai = stats['total_cai'] / stats['count'] if stats['total_cai'] > 0 else 0
            cai_success_rate = stats['cai_success'] / stats['count'] * 100
            



            if avg_cai > 0:


    

    summary_file = Path(config.output_dir) / 'experiment_summary.txt'
    summary_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(summary_file, 'w') as f:






        if failed_count > 0:

        f.write("\n")
        

        for constraint, stats in sorted(constraint_results.items()):
            avg_access = stats['total_accessibility'] / stats['count']
            avg_cai = stats['total_cai'] / stats['count'] if stats['total_cai'] > 0 else 0




    


    
    return success_count == total_experiments

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)