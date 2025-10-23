#!/usr/bin/env python3
"""



"""

import torch
import logging
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner
from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_fixed_variants():

    
    variants = [

        {'variant': '01', 'name': 'STE (α=0,β=1)', 'expected_alpha': 0.0, 'expected_beta': 1.0},

        {'variant': '11', 'name': 'Gumbel+STE (α=1,β=1)', 'expected_alpha': 0.1, 'expected_beta': 1.0},
    ]
    

    logger.info("="*60)
    

    exp_config = UnifiedExperimentConfig()
    exp_config.iterations = 10
    exp_config.enable_cai = False
    

    runner = UnifiedExperimentRunner(exp_config)
    
    # Hook _forward_pass to capture alpha/beta values
    captured_values = {}
    
    original_forward_pass = runner._forward_pass
    def debug_forward_pass(constraint, constraint_type, amino_acid_sequence,
                          alpha, beta, protein_info, utr5_tensor=None, utr3_tensor=None, 
                          enable_discrete_monitoring=True):
        
        # Capture values
        if 'alpha_beta_values' not in captured_values:
            captured_values['alpha_beta_values'] = []
        captured_values['alpha_beta_values'].append({'alpha': alpha, 'beta': beta})
        
        # Call original method
        return original_forward_pass(constraint, constraint_type, amino_acid_sequence,
                                   alpha, beta, protein_info, utr5_tensor, utr3_tensor, 
                                   enable_discrete_monitoring)
    
    runner._forward_pass = debug_forward_pass
    

    for variant_config in variants:
        variant = variant_config['variant']
        name = variant_config['name']
        expected_alpha = variant_config['expected_alpha']
        expected_beta = variant_config['expected_beta']
        


        
        captured_values.clear()
        
        try:
            result = runner.run_single_experiment(
                protein_name='O15263',
                constraint_type='lagrangian',
                variant=variant,
                seed=42,
                show_progress=False,

            )
            

            if captured_values.get('alpha_beta_values'):
                values = captured_values['alpha_beta_values']
                alphas = [v['alpha'] for v in values]
                betas = [v['beta'] for v in values]
                


                

                alpha_is_fixed = len(set(f"{a:.3f}" for a in alphas)) == 1
                beta_is_fixed = len(set(f"{b:.3f}" for b in betas)) == 1
                

                alpha_correct = abs(alphas[0] - expected_alpha) < 0.001
                beta_correct = abs(betas[0] - expected_beta) < 0.001
                
                if alpha_is_fixed and beta_is_fixed:

                else:

                
                if alpha_correct and beta_correct:

                else:

                

            
        except Exception as e:

    
    logger.info(f"\n{'='*60}")







if __name__ == "__main__":
    test_fixed_variants()