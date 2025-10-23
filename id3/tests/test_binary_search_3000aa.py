#!/usr/bin/env python3
"""

"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import time
import torch
import numpy as np

from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging
from id3.optimizers.cai import BinarySearchCAIOptimizer

logger = setup_logging(level='INFO', name='binary_3000aa')


def test_binary_search_3000aa():

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    


    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    weights = np.array([8, 3, 6, 6, 7, 4, 2, 6, 5, 9, 2, 4, 2, 4, 5, 7, 5, 1, 3, 7])
    weights = weights / weights.sum()
    test_sequence = ''.join(np.random.choice(amino_acids, 3000, p=weights))
    


    seq_len = 3000
    target_dist = torch.zeros(seq_len, 6, device=device)
    valid_mask = torch.zeros(seq_len, 6, dtype=torch.bool, device=device)
    
    for pos, aa in enumerate(test_sequence):
        if aa in amino_acids_to_codons:
            num_codons = len(amino_acids_to_codons[aa])
            if num_codons > 0:

                for i in range(min(num_codons, 6)):
                    valid_mask[pos, i] = True
                

                probs = torch.rand(num_codons, device=device)
                probs = probs / probs.sum()
                for i in range(min(num_codons, 6)):
                    target_dist[pos, i] = probs[i]
    


    optimizer = BinarySearchCAIOptimizer(
        species='ecoli_bl21de3',
        device=device,

    )
    


    times = []
    cais = []
    gammas = []
    iterations_list = []
    
    for run in range(5):
        start_time = time.time()
        
        result, metadata = optimizer.optimize(
            pi_accessibility=target_dist,
            target_cai=0.8,
            amino_acid_sequence=test_sequence,
            valid_codon_mask=valid_mask
        )
        
        elapsed_time = time.time() - start_time
        
        times.append(elapsed_time)
        cais.append(metadata.get('final_cai', 0.0))
        gammas.append(metadata.get('gamma', 0.0))
        iterations_list.append(metadata.get('iterations', 0))
        

    

    avg_time = np.mean(times)
    std_time = np.std(times)
    avg_cai = np.mean(cais)
    avg_gamma = np.mean(gammas)
    avg_iterations = np.mean(iterations_list)
    
    logger.info("=" * 80)

    logger.info("=" * 80)






    

    logger.info("")


    epsilon = 1e-4
    theoretical_ops = 3000 * k * np.log(1/epsilon)

    logger.info(f"  n=3000, k≈4, ε=1e-4")

    logger.info(f"  log(1/ε) ≈ {np.log(1/epsilon):.1f}")
    

    logger.info("")


    logger.info(f"  PrecomputedSwitchingSearch: ~2900ms - O(n·k·log(n·k))")

    
    if avg_time < 0.1:  # <100ms

    elif avg_time < 1.0:  # <1s

    else:

    
    return {
        'avg_time_ms': avg_time * 1000,
        'avg_cai': avg_cai,
        'avg_gamma': avg_gamma,
        'avg_iterations': avg_iterations
    }


if __name__ == "__main__":

    results = test_binary_search_3000aa()

