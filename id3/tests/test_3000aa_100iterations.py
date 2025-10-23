#!/usr/bin/env python3
"""


"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import torch
import numpy as np
import time
import random
from datetime import datetime

from id3.optimizers.cai.sado import SADOOptimizer


def test_3000aa_100iterations():

    
    print("\n" + "="*70)

    print("="*70)

    

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    amino_sequence = ''.join(random.choice(amino_acids) for _ in range(3000))

    

    torch.manual_seed(42)
    pi_accessibility = torch.rand(3000, 6)
    pi_accessibility = pi_accessibility / pi_accessibility.sum(dim=1, keepdim=True)
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


    

    iteration_times = []
    cai_values = []
    

    print("-" * 50)
    

    total_start = time.time()
    

    for i in range(100):

        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        

        iter_start = time.time()
        
        distribution, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility.clone(),
            target_cai=0.8,
            use_binary_search=False,

            gamma=0.3
        )
        
        iter_time = time.time() - iter_start
        iteration_times.append(iter_time)
        cai_values.append(metadata['final_cai'])
        

        if (i + 1) % 10 == 0:
            avg_time = np.mean(iteration_times)
            avg_cai = np.mean(cai_values)


    

    total_time = time.time() - total_start
    

    print("\n" + "="*70)

    print("="*70)
    






    






    





    



    first_half_avg = np.mean(iteration_times[:50])
    second_half_avg = np.mean(iteration_times[50:])



    

    print("\n" + "="*70)

    print("="*70)
    
    if total_time < 100:


    elif total_time < 300:


    else:


    

    
    return total_time, iteration_times, cai_values


if __name__ == '__main__':
    try:
        total_time, times, cais = test_3000aa_100iterations()

    except KeyboardInterrupt:

    except Exception as e:

        import traceback
        traceback.print_exc()