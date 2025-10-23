"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='exact_compare')


def simulate_v4_with_while_loop(optimizer, pi_accessibility, target_cai, gamma):
    """



    """
    

    indices = optimizer._initial_optimization(
        pi_accessibility, target_cai, gamma
    )
    

    initial_cai = optimizer._compute_cai_from_indices(indices)
    

    current_cai = initial_cai
    iterations = 0
    
    while current_cai < target_cai and iterations < 100:
        improvements = []
        

        for pos in range(len(indices)):
            if pos >= len(optimizer.codon_choices):
                continue
            
            current_idx = indices[pos]
            choices = optimizer.codon_choices[pos]
            

            for new_idx in range(len(choices)):
                if new_idx <= current_idx:
                    continue
                

                test_indices = indices.copy()
                test_indices[pos] = new_idx
                test_cai = optimizer._compute_cai_from_indices(test_indices)
                
                if test_cai > current_cai:
                    cai_gain = test_cai - current_cai
                    improvements.append((cai_gain, pos, new_idx, test_cai))
        
        if not improvements:
            break
        

        improvements.sort()
        

        selected = None
        for gain, pos, new_idx, new_cai in improvements:
            if new_cai >= target_cai:
                selected = (gain, pos, new_idx, new_cai)
                break
        

        if selected is None and improvements:
            selected = improvements[0]
        
        if selected:
            _, pos, new_idx, new_cai = selected
            indices[pos] = new_idx
            current_cai = new_cai
            iterations += 1
        else:
            break
    
    return indices, initial_cai, current_cai, iterations


def test_exact_comparison():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    


    logger.info("-" * 60)
    

    torch.manual_seed(42)
    pi_test = torch.rand(seq_len, num_codons, device=device)
    pi_test = pi_test * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
    
    target_cai = 0.8
    

    test_gammas = [0.25, 0.3, 0.35, 0.4]
    

    logger.info("")
    
    for gamma in test_gammas:
        logger.info(f"gamma = {gamma}:")
        

        optimizer.reset()
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        current_cai = metadata['final_cai']
        

        optimizer.reset()
        indices_v4, initial_cai, final_cai_v4, iterations = simulate_v4_with_while_loop(
            optimizer, pi_test, target_cai, gamma
        )
        




        
        if current_cai < target_cai:

        if final_cai_v4 >= target_cai and current_cai < target_cai:

        
        logger.info("")
    

    logger.info("="*60)

    logger.info("-" * 60)
    
    gamma = 0.3
    optimizer.reset()
    

    indices = optimizer._initial_optimization(pi_test, target_cai, gamma)
    cai = optimizer._compute_cai_from_indices(indices)
    

    
    if cai < target_cai:






    


    logger.info("-" * 60)
    
    left, right = 0.3, 0.5
    best_gamma = 0.4
    

        mid = (left + right) / 2
        optimizer.reset()
        _, meta = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=mid
        )
        
        if meta['final_cai'] >= target_cai:
            best_gamma = mid
            right = mid
        else:
            left = mid
    
    optimizer.reset()
    _, meta = optimizer.optimize(
        pi_accessibility=pi_test,
        target_cai=target_cai,
        amino_acid_sequence=amino_sequence,
        valid_codon_mask=valid_mask,
        gamma=best_gamma
    )
    

    

    logger.info("\n" + "="*80)

    logger.info("="*80)




    logger.info("")






def main():
    """主函数"""
    test_exact_comparison()


if __name__ == "__main__":
    main()