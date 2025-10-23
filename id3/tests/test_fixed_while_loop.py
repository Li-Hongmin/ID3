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


logger = setup_logging(level='INFO', name='fixed_loop')


def apply_while_loop_fix(optimizer, indices, target_cai, max_iterations=10):
    """


    """
    
    current_cai = optimizer._compute_cai_from_indices(indices)
    iterations = 0
    
    logger.info(f"初始CAI: {current_cai:.4f}, 目标: {target_cai}")
    
    while current_cai < target_cai and iterations < max_iterations:
        improvements = []
        

        for pos in range(len(indices)):
            if pos >= len(optimizer.codon_choices):
                continue
            
            current_idx = indices[pos]
            choices = optimizer.codon_choices[pos]
            

            for new_idx in range(len(choices)):
                if new_idx == current_idx:
                    continue
                

                test_indices = indices.copy()
                test_indices[pos] = new_idx
                test_cai = optimizer._compute_cai_from_indices(test_indices)
                
                if test_cai > current_cai:
                    improvements.append({
                        'pos': pos,
                        'new_idx': new_idx,
                        'new_cai': test_cai,
                        'cai_gain': test_cai - current_cai
                    })
        
        if not improvements:
            logger.info(f"  没有找到改进，停止")
            break
        

        improvements.sort(key=lambda x: x['cai_gain'])
        

        selected = None
        

        for imp in improvements:
            if imp['new_cai'] >= target_cai:

                selected = imp
                logger.info(f"  找到达到目标的改进: 位置{imp['pos']}, "
                          f"CAI {current_cai:.4f} → {imp['new_cai']:.4f}")
                break
        

        if selected is None and improvements:
            selected = improvements[-1]
            logger.info(f"  选择最大改进: 位置{selected['pos']}, "
                      f"CAI {current_cai:.4f} → {selected['new_cai']:.4f}")
        
        if selected:
            indices[selected['pos']] = selected['new_idx']
            current_cai = selected['new_cai']
            iterations += 1
        else:
            break
    
    return indices, current_cai, iterations


def test_fixed_while_loop():

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
    

    torch.manual_seed(42)
    pi_test = torch.rand(seq_len, num_codons, device=device)
    pi_test = pi_test * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
    
    target_cai = 0.8
    

    test_gammas = [0.25, 0.3, 0.35]
    

    logger.info("-" * 60)
    
    for gamma in test_gammas:
        logger.info(f"\ngamma = {gamma}:")
        
        optimizer.reset()
        

        indices = optimizer._initial_optimization(pi_test, target_cai, gamma)
        initial_cai = optimizer._compute_cai_from_indices(indices)
        

        final_indices, final_cai, iterations = apply_while_loop_fix(
            optimizer, indices.copy(), target_cai
        )
        


        
        if final_cai >= target_cai:

        else:

    

    logger.info("\n" + "="*80)

    logger.info("="*80)






def main():
    """主函数"""
    test_fixed_while_loop()


if __name__ == "__main__":
    main()