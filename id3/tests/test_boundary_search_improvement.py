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


logger = setup_logging(level='INFO', name='boundary_improve')


def improved_boundary_search(optimizer, indices, target_cai, max_iterations=50):
    """



    """
    
    current_cai = optimizer._compute_cai_from_indices(indices)
    iterations = 0
    
    logger.info(f"开始边界搜索: 初始CAI={current_cai:.4f}, 目标={target_cai}")
    
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
            logger.info(f"  没有找到改进，停在CAI={current_cai:.4f}")
            break
        

        improvements.sort(key=lambda x: x['cai_gain'])
        

        selected = None
        for imp in improvements:
            if imp['new_cai'] >= target_cai and imp['new_cai'] < target_cai + 0.05:
                selected = imp
                break
        

        if selected is None:
            for imp in improvements:
                if imp['new_cai'] >= target_cai:
                    selected = imp
                    break
        

        if selected is None:
            selected = improvements[-1]
        
        if selected:
            indices[selected['pos']] = selected['new_idx']
            old_cai = current_cai
            current_cai = selected['new_cai']
            iterations += 1
            
            if iterations % 5 == 0:
                logger.info(f"  迭代{iterations}: CAI {old_cai:.4f} → {current_cai:.4f}")
    
    return indices, current_cai, iterations


def test_difficult_cases():

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
    

    optimizer.reset()
    indices = optimizer._initial_optimization(pi_test)
    initial_cai = optimizer._compute_cai_from_indices(indices)
    

    

    final_indices, final_cai, iterations = improved_boundary_search(
        optimizer, indices.copy(), 0.8
    )
    

    
    if final_cai >= 0.8:

    else:

        


        

        max_possible_indices = np.zeros(seq_len, dtype=np.int32)
        for pos in range(seq_len):
            if pos < len(optimizer.codon_choices) and optimizer.codon_choices[pos]:

                max_possible_indices[pos] = 0
        
        max_possible_cai = optimizer._compute_cai_from_indices(max_possible_indices)

        
        if max_possible_cai < 0.8:


    


    logger.info("-" * 60)
    
    success_count = 0
    cai_values = []
    iteration_counts = []
    
    for i in range(100):
        torch.manual_seed(1000 + i)
        pi = torch.rand(seq_len, num_codons, device=device)
        pi = pi * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi[pos] = pi[pos] / pi[pos].sum()
        

        optimizer.reset()
        indices = optimizer._initial_optimization(pi)
        

        final_indices, final_cai, iterations = improved_boundary_search(
            optimizer, indices.copy(), 0.8, max_iterations=50
        )
        
        cai_values.append(final_cai)
        iteration_counts.append(iterations)
        
        if final_cai >= 0.8:
            success_count += 1
    



    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    if success_count >= 95:

    else:

    






def main():
    """主函数"""
    test_difficult_cases()


if __name__ == "__main__":
    main()