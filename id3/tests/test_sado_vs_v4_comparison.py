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


sys.path.append(str(Path(__file__).parent.parent.parent / 
                   'experiments/cai_enhancement_theory/03_SADO_algorithm'))
from algorithms.sado_v4_incremental import IncrementalSADO as SADOOptimizerV4


logger = setup_logging(level='INFO', name='test_comparison')


def compute_cai_from_indices(indices, optimizer):

    cai_product = 1.0
    
    if hasattr(optimizer, 'codon_choices'):

        for pos, idx in enumerate(indices):
            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):
                weight = optimizer.codon_choices[pos][idx]['weight']
                cai_product *= weight
    else:

        for pos, idx in enumerate(indices):
            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):

                weight = optimizer.codon_choices[pos][idx]['weight']
                cai_product *= weight
    
    seq_len = len(indices)
    cai = cai_product ** (1.0 / seq_len) if seq_len > 0 else 0.0
    return cai


def compute_probability(indices, pi_distribution, optimizer):
    """计算序列的概率"""
    prob = 1.0
    
    for pos, idx in enumerate(indices):
        if hasattr(optimizer, 'codon_choices'):

            if pos < len(optimizer.codon_choices) and idx < len(optimizer.codon_choices[pos]):
                choice = optimizer.codon_choices[pos][idx]
                orig_idx = choice.get('original_local_index', idx)
                if orig_idx < pi_distribution.shape[1]:
                    p = pi_distribution[pos, orig_idx].item()
                    prob *= p
        else:

            if pos < pi_distribution.shape[0] and idx < pi_distribution.shape[1]:
                p = pi_distribution[pos, idx].item()
                prob *= p
    
    return prob


def test_comparison():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    


    optimizer_current = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    optimizer_v4 = SADOOptimizerV4(
        sequence=amino_sequence,
        target_cai=0.8
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    


    logger.info("-" * 60)
    
    for test_id in range(5):

        

        torch.manual_seed(42 + test_id)
        np.random.seed(42 + test_id)
        
        pi_test = torch.rand(seq_len, num_codons, device=device)
        pi_test = pi_test * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
        

        for gamma in [0.3, 0.6]:
            logger.info(f"\n  gamma/alpha = {gamma}:")
            

            optimizer_current.reset()
            dist_current, meta_current = optimizer_current.optimize(
                pi_accessibility=pi_test,
                target_cai=0.8,
                amino_acid_sequence=amino_sequence,
                valid_codon_mask=valid_mask,
                gamma=gamma
            )
            
            indices_current = optimizer_current.last_indices
            cai_current = meta_current['final_cai']
            prob_current = compute_probability(indices_current, pi_test, optimizer_current)
            

            optimizer_v4.reset()

            pi_numpy = pi_test.detach().cpu().numpy()
            indices_v4 = optimizer_v4.optimize(pi_numpy, alpha=gamma)
            

            cai_v4 = optimizer_v4._calculate_cai_from_indices(indices_v4)
            


            prob_v4 = 1.0
            for pos in range(len(indices_v4)):
                if pos < pi_numpy.shape[0] and indices_v4[pos] < pi_numpy.shape[1]:
                    prob_v4 *= pi_numpy[pos, indices_v4[pos]]
            



            

            cai_diff = abs(cai_current - cai_v4)
            prob_ratio = prob_current / prob_v4 if prob_v4 > 0 else float('inf')
            

            
            if cai_current > cai_v4 + 0.1:

            elif cai_v4 > cai_current + 0.1:

            
            if prob_v4 > prob_current * 10:

    

    logger.info("\n" + "="*80)

    logger.info("="*80)







def test_detailed_v4_behavior():
    """详细分析v4的行为"""
    logger.info("\n" + "="*80)
    logger.info("详细分析v4版本的优化过程")
    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    

    optimizer_v4 = SADOOptimizerV4(
        sequence=amino_sequence,
        target_cai=0.8
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
    
    torch.manual_seed(42)
    np.random.seed(42)
    
    pi_test = np.random.rand(seq_len, num_codons)

    for pos in range(seq_len):
        aa = amino_sequence[pos]
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])

            if num_valid < num_codons:
                pi_test[pos, num_valid:] = 0
            if pi_test[pos].sum() > 0:
                pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
    

    logger.info("\n不同alpha值下v4的CAI结果：")
    logger.info("-" * 40)
    
    for alpha in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        optimizer_v4.reset()
        indices = optimizer_v4.optimize(pi_test, alpha=alpha)
        cai = optimizer_v4._calculate_cai_from_indices(indices)
        

        prob = 1.0
        for pos in range(len(indices)):
            if indices[pos] < pi_test.shape[1]:
                prob *= pi_test[pos, indices[pos]]
        
        satisfied = "✅" if cai >= 0.8 else "❌"
        logger.info(f"alpha={alpha:.1f}: CAI={cai:.4f}, log_prob={np.log(prob):.2f} {satisfied}")
    
    logger.info("\n观察：v4版本在不同alpha下都能接近目标CAI=0.8")


def main():

    test_comparison()
    test_detailed_v4_behavior()


if __name__ == "__main__":
    main()