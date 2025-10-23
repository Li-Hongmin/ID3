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


logger = setup_logging(level='INFO', name='test_tradeoff')


def compute_sequence_probability(indices, pi_distribution, optimizer):

    prob = 1.0
    
    for pos, idx in enumerate(indices):
        if pos < len(optimizer.codon_choices):
            choice = optimizer.codon_choices[pos][idx]
            orig_idx = choice.get('original_local_index', idx)
            
            if orig_idx < pi_distribution.shape[1]:
                p = pi_distribution[pos, orig_idx].item()
                prob *= p
    
    return prob


def find_boundary_sequence(optimizer, pi_distribution, target_cai=0.8, 
                          amino_sequence=None, valid_mask=None):
    """寻找接近CAI边界的序列"""

    for gamma in [0.05, 0.1, 0.15, 0.2]:
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_distribution,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        
        if metadata['final_cai'] >= target_cai and metadata['final_cai'] < target_cai + 0.05:
            return optimizer.last_indices, metadata['final_cai'], gamma
    

    return optimizer.last_indices, metadata['final_cai'], gamma


def test_tradeoff():

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
    

    logger.info("-" * 60)
    
    results = []
    

    gammas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    for gamma in gammas:

        
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        

        if optimizer.last_indices is not None:
            prob = compute_sequence_probability(
                optimizer.last_indices, pi_test, optimizer
            )
            log_prob = np.log(prob) if prob > 0 else float('-inf')
        else:
            prob = 0
            log_prob = float('-inf')
        
        results.append({
            'gamma': gamma,
            'cai': metadata['final_cai'],
            'prob': prob,
            'log_prob': log_prob,
            'satisfied': metadata['final_cai'] >= 0.8
        })
        
        status = "✅" if metadata['final_cai'] >= 0.8 else "❌"
        logger.info(f"gamma={gamma:.1f}: CAI={metadata['final_cai']:.4f}, "
                   f"log_prob={log_prob:.2f}, P={prob:.2e} {status}")
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    satisfied = [r for r in results if r['satisfied']]
    
    if len(satisfied) >= 2:

        sorted_sat = sorted(satisfied, key=lambda x: x['cai'])
        
        lowest_cai = sorted_sat[0]
        highest_cai = sorted_sat[-1]
        


        logger.info(f"  CAI: {lowest_cai['cai']:.4f}")


        logger.info(f"  gamma: {lowest_cai['gamma']:.1f}")
        

        logger.info(f"  CAI: {highest_cai['cai']:.4f}")


        logger.info(f"  gamma: {highest_cai['gamma']:.1f}")
        

        if lowest_cai['prob'] > 0 and highest_cai['prob'] > 0:
            prob_ratio = lowest_cai['prob'] / highest_cai['prob']

            
            if prob_ratio > 1:



            else:

    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    if satisfied:

        best_prob = max(satisfied, key=lambda x: x['prob'])

        logger.info(f"  gamma: {best_prob['gamma']:.1f}")
        logger.info(f"  CAI: {best_prob['cai']:.4f}")


        

        if best_prob['cai'] > 0.85:


    

    logger.info("\n" + "="*80)

    logger.info("="*80)









def test_gamma_search():
    """测试寻找最优gamma值"""
    logger.info("\n" + "="*80)
    logger.info("寻找最优gamma值（最大化概率同时满足CAI≥0.8）")
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
    

    logger.info("\n对10个随机分布寻找最优gamma：")
    logger.info("-" * 60)
    
    optimal_gammas = []
    
    for i in range(10):
        torch.manual_seed(100 + i)
        pi = torch.rand(seq_len, num_codons, device=device)
        pi = pi * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi[pos] = pi[pos] / pi[pos].sum()
        

        left, right = 0.0, 1.0
        best_gamma = 0.5
        best_prob = 0
        best_cai = 0
        
        for _ in range(10):
            mid = (left + right) / 2
            optimizer.reset()
            
            _, metadata = optimizer.optimize(
                pi_accessibility=pi,
                target_cai=0.8,
                amino_acid_sequence=amino_sequence,
                valid_codon_mask=valid_mask,
                gamma=mid
            )
            
            if optimizer.last_indices is not None:
                prob = compute_sequence_probability(
                    optimizer.last_indices, pi, optimizer
                )
            else:
                prob = 0
            
            cai = metadata['final_cai']
            
            if cai >= 0.8:

                if prob > best_prob:
                    best_gamma = mid
                    best_prob = prob
                    best_cai = cai
                right = mid
            else:

                left = mid
        
        optimal_gammas.append(best_gamma)
        logger.info(f"分布{i+1}: 最优gamma={best_gamma:.3f}, "
                   f"CAI={best_cai:.4f}, P={best_prob:.2e}")
    
    logger.info(f"\n平均最优gamma: {np.mean(optimal_gammas):.3f} "
               f"(±{np.std(optimal_gammas):.3f})")
    
    logger.info("\n建议：")
    logger.info(f"对于arg max P(S|π) s.t. CAI≥0.8的目标，")
    logger.info(f"gamma应该设置在{np.mean(optimal_gammas):.2f}左右，")
    logger.info(f"而不是默认的0.6，以获得更高的条件概率")


def main():

    test_tradeoff()
    test_gamma_search()


if __name__ == "__main__":
    main()