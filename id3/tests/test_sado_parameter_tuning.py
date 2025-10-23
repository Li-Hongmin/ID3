"""



"""

import sys
import torch
import numpy as np
import time
import logging


sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_sado_with_different_gamma():

    

    sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVT"
    target_cai = 0.8
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    logger.info("="*80)

    logger.info("="*80)



    

    gamma_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]
    

    seq_len = len(sequence)
    max_codons = 6
    
    distribution = torch.rand(seq_len, max_codons, device=device)
    valid_mask = torch.zeros(seq_len, max_codons, device=device)
    
    for pos, aa in enumerate(sequence):
        if aa in amino_acids_to_codons:
            num_codons = len(amino_acids_to_codons[aa])
            valid_mask[pos, :num_codons] = 1.0
            distribution[pos] = distribution[pos] * valid_mask[pos]
            if distribution[pos].sum() > 0:
                distribution[pos] = distribution[pos] / distribution[pos].sum()
    
    codon_indices = torch.zeros_like(distribution, dtype=torch.long)
    

    logger.info("-"*60)

    logger.info("-"*60)
    
    best_gamma = None
    best_cai = 0
    
    for gamma in gamma_values:

        operator = CAIEnhancementOperator(
            method='sado',
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=sequence
        )
        

        start_time = time.time()
        result_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=distribution,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_mask,
            codon_indices=codon_indices,
            target_cai=target_cai,

        )
        elapsed_time = (time.time() - start_time) * 1000
        
        final_cai = metadata.get('final_cai', 0.0)
        satisfied = metadata.get('constraint_satisfied', False)
        
        logger.info(f"{gamma:8.2f} | {final_cai:8.4f} | {'✅' if satisfied else '❌':>10} | {elapsed_time:10.2f}")
        

        if satisfied and (best_gamma is None or final_cai < best_cai):
            best_gamma = gamma
            best_cai = final_cai
    
    logger.info("-"*60)
    
    if best_gamma is not None:

    else:

    

    logger.info(f"\n{'='*60}")

    logger.info(f"{'='*60}\n")
    
    if best_gamma is not None:
        test_gamma = best_gamma
    else:

    
    operator = CAIEnhancementOperator(
        method='sado',
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=sequence
    )
    
    cais = []
    satisfied_count = 0
    
    for i in range(10):

        test_dist = distribution + torch.rand_like(distribution) * 0.1
        test_dist = test_dist * valid_mask
        test_dist = test_dist / (test_dist.sum(dim=-1, keepdim=True) + 1e-10)
        
        result_dist, metadata = operator.apply_cai_enhancement(
            pi_accessibility=test_dist,
            amino_acid_sequence=sequence,
            valid_codon_mask=valid_mask,
            codon_indices=codon_indices,
            target_cai=target_cai,
            gamma=test_gamma
        )
        
        final_cai = metadata.get('final_cai', 0.0)
        satisfied = metadata.get('constraint_satisfied', False)
        
        cais.append(final_cai)
        if satisfied:
            satisfied_count += 1
        

    
    avg_cai = np.mean(cais)
    std_cai = np.std(cais)
    



    
    stats = operator.get_statistics()
    if 'diversity' in stats:
        diversity = stats['diversity']


    

    logger.info(f"\n{'='*60}")

    logger.info(f"{'='*60}")
    
    if best_gamma is not None:

    else:


    






if __name__ == "__main__":
    test_sado_with_different_gamma()