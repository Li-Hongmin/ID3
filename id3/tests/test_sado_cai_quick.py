"""



"""

import sys
import torch
import numpy as np
import time
import logging
from typing import Dict, List, Tuple


sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def quick_test_sado():

    

    test_sequences = [



    ]
    
    target_cais = [0.7, 0.75, 0.8, 0.85]
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    logger.info("="*80)

    logger.info("="*80)



    logger.info("")
    

    for seq_idx, sequence in enumerate(test_sequences, 1):
        logger.info(f"\n{'='*60}")

        logger.info(f"{'='*60}")
        

        for target_cai in target_cais:

            logger.info("-"*40)
            

            results = {}
            
            for method in ['binary_search', 'sado']:

                operator = CAIEnhancementOperator(
                    method=method,
                    species='ecoli_bl21de3',
                    device=device,
                    amino_acid_sequence=sequence
                )
                

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
                

                start_time = time.time()
                result_dist, metadata = operator.apply_cai_enhancement(
                    pi_accessibility=distribution,
                    amino_acid_sequence=sequence,
                    valid_codon_mask=valid_mask,
                    codon_indices=codon_indices,
                    target_cai=target_cai
                )

                
                results[method] = {
                    'time': elapsed_time,
                    'final_cai': metadata.get('final_cai', 0.0),
                    'satisfied': metadata.get('constraint_satisfied', False),
                    'metadata': metadata
                }
                
                logger.info(f"{method:15s}: CAI={results[method]['final_cai']:.4f}, "


            

            if 'binary_search' in results and 'sado' in results:
                speedup = results['binary_search']['time'] / results['sado']['time']

    

    logger.info(f"\n{'='*60}")

    logger.info(f"{'='*60}\n")
    

    target_cai = 0.8
    
    for method in ['binary_search', 'sado']:
        operator = CAIEnhancementOperator(
            method=method,
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=sequence
        )
        

        if method == 'sado':
            operator.reset()
        
        seq_len = len(sequence)
        times = []
        cais = []
        
        for i in range(10):

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
            
            start_time = time.time()
            result_dist, metadata = operator.apply_cai_enhancement(
                pi_accessibility=distribution,
                amino_acid_sequence=sequence,
                valid_codon_mask=valid_mask,
                codon_indices=codon_indices,
                target_cai=target_cai
            )
            elapsed_time = (time.time() - start_time) * 1000
            
            times.append(elapsed_time)
            cais.append(metadata.get('final_cai', 0.0))
        
        avg_time = np.mean(times)
        avg_cai = np.mean(cais)
        

        

        if method == 'sado':
            stats = operator.get_statistics()
            if 'diversity' in stats:
                diversity = stats['diversity']


    

    logger.info(f"\n{'='*80}")

    logger.info(f"{'='*80}")






if __name__ == "__main__":
    quick_test_sado()