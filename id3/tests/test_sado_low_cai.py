"""

"""

import sys
sys.path.append('/home/yunqi/ideas/ID3_DeepRaccess_CAI_Paper')

import torch
import numpy as np
import time
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.constants import amino_acids_to_codons

def test_sado_low_cai():

    

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    print("="*60)

    print("="*60)


    print()
    

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
    

    target_cais = [0.5, 0.6, 0.65, 0.7, 0.75, 0.8]
    

    print("-" * 70)
    
    for target_cai in target_cais:
        results = {}
        

        for method in ['binary_search', 'sado']:
            operator = CAIEnhancementOperator(
                method=method,
                species='ecoli_bl21de3', 
                device=device,
                amino_acid_sequence=sequence
            )
            
            start_time = time.time()
            if method == 'sado':
                result_dist, metadata = operator.apply_cai_enhancement(
                    pi_accessibility=distribution,
                    amino_acid_sequence=sequence,
                    valid_codon_mask=valid_mask,
                    codon_indices=codon_indices,
                    target_cai=target_cai,

                )
            else:
                result_dist, metadata = operator.apply_cai_enhancement(
                    pi_accessibility=distribution,
                    amino_acid_sequence=sequence,
                    valid_codon_mask=valid_mask,
                    codon_indices=codon_indices,
                    target_cai=target_cai
                )
            elapsed_time = (time.time() - start_time) * 1000
            
            final_cai = metadata.get('final_cai', 0.0)
            satisfied = metadata.get('constraint_satisfied', False)
            
            results[method] = {
                'time': elapsed_time,
                'cai': final_cai,
                'satisfied': satisfied
            }
        

        for method in ['binary_search', 'sado']:
            r = results[method]
            print(f"{target_cai:8.2f} | {method:>15} | {r['cai']:8.4f} | {'✅' if r['satisfied'] else '❌':>6} | {r['time']:8.1f} | ", end="")
            
            if method == 'sado' and 'binary_search' in results:
                speedup = results['binary_search']['time'] / r['time']
                print(f"{speedup:8.1f}x")
            else:
                print(f"{'':>8}")
        
        print()
    
    print("="*70)






if __name__ == "__main__":
    test_sado_low_cai()