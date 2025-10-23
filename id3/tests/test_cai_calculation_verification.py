#!/usr/bin/env python3
"""

"""

import numpy as np
from pathlib import Path
import json


import sys
sys.path.append(str(Path(__file__).parent.parent.parent))

from id3.utils.constants import amino_acids_to_codons


def verify_cai_calculation():

    

    weights_file = Path("data/codon_references/ecoli_bl21de3_wi_weights_comparison.json")
    with open(weights_file, 'r') as f:
        data = json.load(f)
    
    wi_table = data['wi_table']
    

    rna_wi_table = {}
    for dna_codon, weight in wi_table.items():
        rna_codon = dna_codon.replace('T', 'U')
        rna_wi_table[rna_codon] = weight
    

    print("=" * 60)
    

    for aa, codons in amino_acids_to_codons.items():
        weights = []
        codon_weight_map = {}
        
        for codon in codons:
            if codon in rna_wi_table:
                weight = rna_wi_table[codon]
            else:
                weight = 0.1
            weights.append(weight)
            codon_weight_map[codon] = weight
        
        max_weight = max(weights)

        for codon, weight in codon_weight_map.items():

            print(f"  {codon}: {weight:.4f}{marker}")
    

    test_sequences = {



    }
    
    print("\n" + "=" * 60)

    print("=" * 60)
    
    for name, sequence in test_sequences.items():
        log_weights = []
        
        for aa in sequence:
            if aa not in amino_acids_to_codons:
                continue
            
            codons = amino_acids_to_codons[aa]
            weights = []
            
            for codon in codons:
                if codon in rna_wi_table:
                    weights.append(rna_wi_table[codon])
                else:
                    weights.append(0.1)
            
            max_weight = max(weights)
            log_weights.append(np.log(max(max_weight, 1e-10)))
        

        if log_weights:
            mean_log_weight = np.mean(log_weights)
            theoretical_max_cai = np.exp(mean_log_weight)
        else:
            theoretical_max_cai = 0.0
        
        print(f"\n{name} ({sequence}):")

        

        for i, aa in enumerate(sequence):
            if aa in amino_acids_to_codons:
                codons = amino_acids_to_codons[aa]
                weights = [rna_wi_table.get(c, 0.1) for c in codons]
                max_weight = max(weights)



if __name__ == "__main__":
    verify_cai_calculation()