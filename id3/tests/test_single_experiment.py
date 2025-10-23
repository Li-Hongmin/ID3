#!/usr/bin/env python3
"""


"""

import json
import sys
import time
from pathlib import Path
from datetime import datetime

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.experiments.configs.unified_experiment_config import UnifiedExperimentConfig
from id3.experiments.core.unified_experiment_runner import UnifiedExperimentRunner

def run_single_test():

    
    print("=" * 80)

    print("=" * 80)
    

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"results/{timestamp}_single_lagrangian_test")
    output_dir.mkdir(parents=True, exist_ok=True)
    

    config = UnifiedExperimentConfig(

        constraints=['lagrangian'],


        seeds=1,
        base_seed=42,

        verbose=True,
        output_dir=str(output_dir)
    )
    






    


    start_time = time.time()
    
    runner = UnifiedExperimentRunner(config)
    

    experiment = {
        'protein_name': 'P04637',
        'constraint_type': 'lagrangian',
        'variant': '01',
        'seed': 42
    }
    
    result = runner.run_single_experiment(
        protein_name='P04637',
        constraint_type='lagrangian',
        variant='01',
        seed=42
    )
    
    elapsed_time = time.time() - start_time
    

    print(f"\n" + "=" * 80)

    print("=" * 80)
    
    if result and result.get('status') == 'completed':

        





        




        

        if 'trajectory' in result:
            traj = result['trajectory']
            if 'accessibility' in traj:
                import numpy as np
                best_idx = np.argmin(traj['accessibility'])
                best_value = traj['accessibility'][best_idx]



                

                if 'discrete_sequences' in traj and best_idx < len(traj['discrete_sequences']):
                    best_seq = traj['discrete_sequences'][best_idx]

                    

                    from id3.utils.constants import amino_acids_to_codons
                    
                    def rna_to_aa(rna_seq):
                        codon_to_aa = {}
                        for aa, codons in amino_acids_to_codons.items():
                            for codon in codons:
                                codon_to_aa[codon] = aa
                        
                        aa_seq = []
                        for i in range(0, len(rna_seq) - 2, 3):
                            codon = rna_seq[i:i+3]
                            aa_seq.append(codon_to_aa.get(codon, 'X'))
                        return ''.join(aa_seq)
                    
                    best_aa = rna_to_aa(best_seq)
                    expected_aa = result.get('expected_amino_acids', '')
                    
                    if expected_aa:
                        matches = sum(1 for a, e in zip(best_aa, expected_aa) if a == e)
                        match_rate = matches / len(expected_aa) * 100 if expected_aa else 0

                        
                        if match_rate < 100:

                        else:

            

            if 'discrete_sequences' in traj:
                sequences = traj['discrete_sequences']
                unique_sequences = len(set(sequences))
                uniqueness_rate = unique_sequences / len(sequences) * 100 if sequences else 0




                
                if uniqueness_rate < 10:

                else:

        

        result_file = output_dir / "experiment_result.json"
        with open(result_file, 'w') as f:

            result_summary = {k: v for k, v in result.items() if k != 'trajectory'}
            if 'trajectory' in result:
                result_summary['trajectory_summary'] = {
                    'num_iterations': len(result['trajectory'].get('iterations', [])),
                    'has_discrete_sequences': 'discrete_sequences' in result['trajectory'],
                    'uniqueness_rate': uniqueness_rate if 'uniqueness_rate' in locals() else None
                }
            json.dump(result_summary, f, indent=2)
        

        
    else:

        if result:

    
    print("\n" + "=" * 80)

    print("=" * 80)
    
    if result and result.get('status') == 'completed':

        issues = []
        

        if 'uniqueness_rate' in locals() and uniqueness_rate < 10:

        

        if 'match_rate' in locals() and match_rate < 100:

        
        if issues:

            for issue in issues:
                print(f"  • {issue}")
        else:





if __name__ == "__main__":
    run_single_test()