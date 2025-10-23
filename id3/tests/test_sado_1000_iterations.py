"""



"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
import hashlib
import time
from id3.optimizers.cai.sado import SADOOptimizer
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_1000')


def test_sado_1000_iterations():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    

    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
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
    

    all_sequences = []
    all_hashes = []
    cai_values = []
    

    start_time = time.time()
    

    logger.info("-" * 60)
    

    for i in range(1000):


        torch.manual_seed(42 + i)
        np.random.seed(42 + i)
        

        pi_accessibility = torch.rand(seq_len, num_codons, device=device)
        pi_accessibility = pi_accessibility * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        

        if optimizer.last_indices is not None:
            seq_copy = optimizer.last_indices.copy()
            all_sequences.append(seq_copy)
            

            seq_hash = hashlib.md5(seq_copy.tobytes()).hexdigest()
            all_hashes.append(seq_hash)
            

            cai_values.append(metadata['final_cai'])
        

        if (i + 1) % 100 == 0:
            elapsed = time.time() - start_time
            unique_so_far = len(set(all_hashes[:i+1]))
            repetition_so_far = 1.0 - (unique_so_far / (i + 1))
            avg_cai = np.mean(cai_values[:i+1])
            





    

    total_time = time.time() - start_time
    unique_hashes_set = set(all_hashes)
    total_sequences = len(all_sequences)
    unique_sequences = len(unique_hashes_set)
    repetition_rate = 1.0 - (unique_sequences / total_sequences)
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    







    





    

    if repetition_rate > 0:

        hash_counts = {}
        for h in all_hashes:
            hash_counts[h] = hash_counts.get(h, 0) + 1
        
        repeated_hashes = {h: count for h, count in hash_counts.items() if count > 1}

        

        if repeated_hashes:
            max_repeat = max(repeated_hashes.values())

            

            sorted_repeats = sorted(repeated_hashes.items(), key=lambda x: x[1], reverse=True)

            for i, (h, count) in enumerate(sorted_repeats[:5]):

                first_idx = all_hashes.index(h)
                last_idx = len(all_hashes) - 1 - all_hashes[::-1].index(h)

    




    

    if hasattr(optimizer, 'get_diversity_stats'):
        diversity_stats = optimizer.get_diversity_stats()



    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    if repetition_rate == 0:





    else:

    
    return repetition_rate


def main():
    """主函数"""
    repetition_rate = test_sado_1000_iterations()
    

    if repetition_rate == 0:
        logger.info("\n🎉 测试通过：SADO在1000次迭代中保持了完全的序列唯一性")
        return 0
    else:
        logger.warning(f"\n⚠️ 测试发现重复：{repetition_rate*100:.2f}% 的序列重复率")
        return 1


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)