"""





"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
import hashlib
from id3.optimizers.cai.sado import SADOOptimizer as SADOOriginal
from id3.optimizers.cai.sado_improved import SADOOptimizer as SADOImproved
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_single_exp')


def test_single_experiment(optimizer_class, name):


    logger.info("="*60)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    


    optimizer = optimizer_class(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    sequences_exp1 = []
    

    torch.manual_seed(42)
    np.random.seed(42)
    


        pi_accessibility = torch.ones(seq_len, num_codons, device=device) * 0.5
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
            sequences_exp1.append(optimizer.last_indices.copy())
        
        if (i + 1) % 20 == 0:

    

    unique_hashes_exp1 = set()
    for seq in sequences_exp1:
        seq_hash = hashlib.md5(seq.tobytes()).hexdigest()
        unique_hashes_exp1.add(seq_hash)
    
    repetition_exp1 = 1.0 - (len(unique_hashes_exp1) / len(sequences_exp1))



    


    

    optimizer2 = optimizer_class(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    sequences_exp2 = []
    

    torch.manual_seed(42)
    np.random.seed(42)
    
    for i in range(100):

        pi_accessibility = torch.ones(seq_len, num_codons, device=device) * 0.5
        pi_accessibility = pi_accessibility * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        optimized_dist, metadata = optimizer2.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        
        if optimizer2.last_indices is not None:
            sequences_exp2.append(optimizer2.last_indices.copy())
    

    unique_hashes_exp2 = set()
    for seq in sequences_exp2:
        seq_hash = hashlib.md5(seq.tobytes()).hexdigest()
        unique_hashes_exp2.add(seq_hash)
    
    repetition_exp2 = 1.0 - (len(unique_hashes_exp2) / len(sequences_exp2))



    


    if sequences_exp1 and sequences_exp2:
        first_seq_exp1 = sequences_exp1[0]
        first_seq_exp2 = sequences_exp2[0]
        
        if np.array_equal(first_seq_exp1, first_seq_exp2):

        else:

    

    overlap = unique_hashes_exp1.intersection(unique_hashes_exp2)
    overlap_rate = len(overlap) / len(unique_hashes_exp1) if unique_hashes_exp1 else 0

    
    return repetition_exp1, repetition_exp2, overlap_rate


def main():
    """主测试函数"""
    
    logger.info("\n" + "="*80)
    logger.info("SADO单实验行为测试")
    logger.info("="*80)
    logger.info("验证SADO在单个实验内的重复情况和实验间的独立性")
    

    rep1_orig, rep2_orig, overlap_orig = test_single_experiment(SADOOriginal, "原始SADO")
    

    rep1_imp, rep2_imp, overlap_imp = test_single_experiment(SADOImproved, "改进版SADO（默认）")
    

    logger.info("\n测试 改进版SADO（无全局历史）")
    logger.info("="*60)
    

    SADOImproved.clear_all_global_history()
    

    class SADONoGlobal(SADOImproved):
        def __init__(self, *args, **kwargs):
            kwargs['use_global_history'] = False
            super().__init__(*args, **kwargs)
    
    rep1_no_global, rep2_no_global, overlap_no_global = test_single_experiment(
        SADONoGlobal, "改进版SADO（无全局历史）"
    )
    

    logger.info("\n" + "="*80)
    logger.info("测试总结")
    logger.info("="*80)
    
    logger.info("\n单个实验内重复率:")
    logger.info(f"  原始SADO: 实验1={rep1_orig*100:.1f}%, 实验2={rep2_orig*100:.1f}%")
    logger.info(f"  改进版（默认）: 实验1={rep1_imp*100:.1f}%, 实验2={rep2_imp*100:.1f}%")
    logger.info(f"  改进版（无全局）: 实验1={rep1_no_global*100:.1f}%, 实验2={rep2_no_global*100:.1f}%")
    
    logger.info("\n实验间序列重叠率:")
    logger.info(f"  原始SADO: {overlap_orig*100:.1f}%")
    logger.info(f"  改进版（默认）: {overlap_imp*100:.1f}%")
    logger.info(f"  改进版（无全局）: {overlap_no_global*100:.1f}%")
    
    logger.info("\n分析:")
    logger.info("• 原始SADO: 单实验内不重复，实验间独立")
    logger.info("• 改进版（默认）: 单实验内不重复，实验间也不重复（全局唯一）")
    logger.info("• 改进版（无全局）: 单实验内不重复，实验间独立")
    
    logger.info("\n建议:")
    logger.info("对于论文实验，应该使用无全局历史模式或原始SADO")
    logger.info("这样每个实验都是独立的，符合科学实验的可重现性要求")


if __name__ == "__main__":
    main()