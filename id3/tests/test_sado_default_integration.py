"""







"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.optimizers.cai import DefaultCAIOptimizer, SADOOptimizerImproved
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_integration')


def test_default_optimizer():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    



    

    optimizer = DefaultCAIOptimizer(
        species='ecoli_bl21de3',
        device=torch.device('cpu')
    )
    


    
    return DefaultCAIOptimizer == SADOOptimizerImproved


def test_cai_enhancement_operator():
    """测试CAI增强操作符的默认方法"""
    logger.info("\n" + "="*80)
    logger.info("测试CAI增强操作符")
    logger.info("="*80)
    

    operator = CAIEnhancementOperator()
    
    logger.info(f"默认方法: {operator.method}")
    logger.info(f"优化器类型: {type(operator.optimizer).__name__}")
    

    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    seq_len = len(amino_sequence)
    num_codons = 6
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    pi_accessibility = torch.rand(seq_len, num_codons, device=device)
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    

    pi_accessibility = pi_accessibility * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
    

    discrete_dist, metadata = operator.enhance(
        pi_accessibility=pi_accessibility,
        target_cai=0.8,
        amino_acid_sequence=amino_sequence,
        valid_codon_mask=valid_mask
    )
    
    final_cai = metadata['final_cai']
    method = metadata['method']
    
    logger.info(f"使用方法: {method}")
    logger.info(f"达到的CAI: {final_cai:.4f}")
    logger.info(f"满足目标CAI (≥0.8): {final_cai >= 0.8}")
    
    return method == 'sado' and final_cai >= 0.8


def test_sado_performance_batch():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    operator = CAIEnhancementOperator(method='sado')
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    

    num_tests = 30
    cai_values = []
    all_sequences = []
    

    
    for i in range(num_tests):

        pi_accessibility = torch.rand(seq_len, num_codons, device=device)
        pi_accessibility = pi_accessibility * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
        

        discrete_dist, metadata = operator.enhance(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        
        cai_values.append(metadata['final_cai'])
        

        if hasattr(operator.optimizer, 'last_indices'):
            all_sequences.append(operator.optimizer.last_indices.copy())
        
        if (i + 1) % 10 == 0:

    






    

    satisfied = sum(1 for cai in cai_values if cai >= 0.8)

    

    if all_sequences:
        unique_hashes = set()
        for seq in all_sequences:
            seq_hash = operator.optimizer._hash_sequence(seq)
            unique_hashes.add(seq_hash)
        
        repetition_rate = 1.0 - (len(unique_hashes) / len(all_sequences))

    

    if hasattr(operator.optimizer, 'get_diversity_stats'):
        diversity_stats = operator.optimizer.get_diversity_stats()

    
    success = np.mean(cai_values) >= 0.8 and repetition_rate == 0.0
    
    if success:

    else:

    
    return success


def test_multiple_sessions():
    """测试跨会话的序列唯一性"""
    logger.info("\n" + "="*80)
    logger.info("测试跨会话序列唯一性")
    logger.info("="*80)
    
    amino_sequence = "MKAI"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    seq_len = len(amino_sequence)
    num_codons = 6
    

    valid_mask = torch.zeros(seq_len, num_codons, dtype=torch.bool, device=device)
    for pos, aa in enumerate(amino_sequence):
        if aa in amino_acids_to_codons:
            num_valid = len(amino_acids_to_codons[aa])
            valid_mask[pos, :min(num_valid, num_codons)] = True
    
    all_sequences = []
    

    num_sessions = 3
    sequences_per_session = 5
    
    for session in range(num_sessions):
        logger.info(f"\n会话 {session + 1}/{num_sessions}")
        

        operator = CAIEnhancementOperator()
        
        for i in range(sequences_per_session):

            pi_accessibility = torch.rand(seq_len, num_codons, device=device)
            pi_accessibility = pi_accessibility * valid_mask.float()
            for pos in range(seq_len):
                if valid_mask[pos].any():
                    pi_accessibility[pos] = pi_accessibility[pos] / pi_accessibility[pos].sum()
            

            discrete_dist, metadata = operator.enhance(
                pi_accessibility=pi_accessibility,
                target_cai=0.8,
                amino_acid_sequence=amino_sequence,
                valid_codon_mask=valid_mask
            )
            

            if hasattr(operator.optimizer, 'last_indices'):
                all_sequences.append(operator.optimizer.last_indices.copy())
        
        logger.info(f"  会话生成: {sequences_per_session} 序列")
    

    unique_hashes = set()
    for seq in all_sequences:

        import hashlib
        seq_hash = hashlib.md5(seq.tobytes()).hexdigest()
        unique_hashes.add(seq_hash)
    
    total_sequences = len(all_sequences)
    unique_sequences = len(unique_hashes)
    global_repetition = 1.0 - (unique_sequences / total_sequences)
    
    logger.info("\n全局统计:")
    logger.info(f"  总序列数: {total_sequences}")
    logger.info(f"  唯一序列数: {unique_sequences}")
    logger.info(f"  全局重复率: {global_repetition*100:.1f}%")
    
    if global_repetition == 0:
        logger.info("\n✅ 跨会话唯一性测试通过")
    else:
        logger.warning(f"\n⚠️ 跨会话有 {global_repetition*100:.1f}% 重复")
    
    return global_repetition == 0


def main():

    logger.info("\n" + "="*80)

    logger.info("="*80)

    
    tests = [




    ]
    
    results = []
    for name, test_func in tests:
        try:

            success = test_func()
            results.append((name, success))
        except Exception as e:

            import traceback
            traceback.print_exc()
            results.append((name, False))
    

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:

        logger.info(f"{name:20} {status}")
    

    
    if passed == total:







    else:

    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)