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


logger = setup_logging(level='INFO', name='test_optimal')


def find_optimal_gamma(optimizer, pi_accessibility, target_cai, 
                       amino_sequence, valid_mask, tolerance=0.02):

    
    left, right = 0.0, 1.0
    best_gamma = 0.5
    best_cai = 0.0
    best_diff = float('inf')
    

        mid = (left + right) / 2
        
        optimizer.reset()
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=mid
        )
        
        cai = metadata['final_cai']
        diff = abs(cai - target_cai)
        

        if diff < best_diff:
            best_gamma = mid
            best_cai = cai
            best_diff = diff
        

        if diff < tolerance:
            return best_gamma, best_cai
        

        if cai < target_cai:

            left = mid
        else:

            right = mid
    
    return best_gamma, best_cai


def test_gamma_optimization():
    """测试gamma优化"""
    logger.info("\n" + "="*80)
    logger.info("寻找最优gamma值以实现CAI≈0.8")
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
    final_cais = []
    
    for test_id in range(10):
        torch.manual_seed(200 + test_id)
        pi_test = torch.rand(seq_len, num_codons, device=device)
        pi_test = pi_test * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
        

        best_gamma, best_cai = find_optimal_gamma(
            optimizer, pi_test, 0.8, amino_sequence, valid_mask
        )
        
        optimal_gammas.append(best_gamma)
        final_cais.append(best_cai)
        
        status = "✅" if abs(best_cai - 0.8) < 0.02 else "⚠️"
        logger.info(f"分布{test_id+1}: gamma={best_gamma:.3f} → CAI={best_cai:.4f} {status}")
    

    logger.info("\n" + "="*80)
    logger.info("统计分析")
    logger.info("="*80)
    
    logger.info(f"最优gamma范围: {min(optimal_gammas):.3f} - {max(optimal_gammas):.3f}")
    logger.info(f"平均最优gamma: {np.mean(optimal_gammas):.3f} (±{np.std(optimal_gammas):.3f})")
    logger.info(f"CAI范围: {min(final_cais):.4f} - {max(final_cais):.4f}")
    logger.info(f"平均CAI: {np.mean(final_cais):.4f} (±{np.std(final_cais):.4f})")
    

    logger.info("\n固定分布下不同gamma的效果：")
    logger.info("-" * 60)
    
    torch.manual_seed(42)
    pi_fixed = torch.rand(seq_len, num_codons, device=device)
    pi_fixed = pi_fixed * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_fixed[pos] = pi_fixed[pos] / pi_fixed[pos].sum()
    
    test_gammas = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8]
    
    for gamma in test_gammas:
        optimizer.reset()
        _, metadata = optimizer.optimize(
            pi_accessibility=pi_fixed,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        
        cai = metadata['final_cai']
        distance = abs(cai - 0.8)
        
        if distance < 0.02:
            logger.info(f"gamma={gamma:.2f}: CAI={cai:.4f} ✅ (距离目标 {distance:.4f})")
        else:
            logger.info(f"gamma={gamma:.2f}: CAI={cai:.4f}    (距离目标 {distance:.4f})")
    

    logger.info("\n" + "="*80)
    logger.info("结论和建议")
    logger.info("="*80)
    
    logger.info("1. 当前SADO的gamma默认值0.6过高，导致CAI远超目标0.8")
    logger.info(f"2. 建议将默认gamma调整为 {np.mean(optimal_gammas):.2f} 左右")
    logger.info("3. 或者实现自适应gamma搜索，根据输入分布动态调整")
    logger.info("4. 最理想的是添加boundary search模式，确保CAI刚好满足约束")
    logger.info("")
    logger.info("这样可以实现 arg max P(S|π) s.t. CAI≥0.8 的真正目标：")
    logger.info("- 满足CAI约束（≥0.8）")
    logger.info("- 最大化条件概率（不过度优化CAI）")


def main():

    test_gamma_optimization()


if __name__ == "__main__":
    main()