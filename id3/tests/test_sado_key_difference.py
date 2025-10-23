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


logger = setup_logging(level='INFO', name='test_key_diff')


def simulate_v4_behavior(optimizer, pi_accessibility, target_cai, gamma):

    

    indices = optimizer._initial_optimization(
        pi_accessibility, target_cai, gamma
    )
    

    current_cai = optimizer._compute_cai_from_indices(indices)
    

    

    iterations = 0
    while current_cai < target_cai and iterations < 100:
        improvements = []
        

        for pos in range(len(indices)):
            if pos >= len(optimizer.codon_choices):
                continue
            
            current_idx = indices[pos]
            choices = optimizer.codon_choices[pos]
            

            for new_idx, choice in enumerate(choices):
                if new_idx <= current_idx:
                    continue
                

                test_indices = indices.copy()
                test_indices[pos] = new_idx
                test_cai = optimizer._compute_cai_from_indices(test_indices)
                
                if test_cai > current_cai:
                    cai_gain = test_cai - current_cai
                    improvements.append((cai_gain, pos, new_idx, test_cai))
        
        if not improvements:
            break
        


        

        selected = None
        for gain, pos, new_idx, new_cai in improvements:
            if new_cai >= target_cai:
                selected = (gain, pos, new_idx, new_cai)
                break
        

        if selected is None and improvements:
            selected = improvements[0]
        
        if selected:
            _, pos, new_idx, new_cai = selected
            indices[pos] = new_idx
            current_cai = new_cai
            iterations += 1
        else:
            break
    

    return indices, current_cai


def test_key_difference():
    """测试关键差异"""
    logger.info("\n" + "="*80)
    logger.info("分析当前SADO vs v4的关键差异")
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
    
    logger.info("\n关键差异分析：")
    logger.info("-" * 60)
    logger.info("当前SADO: 一步优化，根据gamma平衡CAI和概率")
    logger.info("v4版本: 两步优化")
    logger.info("  1. 初始优化（类似当前版本）")
    logger.info("  2. while循环逐步提升CAI直到刚好满足约束")
    

    logger.info("\n测试不同概率分布：")
    logger.info("-" * 60)
    
    for test_id in range(5):
        logger.info(f"\n测试 {test_id + 1}:")
        
        torch.manual_seed(100 + test_id)
        pi_test = torch.rand(seq_len, num_codons, device=device)
        pi_test = pi_test * valid_mask.float()
        for pos in range(seq_len):
            if valid_mask[pos].any():
                pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
        

        gamma = 0.3
        target_cai = 0.8
        
        logger.info(f"gamma={gamma}, target_cai={target_cai}")
        

        optimizer.reset()
        dist_current, meta_current = optimizer.optimize(
            pi_accessibility=pi_test,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask,
            gamma=gamma
        )
        logger.info(f"  当前SADO CAI: {meta_current['final_cai']:.4f}")
        

        optimizer.reset()
        indices_v4, cai_v4 = simulate_v4_behavior(
            optimizer, pi_test, target_cai, gamma
        )
        

        diff = meta_current['final_cai'] - cai_v4
        if diff > 0.05:
            logger.info(f"  ⚠️ 当前版本CAI高出 {diff:.3f}")
    

    logger.info("\n" + "="*80)
    logger.info("建议的改进方案")
    logger.info("="*80)
    logger.info("为了实现 arg max P(S|π) s.t. CAI≥0.8 的目标：")
    logger.info("")
    logger.info("方案1: 添加boundary_search模式（推荐）")
    logger.info("  - 在SADO中添加类似v4的while循环")
    logger.info("  - 初始优化后，逐步微调直到刚好满足CAI约束")
    logger.info("  - 这样可以最大化条件概率")
    logger.info("")
    logger.info("方案2: 调整默认gamma值")
    logger.info("  - 将默认gamma从0.6降低到0.3左右")
    logger.info("  - 但这只是近似解决方案")
    logger.info("")
    logger.info("方案3: 自适应gamma搜索")
    logger.info("  - 二分搜索找到刚好满足约束的gamma值")
    logger.info("  - 但需要多次优化，效率较低")


def main():

    test_key_difference()


if __name__ == "__main__":
    main()