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


logger = setup_logging(level='INFO', name='while_analysis')


def analyze_improvement_options(optimizer, indices, target_cai):

    
    current_cai = optimizer._compute_cai_from_indices(indices)

    

    possible_improvements = []
    
    for pos in range(len(indices)):
        if pos >= len(optimizer.codon_choices):
            continue
        
        current_idx = indices[pos]
        choices = optimizer.codon_choices[pos]
        

        current_choice = choices[current_idx] if current_idx < len(choices) else None
        if current_choice:



        

        for new_idx in range(len(choices)):
            if new_idx == current_idx:
                continue
            

            test_indices = indices.copy()
            test_indices[pos] = new_idx
            test_cai = optimizer._compute_cai_from_indices(test_indices)
            
            if test_cai > current_cai:
                improvement = {
                    'pos': pos,
                    'from_idx': current_idx,
                    'to_idx': new_idx,
                    'from_weight': choices[current_idx]['weight'] if current_idx < len(choices) else 0,
                    'to_weight': choices[new_idx]['weight'],
                    'cai_before': current_cai,
                    'cai_after': test_cai,
                    'cai_gain': test_cai - current_cai
                }
                possible_improvements.append(improvement)
    
    return possible_improvements


def test_while_loop_analysis():
    """分析while循环行为"""
    logger.info("\n" + "="*80)
    logger.info("分析v4 while循环为什么没有工作")
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
    

    torch.manual_seed(42)
    pi_test = torch.rand(seq_len, num_codons, device=device)
    pi_test = pi_test * valid_mask.float()
    for pos in range(seq_len):
        if valid_mask[pos].any():
            pi_test[pos] = pi_test[pos] / pi_test[pos].sum()
    

    gamma = 0.3
    target_cai = 0.8
    
    indices = optimizer._initial_optimization(pi_test, target_cai, gamma)
    current_cai = optimizer._compute_cai_from_indices(indices)
    
    logger.info(f"\ngamma={gamma}时的初始优化:")
    logger.info(f"CAI = {current_cai:.4f}")
    

    improvements = analyze_improvement_options(optimizer, indices, target_cai)
    
    logger.info(f"\n找到 {len(improvements)} 个可能的改进")
    
    if improvements:

        improvements.sort(key=lambda x: x['cai_gain'])
        
        logger.info("\n最小的5个改进（最接近目标）：")
        for i, imp in enumerate(improvements[:5]):
            logger.info(f"{i+1}. 位置{imp['pos']}: "
                       f"权重 {imp['from_weight']:.4f} → {imp['to_weight']:.4f}, "
                       f"CAI {imp['cai_before']:.4f} → {imp['cai_after']:.4f} "
                       f"(+{imp['cai_gain']:.4f})")
            
            if imp['cai_after'] >= target_cai:
                logger.info(f"   ✅ 这个改进会达到目标CAI")
        

        suitable_improvements = [imp for imp in improvements 
                                if imp['cai_after'] >= target_cai and 
                                   imp['cai_after'] < target_cai + 0.05]
        
        if suitable_improvements:
            logger.info(f"\n找到 {len(suitable_improvements)} 个接近目标的改进")
        else:
            logger.info("\n⚠️ 没有找到接近目标的改进")
            

            min_improvement = min(improvements, key=lambda x: x['cai_after'])
            logger.info(f"最小的改进会导致CAI = {min_improvement['cai_after']:.4f}")
            
            if min_improvement['cai_after'] > target_cai + 0.1:
                logger.info("问题：任何单个改变都会让CAI跳得太高")
                logger.info("这是离散优化的典型问题 - 无法精确控制CAI值")
    else:
        logger.info("⚠️ 没有找到任何可以提升CAI的改进")
        logger.info("可能原因：")
        logger.info("1. 当前选择已经是在gamma=0.3权重下的最优")
        logger.info("2. 密码子选择已经按权重排序，当前选择可能都是较优的")
    

    logger.info("\n检查密码子选择是否正确排序：")
    for pos in range(min(5, len(optimizer.codon_choices))):
        choices = optimizer.codon_choices[pos]
        weights = [c['weight'] for c in choices]
        is_sorted = all(weights[i] >= weights[i+1] for i in range(len(weights)-1))
        logger.info(f"位置{pos}: 权重={weights[:3]}... {'✅已排序' if is_sorted else '❌未排序'}")
    

    logger.info("\n" + "="*80)
    logger.info("分析总结")
    logger.info("="*80)
    logger.info("1. v4的while循环实际上也无法将CAI从0.7939提升到0.8")
    logger.info("2. 原因可能是：")
    logger.info("   - 任何单个密码子改变都会让CAI跳得太高（>0.92）")
    logger.info("   - 这是离散优化的固有限制")
    logger.info("3. 这解释了为什么v4也停在0.7939")
    logger.info("4. 解决方案：")
    logger.info("   - 接受CAI可能略低于目标（如0.79）")
    logger.info("   - 或微调gamma到0.35-0.36以达到目标")


def main():

    test_while_loop_analysis()


if __name__ == "__main__":
    main()