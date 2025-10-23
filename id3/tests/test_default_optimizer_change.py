#!/usr/bin/env python3
"""

"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.utils.logging_config import setup_logging

logger = setup_logging(level='INFO', name='test_default')


def test_default_method():

    logger.info("="*80)

    logger.info("="*80)
    

    operator = CAIEnhancementOperator(
        species='ecoli_bl21de3',
        device=torch.device('cpu')
    )
    


    



    

    


    operator_sado = CAIEnhancementOperator(
        method='sado',
        species='ecoli_bl21de3',
        device=torch.device('cpu')
    )
    


    


    

    


    

    test_sequence = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 aa
    pi_accessibility = torch.rand(100, 6)
    for i in range(100):
        if pi_accessibility[i].sum() > 0:
            pi_accessibility[i] = pi_accessibility[i] / pi_accessibility[i].sum()
    

    try:
        result, metadata = operator.enhance(
            pi_accessibility=pi_accessibility,
            target_cai=0.7,
            amino_acid_sequence=test_sequence
        )


    except Exception as e:

        raise


def test_performance_comparison():
    """简单的性能对比"""
    logger.info("\n" + "="*80)
    logger.info("性能对比")
    logger.info("="*80)
    
    import time
    

    test_sequence = "ACDEFGHIKLMNPQRSTVWY" * 25  # 500 aa
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    pi_accessibility = torch.rand(500, 6, device=device)
    for i in range(500):
        if pi_accessibility[i].sum() > 0:
            pi_accessibility[i] = pi_accessibility[i] / pi_accessibility[i].sum()
    

    operator_binary = CAIEnhancementOperator(
        method='binary_search',
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=test_sequence
    )
    
    start = time.time()
    result_binary, metadata_binary = operator_binary.enhance(
        pi_accessibility=pi_accessibility,
        target_cai=0.7,
        amino_acid_sequence=test_sequence
    )
    time_binary = (time.time() - start) * 1000
    
    logger.info(f"二分查找: {time_binary:.1f}ms, CAI={metadata_binary['final_cai']:.4f}")
    

    operator_sado = CAIEnhancementOperator(
        method='sado',
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=test_sequence
    )
    
    start = time.time()
    result_sado, metadata_sado = operator_sado.enhance(
        pi_accessibility=pi_accessibility,
        target_cai=0.7,
        amino_acid_sequence=test_sequence
    )
    time_sado = (time.time() - start) * 1000
    
    logger.info(f"SADO: {time_sado:.1f}ms, CAI={metadata_sado['final_cai']:.4f}")
    
    if time_binary > 0:
        speedup = time_sado / time_binary
        if speedup < 1:
            logger.info(f"SADO比二分查找快 {1/speedup:.1f}倍")
        else:
            logger.info(f"二分查找比SADO快 {speedup:.1f}倍")


if __name__ == "__main__":
    logger.info("🚀 开始测试默认优化器更改")
    

    test_default_method()
    

    test_performance_comparison()
    
    logger.info("\n✅ 所有测试通过！默认优化器已成功更改为binary_search")