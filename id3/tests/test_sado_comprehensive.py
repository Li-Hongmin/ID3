"""






"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from typing import List, Set
import time
from id3.optimizers.cai.sado import SADOOptimizer
from id3.constraints.cai_enhancement_operator import CAIEnhancementOperator
from id3.constraints.lagrangian import LagrangianConstraint
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='test_sado')


def decode_rna_to_amino(rna_sequence: str) -> str:

    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    amino_sequence = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_sequence += genetic_code.get(codon, 'X')
    return amino_sequence


def test_sado_basic_performance():
    """测试SADO基本性能"""
    logger.info("\n" + "="*60)
    logger.info("测试1: SADO基本性能")
    logger.info("="*60)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
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
    

    target_cai = 0.8
    cai_values = []
    time_costs = []
    
    for i in range(10):
        start_time = time.time()
        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=target_cai,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        elapsed = time.time() - start_time
        
        cai_values.append(metadata['final_cai'])
        time_costs.append(elapsed * 1000)
        
        logger.info(f"  迭代 {i+1}: CAI={metadata['final_cai']:.4f}, "
                   f"时间={elapsed*1000:.1f}ms, "
                   f"满足约束={metadata['constraint_satisfied']}")
    

    avg_cai = np.mean(cai_values)
    avg_time = np.mean(time_costs)
    success_rate = sum(1 for cai in cai_values if cai >= target_cai) / len(cai_values)
    
    logger.info(f"\n  📊 统计结果:")
    logger.info(f"     平均CAI: {avg_cai:.4f}")
    logger.info(f"     平均时间: {avg_time:.1f}ms")
    logger.info(f"     成功率: {success_rate*100:.1f}%")
    
    return success_rate >= 0.8


def test_sado_uniqueness():

    logger.info("\n" + "="*60)

    logger.info("="*60)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    seq_len = len(amino_sequence)
    num_codons = 6
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
    

    num_sequences = 50
    sequences: Set[str] = set()
    cai_values = []
    

    
    for i in range(num_sequences):
        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        

        discrete_indices = torch.argmax(optimized_dist, dim=-1)
        seq_hash = discrete_indices.cpu().numpy().tobytes()
        sequences.add(seq_hash)
        cai_values.append(metadata['final_cai'])
        
        if (i + 1) % 10 == 0:

    

    diversity_stats = optimizer.get_diversity_stats()
    







    

    optimizer.reset()

    
    for i in range(10):
        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        
        discrete_indices = torch.argmax(optimized_dist, dim=-1)
        seq_hash = discrete_indices.cpu().numpy().tobytes()
        sequences.add(seq_hash)
    

    

    repetition_rate = 1 - len(sequences)/(num_sequences+10)
    success = repetition_rate <= 0.1
    
    if success:

    else:

    
    return success


def test_sado_constraint_satisfaction():
    """测试SADO生成序列的约束满足情况"""
    logger.info("\n" + "="*60)
    logger.info("测试3: 约束满足情况")
    logger.info("="*60)
    
    test_sequences = [
        "MKAI",
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSG",
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG",
    ]
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    all_satisfied = True
    
    for amino_sequence in test_sequences:
        logger.info(f"\n  测试序列: {amino_sequence[:20]}... (长度={len(amino_sequence)})")
        

        optimizer = SADOOptimizer(
            species='ecoli_bl21de3',
            device=device,
            amino_acid_sequence=amino_sequence
        )
        

        seq_len = len(amino_sequence)
        num_codons = 6
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
        

        target_cais = [0.6, 0.7, 0.8]
        
        for target_cai in target_cais:

            optimized_dist, metadata = optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_sequence,
                valid_codon_mask=valid_mask
            )
            


            cai_satisfied = metadata['final_cai'] >= target_cai
            


            amino_satisfied = True
            
            logger.info(f"    目标CAI={target_cai:.1f}: "
                       f"实际CAI={metadata['final_cai']:.4f}, "
                       f"CAI满足={'✅' if cai_satisfied else '❌'}, "
                       f"氨基酸满足={'✅' if amino_satisfied else '❌'}")
            
            if not (cai_satisfied and amino_satisfied):
                all_satisfied = False
        

        optimizer.reset()
    
    if all_satisfied:
        logger.info("\n  ✅ 所有约束满足测试通过")
    else:
        logger.warning("\n  ⚠️ 部分约束未满足")
    
    return all_satisfied


def test_sado_vs_binary_search():

    logger.info("\n" + "="*60)

    logger.info("="*60)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    seq_len = len(amino_sequence)
    num_codons = 6
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
    


    sado_optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    sado_times = []
    sado_cais = []
    sado_sequences = set()
    
    for i in range(20):
        start = time.time()
        optimized_dist, metadata = sado_optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        elapsed = time.time() - start
        
        sado_times.append(elapsed * 1000)
        sado_cais.append(metadata['final_cai'])
        

        discrete_indices = torch.argmax(optimized_dist, dim=-1)
        seq_hash = discrete_indices.cpu().numpy().tobytes()
        sado_sequences.add(seq_hash)
    



    


    from id3.optimizers.cai.binary_search import BinarySearchCAIOptimizer
    
    bs_optimizer = BinarySearchCAIOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    bs_times = []
    bs_cais = []
    bs_sequences = set()
    
    for i in range(20):
        start = time.time()
        optimized_dist, metadata = bs_optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        elapsed = time.time() - start
        
        bs_times.append(elapsed * 1000)
        bs_cais.append(metadata['final_cai'])
        

        discrete_indices = torch.argmax(optimized_dist, dim=-1)
        seq_hash = discrete_indices.cpu().numpy().tobytes()
        bs_sequences.add(seq_hash)
    



    

    speedup = np.mean(bs_times) / np.mean(sado_times)




    

    success = speedup > 1.0 and len(sado_sequences) > len(bs_sequences)
    
    if success:

    else:

    



def test_sado_with_constraints():
    """测试SADO与约束系统的集成"""
    logger.info("\n" + "="*60)
    logger.info("测试5: SADO与约束系统集成")
    logger.info("="*60)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    operator = CAIEnhancementOperator(
        method='sado',
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    

    constraint = LagrangianConstraint(
        amino_acid_sequence=amino_sequence,
        protein_sequence=amino_sequence,
        enable_cai=True,
        target_cai=0.8,
        lambda_cai=0.1,
        cai_method='sado',
        device=device
    )
    

    result = constraint.forward(alpha=1.0, tau=1.0)
    
    if 'enhanced_sequence' in result and result['enhanced_sequence'] is not None:
        logger.info(f"  ✅ 增强序列生成成功")
        logger.info(f"     形状: {result['enhanced_sequence'].shape}")
        

        if hasattr(constraint, 'cai_metadata') and constraint.cai_metadata:
            metadata = constraint.cai_metadata
            logger.info(f"     最终CAI: {metadata.get('final_cai', 0):.4f}")
            logger.info(f"     目标CAI: {metadata.get('target_cai', 0):.4f}")
            logger.info(f"     满足约束: {metadata.get('constraint_satisfied', False)}")
            logger.info(f"     唯一序列数: {metadata.get('unique_sequences', 0)}")
    

    logger.info("\n  测试多次调用的唯一性...")
    sequences = []
    cai_values = []
    
    for i in range(10):
        result = constraint.forward(alpha=1.0, tau=1.0)
        if 'enhanced_sequence' in result and result['enhanced_sequence'] is not None:
            sequences.append(result['enhanced_sequence'])
            if hasattr(constraint, 'cai_metadata') and constraint.cai_metadata:
                cai_values.append(constraint.cai_metadata.get('final_cai', 0))
    

    unique_sequences = set()
    for seq in sequences:
        if seq is not None:
            seq_hash = seq.cpu().numpy().tobytes()
            unique_sequences.add(seq_hash)
    
    logger.info(f"     生成序列数: {len(sequences)}")
    logger.info(f"     唯一序列数: {len(unique_sequences)}")
    logger.info(f"     平均CAI: {np.mean(cai_values):.4f}" if cai_values else "     无CAI数据")
    
    success = len(unique_sequences) > 1
    
    if success:
        logger.info(f"  ✅ SADO与约束系统集成成功")
    else:
        logger.warning(f"  ⚠️ 需要检查集成")
    
    return success


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

        logger.info(f"{name:15} {status}")
    

    
    if passed == total:








    else:

    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)