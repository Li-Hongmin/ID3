"""







"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))

import torch
import numpy as np
from id3.optimizers.cai.sado import SADOOptimizer
from id3.cai.unified_calculator import UnifiedCAICalculator
from id3.utils.constants import amino_acids_to_codons
from id3.utils.logging_config import setup_logging


logger = setup_logging(level='INFO', name='verify_sado')


def onehot_to_rna_sequence(onehot_sequence: torch.Tensor, amino_acid_sequence: str) -> str:
    """

    
    Args:


        
    Returns:

    """
    rna_sequence = ""
    
    for pos, aa in enumerate(amino_acid_sequence):
        if aa in amino_acids_to_codons:
            codons = amino_acids_to_codons[aa]
            

            if pos < len(onehot_sequence):
                selected_idx = torch.argmax(onehot_sequence[pos]).item()
                

                if selected_idx < len(codons):
                    selected_codon = codons[selected_idx]
                    rna_sequence += selected_codon
                else:

                    rna_sequence += codons[0]
                    logger.warning(f"位置 {pos}: 索引 {selected_idx} 超出密码子范围，使用默认密码子")
            else:

                rna_sequence += codons[0]
                logger.warning(f"位置 {pos} 超出序列范围")
        else:
            logger.warning(f"未知氨基酸: {aa}")
            rna_sequence += "NNN"
    
    return rna_sequence


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


def test_sado_sequence_verification():
    """测试SADO序列转换和CAI验证"""
    logger.info("\n" + "="*80)
    logger.info("SADO序列转换和CAI验证测试")
    logger.info("="*80)
    

    test_sequences = [
        "MKAI",
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSG",
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG",
    ]
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    cai_calculator = UnifiedCAICalculator(species='ecoli_bl21de3', device=device)
    
    all_matches = True
    
    for amino_sequence in test_sequences:
        logger.info(f"\n测试序列: {amino_sequence[:20]}... (长度={len(amino_sequence)})")
        logger.info("-" * 60)
        

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
        

        target_cais = [0.6, 0.7, 0.8, 0.9]
        
        for target_cai in target_cais:
            logger.info(f"\n  目标CAI: {target_cai}")
            

            optimized_dist, metadata = optimizer.optimize(
                pi_accessibility=pi_accessibility,
                target_cai=target_cai,
                amino_acid_sequence=amino_sequence,
                valid_codon_mask=valid_mask
            )
            
            sado_reported_cai = metadata['final_cai']
            logger.info(f"  1. SADO报告的CAI: {sado_reported_cai:.6f}")
            

            if optimized_dist.dim() == 3:
                optimized_dist = optimized_dist.squeeze(0)
            
            rna_sequence = onehot_to_rna_sequence(optimized_dist, amino_sequence)
            logger.info(f"  2. 生成的RNA序列: {rna_sequence[:30]}...")
            

            decoded_amino = decode_rna_to_amino(rna_sequence)
            logger.info(f"  3. 解码的氨基酸序列: {decoded_amino[:20]}...")
            

            amino_match = decoded_amino == amino_sequence
            if amino_match:
                logger.info(f"  ✅ 氨基酸序列匹配")
            else:
                logger.error(f"  ❌ 氨基酸序列不匹配!")
                logger.error(f"     原始: {amino_sequence}")
                logger.error(f"     解码: {decoded_amino}")
                all_matches = False
            

            actual_cai = cai_calculator.compute_cai(rna_sequence, method='standard')
            logger.info(f"  4. RNA序列计算的CAI: {actual_cai:.6f}")
            

            cai_diff = abs(sado_reported_cai - actual_cai)
            logger.info(f"  5. CAI差异: {cai_diff:.6f}")
            

            if cai_diff < 0.01:
                logger.info(f"  ✅ CAI值一致（差异 < 1%）")
            else:
                logger.warning(f"  ⚠️ CAI值有差异（差异 = {cai_diff:.4f}）")
                if cai_diff > 0.05:
                    all_matches = False
            

            if actual_cai >= target_cai:
                logger.info(f"  ✅ 满足目标CAI（{actual_cai:.4f} >= {target_cai}）")
            else:
                logger.warning(f"  ⚠️ 未达到目标CAI（{actual_cai:.4f} < {target_cai}）")
        

        optimizer.reset()
    
    return all_matches


def test_sado_batch_verification():

    logger.info("\n" + "="*80)

    logger.info("="*80)
    
    amino_sequence = "MSKGEELFTGVVPILVELDGDVNGHKFSVSG"
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    

    optimizer = SADOOptimizer(
        species='ecoli_bl21de3',
        device=device,
        amino_acid_sequence=amino_sequence
    )
    
    cai_calculator = UnifiedCAICalculator(species='ecoli_bl21de3', device=device)
    

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
    

    num_tests = 20
    cai_diffs = []
    amino_matches = 0
    

    
    for i in range(num_tests):

        optimized_dist, metadata = optimizer.optimize(
            pi_accessibility=pi_accessibility,
            target_cai=0.8,
            amino_acid_sequence=amino_sequence,
            valid_codon_mask=valid_mask
        )
        
        sado_cai = metadata['final_cai']
        

        if optimized_dist.dim() == 3:
            optimized_dist = optimized_dist.squeeze(0)
        
        rna_sequence = onehot_to_rna_sequence(optimized_dist, amino_sequence)
        

        decoded_amino = decode_rna_to_amino(rna_sequence)
        if decoded_amino == amino_sequence:
            amino_matches += 1
        

        actual_cai = cai_calculator.compute_cai(rna_sequence, method='standard')
        

        diff = abs(sado_cai - actual_cai)
        cai_diffs.append(diff)
        
        if (i + 1) % 5 == 0:

    







    

    success = (amino_matches == num_tests) and (np.mean(cai_diffs) < 0.01)
    
    if success:

    else:

    
    return success


def main():
    """主测试函数"""
    logger.info("\n" + "="*80)
    logger.info("SADO序列验证测试套件")
    logger.info("="*80)
    logger.info("验证one-hot序列转换为RNA后的CAI计算一致性")
    
    tests = [
        ("序列转换和CAI验证", test_sado_sequence_verification),
        ("批量CAI准确性验证", test_sado_batch_verification),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            logger.error(f"{name}测试失败: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    

    logger.info("\n" + "="*80)
    logger.info("测试总结")
    logger.info("="*80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:
        status = "✅ 通过" if success else "❌ 失败"
        logger.info(f"{name:25} {status}")
    
    logger.info(f"\n总计: {passed}/{total} 测试通过")
    
    if passed == total:
        logger.info("\n🎉 所有验证通过！")
        logger.info("\n验证结果：")
        logger.info("1. ✅ One-hot序列可以正确转换为RNA序列")
        logger.info("2. ✅ RNA序列正确编码目标氨基酸")
        logger.info("3. ✅ 转换后的CAI值与SADO报告值一致")
        logger.info("4. ✅ SADO的CAI计算是准确的")
    else:
        logger.warning("\n⚠️ 部分验证失败，需要检查")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)