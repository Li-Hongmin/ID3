#!/usr/bin/env python3
"""


"""

import torch
import sys
from pathlib import Path

# Add project path
sys.path.append(str(Path(__file__).parent))

from id3.utils.constraint_satisfied_argmax import get_constraint_satisfied_argmax
from id3.utils.constants import amino_acids_to_codons, NUCLEOTIDES

def rna_to_amino_acids(rna_seq):

    codon_to_aa = {}
    for aa, codons in amino_acids_to_codons.items():
        for codon in codons:

    
    aa_seq = []
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        if codon in codon_to_aa:
            aa_seq.append(codon_to_aa[codon])
        else:
            aa_seq.append('X')
    return ''.join(aa_seq)

def test_batch_independence():
    """测试批次独立性（修复前会失败）"""
    print("=" * 80)
    print("测试1：批次独立性")
    print("=" * 80)
    
    amino_acid_sequence = "MEEPQSD"
    batch_size = 5
    seq_len = len(amino_acid_sequence) * 3
    

    torch.manual_seed(42)
    rna_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)
    

    result = get_constraint_satisfied_argmax(rna_probs, amino_acid_sequence)
    

    sequences = []
    for b in range(batch_size):
        indices = torch.argmax(result[b], dim=-1)
        rna_seq = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
        sequences.append(rna_seq)
    

    unique_sequences = set(sequences)
    diversity_rate = len(unique_sequences) / len(sequences)
    
    print(f"批次大小: {batch_size}")
    print(f"唯一序列数: {len(unique_sequences)}")
    print(f"多样性率: {diversity_rate:.1%}")
    

    for i, seq in enumerate(sequences[:3]):
        print(f"  批次{i}: {seq}")
    

    if diversity_rate > 0.5:
        print("✅ 通过：批次独立选择密码子")
    else:
        print("❌ 失败：批次使用相同密码子（批平均问题）")
    
    return diversity_rate > 0.5

def test_constraint_satisfaction():

    print("\n" + "=" * 80)

    print("=" * 80)
    
    amino_acid_sequence = "MEEPQSDPSVEP"
    batch_size = 10
    seq_len = len(amino_acid_sequence) * 3
    

    torch.manual_seed(123)
    rna_probs = torch.softmax(torch.randn(batch_size, seq_len, 4), dim=-1)
    

    result = get_constraint_satisfied_argmax(rna_probs, amino_acid_sequence)
    

    all_match = True
    for b in range(batch_size):
        indices = torch.argmax(result[b], dim=-1)
        rna_seq = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
        aa_seq = rna_to_amino_acids(rna_seq)
        
        if aa_seq != amino_acid_sequence:

            all_match = False
    
    if all_match:

    else:

    
    return all_match

def test_deterministic_with_same_input():
    """测试相同输入的确定性"""
    print("\n" + "=" * 80)
    print("测试3：相同输入的确定性")
    print("=" * 80)
    
    amino_acid_sequence = "MGK"
    batch_size = 3
    seq_len = len(amino_acid_sequence) * 3
    

    torch.manual_seed(42)
    base_prob = torch.softmax(torch.randn(1, seq_len, 4), dim=-1)
    rna_probs = base_prob.repeat(batch_size, 1, 1)
    

    result = get_constraint_satisfied_argmax(rna_probs, amino_acid_sequence)
    

    sequences = []
    for b in range(batch_size):
        indices = torch.argmax(result[b], dim=-1)
        rna_seq = ''.join([NUCLEOTIDES[idx] for idx in indices.cpu().numpy()])
        sequences.append(rna_seq)
    

    if len(set(sequences)) == 1:
        print(f"✅ 通过：相同输入产生相同输出")
        print(f"  序列: {sequences[0]}")
    else:
        print(f"❌ 失败：相同输入产生不同输出")
    
    return len(set(sequences)) == 1

def main():
    print("\n" + "=" * 80)
    print("🔬 Lagrangian约束离散化修复验证")
    print("=" * 80)
    
    tests_passed = 0
    total_tests = 3
    

    if test_batch_independence():
        tests_passed += 1
    
    if test_constraint_satisfaction():
        tests_passed += 1
    
    if test_deterministic_with_same_input():
        tests_passed += 1
    

    print("\n" + "=" * 80)
    print("📊 测试总结")
    print("=" * 80)
    print(f"通过测试: {tests_passed}/{total_tests}")
    
    if tests_passed == total_tests:
        print("✅ 所有测试通过！离散化修复成功。")
        print("\n建议：")
        print("1. 重新运行Lagrangian实验（特别是variant 01）")
        print("2. 验证序列多样性是否改善")
        print("3. 检查最佳结果是否满足约束")
    else:
        print("❌ 某些测试失败，需要进一步调试。")
    
    return tests_passed == total_tests

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)