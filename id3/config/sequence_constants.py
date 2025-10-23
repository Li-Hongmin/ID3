#!/usr/bin/env python3
"""
序列常量定义模块

定义所有序列相关的常量，包括密码子、UTR序列等
通过UTR加载器动态获取UTR模板，避免硬编码
"""

from typing import Dict, Set, Optional
from id3.config.utr_loader import get_utr_loader


# =============================================================================
# 密码子常量
# =============================================================================

# 起始密码子
START_CODON = 'ATG'
START_CODON_RNA = 'AUG'

# 终止密码子 (DNA)
STOP_CODONS_DNA = {'TAA', 'TAG', 'TGA'}

# 终止密码子 (RNA)  
STOP_CODONS_RNA = {'UAA', 'UAG', 'UGA'}

# 非同义密码子（只有一个密码子编码的氨基酸）
NON_SYNONYMOUS_CODONS = {'ATG', 'TGG'}  # Met, Trp

# 标准遗传密码表 (DNA)
GENETIC_CODE_DNA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

# 标准遗传密码表 (RNA)
GENETIC_CODE_RNA = {
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


# =============================================================================
# UTR序列获取函数
# =============================================================================

def get_default_utr5() -> str:
    """获取默认5'UTR序列
    
    从UTR模板文件动态加载，避免硬编码
    
    Returns:
        str: 5'UTR序列
    """
    loader = get_utr_loader()
    return loader.get_default_utr5()


def get_default_utr3() -> str:
    """获取默认3'UTR序列
    
    从UTR模板文件动态加载，避免硬编码
    
    Returns:
        str: 3'UTR序列
    """
    loader = get_utr_loader()
    return loader.get_default_utr3()


def get_utr_sequences() -> Dict[str, str]:
    """获取所有默认UTR序列
    
    Returns:
        Dict[str, str]: 包含utr5和utr3的字典
    """
    loader = get_utr_loader()
    return {
        'utr5': loader.get_default_utr5(),
        'utr3': loader.get_default_utr3()
    }


def get_utr_lengths() -> Dict[str, int]:
    """获取UTR序列长度信息
    
    Returns:
        Dict[str, int]: UTR长度信息
    """
    loader = get_utr_loader()
    info = loader.get_utr_info()
    return {
        'utr5_length': info['utr5_length'],
        'utr3_length': info['utr3_length'],
        'total_utr_length': info['total_utr_length']
    }


# =============================================================================
# 序列工具函数
# =============================================================================

def dna_to_rna(sequence: str) -> str:
    """将DNA序列转换为RNA序列
    
    Args:
        sequence: DNA序列
        
    Returns:
        str: RNA序列
    """
    return sequence.replace('T', 'U')


def rna_to_dna(sequence: str) -> str:
    """将RNA序列转换为DNA序列
    
    Args:
        sequence: RNA序列
        
    Returns:
        str: DNA序列
    """
    return sequence.replace('U', 'T')


def translate_dna(sequence: str) -> str:
    """翻译DNA序列为氨基酸序列
    
    Args:
        sequence: DNA序列
        
    Returns:
        str: 氨基酸序列
    """
    if len(sequence) % 3 != 0:
        raise ValueError("DNA序列长度必须是3的倍数")
    
    amino_acids = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE_DNA.get(codon, 'X')  # X表示未知氨基酸
        amino_acids.append(amino_acid)
    
    return ''.join(amino_acids)


def translate_rna(sequence: str) -> str:
    """翻译RNA序列为氨基酸序列
    
    Args:
        sequence: RNA序列
        
    Returns:
        str: 氨基酸序列
    """
    if len(sequence) % 3 != 0:
        raise ValueError("RNA序列长度必须是3的倍数")
    
    amino_acids = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE_RNA.get(codon, 'X')  # X表示未知氨基酸
        amino_acids.append(amino_acid)
    
    return ''.join(amino_acids)


def is_start_codon(codon: str, rna: bool = False) -> bool:
    """检查是否为起始密码子
    
    Args:
        codon: 密码子序列
        rna: 是否为RNA序列
        
    Returns:
        bool: 是否为起始密码子
    """
    if rna:
        return codon.upper() == START_CODON_RNA
    else:
        return codon.upper() == START_CODON


def is_stop_codon(codon: str, rna: bool = False) -> bool:
    """检查是否为终止密码子
    
    Args:
        codon: 密码子序列
        rna: 是否为RNA序列
        
    Returns:
        bool: 是否为终止密码子
    """
    codon = codon.upper()
    if rna:
        return codon in STOP_CODONS_RNA
    else:
        return codon in STOP_CODONS_DNA


def validate_nucleotide_sequence(sequence: str, rna: bool = False) -> bool:
    """验证核苷酸序列的有效性
    
    Args:
        sequence: 核苷酸序列
        rna: 是否为RNA序列
        
    Returns:
        bool: 序列是否有效
    """
    sequence = sequence.upper()
    if rna:
        valid_nucleotides = set('AUCG')
    else:
        valid_nucleotides = set('ATCG')
    
    return set(sequence).issubset(valid_nucleotides)


# =============================================================================
# 配置验证函数
# =============================================================================

def validate_sequence_constants() -> Dict[str, bool]:
    """验证所有序列常量的有效性
    
    Returns:
        Dict[str, bool]: 验证结果
    """
    results = {}
    
    try:
        # 验证UTR序列
        utr5 = get_default_utr5()
        utr3 = get_default_utr3()
        
        results['utr5_valid'] = validate_nucleotide_sequence(utr5) and len(utr5) > 0
        results['utr3_valid'] = validate_nucleotide_sequence(utr3) and len(utr3) > 0
        
        # 验证密码子常量
        results['start_codon_valid'] = len(START_CODON) == 3 and validate_nucleotide_sequence(START_CODON)
        results['stop_codons_valid'] = all(len(codon) == 3 and validate_nucleotide_sequence(codon) 
                                         for codon in STOP_CODONS_DNA)
        
        # 验证遗传密码表
        results['genetic_code_complete'] = len(GENETIC_CODE_DNA) == 64
        
        results['overall_valid'] = all(results.values())
        
    except Exception as e:
        results['error'] = str(e)
        results['overall_valid'] = False
    
    return results


# =============================================================================
# 兼容性函数 (为向后兼容保留)
# =============================================================================

def get_default_5utr() -> str:
    """获取默认5'UTR序列 (兼容性函数)"""
    return get_default_utr5()


def get_default_3utr() -> str:
    """获取默认3'UTR序列 (兼容性函数)"""
    return get_default_utr3()


if __name__ == "__main__":
    # 测试代码
    print("序列常量模块测试:")
    print(f"起始密码子: {START_CODON}")
    print(f"终止密码子: {STOP_CODONS_DNA}")
    
    try:
        utr5 = get_default_utr5()
        utr3 = get_default_utr3()
        
        print(f"5'UTR长度: {len(utr5)}")
        print(f"3'UTR长度: {len(utr3)}")
        
        lengths = get_utr_lengths()
        print(f"UTR长度信息: {lengths}")
        
        validation = validate_sequence_constants()
        print(f"序列验证结果: {validation}")
        
    except Exception as e:
        print(f"测试失败: {e}")