#!/usr/bin/env python3
"""
计算E.coli BL21(DE3)密码子使用频率

从真实的BL21(DE3)基因组序列中提取高表达基因，
计算加权的密码子使用频率，用于CAI优化。
"""

import json
import re
from typing import Dict, List, Tuple
from collections import defaultdict, Counter
from pathlib import Path
import argparse

# 遗传密码表 (DNA -> Amino Acid)
GENETIC_CODE = {
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


class BL21DE3CodonUsageCalculator:
    """BL21(DE3)密码子使用频率计算器"""
    
    def __init__(self, reference_genes_dir: str):
        """
        初始化计算器
        
        Args:
            reference_genes_dir: 参考基因目录路径
        """
        self.reference_genes_dir = Path(reference_genes_dir)
        self.highly_expressed_genes = self._load_highly_expressed_genes()
        self.expression_weights = self._load_expression_weights()
        
    def _load_highly_expressed_genes(self) -> Dict:
        """加载高表达基因列表"""
        genes_file = self.reference_genes_dir / "highly_expressed_genes.json"
        with open(genes_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    def _load_expression_weights(self) -> Dict[str, float]:
        """加载基因表达权重"""
        weights = {}
        gene_data = self.highly_expressed_genes
        
        for weight_category in gene_data["expression_weight_assignments"]["assignments"]:
            weight_info = gene_data["expression_weight_assignments"]["assignments"][weight_category]
            weight_value = weight_info["weight"]
            
            for gene in weight_info["genes"]:
                weights[gene] = weight_value
                
        return weights
    
    def parse_genbank_file(self, genbank_file: str, fasta_file: str) -> Dict[str, str]:
        """
        解析GenBank文件，提取CDS序列
        
        Args:
            genbank_file: GenBank文件路径
            fasta_file: 对应的FASTA基因组文件
            
        Returns:
            gene_sequences: {gene_name: cds_sequence}
        """
        gene_sequences = {}
        
        try:
            # 读取FASTA文件获取完整基因组序列
            with open(fasta_file, 'r', encoding='utf-8') as f:
                fasta_content = f.read()
            
            # 提取基因组序列（去除标题行）
            genome_sequence = ''.join(line.strip() for line in fasta_content.split('\n') if not line.startswith('>'))
            
            # 读取GenBank文件
            with open(genbank_file, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # 改进的CDS提取模式
            cds_pattern = r'     CDS\s+(?:complement\()?(\d+)\.\.(\d+)\)?.*?(?:/gene="([^"]+)".*?)?(?=/\w+|     \w|\nORIGIN|\n\n|\nLOCUS|\Z)'
            cds_matches = re.findall(cds_pattern, content, re.MULTILINE | re.DOTALL)
            
            for match in cds_matches:
                start, end, gene_name = match
                start_pos = int(start) - 1  # 转为0-based索引
                end_pos = int(end)
                
                if gene_name and gene_name in self.expression_weights:
                    # 提取CDS序列
                    cds_sequence = genome_sequence[start_pos:end_pos]
                    gene_sequences[gene_name] = cds_sequence
                        
        except Exception as e:
            print(f"解析GenBank文件时出错: {e}")
            print("注意: 此处需要真实的GenBank文件来提取CDS序列")
            
        return gene_sequences
    
    def extract_codons_from_sequence(self, sequence: str) -> List[str]:
        """
        从DNA序列中提取密码子
        
        Args:
            sequence: DNA序列
            
        Returns:
            codons: 密码子列表
        """
        # 确保序列长度是3的倍数
        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]
        
        codons = []
        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].upper()
            if len(codon) == 3 and all(base in 'ATCG' for base in codon):
                codons.append(codon)
                
        return codons
    
    def calculate_weighted_codon_frequencies(self, gene_sequences: Dict[str, str]) -> Dict[str, float]:
        """
        计算加权的密码子使用频率
        
        Args:
            gene_sequences: {gene_name: sequence}
            
        Returns:
            codon_frequencies: {codon: frequency_per_1000}
        """
        weighted_codon_counts = defaultdict(float)
        total_weighted_codons = 0
        
        for gene_name, sequence in gene_sequences.items():
            if gene_name not in self.expression_weights:
                continue
                
            weight = self.expression_weights[gene_name]
            codons = self.extract_codons_from_sequence(sequence)
            
            # 应用权重到密码子计数
            for codon in codons:
                if codon in GENETIC_CODE and GENETIC_CODE[codon] != '*':  # 排除终止密码子
                    weighted_codon_counts[codon] += weight
                    total_weighted_codons += weight
        
        # 转换为每1000个密码子的频率
        codon_frequencies = {}
        if total_weighted_codons > 0:
            for codon, count in weighted_codon_counts.items():
                codon_frequencies[codon] = (count / total_weighted_codons) * 1000
        
        return codon_frequencies
    
    def calculate_relative_adaptiveness(self, codon_frequencies: Dict[str, float]) -> Dict[str, float]:
        """
        计算相对适应性权重 (w_i = f_i / max(f_j))
        
        Args:
            codon_frequencies: 密码子频率
            
        Returns:
            relative_weights: 相对适应性权重
        """
        # 按氨基酸分组
        aa_groups = defaultdict(list)
        for codon, freq in codon_frequencies.items():
            if codon in GENETIC_CODE and GENETIC_CODE[codon] != '*':
                aa = GENETIC_CODE[codon]
                aa_groups[aa].append((codon, freq))
        
        # 计算每个氨基酸内的相对权重
        relative_weights = {}
        for aa, codon_freqs in aa_groups.items():
            if codon_freqs:
                max_freq = max(freq for _, freq in codon_freqs)
                for codon, freq in codon_freqs:
                    relative_weights[codon] = freq / max_freq if max_freq > 0 else 0
        
        return relative_weights
    
    def generate_output_files(self, codon_frequencies: Dict[str, float], 
                            relative_weights: Dict[str, float], 
                            output_dir: str):
        """
        生成输出文件
        
        Args:
            codon_frequencies: 密码子频率
            relative_weights: 相对适应性权重  
            output_dir: 输出目录
        """
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        # 1. BL21(DE3)密码子使用频率文件
        frequency_data = {
            "metadata": {
                "source": "E.coli BL21(DE3) highly expressed genes",
                "calculation_method": "Expression-weighted codon counting",
                "reference_genes": len(self.expression_weights),
                "date_generated": "2025-01-15",
                "genome_source": "NCBI RefSeq NC_012971.2"
            },
            "raw_frequencies": codon_frequencies,
            "relative_adaptiveness": relative_weights,
            "rare_codons": self._identify_rare_codons(codon_frequencies),
            "preferred_codons": self._identify_preferred_codons(codon_frequencies)
        }
        
        freq_file = output_path / "bl21de3_codon_frequencies.json"
        with open(freq_file, 'w', encoding='utf-8') as f:
            json.dump(frequency_data, f, indent=2, ensure_ascii=False)
        
        print(f"生成频率文件: {freq_file}")
        
        # 2. 与K12对比文件
        self._generate_k12_comparison(codon_frequencies, output_path)
        
        # 3. 稀有密码子分析文件
        self._generate_rare_codon_analysis(codon_frequencies, output_path)
    
    def _identify_rare_codons(self, frequencies: Dict[str, float], threshold: float = 5.0) -> Dict[str, float]:
        """识别稀有密码子 (< 5/1000)"""
        return {codon: freq for codon, freq in frequencies.items() if freq < threshold}
    
    def _identify_preferred_codons(self, frequencies: Dict[str, float], threshold: float = 25.0) -> Dict[str, float]:
        """识别偏好密码子 (> 25/1000)"""
        return {codon: freq for codon, freq in frequencies.items() if freq > threshold}
    
    def _generate_k12_comparison(self, bl21_frequencies: Dict[str, float], output_path: Path):
        """生成与K12的对比数据"""
        # 加载K12数据 (来自已有的Kazusa数据)
        k12_file = Path("data/codon_references/ecoli_k12_kazusa.json")
        
        if k12_file.exists():
            with open(k12_file, 'r', encoding='utf-8') as f:
                k12_data = json.load(f)
            
            k12_frequencies = k12_data["converted_data_dna_codons"]
            
            # 计算差异
            comparison = {
                "metadata": {
                    "bl21_source": "Calculated from BL21(DE3) highly expressed genes",
                    "k12_source": "Kazusa Codon Usage Database",
                    "comparison_date": "2025-01-15"
                },
                "frequency_differences": {},
                "ratio_bl21_to_k12": {}
            }
            
            for codon in set(bl21_frequencies.keys()) | set(k12_frequencies.keys()):
                bl21_freq = bl21_frequencies.get(codon, 0)
                k12_freq = k12_frequencies.get(codon, 0)
                
                comparison["frequency_differences"][codon] = bl21_freq - k12_freq
                if k12_freq > 0:
                    comparison["ratio_bl21_to_k12"][codon] = bl21_freq / k12_freq
            
            comp_file = output_path / "comparison_with_k12.json"
            with open(comp_file, 'w', encoding='utf-8') as f:
                json.dump(comparison, f, indent=2, ensure_ascii=False)
            
            print(f"生成K12对比文件: {comp_file}")
    
    def _generate_rare_codon_analysis(self, frequencies: Dict[str, float], output_path: Path):
        """生成稀有密码子分析"""
        rare_codons = self._identify_rare_codons(frequencies)
        
        analysis = {
            "metadata": {
                "threshold": "< 5 per 1000 codons",
                "total_rare_codons": len(rare_codons),
                "analysis_date": "2025-01-15"
            },
            "rare_codons": rare_codons,
            "amino_acid_analysis": {}
        }
        
        # 按氨基酸分析稀有密码子
        for codon, freq in rare_codons.items():
            if codon in GENETIC_CODE:
                aa = GENETIC_CODE[codon]
                if aa not in analysis["amino_acid_analysis"]:
                    analysis["amino_acid_analysis"][aa] = []
                analysis["amino_acid_analysis"][aa].append({"codon": codon, "frequency": freq})
        
        rare_file = output_path / "rare_codons_analysis.json"
        with open(rare_file, 'w', encoding='utf-8') as f:
            json.dump(analysis, f, indent=2, ensure_ascii=False)
        
        print(f"生成稀有密码子分析文件: {rare_file}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="计算BL21(DE3)密码子使用频率")
    parser.add_argument("--genbank", required=True, help="BL21(DE3) GenBank文件路径")
    parser.add_argument("--fasta", help="对应的FASTA基因组文件路径（如未指定，将自动推导）")
    parser.add_argument("--reference-genes", default="data/codon_references/ecoli_bl21de3_reference/reference_genes", 
                       help="参考基因目录")
    parser.add_argument("--output", default="data/codon_references/ecoli_bl21de3_reference/codon_usage",
                       help="输出目录")
    
    args = parser.parse_args()
    
    # 创建计算器实例
    calculator = BL21DE3CodonUsageCalculator(args.reference_genes)
    
    print("开始计算BL21(DE3)密码子使用频率...")
    print(f"参考基因数量: {len(calculator.expression_weights)}")
    
    # 解析GenBank文件和FASTA文件
    print("解析GenBank文件...")
    fasta_file = args.fasta if args.fasta else args.genbank.replace('.gbk', '.fasta')
    gene_sequences = calculator.parse_genbank_file(args.genbank, fasta_file)
    print(f"成功提取基因序列: {len(gene_sequences)}")
    
    if not gene_sequences:
        print("警告: 未能提取到基因序列，请检查GenBank文件格式")
        print("注意: 需要先下载真实的BL21(DE3) GenBank文件")
        return
    
    # 计算密码子频率
    print("计算加权密码子频率...")
    codon_frequencies = calculator.calculate_weighted_codon_frequencies(gene_sequences)
    print(f"计算得到 {len(codon_frequencies)} 个密码子的频率")
    
    # 计算相对适应性权重
    print("计算相对适应性权重...")
    relative_weights = calculator.calculate_relative_adaptiveness(codon_frequencies)
    
    # 生成输出文件
    print("生成输出文件...")
    calculator.generate_output_files(codon_frequencies, relative_weights, args.output)
    
    print("密码子使用频率计算完成!")
    print(f"输出文件保存在: {args.output}")


if __name__ == "__main__":
    main()