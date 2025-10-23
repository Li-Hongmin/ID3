# E.coli CAI参考数据构建过程

此目录完整记录了E.coli K12和BL21(DE3) CAI参考序列的构建过程。

## 构建流程概览

```
下载基因组 → 定义高表达基因 → 提取CDS序列 → 生成JSON文件
     ↓              ↓              ↓            ↓
  downloads/    gene_lists/    scripts/    ../ecoli_*.json
```

## 步骤1: 基因组数据下载

### 下载命令
```bash
# K12 MG1655基因组
curl -o downloads/NC_000913.3.gbk \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text"

# BL21(DE3)基因组 (GenBank格式)
curl -o downloads/NC_012971.2.gbk \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012971.2&rettype=gbwithparts&retmode=text"

# BL21(DE3)基因组 (FASTA格式)
curl -o downloads/NC_012971.2.fasta \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012971.2&rettype=fasta&retmode=text"
```

### 下载的文件
- `downloads/NC_000913.3.gbk` (11.9MB) - E.coli K12 MG1655完整基因组
- `downloads/NC_012971.2.gbk` (9.9MB) - E.coli BL21(DE3)完整基因组  
- `downloads/NC_012971.2.fasta` (4.6MB) - BL21(DE3) FASTA序列

## 步骤2: 高表达基因定义

### 基因分类标准
基于生物学功能和表达水平，选择以下类别的基因：

#### K12系统 (28个基因)
- **核糖体蛋白**: rpsA-rpsU, rplA-rplY 
- **翻译因子**: fusA, tufA, tufB, tsf, infA, infB, infC, prfA, prfB, prfC
- **糖酵解酶**: pgi, pfkA, pfkB, fbaA, tpiA, gapA, pgk, gpmA, eno, pykA, pykF

#### BL21(DE3)系统 (29个基因，按表达量加权)
- **很高表达 (权重1.0)**: 核糖体蛋白 
- **高表达 (权重0.8)**: 翻译因子
- **中高表达 (权重0.6)**: 糖酵解酶

### 基因列表文件
- `gene_lists/highly_expressed_genes.json` - BL21(DE3)完整基因列表和权重
- `gene_lists/ribosomal_proteins.json` - 核糖体蛋白基因
- `gene_lists/translation_factors.json` - 翻译因子基因  
- `gene_lists/glycolysis_enzymes.json` - 糖酵解酶基因

## 步骤3: CDS序列提取

### 处理脚本
- `scripts/calculate_bl21de3_codon_usage.py` - 主要处理脚本

### 脚本功能
1. **GenBank解析**: 从.gbk文件中提取基因注释信息
2. **序列匹配**: 根据基因名匹配高表达基因列表
3. **CDS提取**: 从基因组序列中提取完整编码序列
4. **格式转换**: 生成标准化JSON格式输出

### 执行命令
```bash
cd scripts/
python calculate_bl21de3_codon_usage.py \\
    --genbank ../downloads/NC_012971.2.gbk \\
    --fasta ../downloads/NC_012971.2.fasta \\
    --reference-genes ../gene_lists \\
    --output ../intermediate_data
```

## 步骤4: 最终JSON文件生成

### 输出文件结构
```json
{
  "metadata": {
    "description": "E.coli [K12|BL21(DE3)] highly expressed gene CDS sequences",
    "source_genome": "NC_000913.3|NC_012971.2", 
    "total_sequences": 28|29,
    "extraction_date": "2025-07-24"
  },
  "sequences": {
    "gene_name": "ATGAAAAAACTGCTGGTG...",
    ...
  }
}
```

### 最终文件
- `../ecoli_k12_reference_sequences.json` - K12参考序列
- `../ecoli_bl21de3_reference_sequences.json` - BL21(DE3)参考序列

## 数据验证

### 质量检查
1. **序列完整性**: 所有CDS以ATG开始，长度为3的倍数
2. **基因覆盖率**: K12 (28/28), BL21(DE3) (29/29)  
3. **生物学验证**: 与Benjamin Lee CAI库兼容测试通过

### 文献参考
- Sharp & Li (1987) - CAI计算方法
- Nucleic Acids Research (2017) - BL21(DE3)基因组分析
- NCBI RefSeq - 官方基因组注释

## 重现完整构建过程

```bash
# 1. 下载基因组数据
bash download_genomes.sh

# 2. 提取CDS序列  
python scripts/calculate_bl21de3_codon_usage.py --genbank downloads/NC_012971.2.gbk

# 3. 复制到最终位置
cp intermediate_data/reference_sequences.json ../ecoli_bl21de3_reference_sequences.json
```

## 目录结构
```
build_process/
├── README.md                          # 本文件
├── downloads/                         # 原始基因组数据
│   ├── NC_000913.3.gbk               # K12基因组 (11.9MB)
│   ├── NC_012971.2.gbk               # BL21基因组 (9.9MB)  
│   └── NC_012971.2.fasta             # BL21 FASTA (4.6MB)
├── gene_lists/                        # 高表达基因定义
│   ├── highly_expressed_genes.json   # 完整基因列表
│   ├── ribosomal_proteins.json       # 核糖体蛋白
│   ├── translation_factors.json      # 翻译因子
│   └── glycolysis_enzymes.json       # 糖酵解酶
├── scripts/                           # 处理脚本
│   └── calculate_bl21de3_codon_usage.py
└── intermediate_data/                 # 中间处理结果
    └── (生成的临时文件)
```

---
**创建日期**: 2025-07-24  
**数据来源**: NCBI RefSeq, KEGG, UniProt  
**验证状态**: ✅ 已通过Benjamin Lee CAI库验证