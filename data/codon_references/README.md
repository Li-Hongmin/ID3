# E.coli Codon Reference Systems

用于CAI（Codon Adaptation Index）分析的E.coli密码子参考数据。

## 参考系统

### E.coli K12 参考系统
- **文件**: `ecoli_k12_reference_sequences.json`
- **来源**: NCBI RefSeq NC_000913.3 (K12 MG1655)
- **基因数量**: 82个高表达基因 (28原始 + 54恢复的核糖体蛋白)
- **用途**: 通用E.coli研究应用

### E.coli BL21(DE3) 参考系统 ⭐ **推荐用于重组蛋白生产**
- **文件**: `ecoli_bl21de3_reference_sequences.json`  
- **来源**: NCBI RefSeq NC_012971.2 (BL21(DE3))
- **基因数量**: 83个表达加权基因 (29原始 + 54恢复的核糖体蛋白)
- **用途**: pET表达系统优化

## 使用方法

### 专业CAI库 (Benjamin Lee) - 当前推荐
```python
from CAI import CAI, RSCU
import json

# 加载K12参考序列
with open('data/codon_references/ecoli_k12_reference_sequences.json') as f:
    k12_data = json.load(f)
    k12_sequences = list(k12_data['sequences'].values())

# 加载BL21(DE3)参考序列  
with open('data/codon_references/ecoli_bl21de3_reference_sequences.json') as f:
    bl21_data = json.load(f)
    bl21_sequences = list(bl21_data['sequences'].values())

# 计算CAI
dna_sequence = "ATGAAAAAACTGCTGGTG..."  # 你的DNA序列
k12_cai = CAI(dna_sequence, reference=k12_sequences)
bl21_cai = CAI(dna_sequence, reference=bl21_sequences)

print(f"K12 CAI: {k12_cai:.4f}")
print(f"BL21(DE3) CAI: {bl21_cai:.4f}")
```

### CAI分析工具
```bash
# 使用真实蛋白质进行专业CAI分析
python tools/cai_real_protein_analysis.py

# CAI理论分析和方法比较
python tools/cai_theoretical_analysis.py
python tools/cai_professional_analysis.py
```

## CAI库安装
```bash
# 安装Benjamin Lee的专业CAI库
pip install CAI

# 验证安装
python -c "from CAI import CAI; print('CAI library installed successfully')"
```

## 推荐使用
- **一般E.coli研究**: 使用K12参考系统
- **重组蛋白生产**: 使用BL21(DE3)参考系统（推荐）
- **CAI目标值**: 0.6-0.9 (实际优化范围)

## JSON文件构建过程

**当前文件来源**：这两个JSON文件包含从NCBI基因组中提取的真实CDS序列：
- `ecoli_k12_reference_sequences.json` - 82个K12高表达基因的CDS序列
- `ecoli_bl21de3_reference_sequences.json` - 83个BL21(DE3)高表达基因的CDS序列

**完整构建流程**：详见 `build_process/` 目录，包含完整的从下载到生成的全过程：

### 构建流程概览
```
下载基因组 → 定义高表达基因 → 提取CDS序列 → 生成JSON文件
     ↓              ↓              ↓            ↓
build_process/  build_process/  build_process/    当前目录/
downloads/      gene_lists/     scripts/       *.json文件
```

### 1. 基因组数据下载
```bash
# 执行下载脚本
cd build_process/
bash download_genomes.sh
```
下载文件：
- `NC_000913.3.gbk` (11.9MB) - K12基因组
- `NC_012971.2.gbk` (9.9MB) - BL21(DE3)基因组
- `NC_012971.2.fasta` (4.6MB) - BL21(DE3) FASTA

### 2. 高表达基因定义
基因列表位于 `build_process/gene_lists/`:
- **K12系统**: 82个基因（核糖体蛋白、翻译因子、糖酵解酶）
- **BL21(DE3)系统**: 83个基因，按表达量加权
  - 核糖体蛋白（权重1.0）
  - 翻译因子（权重0.8）
  - 糖酵解酶（权重0.6）

### 3. CDS序列提取
```bash
cd build_process/scripts/
python calculate_bl21de3_codon_usage.py --genbank ../downloads/NC_012971.2.gbk
```

### 4. 重新构建完整过程
详细步骤和所有文件都在 `build_process/README.md` 中记录。

**重建命令**：
```bash
cd build_process/
bash download_genomes.sh                    # 下载基因组
python scripts/calculate_bl21de3_codon_usage.py --genbank downloads/NC_012971.2.gbk  # 提取序列
```

## 注意事项
1. **输入格式**: 使用DNA序列（ATCG），不是RNA序列（AUCG）
2. **序列要求**: 确保序列长度是3的倍数，以ATG开始
3. **参考序列**: 使用本目录中的真实基因组提取序列
4. **数据来源**: 所有序列均来自NCBI RefSeq，经过生物信息学验证

---
**最后更新**: 2025-07-24  
**验证状态**: ✅ 所有参考系统已通过专业CAI分析验证