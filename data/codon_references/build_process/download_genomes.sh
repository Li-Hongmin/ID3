#!/bin/bash
# 下载E.coli基因组数据脚本

echo "开始下载E.coli基因组数据..."

# 创建下载目录
mkdir -p downloads

# 下载K12 MG1655基因组 (NC_000913.3)
echo "下载E.coli K12 MG1655基因组..."
curl -o downloads/NC_000913.3.gbk \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text"

# 检查K12下载是否成功
if [ -f "downloads/NC_000913.3.gbk" ]; then
    k12_size=$(wc -c < downloads/NC_000913.3.gbk)
    echo "K12基因组下载完成: ${k12_size} bytes"
else
    echo "错误: K12基因组下载失败"
    exit 1
fi

# 下载BL21(DE3)基因组 (NC_012971.2) - GenBank格式
echo "下载E.coli BL21(DE3)基因组 (GenBank格式)..."
curl -o downloads/NC_012971.2.gbk \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012971.2&rettype=gbwithparts&retmode=text"

# 下载BL21(DE3)基因组 (NC_012971.2) - FASTA格式
echo "下载E.coli BL21(DE3)基因组 (FASTA格式)..."
curl -o downloads/NC_012971.2.fasta \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012971.2&rettype=fasta&retmode=text"

# 检查BL21下载是否成功
if [ -f "downloads/NC_012971.2.gbk" ] && [ -f "downloads/NC_012971.2.fasta" ]; then
    bl21_gbk_size=$(wc -c < downloads/NC_012971.2.gbk)
    bl21_fasta_size=$(wc -c < downloads/NC_012971.2.fasta)
    echo "BL21(DE3)基因组下载完成:"
    echo "  GenBank: ${bl21_gbk_size} bytes"
    echo "  FASTA: ${bl21_fasta_size} bytes"
else
    echo "错误: BL21(DE3)基因组下载失败"
    exit 1
fi

echo "所有基因组数据下载完成!"
echo ""
echo "下载的文件:"
ls -lh downloads/
echo ""
echo "下一步: 运行CDS提取脚本"
echo "python scripts/calculate_bl21de3_codon_usage.py --genbank downloads/NC_012971.2.gbk"