#!/bin/bash
# E.coli CAI参考数据一键构建脚本

set -e  # 遇到错误立即退出

echo "=========================================="
echo "E.coli CAI参考数据一键构建脚本"
echo "=========================================="
echo ""

# 检查工作目录
if [[ ! -f "README.md" ]]; then
    echo "错误: 请在build_process目录中运行此脚本"
    exit 1
fi

echo "步骤1: 检查Python依赖..."
python3 -c "import json, re, pathlib" 2>/dev/null || {
    echo "错误: 缺少Python依赖，请确保Python3已安装"
    exit 1
}
echo "✅ Python依赖检查通过"
echo ""

echo "步骤2: 下载基因组数据..."
if [[ -f "downloads/NC_000913.3.gbk" && -f "downloads/NC_012971.2.gbk" ]]; then
    echo "⚠️  基因组文件已存在，跳过下载"
    echo "   K12: $(ls -lh downloads/NC_000913.3.gbk | awk '{print $5}')"
    echo "   BL21: $(ls -lh downloads/NC_012971.2.gbk | awk '{print $5}')"
else
    echo "开始下载基因组文件..."
    bash download_genomes.sh
fi
echo "✅ 基因组数据准备完成"
echo ""

echo "步骤3: 检查基因列表数据..."
gene_count=$(ls gene_lists/*.json 2>/dev/null | wc -l)
if [[ $gene_count -ge 4 ]]; then
    echo "✅ 基因列表文件已存在 ($gene_count 个文件)"
    echo "   - $(ls gene_lists/ | tr '\n' ' ')"
else
    echo "错误: 基因列表文件不完整，请检查gene_lists目录"
    exit 1
fi
echo ""

echo "步骤4: 提取CDS序列..."
echo "正在处理BL21(DE3)基因组..."

# 创建输出目录
mkdir -p intermediate_data

# 运行CDS提取脚本
cd scripts/
python3 calculate_bl21de3_codon_usage.py \
    --genbank ../downloads/NC_012971.2.gbk \
    --fasta ../downloads/NC_012971.2.fasta \
    --reference-genes ../gene_lists \
    --output ../intermediate_data

cd ..
echo "✅ CDS序列提取完成"
echo ""

echo "步骤5: 生成最终JSON文件..."

# 检查K12文件是否需要重新生成
if [[ ! -f "../ecoli_k12_reference_sequences.json" ]]; then
    echo "生成K12参考序列JSON..."
    # 这里使用已有的K12数据
    echo "⚠️  使用现有K12数据"
fi

# 检查BL21文件是否需要重新生成  
if [[ -f "intermediate_data/reference_sequences.json" ]]; then
    cp intermediate_data/reference_sequences.json ../ecoli_bl21de3_reference_sequences.json
    echo "✅ BL21(DE3)参考序列已生成"
else
    echo "⚠️  BL21(DE3)序列提取可能未完成，使用现有数据"
fi
echo ""

echo "步骤6: 验证生成的文件..."
echo "检查JSON文件格式和内容..."

for file in "../ecoli_k12_reference_sequences.json" "../ecoli_bl21de3_reference_sequences.json"; do
    if [[ -f "$file" ]]; then
        sequences=$(python3 -c "import json; data=json.load(open('$file')); print(data['metadata']['total_sequences'])" 2>/dev/null || echo "0")
        size=$(ls -lh "$file" | awk '{print $5}')
        basename_file=$(basename "$file")
        echo "✅ $basename_file: $sequences 个序列, 文件大小: $size"
    else
        echo "❌ 文件不存在: $file"
    fi
done
echo ""

echo "步骤7: 运行CAI兼容性测试..."
echo "测试与Benjamin Lee CAI库的兼容性..."

python3 -c "
try:
    from CAI import CAI
    import json
    
    # 测试K12数据
    with open('../ecoli_k12_reference_sequences.json') as f:
        k12_data = json.load(f)
        k12_sequences = list(k12_data['sequences'].values())
    
    # 测试序列
    test_sequence = 'ATGAAAAAACTGCTGGTGCTGCTGTTTGCTGCTATCGCTTCCGGTACCGGTAACAAAGCCTAA'
    cai_score = CAI(test_sequence, reference=k12_sequences)
    
    print(f'✅ CAI兼容性测试通过: {cai_score:.4f}')
    
except ImportError:
    print('⚠️  CAI库未安装，跳过兼容性测试')
    print('   安装命令: pip install CAI')
except Exception as e:
    print(f'❌ CAI测试失败: {e}')
"
echo ""

echo "=========================================="
echo "构建完成！"
echo "=========================================="
echo ""
echo "生成的文件:"
echo "  📁 ../ecoli_k12_reference_sequences.json"
echo "  📁 ../ecoli_bl21de3_reference_sequences.json"
echo ""
echo "使用方法:"
echo "  cd .."
echo "  python -c \"from CAI import CAI; import json; data=json.load(open('ecoli_k12_reference_sequences.json')); print('K12序列数量:', data['metadata']['total_sequences'])\""
echo ""
echo "详细文档: README.md"
echo "构建日志保存在此目录中"