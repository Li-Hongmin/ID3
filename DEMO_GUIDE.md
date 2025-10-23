# ID3 Framework Usage Guide

## 两种使用方式

ID3框架提供两种工具：

### 1. `demo.py` - 交互式Demo
**用途**: 快速测试、学习、单个蛋白优化
**特点**:
- ✅ 可视化进度
- ✅ 即时反馈
- ✅ 简单易用

### 2. `run_unified_experiment.py` - 系统化实验
**用途**: 批量实验、论文结果重现、统计分析
**特点**:
- ✅ 批量处理多个蛋白
- ✅ 多随机种子统计
- ✅ 详细结果追踪
- ✅ 12×12实验矩阵

---

## 快速开始

### Demo方式
```bash
# 1. 安装依赖
pip install -r requirements.txt

# 2. 运行demo（首次会自动提示安装DeepRaccess）
python demo.py
```

### 实验方式
```bash
# 快速测试
python run_unified_experiment.py --preset quick-test --device cpu

# 完整实验（论文级别）
python run_unified_experiment.py --preset full-12x12 --device cpu
```

---

## 3种约束机制

所有约束机制都支持与DeepRaccess联合优化：

### 1. Lagrangian (拉格朗日) - 默认，推荐
```bash
python demo.py --constraint lagrangian --iterations 50
```
- 软惩罚优化
- 自适应约束强度
- 论文主要使用的方法

### 2. AMS (Amino Matching Softmax)
```bash
python demo.py --constraint amino_matching --iterations 50
```
- Softmax匹配
- 自然概率约束

### 3. CPC (Codon Profile Constraint)
```bash
python demo.py --constraint codon_profile --iterations 50
```
- 保持密码子使用分布
- 适合微调

---

## 重要参数

### 必需输入
- `--protein` 或 `--protein-file`: 蛋白质序列

### 优化控制
- `--constraint`: 约束机制（lagrangian/amino_matching/codon_profile）
- `--iterations`: 迭代次数（默认20，生产环境建议100-1000）
- `--learning-rate`: 学习率（默认0.01）

### CAI参数
- `--cai-target`: 目标CAI值（默认0.8）
- `--cai-weight`: CAI权重（默认0.1）

### UTR序列
- `--utr5-file`: 自定义5' UTR（默认从data/utr_templates/加载）
- `--utr3-file`: 自定义3' UTR（默认从data/utr_templates/加载）

**为什么UTR重要？**
- RNA可及性预测需要完整的mRNA序列（UTR5 + CDS + UTR3）
- ATG起始密码子附近的可及性对翻译效率至关重要
- 默认UTR来自经典的T7启动子和终止子

---

## 使用示例

### 示例1: 基本优化
```bash
python demo.py --protein MSKGEELFTGVVPILVELDGDVNGHKFSVSGEG --iterations 50
```

输出:
```
✓ RNA accessibility optimized (DeepRaccess)
✓ CAI optimized (target: 0.8, achieved: 0.895)
✓ Amino acid constraints maintained
✓ Final accessibility score: 1.344
```

### 示例2: 从FASTA文件
```bash
python demo.py --protein-file data/proteins/P04637.fasta --iterations 100
```

### 示例3: 对比3种约束机制
```bash
for constraint in lagrangian amino_matching codon_profile; do
    python demo.py --constraint $constraint --iterations 50 --output result_$constraint.fasta
done
```

### 示例4: 自定义UTR
```bash
# 使用自己的UTR序列
python demo.py --protein MSKGEELFT \
               --utr5-file my_custom_5utr.txt \
               --utr3-file my_custom_3utr.txt \
               --iterations 100
```

---

## 优化目标

Demo优化以下目标的联合损失：

```
L_total = L_accessibility + λ_constraint·C_amino + λ_CAI·L_CAI

其中:
- L_accessibility: DeepRaccess预测的RNA可及性（越低越好）
- C_amino: 氨基酸约束惩罚（确保正确翻译）
- L_CAI: CAI损失（接近目标CAI值）
```

---

## 输出说明

Demo会显示：
1. **优化过程**: 进度条显示 loss/acc/cai
2. **最佳结果**: 最低损失时的可及性和CAI
3. **最终评估**: 离散序列的实际指标
4. **约束验证**: 确认氨基酸序列正确

保存的文件包含：
```
>MSKGEELFT...
# Accessibility: 1.3444
# CAI: 0.8953
AUGAGCAAAGGCGAAGAACUGUUCACU...
```

---

## 常见问题

### Q: 需要多少次迭代？
- 快速测试: 20-50次（几秒）
- 正常优化: 100-200次（几分钟）
- 论文结果: 1000次（30-45分钟）

### Q: 哪个约束机制最好？
- 不确定 → 用Lagrangian（默认）
- 需要稳定性 → Lagrangian
- 对比实验 → 三种都试试

### Q: 可及性分数是什么意思？
- DeepRaccess输出的自由能估计
- 越低 = RNA越accessible（核糖体更容易结合）
- 典型范围: -5 到 +10

### Q: 为什么需要UTR？
- DeepRaccess需要完整mRNA序列来预测可及性
- ATG起始密码子附近的二级结构对翻译很重要
- 论文使用-19到+15位置的35nt窗口

---

## 技术细节

### 梯度流
```
约束参数 → 软概率 → DeepRaccess → 可及性 → 总损失 → 反向传播
```

### 完整mRNA结构
```
5' UTR (70nt) + CDS (protein_len×3) + 3' UTR (63nt)
         ↓
    ATG position
    (起始密码子)
         ↓
   优化这个区域的可及性
```

### ATG窗口优化
- 提取ATG位置-19到+15的35nt窗口
- 基于Terai & Asai (2020)的研究
- 这个区域的可及性对核糖体结合最关键
