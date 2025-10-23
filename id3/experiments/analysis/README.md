# ID3实验分析脚本目录

## 📋 脚本分类总览

### 🎯 论文图表生成脚本 (核心)
这些脚本用于生成论文中的主要表格和图表：

#### 表格生成
- **`generate_tab1_access_comparison.py`** - 生成Tab1: Access-only实验结果对比表
- **`generate_tab2_access_cai_comparison.py`** - 生成Tab2: Access+CAI实验结果对比表

#### 收敛分析图表
- **`generate_fig_convergence_access_only.py`** - 生成Access-only收敛分析图表
  - 输出: `fig_convergence_analysis_access_only.png/.pdf`
- **`generate_fig_convergence_cai.py`** - 生成Access+CAI收敛分析图表  
  - 输出: `fig_convergence_analysis_cai.png/.pdf`

#### 性能分析图表
- **`generate_convergence_performance_fast.py`** - 生成性能对比图表
  - 输出: `fig_convergence_performance.png/.pdf`

### 📊 数据处理脚本
这些脚本用于处理和提取实验数据：

#### 轨迹数据提取
- **`extract_access_only_trajectories.py`** - 并行提取Access-only实验轨迹数据
  - 输出: `trajectory_access_only_1000steps.csv`
- **`extract_cai_trajectories.py`** - 并行提取CAI实验轨迹数据
  - 输出: `trajectory_cai_1000steps.csv`

#### 数据探索
- **`explore_trajectory_data.py`** - 探索轨迹数据结构和格式

### 🔍 验证和质控脚本
这些脚本用于验证实验结果的正确性：

- **`verify_l11_cai.py`** - 验证L11配置的CAI数值和约束满足
- **`experiment_validator.py`** - 实验结果验证工具
- **`experiment_comparator.py`** - 实验结果对比工具

### 🛠️ 核心分析工具
底层分析工具和解析器：

- **`fast_analyzer.py`** - 快速实验结果分析器
- **`result_parser.py`** - 实验结果解析器
- **`unified_experiment_analyzer.py`** - 统一实验分析器
- **`performance_tables.py`** - 性能表格生成工具
- **`visualization.py`** - 可视化工具集

### 📁 已废弃脚本 (`deprecated/`)
这些脚本已被更好的版本替代或不再使用：

- 各种调试脚本和过时的分析工具
- 早期版本的收敛分析脚本
- 实验性的插值分析脚本

## 🚀 快速使用指南

### 生成论文图表
```bash
# 生成表格
python generate_tab1_access_comparison.py
python generate_tab2_access_cai_comparison.py

# 生成收敛分析图
python generate_fig_convergence_access_only.py
python generate_fig_convergence_cai.py

# 生成性能对比图
python generate_convergence_performance_fast.py
```

### 提取轨迹数据
```bash
# 提取Access-only轨迹
python extract_access_only_trajectories.py

# 提取CAI轨迹
python extract_cai_trajectories.py
```

### 验证实验结果
```bash
# 验证L11实验
python verify_l11_cai.py

# 一般验证
python experiment_validator.py
```

## 📊 输出文件说明

### 表格文件 (保存到 `/paper_experiment_results/tables/`)
- `data_tab1_access_comparison.csv` - Tab1原始数据
- `data_tab2_access_cai_comparison.json` - Tab2原始数据
- `tab2_access_cai_comparison.tex` - Tab2 LaTeX格式

### 图表文件 (保存到 `/paper_experiment_results/figures/`)
- `fig_convergence_analysis_access_only.png/.pdf` - Access-only收敛图
- `fig_convergence_analysis_cai.png/.pdf` - CAI收敛图
- `fig_convergence_performance.png/.pdf` - 性能对比图

### 轨迹数据文件 (保存到 `/paper_experiment_results/figures/`)
- `trajectory_access_only_1000steps.csv` - Access-only轨迹数据
- `trajectory_cai_1000steps.csv` - CAI轨迹数据

## 🔗 依赖关系

### 数据源依赖
- Access-only实验: `/paper_experiment_results/20250910_174703_unified_access_experiments/`
- CAI实验: `/paper_experiment_results/merged_cai_experiments/`

### 脚本间依赖
- 图表生成脚本依赖对应的轨迹提取脚本
- 某些分析脚本依赖 `result_parser.py` 和 `unified_experiment_analyzer.py`

## 🎯 论文使用说明

1. **Tab1和Tab2**: 直接运行对应的生成脚本即可
2. **收敛分析图**: 需要先运行轨迹提取脚本，再运行图表生成脚本
3. **所有脚本都命名直观**: 文件名直接反映其功能和输出

## ⚠️ 注意事项

- 所有脚本都假设从项目根目录运行
- 轨迹提取需要较大内存(处理144,000行数据)
- 图表生成需要matplotlib和seaborn
- 建议按照依赖顺序运行脚本