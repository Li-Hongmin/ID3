# ID3分析脚本使用优先级指南

## 🚀 论文发表级脚本 (最高优先级)
这些脚本直接用于论文图表生成，必须保持完美运行状态：

### 表格生成 (Tab1-2)
```bash
# Priority 1: 论文表格
python generate_tab1_access_comparison.py     # → Tab1: Access-only对比
python generate_tab2_access_cai_comparison.py # → Tab2: Access+CAI对比
```

### 收敛分析图表 (Fig)
```bash
# Priority 1: 论文主图
python extract_access_only_trajectories.py    # → 提取Access-only轨迹
python generate_fig_convergence_access_only.py # → Access-only收敛图

python extract_cai_trajectories.py            # → 提取CAI轨迹  
python generate_fig_convergence_cai.py        # → CAI收敛图

python generate_convergence_performance_fast.py # → 性能对比图
```

## 🔍 质量控制脚本 (高优先级)
用于验证实验数据的正确性：

```bash
# Priority 2: 数据验证
python verify_l11_cai.py                      # → L11配置验证
python experiment_validator.py                # → 通用实验验证
python explore_trajectory_data.py             # → 数据结构探索
```

## 🛠️ 开发支持脚本 (中优先级)
用于数据处理和分析的底层工具：

```bash
# Priority 3: 分析工具
python fast_analyzer.py                       # → 快速分析
python unified_experiment_analyzer.py         # → 统一分析器
python result_parser.py                       # → 结果解析
python performance_tables.py                  # → 性能表格工具
```

## 📊 可选分析脚本 (低优先级)
用于深入分析和比较的工具：

```bash
# Priority 4: 扩展分析
python experiment_comparator.py               # → 实验对比
python visualization.py                       # → 可视化工具
python core/sequence_uniqueness_analyzer.py   # → 序列唯一性分析
python core/diversity_metrics.py              # → 多样性指标
```

## 🗂️ 支持模块 (库文件)
这些是被其他脚本调用的支持模块：

```bash
# Priority 5: 支持库
python core/base_analyzer.py                  # → 基础分析类
python utils/sequence_utils.py                # → 序列工具函数
python configs/analysis_config.py             # → 配置文件
python visualizers/uniqueness_visualizer.py   # → 唯一性可视化
```

## ⚡ 快速命令执行顺序

### 论文图表完整生成流程：
```bash
# 1. 生成表格（独立运行）
python generate_tab1_access_comparison.py
python generate_tab2_access_cai_comparison.py

# 2. 生成收敛分析图（有依赖关系）
python extract_access_only_trajectories.py
python generate_fig_convergence_access_only.py

python extract_cai_trajectories.py  
python generate_fig_convergence_cai.py

# 3. 生成性能对比图（独立运行）
python generate_convergence_performance_fast.py

# 4. 验证结果
python verify_l11_cai.py
```

### 一键执行脚本（建议创建）：
```bash
#!/bin/bash
# generate_all_figures.sh

echo "🚀 开始生成所有论文图表..."

echo "📊 生成表格..."
python generate_tab1_access_comparison.py
python generate_tab2_access_cai_comparison.py

echo "📈 提取轨迹数据..."
python extract_access_only_trajectories.py &
python extract_cai_trajectories.py &
wait

echo "🎨 生成收敛分析图..."
python generate_fig_convergence_access_only.py &
python generate_fig_convergence_cai.py &
python generate_convergence_performance_fast.py &
wait

echo "✅ 验证结果..."
python verify_l11_cai.py

echo "🎯 所有图表生成完成！"
```

## 📋 脚本状态总结

- **总脚本数**: 27个活跃脚本
- **已废弃**: 13个脚本移至 `deprecated/` 目录
- **论文核心**: 8个脚本直接用于论文图表
- **支持工具**: 19个支持和分析脚本
- **完整文档**: `README.md` 提供详细说明

## 🎯 维护建议

1. **优先维护**: 论文级脚本必须保持100%可用性
2. **定期测试**: 每次数据更新后测试核心脚本
3. **版本控制**: 重要修改前备份脚本
4. **命名规范**: 保持直观的命名约定
5. **依赖管理**: 明确脚本间的依赖关系

## ⚠️ 重要提醒

- 所有脚本从项目根目录运行
- 轨迹提取脚本需要大内存
- 图表生成需要完整的Python科学计算环境
- 定期清理临时文件和过时脚本