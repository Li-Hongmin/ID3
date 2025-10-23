# ID3 Utilities Module

## 📋 概述

工具模块提供项目所需的通用功能函数，包括数据处理、数值计算、序列操作等辅助功能。

## 📁 模块内容

- **`constants.py`** - 项目常量定义
- **`functions.py`** - 通用功能函数
- **`functions_vectorized.py`** - 向量化计算函数
- **`deepraccess_wrapper.py`** - DeepRaccess模型封装
- **`dimension_converter.py`** - 维度转换工具
- **`error_handling.py`** - 错误处理工具
- **`gradient_utils.py`** - 梯度计算工具
- **`memory_utils.py`** - 内存管理工具
- **`numerical_stability.py`** - 数值稳定性工具
- **`serialization.py`** - 序列化工具
- **`utils.py`** - 其他通用工具

## 🔧 主要功能分类

### 序列处理
- 氨基酸与密码子转换
- 序列验证和格式化
- RNA结构分析

### 数值计算
- 向量化操作
- 梯度计算
- 数值稳定性处理

### 系统工具
- 内存管理
- 错误处理
- 数据序列化

## 📚 使用示例

```python
from id3.utils import validate_sequence, convert_to_rna
from id3.utils.gradient_utils import safe_gradient

# 序列验证
is_valid = validate_sequence("ATGAAACCC")

# RNA转换
rna_seq = convert_to_rna(dna_sequence)

# 安全梯度计算
grad = safe_gradient(loss, parameters)
```

## ⚠️ 开发规范

- 新增工具函数需要单元测试
- 保持函数功能单一和职责明确
- 向量化操作优先于循环实现
- 关注内存使用和性能优化

---
**最后更新**: 2025-08-17