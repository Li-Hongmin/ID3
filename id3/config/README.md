# ID3 Configuration Module

## 📋 概述

配置模块管理项目的全局配置、常量定义和序列模板，确保配置的统一性和可维护性。

## 📁 模块内容

- **`id3_config.json`** - 主配置文件
- **`sequence_constants.py`** - 序列相关常量
- **`unified_config.py`** - 统一配置接口
- **`utr_loader.py`** - UTR模板加载器

## 🔧 配置管理

### 主要配置项
- 模型参数设置
- 序列长度限制
- 优化算法参数
- 输出格式配置

### 使用方法
```python
from id3.config import get_config, SEQUENCE_CONSTANTS

# 获取配置
config = get_config()

# 使用常量
max_length = SEQUENCE_CONSTANTS.MAX_SEQUENCE_LENGTH
```

## ⚠️ 维护原则

- 配置变更需更新相关文档
- 兼容性变更需要版本标记
- 敏感配置使用环境变量

---
**最后更新**: 2025-08-17