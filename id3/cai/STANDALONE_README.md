# Standalone CAI Discretization Module

## 概述 / Overview

This directory contains standalone, highly maintainable CAI (Codon Adaptation Index) optimization modules that follow Unix philosophy - each module does one thing well and can be composed with other tools.

这个目录包含独立的、高度可维护的CAI优化模块，遵循Unix哲学 - 每个模块只做好一件事，可以与其他工具组合使用。

## 核心模块 / Core Modules

### 1. `cai_discretizer.py`
**Purpose**: Convert continuous RNA probabilities to discrete sequences with CAI optimization  
**目的**: 将连续RNA概率转换为具有CAI优化的离散序列

- No dependencies on ID3 framework / 不依赖ID3框架
- Binary search for globally optimal CAI / 二分查找全局最优CAI
- Maintains gradient flow via straight-through estimator / 通过直通估计器保持梯度流
- Clean, documented interface / 清晰的文档化接口

### 2. `binary_search_optimizer.py`
**Purpose**: Core binary search algorithm for CAI optimization  
**目的**: CAI优化的核心二分查找算法

- Mathematical foundation: α ∈ [0,1] where S(α) = (1-α)·π + α·w
- Guarantees optimal trade-off between probability and CAI
- Precision control with configurable tolerance

### 3. `accessibility_optimizer.py`
**Purpose**: Example integration of CAI discretization with accessibility optimization  
**目的**: CAI离散化与可访问性优化的集成示例

- Demonstrates gradient-based optimization
- Uses CAI discretizer for sequence generation
- Mock accessibility model for testing

### 4. `minimal_example.py`
**Purpose**: Minimal working example without any framework dependencies  
**目的**: 不依赖任何框架的最小工作示例

- Pure PyTorch implementation
- Shows essential workflow
- Easy to understand and modify

## 设计原则 / Design Principles

1. **模块化 / Modularity**
   - Each module is self-contained
   - Clear interfaces between modules
   - No hidden dependencies

2. **可维护性 / Maintainability**
   - Comprehensive documentation in code headers
   - Clear purpose and interface descriptions
   - Type hints for all functions

3. **Unix哲学 / Unix Philosophy**
   - Do one thing well
   - Compose with other tools
   - Text/tensor in, text/tensor out

## 使用方法 / Usage

### Basic Discretization
```python
from cai_discretizer import discretize_with_cai
import torch

# Create RNA probabilities
probs = torch.rand(36, 4)  # 12 amino acids * 3
probs = torch.softmax(probs, dim=-1)

# Discretize with CAI target
discrete_seq, metadata = discretize_with_cai(
    probs, 
    amino_sequence="MSKGEELFTGVV",
    cai_target=0.8
)

print(f"CAI achieved: {metadata['mean_cai']:.3f}")
```

### With Accessibility Optimization
```python
from accessibility_optimizer import AccessibilityOptimizer

optimizer = AccessibilityOptimizer(device='cpu')
result = optimizer.optimize(
    amino_sequence="MSKGEELFTGVV",
    utr5="AGATCT" * 10,
    utr3="TAATAA" * 10,
    cai_target=0.8,
    iterations=100
)

print(f"Accessibility: {result['accessibility']:.3f} kcal/mol")
print(f"CAI: {result['cai']:.3f}")
```

### Minimal Example
```python
# Run the minimal example
python minimal_example.py
```

## 接口说明 / Interface Documentation

### `discretize_with_cai()`

**输入 / Input**:
- `probabilities`: torch.Tensor - RNA probability distribution [seq_len, 4]
- `amino_sequence`: str - Target amino acid sequence
- `cai_target`: float - Target CAI value (0.0-1.0)
- `return_scores`: bool - If True, return (sequence_str, cai_score)
- `device`: str - Computation device ('cpu' or 'cuda')

**输出 / Output**:
- If `return_scores=False`: (discrete_tensor, metadata_dict)
- If `return_scores=True`: (sequence_string, cai_score)

### `BinarySearchCAIOptimizer.optimize()`

**输入 / Input**:
- `amino_acid_sequence`: str - Target amino acids
- `codon_probabilities`: dict - Probability distribution over codons

**输出 / Output**:
- `sequence`: str - Optimized DNA sequence
- `cai`: float - Achieved CAI score
- `details`: dict - Optimization details

## 数学原理 / Mathematical Foundation

The binary search finds the minimum α such that:
```
CAI(S(α)) ≥ τ
where S(α) = (1-α)·π + α·w
```

- π: probability-optimal sequence (argmax of probabilities)
- w: CAI-optimal sequence (highest CAI weights)
- τ: target CAI threshold
- α: interpolation parameter ∈ [0,1]

## 性能特征 / Performance Characteristics

- **时间复杂度 / Time Complexity**: O(n·log(1/ε)) where n = sequence length, ε = precision
- **空间复杂度 / Space Complexity**: O(n)
- **收敛保证 / Convergence**: Guaranteed by binary search properties
- **梯度流 / Gradient Flow**: Maintained through straight-through estimator

## 测试 / Testing

### Unit Tests
```bash
# Test CAI discretizer
python -c "from cai_discretizer import discretize_with_cai; print('Import successful')"

# Test binary search optimizer  
python test_cai_integration.py
```

### Integration Test
```bash
# Run accessibility optimization demo
python -m accessibility_optimizer

# Run minimal example
python minimal_example.py
```

## 维护建议 / Maintenance Guidelines

1. **添加功能 / Adding Features**
   - Keep modules independent
   - Document new interfaces clearly
   - Maintain backward compatibility

2. **优化性能 / Optimizing Performance**
   - Profile before optimizing
   - Maintain correctness tests
   - Document performance changes

3. **调试 / Debugging**
   - Use minimal_example.py for isolation
   - Check device compatibility (CPU vs CUDA)
   - Verify amino acid sequence lengths

## 常见问题 / Common Issues

### Device Mismatch
**Problem**: RuntimeError about tensors on different devices  
**Solution**: Ensure all tensors use same device, pass `device='cpu'` or `device='cuda'`

### Sequence Length Mismatch
**Problem**: ValueError about sequence length  
**Solution**: Check that RNA length = amino_acid_length × 3

### CAI Target Not Met
**Problem**: Cannot achieve target CAI  
**Solution**: Lower target or check codon usage table for species

## 未来工作 / Future Work

- [ ] Add real DeepRaccess model integration
- [ ] Support batch processing optimization
- [ ] Add more species CAI tables
- [ ] Implement adaptive CAI targets
- [ ] Add visualization tools

## 作者 / Authors

ID3-DeepRaccess Team

## 许可 / License

See project LICENSE file

---

最后更新 / Last Updated: 2025-01-09