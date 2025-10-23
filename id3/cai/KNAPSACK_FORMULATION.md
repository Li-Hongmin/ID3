# CAI优化的背包问题形式化

## 1. 问题定义

### 1.1 原始CAI优化问题

给定氨基酸序列 **S** = (a₁, a₂, ..., aₙ)，每个氨基酸 aᵢ 可由密码子集合 Cᵢ = {cᵢ₁, cᵢ₂, ..., cᵢₖ} 编码。

**目标函数：**
$$\max \prod_{i=1}^{n} P(c_i) \quad \text{where } c_i \in C_i$$

**约束条件：**
$$\left(\prod_{i=1}^{n} w(c_i)\right)^{\frac{1}{n}} \geq \tau_{CAI}$$

其中：
- P(cᵢ): 密码子 cᵢ 的选择概率
- w(cᵢ): 密码子 cᵢ 的CAI权重，w(cᵢ) ∈ (0, 1]
- τ_CAI: 目标CAI阈值

### 1.2 对数域转换

取对数将乘积转为求和：

**目标函数：**
$$\max \sum_{i=1}^{n} \log P(c_i)$$

**约束条件：**
$$\frac{1}{n} \sum_{i=1}^{n} \log w(c_i) \geq \log \tau_{CAI}$$

## 2. 背包问题映射

### 2.1 Retreat策略形式化

从最优CAI配置开始（∀i: cᵢ = argmax_{c∈Cᵢ} w(c)），通过替换降低CAI以提高概率。

定义替换变量：
- xᵢⱼ ∈ {0,1}: 是否在位置 i 选择密码子 j
- 约束：∑ⱼ xᵢⱼ = 1, ∀i

### 2.2 标准背包形式

令 c*ᵢ = argmax_{c∈Cᵢ} w(c) 为位置 i 的最优CAI密码子。

**价值函数（概率增益）：**
$$v_{ij} = \log P(c_{ij}) - \log P(c_i^*)$$

**权重函数（CAI损失）：**
$$w_{ij} = \log w(c_i^*) - \log w(c_{ij})$$

**背包问题：**
$$\max \sum_{i,j} v_{ij} \cdot x_{ij}$$

$$\text{s.t.} \quad \sum_{i,j} w_{ij} \cdot x_{ij} \leq n(\log 1.0 - \log \tau_{CAI})$$

$$\sum_j x_{ij} = 1, \quad \forall i$$

## 3. 贪心算法分析

### 3.1 性价比定义

定义性价比（价值密度）：
$$\rho_{ij} = \frac{v_{ij}}{w_{ij}} = \frac{\log P(c_{ij}) - \log P(c_i^*)}{\log w(c_i^*) - \log w(c_{ij})}$$

### 3.2 Retreat贪心策略

**算法步骤：**
1. 初始化：X⁽⁰⁾ = {c*₁, c*₂, ..., c*ₙ}
2. 迭代 t = 1, 2, ...：
   - 计算当前CAI：CAI⁽ᵗ⁾ = (∏ᵢ w(cᵢ⁽ᵗ⁾))^(1/n)
   - 若 CAI⁽ᵗ⁾ ≤ τ_CAI，终止
   - 选择最优替换：(i*, j*) = argmax_{i,j} ρᵢⱼ
   - 更新：cᵢ*⁽ᵗ⁺¹⁾ = cᵢ*ⱼ*

### 3.3 时间复杂度

- 每次迭代：O(nk) 计算所有性价比
- 总迭代次数：O(m)，m ≤ n
- 总复杂度：O(n²k)

## 4. 理论性质

### 4.1 分数背包性质

**定理1（贪心最优性）：** 
对于分数背包问题，按价值密度 ρ 降序选择物品得到最优解。

**证明概要：**
设最优解中存在物品 i, j 满足 ρᵢ > ρⱼ，且 xⱼ > 0 而 xᵢ < 1。
交换量 ε = min(xⱼ, 1-xᵢ) 可严格改进目标值，矛盾。

### 4.2 CAI约束的特殊结构

**引理1（单调性）：**
$$\frac{\partial \text{CAI}}{\partial w(c_i)} > 0, \quad \forall i$$

**引理2（凸性）：**
log CAI 关于 log w(cᵢ) 是线性的：
$$\log \text{CAI} = \frac{1}{n} \sum_{i=1}^{n} \log w(c_i)$$

### 4.3 近似比分析

**定理2（近似保证）：**
Retreat算法的解满足：
$$\text{OPT} - \text{ALG} \leq v_{\max}$$

其中 v_max = max_{i,j} vᵢⱼ 是单个替换的最大价值。

## 5. 增量CAI计算

### 5.1 O(1)更新公式

替换位置 i 的密码子从 c_old 到 c_new：

$$\text{CAI}_{\text{new}} = \text{CAI}_{\text{old}} \cdot \left(\frac{w(c_{\text{new}})}{w(c_{\text{old}})}\right)^{\frac{1}{n}}$$

### 5.2 数值稳定性

对数域计算避免数值下溢：
$$\log \text{CAI}_{\text{new}} = \log \text{CAI}_{\text{old}} + \frac{1}{n}(\log w(c_{\text{new}}) - \log w(c_{\text{old}}))$$

## 6. 扩展形式

### 6.1 多目标优化

$$\min_{\mathbf{x}} \left[ -\sum_{i,j} \log P(c_{ij}) x_{ij}, \quad d(\mathbf{x}, \mathbf{x}^{(0)}) \right]$$

其中 d(·,·) 是序列距离度量。

### 6.2 带正则化的形式

$$\max \sum_{i,j} \log P(c_{ij}) x_{ij} - \lambda \sum_{i,j} x_{ij} \log x_{ij}$$

熵正则项促进密码子使用多样性。

## 7. 计算复杂度对比

| 算法 | 时间复杂度 | 空间复杂度 | 解质量 |
|-----|-----------|-----------|--------|
| Retreat贪心 | O(n²k) | O(nk) | 1-近似* |
| 标准DP | O(nWk) | O(nW) | 最优 |
| FPTAS | O(n²k/ε) | O(n/ε) | (1+ε)-近似 |

*对于分数背包松弛

## 8. 结论

CAI优化问题本质是一个特殊的背包问题：
- **目标**：最大化概率收益 ∑ vᵢⱼxᵢⱼ
- **约束**：CAI预算限制 ∑ wᵢⱼxᵢⱼ ≤ B
- **特性**：具有分数背包的连续松弛性质
- **算法**：贪心策略提供高质量近似解

Retreat算法利用了问题的特殊结构（对数线性、单调性），在实践中表现优异。