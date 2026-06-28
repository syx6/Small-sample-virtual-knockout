# 虚拟敲除方法背景、主流路线和当前方法选择

## 1. 什么是虚拟敲除

虚拟敲除指的是：

```text
在不真正做 CRISPR/siRNA/药物扰动实验的情况下，
用计算模型预测某个基因被敲除后，细胞状态会如何变化。
```

在单细胞场景里，虚拟敲除通常希望回答：

- 敲除某个基因后，细胞会不会从一个状态转到另一个状态？
- 哪些 pathway/program/protein 会升高或降低？
- 哪些细胞亚群对这个 KO 更敏感？
- 单基因 KO 和双基因 KO 的效应是否不同？
- 能不能用已有 perturbation 数据预测新的 KO？

简单说：

```text
真实敲除 = 在实验里真的敲基因
虚拟敲除 = 在计算模型里模拟敲基因
```

## 2. 为什么要做虚拟敲除

真实 CRISPR/perturbation 实验很有价值，但成本高、周期长，而且不能无限组合。

虚拟敲除的价值在于：

1. **降低实验成本**

   可以先用模型筛选候选基因，再挑最有希望的做真实实验。

2. **支持大规模假设生成**

   真实实验很难一次做所有基因、所有双基因组合、所有细胞类型。虚拟 KO 可以先做计算优先级排序。

3. **解释细胞状态变化**

   不只看某个基因表达变不变，而是看 pathway/program/protein 层面的状态迁移。

4. **辅助多基因组合设计**

   双基因 KO 数量随基因数平方增长，真实实验成本很高。虚拟 KO 可以帮助缩小组合空间。

5. **迁移到普通 10X 数据**

   大多数实验室有普通 scRNA-seq 数据，但没有 perturb-seq。理想情况下，可以从公开 perturbation reference 学到 KO 方向，再应用到自己的普通 10X 细胞群。

## 3. 虚拟敲除有哪些主流方法

### 3.1 差异表达 / pathway enrichment 型

最简单路线是：

```text
真实 KO vs control
-> 差异表达基因
-> pathway enrichment
```

优点：

- 简单、可解释。
- 对小数据也能做。
- 生物学读者容易理解。

缺点：

- 只能分析已经做过的 KO。
- 不能很好预测未见 KO。
- 通常不是 cell-level 生成。

### 3.2 单基因效应相加模型

双基因 KO 常见 baseline：

```text
effect(A+B) = effect(A) + effect(B)
```

优点：

- 简单。
- 对部分组合有效。
- 可以作为双敲 baseline。

缺点：

- 不能表达非线性组合效应。
- 对 MAPK/TGFB、Pro-growth 等上下文依赖程序容易失败。

### 3.3 回归模型 / latent model

例如 ridge regression、PLS、linear model、latent factor model：

```text
KO gene features -> cell-state delta
```

优点：

- 小样本下相对稳。
- 可解释。
- 容易接系统先验。

缺点：

- 表达复杂非线性能力有限。
- 对细胞分布形状模拟较弱。

### 3.4 基于系统先验的模型

把 KO gene 映射到：

- Reactome pathway
- MSigDB gene set
- TF-target network
- PPI network

再预测状态变化。

优点：

- 更可解释。
- 对 unseen gene 或 gene pair 有一定泛化能力。
- 适合小样本。

缺点：

- 依赖先验质量。
- 先验缺失的机制难以预测。
- 对强非线性组合仍然有限。

### 3.5 VAE / cVAE / scGen 类方法

典型思想：

```text
把细胞编码到 latent space
学习 control -> perturbation 的 latent shift
再解码成 perturbation cell
```

优点：

- 可以做 cell-level 生成。
- 能模拟一定分布变化。
- 在数据量足够时效果好。

缺点：

- 小样本容易不稳。
- 需要较多 perturbation 和细胞数。
- 解释性较弱。

### 3.6 diffusion / flow matching / conditional generative model

更复杂的生成模型：

```text
condition = KO gene / cell state / modality
model learns perturbation distribution
```

优点：

- 理论上能生成复杂分布。
- 可以建模不确定性和多峰分布。

缺点：

- 训练数据要求高。
- 计算资源要求高。
- 小样本下容易学噪声或方向跑偏。
- 如果缺少 hard constraint，生成结果可能看起来像细胞，但不是正确 KO 方向。

## 4. 为什么我们现在选择 residual/PLS + 系统先验

我们的目标不是做最大规模预训练模型，而是解决：

```text
小样本、多模态、可解释、多基因 KO 的虚拟敲除问题。
```

因此我们选择：

```text
prior-constrained residual/PLS virtual KO
```

核心公式：

```text
virtual KO cell = input/control cell + predicted KO delta
```

`predicted KO delta` 由系统先验和训练 KO 学到。

选择这个方法的原因：

1. **小样本更稳**

   我们实际比较过 VAE、diffusion、flow matching，当前小样本下都不如 residual baseline。

2. **方向可控**

   虚拟 KO 最重要的是方向。先把 KO delta 方向约束住，比自由生成更安全。

3. **可解释**

   输出是 pathway/program/protein 的变化，而不是难解释的 latent vector。

4. **容易接多模态**

   RNA 可以转 pathway/program score；protein 可以直接进入 phenotype score；ATAC 可以进入 regulatory score。

5. **适合普通 10X 扩展**

   普通 10X 没有 KO 标签时，可以先转成相同 state space，再应用 reference perturbation model。

## 5. 我们的方法具体为了解决什么问题

### 5.1 小样本问题

很多实验室没有大规模 perturb-seq 数据，也没有足够算力训练大型生成模型。

我们的方法把问题压缩到：

```text
低维、可解释的 pathway/program/protein state
```

这样小样本下更可行。

### 5.2 多模态问题

真实细胞状态不只由 RNA 决定。

我们希望支持：

```text
RNA -> pathway/program score
protein/ADT -> surface phenotype score
ATAC/gene activity -> regulatory score
```

当前实际跑通：

- RNA-only
- RNA + ADT protein

下一步需要补：

- RNA + ATAC / multiome 示例

### 5.3 普通 10X 应用问题

大多数用户只有普通 10X scRNA-seq，没有 perturbation 标签。

因此我们新增了 reference model 流程：

```text
train-reference:
  从 perturbation 数据训练 KO delta model

apply-reference:
  把 KO delta model 应用到普通 10X 细胞
```

这一步解决的是：

```text
没有 KO 标签的数据也可以获得 virtual KO prediction
```

但要注意：

```text
没有真实 KO，就不能在这个数据集内部评估准确性。
```

### 5.4 双基因 KO 问题

双基因 KO 不是简单相加。

当前方法支持：

```text
--target-kos CEBPB+CEBPA
```

目前已经新增了一个基础交互修正模型：

```text
single-gene additive prediction
+ prior-based gene-gene interaction residual
= optimized double-KO prediction
```

在 Norman 52 个双敲组合上，平均 R2 从 `0.008` 提升到 `0.617`，说明显式交互特征对双敲很有帮助。

下一步需要把它进一步集成到主接口，并增强：

- gene-gene interaction feature
- pathway interaction prior
- constrained nonlinear correction

## 6. 当前方法优势

| 优势 | 说明 |
|---|---|
| 小样本友好 | 不依赖大规模预训练 |
| 可解释 | 输出 pathway/program/protein 变化 |
| 方向受约束 | 不让生成模型自由跑偏 |
| 支持单敲/双敲输入 | `STAT1` 或 `STAT1+JAK2` |
| 可接多模态 | RNA + protein 已跑通，ATAC 接口预留 |
| 可接普通 10X | 通过 reference model 应用到无标签细胞 |
| 可视化清楚 | heatmap、UMAP、AUC/ROC、PCA |

## 7. 当前方法劣势

| 劣势 | 说明 |
|---|---|
| 分布形状模拟有限 | residual/PLS 主要移动状态，不擅长复杂多峰分布 |
| 双敲非线性仍需增强 | interaction residual 已改善，但还没完全主接口化 |
| 依赖先验 | Reactome/MSigDB/TF/PPI 缺失的机制难预测 |
| 普通 10X 无法内部评估 | 没有真实 KO 标签，就没有 ground truth |
| AUC 需要谨慎解释 | feature 数少时 AUC 容易过高 |
| ATAC/multiome 还需真实示例 | 当前主要多模态证据来自 RNA + ADT protein |

## 8. 当前证据总结

### 单敲

Papalexi `STAT1`:

| 指标 | 值 |
|---|---:|
| mean distribution improvement | 0.144 |
| improved fraction | 0.821 |
| direction cosine | 0.879 |
| AUC | 0.878 |

说明单敲方向和强响应排序较好。

### 双敲

Norman `CEBPB+CEBPA`:

| 指标 | 值 |
|---|---:|
| mean distribution improvement | 0.241 |
| improved fraction | 0.600 |
| direction cosine | 0.960 |
| AUC | 1.000 |

说明部分双敲可以预测得很好。

但 52 个 Norman 双敲组合整体：

| 指标 | 值 |
|---|---:|
| mean MAE | 0.131 |
| mean R2 | 0.274 |
| mean ROC-AUC | 0.682 |
| mean PR-AUC | 0.704 |

新增 interaction residual 后，双敲整体明显改善：

| 指标 | single-gene additive | interaction residual |
|---|---:|---:|
| mean MAE | 0.150 | 0.076 |
| mean R2 | 0.008 | 0.617 |
| mean ROC-AUC | 0.707 | 0.894 |
| mean PR-AUC | 0.682 | 0.845 |

说明强约束基础上的交互修正比直接自由生成更适合当前双敲优化。

## 9. 最终定位

当前方法最准确的定位是：

```text
一个面向小样本、多模态 perturbation 数据的
prior-constrained, pathway/program-level, cell-level virtual knockout framework。
```

它最适合：

- 有少量 perturbation 数据
- 想做单基因或双基因 KO 初筛
- 关注 pathway/program/protein 层面变化
- 希望把 reference perturbation model 应用到普通 10X 数据

它不适合被过度宣传为：

- 已经完全解决所有普通 10X 无标签虚拟 KO
- 已经完全解决双敲非线性
- 已经能像大型 diffusion 模型一样模拟复杂细胞分布

## 10. 下一步开发重点

1. **完善 reference model 流程**

   当前已新增 `train-reference` 和 `apply-reference`，下一步要增加模型版本、批量 KO、cell type 分层输出。

2. **真正 multiome 示例**

   加入 RNA + ATAC 或 RNA + protein + ATAC 公开数据。

3. **双敲交互增强**

   基础版 interaction residual 已完成。下一步要集成进主 CLI/reference model，并支持批量双敲。

4. **MAPK/TGFB 专项修正**

   针对失败 program 做先验增强和非线性校正。

5. **普通 10X 应用报告**

   对无标签数据输出 prediction-only report，明确不报告真实准确率。
