# 当前方法支持范围和限制说明

## 1. 先说结论

你指出的问题是对的。当前版本应该这样定位：

```text
已经跑通：
1. 原始 scRNA-seq / h5ad 输入
2. RNA -> pathway/program score 自动转换
3. Papalexi RNA + ADT protein 的多模态示例
4. 有 KO 标签的 perturbation 数据上的单敲/双敲评估
5. 普通 10X 无 KO 标签数据的状态转换

还没有完全解决：
1. 普通 10X 无 KO 标签数据上的真实准确性评估
2. 普通 10X + 外部 reference perturbation model 的批量化、分 cell type 和模型版本管理
3. 双敲交互模型的主接口集成和跨数据验证
4. MAPK_TGFB 这类 program 的进一步稳定预测
5. 自由 diffusion/VAE 生成复杂单细胞分布
```

所以当前方法不是“已经能处理所有普通 10X 并自动给出真实准确率”的最终软件，而是一个小样本、强约束、可解释的虚拟 KO 框架原型。

## 2. 当前模型学过多少 KO，没学过的基因怎么预测

这里必须区分两个概念：

```text
接入过的真实扰动数据覆盖范围
≠
某一次模型训练实际学过的 KO 基因数
```

当前本地接入和整理过的 perturbation 数据大致如下：

| 数据集 | 模态 | 真实扰动标签 |
|---|---|---:|
| Papalexi ECCITE-seq | RNA + ADT protein | 25 个单基因 KO |
| Norman Perturb-seq | RNA | 25 个单基因 KO + 52 个双基因 KO |
| HMPCITE-seq | RNA + ADT + guide-derived labels | 11 个基因组成的 66 个单敲/双敲条件 |
| Liscovitch ATAC / scPerturb | ATAC gene activity + chromVAR + peak features | 21 个单基因 KO |
| Datlinger CRISPR | RNA | 约 20 个基因，40 个 guide/扰动标签 |
| Dixit TF perturbation | RNA | 10 个 TF KO |

把这些数据集合并看，当前接入过的是百级别 KO/perturbation gene 覆盖，大约 129 个基因名。这个数字只是“已经接入过的数据覆盖范围”，不是说某一个最终模型已经在一次训练里学遍了 129 个基因。不同数据集的物种、命名大小写、状态空间和模态不同，不能简单粗暴地当作同一个全基因组预训练模型。

### 2.1 单次训练实际学什么

一次模型训练通常发生在一个具体 reference 数据集里。例如：

```text
Papalexi:
  用部分真实单基因 KO 学 KO direction
  再预测留出的 STAT1/JAK2/IFNGR2/IRF1 等 KO

Norman:
  用单基因 KO 学方向
  再预测双基因 KO，或评估 seen/unseen gene pair

HMPCITE:
  用 RNA + ADT + GDO 标签评估真实多模态 double-KO

ATAC:
  用 gene activity / chromVAR / peak features 测试 regulatory-state KO prediction
```

因此当前方法不是“记住所有基因的 KO 效果”，而是学习：

```text
KO gene 的系统先验特征
-> KO 后 pathway / protein / ATAC 状态变化方向
```

### 2.2 没学过的基因如何预测

如果目标基因在训练 KO 中没出现过，模型会把它转成系统先验向量：

```text
target gene
-> Reactome/MSigDB pathway membership
-> TF-target network membership
-> PPI network neighborhood
-> motif / chromVAR / ATAC regulatory prior
-> KO prior vector
```

模型已经从训练 KO 学过：

```text
KO prior vector -> state delta
```

所以未见基因的虚拟 KO 是：

```text
unseen gene 的网络位置
-> 预测它可能影响哪些 pathway/protein/regulatory state
-> 得到 predicted KO delta
-> control cell + predicted KO delta
-> virtual KO cell
```

这属于 prior-based zero-shot / few-shot 外推。它不是凭空预测，而是借助已知生物网络和已见 KO 的响应模式。

### 2.3 未见基因什么时候可信，什么时候不可信

相对可信的情况：

```text
目标基因出现在 Reactome/MSigDB/TF-target/PPI/motif 先验中
目标基因与训练过的 KO 基因处在相似 pathway 或 network module
当前细胞类型中该基因或其下游通路有表达/活性
多模态信息提供了一致信号，例如 RNA + protein 或 RNA + ATAC 都支持同一方向
```

风险较高的情况：

```text
目标基因没有真实 KO 训练样本
目标基因在系统先验中覆盖很弱
目标基因在该细胞类型中表达很低或调控关系未知
KO 效应主要来自当前先验没有覆盖的机制
双基因 KO 里两个基因都没见过，且没有相似组合可参考
```

因此对 unseen gene 应该输出置信度提示。更准确的对外说法是：

```text
模型可以对未见基因做基于系统先验的外推，
适合候选筛选和方向判断；
但如果目标基因缺少先验覆盖或训练相似性，
结果应标记为低置信度，不能当作已验证准确预测。
```

### 2.4 seen gene / unseen gene 的定义

| 术语 | 含义 |
|---|---|
| seen gene | 该基因在训练集中出现过真实单基因扰动 |
| unseen gene | 该基因没有出现在训练单基因扰动里，只能靠系统先验外推 |
| seen combo | 双敲组合中的两个基因都有单基因训练信息 |
| unseen combo | 双敲组合中至少一个基因没有单基因训练信息 |

双敲预测可以分三类看：

| 双敲类型 | 可靠性 |
|---|---|
| 两个基因都 seen | 相对最可靠 |
| 一个 seen、一个 unseen | 中等，取决于 unseen gene 的先验覆盖 |
| 两个基因都 unseen | 风险较高，只适合初筛 |

## 3. 当前输入到底有没有多模态

有，但要分清楚“实际示例”和“接口能力”。

### 3.1 已经实际跑通的多模态

Papalexi ECCITE-seq 示例是多模态输入：

```text
RNA expression matrix: adata.X
ADT protein matrix: adata.obsm["protein"]
KO label: adata.obs["ko_target"]
```

运行命令中用了：

```powershell
--protein-obsm protein
```

因此这个示例不是 RNA-only。它的内部状态是：

```text
RNA-derived pathway score + ADT protein score
```

对应结果目录：

```text
results/software_interface_single_gene_demo
results/software_interface_raw_papalexi
```

### 3.2 当前主要结果里哪些是 RNA-only

下面这些主要是 RNA 或 RNA-derived program：

- Norman Perturb-seq: RNA perturb-seq -> gene program score
- Datlinger CRISPR RNA: RNA-only -> pathway/program score
- Dixit Perturb-seq RNA: RNA-only -> pathway/program score

所以现在的结果并不是全都多模态。真正多模态的主要是 Papalexi RNA + ADT protein。

### 3.3 ATAC / multiome 目前是什么状态

接口层面预留了 ATAC/gene activity 输入，但还没有完成一个正式 multiome 示例。

也就是说：

```text
protein/ADT: 已经实际跑通
RNA-only: 已经实际跑通
ATAC/multiome: 接口预留，仍需要真实数据示例验证
```

## 4. 普通 10X 没有 KO 标签时支持到什么程度

普通 10X 无 KO 标签时，目前支持的是“状态转换”，不是完整评估。

可以做：

```powershell
python -m vkx.cli score \
  --input-h5ad your_10x_data.h5ad \
  --prior-dir data/priors \
  --out-dir results/your_10x_state_scores
```

这会生成：

```text
derived_state_scores.csv
derived_state_manifest.csv
```

含义是：

```text
普通 10X scRNA-seq
-> 自动转换成 pathway/program cell-state table
```

不能做：

```text
不能在没有真实 KO 标签的数据内部计算：
- true vs virtual KO heatmap
- AUC
- UMAP improvement
- distribution improvement
```

原因很简单：没有真实 KO 细胞，就没有 ground truth。

## 5. 普通 10X 怎么真正用于虚拟 KO

合理路线应该是两步：

### 5.1 先从 perturbation/reference 数据训练 KO 方向

例如：

```text
Papalexi / Norman / Datlinger / Dixit / 其他 CRISPR perturbation 数据
-> 学到 KO gene -> pathway/program delta
```

### 5.2 再把这个 KO delta 应用到普通 10X 细胞

```text
普通 10X control-like cells
+ reference model 预测出的 KO delta
= 普通 10X 上的 virtual KO cells
```

这一步已经有基础接口：

```text
train-reference
apply-reference
```

也就是保存一个 reference KO model，再加载它应用到没有 KO 标签的普通 10X 数据。

当前示例：

```text
results/reference_models/papalexi_rna_protein_reference.pkl
results/reference_apply_stat1_demo
```

仍然需要继续增强：

- 批量 KO 输出
- 按 cell type 分层应用和可视化
- reference model 版本记录
- 不同数据集之间的 feature scaling / batch 校正

## 6. 双敲非线性效应

双敲已经有了新的 interaction residual 优化，但还没有完全产品化到主接口。

当前支持双敲输入：

```powershell
--target-kos CEBPB+CEBPA
```

并且 `CEBPB+CEBPA` 示例效果很好：

| 指标 | 值 |
|---|---:|
| mean distribution improvement | 0.241 |
| improved fraction | 0.600 |
| direction cosine | 0.960 |
| AUC | 1.000 |

原始 system prior / additive 版本的 52 个 Norman 双基因组合整体评估说明：

- `all_genes_seen` 比 `has_unseen_gene` 明显更好。
- 双敲不是简单单基因相加。
- MAPK/TGFB 和 Pro-growth 这类 program 有明显非线性。

新增 interaction residual 后，双敲整体明显改善：

| model | subset | mean R2 | mean ROC-AUC |
|---|---|---:|---:|
| single_gene_additive | all_combos | 0.008 | 0.707 |
| interaction_residual | all_combos | 0.617 | 0.894 |
| single_gene_additive | has_unseen_gene | -0.101 | 0.638 |
| interaction_residual | has_unseen_gene | 0.697 | 0.895 |

所以当前结论应该是：

```text
双敲可以支持和展示；interaction residual 已显著改善非线性组合效应；
下一步是把它集成进主 CLI/reference model，并做跨数据验证。
```

## 7. MAPK_TGFB 为什么困难

Norman 双敲整体评估中，`MAPK_TGFB` 表现不稳定：

| subset | R2 | ROC-AUC |
|---|---:|---:|
| all_combos | -0.139 | 0.536 |
| all_genes_seen | -0.608 | 0.909 |
| has_unseen_gene | -0.091 | 0.462 |

解释：

- R2 差，说明连续幅度预测不好。
- AUC 在 seen combo 中还可以，但 unseen combo 变差。
- 这说明 MAPK/TGFB 的强响应有时能识别，但具体变化幅度和跨组合泛化不稳定。

可能原因：

1. MAPK/TGFB 是高度上下文依赖的信号通路。
2. 双基因组合可能有非加性调控。
3. 现有 PLS/residual 模型偏线性，难以表达复杂交互。
4. 小样本下直接上自由深度模型又容易过拟合或方向跑偏。

下一步需要的是“受约束的交互模型”，不是直接让 diffusion/VAE 自由生成。

## 8. 当前生成模型定位

当前最终采用的是：

```text
hard-constrained residual / PLS baseline
```

不是自由 diffusion/VAE/flow matching。

原因是我们已经比较过：

| 模型 | mean distribution improvement |
|---|---:|
| Residual baseline | 0.174 |
| Conditional VAE | -0.189 |
| Diffusion | -0.209 |
| Flow matching | -0.577 |
| Guided Diffusion | -0.769 |

在当前小样本设置下，自由生成模型明显不如 residual baseline。

所以现在的选择是有意的：

```text
先让 residual/PLS 决定 KO 方向；
生成模型最多只能学习方向附近的不确定性；
不能让生成模型自由决定 KO 往哪里走。
```

局限是：

- 复杂单细胞分布形状模拟能力有限。
- 细胞亚群比例变化模拟不足。
- 对强非线性双敲仍然不够。

## 9. AUC 怎么解释

AUC 是有用的，但要谨慎。

它回答的是：

```text
模型能不能把强响应 feature 排在前面？
```

它不等于：

```text
模型已经精确预测了所有 feature 的变化幅度。
```

尤其在双敲 `CEBPB+CEBPA` 示例里，只有 5 个 gene program：

```text
AUC = 1.000
```

这个结果只能说明：

```text
在这 5 个 program 里，强响应排序是对的。
```

不能解释成大规模泛化能力已经完美。

## 10. 当前应该怎么对外描述

比较准确的说法：

```text
这是一个面向小样本 perturbation 数据的 prior-constrained virtual KO 框架。
它支持原始 scRNA-seq 输入，并能合并 protein/ADT 等多模态信息；
普通 10X 无 KO 标签数据目前支持状态转换，也已经有基础版 reference model 应用接口；但没有真实 KO 标签时仍不能在该数据内部评估准确性。
当前方法对单敲更稳；双敲在 interaction residual 后明显改善，但仍需主接口集成和跨数据验证。
```

不应该说：

```text
已经完全解决普通 10X 无标签虚拟 KO。
已经完全解决多模态 multiome。
已经完全产品化解决所有双敲非线性。
AUC 高就说明模型全面准确。
```

## 11. 下一步最应该做什么

优先级建议：

1. **reference model 接口增强**

   基础版 `train-reference` 和 `apply-reference` 已完成。下一步增强模型版本、批量 KO、cell type 分层输出和跨数据校正。

2. **真正 multiome 示例**

   至少加入一个 RNA + ATAC 或 RNA + protein + ATAC 的公开数据示例。

3. **双敲交互模型**

   在 residual/PLS hard constraint 上加入 gene-gene interaction features，而不是自由深度生成。

4. **MAPK/TGFB 专项诊断**

   单独分析 MAPK/TGFB 失败组合，判断是先验缺失、状态表示不够，还是模型线性假设不够。
