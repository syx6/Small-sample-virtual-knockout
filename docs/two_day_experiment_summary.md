# 小样本多模态虚拟敲除方法：两天实验总结

## 1. 这两天我们到底在做什么

目标是开发一个适合小样本、多模态、可解释的虚拟基因敲除方法。这个方法希望解决三个问题：

- 支持单基因敲除和多基因敲除。
- 不依赖大规模预训练和高算力。
- 输出结果不能只是抽象指标，而要能让人看懂：敲了什么基因，细胞状态往哪里变，哪些 pathway/program/protein 发生变化，虚拟结果和真实 KO 有多像。

最后形成的主线方法是：

```text
输入细胞状态
  RNA-only: pathway/program score
  multimodal/multiome: pathway score + protein/ATAC score

输入 KO 条件
  KO gene / KO gene pair / KO gene set
  + Reactome / MSigDB / TF-target / PPI 等系统网络先验

预测 KO 方向
  residual / PLS baseline 学习 control -> KO 的平均状态变化方向

生成虚拟 KO 单细胞
  virtual KO cell = control cell + constrained KO delta

幅度校准
  对 RNA-only pathway/program score 做非负倍率校准，让变化幅度更接近真实 KO
```

核心思想不是从零训练一个自由生成模型，而是先用稳定的、可解释的 baseline 固定 KO 方向，再在这个方向附近生成单细胞状态。对小样本来说，这比直接训练 VAE、flow matching 或 diffusion 更稳。

## 2. 使用的数据和输入输出

这两天实际测试了多种 perturbation 数据：

| 数据集 | 类型 | 输入状态表示 | 测试内容 |
|---|---|---|---|
| Papalexi ECCITE-seq | RNA + ADT protein | pathway score + protein score | 单基因 KO，多模态虚拟敲除 |
| Norman Perturb-seq | RNA perturb-seq | gene program score | 多基因 KO |
| Datlinger CRISPR RNA | RNA CRISPR screen | pathway/program score | RNA-only 单基因 KO |
| Dixit Perturb-seq RNA | RNA perturb-seq | pathway/program score | RNA-only TF KO |

方法输入：

- 单细胞表达矩阵，或已经整理好的 pathway/program/protein/ATAC score。
- 每个细胞对应的扰动标签，例如 `control`, `STAT1`, `CEBPB+CEBPA`。
- 系统先验网络，例如 Reactome、MSigDB、TF-target、PPI。
- 需要预测的 KO 目标，可以是单基因，也可以是多基因组合。

方法输出：

- 虚拟 KO 后的单细胞状态表。
- 每个 KO 的 pathway/program/protein 变化方向。
- 真实 KO vs 虚拟 KO 的 heatmap。
- KO 前后 UMAP 图。
- 指标表，包括分布距离改善、方向一致性、幅度误差。

## 3. 第一阶段：pathway / protein / gene-prior baseline

最早的实验是把 Papalexi 数据整理成 pathway score 和 protein score，然后尝试用 gene prior / pathway prior 去预测 KO 后的状态变化。

这一阶段确认了两件事：

- 单纯看原始 gene expression 太难解释，换成 pathway score 更适合展示 KO 效果。
- protein 模态能提供 RNA-only 没有的信息，尤其在免疫检查点、干扰素响应相关 KO 上更有帮助。

因此后面我们把状态表示统一成更解释性的形式：

- RNA-only 不再用 SVD，而是 pathway/program score。
- RNA + protein 用 pathway score + protein score。
- 将来 multiome 可以对应 pathway score + ATAC regulatory score。

## 4. 第二阶段：AUC、R2、heatmap 和可视化修正

中间我们调试过 AUC 和 R2 的展示方式。

结论是：

- AUC 应该用 ROC curve 展示，而不是只用柱状图。
- R2 图必须清楚说明比较对象，否则很容易让人误解。
- 对这个任务来说，最直观的可视化不是单独的 R2，而是：
  - 真实 KO change heatmap
  - 虚拟 KO change heatmap
  - prediction error heatmap
  - KO 前后 UMAP

因此后续主图改成了 heatmap + UMAP。

重要图：

- `results/figures/papalexi_roc_curves.png`
- `results/figures/multi_dataset_true_vs_virtual_agreement_heatmap.png`
- `results/figures/multi_dataset_virtual_ko_umap_examples.png`

## 5. 第三阶段：真正进入 cell-level 虚拟敲除

我们从只预测平均 pathway/protein 变化，进一步扩展到单细胞层面的虚拟 KO。

方法变成：

```text
从 control cells 取真实单细胞
+ 预测出来的 KO delta
= 虚拟 KO 单细胞
```

这一步可以直接画出：

- KO 前 control cells 的 UMAP。
- 虚拟 KO 后 cells 的 UMAP。
- 真实 KO cells 的 UMAP。
- 虚拟 KO 是否往真实 KO 云团移动。

单细胞 UMAP 结果：

| 数据 | KO | UMAP centroid improvement |
|---|---:|---:|
| Papalexi single-gene KO | IRF1 | 0.676 |
| Papalexi single-gene KO | STAT1 | 0.297 |
| Papalexi single-gene KO | JAK2 | 0.076 |
| Papalexi single-gene KO | IFNGR2 | -0.396 |
| Norman multi-gene KO | CEBPB+CEBPA | 0.902 |
| Norman multi-gene KO | AHR+KLF1 | 0.280 |
| Norman multi-gene KO | MAPK1+TGFBR2 | 0.093 |
| Norman multi-gene KO | CBL+UBASH3B | 0.054 |

解释：

- 多数 KO 的虚拟细胞状态确实往真实 KO 方向移动。
- 不是所有 KO 都成功，IFNGR2 是明显失败例子。
- 多基因 KO 里 CEBPB+CEBPA 效果很好，说明方法有支持多基因 KO 的潜力。

## 6. 第四阶段：比较 VAE / flow matching / diffusion

我们测试了把 cell-level 生成模型换成更复杂的模型：

| 模型 | mean distribution improvement | improved fraction |
|---|---:|---:|
| Residual baseline | 0.174 | 0.778 |
| Conditional VAE | -0.189 | 0.333 |
| Diffusion | -0.209 | 0.500 |
| Guided Conditional VAE | -0.566 | 0.417 |
| Flow matching | -0.577 | 0.222 |
| Guided Flow matching | -0.766 | 0.389 |
| Guided Diffusion | -0.769 | 0.361 |

这是一个很关键的结果。

直接训练复杂生成模型并没有变好，反而明显变差。原因大概率是：

- 小样本里每个 KO 的真实细胞数有限。
- 条件生成模型需要更多 perturbation、多细胞状态和更强训练信号。
- 如果模型自由度太高，会学到分布噪声，而不是 KO 方向。
- soft guidance 不够，生成模型仍然可以偏离真实 KO 方向。

因此我们确定了后续策略：

```text
residual/PLS baseline 是 hard constraint
生成模型只学方向附近的不确定性范围
不能让生成模型自由决定 KO 往哪里走
```

## 7. 第五阶段：Reactome / MSigDB / TF-target / PPI 系统先验

我们引入系统网络先验之后，比较了几种方向预测方法：

| 模型 | mean distribution improvement | improved fraction |
|---|---:|---:|
| Cross-data PLS | 0.138 | 0.694 |
| PLS | 0.137 | 0.694 |
| KNN-3 | 0.088 | 0.667 |
| Constrained ensemble | 0.088 | 0.722 |
| KNN-7 | 0.070 | 0.667 |
| Zero shrink | 0.000 | 0.000 |

结论：

- PLS 已经是很强的稳定 baseline。
- 跨数据 PLS 略有帮助，但没有大幅超过普通 PLS。
- KNN 和 ensemble 容易过度收缩。
- 系统先验的价值主要是让 KO 条件变得可泛化、可解释，而不是让复杂模型自动变强。

## 8. 第六阶段：hard constraint uncertainty

我们尝试引入更多 perturbation 数据，让生成模型只学习不确定性范围，而不是改变方向。

结果：

| 模型 | mean distribution improvement | improved fraction |
|---|---:|---:|
| hard_mean | 0.061 | 0.694 |
| hard_uncertainty_samples | 0.061 | 0.694 |

外部 perturbation 数据估计出的噪声幅度很小，最后选择的 noise scale 接近 0。

解释：

- 现阶段加入随机不确定性并没有提升平均预测效果。
- 但是它可以作为输出可信区间使用。
- 更合理的定位是：输出 uncertainty band，而不是让噪声参与改变 KO 主方向。

## 9. 第七阶段：多数据集、多模态验证

我们把同一套方法放到四类数据上测试。

未做幅度校准前：

| 数据集 | 输入模态 | 状态表示 | mean distribution improvement | improved fraction |
|---|---|---|---:|---:|
| Papalexi ECCITE-seq | RNA + ADT protein | pathway/protein scores | 0.137 | 0.750 |
| Norman Perturb-seq | RNA perturb-seq | gene program scores | -0.233 | 0.500 |
| Dixit Perturb-seq RNA | RNA CRISPR screen | pathway/program scores | -0.114 | 0.232 |
| Datlinger CRISPR RNA | RNA CRISPR screen | pathway/program scores | -0.196 | 0.179 |

这个结果说明：

- 多模态 Papalexi 表现最好。
- RNA-only 用 pathway/program score 后，方向可解释，但分布距离仍然不够好。
- Norman 多基因 KO 是混合结果，有些组合好，有些组合差。

## 10. 真实 KO vs 虚拟 KO heatmap 结果

方向一致性是判断模型是否“敲对方向”的关键。

| 数据集 | KO | direction cosine | mean abs delta error |
|---|---|---:|---:|
| Papalexi | STAT1 | 0.931 | 0.196 |
| Papalexi | JAK2 | 0.600 | 0.185 |
| Papalexi | IFNGR2 | 0.966 | 0.102 |
| Papalexi | IRF1 | 0.189 | 0.101 |
| Norman | AHR+KLF1 | 0.694 | 0.175 |
| Norman | CEBPB+CEBPA | 0.942 | 0.115 |
| Norman | MAPK1+TGFBR2 | 0.134 | 0.172 |
| Norman | CBL+UBASH3B | 0.712 | 0.225 |
| Datlinger | LAT | 0.877 | 0.052 |
| Datlinger | LCK | 0.290 | 0.212 |
| Datlinger | JUND | 0.863 | 0.044 |
| Datlinger | FOS | 0.944 | 0.048 |
| Dixit | ELF1 | 0.994 | 0.015 |
| Dixit | CREB1 | -0.413 | 0.068 |
| Dixit | ELK1 | 0.969 | 0.016 |
| Dixit | GABPA | 0.936 | 0.024 |

解释：

- 很多 RNA-only KO 的方向其实很好，例如 Dixit 的 ELF1、ELK1、GABPA，Datlinger 的 FOS、JUND、LAT。
- 但是分布距离仍然不好，说明主要问题不是方向，而是变化幅度和单细胞分布形状。
- CREB1、LCK、IRF1、MAPK1+TGFBR2 是明显困难或失败例子。

## 11. 第八阶段：pathway/program score 幅度校准

我们专门做了幅度校准，让 RNA-only 数据不仅方向对，而且变化幅度更接近真实 KO。

校准规则：

```text
virtual KO state = control state + alpha * predicted KO delta
```

其中 alpha 限制在 `0.15 到 2.5`，所以：

- 不能把方向反过来。
- 不能把虚拟 KO 缩成完全没变化。
- 只能调节 KO 效应大小。

校准后结果：

| 数据集 | 校准方式 | 分布改善 before | 分布改善 after | 幅度误差 before | 幅度误差 after |
|---|---|---:|---:|---:|---:|
| Datlinger | feature scale | -0.170 | -0.059 | 0.059 | 0.052 |
| Dixit | global scale | -0.172 | -0.024 | 0.044 | 0.035 |
| Norman | feature scale | -0.766 | -0.643 | 0.180 | 0.171 |
| Papalexi | feature scale | 0.004 | -0.046 | 0.145 | 0.175 |

解释：

- RNA-only 的 Datlinger 和 Dixit 明显改善。
- Norman 也有小幅改善，但仍然受多基因方向错误影响。
- Papalexi 的多模态输入本来已经较稳，强行做 feature-level 幅度校准反而会变差。

因此建议：

```text
RNA-only: pathway/program score + non-negative magnitude calibration
RNA + protein / multiome: 默认不强制校准，作为可选后处理
```

## 12. 方法现在的适用性判断

### 最适合的场景

- 小样本 perturbation 数据。
- 有一部分真实 KO 可以作为训练集。
- 研究者关心 pathway/program/protein 层面的状态变化，而不是逐基因表达的精细还原。
- 希望支持单基因和多基因 KO。
- 数据可以提供多模态信息，或者至少可以从 RNA 转成 pathway/program score。

### 比较有优势的地方

- 可解释：能说清楚哪个 KO 引起哪些 pathway/program/protein 改变。
- 小样本友好：不需要大规模预训练。
- 多模态友好：RNA、protein、ATAC 都可以转成 state score 后进入同一框架。
- 多基因可扩展：KO 条件通过 gene set / network prior 表示，不限于单基因。
- 可视化清楚：heatmap 和 UMAP 能直接展示虚拟 KO 是否接近真实 KO。

### 当前不适合的场景

- 完全没有训练 KO 的 zero-shot 任务。
- 想预测每个基因表达值的精细变化。
- KO 真实效应主要由未知调控机制决定，而不在现有 Reactome/MSigDB/TF-target/PPI 先验里。
- 需要精确模拟细胞亚群比例变化或复杂分布形状。
- 方向本身已经预测错的 KO，幅度校准救不了。

## 13. 为什么小样本可以做，但不能盲目乐观

这两天的实验给了一个比较清楚的答案：

小样本能做，不是因为新算法神奇地替代了数据量，而是因为我们把任务缩小到更稳定、可解释的层面：

- 不直接预测全基因表达矩阵。
- 不让深度生成模型自由学习 KO 分布。
- 把细胞状态压缩成 pathway/program/protein/ATAC score。
- 用系统网络先验约束 KO 条件。
- 用 residual/PLS 固定方向。
- 必要时只校准幅度。

也就是说：

```text
小样本可行 = 多模态信息 + 系统先验 + 低维可解释状态 + 强约束模型
```

不是单靠某一个“小样本算法”就能解决。

## 14. 为什么多模态仍然重要

Papalexi 的表现整体最好，说明多模态确实有价值。

原因是：

- RNA 反映 transcriptional response。
- protein 反映更接近功能表型的状态。
- ATAC 将来可以反映 regulatory accessibility。
- 多模态可以让 KO 后的状态变化不只停留在转录层面。

但多模态不是万能的。更准确的说法是：

```text
模态越多，状态表征越充分；
但样本太少时，模型仍然必须强约束。
```

## 15. 当前推荐的最终方法版本

推荐把方法定义为：

**Prior-constrained pathway/program residual virtual knockout**

流程：

1. 把输入单细胞数据转成可解释状态。
   - RNA-only: pathway/program score。
   - RNA + protein: pathway score + protein score。
   - Multiome: pathway score + ATAC regulatory score。
2. 用训练 KO 计算真实 KO delta。
3. 用 Reactome/MSigDB/TF-target/PPI 构建 KO 条件向量。
4. 用 PLS/residual baseline 预测未见 KO 的 delta。
5. 对 control cells 加上 delta，生成虚拟 KO cells。
6. RNA-only 可加非负幅度校准。
7. 输出 heatmap、UMAP、方向一致性、幅度误差、分布距离。

一句话解释给别人：

```text
我们不是让模型凭空生成 KO 细胞，而是用生物网络先验预测 KO 应该把细胞状态推向哪里，
再把真实 control cells 沿这个方向移动，得到可解释的虚拟 KO 单细胞状态。
```

## 16. 最重要的图和文件

综合报告相关：

- `docs/final_readable_report.md`
- `docs/multi_dataset_virtual_ko_demo.md`
- `docs/multi_dataset_true_vs_virtual_heatmaps.md`
- `docs/pathway_magnitude_calibration.md`

关键图：

- `results/figures/multi_dataset_virtual_ko_effect_summary.png`
- `results/figures/multi_dataset_virtual_ko_umap_examples.png`
- `results/figures/multi_dataset_true_vs_virtual_agreement_heatmap.png`
- `results/figures/pathway_magnitude_calibration_summary.png`
- `results/figures/pathway_magnitude_calibration_delta_error_heatmap.png`
- `results/figures/pathway_magnitude_calibration_heatmap_datlinger_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_dixit_global_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_norman_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_papalexi_feature_scale.png`

关键结果表：

- `results/multi_dataset_virtual_ko_summary.csv`
- `results/multi_dataset_true_vs_virtual_heatmap_summary.csv`
- `results/pathway_magnitude_calibration_summary.csv`
- `results/pathway_magnitude_calibration_factors.csv`
- `results/papalexi_cell_level_deep_generator_summary.csv`

## 17. 下一步建议

最值得继续做的是三件事：

1. **把方法整理成可复用软件接口**

   输入固定为：单细胞矩阵、扰动标签、模态信息、先验网络、待预测 KO。

   输出固定为：虚拟 KO 单细胞表、KO effect heatmap、UMAP、指标表。

2. **扩展 multiome / ATAC score**

   当前 protein 已经验证有价值。下一步应该加入 ATAC regulatory program score，让方法真正覆盖 RNA + protein + ATAC。

3. **更严格评估多基因 KO**

   Norman 结果说明多基因 KO 可行，但组合效应不是简单相加。下一步应该专门比较：

   - 单基因 delta 相加。
   - PLS/system prior 预测组合 delta。
   - 带 pathway interaction 的组合模型。

## 18. 总体结论

这两天的结果支持这个方向继续做，但也说明它应该被定位成一个“强约束、可解释、小样本虚拟敲除方法”，而不是一个自由深度生成模型。

最可靠的结论是：

- 多模态数据确实提升了状态表征质量。
- pathway/program score 比 SVD 更适合这个方法，因为可解释。
- residual/PLS baseline 比 VAE、flow matching、diffusion 更适合当前小样本设置。
- 系统网络先验是必要的，它让 KO 条件可泛化。
- RNA-only 可以用幅度校准改善结果。
- 多基因 KO 有希望，但需要更专门的组合效应建模。

当前方法的合理定位：

```text
一个面向小样本、多模态 perturbation 数据的
prior-constrained, pathway-level, cell-level virtual knockout framework。
```

