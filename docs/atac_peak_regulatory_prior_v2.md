# ATAC peak-level regulatory prior v2

这一步把 ATAC peak selection 从“按全局方差/可及细胞数筛选”升级为“调控先验综合打分”。

## 1. 为什么要增强 peak-level prior

ATAC peak 很稀疏，单纯按方差或可及细胞数筛选容易选到：

- 很活跃但和目标 KO 无关的 peak；
- 全局开放 peak；
- 技术上稳定但生物解释弱的 peak。

虚拟 KO 更需要回答：

```text
敲掉这个基因后，哪些调控区域可能变化？
这些 peak 是否在 target gene locus、marker peak、KO effect 或 TF/motif 网络里有支持？
```

## 2. 新的 peak 打分

每个 peak 现在会计算 5 个分数：

| 分数 | 含义 |
|---|---|
| `locus_score` | peak 是否靠近 target gene 或其网络邻居，promoter/intronic/distal 有不同权重 |
| `marker_score` | 是否来自 `markerpeak_target` 中与 target gene 相关的 marker peak |
| `ko_effect_score` | target KO 细胞和 control 细胞之间 peak accessibility 的实际差异 |
| `accessibility_score` | peak 在多少细胞中可及，避免选过稀疏 peak |
| `motif_tf_score` | target gene / TF-target / PPI / pathway 网络是否支持该 peak 的 nearest gene |

总分：

```text
total_score =
  0.30 * locus_score
+ 0.23 * KO effect score
+ 0.20 * marker score
+ 0.15 * motif/TF score
+ 0.12 * accessibility score
```

同时保存 `selection_reason`：

```text
target_or_network_locus
global_accessible+marker_peak
global_accessible+marker_peak+ko_effect
...
```

## 3. 新输出

增强脚本：

```text
scripts/37_prepare_scperturb_atac_with_peak_features.py
```

输出：

```text
data/scperturb_atac/liscovitch_k562_gene_activity_chromvar_peaks.h5ad
data/scperturb_atac/liscovitch_k562_selected_peak_metadata.csv
data/scperturb_atac/liscovitch_k562_peak_regulatory_prior_scores_top1000.csv
```

`liscovitch_k562_selected_peak_metadata.csv` 现在包含：

- peak 坐标
- nearest gene
- peak type
- locus / marker / KO effect / accessibility / motif-TF 分数
- total score
- selection reason

## 4. KDM6A 结果

新运行目录：

```text
results/scperturb_atac_regulatory_peak_prior_kdm6a
```

结果：

| 指标 | 新 regulatory peak prior |
|---|---:|
| ROC-AUC | 0.674 |
| direction cosine | 0.771 |
| mean abs delta error | 0.061 |
| improved fraction | 0.375 |
| mean distribution improvement | -0.176 |

和上一版 peak model 相比：

- AUC 从约 0.658 提高到 0.674。
- direction cosine 从约 0.684 提高到 0.771。
- 说明强响应排序和整体方向更好。
- 但 distribution improvement 仍为负，说明 peak-level 单细胞分布形状仍然困难。

## 5. 解决单细胞分布形状：shape calibration

ATAC peak 的难点不是只有“平均方向”。很多 peak 的真实 KO 效应表现为：

- 可及细胞比例变化；
- 方差变化；
- 分布尾部变化；
- zero-inflated / sparse signal 的形状变化。

原来的虚拟 KO 是：

```text
virtual cell = control cell + predicted delta
```

这会保留 control 的分布形状，只移动平均值。因此当真实 KO 改变 peak 分布形状时，Wasserstein distance 仍然可能很差。

因此我们把 shape calibration 分成两档：

```text
--shape-calibrate variance
--shape-calibrate quantile
```

### 5.1 variance shape calibration

`variance` 只从训练 KO 学习每个 feature 的：

```text
KO std / control std
```

然后生成：

```text
virtual cell =
  control mean
  + shape_alpha * (control cell - control mean)
  + predicted KO delta
```

这仍然是 hard-constrained：

- KO 主方向仍由 residual/PLS + prior 决定；
- shape calibration 只轻微拉伸或压缩单细胞分布；
- 不使用 held-out KDM6A 的真实分布，因此不是结果泄漏。

### 5.2 zero-inflated / quantile shape calibration

`quantile` 是更适合 ATAC peak 的版本。它不只校准方差，还从训练 KO 学习：

- 每个 peak 的开放比例，也就是多少细胞中该 peak 是开放/非零的；
- 每个 feature 的分位数形状，例如 1%、5%、10%、50%、90%、99% 分位数；
- 稀疏非负 peak 的 closed/open 结构。

生成虚拟细胞时，它先保留 residual/PLS 预测出的 KO 平均方向，再把 control 细胞的单细胞排序映射到训练 KO 学到的分位数形状。对原始非负稀疏 peak/count，还会显式控制开放细胞比例。

一个重要保护规则是：

```text
如果 peak 已经被中心化/标准化成 peak state score，
软件只记录 open fraction，并使用 quantile shape calibration；
不会把 closed cells 强行置零。
```

原因是中心化后的 peak score 已经不是原始 0/1 或 count 矩阵。我们实测发现，对这种 score 强行 hard-zero 会显著破坏 Wasserstein distribution improvement。因此当前版本会区分“原始非负 peak 矩阵”和“中心化 peak state score”。

简化理解：

```text
1. residual/PLS 预测 KO 会让哪些状态特征上升或下降；
2. variance/quantile shape calibration 决定这些变化在单细胞之间怎么分布；
3. quantile 额外让 sparse peak 的“开/关比例”更像真实 KO 后的形状。
```

这一步仍然不使用 held-out KDM6A 的真实分布，而是只从训练 KO 学习通用形状，因此不是对测试结果的泄漏。

### Shape-calibrated KDM6A 结果

新运行目录：

```text
results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a
```

对比：

| 版本 | distribution improvement | improved fraction | direction | AUC |
|---|---:|---:|---:|---:|
| old peak prior | -0.159 | 0.378 | 0.684 | 0.658 |
| regulatory peak prior | -0.176 | 0.375 | 0.771 | 0.674 |
| regulatory peak prior + variance shape | -0.059 | 0.513 | 0.771 | 0.674 |
| regulatory peak prior + quantile shape | 0.166 | 0.788 | 0.771 | 0.674 |

解释：

- regulatory peak prior 主要改善方向和强响应排序。
- variance shape calibration 明显改善单细胞分布距离，但平均改善仍未转正。
- quantile shape calibration 进一步把 distribution improvement 提到正值，说明虚拟 KO 的 peak 单细胞分布整体比 control 更接近真实 KO。
- AUC、direction 和 MAE 基本不变是合理的，因为 shape calibration 主要改变单细胞分布形状，不改变 KO 平均方向和强响应排序。

## 6. 新 peak-level 图

默认 `run` 会输出：

```text
results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a/05_atac_peak_level_changes.png
```

新版图包含三部分：

1. 真实 KO vs 虚拟 KO 的 peak delta barplot。
2. true / virtual / error heatmap。
3. 真实 peak delta 与虚拟 peak delta 的方向一致性散点图。

用户版总览图：

```text
results/user_facing_figures/19_atac_peak_level_visualization.png
```

## 7. 结论

这次增强让 ATAC peak 选择更有调控解释性，而不是只看方差。结果说明：

```text
peak-level prior 能改善强响应排序和方向一致性；
quantile shape calibration 能进一步改善单细胞分布形状；
但 hard-zero open/closed 约束只适合原始非负 peak/count，不适合中心化 peak state score。
```

下一步应该继续加入：

- motif-to-peak annotation；
- promoter/enhancer 更精细分类；
- peak-gene linkage；
- batch-aware peak normalization。
