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

## 5. 解决单细胞分布形状：variance shape calibration

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

现在新增：

```text
--shape-calibrate variance
```

它只从训练 KO 学习每个 feature 的：

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

### Shape-calibrated KDM6A 结果

新运行目录：

```text
results/scperturb_atac_regulatory_peak_prior_shape_kdm6a
```

对比：

| 版本 | distribution improvement | improved fraction | direction | AUC |
|---|---:|---:|---:|---:|
| old peak prior | -0.159 | 0.378 | 0.684 | 0.658 |
| regulatory peak prior | -0.176 | 0.375 | 0.771 | 0.674 |
| regulatory peak prior + shape | -0.059 | 0.513 | 0.771 | 0.674 |

解释：

- regulatory peak prior 主要改善方向和强响应排序。
- variance shape calibration 明显改善单细胞分布距离。
- distribution improvement 仍未转正，说明 ATAC peak 形状还没有完全解决，但已经从“明显更差”变成“接近 control/true 的边界”。

## 6. 新 peak-level 图

默认 `run` 会输出：

```text
results/scperturb_atac_regulatory_peak_prior_kdm6a/05_atac_peak_level_changes.png
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
variance shape calibration 能改善单细胞分布形状；
但 ATAC peak 分布仍比 RNA/ADT 难预测。
```

下一步应该继续加入：

- motif-to-peak annotation；
- promoter/enhancer 更精细分类；
- peak-gene linkage；
- batch-aware peak normalization。
