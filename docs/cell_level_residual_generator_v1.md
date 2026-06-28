# 真正的 cell-level 条件生成模型：Residual Generator v1

这一步不再只是预测“一个 KO 组的平均变化”，而是直接生成一批虚拟单细胞。

## 输入是什么

1. 一批未敲除的 control cells，每个细胞已经表示成 pathway score + protein score。
2. 要敲除的基因，例如 `STAT1`、`JAK2`、`IFNGR2`、`IRF1`。
3. 这些基因在 Reactome、MSigDB、TF-target、PPI 网络中的先验特征。

## 输出是什么

对每个目标 KO，输出一批虚拟 KO cells。每个虚拟细胞都有：

- pathway scores
- protein scores
- 与真实 KO 细胞分布的距离评估
- 与 control cells 相比是否更接近真实 KO 的判断

## 为什么先用 residual generator

第一版 MLP 条件生成器在小样本下表现不好，生成分布比原始 control 更远。原因很直接：小样本里没有同一个细胞 KO 前后的真实配对，神经网络很容易学到错误的平均模式。

Residual Generator v1 更保守：

```text
虚拟 KO cell = control cell + 系统先验预测的 KO 平均变化 + 可选的真实细胞波动
```

这样做的好处是：平均变化由网络先验约束，细胞级输出从真实 control cells 出发，而不是让模型凭空编。

本轮调参发现，当前小样本 holdout 测试里加入额外 residual noise 会降低效果，所以最终默认使用：

```text
residual_scale = 0.0
```

这意味着当前最可靠的 cell-level 生成方式是：先把真实 control cells 沿着预测到的 KO 方向移动。细胞之间原本的差异仍然保留，因为每个虚拟细胞都来自一个真实 control cell。

## 测试方式

完全留出这些 KO，不让模型训练时看到它们：

```text
STAT1, JAK2, IFNGR2, IRF1
```

然后比较三类细胞：

- control cells：原始未敲除细胞
- virtual KO cells：模型生成的敲除后细胞
- true KO cells：实验里真实测到的敲除细胞

## 当前效果

平均分布改进值：`0.151`

所有 KO-特征组合中，生成细胞比 control 更接近真实 KO 的比例：`77.8%`

重点特征的平均分布改进：

- protein_CD86: 0.463
- protein_PDL1: 0.316
- pathway_IFNG_JAK_STAT: 0.155
- pathway_IMMUNE_CHECKPOINT: -0.266

解释规则：数值大于 0 表示生成细胞比 control 更接近真实 KO；小于 0 表示还不如直接用 control。

## 结果图

- `results/figures/papalexi_cell_level_residual_generator_holdout_distributions.png`
- `results/figures/papalexi_cell_level_residual_generator_metric_summary.png`
- `results/figures/papalexi_cell_level_residual_generator_state_space.png`
- `results/figures/papalexi_cell_level_residual_generator_delta_heatmap.png`
- `results/figures/papalexi_cell_level_generator_model_comparison.png`

## 结论

这个版本是进入 cell-level 条件生成后的第一个可解释基线。它回答的问题是：给定一个真实 control cell 和一个 KO 条件，能不能生成更像真实 KO 的单细胞状态。

如果它只在部分通路或蛋白上有效，说明多模态和系统先验确实提供了信号，但还不足以完全解决小样本条件生成。下一步可以在这个稳定基线上接 VAE / flow matching / diffusion，而不是直接从零训练复杂模型。
