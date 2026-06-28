# Cell-level 条件生成模型第一版

这一版开始从 KO 组平均预测，升级到单细胞条件生成。

## 模型输入

```text
单个 control 细胞的 pathway/protein state
+ KO gene 的 Reactome/MSigDB/TF-target/PPI 先验特征
```

## 模型输出

```text
该细胞在 KO 条件下的 pathway/protein state
```

## 训练方式

Papalexi 数据中没有同一个细胞 KO 前后的真实配对，因此训练时使用随机配对：

```text
control cell state + KO condition -> sampled true KO cell state
```

这使模型学习条件分布，而不是单个细胞的精确一一对应。

## 测试方式

完全留出这些 KO，不参与训练：

```text
STAT1, JAK2, IFNGR2, IRF1
```

然后用 control 细胞生成这些 KO 的虚拟单细胞状态，并与真实 KO 细胞分布比较。

## 输出图

- `results/figures/papalexi_cell_level_generator_holdout_distributions.png`
- `results/figures/papalexi_cell_level_generator_state_space.png`
- `results/figures/papalexi_cell_level_generator_metric_summary.png`

## 如何解释

如果 virtual KO cells 的分布比 control cells 更接近 true KO cells，说明模型学到了 KO 条件下的单细胞状态移动。

当前模型是第一版轻量 MLP 条件生成器，还不是最终的 VAE / flow matching / diffusion 模型。
