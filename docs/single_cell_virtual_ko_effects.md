# 单细胞层面的虚拟敲除效果

之前的多数结果是 perturbation-level，也就是把同一个 KO 下的细胞聚合后看平均变化。为了展示单细胞效果，现在补充了单细胞分布图。

## 图怎么看

文件：

- `results/figures/papalexi_single_cell_virtual_ko_distributions.png`
- `results/figures/papalexi_single_cell_state_space_pca.png`
- `results/figures/norman_single_cell_virtual_combo_distributions.png`

每个小图比较三类细胞状态：

- `control`：真实 control 细胞的单细胞分布。
- `true KO cells`：真实 KO 细胞的单细胞分布。
- `virtual KO cells`：把模型预测的 KO delta 加到 control 细胞上，得到的虚拟 KO 单细胞分布。

如果 `virtual KO cells` 和 `true KO cells` 的分布接近，说明模型在单细胞状态空间里预测得较好。

## 当前结论

Papalexi 中，STAT1/JAK2/IFNGR2 相关 KO 对 IFNG-JAK-STAT pathway 和 PDL1 protein 的单细胞分布有明显移动，虚拟 KO 能捕捉主要方向，但幅度有时偏低。

Norman 中，部分双基因组合如 `AHR+KLF1`、`CEBPB+CEBPA` 的主要程序变化方向能被捕捉；但 `MAPK1+TGFBR2` 等组合仍预测不足。

## 重要限制

当前 virtual KO cells 是一个简化版本：

```text
virtual single-cell score = control single-cell score + predicted perturbation delta
```

也就是说，它展示的是“模型预测的 KO 后单细胞状态分布”，但还不是最终的 cell-level 条件生成模型。下一步需要训练真正的 cell-state-conditioned model，让每个 control 细胞都有自己的条件化 KO 响应。
