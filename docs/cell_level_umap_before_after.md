# Cell-level UMAP：敲除前后单细胞状态移动

这组图专门回答一个直观问题：

```text
敲除前的 control cells 在哪里？
虚拟敲除后的 cells 移动到哪里？
真实实验 KO cells 在哪里？
```

## 图怎么看

每个 KO 都有两列：

- 左列：敲除前，只显示 control cells。
- 右列：敲除后，显示 virtual KO cells 和 true KO cells，同时淡灰色保留 control cells 作为参照。

箭头含义：

- 橙色实线箭头：control 质心 -> virtual KO 质心。
- 绿色虚线箭头：control 质心 -> true KO 质心。

如果橙色箭头和绿色箭头方向接近，并且橙色点云靠近绿色点云，说明虚拟敲除捕捉到了单细胞状态变化。

## 单基因敲除

图：`results/figures/papalexi_cell_level_umap_single_gene_before_after.png`

使用 Papalexi 多模态数据：pathway score + protein score。

平均 UMAP 质心改进：`0.163`

## 多基因组合敲除

图：`results/figures/norman_cell_level_umap_multi_gene_before_after.png`

使用 Norman 组合扰动数据：gene program score。

平均 UMAP 质心改进：`0.332`

## 辅助解释表

注意：UMAP 质心距离只用于解释可视化，不作为正式性能指标。正式评价仍然看前面的 Wasserstein、ROC-AUC、R2。

| dataset | KO | control->true | virtual->true | improvement |
|---|---:|---:|---:|---:|
| Papalexi single-gene KO | IFNGR2 | 1.366 | 1.907 | -0.396 |
| Papalexi single-gene KO | IRF1 | 1.484 | 0.481 | 0.676 |
| Papalexi single-gene KO | JAK2 | 0.761 | 0.703 | 0.076 |
| Papalexi single-gene KO | STAT1 | 1.973 | 1.387 | 0.297 |
| Norman multi-gene KO | AHR+KLF1 | 1.391 | 1.001 | 0.28 |
| Norman multi-gene KO | CBL+UBASH3B | 3.127 | 2.957 | 0.054 |
| Norman multi-gene KO | CEBPB+CEBPA | 4.689 | 0.462 | 0.902 |
| Norman multi-gene KO | MAPK1+TGFBR2 | 2.512 | 2.279 | 0.093 |

## 当前结论

这两张 UMAP 图能让人直接看到：模型不是只输出一个分数，而是在 cell-level 状态空间里把 control cells 移向 KO-like 状态。

但也要诚实：如果某些 KO 的橙色云团没有靠近绿色云团，说明当前稳定基线还不够，需要在这个基础上继续接 VAE / flow matching / diffusion，而不是直接从零训练复杂模型。
