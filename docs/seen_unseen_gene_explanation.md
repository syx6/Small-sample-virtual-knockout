# seen gene / unseen gene 是什么意思

在 Norman 2019 的多基因组合实验里，我们故意做了一个外推测试：

```text
训练集：单基因扰动
测试集：双基因组合扰动
```

因此：

- **seen gene**：这个基因在训练集中出现过单基因扰动。例如训练里有 `KLF1` 单基因扰动，那么在 `AHR+KLF1` 里，`KLF1` 是 seen gene。
- **unseen gene**：这个基因没有单基因训练样本，只在双基因组合里出现。例如训练里没有 `AHR` 单基因扰动，那么 `AHR+KLF1` 这个组合属于含 unseen gene 的组合。
- **seen combo**：双基因组合里的两个基因都在单基因训练集中出现过。
- **unseen combo**：双基因组合里至少一个基因没有单基因训练样本。

这不是说生物学上“没见过这个基因”，而是说模型训练时没有见过这个基因的单独扰动效果。

## 为什么要这样分

因为多基因虚拟敲除最难的不是“见过的基因怎么相加”，而是：

1. 没有某个基因的单基因 KO 数据时，能不能靠 pathway/TF/PPI 先验外推？
2. 两个基因组合时，效应是否仍然接近加和？
3. 哪些组合出现非线性，模型预测会失败？

## 具体结果怎么看

推荐看两个文件：

- `results/norman_combo_explanation_table.csv`
- `results/figures/norman_combo_true_vs_pred_heatmap.png`
- `results/figures/norman_representative_combo_bars.png`

这些文件会告诉你：

```text
敲了什么基因组合
真实引起了哪个通路/程序变化
模型预测了什么变化
预测和真实差多少
这个组合命中了哪些系统先验 term
```
