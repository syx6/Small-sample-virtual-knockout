# True KO vs virtual KO heatmap visualization

这组图专门展示：

```text
真实敲除后状态变化
vs
虚拟敲除预测状态变化
vs
预测误差
```

每个数据集都有一张三栏 heatmap：

- Real KO change：真实 KO 相对 control 的变化。
- Virtual KO change：虚拟 KO 相对 control 的变化。
- Prediction error：virtual - true。

另外有一张总览图：

- `results/figures/multi_dataset_true_vs_virtual_agreement_heatmap.png`

单数据集 heatmap：

- `results/figures/multi_dataset_heatmap_papalexi_true_virtual_error.png`
- `results/figures/multi_dataset_heatmap_norman_true_virtual_error.png`
- `results/figures/multi_dataset_heatmap_datlinger_true_virtual_error.png`
- `results/figures/multi_dataset_heatmap_dixit_true_virtual_error.png`

读图规则：

- 左右两张热图颜色方向一致，说明模型预测到了正确变化方向。
- 误差热图颜色越浅，说明虚拟 KO 越接近真实 KO。
- Direction agreement 越接近 1，说明整体变化方向越一致。
