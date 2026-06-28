# Cell-level 条件生成模型比较：Residual / VAE / Flow / Diffusion

这一步把三类生成模型放到同一个小样本多模态虚拟敲除框架里比较，并额外测试了“接在稳定 residual baseline 上”的 guided 版本。

## 比较对象

- Residual baseline：control cell + 系统先验预测的 KO 平均移动。
- Conditional VAE：直接学习 KO 条件下的潜在分布。
- Flow matching：直接学习从 control cell 到 KO cell 的移动速度场。
- Diffusion：直接从噪声去噪生成 KO-like cells。
- Guided VAE / Guided Flow / Guided Diffusion：先用 residual baseline 得到 KO anchor，再让深度模型学习剩余细胞级修正。

## 公平测试方式

完全留出这些 KO，不参与训练：

```text
STAT1, JAK2, IFNGR2, IRF1
```

所有模型都使用同一组 pathway/protein state 和同一套 Reactome/MSigDB/TF-target/PPI 条件先验。

## 当前结果

- Residual baseline: 平均分布改进 0.174, 改进比例 77.8%
- Conditional VAE: 平均分布改进 -0.189, 改进比例 33.3%
- Diffusion: 平均分布改进 -0.209, 改进比例 50.0%
- Guided Conditional VAE: 平均分布改进 -0.566, 改进比例 41.7%
- Flow matching: 平均分布改进 -0.577, 改进比例 22.2%
- Guided Flow matching: 平均分布改进 -0.766, 改进比例 38.9%
- Guided Diffusion: 平均分布改进 -0.769, 改进比例 36.1%

解释：平均分布改进大于 0 表示生成细胞比 control 更接近真实 KO；改进比例表示多少个 KO-特征组合是有帮助的。

## 图

- `results/figures/papalexi_cell_level_generator_deep_model_comparison.png`
- `results/figures/papalexi_cell_level_generator_model_by_ko_heatmap.png`
- `results/figures/papalexi_cell_level_generator_deep_model_distributions.png`

## 当前结论

当前小样本条件下，复杂生成模型没有超过稳定 residual baseline。直接版 VAE、flow matching、diffusion 多数为负；guided 版本也没有改善，说明现在的细胞级残差还没有足够稳定的规律可供小模型学习。

这个结果支持一个重要开发判断：下一步不应该盲目加深模型，而应该先增强约束，包括更强的 KO 方向先验、更多训练 KO、跨数据预训练，或者把 residual baseline 作为 hard constraint，只允许生成模型学习低幅度不确定性。
