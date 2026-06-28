# 增强 KO 方向先验：constrained direction prior

这一步不再加深生成模型，而是先把“KO 后平均状态往哪个方向移动”预测得更稳。

## 做了什么

对每个 held-out KO，模型融合了五类方向预测：

- `pls`：只用 Papalexi 训练 KO 的系统先验 PLS。
- `cross_pls`：先用 Norman 数据预训练一个 perturbation response embedding，再用于 Papalexi。
- `knn3`：找系统先验最相似的 3 个训练 KO，做加权平均。
- `knn7`：找系统先验最相似的 7 个训练 KO，做加权平均。
- `zero`：收缩到 0，防止方向预测过猛。

融合权重不是手调的，而是在 Papalexi 训练 KO 上做 leave-one-KO-out 自动选择。

## 自动选择的权重

- pls: 0.20
- knn3: 0.40
- zero: 0.40

## 训练 KO 留一验证

- pls: LOO MAE 0.064, ensemble weight 0.20
- cross_pls: LOO MAE 0.064, ensemble weight 0.00
- knn3: LOO MAE 0.063, ensemble weight 0.40
- knn7: LOO MAE 0.061, ensemble weight 0.00
- zero: LOO MAE 0.061, ensemble weight 0.40
- ensemble: LOO MAE 0.057, ensemble weight 1.00

## Held-out KO 测试结果

- cross_pls: 平均分布改进 0.138, 改进比例 69.4%
- pls: 平均分布改进 0.137, 改进比例 69.4%
- knn3: 平均分布改进 0.088, 改进比例 66.7%
- constrained_ensemble: 平均分布改进 0.088, 改进比例 72.2%
- knn7: 平均分布改进 0.070, 改进比例 66.7%
- zero: 平均分布改进 0.000, 改进比例 0.0%

## 图

- `results/figures/papalexi_constrained_direction_prior_summary.png`
- `results/figures/papalexi_constrained_direction_prior_by_ko_heatmap.png`
- `results/figures/papalexi_constrained_direction_prior_distributions.png`

## 当前结论

这个实验检验了三件事：更强 KO 方向先验、更多训练 KO 的留一调权、以及 Norman -> Papalexi 的跨数据预训练 embedding。

如果 constrained ensemble 超过普通 PLS，说明方向先验增强有效；如果没有超过，说明当前限制主要来自跨数据状态不一致或训练 KO 数量仍然不足。无论哪种结果，都比继续盲目加深 VAE/flow/diffusion 更有信息量。
