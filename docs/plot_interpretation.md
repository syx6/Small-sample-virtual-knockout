# 图该怎么看

## 1. ROC 曲线

文件：

- `results/figures/papalexi_roc_curves.png`

用途：

- 判断模型能否识别“强响应 KO”。
- 曲线越靠近左上角越好。
- AUC 越接近 1 越好，0.5 接近随机。

当前结论：

- IFNG pathway decrease：AUC 约 0.95。
- PDL1 protein decrease：AUC 约 0.96-1.00。
- CD86 protein strong change：AUC 约 0.95-0.98。

这说明当前模型在“筛出强响应 KO”这个任务上表现好。

## 2. 真实值 vs 预测值散点图

文件：

- `results/figures/norman_true_vs_pred_system_prior.png`

用途：

- 判断连续效应幅度预测是否准确。
- 点越靠近虚线 y=x，预测越准。
- 如果点只大致有趋势但偏离对角线，说明模型能抓方向或排序，但幅度还不准。

当前结论：

- Granulocyte/apoptosis 预测较好。
- Erythroid 有一定趋势。
- MAPK/TGFB 和 Pro-growth 仍然较弱。

## 3. R2 improvement 图

文件：

- `results/figures/norman_r2_improvement_system_prior.png`

用途：

- 只回答一个问题：系统先验是否比 additive baseline 更好。
- 正值表示系统先验提升。
- 负值表示系统先验反而不如简单加和。

当前结论：

- 系统先验明显改善 Pro-growth、MAPK/TGFB、Granulocyte/apoptosis。
- 对 Erythroid 和 Pioneer TF 不一定提升，说明原始 additive baseline 已经较强或系统先验引入了噪声。

## 总结

最标准、最不容易误解的展示组合是：

1. ROC 曲线：说明强响应筛选能力。
2. 真实值 vs 预测值散点：说明连续效应预测能力。
3. R2 improvement：说明系统先验是否真的有帮助。
4. 具体 KO 组合热图/条形图：说明敲了什么基因，真实引起什么程序变化，模型预测了什么。

## 4. 具体 KO 组合解释图

文件：

- `results/figures/norman_combo_true_vs_pred_heatmap.png`
- `results/figures/norman_representative_combo_bars.png`
- `results/norman_combo_explanation_table.csv`

用途：

- 每一行是一个双基因 KO，例如 `AHR+KLF1`。
- 左侧/绿色是真实双基因扰动相对 control 的程序分数变化。
- 右侧/橙色是模型预测的程序分数变化。
- 正值表示该程序增强，负值表示该程序降低。

这类图比单纯 R2 更适合回答生物问题：

```text
敲了什么基因？
真实改变了什么通路/程序？
模型预测对了吗？
错在哪里？
```
