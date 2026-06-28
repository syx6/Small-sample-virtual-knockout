# 用户视角的方法效果与输入输出

## 一句话结论

当前原型**不是已经能精确预测所有 KO 后表达变化的最终模型**。它目前最可靠的能力是：

1. 识别哪些 KO 会产生强通路/蛋白响应。
2. 用 pathway + protein 输出解释 KO 后细胞状态。
3. 在部分双基因组合上做可解释外推，系统先验对未见基因组合有帮助。

还不够强的地方是：

1. 连续变化幅度的精确预测仍一般。
2. MAPK/TGFB、Pro-growth 这类非线性组合效应仍难。
3. 当前模型还是 perturbation-level 聚合模型，不是最终 cell-level 条件生成模型。

## 当前效果判断

| 问题 | 指标 | 判断 | 含义 |
| --- | --- | --- | --- |
| 能不能识别强响应 KO? | IFNG decrease AUC=0.95; PDL1 decrease AUC=1.00 | 好 | 适合先回答哪些 KO 可能产生明显通路/蛋白响应。 |
| 能不能精确预测连续变化幅度? | PDL1 R2=0.68; IFNG pathway R2=0.23 | 中等 | 蛋白表型较稳，RNA 通路幅度预测还需要更好的 cell-level 模型。 |
| 能不能外推到双基因组合 KO? | unseen-combo mean R2 gain=0.52 | 部分有效 | 系统先验能改善未见基因组合，但组合效应还不能完全解决。 |
| 小样本多模态思路是否成立? | pathway-protein coupling r≈0.70 | 成立 | 多模态表型给通路状态提供了可解释约束，是当前方法最强证据。 |

## 正确看图方式

AUC 的正式展示应使用 ROC 曲线，而不是只看柱状图。柱状图只能作为 AUC 数值摘要。

推荐正式使用：

- `results/figures/papalexi_roc_curves.png`：强响应识别的 ROC 曲线。
- `results/figures/norman_true_vs_pred_system_prior.png`：双基因组合 KO 的真实值 vs 预测值。
- `results/figures/norman_r2_improvement_system_prior.png`：系统先验相对 additive baseline 的 R2 改善。

R2 可以用柱状图展示不同模型/任务的数值，但它不直观。更清楚的方式是看真实值-预测值散点：点越靠近对角线，连续效应预测越好；点偏离越大，说明虽然可能能排序，但幅度预测不准。

## 给别人使用时的输入

最小输入：

- 单细胞 RNA 表达矩阵，格式可以是 `.h5ad` / `.h5mu` / count matrix。
- 每个细胞的 perturbation 标签，例如 `STAT1`、`JAK2`、`KLF1+BAK1`、`ctrl`。
- control/negative-control 细胞标记。

推荐输入：

- RNA + protein/CITE-seq。
- RNA + ATAC 或 gene activity。
- cell type / sample / batch metadata。
- KO gene list，支持单基因和多基因组合。

系统先验输入：

- Reactome / MSigDB pathway gene sets。
- TF-target gene sets。
- PPI hub 或 gene network。

## 给别人使用时的输出

主要输出：

- 每个 KO 或组合 KO 的 pathway activity delta。
- 可选 protein phenotype delta，例如 PDL1、CD86。
- 强响应 KO 排名。
- 多基因组合 KO 的预测表。
- 命中的 pathway / TF / PPI prior terms，用于解释。

评估输出：

- 回归：MAE、R2。
- 排序/筛选：ROC-AUC、PR-AUC。
- 可视化：效果判断图、输入输出图、强响应图、组合 KO 外推图。

## 推荐对外表述

这个方法目前应定位为：

> 一个面向小样本多模态单细胞数据的机制先验驱动虚拟敲除框架，用于预测和解释 KO 后的通路/蛋白状态变化，特别适合强响应筛选和多基因组合扰动的初步外推。

不建议现在声称：

- 可以精确重建 KO 后全转录组。
- 可以可靠预测所有未知多基因组合。
- 已经优于大型预训练生成模型。
