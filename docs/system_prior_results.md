# 系统先验优化结果

## 使用的先验库

通过 Enrichr API 下载并本地保存为 GMT：

- `Reactome_Pathways_2024`
- `MSigDB_Hallmark_2020`
- `ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X`
- `PPI_Hub_Proteins`

本地路径：

- `data/priors/reactome.gmt`
- `data/priors/hallmark.gmt`
- `data/priors/tf_target.gmt`
- `data/priors/ppi_hub.gmt`

## 方法

对每个 KO 或组合 KO，将基因映射到系统先验 term：

```text
KO genes -> Reactome / Hallmark / TF-target / PPI hub features
```

然后在 Norman 2019 中做：

```text
single-gene perturbation training -> double-gene perturbation testing
```

## Norman 2019 结果

系统先验模型选中了 449 个与扰动基因相关的 prior terms。

在全体 52 个双基因组合上：

- Granulocyte/apoptosis: additive R2 ≈ 0.55，system prior R2 ≈ 0.75。
- MAPK/TGFB: additive R2 ≈ -0.64，system prior R2 ≈ -0.14。
- Pro-growth: additive R2 ≈ -0.89，system prior R2 ≈ -0.08。
- Erythroid: additive R2 ≈ 0.51，system prior R2 ≈ 0.39。
- Pioneer TF: additive R2 ≈ 0.50，system prior R2 ≈ 0.45。

在含未见基因的 40 个组合上：

- Granulocyte/apoptosis: additive R2 ≈ 0.20，system prior R2 ≈ 0.53。
- MAPK/TGFB: additive R2 ≈ -0.49，system prior R2 ≈ -0.09。
- Pro-growth: additive R2 ≈ -1.03，system prior R2 ≈ -0.18。

结论：系统先验明显改善了部分未见基因组合的外推，但对原本 additive 已经很强的程序不一定继续提升。

## Papalexi 多模态 AUC

强响应识别任务中：

- IFNG-JAK-STAT 下降：PLS ROC-AUC ≈ 0.95。
- PDL1 蛋白下降：Ridge ROC-AUC = 1.00。
- CD86 蛋白强响应：PLS ROC-AUC ≈ 0.98。

这说明多模态读数更适合用于“识别强响应 KO”的排序任务，而不仅仅是连续值回归。

## 图文件

- `results/figures/papalexi_auc_strong_response.png`
- `results/figures/papalexi_pathway_protein_correlation.png`
- `results/figures/norman_system_prior_r2.png`
- `results/figures/norman_system_prior_scatter.png`
- `results/figures/norman_prior_term_hits.png`

## 解释

系统先验的价值主要体现在：

1. 未见基因可以通过共享 pathway/TF/PPI term 获得特征。
2. 双基因组合不再只能做简单加和。
3. 可解释性更强，可以追踪某个组合 KO 命中了哪些 pathway、TF 或 PPI hub。

当前限制：

1. Norman 轻量版只有 2000 个基因，程序分数仍受覆盖限制。
2. 先验 term 仍是二值/重叠特征，没有利用边方向、激活/抑制关系和网络距离。
3. 训练仍是 perturbation-level 聚合，下一步应升级到 cell-state-conditioned model。
