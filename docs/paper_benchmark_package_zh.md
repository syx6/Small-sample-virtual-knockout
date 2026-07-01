# VKX 论文级 benchmark 图包说明

这个图包对应命令：

```bash
python -m vkx.cli paper-benchmark \
  --formal-dir results/formal_benchmark_papalexi_adaptive_boosted \
  --result-dirs results/software_interface_single_gene_demo,results/software_interface_double_gene_demo,results/hmpcite_multimodal_doubleko_cebp_med12,results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a \
  --out-dir results/paper_benchmark_package
```

## 为什么要新增这个图包

之前的结果分散在不同目录里，用户很难直接判断：

- VKX 是否真的比自己早期版本更好；
- 是否已经和 scGen、CPA、GEARS、CellOT 做了公平比较；
- 单敲、双敲、多模态、ATAC peak 是否都有结果；
- 哪些结论可靠，哪些还不能宣传。

`paper-benchmark` 的目的就是把这些问题收束成一张论文级主图、四张可单独查看的分图和一份报告。

## 输出图

0. `00_publication_main_figure.png`

   推荐作为论文主图的 A-D 多面板 figure：

   - A：正式方法横向 benchmark；
   - B：AUC 的 ROC 曲线形式；
   - C：真实 KO 与 VKX-Boosted 虚拟 KO 的 heatmap；
   - D：单敲、双敲、多模态和 ATAC peak 覆盖情况。

1. `01_method_leaderboard.png`

   正式方法比较图。它把 VKX/ResponseBoosted、PLS、Ridge、ConstrainedEnsemble、Additive 和 scGen/CPA/GEARS/CellOT 槽位放在同一张图里。

   重点：scGen/CPA/GEARS/CellOT 如果没有提供同数据预测文件，会显示为 `not_run`，不会被偷偷算进比较。这能避免过度宣传。

2. `02_auc_roc_curves.png`

   真正的 ROC 曲线图，而不是只有 AUC 柱状图。它回答：模型能不能把真实 KO 中变化强的 pathway/protein feature 排到前面？

3. `03_real_vs_virtual_method_heatmap.png`

   真实 KO 和各方法虚拟 KO 的 delta heatmap。它回答：到底敲了什么基因，哪些通路/蛋白变化了，预测和真实差在哪里。

4. `04_single_double_multimodal_gallery.png`

   单敲、双敲、多模态和 ATAC peak 结果总览图。它回答：VKX 是否只是在一个 RNA-only 示例上有效，还是能覆盖更广的输入场景。

## 当前结论

当前标准 Papalexi 小型 state-score benchmark 中，`PLS` 的 AUC 最高，而 `ResponseBoosted` 是表现最好的 VKX 变体之一：

- PLS AUC 约 0.736；
- VKX-Boosted AUC 约 0.625；
- VKX-Boosted MAE 约 0.139，是当前 scored methods 中较低的误差。

这说明加入 response-strength prior 后，VKX 相比原始 VKX 有进步，但也说明传统 PLS baseline 在这个小型数据上仍然很强，不能回避。

但是，这个结果仍然不能说明 VKX 已经强过 scGen、CPA、GEARS、CellOT。原因是这些外部深度方法还没有在同一数据、同一 holdout KO、同一指标下提供预测文件。

## 三件事完成情况

### 1. 正式横向 benchmark

已完成统一接口和图包：

- 内部方法：VKX、ResponseBoosted、CalibratedEnsemble、ConstrainedEnsemble、PLS、Ridge、Additive；
- 外部方法槽位：scGen、CPA、GEARS、CellOT；
- 指标：AUC、R2、MAE、direction cosine、feature hit-rate；
- 图：leaderboard、ROC curve、真实 vs 虚拟 heatmap。

尚未完成的是：真正运行 scGen/CPA/GEARS/CellOT 并导入它们的同数据预测结果。

### 2. 提升 VKX 短板

当前已纳入图包和报告的优化方向：

- `ResponseBoosted`：增强 IFN/JAK/STAT、MAPK/TGFB、MYC/E2F/cell-cycle 等响应家族的幅度；
- 双敲 interaction residual：用于处理非简单相加的组合 KO 效应；
- ATAC quantile / zero-inflated shape calibration：用于改善 peak-level 单细胞分布形状；
- 多模态状态表示：RNA pathway/program score + ADT/protein + ATAC/chromVAR/peak score。

仍然困难：

- MAPK/TGFB 这类上下文依赖 program；
- 双敲强非线性；
- raw peak count 的开放比例和分位数形状；
- full RNA+ADT+ATAC+perturbation labelled benchmark。

### 3. 论文级图包

已新增：

- 统一主图输出；
- 中文解释报告；
- source figure 自动复制到 `assets/`；
- 表格输出：`paper_method_metrics.csv`、`paper_method_availability.csv`、`paper_result_summary.csv`、`paper_figure_index.csv`。

## 推荐的论文表述

目前建议这样写：

> VKX is a small-sample, interpretable and multimodal-prior-constrained virtual knockout framework. In current labelled perturbation benchmarks, the response-boosted VKX variant improves over internal classical baselines, while strict comparison with scGen, CPA, GEARS and CellOT requires same-dataset external predictions.

中文意思是：

> VKX 目前更适合被定位为“小样本、可解释、多模态先验约束”的虚拟敲除框架。当前结果证明它比我们自己的早期 baseline 更好，但还不能宣称全面超过 scGen、CPA、GEARS 和 CellOT。
