# 当前 VKX 方法还缺什么，以及下一步怎么补齐

更新时间：2026-07-01

这份文档把当前模型方法的能力、缺口和下一步优先级整理成可执行清单。重点是让 VKX 不只是能跑出结果，而是能让用户看懂、能和已有方法比较、能明确知道什么情况下结果可信。

## 1. 已经具备的核心能力

VKX 当前已经支持：

- 单基因 KO 和双基因 KO。
- 小样本 perturbation 数据上的 prior-constrained residual/PLS baseline。
- 双 KO interaction residual，并已接入 reference model。
- batch covariate 显式建模，用 control-centered batch correction 降低批次混杂。
- RNA-only 输入自动转成 pathway/program score。
- RNA+ADT、RNA+ATAC、RNA+ADT+ATAC 等多模态输入接口。
- ATAC peak-level feature selection，支持 peak-gene linkage、motif-to-peak、marker peak 和 regulatory prior score。
- chromVAR/TF-target/motif/PPI/pathway 等系统先验。
- 普通 10X 无 KO 标签数据的 prediction-only reference application。
- hard-constrained uncertainty band，保守模拟 KO 主方向附近的不确定性范围。
- AUC/ROC、真实 vs 虚拟 heatmap、UMAP/PCA、KO summary card、figure package。
- 新增 `diagnose-results`，自动解释哪些 KO/feature 可信、哪些需要谨慎以及原因。

## 2. 还缺什么

### 2.1 缺正式横向 benchmark

当前 VKX 已有 `method-comparison` 接口，但还需要在同一批数据上系统比较：

- ridge/PLS/additive baseline
- scGen
- CPA
- GEARS
- CellOT
- diffusion/flow matching 类生成模型

比较时不能只看一个 AUC，而应该同时展示：

- AUC/ROC：能否识别强响应特征。
- MAE：变化幅度是否接近真实 KO。
- direction cosine / direction match：方向是否正确。
- R2：整体 delta 拟合程度。
- feature hit-rate：top 改变特征是否命中。
- prediction-only transfer confidence：普通 10X 应用时的迁移可靠性。

### 2.2 缺真正 full trimodal labelled benchmark

目前公开数据中，比较确定的是：

- RNA+ADT+perturbation labelled benchmark。
- RNA+ATAC+perturbation labelled benchmark。
- RNA+ADT+ATAC 但无 genetic perturbation label 的 DOGMA/TEA-seq 类数据。

真正公开、同一批细胞同时具备 RNA+ADT+ATAC+genetic perturbation label 的 benchmark 仍未确认。因此当前不能把 DOGMA/TEA-seq 当作准确性验证，只能做三模态输入兼容和 reference application。

### 2.3 神经生成模型仍是保守入口

当前 `vae`、`flow`、`diffusion` 入口仍然以 hard-constrained residual uncertainty 为主。这样做适合小样本，但不是完整自由生成模型。

下一步需要在有更多同类型 perturbation 数据后训练轻量 conditional VAE / flow matching / diffusion，并保持约束：

```text
virtual KO mean direction = residual/PLS or interaction residual hard constraint
neural generator = only local uncertainty around the constrained direction
```

### 2.4 MAPK/TGFB 等非线性 program 仍困难

MAPK/TGFB 类 program 往往涉及强非线性、反馈调控、细胞类型依赖和时间动态。当前线性 residual/PLS baseline 能提供稳定方向，但复杂幅度和组合效应可能不足。

下一步需要：

- 增强 pathway/TF/PPI/motif prior。
- 对 MAPK/TGFB 增加专项非线性校正项。
- 在 failure diagnosis 中标出这些 program，避免用户误读。

### 2.5 可视化还需要论文级统一排版

当前软件能自动输出很多图，但用于论文或给别人展示时，还需要统一成更清晰的 panel：

- Panel A：方法流程图。
- Panel B：输入模态矩阵和状态表示。
- Panel C：真实 vs 虚拟 heatmap。
- Panel D：AUC/ROC 曲线。
- Panel E：UMAP/PCA 单细胞状态移动。
- Panel F：KO summary card。
- Panel G：failure diagnosis overview。
- Panel H：ATAC peak-level locus/peak/motif 图。

## 3. 新增 failure diagnosis 的意义

用户之前的问题是：图很多，但不知道方法到底好不好。新增诊断模块就是为了解决这个问题。

命令：

```bash
python -m vkx.cli diagnose-results \
  --delta-csv results/labelled_virtual_ko/delta_table.csv \
  --manifest-csv results/labelled_virtual_ko/derived_state_manifest.csv \
  --out-dir results/labelled_virtual_ko/diagnosis
```

输出：

- `ko_failure_diagnosis.csv`
- `feature_failure_diagnosis.csv`
- `failure_diagnosis_report.md`
- `01_failure_diagnosis_overview.png`
- `02_feature_error_heatmap.png`

它会自动标记：

- `ok`: 方向和幅度基本可信。
- `warning`: 误差偏大、真实变化接近 0 或可能 overcall。
- `high_risk`: 方向错误或漏掉真实强效应。
- `prediction_only`: 无真实 KO 标签，不能当准确性评估。

常见原因标签：

- `direction_mismatch`
- `large_error`
- `missed_real_effect`
- `overcalled_effect`
- `near_zero_real_effect`
- `sparse_atac_shape_or_prior_issue`
- `difficult_mapk_tgfb_program`
- `prediction_only_no_truth`

## 3.1 新增一键用户可读结果报告

为了避免用户在多个输出目录之间来回找图，新增：

```bash
python -m vkx.cli summarize-result \
  --result-dir results/labelled_virtual_ko
```

它会自动检查结果目录中是否存在：

- `delta_table.csv`
- `predicted_ko_delta.csv`
- `derived_state_manifest.csv`
- `auc_summary.csv`
- `transfer_confidence.csv`
- PNG 图

然后统一生成：

- `readable_result_report/user_readable_result_report.md`
- `readable_result_report/ko_cards/`
- `readable_result_report/diagnosis/`
- `readable_result_report/figure_package/`

这个总报告会首先说明当前结果属于哪种模式：

- labelled evaluation：可以解释真实准确率、AUC、MAE/R2 和真实 vs 虚拟 heatmap。
- prediction-only：只能解释预测状态变化和迁移置信度，不能当作真实准确率。
- double-KO evaluation：可以解释 additive vs interaction residual。

因此，后续用户只需要打开 `user_readable_result_report.md`，就能看到：

1. 这个结果能不能报告准确性。
2. 每个 KO 的一句话结论。
3. 哪些 feature 高风险。
4. KO summary card。
5. failure diagnosis 图。
6. 完整 figure package 入口。

## 4. 下一步推荐优先级

1. 先做正式横向 benchmark：VKX vs ridge/PLS/additive/scGen/CPA/GEARS/CellOT。
2. 把 benchmark 结果统一进 figure package，形成一套论文级主图。
3. 继续搜索并人工确认 full RNA+ADT+ATAC+genetic perturbation labelled benchmark；没有确认前保持 `not_confirmed_public`。
4. 针对 MAPK/TGFB 和 ATAC peak sparse shape 做专项增强。
5. 在更多 perturbation 数据上训练轻量 neural generator，但保持 hard constraint。

## 5. 对外表述建议

VKX 最稳妥的定位是：

> VKX is a small-sample, multimodal, prior-constrained virtual knockout framework that predicts interpretable cell-state shifts and provides explicit reliability diagnostics. It is designed for settings where labelled perturbation data and compute are limited, and where users need transparent pathway/protein/ATAC-level outputs rather than a black-box free generator.

中文表述：

> VKX 是一种面向小样本多模态单细胞数据的先验约束虚拟敲除方法。它不依赖大规模预训练和高算力，而是利用真实 perturbation 数据、系统生物学网络先验和多组学状态表示，预测单基因或双基因敲除后的细胞状态变化，并输出可解释的图、表和可靠性诊断。
