# 当前 VKX 方法状态、正式 benchmark 与下一步

更新时间：2026-07-01

这份文档整理当前 VKX 的最新软件能力、已经补齐的部分、仍然不能夸大的边界，以及下一步最重要的实验。

## 1. 当前已经具备的能力

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
- hard-constrained uncertainty band。
- AUC/ROC、真实 vs 虚拟 heatmap、UMAP/PCA、KO summary card、figure package。
- `diagnose-results`：解释哪些 KO/feature 可信、哪些需要谨慎以及原因。
- `summarize-result`：把结果目录整理成用户可读总报告。
- `formal-benchmark`：正式横向比较 VKX、PLS、ridge、additive，并保留 scGen/CPA/GEARS/CellOT 外部预测槽位。
- `train-hard-generator`：在稳定 VKX baseline 上训练 hard-constrained residual generator。

## 2. 正式横向 benchmark

新增命令：

```bash
python -m vkx.cli formal-benchmark \
  --state-csv results/labelled_virtual_ko/derived_state_scores.csv \
  --ko-col ko_target \
  --target-kos STAT1,JAK2,STAT1+JAK2 \
  --prior-dir data/priors \
  --methods vkx,pls,ridge,additive,scgen,cpa,gears,cellot \
  --out-dir results/formal_method_benchmark
```

固定输出：

- `formal_benchmark_predictions.csv`
- `formal_benchmark_truth.csv`
- `formal_benchmark_metrics.csv`
- `method_metric_comparison.csv`
- `method_availability.csv`
- `formal_benchmark_report.md`
- `01_formal_benchmark_metric_panel.png`
- `02_formal_benchmark_delta_heatmap.png`
- `03_method_availability.png`
- `figure_package/figure_package_report.md`

### 如何解释 scGen/CPA/GEARS/CellOT

这些方法被作为正式外部比较槽位接入。VKX 不会假装已经运行它们。

如果没有外部预测文件，它们会在 `method_availability.csv` 中标为 `not_run`。这是正确边界，不是失败。

如果已经单独运行这些方法，可以整理成：

```text
method,ko_target,pred_delta_<feature1>,pred_delta_<feature2>,...
scGen,STAT1,...
CPA,STAT1,...
GEARS,STAT1,...
CellOT,STAT1,...
```

然后通过：

```bash
--external-predictions-csv external_method_predictions.csv
```

导入同一套 benchmark 指标和图。

### benchmark 指标

正式 benchmark 不只看 AUC，而是同时输出：

- AUC/ROC：能否识别强响应特征。
- MAE：变化幅度是否接近真实 KO。
- R2：整体 delta 拟合程度。
- direction cosine：KO 方向是否正确。
- feature hit-rate：top 改变特征是否命中。

## 3. 论文级主图整合

`formal-benchmark` 会自动调用 figure package，把主图收进：

```text
results/formal_method_benchmark/figure_package/
```

当前建议论文主图组织为：

- Panel A：方法流程图。
- Panel B：输入模态和状态表示。
- Panel C：正式 benchmark metric panel。
- Panel D：真实 vs 预测 delta heatmap。
- Panel E：方法可用性和外部方法边界。
- Panel F：单细胞状态移动 UMAP/PCA。
- Panel G：failure diagnosis overview。
- Panel H：ATAC peak-level locus/peak/motif 图。

## 4. Hard-constrained residual generator

新增命令：

```bash
python -m vkx.cli train-hard-generator \
  --state-csvs results/dataset1/derived_state_scores.csv,results/dataset2/derived_state_scores.csv \
  --ko-col ko_target \
  --target-kos STAT1,STAT1+JAK2 \
  --prior-dir data/priors \
  --samples-per-ko 300 \
  --max-residual-fraction 0.35 \
  --out-dir results/hard_constrained_generator
```

生成细胞遵循：

```text
virtual cell = control cell + VKX baseline KO delta + bounded learned residual
```

也就是说，generator 不允许自由改变 KO 主方向。它只学习真实 perturbation 细胞围绕 KO 平均方向的残差形状。

如果 PyTorch 可用，训练轻量 residual VAE；如果 PyTorch 不可用，自动退回 PCA residual sampler，并在 `hard_generator_report.md` 中说明。

固定输出：

- `generator_residual_bank.csv`
- `hard_generator_virtual_cells.csv`
- `hard_generator_intervals.csv`
- `hard_generator_metrics.csv`
- `hard_generator_report.md`
- `01_hard_generator_metric_panel.png`
- `02_hard_generator_intervals.png`
- `figure_package/figure_package_report.md`

## 5. 仍然不能夸大的边界

### 5.1 full trimodal labelled benchmark 仍未确认

目前公开数据中比较确定的是：

- RNA+ADT+perturbation labelled benchmark。
- RNA+ATAC+perturbation labelled benchmark。
- RNA+ADT+ATAC 但无 genetic perturbation label 的 DOGMA/TEA-seq 类数据。

真正公开、同一批细胞同时具备 RNA+ADT+ATAC+genetic perturbation label 的 benchmark 仍未确认。因此当前不能把 DOGMA/TEA-seq 当作准确性验证，只能做三模态输入兼容和 reference application。

### 5.2 external deep methods 需要真实预测文件

scGen、CPA、GEARS、CellOT 已经有比较接口，但没有外部预测 CSV 时不能声称已经完成真实比较。

### 5.3 MAPK/TGFB 等非线性 program 仍困难

MAPK/TGFB 类 program 往往涉及强非线性、反馈调控、细胞类型依赖和时间动态。当前模型能提供稳定方向，但复杂幅度和组合效应仍可能不足。

下一步需要：

- 增强 pathway/TF/PPI/motif prior。
- 对 MAPK/TGFB 增加专项非线性校正项。
- 在 failure diagnosis 中继续标出这些 program，避免用户误读。

## 6. 下一步优先级

1. 用 `formal-benchmark` 在 Papalexi、Norman、HMPCITE、ATAC 等统一 state CSV 上跑 VKX vs PLS/ridge/additive。
2. 单独运行或导入 scGen/CPA/GEARS/CellOT 的外部预测 CSV，补齐真实横向比较。
3. 用 `train-hard-generator` 合并更多 perturbation 数据训练 residual generator，并比较 residual VAE 与 PCA fallback。
4. 继续搜索并人工确认 full RNA+ADT+ATAC+genetic perturbation labelled benchmark；没有确认前保持 `not_confirmed_public`。
5. 针对 MAPK/TGFB 和 ATAC peak sparse shape 做专项增强。

## 7. 对外表述建议

英文：

> VKX is a small-sample, multimodal, prior-constrained virtual knockout framework that predicts interpretable cell-state shifts and provides explicit reliability diagnostics. It is designed for settings where labelled perturbation data and compute are limited, and where users need transparent pathway/protein/ATAC-level outputs rather than a black-box free generator.

中文：

> VKX 是一种面向小样本多模态单细胞数据的先验约束虚拟敲除方法。它不依赖大规模预训练和高算力，而是利用真实 perturbation 数据、系统生物学网络先验和多组学状态表示，预测单基因或双基因敲除后的细胞状态变化，并输出可解释的图、表和可靠性诊断。
