# Small-sample Multimodal Virtual Knockout (VKX)

VKX 是一个面向小样本单细胞 perturbation 数据的虚拟敲除方法和软件原型。它的目标不是用大模型自由生成一个看起来像 KO 的细胞，而是在可解释的 pathway/program/protein/ATAC 状态空间中，利用真实 KO 数据、系统网络先验和多模态信息，预测一个基因或两个基因敲除后细胞状态会怎样改变。

## 一句话概括

用户输入普通单细胞矩阵或多组学矩阵，软件内部自动转成可解释的状态表示，然后输出虚拟 KO 之后的 pathway/protein/ATAC 变化、单细胞状态移动图、真实 vs 虚拟 heatmap、ROC/AUC 曲线、KO 总结卡片、失败原因诊断和图文报告。

## 为什么这个方法适合小样本？

VKX 采用 hard-constrained residual/PLS baseline：

- KO 平均方向由真实 perturbation 数据和 Reactome/MSigDB/TF-target/PPI/motif/peak-gene 等先验约束。
- 多模态信息作为额外状态特征加入，例如 RNA pathway score、ADT protein、ATAC gene activity、chromVAR motif activity、peak-level accessibility。
- 双基因 KO 优先使用 interaction residual；数据不足时回退到稳定的 additive/prior-constrained baseline。
- 轻量生成模块只学习 hard constraint 附近的不确定性范围，不让 VAE/flow/diffusion 自由改变 KO 主方向。

这意味着它更保守、更容易解释，也更适合小样本；代价是复杂分布形状和强非线性效应不一定能完全模拟。

## 支持什么输入？

### 1. 带 KO 标签的 perturbation / CRISPR / Perturb-seq 数据

推荐 h5ad：

```text
adata.X                 cells x genes RNA matrix
adata.var_names         gene symbols
adata.obs["ko_target"]  control / STAT1 / STAT1+JAK2 等 KO 标签
adata.obs["cell_type"]  可选
adata.obs["batch"]      可选
adata.obsm["protein"]   可选，ADT/CITE-seq protein
adata.obsm["atac"]      可选，ATAC gene activity
adata.obsm["chromvar"]  可选，TF/motif activity
adata.obsm["peak"]      可选，raw peak count / peak accessibility
```

有真实 KO 标签时，可以做准确性评估，输出 AUC、R2/MAE、真实 vs 虚拟 heatmap 和 UMAP/PCA 状态移动图。

### 2. 普通 10X 单细胞数据

没有 KO 标签也支持，但只能做 prediction-only reference application：

- 可以预测“如果敲 STAT1 或 STAT1+JAK2，细胞状态可能往哪里移动”。
- 不能在该数据内部报告真实准确率、AUC、R2 或 MAE。
- 软件会输出 prior coverage、transfer confidence、uncertainty band、虚拟状态图和 KO 总结卡片。

### 3. 多模态数据

支持 RNA-only、RNA+ADT、RNA+ATAC、RNA+ADT+ATAC 输入。当前公开 labelled benchmark 中已经比较明确的是：

- RNA+ADT+perturbation：例如 ECCITE/Perturb-CITE 类数据。
- RNA+ATAC+perturbation：例如 multiome perturbation 类数据。
- RNA+ADT+ATAC 但无 genetic perturbation 标签：例如 DOGMA/TEA-seq，更适合 prediction-only/reference application。

真正公开、同一批细胞同时具备 RNA+ADT+ATAC+genetic perturbation 标签的数据集，目前在本项目 registry 中仍标记为 `not_confirmed_public`，不能假装已经完成 full trimodal labelled benchmark。

## 最常用命令

### 导入 10X 或 Seurat/AnnData 数据

```bash
python -m vkx.cli import-data \
  --input path/to/filtered_feature_bc_matrix \
  --format 10x_mtx \
  --metadata-csv metadata.csv \
  --cell-id-col cell_id \
  --out-dir results/import_10x
```

### 有 KO 标签时直接评估虚拟敲除

```bash
python -m vkx.cli run \
  --input-h5ad results/import_10x/imported_data.h5ad \
  --ko-col ko_target \
  --target-kos STAT1,STAT1+JAK2 \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv peak_annotation.csv \
  --shape-calibrate quantile \
  --out-dir results/labelled_virtual_ko
```

### 训练 reference model

```bash
python -m vkx.cli train-reference \
  --input-h5ad perturbation_reference.h5ad \
  --ko-col ko_target \
  --batch-col donor \
  --interaction-mode auto \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv peak_annotation.csv \
  --output-model results/reference_models/vkx_reference.pkl
```

### 应用到普通 10X 或无 KO 标签数据

```bash
python -m vkx.cli apply-reference \
  --reference-model results/reference_models/vkx_reference.pkl \
  --input-h5ad ordinary_10x_or_multiome.h5ad \
  --target-kos STAT1,STAT1+JAK2 \
  --cell-type-col cell_type \
  --batch-col donor \
  --uncertainty-method hard-residual \
  --out-dir results/prediction_only_STAT1
```

### 结果诊断：告诉用户哪里可信、哪里要谨慎

```bash
python -m vkx.cli diagnose-results \
  --delta-csv results/labelled_virtual_ko/delta_table.csv \
  --manifest-csv results/labelled_virtual_ko/derived_state_manifest.csv \
  --out-dir results/labelled_virtual_ko/diagnosis
```

`run`、`fit` 和 `apply-reference` 会自动尝试生成诊断结果；这个命令适合对已有结果重新诊断。

### 一键整理成用户可读报告

```bash
python -m vkx.cli summarize-result \
  --result-dir results/labelled_virtual_ko
```

这个命令会把一个结果目录自动整理成：

- `readable_result_report/user_readable_result_report.md`
- `readable_result_report/ko_cards/`
- `readable_result_report/diagnosis/`
- `readable_result_report/figure_package/`

适合把结果发给别人看，也适合写论文或组会汇报前快速检查结果是否可信。

## 主要输出怎么看？

每个结果目录通常包含：

- `01_summary_dashboard.png`: 总体表现，包括方向一致性、误差和改进比例。
- `02_true_vs_virtual_heatmap.png`: 最重要的图，直接比较真实 KO delta、虚拟 KO delta 和误差。
- `03_cell_state_umap.png`: 单细胞状态空间中 control、virtual KO、true KO 的位置。
- `04_auc_strong_response_roc.png`: AUC 曲线，不是柱状图，用来判断能否识别强响应特征。
- `ko_cards/ko_card_<KO>.png`: 每个 KO 的用户可读总结卡片。
- `diagnosis/01_failure_diagnosis_overview.png`: 哪些 KO/特征风险高。
- `diagnosis/02_feature_error_heatmap.png`: feature-level 误差热图。
- `figure_package/figure_package_report.md`: 自动整理好的图文报告。
- `readable_result_report/user_readable_result_report.md`: 一键汇总后的用户可读总报告。

## 怎么做横向比较？

```bash
python -m vkx.cli method-comparison \
  --metric-csv results/vkx_metrics.csv,results/ridge_metrics.csv,results/gears_metrics.csv \
  --out-dir results/method_comparison
```

VKX 的定位是小样本、多模态、先验约束、可解释 baseline。它应该和 ridge/PLS、scGen、CPA、GEARS、CellOT、diffusion/flow 类方法比较，但结论需要限定在具体 benchmark 和数据规模内。

## 当前还缺什么？

1. 需要更正式的横向 benchmark，把 VKX 与 ridge/PLS、scGen、CPA、GEARS、CellOT 等方法在同一批数据上比较。
2. 真正公开 RNA+ADT+ATAC+genetic perturbation labelled benchmark 仍未确认；找到后才能升级为 full trimodal labelled benchmark。
3. 当前 VAE/flow/diffusion 入口仍以 hard-constrained residual uncertainty 为主，完整神经生成模型需要更多同类型 perturbation 数据后再训练。
4. MAPK/TGFB 等强非线性 program 仍是困难案例，需要更强 pathway/TF/PPI/motif 先验和非线性修正。
5. 可视化已经能自动生成，但论文级多 panel figure 还需要针对最终 benchmark 再统一排版。

## 结果解释原则

- 有 KO 标签：可以报告真实准确性，例如 AUC、MAE、R2、方向一致性、真实 vs 虚拟 heatmap。
- 没有 KO 标签：只能报告预测状态变化、先验覆盖、迁移置信度和不确定性范围，不能报告真实准确性。
- 多模态输入通常会让预测更可靠，但前提是模态质量好、与 KO 机制相关，并且有合理的 feature selection 和先验约束。
- AUC 如果接近 1，需要检查特征数、阈值和正负样本数量；特征很少时 AUC 容易显得过于理想，必须谨慎解释。
