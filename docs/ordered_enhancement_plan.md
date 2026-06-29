# 按优先级增强路线

这份文档按当前开发顺序整理后续增强，每一步都要求输出直观图，而不是只输出表格。

## 1. Seurat / 10x 直接导入

已新增命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli import-data `
  --input path\to\filtered_feature_bc_matrix `
  --format 10x_mtx `
  --out-dir results\import_10x_demo
```

支持格式：

| 格式 | 参数 |
|---|---|
| AnnData | `--format h5ad` |
| 10x mtx 文件夹 | `--format 10x_mtx` |
| 10x filtered_feature_bc_matrix.h5 | `--format 10x_h5` |
| h5Seurat 基础读取 | `--format h5seurat` |

如果有 KO 标签、cell type 或 batch 信息，可以用 metadata CSV 合并：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli import-data `
  --input path\to\filtered_feature_bc_matrix.h5 `
  --format 10x_h5 `
  --metadata-csv metadata.csv `
  --cell-id-col cell_id `
  --out-dir results\import_10x_with_metadata
```

固定输出：

```text
imported_data.h5ad
input_summary.csv
input_overview.png
import_report.md
```

`input_overview.png` 会展示：

- 细胞数
- RNA feature 数
- protein / peak / guide 等 obsm 模态
- metadata 中主要标签分布

## 2. 真正 RNA+ADT+ATAC 且带 perturbation 标签的公开 benchmark

当前原则：

```text
只有同时具备多模态矩阵和 perturbation/guide 标签的数据，才能作为准确性 benchmark。
无 perturbation 标签的 DOGMA/TEA-seq/Multiome 只能做 prediction-only application。
```

候选数据：

| 数据 | 模态 | 标签 | 当前用途 |
|---|---|---|---|
| CAT-ATAC | RNA + ATAC + CRISPR guide | 有 perturbation identity | 真 multiome perturbation benchmark 候选 |
| scPerturb ATAC files | ATAC gene_scores + chromVAR + peak_bc + perturbation | 有 perturbation labels | ATAC regulatory benchmark |
| HMPCITE-seq | RNA + ADT + guide labels | 有 perturbation labels | RNA+ADT double-KO benchmark |
| DOGMA/TEA-seq | RNA + protein + ATAC | 通常无 perturbation labels | prediction-only 三模态兼容测试 |

来源：

- CAT-ATAC：Cell Reports Methods / PubMed 页面说明该方法在 10x Multiome 基础上加入 CRISPR gRNA capture，可同时获得 RNA、ATAC 和 perturbation assignment。
  https://pubmed.ncbi.nlm.nih.gov/41218606/
- scPerturb ATAC files：Figshare 页面说明提供 ChromVar、LSI_embedding、gene_scores、markerpeak_target、peak_bc 等 ATAC feature matrices。
  https://plus.figshare.com/articles/dataset/scPerturb_Single-Cell_Perturbation_Data_ATAC_files/24160968
- scPerturb resource：Nature Methods 页面说明 scPerturb 收集了包含 transcriptomics、proteomics、epigenomics readouts 的单细胞 perturbation-response 数据集。
  https://www.nature.com/articles/s41592-023-02144-y

下一步接入规则：

1. 如果 RNA 和 ATAC 是分开的矩阵，先用 `assemble-multiome` 组装成 h5ad：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli assemble-multiome `
  --rna-input path\to\gex_matrix `
  --rna-format 10x_mtx `
  --atac-input path\to\atac_matrix `
  --atac-format 10x_mtx `
  --metadata-csv metadata_with_ko.csv `
  --cell-id-col cell_id `
  --ko-col ko_target `
  --max-atac-features 500 `
  --output-h5ad data\public_multiome_perturbation.h5ad `
  --out-dir results\public_multiome_assembly
```

固定输出：

```text
multiome_assembly_summary.csv
multiome_assembly_overview.png
multiome_assembly_report.md
```

2. 如果数据已经是整合好的 h5ad，可以跳过组装。
3. h5ad 内部推荐结构：

```text
adata.X                  RNA matrix
adata.var_names          gene symbols
adata.obs["ko_target"]   perturbation labels
adata.obsm["protein"]    optional ADT/protein
adata.obsm["atac"]       ATAC/gene activity/selected peak matrix
adata.obsm["chromvar"]   optional motif activity
adata.obsm["peak"]       optional peak-level matrix
```

4. 检查 benchmark readiness：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli validate-benchmark `
  --input-h5ad data\public_multiome_perturbation.h5ad `
  --ko-col ko_target `
  --extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak `
  --out-dir results\public_trimodal_readiness
```

5. 如果 readiness 是 `ok`，再运行正式 benchmark：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\public_multiome_perturbation.h5ad `
  --ko-col ko_target `
  --target-kos GENE1,GENE2,GENE1+GENE2 `
  --prior-dir data\priors `
  --extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak `
  --max-extra-features-per-obsm 100 `
  --extra-feature-selection hybrid `
  --out-dir results\public_trimodal_benchmark
```

固定输出：

```text
benchmark_readiness.csv
benchmark_label_counts.csv
benchmark_overview.png
benchmark_modalities.png
benchmark_readiness_report.md
```

必须输出图：

- `01_summary_dashboard.png`
- `02_true_vs_virtual_heatmap.png`
- `03_cell_state_umap.png`
- `04_auc_strong_response_roc.png`
- `05_atac_peak_level_changes.png`，如果有 peak 特征

## 3. ATAC peak-level 调控先验继续增强

当前 v2 已完成：

- `peak:peak` obsm 输入
- locus-aware peak selection
- `05_atac_peak_level_changes.png`
- motif/TF-target weighted prior
- target locus + marker peak + KO effect + accessibility + motif/TF prior 综合打分
- `selection_reason` 和 peak regulatory score 输出
- variance shape calibration
- zero-inflated / quantile shape calibration，用于记录 peak 开放比例并校准分位数形状；对原始非负 peak/count 才启用 hard-zero open/closed 约束

新增输出：

```text
data/scperturb_atac/liscovitch_k562_selected_peak_metadata.csv
data/scperturb_atac/liscovitch_k562_peak_regulatory_prior_scores_top1000.csv
results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a/05_atac_peak_level_changes.png
```

结果：

```text
KDM6A regulatory peak prior + quantile shape calibration:
ROC-AUC = 0.674
direction cosine = 0.771
distribution improvement = 0.166
improved fraction = 0.788
```

下一步继续增强：

- motif-to-peak annotation。
- promoter/enhancer 更精细分类。
- peak-gene linkage。
- batch-aware peak normalization。

## 4. 多基因 KO 组合效应主接口集成

当前已有：

- `double-interaction` 独立命令
- interaction residual 评估
- seen/unseen gene pair 分层

下一步：

- 在 `run` 和 `apply-reference` 中自动识别双敲。
- 有真实 double-KO 标签时自动调用 interaction benchmark。
- 无真实标签时输出 prediction-only double-KO interpretation。

必须输出图：

- additive vs interaction MAE/R2/AUC 图
- 多基因组合 true vs virtual heatmap
- seen/unseen gene pair 分层图

## 5. Batch covariate 显式建模

当前 transfer confidence 会检查输入细胞是否偏离 reference 分布，但还没有显式 batch covariate。

下一步：

- `run` / `train-reference` 增加 `--batch-col`。
- 在计算 KO delta 前做 batch-aware centering。
- 输出 batch composition 图。
- 在 prediction-only report 里说明输入 batch 是否偏离 reference。

必须输出图：

- batch composition barplot
- batch-wise state shift heatmap
- batch-aware vs uncorrected comparison
