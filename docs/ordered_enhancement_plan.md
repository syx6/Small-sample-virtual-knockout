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

1. 用 `import-data` 或专门脚本转成 h5ad。
2. RNA 放入 `adata.X`。
3. ADT/protein 放入 `adata.obsm["protein"]`。
4. ATAC gene activity / chromVAR / peak score 放入 `adata.obsm["atac"]`、`adata.obsm["chromvar"]`、`adata.obsm["peak"]`。
5. perturbation 标签放入 `adata.obs["ko_target"]`。
6. 运行：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\public_trimodal_perturbation.h5ad `
  --ko-col ko_target `
  --target-kos GENE1,GENE2,GENE1+GENE2 `
  --prior-dir data\priors `
  --extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak `
  --max-extra-features-per-obsm 100 `
  --extra-feature-selection hybrid `
  --out-dir results\public_trimodal_benchmark
```

必须输出图：

- `01_summary_dashboard.png`
- `02_true_vs_virtual_heatmap.png`
- `03_cell_state_umap.png`
- `04_auc_strong_response_roc.png`
- `05_atac_peak_level_changes.png`，如果有 peak 特征

## 3. ATAC peak-level 调控先验继续增强

当前已有：

- `peak:peak` obsm 输入
- locus-aware peak selection
- `05_atac_peak_level_changes.png`
- motif/TF-target weighted prior

下一步增强：

- peak 按 target gene locus、TF motif、markerpeak_target、KO effect 共同打分。
- 输出每个 peak 的 locus 类型：promoter、intronic、distal、global variable。
- peak 图固定展示真实 KO、虚拟 KO、误差和单细胞分布。

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
