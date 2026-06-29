# 虚拟敲除软件接口说明：用户输入原始矩阵，软件内部自动转成通路状态

## 1. 先澄清一个关键点

用户不需要自己准备 pathway/program score。

真正给用户的输入应该是常见单细胞文件：

- 单细胞 RNA 表达矩阵
- CITE-seq / ADT protein 矩阵
- scATAC / gene activity / peak score 矩阵
- 每个细胞对应的 KO / perturbation 标签

`pathway/program score` 是软件内部自动生成的中间表示，不是要求用户手工提供的输入。

也就是说，软件逻辑应该是：

```text
用户输入：原始单细胞/组学矩阵 + KO 标签
        ↓
软件内部：RNA -> pathway/program score
        ↓
软件内部：protein/ATAC -> phenotype/regulatory score
        ↓
模型输入：解释性 state score
        ↓
输出：虚拟 KO 单细胞状态 + heatmap + UMAP + AUC 曲线
```

## 2. 推荐用户输入格式

### 格式 0：Seurat / 10x 先导入

如果用户手里不是 h5ad，而是 10x 或 Seurat 文件，可以先运行：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli import-data `
  --input path\to\filtered_feature_bc_matrix `
  --format 10x_mtx `
  --out-dir results\import_10x_demo
```

可选格式：

```text
h5ad
10x_mtx
10x_h5
h5seurat
```

导入结果：

```text
imported_data.h5ad
input_summary.csv
input_overview.png
import_report.md
```

其中 `input_overview.png` 是第一张必须给用户看的图：它告诉用户数据里有多少细胞、多少 RNA features、是否检测到 protein/peak/guide 模态，以及 metadata 标签分布。

### 格式 A：AnnData h5ad

这是最推荐的格式。

要求：

- `adata.X`: cells x genes 的 RNA 表达矩阵
- `adata.var_names`: gene symbols，例如 `STAT1`, `JAK2`, `IFNGR2`
- `adata.obs[ko_col]`: 每个细胞的 KO 标签
- 可选：`adata.obsm["protein"]` 或其他 obsm 矩阵，存放 ADT/protein/ATAC/gene activity

例子：

```text
data/papalexi_small_pathway.h5ad
```

### 格式 B：CSV 原始矩阵

也可以输入 CSV。

要求：

| 列 | 说明 |
|---|---|
| `ko_target` | 每个细胞的 KO 标签 |
| gene columns | 每个基因一列，列名是 gene symbol，值是表达量 |

示意：

| ko_target | STAT1 | JAK2 | IRF1 | ... |
|---|---:|---:|---:|---:|
| control | 0.1 | 1.3 | 0.4 | ... |
| STAT1 | 0.0 | 1.2 | 0.2 | ... |

## 3. 软件内部自动做什么

用户输入原始 RNA 后，软件自动：

1. 读取 gene expression matrix。
2. 读取 Reactome / MSigDB / TF-target / PPI 等 GMT 先验。
3. 找出和当前数据 gene symbols 有交集的 pathway/program。
4. 对每个细胞计算 pathway/program score。
5. 如果提供 protein/ATAC 矩阵，就把它一起加入 state representation。
6. 用这些解释性 state score 训练 prior-constrained residual/PLS 虚拟 KO 模型。

因此，文档里出现的 pathway/program score 指的是：

```text
模型内部状态表示，不是用户必须准备的原始输入。
```

## 4. 从原始 h5ad 一键运行

示例命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --ko-col ko_target `
  --target-kos STAT1 `
  --prior-dir data\priors `
  --out-dir results\software_interface_raw_papalexi `
  --dataset-name "Papalexi ECCITE-seq" `
  --modality "raw RNA matrix + ADT protein obsm" `
  --representation "auto-derived pathway/protein scores" `
  --protein-obsm protein `
  --calibrate none `
  --max-pathways 24
```

这条命令的含义是：

- 输入的是 `.h5ad` 原始矩阵。
- KO 标签在 `adata.obs["ko_target"]`。
- RNA pathway score 由软件自动计算。
- protein 矩阵从 `adata.obsm["protein"]` 读取。
- 默认虚拟敲除一个基因：`STAT1`。

如果想虚拟敲除两个基因，写成：

```powershell
--target-kos STAT1+JAK2
```

如果是在有真实 KO 标签的数据里做批量评估，才写成：

```powershell
--target-kos STAT1,JAK2,IFNGR2,IRF1
```

## 5. 从原始 CSV 一键运行

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-csv your_raw_matrix.csv `
  --ko-col ko_target `
  --target-kos STAT1 `
  --prior-dir data\priors `
  --out-dir results\your_dataset_result `
  --dataset-name "Your dataset" `
  --modality "raw RNA matrix" `
  --representation "auto-derived pathway/program scores" `
  --calibrate auto
```

RNA-only 数据推荐 `--calibrate auto`，因为 RNA-only 常见问题是方向还可以，但变化幅度需要校准。

RNA + protein / multiome 数据建议先用 `--calibrate none`，因为多模态状态本身已经提供了更强约束。

## 6. 普通 10X 单细胞数据怎么支持

大多数用户手里是普通 10X scRNA-seq，不一定有 perturb-seq 标签。这个场景要分清楚：

### 6.1 只有普通 10X，没有 KO 标签

可以支持“状态准备”，也就是把原始 RNA 矩阵自动转成 pathway/program score：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli score `
  --input-h5ad your_10x_data.h5ad `
  --prior-dir data\priors `
  --out-dir results\your_10x_state_scores
```

输出：

- `derived_state_scores.csv`
- `derived_state_manifest.csv`

这一步适合把普通 10X 数据变成后续虚拟 KO 可用的细胞状态。

但是，如果完全没有真实 KO/perturbation 标签，软件不能在这个数据集内部评估“虚拟 KO 是否接近真实 KO”，因为没有真实 KO 作为答案。

### 6.2 有 KO 标签的 10X

有 KO 标签的 10X 是支持的。典型文件结构是：

```text
filtered_feature_bc_matrix/
├── matrix.mtx.gz
├── features.tsv.gz
└── barcodes.tsv.gz

metadata.csv
```

`metadata.csv` 用 10X barcode 对齐每个细胞的 KO 标签：

```csv
cell_id,ko_target,cell_type,batch
AAACCCAAGAAACACT-1,control,T_cell,batch1
AAACCCAAGAAACCAT-1,STAT1,T_cell,batch1
AAACCCAAGAAAGTGG-1,JAK2,Mono,batch2
AAACCCAAGAAATCCA-1,STAT1+JAK2,Mono,batch2
```

关键列：

| 列名 | 作用 |
|---|---|
| `cell_id` | 10X barcode，必须能和 `barcodes.tsv.gz` 对齐 |
| `ko_target` | 每个细胞的扰动标签，例如 `control`、`STAT1`、`STAT1+JAK2` |
| `cell_type` | 可选，用于 cell type 分层输出 |
| `batch` | 可选，用于后续 batch covariate 分析 |

导入：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli import-data `
  --input path\to\filtered_feature_bc_matrix `
  --format 10x_mtx `
  --metadata-csv metadata.csv `
  --cell-id-col cell_id `
  --out-dir results\import_10x_ko
```

评估真实 KO：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad results\import_10x_ko\imported_data.h5ad `
  --ko-col ko_target `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --prior-dir data\priors `
  --out-dir results\tenx_ko_virtual_eval
```

如果 `ko_target` 里有 control 和真实 KO 标签，这就是 evaluation 模式，可以报告：

- true vs virtual heatmap
- UMAP
- ROC/AUC 曲线
- distribution improvement
- direction cosine
- MAE/R2 等误差指标

### 6.3 普通 10X 作为待预测细胞群

这是合理应用场景：

```text
普通 10X control cells
+ 外部 perturbation/reference model 学到的 KO delta
= 这个 10X 数据上的虚拟 KO 细胞状态
```

也就是说，普通 10X 可以作为“要被虚拟敲除的细胞群”。但 KO 方向需要来自：

- 同实验里的少量 perturbation 数据，或者
- 外部 perturb-seq / CRISPR reference 数据，或者
- 以后保存好的 reference virtual KO model。

当前版本已经支持普通 10X 的自动状态转换，也已经把 reference model 的训练、保存、加载和应用做成正式接口。

reference model v2 还新增了：

- 模型 metadata JSON，记录训练 KO、训练基因、状态特征和支持能力。
- `inspect-reference`，应用前检查 reference model 和目标 KO 的 prior 覆盖。
- 批量 KO 输入，例如 `STAT1,JAK2,STAT1+JAK2`。
- prediction-only 报告，明确普通 10X/无标签 multiome 不能报告真实准确率。
- 可选 cell type 分层输出。

训练 reference model：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli train-reference `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --ko-col ko_target `
  --prior-dir data\priors `
  --output-model results\reference_models\papalexi_rna_protein_reference.pkl `
  --dataset-name "Papalexi ECCITE-seq reference" `
  --protein-obsm protein `
  --max-pathways 24
```

检查 reference model：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli inspect-reference `
  --reference-model results\reference_models\papalexi_rna_protein_reference.pkl `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --out-dir results\papalexi_reference_inspection
```

应用到普通 h5ad 细胞：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\papalexi_rna_protein_reference.pkl `
  --input-h5ad your_10x_data.h5ad `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --cell-type-col cell_type `
  --out-dir results\your_10x_virtual_ko_batch
```

输出包括：

- `applied_virtual_cells.csv`
- `predicted_ko_delta.csv`
- `target_interpretation.csv`
- `prior_coverage.csv`
- `transfer_confidence.csv`
- `prediction_only_report.md`
- `01_predicted_ko_delta_heatmap.png`
- `02_input_vs_virtual_pca.png`
- `03_transfer_confidence.png`
- `05_prior_coverage.png`
- `04_cell_type_predicted_delta_heatmap.png`，如果提供 `--cell-type-col`
- `apply_report.md`

注意：`apply-reference` 是 prediction-only 模式。如果输入数据没有真实 KO 标签，不会输出 AUC、distribution improvement 或 true-vs-virtual heatmap。

## 7. 输出文件

### 7.1 Benchmark readiness 输出

如果 RNA 和 ATAC 是分开的矩阵，先组装：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli assemble-multiome `
  --rna-input path\to\rna_matrix `
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

输出：

| 文件 | 作用 |
|---|---|
| `multiome_assembly_summary.csv` | RNA/ATAC 共享细胞数、特征数和 KO 标签数 |
| `multiome_assembly_overview.png` | 直观显示 RNA/ATAC 模态和 KO 标签分布 |
| `multiome_assembly_report.md` | 说明组装结果和下一步命令 |

在正式运行多模态 benchmark 前，推荐先运行：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli validate-benchmark `
  --input-h5ad your_multimodal_perturbation.h5ad `
  --ko-col ko_target `
  --extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak `
  --out-dir results\benchmark_readiness
```

输出：

| 文件 | 作用 |
|---|---|
| `benchmark_readiness.csv` | 检查细胞数、RNA features、KO 标签、control、extra modality 是否存在 |
| `benchmark_label_counts.csv` | 每个 KO/perturbation 标签有多少细胞 |
| `benchmark_overview.png` | 一页图显示 benchmark 是否合格和主要 KO 标签 |
| `benchmark_modalities.png` | 显示 RNA/protein/ATAC/chromVAR/peak 等模态特征数 |
| `benchmark_readiness_report.md` | 自动解释这个数据是否适合作为真实 benchmark |

如果 readiness 是 `ok`，说明它有真实 KO 标签和至少一种额外模态，可以继续跑 `run` 并解释 AUC/heatmap/UMAP。  
如果是 `partial`，说明可以运行部分功能，但不能宣传成完整多模态 perturbation benchmark。

运行后会生成：

| 文件 | 作用 |
|---|---|
| `derived_state_scores.csv` | 软件自动从原始矩阵派生出的 pathway/protein/ATAC 状态表 |
| `derived_state_manifest.csv` | 每个 state feature 来自哪里，例如 RNA pathway 或 protein obsm |
| `report.md` | 自动文字报告 |
| `summary.csv` | 总体效果指标 |
| `metrics.csv` | 每个 KO、每个 feature 的分布距离 |
| `delta_table.csv` | 每个 KO 的真实变化和虚拟变化 |
| `virtual_cells.csv` | 虚拟 KO 单细胞状态 |
| `auc_summary.csv` | 强响应识别 AUC |
| `calibration.csv` | 每个 feature 的校准倍率 |
| `umap_cells.csv` | UMAP 坐标 |

## 8. 默认主图

### 01_summary_dashboard.png

一页总览：

- `Distribution improvement > 0`: 虚拟 KO 比 control 更接近真实 KO。
- `Improved features`: 有多少比例的状态特征被改善。
- `Direction cosine`: KO 方向是否预测对。
- `Magnitude error`: KO 变化幅度误差。

### 02_true_vs_virtual_heatmap.png

最重要的解释图：

```text
左：真实 KO 相对 control 的变化
中：虚拟 KO 相对 control 的变化
右：虚拟 - 真实 的误差
```

如果左图和中图颜色方向一致，说明 KO 方向预测对了。
如果右图颜色浅，说明变化幅度也接近真实 KO。

### 03_cell_state_umap.png

展示单细胞状态移动：

```text
control cells -> virtual KO cells -> true KO cells
```

理想情况是 virtual KO 细胞云团从 control 向 true KO 靠近。

### 04_auc_strong_response_roc.png

这个图回答一个简单问题：

```text
模型能不能识别哪些 pathway/protein/program 会发生强 KO 响应？
```

AUC 越高，说明强响应排序能力越好。

## 9. 原始 h5ad 批量评估示例结果

示例输出目录：

```text
results/software_interface_raw_papalexi
```

结果：

| 指标 | 值 |
|---|---:|
| mean distribution improvement | 0.111 |
| improved fraction | 0.759 |
| mean direction cosine | 0.679 |
| mean abs delta error | 0.317 |
| strong-response ROC-AUC | 0.802 |

解释：

- 输入是原始 RNA 矩阵和 protein obsm，不是用户手工准备的 pathway 表。
- 软件自动生成 28 个状态特征，其中包括 RNA pathway score 和 protein score。
- 约 75.9% 的状态特征在虚拟 KO 后比 control 更接近真实 KO。
- AUC 为 0.802，说明模型能较好识别强响应特征。

## 10. 什么时候还可以直接输入 state score table

如果高级用户已经有自己计算好的 pathway/program/protein/ATAC score，也可以用旧的 `fit` 子命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli fit `
  --state-csv your_state_score_table.csv `
  --ko-col ko_target `
  --target-kos STAT1 `
  --prior-dir data\priors `
  --out-dir results\your_state_score_result
```

但这不是普通用户的主要入口。普通用户应该使用：

```text
python -m vkx.cli run
```

## 11. 当前接口边界

已经支持：

- 原始 h5ad RNA 矩阵输入。
- 原始 CSV gene expression 矩阵输入。
- 普通 10X h5ad 的 pathway/program state score 自动转换。
- 可选 protein/ADT obsm 输入。
- 可选 ATAC/gene activity/chromVAR/peak obsm 输入。
- 自动 RNA pathway/program score。
- 自动结果图：summary、heatmap、UMAP、AUC。
- reference model 训练和应用。
- 批量单敲/双敲 prediction-only 应用。
- cell type 分层 prediction-only 输出。

下一步需要继续扩展：

- 多基因 KO 组合效应的主接口深度集成。
- 真正 RNA+ADT+ATAC 且带 perturbation 标签的公开 benchmark。
- ATAC peak-level 调控先验继续增强。
- batch covariate 的显式建模。
- 从 Seurat/10x 等格式直接导入。
