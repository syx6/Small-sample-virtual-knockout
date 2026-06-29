# Reference model v2 工作流

这份文档对应当前新增的可复用软件接口，重点解决五件事：

1. reference model 版本管理
2. 批量 KO 应用
3. cell type 分层输出
4. 普通 10X / DOGMA / TEA-seq / multiome 的 prediction-only 报告
5. 双敲 interaction、motif/TF-target 加权和 MAPK/TGFB 专项先验

## 1. 训练 reference model

从带真实 perturbation 标签的数据训练 reference model：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli train-reference `
  --input-h5ad data\hmpcite_gse243244\hmpcite_perturbation_rna_adt_doubleko.h5ad `
  --ko-col ko_target `
  --prior-dir data\priors `
  --output-model results\reference_models\hmpcite_rna_adt_reference.pkl `
  --dataset-name "HMPCITE-seq RNA+ADT perturbation reference" `
  --extra-obsm protein:protein `
  --max-pathways 40
```

输出：

```text
hmpcite_rna_adt_reference.pkl
hmpcite_rna_adt_reference.pkl.metadata.json
```

metadata JSON 会记录：

- reference 数据集名称
- 模型版本
- 训练 KO 标签数
- 训练基因数
- 状态特征数
- 是否支持 batch KO、double KO、cell type stratification、prediction-only application

## 2. 批量应用到普通 10X 或 multiome 数据

应用前可以先检查 reference model：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli inspect-reference `
  --reference-model results\reference_models\hmpcite_rna_adt_reference.pkl `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --out-dir results\hmpcite_reference_inspection
```

输出：

```text
reference_summary.csv
reference_prior_libraries.csv
reference_state_features.csv
reference_training_genes.csv
reference_target_prior_coverage.csv
reference_inspection_report.md
```

这个步骤回答：

```text
这个 reference model 学过多少 KO？
有哪些状态特征？
用了哪些 prior library？
我要预测的 KO 在 pathway/TF/PPI/motif 先验里覆盖好不好？
```

一个命令可以同时预测多个单敲和双敲：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\hmpcite_rna_adt_reference.pkl `
  --input-h5ad your_10x_or_multiome.h5ad `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --out-dir results\your_prediction_only_virtual_ko `
  --max-cells 1000
```

输出：

```text
applied_virtual_cells.csv
predicted_ko_delta.csv
target_interpretation.csv
transfer_confidence.csv
prior_coverage.csv
apply_report.md
prediction_only_report.md
01_predicted_ko_delta_heatmap.png
02_input_vs_virtual_pca.png
03_transfer_confidence.png
05_prior_coverage.png
```

解释：

- `predicted_ko_delta.csv`：每个 KO 的预测状态变化。
- `applied_virtual_cells.csv`：输入细胞和虚拟 KO 细胞。
- `target_interpretation.csv`：说明每个 target 是单敲、双敲还是探索性多基因 KO。
- `transfer_confidence.csv`：说明 KO gene 是否在 reference 中见过，输入细胞是否偏离 reference 分布。
- `prior_coverage.csv`：说明每个 KO 命中了多少 pathway/TF/PPI/motif prior。

如果某个 KO 是 unseen gene，而且 `prior_coverage.csv` 里 `n_prior_terms_hit` 接近 0，那么这个预测应标记为低置信度。

## 3. Cell type 分层输出

如果 h5ad 里有 cell type 注释，例如 `adata.obs["cell_type"]`：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\hmpcite_rna_adt_reference.pkl `
  --input-h5ad your_10x_or_multiome.h5ad `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --cell-type-col cell_type `
  --out-dir results\your_celltype_virtual_ko
```

额外输出：

```text
cell_type_predicted_delta.csv
04_cell_type_predicted_delta_heatmap.png
```

这个图回答的问题是：

```text
同一个虚拟 KO 在不同 cell type 里，预测影响哪些 pathway/protein/ATAC state？
```

注意：如果输入数据没有真实 KO 标签，这仍然是 prediction-only，不是准确性验证。

## 4. 真正 multiome 示例怎么运行

如果数据是 RNA + ADT + ATAC 且带 perturbation 标签，可以直接评估：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad your_rna_adt_atac_perturbation.h5ad `
  --ko-col ko_target `
  --target-kos GENE1+GENE2 `
  --prior-dir data\priors `
  --out-dir results\trimodal_true_ko_eval `
  --extra-obsm protein:protein,atac:atac,chromvar:tf `
  --max-extra-features-per-obsm 100 `
  --extra-feature-selection hybrid
```

如果数据是 DOGMA/TEA-seq 这类 RNA + protein + ATAC，但没有 perturbation 标签：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\hmpcite_rna_adt_reference.pkl `
  --input-h5ad your_unlabeled_dogma_or_teaseq.h5ad `
  --target-kos STAT1,JAK2,STAT1+JAK2 `
  --cell-type-col cell_type `
  --out-dir results\dogma_teaseq_prediction_only
```

这时输出的是虚拟状态变化报告，不报告 AUC/R2/MAE。

## 5. 双敲 interaction 怎么用

有真实 double-KO 标签时，用专门命令评估 interaction residual：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli double-interaction `
  --delta-csv results\your_doubleko_delta_table.csv `
  --ko-col ko_genes `
  --n-ko-col n_ko_genes `
  --target-prefix delta_program_ `
  --prior-dir data\priors `
  --out-dir results\your_double_interaction_eval
```

当前 interaction residual 已加入：

- gene-gene shared prior terms
- weighted TF-target/motif prior
- MAPK/TGFB/SMAD/RAS/ERK/JNK/P38 等专项机制加权

它适合回答：

```text
双敲是不是比两个单敲简单相加更接近真实 KO？
```

## 6. MAPK/TGFB 专项修正

MAPK/TGFB 这类 program 难预测，主要原因是：

- 通路高度上下文依赖。
- 双基因组合经常有非线性。
- 单纯 PLS/residual 容易低估复杂交互。

当前优化包括：

- pathway 选择时优先保留 MAPK/TGFB/SMAD/MAPK cascade 相关 term。
- prior weight 中对 MAPK/TGFB/SMAD/RAS/ERK/JNK/P38 term 加权。
- double-interaction 中对这些 term 的 shared pathway interaction 加权。

这不会保证 MAPK/TGFB 一定变好，但会让模型在这些失败 program 上有更强的生物学约束。

## 7. 对普通用户怎么解释

普通 10X 或无标签 multiome 的输出应该这样解释：

```text
这是模型根据 perturbation reference 推断出的虚拟 KO 状态变化。
因为输入数据没有真实 KO 标签，所以不能在这个数据内部证明准确率。
请主要看 predicted delta heatmap、cell type 分层变化和 transfer confidence。
```

不要写成：

```text
模型已经在你的普通 10X 数据里验证了 KO 准确率。
```
