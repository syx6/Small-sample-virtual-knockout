# 五点优化后的当前状态

## 1. 结果模式分层

现在每个主要输出目录都会生成 `analysis_mode.md`，用于告诉用户当前结果属于哪一种模式：

- `evaluation`：输入数据有真实 KO 标签，可以解释 heatmap、UMAP、AUC 和误差指标。
- `prediction_only`：把 reference model 应用到普通 10X 或无标签细胞，只能解释预测状态变化，不能解释为真实准确率。
- `state_scoring_only`：只把普通单细胞矩阵转成 pathway/program/protein state score。
- `double_ko_evaluation`：有真实双敲标签，可以比较 additive baseline 和 interaction residual model。

这样可以避免普通 10X 无标签数据被误读成已经完成真实准确性验证。

## 2. 普通 10X 输入

普通用户可以输入 `.h5ad` 单细胞矩阵：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli score `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --prior-dir data\priors `
  --out-dir results\ordinary_10x_like_input_demo `
  --protein-obsm protein `
  --max-pathways 24
```

输出：

- `derived_state_scores.csv`
- `derived_state_manifest.csv`
- `analysis_mode.md`

注意：如果没有 KO 标签，这一步不能算 AUC/R2/MAE，只是把普通细胞准备成虚拟敲除可用的状态表示。

## 3. 多模态单敲/双敲应用

已训练新的 reference model：

```text
results/reference_models/papalexi_rna_protein_reference_v2.pkl
```

它来自 Papalexi RNA+ADT 多模态 perturbation 数据，包含 RNA pathway score 和 ADT/protein features。

应用到输入细胞并同时预测单敲和双敲：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\papalexi_rna_protein_reference_v2.pkl `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --target-kos STAT1,STAT1+JAK2 `
  --out-dir results\multimodal_single_double_apply_v2_demo `
  --max-cells 500
```

输出：

- `predicted_ko_delta.csv`
- `applied_virtual_cells.csv`
- `transfer_confidence.csv`
- `01_predicted_ko_delta_heatmap.png`
- `02_input_vs_virtual_pca.png`
- `03_transfer_confidence.png`
- `analysis_mode.md`

这个结果说明软件支持多模态输入，也支持一个或两个基因的虚拟敲除请求。

但是：Papalexi 当前公开小样本例子没有真实双敲标签，所以 `STAT1+JAK2` 在这里属于 prediction-only 应用，不能在该数据内部报告真实 AUC/R2/MAE。

## 4. Reference transfer 置信度

`apply-reference` 现在会输出 `transfer_confidence.csv`。

当前示例：

```text
STAT1       high
STAT1+JAK2  high
```

置信度依据：

- 输入细胞有多少比例超出 reference 训练分布的 q95 距离。
- 请求敲除的基因是否出现在 reference 训练 KO 标签中。

这不是准确率，而是“这个 reference model 是否适合迁移到当前输入数据”的提示。

## 5. 双敲 interaction model 正式接口

双敲非线性模型已从临时脚本整理成正式 CLI：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli double-interaction `
  --delta-csv results\norman_program_delta.csv `
  --ko-col ko_genes `
  --n-ko-col n_ko_genes `
  --target-prefix delta_program_ `
  --prior-dir data\priors `
  --out-dir results\software_interface_double_interaction_cli_demo
```

输出：

- `double_interaction_metrics.csv`
- `double_interaction_predictions.csv`
- `double_interaction_metrics.png`
- `double_interaction_report.md`
- `analysis_mode.md`

Norman 52 个双敲组合的结果：

```text
single_gene_additive: mean MAE 0.150, mean R2 0.008, mean ROC-AUC 0.707
interaction_residual: mean MAE 0.076, mean R2 0.617, mean ROC-AUC 0.894
```

说明加入系统网络先验的 interaction residual 后，双敲非线性预测明显改善。

## 仍然缺什么

目前仍然缺一个真正理想的公开小样本数据：

```text
同一批细胞同时具备：
RNA + ADT 或 RNA + ATAC 多模态
真实 single KO 标签
真实 double KO 标签
```

没有这种数据时，我们可以：

- 在 Papalexi 这类 RNA+ADT 数据上验证多模态单敲。
- 在 Norman 这类 RNA-only perturb-seq 上验证真实双敲。
- 在多模态数据上应用双敲预测，但不能声称已经在该数据中真实验证双敲准确性。

下一步真正提升说服力的关键，是接入一个多模态 perturbation double-KO 数据，或由用户提供自己的 RNA+ADT/RNA+ATAC perturbation 数据。
