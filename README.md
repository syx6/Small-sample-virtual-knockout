# 小样本多模态虚拟敲除

这是一个面向小样本 perturbation / CRISPR 单细胞数据的虚拟敲除原型。

关键设计原则：

- 用户输入普通单细胞或多组学矩阵，不需要手工准备通路分数。
- 软件内部自动把 RNA 转成 pathway/program score。
- protein/ADT/ATAC/gene activity 可以作为额外模态加入。
- 模型用系统先验和 residual/PLS 约束 KO 方向，避免小样本下自由生成模型跑偏。
- 输出清晰可读的 heatmap、UMAP 和 AUC/ROC 曲线。

## 推荐输入

普通用户推荐直接输入 `.h5ad`：

```text
adata.X              cells x genes RNA matrix
adata.var_names      gene symbols
adata.obs[ko_col]    每个细胞的 KO / perturbation 标签
adata.obsm[...]      可选 protein/ADT/ATAC/gene activity 矩阵
```

也支持 CSV 原始表达矩阵：

```text
ko_target,STAT1,JAK2,IRF1,...
control,0.1,1.3,0.4,...
STAT1,0.0,1.2,0.2,...
```

`pathway/program score` 是软件内部自动生成的模型状态表示，不是用户必须提供的输入。

## 一键运行 h5ad 示例

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

## 输出

运行后会生成：

- `derived_state_scores.csv`: 软件自动派生出的 pathway/protein/ATAC 状态表
- `derived_state_manifest.csv`: 每个状态特征来自哪里
- `summary.csv`: 总体效果
- `metrics.csv`: 每个 KO、每个 feature 的分布距离
- `delta_table.csv`: 真实 KO 变化 vs 虚拟 KO 变化
- `virtual_cells.csv`: 虚拟 KO 单细胞状态
- `auc_summary.csv`: 强响应识别 AUC
- `report.md`: 自动报告

默认主图：

- `01_summary_dashboard.png`
- `02_true_vs_virtual_heatmap.png`
- `03_cell_state_umap.png`
- `04_auc_strong_response_roc.png`

## 当前 h5ad 示例结果

输出目录：

```text
results/software_interface_raw_papalexi
```

结果：

这个命令默认展示“虚拟敲除一个基因”的使用方式。评估多个 KO 时才需要写成逗号分隔，例如 `STAT1,JAK2`。

## 模型学过多少 KO？没学过的基因怎么预测？

当前原型接入过多个真实 perturbation 数据集，包括 Papalexi、Norman、HMPCITE、Liscovitch ATAC、Datlinger 和 Dixit。合并看，本地实验覆盖了百级别扰动基因，大约 129 个基因名。

但这个数字不等于“某一个模型已经一次性学过 129 个基因”。每次训练通常发生在一个具体 reference 数据集里，例如用 Papalexi 的部分单基因 KO 学方向，再预测留出的 KO；或用 Norman 的单基因 KO 预测双基因组合。

对于训练中没见过的基因，模型不会使用 one-hot 记忆，而是构建系统先验表示：

```text
target gene
-> Reactome/MSigDB pathway
-> TF-target network
-> PPI neighborhood
-> motif / chromVAR / ATAC prior
-> KO prior vector
-> predicted KO delta
```

因此 unseen gene prediction 本质上是基于生物网络先验的 zero-shot / few-shot 外推。它适合做候选筛选和方向判断；如果目标基因缺少 pathway、TF-target、PPI 或 motif 先验覆盖，软件应该把结果标记为低置信度。

## 普通 10X 单细胞数据

如果用户手里只是普通 10X scRNA-seq，没有 KO/perturbation 标签，可以先运行：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli score `
  --input-h5ad your_10x_data.h5ad `
  --prior-dir data\priors `
  --out-dir results\your_10x_state_scores
```

这一步会自动生成：

- `derived_state_scores.csv`
- `derived_state_manifest.csv`

注意：没有真实 KO 标签时，软件可以准备细胞状态并用于后续虚拟 KO 应用，但不能计算真实 KO 对照、AUC、UMAP 改善等评估指标。要评估“好不好”，必须有 perturb-seq / CRISPR / 药物扰动等真实 perturbation 标签。

## Reference model 应用到普通 10X

从 perturbation 数据训练 reference KO 模型：

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

把 reference model 应用到普通 h5ad 细胞：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli apply-reference `
  --reference-model results\reference_models\papalexi_rna_protein_reference.pkl `
  --input-h5ad your_10x_data.h5ad `
  --target-kos STAT1 `
  --out-dir results\your_10x_virtual_STAT1
```

示例输出：

```text
results/reference_apply_stat1_demo
```

注意：这是 prediction-only 模式。如果输入普通 10X 没有真实 KO 标签，不会报告真实准确率、AUC 或 distribution improvement。

## 当前 h5ad 多 KO 评估示例

如果要复现实验评估，可以一次评估多个真实 KO：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --ko-col ko_target `
  --target-kos STAT1,JAK2,IFNGR2,IRF1 `
  --prior-dir data\priors `
  --out-dir results\software_interface_raw_papalexi `
  --dataset-name "Papalexi ECCITE-seq" `
  --modality "raw RNA matrix + ADT protein obsm" `
  --representation "auto-derived pathway/protein scores" `
  --protein-obsm protein `
  --calibrate none `
  --max-pathways 24
```

示例结果：

- mean distribution improvement: `0.111`
- improved fraction: `0.759`
- mean direction cosine: `0.679`
- strong-response ROC-AUC: `0.802`

## 文档

- `docs/software_interface.md`: 用户输入、运行方式和输出解释
- `docs/reference_model_workflow_v2.md`: reference model v2、批量 KO、cell type 分层、prediction-only 10X/multiome 报告
- `docs/single_double_knockout_results.md`: 单基因敲除和双基因敲除结果整理
- `docs/current_scope_and_limitations.md`: 当前支持范围和限制说明
- `docs/virtual_knockout_method_background.md`: 虚拟敲除背景、主流方法和当前方法选择
- `docs/norman_double_interaction_model.md`: 双敲 interaction residual 优化结果
- `docs/two_day_experiment_summary.md`: 两天实验和调试总结
- `docs/pathway_magnitude_calibration.md`: RNA-only 幅度校准结果
- `docs/multi_dataset_virtual_ko_demo.md`: 多数据集虚拟 KO 示例

## 开发结构

- `vkx/preprocess.py`: 原始矩阵到 pathway/program/protein 状态表
- `vkx/core.py`: prior-constrained virtual KO 模型
- `vkx/visualization.py`: heatmap、UMAP、AUC、summary dashboard
- `vkx/cli.py`: 命令行入口

下一步重点是在这个接口上扩展和评估多基因 KO。
