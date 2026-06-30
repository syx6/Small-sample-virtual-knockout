# 2026-06-30 模型增强说明：reference、双敲、batch、三模态与 ATAC peak

本次更新把之前分散在实验脚本里的几个能力进一步接入主软件接口。核心目标是：用户不需要理解内部脚本，也能用统一命令完成 reference model 训练、普通 10X/多模态数据应用、双敲预测和可视化输出。

## 1. Interaction residual 已接入 reference model

以前 double-KO interaction residual 只能通过 `double-interaction` 单独评估。现在 `train-reference` 会在训练集中同时存在足够单敲和双敲标签时，自动训练一个双敲交互残差模型。

训练逻辑：

- 单敲均值效应先作为 additive baseline。
- 双敲真实效应减去 additive baseline，得到 interaction residual。
- residual 再由 KO gene、pathway/TF/PPI/motif prior 和 pair-wise prior features 预测。
- `apply-reference` 遇到双敲 `GENE1+GENE2` 时，会优先使用 interaction residual；如果训练数据不够，则回退到 prior-constrained PLS。

输出解释：

- `predicted_ko_delta.csv` 新增 `prediction_source`。
- `prediction_source = interaction_residual` 表示该双敲使用了训练好的交互残差模型。
- `prediction_source = prior_constrained_pls` 表示训练数据不足或目标不是双敲，使用稳定 baseline。

推荐训练命令：

```bash
python -m vkx.cli train-reference \
  --state-csv path/to/state_scores_with_ko.csv \
  --ko-col ko_target \
  --prior-dir data/priors \
  --output-model results/reference_model.pkl \
  --interaction-mode auto
```

强制要求必须训练 interaction residual：

```bash
python -m vkx.cli train-reference ... --interaction-mode on
```

如果数据里没有至少 3 个单敲和 3 个双敲标签，`--interaction-mode on` 会报错；`auto` 会自动回退并在 metadata 里写明原因。

## 2. Batch covariate 显式建模

新增 `--batch-col`。它适合 donor、sample、batch、library、replicate 等批次变量。

目前采用小样本更稳的 control-centered batch correction：

1. 在每个 batch 内计算 control cells 的平均状态。
2. 计算全局 control cells 的平均状态。
3. 用二者差值校正该 batch 的所有细胞。
4. KO effect 在校正后的 state space 中学习。

这样做的目的不是做复杂 batch integration，而是尽量防止模型把 donor/sample 差异误认为 KO 效应。

训练：

```bash
python -m vkx.cli train-reference \
  --state-csv path/to/state_scores_with_ko.csv \
  --ko-col ko_target \
  --batch-col donor \
  --output-model results/reference_model.pkl
```

应用：

```bash
python -m vkx.cli apply-reference \
  --reference-model results/reference_model.pkl \
  --state-csv path/to/ordinary_10x_state_scores.csv \
  --target-kos STAT1,STAT1+JAK2 \
  --batch-col donor \
  --out-dir results/reference_apply
```

新增输出：

- `batch_composition.csv`
- `06_batch_composition.png`
- reference metadata 中的 `batch_summary`

## 3. Hard-constrained uncertainty 输出

新增 `--uncertainty-method` 和 `--uncertainty-scale`。

当前默认策略不是自由 VAE / flow / diffusion，而是 hard-constrained residual band：

- KO 平均方向仍由 residual/PLS baseline 固定。
- 训练 KO 细胞围绕各自 KO 均值的残差宽度，用来估计单细胞不确定性。
- 输出区间只表示“在 hard constraint 附近可能波动多大”，不允许生成器自由改变 KO 方向。

应用命令：

```bash
python -m vkx.cli apply-reference \
  --reference-model results/reference_model.pkl \
  --state-csv path/to/ordinary_10x_state_scores.csv \
  --target-kos STAT1+JAK2 \
  --uncertainty-method hard-residual \
  --uncertainty-scale 0.25 \
  --out-dir results/reference_apply
```

输出：

- `uncertainty_intervals.csv`
- `07_uncertainty_intervals.png`

说明：CLI 也预留了 `vae`、`flow`、`diffusion` 入口，但当前它们仍被限制在 hard residual anchor 框架内。后续真正接入轻量 neural generator 时，仍会保持“baseline 定方向，generator 只学方向附近不确定性”的原则。

## 4. RNA + ADT + ATAC 输入方式

用户输入仍然可以是普通单细胞矩阵，不需要自己提前算 pathway。

推荐 h5ad 结构：

```text
adata.X                 RNA gene matrix
adata.obs["ko_target"]  perturbation / KO label，可选
adata.obs["batch"]      batch/sample/donor，可选
adata.obs["cell_type"]  cell type，可选
adata.obsm["protein"]   ADT / CITE-seq protein，可选
adata.obsm["atac"]      ATAC gene activity 或 peak-level features，可选
adata.obsm["chromvar"]  chromVAR motif activity，可选
adata.obsm["peak"]      raw peak count / peak accessibility，可选
```

带 KO 标签时可以训练或 benchmark：

```bash
python -m vkx.cli train-reference \
  --input-h5ad perturb_rna_adt_atac.h5ad \
  --ko-col ko_target \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv peak_annotation.csv \
  --max-extra-features-per-obsm 500 \
  --output-model results/trimodal_reference.pkl
```

无 KO 标签的普通 10X / DOGMA / TEA-seq 只能做 reference application：

```bash
python -m vkx.cli apply-reference \
  --reference-model results/trimodal_reference.pkl \
  --input-h5ad ordinary_multiome.h5ad \
  --target-kos STAT1,STAT1+JAK2 \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --out-dir results/prediction_only_multiome
```

注意：无 KO 标签时，软件不会报告真实准确率、AUC、R2 或 MAE；只能报告 predicted state shift、PCA/UMAP 类状态变化图、prior coverage、transfer confidence 和 uncertainty band。

## 5. 公开 benchmark 边界

目前公开数据可以分成四类：

| 数据类型 | perturbation label | RNA | ADT/protein | ATAC | 用途 |
|---|---:|---:|---:|---:|---|
| Multiome Perturb-seq / MultiPerturb-seq | yes | yes | no | yes | RNA+ATAC labelled benchmark |
| Perturb-CITE-seq / ECCITE-seq | yes | yes | yes | no | RNA+ADT labelled benchmark |
| Perturb-ATAC / scPerturb ATAC | yes | partial/no | no | yes | ATAC regulatory benchmark |
| DOGMA-seq / TEA-seq | usually no | yes | yes | yes | 三模态输入兼容和 reference application |

因此当前不能把 DOGMA/TEA-seq 当作准确性验证，也不能声称已经有公开 RNA+ADT+ATAC+CRISPR 同细胞完整 benchmark。更详细的公开数据清单见：

```text
docs/public_multimodal_perturbation_benchmark_registry_2026_06_30.md
```

软件也提供 registry 命令：

```bash
python -m vkx.cli benchmark-registry \
  --out-dir results/public_benchmark_registry
```

并且 `validate-benchmark` 会输出更具体的模式：

- `labelled_rna_atac_benchmark`
- `labelled_rna_adt_benchmark`
- `labelled_rna_adt_atac_benchmark`
- `trimodal_prediction_only`
- `prediction_only_or_incomplete`

只有 labelled benchmark 才能报告真实准确率、AUC、R2 和 MAE；`trimodal_prediction_only` 只能看预测状态变化、prior coverage、transfer confidence 和 uncertainty band。

## 6. ATAC raw peak count 与 peak annotation

新增 `extra-feature-selection=atac_peak`。

它比单纯 variance 更适合 raw peak/count：

- 结合 peak 方差；
- 结合 control vs KO 的 peak effect；
- 结合开放比例，避免只选全局开放或几乎全关闭的 peak；
- 结合外部 peak annotation 权重，例如 peak-gene linkage、motif-to-peak、marker peak 和 regulatory prior。

新增 `--extra-feature-metadata-csv` 后，可以直接接收外部 peak annotation 表。推荐列名：

```text
obsm_key
feature_index
feature_name 或 peak
target_gene
peak_gene_link_score
motif_to_peak_score
marker_score
regulatory_prior_score
```

这些列会写入 `derived_state_manifest.csv`。其中 `peak_gene_link_score`、`motif_to_peak_score`、`marker_score` 和 `regulatory_prior_score` 会参与 `atac_peak` feature selection。

如果用户还没有 annotation 表，可以先用 `annotate-peaks` 自动生成一个基础版本：

```bash
python -m vkx.cli build-peak-annotation \
  --input-h5ad perturb_multiome.h5ad \
  --obsm-key peak \
  --gtf gencode.annotation.gtf \
  --raw-motif-hits-csv raw_motif_hits.csv \
  --raw-marker-peaks-csv raw_marker_peaks.csv \
  --target-genes KDM6A,STAT1 \
  --out-csv results/peak_annotation.csv
```

如果想分步检查，也可以这样运行：

```bash
python -m vkx.cli make-gene-tss \
  --gtf gencode.annotation.gtf \
  --out-csv results/gene_tss.csv

python -m vkx.cli standardize-peak-scores \
  --input-csv raw_motif_hits.csv \
  --table-type motif \
  --out-csv results/motif_to_peak.csv

python -m vkx.cli standardize-peak-scores \
  --input-csv raw_marker_peaks.csv \
  --table-type marker \
  --out-csv results/marker_peaks.csv

python -m vkx.cli annotate-peaks \
  --input-h5ad perturb_multiome.h5ad \
  --obsm-key peak \
  --gene-tss-csv results/gene_tss.csv \
  --motif-hits-csv results/motif_to_peak.csv \
  --marker-peaks-csv results/marker_peaks.csv \
  --target-genes KDM6A,STAT1 \
  --out-csv results/peak_annotation.csv
```

其中 `gene_tss.csv` 至少包含 `gene,chrom,tss`，`motif_to_peak.csv` 和 `marker_peaks.csv` 至少包含 `peak` 或 `feature_name`。如果只有 peak 名称、没有 h5ad，也可以用 `--feature-names-csv peak_names.csv`。

`annotate-peaks` 和 `build-peak-annotation` 会额外输出：

- `peak_annotation.report.md`
- `peak_annotation.summary.png`

这张 summary 图会展示 top regulatory-prior peaks 以及 peak-gene、motif、marker、locus 四类证据的贡献。

推荐：

```bash
python -m vkx.cli run \
  --input-h5ad perturb_multiome.h5ad \
  --ko-col ko_target \
  --extra-obsm peak:peak,chromvar:tf \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv peak_annotation.csv \
  --shape-calibrate quantile \
  --target-kos KDM6A \
  --out-dir results/atac_peak_virtual_ko
```

如果 peak 是原始非负 sparse peak/count，建议配合：

- `--extra-feature-selection atac_peak`
- `--shape-calibrate quantile`

这样可以同时处理平均方向、开放比例和分位数形状。

## 7. 仍需继续补齐

- 真正公开 RNA+ADT+ATAC 且带 perturbation 标签的数据集仍未确认；找到后再升级为 full trimodal labelled benchmark。
- motif-to-peak annotation 和 peak-gene linkage 已经可以通过 `--extra-feature-metadata-csv` 输入并写入 feature metadata；下一步是开发自动生成 annotation 表的辅助脚本。
- 已新增 `make-gene-tss`、`standardize-peak-scores`、`annotate-peaks` 和一键式 `build-peak-annotation`，可以从 GTF、motif hits、marker peaks 自动生成基础 peak annotation 表与 QC 图；下一步可以继续接 Ensembl/GENCODE 自动下载和更具体的 motif scanner 输出格式转换。
- VAE / flow matching / diffusion 入口已保留，但当前仍采用 hard residual uncertainty band；真正 neural generator 需要在有更多同类型 perturbation 数据后再训练。
