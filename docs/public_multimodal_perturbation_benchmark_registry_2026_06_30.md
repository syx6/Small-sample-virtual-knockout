# Public multimodal perturbation benchmark registry, checked 2026-06-30

本文件用于明确：哪些公开数据可以作为真实准确性 benchmark，哪些只能作为输入兼容或 reference application。这样后续报告不会把无标签多组学数据误解释为验证结果。

## 当前结论

截至 2026-06-30，尚未确认存在一个公开、可直接作为 benchmark 的 **同细胞 RNA + ADT/protein + ATAC + genetic perturbation label** 数据集。

也就是说，当前不能宣称本项目已经完成 full trimodal labelled benchmark。软件现在明确区分：

- RNA+ATAC+guide labelled benchmark；
- RNA+ADT+guide labelled benchmark；
- ATAC perturbation benchmark；
- RNA+ADT+ATAC 无标签 reference application；
- full RNA+ADT+ATAC+perturbation benchmark：`not_confirmed_public`。

## 判定标准

一个数据集只有同时满足下面 4 点，才可以升级为 full trimodal labelled benchmark：

1. 同一细胞或同一可配对 cell barcode 有 RNA 表达矩阵；
2. 同一细胞或同一可配对 cell barcode 有 ADT/protein 矩阵；
3. 同一细胞或同一可配对 cell barcode 有 ATAC / peak / chromatin accessibility 矩阵；
4. 同一细胞或同一可配对 cell barcode 有真实 perturbation label，例如 sgRNA、guide、KO target 或明确 perturbation target。

如果只有 RNA+ATAC+guide，不能叫 RNA+ADT+ATAC benchmark。  
如果只有 RNA+ADT+guide，不能验证 ATAC peak-level prediction。  
如果有 RNA+ADT+ATAC 但没有 KO/guide label，只能做 prediction-only reference application。

## 当前公开数据分类

| 数据类型 | perturbation label | RNA | ADT/protein | ATAC | 当前用途 |
|---|---:|---:|---:|---:|---|
| MultiPerturb-seq / Multiome Perturb-seq | yes | yes | no | yes | RNA+ATAC labelled benchmark |
| Perturb-CITE-seq / ECCITE-seq | yes | yes | yes | no | RNA+ADT labelled benchmark |
| Perturb-ATAC / scPerturb ATAC | yes | partial/no | no | yes | ATAC-only 或 ATAC regulatory benchmark |
| DOGMA-seq / TEA-seq | usually no genetic KO label | yes | yes | yes | 三模态输入兼容与 reference application |
| Full RNA+ADT+ATAC+genetic perturbation benchmark | not confirmed | yes | yes | yes | 暂不启用为准确性 benchmark |

## 1. MultiPerturb-seq / Multiome Perturb-seq

用途：真正的 RNA+ATAC+CRISPR perturbation benchmark。

可以验证：

- RNA pathway/program delta；
- ATAC gene activity / peak delta；
- chromVAR / motif activity；
- perturbation label 下的真实 vs 虚拟 KO。

限制：

- 没有 ADT/protein，因此不是 full RNA+ADT+ATAC benchmark。

软件建议：

```bash
python -m vkx.cli validate-benchmark \
  --input-h5ad multiperturb_processed.h5ad \
  --ko-col ko_target \
  --extra-obsm atac:atac,chromvar:tf,peak:peak \
  --out-dir results/multiperturb_readiness
```

## 2. Perturb-CITE-seq / ECCITE-seq

用途：RNA+ADT perturbation benchmark。

可以验证：

- RNA pathway/program delta；
- ADT/protein delta；
- protein-level KO response。

限制：

- 没有 ATAC，因此不能验证 peak-level regulatory prediction。

## 3. Perturb-ATAC / scPerturb ATAC

用途：ATAC regulatory benchmark。

可以验证：

- ATAC gene activity；
- chromVAR/motif activity；
- peak-level regulatory prior；
- sparse peak shape calibration。

限制：

- 通常不能作为 full RNA+ADT+ATAC benchmark。

## 4. DOGMA-seq / TEA-seq

用途：三模态输入兼容测试与 reference application。

可以输出：

- pathway + protein + ATAC state score；
- virtual KO predicted state shift；
- input vs virtual PCA/UMAP；
- transfer confidence；
- prior coverage；
- uncertainty band。

不能输出：

- 真实准确率；
- 该数据内部 ROC-AUC；
- 该数据内部 R2/MAE；
- 真实 KO vs 虚拟 KO heatmap。

原因：普通 DOGMA/TEA-seq 数据通常没有 genetic KO / guide perturbation label。

## 5. Full Trimodal Labelled Benchmark 状态

当前状态：`not_confirmed_public`。

在找到并人工确认公开数据同时具有 RNA、ADT/protein、ATAC 和 perturbation labels 之前，软件不会把任何数据自动标记为 full trimodal labelled benchmark。

如果后续找到候选数据，应先检查：

```text
adata.X or layer["rna"]              RNA matrix
adata.obsm["protein"] or ADT matrix  ADT/protein
adata.obsm["peak"] / "atac"          ATAC/peak/chromatin
adata.obs["ko_target"] or guide col  perturbation target
```

确认后再运行：

```bash
python -m vkx.cli validate-benchmark \
  --input-h5ad candidate_full_trimodal.h5ad \
  --ko-col ko_target \
  --extra-obsm protein:protein,peak:peak,chromvar:tf \
  --out-dir results/full_trimodal_readiness
```

只有 `benchmark_mode = labelled_rna_adt_atac_benchmark` 时，才能报告 full trimodal labelled accuracy。

## 参考来源

- Multiome Perturb-seq / Nature Methods 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11291144/
- MultiPerturb-seq GEO GSE278910: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278910
- MultiPerturb-seq method repository: https://gitlab.com/sanjanalab/mps
- Perturb-ATAC GEO GSE116249: https://www.omicsdi.org/dataset/geo/GSE116249
- ECCITE-seq / pooled CRISPR + RNA + protein: https://pmc.ncbi.nlm.nih.gov/articles/PMC8011839/
- DOGMA-seq / GSE200417 registry: https://www.omicsdi.org/dataset/biostudies-literature/S-EPMC9219143
