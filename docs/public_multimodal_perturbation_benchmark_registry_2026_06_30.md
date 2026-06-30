# Public multimodal perturbation benchmark registry, checked 2026-06-30

本文件用于明确：哪些公开数据可以作为真实准确性 benchmark，哪些只能作为输入兼容或 reference application。这样后续报告不会把无标签多组学数据误解释为验证结果。

## 当前结论

截至 2026-06-30，我能确认的公开数据类型如下：

| 数据类型 | 是否有 perturbation label | 是否有 RNA | 是否有 ADT/protein | 是否有 ATAC | 可用于本项目哪类实验 |
|---|---:|---:|---:|---:|---|
| Multiome Perturb-seq / MultiPerturb-seq | yes | yes | no | yes | RNA+ATAC perturbation benchmark |
| Perturb-ATAC / Spear-ATAC | yes | no 或不完整 | no | yes | ATAC-only 或 ATAC regulatory benchmark |
| Perturb-CITE-seq / ECCITE-seq | yes | yes | yes | no | RNA+ADT perturbation benchmark |
| DOGMA-seq / TEA-seq | usually no KO label | yes | yes | yes | 三模态输入兼容与 reference application，不作为准确性 benchmark |
| RNA+ADT+ATAC+CRISPR 同细胞公开 benchmark | not confirmed | yes | yes | yes | 暂不宣称已有准确性 benchmark |

## 可优先接入的数据

### 1. MultiPerturb-seq / Multiome Perturb-seq

用途：真正的 RNA+ATAC+CRISPR perturbation benchmark。

可验证：

- RNA pathway/program delta；
- ATAC gene activity / peak delta；
- chromVAR / motif activity；
- perturbation label 下的真实 vs 虚拟 KO。

公开入口：

- GEO: GSE278910, MultiPerturb-seq data consists of single-cell ATAC-seq, RNA-seq, and CRISPR guide RNA capture.
- Paper: Multiome Perturb-seq / MultiPerturb-seq extends single-cell CRISPR screens to simultaneously measure gene expression and chromatin accessibility.

软件建议：

```bash
python -m vkx.cli validate-benchmark \
  --input-h5ad multiperturb_processed.h5ad \
  --ko-col ko_target \
  --extra-obsm atac:atac,chromvar:tf,peak:peak \
  --out-dir results/multiperturb_readiness
```

```bash
python -m vkx.cli train-reference \
  --input-h5ad multiperturb_processed.h5ad \
  --ko-col ko_target \
  --extra-obsm atac:atac,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv peak_annotation.csv \
  --interaction-mode auto \
  --batch-col batch \
  --output-model results/multiperturb_reference.pkl
```

### 2. Perturb-CITE-seq / ECCITE-seq

用途：RNA+ADT perturbation benchmark。

可验证：

- RNA pathway/program delta；
- ADT/protein delta；
- protein-level KO response。

限制：

- 没有 ATAC，因此不能验证 peak-level regulatory prediction。

### 3. DOGMA-seq / TEA-seq

用途：三模态输入兼容测试与 reference application。

可输出：

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

原因：普通 DOGMA/TEA-seq 数据通常没有 KO perturbation label。

## 下一步接入策略

1. 优先下载并整理 MultiPerturb-seq / Multiome Perturb-seq，作为 RNA+ATAC labelled benchmark。
2. 保留 Perturb-CITE/ECCITE 作为 RNA+ADT benchmark。
3. DOGMA/TEA-seq 只作为三模态输入兼容示例。
4. 如果后续找到真正 RNA+ADT+ATAC+CRISPR 同细胞公开数据，再升级为 full trimodal labelled benchmark。

## 参考来源

- Multiome Perturb-seq / Nature Methods 2024: https://pmc.ncbi.nlm.nih.gov/articles/PMC11291144/
- MultiPerturb-seq GEO GSE278910: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278910
- MultiPerturb-seq method repository: https://gitlab.com/sanjanalab/mps
- Perturb-ATAC GEO GSE116249: https://www.omicsdi.org/dataset/geo/GSE116249
- ECCITE-seq / pooled CRISPR + RNA + protein: https://pmc.ncbi.nlm.nih.gov/articles/PMC8011839/
- DOGMA-seq benchmark / GSE200417 registry: https://www.omicsdi.org/dataset/biostudies-literature/S-EPMC9219143
