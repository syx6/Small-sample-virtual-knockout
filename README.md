# Small-sample Multimodal Virtual Knockout (VKX)

VKX is an interpretable virtual knockout framework for small-sample single-cell and multimodal perturbation data. It predicts how cell states may change after single-gene or double-gene knockout by combining interpretable state scores, system-level regulatory priors, and a hard-constrained residual/PLS baseline.

中文发布版主文档：[`docs/VKX_publication_package_zh.md`](docs/VKX_publication_package_zh.md)

## What VKX Does

VKX takes ordinary single-cell or multimodal matrices as input:

- RNA expression matrix
- optional ADT / CITE-seq protein matrix
- optional ATAC gene activity, chromVAR motif activity, or peak-level accessibility
- optional perturbation labels for labelled benchmark
- optional batch, donor, and cell type metadata

It outputs:

- predicted virtual KO state changes
- virtual KO cells
- real vs virtual KO heatmaps
- ROC/AUC curves
- before/after UMAP or PCA views
- single-KO vs double-KO response maps
- ATAC peak locus tracks
- prediction-only reports for ordinary 10X data without KO labels

## Core Idea

VKX does not train a fully free generative model from scratch. Instead, it first learns a stable and interpretable KO direction, then generates virtual cells near that hard constraint:

```text
virtual KO state = control state + predicted KO delta + bounded residual
```

This design is intended for small-sample settings where large perturbation pretraining is not available.

## Model Schematic

![VKX model schematic](docs/assets/draft_figures/01_vkx_model_algorithm_schematic.png)

## Current Results

| Dataset / task | Modality | KO type | AUC | Direction | MAE |
|---|---|---:|---:|---:|---:|
| Papalexi ECCITE-seq | RNA + ADT | single KO | 0.878 | 0.879 | 0.424 |
| Norman Perturb-seq | RNA program | double KO | 1.000 | 0.960 | 0.127 |
| HMPCITE-seq GSE243244 | RNA + ADT + GDO-derived | double KO | 0.978 | 0.976 | 0.114 |
| scPerturb ATAC K562 | ATAC peak + prior | peak-level | 0.674 | 0.771 | 0.061 |

The Norman AUC is based on only five program-level features and should be interpreted together with heatmaps, MAE, and response maps.

## Key Figures

![Publication main figure](docs/assets/draft_figures/00_publication_main_figure.png)

![ROC AUC curves](docs/assets/draft_figures/02_auc_roc_curves.png)

![Real vs virtual heatmap](docs/assets/draft_figures/03_real_vs_virtual_method_heatmap.png)

![Before after UMAP](docs/assets/draft_figures/07_before_after_umap_panel.png)

## Quick Start

### Labelled perturbation data

```bash
python -m vkx.cli run \
  --input-h5ad labelled_perturbation.h5ad \
  --ko-col ko_target \
  --target-kos STAT1,STAT1+JAK2 \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --out-dir results/labelled_virtual_ko
```

### Train a reference model

```bash
python -m vkx.cli train-reference \
  --input-h5ad perturbation_reference.h5ad \
  --ko-col ko_target \
  --batch-col donor \
  --interaction-mode auto \
  --prior-dir data/priors \
  --output-model results/reference_models/vkx_reference.pkl
```

### Apply to ordinary 10X or unlabelled data

```bash
python -m vkx.cli apply-reference \
  --reference-model results/reference_models/vkx_reference.pkl \
  --input-h5ad ordinary_10x_or_multiome.h5ad \
  --target-kos STAT1,STAT1+JAK2 \
  --cell-type-col cell_type \
  --batch-col donor \
  --out-dir results/prediction_only_STAT1
```

For data without KO labels, VKX only produces a prediction-only report. It must not report real AUC, MAE, R2, or true accuracy inside that dataset.

## Publication Package

The current release-ready write-up is:

- [`docs/VKX_publication_package_zh.md`](docs/VKX_publication_package_zh.md): publication-style Chinese method and results package
- [`docs/vkx_method_results_draft_zh.md`](docs/vkx_method_results_draft_zh.md): longer development draft with historical results
- [`docs/vkx_method_results_draft_zh.html`](docs/vkx_method_results_draft_zh.html): local HTML reading version

## Current Boundaries

- VKX supports single KO, double KO, reference model application, batch/cell-type metadata, and multimodal extra-obsm inputs.
- The current neural generator is hard-constrained residual uncertainty modelling, not a fully free diffusion/VAE/flow model.
- A true public RNA+ADT+ATAC+genetic perturbation labelled benchmark still needs confirmation.
- scGen, CPA, GEARS, and CellOT need a fully reproducible same-split benchmark before making claims of superiority.

