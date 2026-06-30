from __future__ import annotations

from pathlib import Path


def write_workflow_template(mode: str, out_dir: str | Path) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mode = mode.lower().replace("_", "-")
    if mode == "labelled-benchmark":
        text = LABELLED_BENCHMARK
    elif mode == "prediction-only":
        text = PREDICTION_ONLY
    elif mode == "atac-peak":
        text = ATAC_PEAK
    elif mode == "reference":
        text = REFERENCE
    else:
        text = ALL_WORKFLOWS
    path = out / "workflow_template.md"
    path.write_text(text, encoding="utf-8")
    return path


LABELLED_BENCHMARK = r"""# Workflow: labelled perturbation benchmark

Use this when the input has real KO / perturbation labels.

```bash
python -m vkx.cli validate-benchmark \
  --input-h5ad perturbation_data.h5ad \
  --ko-col ko_target \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --out-dir results/benchmark_readiness

python -m vkx.cli run \
  --input-h5ad perturbation_data.h5ad \
  --ko-col ko_target \
  --target-kos STAT1,STAT1+JAK2 \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --shape-calibrate quantile \
  --out-dir results/labelled_virtual_ko
```

Interpretation: AUC, R2, MAE, true-vs-virtual heatmaps and state movement plots are valid only if labels are real perturbations.
"""


REFERENCE = r"""# Workflow: train and apply a reference model

Use this when one labelled perturbation dataset trains the model and another ordinary dataset receives virtual KO.

```bash
python -m vkx.cli train-reference \
  --input-h5ad labelled_reference.h5ad \
  --ko-col ko_target \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --extra-feature-selection atac_peak \
  --interaction-mode auto \
  --batch-col batch \
  --output-model results/reference_model.pkl

python -m vkx.cli apply-reference \
  --reference-model results/reference_model.pkl \
  --input-h5ad ordinary_cells.h5ad \
  --target-kos STAT1,STAT1+JAK2 \
  --cell-type-col cell_type \
  --batch-col batch \
  --uncertainty-method hard-residual \
  --uncertainty-samples-per-ko 250 \
  --out-dir results/reference_application
```

Interpretation: The application dataset is prediction-only unless it also contains real KO labels for evaluation.
"""


PREDICTION_ONLY = r"""# Workflow: ordinary 10X / unlabeled multiome prediction-only application

Use this when the input has no perturbation labels.

```bash
python -m vkx.cli score \
  --input-h5ad ordinary_10x_or_multiome.h5ad \
  --prior-dir data/priors \
  --extra-obsm protein:protein,chromvar:tf,peak:peak \
  --out-dir results/state_scores_only

python -m vkx.cli apply-reference \
  --reference-model results/reference_model.pkl \
  --state-csv results/state_scores_only/derived_state_scores.csv \
  --target-kos STAT1,STAT1+JAK2 \
  --uncertainty-method hard-residual \
  --out-dir results/prediction_only_virtual_ko
```

Interpretation: do not report true accuracy, AUC, R2 or MAE on unlabeled data.
"""


ATAC_PEAK = r"""# Workflow: ATAC peak-level regulatory prior

Use this before `run` or `train-reference` when peak-level ATAC features are available.

```bash
python -m vkx.cli build-peak-annotation \
  --input-h5ad perturb_multiome.h5ad \
  --obsm-key peak \
  --gtf gencode.annotation.gtf \
  --raw-motif-hits-csv raw_motif_hits.csv \
  --raw-marker-peaks-csv raw_marker_peaks.csv \
  --target-genes KDM6A,STAT1 \
  --out-csv results/peak_annotation.csv

python -m vkx.cli run \
  --input-h5ad perturb_multiome.h5ad \
  --ko-col ko_target \
  --target-kos KDM6A \
  --extra-obsm peak:peak,chromvar:tf \
  --extra-feature-selection atac_peak \
  --extra-feature-metadata-csv results/peak_annotation.csv \
  --shape-calibrate quantile \
  --out-dir results/atac_peak_virtual_ko
```

Interpretation: peak-level distribution plots require sparse peak/count features; centered peak scores should not be interpreted as raw open/closed counts.
"""


ALL_WORKFLOWS = "\n\n".join([LABELLED_BENCHMARK, REFERENCE, PREDICTION_ONLY, ATAC_PEAK])
