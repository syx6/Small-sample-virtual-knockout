from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PUBLIC_BENCHMARKS = [
    {
        "dataset": "MultiPerturb-seq / Multiome Perturb-seq",
        "status": "labelled_benchmark",
        "rna": True,
        "adt": False,
        "atac": True,
        "perturbation_labels": True,
        "recommended_use": "RNA+ATAC labelled perturbation benchmark",
        "accession": "GSE278910",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278910",
        "notes": "Use for RNA pathway/program and ATAC/chromVAR/peak validation. It is not an RNA+ADT+ATAC benchmark.",
    },
    {
        "dataset": "Perturb-CITE-seq / ECCITE-seq",
        "status": "labelled_benchmark",
        "rna": True,
        "adt": True,
        "atac": False,
        "perturbation_labels": True,
        "recommended_use": "RNA+ADT labelled perturbation benchmark",
        "accession": "varies by study",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8011839/",
        "notes": "Use for RNA and protein/ADT validation. It cannot validate ATAC peak-level effects.",
    },
    {
        "dataset": "Perturb-ATAC / scPerturb ATAC",
        "status": "labelled_benchmark",
        "rna": False,
        "adt": False,
        "atac": True,
        "perturbation_labels": True,
        "recommended_use": "ATAC regulatory benchmark",
        "accession": "GSE116249 and related processed files",
        "url": "https://www.omicsdi.org/dataset/geo/GSE116249",
        "notes": "Use for ATAC gene activity, chromVAR, and peak-level regulatory prior tests.",
    },
    {
        "dataset": "DOGMA-seq / TEA-seq",
        "status": "prediction_only_reference_application",
        "rna": True,
        "adt": True,
        "atac": True,
        "perturbation_labels": False,
        "recommended_use": "Three-modality input compatibility and reference application",
        "accession": "e.g. GSE200417 registry",
        "url": "https://www.omicsdi.org/dataset/biostudies-literature/S-EPMC9219143",
        "notes": "Do not report true accuracy, AUC, R2, or MAE unless real perturbation labels are present.",
    },
    {
        "dataset": "Full RNA+ADT+ATAC+perturbation benchmark",
        "status": "not_confirmed_public",
        "rna": True,
        "adt": True,
        "atac": True,
        "perturbation_labels": True,
        "recommended_use": "Not yet enabled as a full trimodal labelled benchmark",
        "accession": "not confirmed",
        "url": "not confirmed",
        "notes": "As of the current registry, no public dataset has been confirmed to provide RNA, ADT/protein, ATAC, and real genetic perturbation labels in the same benchmark-ready cells. Upgrade only after manual verification of all four requirements.",
    },
]


def _control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("nt|control|non|safe|neg")


def _parse_extra_obsm_specs(value: str | None) -> list[tuple[str, str]]:
    if not value:
        return []
    specs = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            key, prefix = [part.strip() for part in item.split(":", 1)]
        else:
            key, prefix = item, item
        specs.append((key, prefix))
    return specs


def public_benchmark_registry(out_dir: str | Path | None = None) -> pd.DataFrame:
    table = pd.DataFrame(PUBLIC_BENCHMARKS)
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        table.to_csv(out / "public_benchmark_registry.csv", index=False)
        _write_public_benchmark_report(table, out)
    return table


def _write_public_benchmark_report(table: pd.DataFrame, out_dir: Path) -> None:
    text = f"""# Public Multimodal Perturbation Benchmark Registry

This registry separates true labelled benchmarks from prediction-only multiome application data.

## Practical Rule

- Labelled RNA+ATAC perturbation benchmark: use MultiPerturb-seq / Multiome Perturb-seq.
- Labelled RNA+ADT perturbation benchmark: use Perturb-CITE-seq / ECCITE-seq.
- ATAC regulatory benchmark: use Perturb-ATAC / scPerturb ATAC.
- RNA+ADT+ATAC without perturbation labels: use only for input compatibility and reference application.
- Full RNA+ADT+ATAC+genetic perturbation benchmark: not confirmed yet; do not claim full trimodal labelled validation until all four requirements are verified in the same cells.

## Registry

{table.to_string(index=False)}

## Interpretation Boundary

Do not report dataset-specific ROC-AUC, R2, MAE, or true-vs-virtual heatmaps on data without real perturbation labels.
For DOGMA-seq / TEA-seq style data, report predicted state shifts, transfer confidence, prior coverage, and uncertainty bands only.
"""
    (out_dir / "public_benchmark_registry_report.md").write_text(text, encoding="utf-8")


def validate_multimodal_benchmark(
    input_h5ad: str | Path,
    ko_col: str,
    extra_obsm: str | None,
    out_dir: str | Path,
) -> dict:
    import anndata as ad

    adata = ad.read_h5ad(input_h5ad)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    specs = _parse_extra_obsm_specs(extra_obsm)
    summary = _benchmark_summary(adata, ko_col, specs)
    labels = _label_summary(adata, ko_col)
    summary.to_csv(out / "benchmark_readiness.csv", index=False)
    labels.to_csv(out / "benchmark_label_counts.csv", index=False)
    _plot_benchmark_overview(summary, labels, out)
    _plot_modality_label_matrix(adata, ko_col, specs, out)
    _write_benchmark_report(summary, labels, ko_col, specs, out)
    return {"summary": summary, "labels": labels}


def _benchmark_summary(adata: ad.AnnData, ko_col: str, specs: list[tuple[str, str]]) -> pd.DataFrame:
    rows = [
        {"check": "cells", "status": "ok", "value": adata.n_obs, "interpretation": "number of cells"},
        {"check": "rna_features", "status": "ok" if adata.n_vars > 0 else "fail", "value": adata.n_vars, "interpretation": "RNA/gene activity features in adata.X"},
    ]
    if ko_col in adata.obs:
        labels = adata.obs[ko_col].astype(str)
        ctrl_mask = _control_mask(labels)
        ctrl = int(ctrl_mask.sum())
        ko = int((~ctrl_mask).sum())
        status = "ok" if ctrl >= 3 and ko >= 3 else "fail"
        rows.append({"check": "ko_column", "status": status, "value": ko_col, "interpretation": f"control cells={ctrl}; perturbed cells={ko}"})
        rows.append({"check": "ko_labels", "status": "ok" if labels.nunique() >= 2 else "fail", "value": int(labels.nunique()), "interpretation": "number of perturbation labels"})
    else:
        rows.append({"check": "ko_column", "status": "fail", "value": ko_col, "interpretation": "KO column not found in adata.obs"})
    for key, prefix in specs:
        if key in adata.obsm:
            value = adata.obsm[key]
            n_features = int(value.shape[1]) if hasattr(value, "shape") and len(value.shape) > 1 else 1
            rows.append({"check": f"obsm:{key}", "status": "ok", "value": n_features, "interpretation": f"extra modality mapped to prefix {prefix}"})
        else:
            rows.append({"check": f"obsm:{key}", "status": "missing", "value": 0, "interpretation": "requested extra modality not found"})
    has_atac = any(key in adata.obsm for key, prefix in specs if prefix.lower() in {"atac", "tf", "peak", "chromvar"})
    has_protein = any(key in adata.obsm for key, prefix in specs if prefix.lower() in {"protein", "adt"})
    has_labels = ko_col in adata.obs and labels.nunique() >= 2 if ko_col in adata.obs else False
    if has_labels and has_atac and has_protein:
        mode = "labelled_rna_adt_atac_benchmark"
        status = "ok"
    elif has_labels and has_atac:
        mode = "labelled_rna_atac_benchmark"
        status = "ok"
    elif has_labels and has_protein:
        mode = "labelled_rna_adt_benchmark"
        status = "ok"
    elif has_atac and has_protein:
        mode = "trimodal_prediction_only"
        status = "partial"
    elif has_labels:
        mode = "labelled_rna_only_benchmark"
        status = "partial"
    else:
        mode = "prediction_only_or_incomplete"
        status = "partial"
    rows.append({"check": "benchmark_mode", "status": status, "value": mode, "interpretation": _mode_interpretation(mode)})
    return pd.DataFrame(rows)


def _mode_interpretation(mode: str) -> str:
    interpretations = {
        "labelled_rna_adt_atac_benchmark": "KO labels plus RNA, ADT/protein, and ATAC were detected; this can support full trimodal validation if labels are real perturbations.",
        "labelled_rna_atac_benchmark": "KO labels plus RNA and ATAC were detected; this is a real RNA+ATAC perturbation benchmark, not an ADT benchmark.",
        "labelled_rna_adt_benchmark": "KO labels plus RNA and ADT/protein were detected; this is a real RNA+ADT benchmark, not an ATAC benchmark.",
        "trimodal_prediction_only": "RNA, ADT/protein, and ATAC were detected but KO labels were not; use prediction-only reports, not accuracy metrics.",
        "labelled_rna_only_benchmark": "KO labels and RNA were detected but requested non-RNA modalities are missing.",
        "prediction_only_or_incomplete": "Insufficient labels or requested modalities for benchmark accuracy.",
    }
    return interpretations.get(mode, "Unknown benchmark mode.")


def _label_summary(adata: ad.AnnData, ko_col: str) -> pd.DataFrame:
    if ko_col not in adata.obs:
        return pd.DataFrame(columns=["ko_target", "n_cells", "is_control"])
    labels = adata.obs[ko_col].astype(str)
    counts = labels.value_counts().reset_index()
    counts.columns = ["ko_target", "n_cells"]
    counts["is_control"] = _control_mask(counts["ko_target"])
    return counts


def _plot_benchmark_overview(summary: pd.DataFrame, labels: pd.DataFrame, out_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    from .visualization import setup_plot

    setup_plot()
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    status_color = {"ok": "#009E73", "partial": "#E69F00", "missing": "#999999", "fail": "#D55E00"}
    plot = summary.copy()
    plot["score"] = plot["status"].map({"ok": 1.0, "partial": 0.65, "missing": 0.35, "fail": 0.1}).fillna(0.4)
    axes[0].barh(plot["check"], plot["score"], color=[status_color.get(x, "#999999") for x in plot["status"]])
    axes[0].set_xlim(0, 1.05)
    axes[0].set_xlabel("Readiness")
    axes[0].set_title("Benchmark Readiness")
    for y, (_, row) in enumerate(plot.iterrows()):
        axes[0].text(0.03, y, row["status"], va="center", color="white", fontweight="bold")
    if labels.empty:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No KO labels detected", ha="center", va="center", fontsize=13)
    else:
        top = labels.head(15).copy()
        sns.barplot(data=top, x="n_cells", y="ko_target", hue="is_control", dodge=False, palette={True: "#4C78A8", False: "#F58518"}, ax=axes[1])
        axes[1].set_title("Top Perturbation Labels")
        axes[1].set_xlabel("Cells")
        axes[1].set_ylabel("")
        axes[1].legend(title="control")
    fig.savefig(out_dir / "benchmark_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_modality_label_matrix(adata: ad.AnnData, ko_col: str, specs: list[tuple[str, str]], out_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    from .visualization import setup_plot

    setup_plot()
    modalities = [{"modality": "RNA / X", "features": adata.n_vars}]
    for key, prefix in specs:
        if key in adata.obsm:
            matrix = adata.obsm[key]
            modalities.append({"modality": f"{prefix} ({key})", "features": matrix.shape[1] if len(matrix.shape) > 1 else 1})
    mod = pd.DataFrame(modalities)
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.45 * len(mod) + 2)), constrained_layout=True)
    sns.barplot(data=mod, x="features", y="modality", color="#4C78A8", ax=ax)
    ax.set_title("Detected Benchmark Modalities")
    ax.set_xlabel("Feature count")
    ax.set_ylabel("")
    if ko_col in adata.obs:
        ax.text(0.98, 0.05, f"KO labels: {adata.obs[ko_col].nunique()}", transform=ax.transAxes, ha="right", va="bottom", fontsize=11)
    fig.savefig(out_dir / "benchmark_modalities.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_benchmark_report(summary: pd.DataFrame, labels: pd.DataFrame, ko_col: str, specs: list[tuple[str, str]], out_dir: Path) -> None:
    status = summary.loc[summary["check"] == "benchmark_mode", "status"].iloc[0]
    mode = summary.loc[summary["check"] == "benchmark_mode", "value"].iloc[0]
    interpretation = summary.loc[summary["check"] == "benchmark_mode", "interpretation"].iloc[0]
    accuracy_allowed = "yes" if str(mode).startswith("labelled_") and status == "ok" else "no"
    text = f"""# Multimodal Perturbation Benchmark Readiness Report

## Verdict

Benchmark mode: `{status}`

Detailed mode: `{mode}`

Interpretation: {interpretation}

Can report true-vs-virtual accuracy on this dataset: `{accuracy_allowed}`

If this value is `no`, use prediction-only plots and do not report dataset-specific ROC-AUC, R2, MAE, or true-vs-virtual heatmaps.

## KO Label Column

`{ko_col}`

## Requested Extra Modalities

{specs}

## Summary

{summary.to_string(index=False)}

## Top Labels

{labels.head(20).to_string(index=False) if not labels.empty else 'No labels detected.'}

## Required Figures

- `benchmark_overview.png`
- `benchmark_modalities.png`

## Next Command

If this is a labelled benchmark, run `python -m vkx.cli run` with the same `--ko-col` and `--extra-obsm` settings to generate heatmap, UMAP and ROC/AUC figures.

If this is prediction-only, use `train-reference` on a labelled reference dataset and then `apply-reference` on this input.
"""
    (out_dir / "benchmark_readiness_report.md").write_text(text, encoding="utf-8")
