from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .core import control_mask
from .preprocess import parse_extra_obsm_specs
from .visualization import setup_plot


def validate_multimodal_benchmark(
    input_h5ad: str | Path,
    ko_col: str,
    extra_obsm: str | None,
    out_dir: str | Path,
) -> dict:
    adata = ad.read_h5ad(input_h5ad)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    specs = parse_extra_obsm_specs(extra_obsm)
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
        ctrl = int(control_mask(labels).sum())
        ko = int((~control_mask(labels)).sum())
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
    benchmark_status = "ok" if ko_col in adata.obs and (has_atac or has_protein) else "partial"
    rows.append({"check": "benchmark_mode", "status": benchmark_status, "value": benchmark_status, "interpretation": "ok means labelled multimodal perturbation benchmark; partial may still run but is not full trimodal validation"})
    return pd.DataFrame(rows)


def _label_summary(adata: ad.AnnData, ko_col: str) -> pd.DataFrame:
    if ko_col not in adata.obs:
        return pd.DataFrame(columns=["ko_target", "n_cells", "is_control"])
    labels = adata.obs[ko_col].astype(str)
    counts = labels.value_counts().reset_index()
    counts.columns = ["ko_target", "n_cells"]
    counts["is_control"] = control_mask(counts["ko_target"])
    return counts


def _plot_benchmark_overview(summary: pd.DataFrame, labels: pd.DataFrame, out_dir: Path) -> None:
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
    text = f"""# Multimodal Perturbation Benchmark Readiness Report

## Verdict

Benchmark mode: `{status}`

`ok` means the file contains KO labels and at least one requested non-RNA modality. `partial` means the data may still run, but should not be described as a full multimodal perturbation benchmark.

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

If benchmark mode is `ok`, run `python -m vkx.cli run` with the same `--ko-col` and `--extra-obsm` settings to generate heatmap, UMAP and ROC/AUC figures.
"""
    (out_dir / "benchmark_readiness_report.md").write_text(text, encoding="utf-8")
