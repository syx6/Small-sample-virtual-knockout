from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse

from .importers import read_10x_h5, read_10x_mtx
from .visualization import setup_plot


def _load_matrix(path: str | Path, fmt: str) -> ad.AnnData:
    fmt = fmt.lower()
    if fmt == "h5ad":
        return ad.read_h5ad(path)
    if fmt == "10x_mtx":
        return read_10x_mtx(path)
    if fmt == "10x_h5":
        return read_10x_h5(path)
    raise ValueError("format must be h5ad, 10x_mtx, or 10x_h5")


def _select_atac_features(adata: ad.AnnData, max_features: int | None) -> tuple[np.ndarray, list[str]]:
    n_features = adata.n_vars
    if max_features is None or max_features <= 0 or n_features <= max_features:
        keep = np.arange(n_features)
    else:
        x = adata.X
        if sparse.issparse(x):
            nonzero = np.asarray((x > 0).mean(axis=0)).reshape(-1)
            mean = np.asarray(x.mean(axis=0)).reshape(-1)
            score = nonzero * np.log1p(mean)
        else:
            arr = np.asarray(x, dtype=float)
            score = np.nanvar(arr, axis=0) * (np.nanmean(arr > 0, axis=0) + 1e-6)
        keep = np.sort(np.argsort(np.nan_to_num(score, nan=-np.inf))[::-1][:max_features])
    names = [str(adata.var_names[i]) for i in keep]
    return keep, names


def assemble_multiome_benchmark(
    rna_input: str | Path,
    rna_format: str,
    atac_input: str | Path,
    atac_format: str,
    metadata_csv: str | Path,
    output_h5ad: str | Path,
    out_dir: str | Path,
    cell_id_col: str = "cell_id",
    ko_col: str = "ko_target",
    max_atac_features: int | None = 500,
) -> ad.AnnData:
    rna = _load_matrix(rna_input, rna_format)
    atac = _load_matrix(atac_input, atac_format)
    metadata = pd.read_csv(metadata_csv)
    if cell_id_col not in metadata.columns:
        raise ValueError(f"Metadata CSV must contain '{cell_id_col}'.")
    if ko_col not in metadata.columns:
        raise ValueError(f"Metadata CSV must contain '{ko_col}'.")
    shared = pd.Index(rna.obs_names.astype(str)).intersection(pd.Index(atac.obs_names.astype(str))).intersection(pd.Index(metadata[cell_id_col].astype(str)))
    if len(shared) < 10:
        raise ValueError("Fewer than 10 shared barcodes across RNA, ATAC and metadata. Check barcode suffixes and cell_id_col.")
    shared = shared.sort_values()
    rna = rna[shared, :].copy()
    atac = atac[shared, :].copy()
    metadata = metadata.set_index(cell_id_col).loc[shared]
    keep, atac_names = _select_atac_features(atac, max_atac_features)
    atac_matrix = atac[:, keep].X
    out = rna.copy()
    out.obs = pd.concat([out.obs, metadata], axis=1)
    out.obsm["atac"] = atac_matrix.toarray() if sparse.issparse(atac_matrix) else np.asarray(atac_matrix)
    out.uns["atac_names"] = atac_names
    out.uns["multiome_assembly"] = {
        "rna_input": str(rna_input),
        "atac_input": str(atac_input),
        "metadata_csv": str(metadata_csv),
        "max_atac_features": max_atac_features,
        "shared_cells": int(len(shared)),
    }
    output_h5ad = Path(output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(output_h5ad)
    write_multiome_assembly_outputs(out, out_dir, ko_col=ko_col)
    return out


def write_multiome_assembly_outputs(adata: ad.AnnData, out_dir: str | Path, ko_col: str) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(
        [
            {"item": "shared_cells", "value": adata.n_obs},
            {"item": "rna_features", "value": adata.n_vars},
            {"item": "atac_features_in_obsm", "value": adata.obsm["atac"].shape[1] if "atac" in adata.obsm else 0},
            {"item": "ko_labels", "value": adata.obs[ko_col].nunique() if ko_col in adata.obs else 0},
        ]
    )
    summary.to_csv(out_dir / "multiome_assembly_summary.csv", index=False)
    _plot_multiome_assembly(adata, ko_col, out_dir)
    _write_multiome_report(summary, ko_col, out_dir)


def _plot_multiome_assembly(adata: ad.AnnData, ko_col: str, out_dir: Path) -> None:
    setup_plot()
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    mod = pd.DataFrame(
        [
            {"modality": "RNA genes", "features": adata.n_vars},
            {"modality": "ATAC selected features", "features": adata.obsm["atac"].shape[1] if "atac" in adata.obsm else 0},
        ]
    )
    sns.barplot(data=mod, x="features", y="modality", color="#4C78A8", ax=axes[0])
    axes[0].set_title("Assembled Multiome Modalities")
    axes[0].set_xlabel("Features")
    axes[0].set_ylabel("")
    if ko_col in adata.obs:
        counts = adata.obs[ko_col].astype(str).value_counts().head(15).reset_index()
        counts.columns = ["ko_target", "cells"]
        sns.barplot(data=counts, x="cells", y="ko_target", color="#F58518", ax=axes[1])
        axes[1].set_title("Top KO Labels")
        axes[1].set_xlabel("Cells")
        axes[1].set_ylabel("")
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, f"KO column not found: {ko_col}", ha="center", va="center", fontsize=12)
    fig.suptitle(f"Assembled RNA+ATAC benchmark: {adata.n_obs:,} shared cells", fontsize=15)
    fig.savefig(out_dir / "multiome_assembly_overview.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_multiome_report(summary: pd.DataFrame, ko_col: str, out_dir: Path) -> None:
    text = f"""# Multiome Benchmark Assembly Report

This command assembled RNA and ATAC matrices by shared cell barcodes and attached perturbation metadata.

## Summary

{summary.to_string(index=False)}

## KO Column

`{ko_col}`

## Main Figure

- `multiome_assembly_overview.png`

## Next Step

Run `validate-benchmark` on the generated h5ad. If readiness is `ok`, run the virtual KO benchmark with `--extra-obsm atac:atac`.
"""
    (out_dir / "multiome_assembly_report.md").write_text(text, encoding="utf-8")
