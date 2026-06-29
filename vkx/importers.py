from __future__ import annotations

from pathlib import Path
import gzip

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.io import mmread

from .visualization import setup_plot


def _read_table_auto(path: Path, header: None | int = None) -> pd.DataFrame:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return pd.read_csv(handle, sep="\t", header=header)
    return pd.read_csv(path, sep="\t", header=header)


def _find_existing(folder: Path, names: list[str]) -> Path:
    for name in names:
        candidate = folder / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find any of: {', '.join(names)} in {folder}")


def read_10x_mtx(folder: str | Path) -> ad.AnnData:
    folder = Path(folder)
    matrix_path = _find_existing(folder, ["matrix.mtx.gz", "matrix.mtx"])
    features_path = _find_existing(folder, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
    barcodes_path = _find_existing(folder, ["barcodes.tsv.gz", "barcodes.tsv"])
    matrix = mmread(matrix_path).tocsr().T
    features = _read_table_auto(features_path, header=None)
    barcodes = _read_table_auto(barcodes_path, header=None)
    gene_symbols = features.iloc[:, 1].astype(str) if features.shape[1] > 1 else features.iloc[:, 0].astype(str)
    feature_types = features.iloc[:, 2].astype(str) if features.shape[1] > 2 else pd.Series(["Gene Expression"] * len(features))
    adata = ad.AnnData(matrix)
    adata.obs_names = barcodes.iloc[:, 0].astype(str).to_list()
    adata.var_names = gene_symbols.to_list()
    adata.var["feature_id"] = features.iloc[:, 0].astype(str).to_numpy()
    adata.var["feature_type"] = feature_types.to_numpy()
    adata.var_names_make_unique()
    return _split_feature_types(adata)


def read_10x_h5(path: str | Path) -> ad.AnnData:
    path = Path(path)
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        shape = tuple(group["shape"][:])
        matrix = sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
        barcodes = [x.decode() if isinstance(x, bytes) else str(x) for x in group["barcodes"][:]]
        features = group["features"]
        names = [x.decode() if isinstance(x, bytes) else str(x) for x in features["name"][:]]
        ids = [x.decode() if isinstance(x, bytes) else str(x) for x in features["id"][:]]
        if "feature_type" in features:
            types = [x.decode() if isinstance(x, bytes) else str(x) for x in features["feature_type"][:]]
        else:
            types = ["Gene Expression"] * len(names)
    adata = ad.AnnData(matrix)
    adata.obs_names = barcodes
    adata.var_names = names
    adata.var["feature_id"] = ids
    adata.var["feature_type"] = types
    adata.var_names_make_unique()
    return _split_feature_types(adata)


def read_h5seurat_basic(path: str | Path, assay: str = "RNA") -> ad.AnnData:
    path = Path(path)
    with h5py.File(path, "r") as handle:
        assay_group = handle.get(f"assays/{assay}")
        if assay_group is None:
            raise ValueError(f"h5Seurat assay '{assay}' was not found. Convert with SeuratDisk to h5ad if this file uses a newer layout.")
        counts = assay_group.get("counts") or assay_group.get("data")
        if counts is None:
            raise ValueError("Could not find counts/data in the h5Seurat assay. Convert to h5ad with SeuratDisk for this layout.")
        if all(key in counts for key in ["data", "indices", "indptr"]):
            data = counts["data"][:]
            indices = counts["indices"][:]
            indptr = counts["indptr"][:]
            shape = tuple(counts["dims"][:] if "dims" in counts else counts["shape"][:])
            matrix = sparse.csc_matrix((data, indices, indptr), shape=shape).T.tocsr()
        else:
            matrix = np.asarray(counts).T
        features = assay_group.get("features")
        cells = handle.get("cell.names") or handle.get("cells")
        var_names = [x.decode() if isinstance(x, bytes) else str(x) for x in features[:]] if features is not None else [f"gene_{i+1}" for i in range(matrix.shape[1])]
        obs_names = [x.decode() if isinstance(x, bytes) else str(x) for x in cells[:]] if cells is not None else [f"cell_{i+1}" for i in range(matrix.shape[0])]
    adata = ad.AnnData(matrix)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.var_names_make_unique()
    return adata


def _split_feature_types(adata: ad.AnnData) -> ad.AnnData:
    if "feature_type" not in adata.var:
        return adata
    feature_type = adata.var["feature_type"].astype(str)
    gene_mask = feature_type.str.lower().eq("gene expression")
    if not gene_mask.any():
        return adata
    primary = adata[:, gene_mask.to_numpy()].copy()
    for label, prefix in [("Antibody Capture", "protein"), ("Peaks", "peak"), ("CRISPR Guide Capture", "guide")]:
        mask = feature_type.eq(label).to_numpy()
        if mask.any():
            matrix = adata[:, mask].X
            primary.obsm[prefix] = matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)
            primary.uns[f"{prefix}_names"] = adata.var_names[mask].astype(str).to_list()
    return primary


def import_single_cell_data(
    input_path: str | Path,
    input_format: str,
    out_h5ad: str | Path,
    ko_metadata_csv: str | Path | None = None,
    cell_id_col: str = "cell_id",
) -> ad.AnnData:
    input_path = Path(input_path)
    input_format = input_format.lower()
    if input_format == "h5ad":
        adata = ad.read_h5ad(input_path)
    elif input_format == "10x_mtx":
        adata = read_10x_mtx(input_path)
    elif input_format == "10x_h5":
        adata = read_10x_h5(input_path)
    elif input_format == "h5seurat":
        adata = read_h5seurat_basic(input_path)
    else:
        raise ValueError("input_format must be one of: h5ad, 10x_mtx, 10x_h5, h5seurat")
    if ko_metadata_csv:
        meta = pd.read_csv(ko_metadata_csv)
        if cell_id_col not in meta.columns:
            raise ValueError(f"Metadata file must contain '{cell_id_col}'.")
        meta = meta.set_index(cell_id_col)
        aligned = meta.reindex(adata.obs_names)
        for col in aligned.columns:
            adata.obs[col] = aligned[col].astype("object").to_numpy()
    out_h5ad = Path(out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    return adata


def write_import_outputs(adata: ad.AnnData, out_dir: str | Path, input_format: str) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = _input_summary(adata, input_format)
    summary.to_csv(out_dir / "input_summary.csv", index=False)
    _plot_input_overview(adata, summary, out_dir)
    _write_import_report(summary, out_dir)


def _input_summary(adata: ad.AnnData, input_format: str) -> pd.DataFrame:
    rows = [
        {"item": "input_format", "value": input_format},
        {"item": "cells", "value": int(adata.n_obs)},
        {"item": "rna_features", "value": int(adata.n_vars)},
    ]
    for key in adata.obsm.keys():
        value = adata.obsm[key]
        rows.append({"item": f"obsm:{key}_features", "value": int(value.shape[1]) if hasattr(value, "shape") and len(value.shape) > 1 else 1})
    for col in adata.obs.columns:
        n = adata.obs[col].nunique(dropna=True)
        if n <= 50:
            rows.append({"item": f"obs:{col}_labels", "value": int(n)})
    return pd.DataFrame(rows)


def _plot_input_overview(adata: ad.AnnData, summary: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    obs_counts = []
    for col in adata.obs.columns:
        n = adata.obs[col].nunique(dropna=True)
        if 1 < n <= 30:
            top = adata.obs[col].astype(str).value_counts().head(10)
            obs_counts.extend({"column": col, "label": idx, "cells": value} for idx, value in top.items())
    modality_rows = []
    modality_rows.append({"modality": "RNA genes", "features": adata.n_vars})
    for key in adata.obsm.keys():
        value = adata.obsm[key]
        modality_rows.append({"modality": f"obsm:{key}", "features": value.shape[1] if hasattr(value, "shape") and len(value.shape) > 1 else 1})
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)
    mod = pd.DataFrame(modality_rows)
    sns.barplot(data=mod, x="features", y="modality", color="#4C78A8", ax=axes[0])
    axes[0].set_title("Imported Modalities")
    axes[0].set_xlabel("Feature count")
    axes[0].set_ylabel("")
    if obs_counts:
        obs = pd.DataFrame(obs_counts)
        obs["name"] = obs["column"] + ": " + obs["label"]
        sns.barplot(data=obs.head(18), x="cells", y="name", color="#F58518", ax=axes[1])
        axes[1].set_title("Top Metadata Labels")
        axes[1].set_xlabel("Cells")
        axes[1].set_ylabel("")
    else:
        axes[1].axis("off")
        axes[1].text(0.5, 0.5, "No categorical metadata detected", ha="center", va="center", fontsize=13)
    fig.suptitle(f"Input overview: {adata.n_obs:,} cells", fontsize=15)
    fig.savefig(out_dir / "input_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _write_import_report(summary: pd.DataFrame, out_dir: Path) -> None:
    text = f"""# Input Import Report

The input file was converted into an AnnData `.h5ad` file for the virtual knockout workflow.

## Summary

{summary.to_string(index=False)}

## Main Figure

- `input_overview.png`: shows detected modalities and top metadata labels.

## Next Step

Use the generated h5ad as input to `score`, `run`, `train-reference`, or `apply-reference`.
"""
    (out_dir / "import_report.md").write_text(text, encoding="utf-8")
