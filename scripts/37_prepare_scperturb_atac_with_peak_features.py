from __future__ import annotations

from pathlib import Path
import gzip

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "scperturb_atac" / "Liscovitch-BrauerSanjana2021_K562_1"
INPUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity_chromvar.h5ad"
OUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity_chromvar_peaks.h5ad"
TARGET_GENES = {"KDM6A"}


def peak_label(row: pd.Series) -> str:
    return f"{row['chromosome']}:{int(row['start'])}-{int(row['end'])}|{row.get('nearestGene', 'NA')}|{row.get('peakType', 'peak')}"


def zscore_columns(values: np.ndarray) -> np.ndarray:
    values = np.log1p(values.astype(np.float32))
    mean = np.nanmean(values, axis=0)
    std = np.nanstd(values, axis=0)
    std = np.where(std < 1e-6, 1.0, std)
    return ((values - mean.reshape(1, -1)) / std.reshape(1, -1)).astype(np.float32)


def peak_coord(row: pd.Series) -> str:
    return f"{row['chromosome']}:{int(row['start'])}-{int(row['end'])}"


def markerpeak_ranked_indices(base: Path, peak_var: pd.DataFrame, target_genes: set[str], limit: int = 48) -> list[int]:
    marker_dir = base / "markerpeak_target"
    if not marker_dir.exists():
        return []
    marker_obs = pd.read_csv(marker_dir / "obs.csv", index_col=0)
    marker_var = pd.read_csv(marker_dir / "var.csv", index_col=0)
    with gzip.open(marker_dir / "counts.mtx.gz", "rb") as handle:
        marker_matrix = mmread(handle)
    if sparse.issparse(marker_matrix):
        marker_matrix = marker_matrix.tocsr()
    else:
        marker_matrix = sparse.csr_matrix(np.asarray(marker_matrix, dtype=np.float32))
    coord_to_peak_idx = {peak_coord(row): int(idx) for idx, row in peak_var.iterrows()}
    selected = []
    for target in target_genes:
        matches = [i for i, name in enumerate(marker_obs.index.astype(str)) if name.upper() == target.upper()]
        if not matches:
            continue
        row = np.asarray(marker_matrix[matches[0]].todense()).reshape(-1)
        order = np.argsort(np.abs(row))[::-1]
        for marker_idx in order:
            coord = str(marker_var.index[marker_idx])
            peak_idx = coord_to_peak_idx.get(coord)
            if peak_idx is not None and peak_idx not in selected:
                selected.append(peak_idx)
            if len(selected) >= limit:
                return selected
    return selected


def main() -> None:
    adata = ad.read_h5ad(INPUT)
    peak_obs = pd.read_csv(BASE / "peak_bc" / "obs.csv")
    peak_var = pd.read_csv(BASE / "peak_bc" / "var.csv")
    with gzip.open(BASE / "peak_bc" / "counts.mtx.gz", "rb") as handle:
        peak_matrix = mmread(handle)
    if sparse.issparse(peak_matrix):
        peak_matrix = peak_matrix.tocsr().astype(np.float32)
    else:
        peak_matrix = sparse.csr_matrix(np.asarray(peak_matrix, dtype=np.float32))

    peak_obs_names = peak_obs["cell_barcode"].astype(str).to_numpy()
    pos = pd.Series(np.arange(len(peak_obs_names)), index=peak_obs_names)
    keep_cells = pos.reindex(adata.obs_names.astype(str)).to_numpy()
    if np.isnan(keep_cells).any():
        missing = int(np.isnan(keep_cells).sum())
        raise ValueError(f"{missing} cells from the gene-activity h5ad were not found in peak_bc obs.csv.")
    peak_matrix = peak_matrix[keep_cells.astype(int)]

    target_mask = peak_var["nearestGene"].astype(str).str.upper().isin(TARGET_GENES)
    target_idx = list(np.flatnonzero(target_mask.to_numpy()))
    marker_idx = markerpeak_ranked_indices(BASE, peak_var, TARGET_GENES, limit=48)
    ranked = peak_var.assign(score=peak_var["ncells"].fillna(0).astype(float)).sort_values("score", ascending=False)
    global_idx = [int(i) for i in ranked.index[:80]]
    selected = []
    for idx in [*target_idx, *marker_idx, *global_idx]:
        if idx not in selected:
            selected.append(idx)
        if len(selected) >= 128:
            break

    selected_var = peak_var.iloc[selected].reset_index(drop=True)
    peak_values = peak_matrix[:, selected].toarray()
    adata.obsm["peak"] = zscore_columns(peak_values)
    adata.uns["peak_names"] = [peak_label(row) for _, row in selected_var.iterrows()]
    adata.uns["peak_metadata"] = selected_var.to_dict(orient="list")
    adata.uns["peak_source"] = "peak_bc selected peaks from scPerturb ATAC Liscovitch-Brauer/Sanjana 2021"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")
    selected_var.to_csv(OUT.with_name("liscovitch_k562_selected_peak_metadata.csv"), index=False)
    print(adata)
    print(f"Saved {OUT}")
    print(f"Selected peak features: {len(selected)}")
    print(f"  target-locus peaks: {len(target_idx)}")
    print(f"  markerpeak-target candidates: {len(marker_idx)}")
    print(selected_var[["chromosome", "start", "end", "nearestGene", "peakType", "ncells"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
