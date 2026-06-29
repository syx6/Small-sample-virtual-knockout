from __future__ import annotations

from pathlib import Path
import gzip
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from vkx.core import control_mask, parse_gmt
BASE = ROOT / "data" / "scperturb_atac" / "Liscovitch-BrauerSanjana2021_K562_1"
INPUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity_chromvar.h5ad"
OUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity_chromvar_peaks.h5ad"
TARGET_GENES = {"KDM6A"}
MAX_PEAKS = 160


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


def markerpeak_scores(base: Path, peak_var: pd.DataFrame, target_genes: set[str]) -> pd.Series:
    marker_dir = base / "markerpeak_target"
    scores = pd.Series(0.0, index=peak_var.index)
    if not marker_dir.exists():
        return scores
    marker_obs = pd.read_csv(marker_dir / "obs.csv", index_col=0)
    marker_var = pd.read_csv(marker_dir / "var.csv", index_col=0)
    with gzip.open(marker_dir / "counts.mtx.gz", "rb") as handle:
        marker_matrix = mmread(handle)
    marker_matrix = marker_matrix.tocsr() if sparse.issparse(marker_matrix) else sparse.csr_matrix(np.asarray(marker_matrix, dtype=np.float32))
    coord_to_peak_idx = {peak_coord(row): int(idx) for idx, row in peak_var.iterrows()}
    for target in target_genes:
        matches = [i for i, name in enumerate(marker_obs.index.astype(str)) if name.upper() == target.upper()]
        if not matches:
            continue
        row = np.asarray(marker_matrix[matches[0]].todense()).reshape(-1)
        abs_row = np.abs(row)
        if abs_row.max() > 0:
            abs_row = abs_row / abs_row.max()
        for marker_idx, value in enumerate(abs_row):
            peak_idx = coord_to_peak_idx.get(str(marker_var.index[marker_idx]))
            if peak_idx is not None:
                scores.loc[peak_idx] = max(float(scores.loc[peak_idx]), float(value))
    return scores


def scaled(values: pd.Series) -> pd.Series:
    arr = values.astype(float).replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy()
    if np.nanmax(arr) <= np.nanmin(arr):
        return pd.Series(np.zeros_like(arr, dtype=float), index=values.index)
    return pd.Series((arr - np.nanmin(arr)) / (np.nanmax(arr) - np.nanmin(arr)), index=values.index)


def locus_scores(peak_var: pd.DataFrame, target_genes: set[str], network_genes: set[str]) -> pd.Series:
    gene = peak_var["nearestGene"].astype(str).str.upper()
    peak_type = peak_var["peakType"].astype(str).str.upper()
    dist = peak_var.get("distToTSS", peak_var.get("distToGeneStart", pd.Series(np.inf, index=peak_var.index))).abs().astype(float)
    target = gene.isin(target_genes)
    network = gene.isin(network_genes)
    type_score = pd.Series(0.0, index=peak_var.index)
    type_score.loc[peak_type.str.contains("PROMOTER", na=False)] = 1.0
    type_score.loc[peak_type.str.contains("EXONIC|INTRONIC", na=False)] = 0.75
    type_score.loc[peak_type.str.contains("DISTAL", na=False)] = 0.45
    distance_score = 1.0 / (1.0 + dist / 50_000.0)
    return (target.astype(float) * (2.2 + type_score + distance_score) + network.astype(float) * (0.6 + 0.4 * type_score)).astype(float)


def target_network_genes(prior_dir: Path, target_genes: set[str], limit: int = 500) -> set[str]:
    network: set[str] = set()
    for path in sorted(prior_dir.glob("*.gmt")):
        if path.stem not in {"tf_target", "ppi_hub", "reactome", "hallmark"}:
            continue
        for term, genes in parse_gmt(path, include_term_gene=path.stem in {"tf_target", "ppi_hub"}):
            if genes & target_genes:
                network.update(genes)
            if len(network) >= limit:
                break
        if len(network) >= limit:
            break
    return {gene.upper() for gene in network}


def motif_tf_scores(peak_var: pd.DataFrame, target_genes: set[str], network_genes: set[str]) -> pd.Series:
    gene = peak_var["nearestGene"].astype(str).str.upper()
    score = pd.Series(0.0, index=peak_var.index)
    score.loc[gene.isin(target_genes)] = 1.0
    score.loc[gene.isin(network_genes)] = np.maximum(score.loc[gene.isin(network_genes)], 0.45)
    return score


def ko_effect_scores(peak_matrix: sparse.csr_matrix, labels: pd.Series, target_genes: set[str]) -> pd.Series:
    labels = labels.astype(str).reset_index(drop=True)
    ctrl = control_mask(labels).to_numpy()
    target_mask = labels.str.upper().isin(target_genes).to_numpy()
    if ctrl.sum() < 3 or target_mask.sum() < 3:
        return pd.Series(np.zeros(peak_matrix.shape[1], dtype=float))
    ctrl_mean = np.asarray(peak_matrix[ctrl].mean(axis=0)).reshape(-1)
    ko_mean = np.asarray(peak_matrix[target_mask].mean(axis=0)).reshape(-1)
    effect = np.abs(np.log1p(ko_mean) - np.log1p(ctrl_mean))
    return pd.Series(effect)


def peak_selection_table(peak_var: pd.DataFrame, peak_matrix: sparse.csr_matrix, labels: pd.Series, target_genes: set[str]) -> pd.DataFrame:
    network = target_network_genes(ROOT / "data" / "priors", {gene.upper() for gene in target_genes})
    table = peak_var.copy()
    table["peak_coord"] = [peak_coord(row) for _, row in table.iterrows()]
    table["locus_score"] = locus_scores(table, {gene.upper() for gene in target_genes}, network)
    table["marker_score"] = markerpeak_scores(BASE, table, target_genes)
    table["ko_effect_score"] = ko_effect_scores(peak_matrix, labels, {gene.upper() for gene in target_genes}).to_numpy()
    table["accessibility_score"] = scaled(table["ncells"].fillna(0).astype(float))
    table["motif_tf_score"] = motif_tf_scores(table, {gene.upper() for gene in target_genes}, network)
    table["total_score"] = (
        0.30 * scaled(table["locus_score"])
        + 0.23 * scaled(table["ko_effect_score"])
        + 0.20 * scaled(table["marker_score"])
        + 0.15 * scaled(table["motif_tf_score"])
        + 0.12 * scaled(table["accessibility_score"])
    )
    table["selection_reason"] = "global_accessible"
    table.loc[table["locus_score"] > 0, "selection_reason"] = "target_or_network_locus"
    table.loc[table["marker_score"] > 0, "selection_reason"] = table.loc[table["marker_score"] > 0, "selection_reason"] + "+marker_peak"
    table.loc[table["ko_effect_score"] >= table["ko_effect_score"].quantile(0.98), "selection_reason"] = table.loc[
        table["ko_effect_score"] >= table["ko_effect_score"].quantile(0.98), "selection_reason"
    ] + "+ko_effect"
    return table


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

    labels = adata.obs["ko_target"].astype(str) if "ko_target" in adata.obs else adata.obs["perturbation"].astype(str)
    scoring = peak_selection_table(peak_var, peak_matrix, labels, TARGET_GENES)
    gene_upper = scoring["nearestGene"].astype(str).str.upper()
    target_idx = list(scoring.index[gene_upper.isin({gene.upper() for gene in TARGET_GENES})])
    target_idx = list(scoring.loc[target_idx].sort_values("total_score", ascending=False).head(64).index)
    marker_idx = list(scoring.sort_values("marker_score", ascending=False).head(64).index)
    ko_idx = list(scoring.sort_values("ko_effect_score", ascending=False).head(48).index)
    ranked_idx = list(scoring.sort_values("total_score", ascending=False).index)
    selected = []
    for idx in [*target_idx, *marker_idx, *ko_idx, *ranked_idx]:
        idx = int(idx)
        if idx not in selected:
            selected.append(idx)
        if len(selected) >= MAX_PEAKS:
            break

    selected_var = scoring.iloc[selected].reset_index(drop=True)
    peak_values = peak_matrix[:, selected].toarray()
    adata.obsm["peak"] = zscore_columns(peak_values)
    adata.uns["peak_names"] = [peak_label(row) for _, row in selected_var.iterrows()]
    adata.uns["peak_metadata"] = selected_var.to_dict(orient="list")
    adata.uns["peak_source"] = "peak_bc selected peaks from scPerturb ATAC Liscovitch-Brauer/Sanjana 2021"

    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")
    selected_var.to_csv(OUT.with_name("liscovitch_k562_selected_peak_metadata.csv"), index=False)
    scoring.sort_values("total_score", ascending=False).head(1000).to_csv(
        OUT.with_name("liscovitch_k562_peak_regulatory_prior_scores_top1000.csv"),
        index=False,
    )
    print(adata)
    print(f"Saved {OUT}")
    print(f"Selected peak features: {len(selected)}")
    print(f"  target/network locus peaks: {len(target_idx)}")
    print(f"  markerpeak-target candidates: {len(marker_idx)}")
    print(f"  KO-effect candidates: {len(ko_idx)}")
    print(selected_var[["chromosome", "start", "end", "nearestGene", "peakType", "selection_reason", "total_score"]].head(12).to_string(index=False))


if __name__ == "__main__":
    main()
