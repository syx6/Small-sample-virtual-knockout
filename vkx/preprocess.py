from __future__ import annotations

from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .core import GENE_RE, control_mask, parse_gmt


def _safe_label(text: str, prefix: str) -> str:
    label = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return f"{prefix}_{label[:52]}"


def _matrix_to_float(x):
    if sparse.issparse(x):
        return x.astype(np.float32)
    return np.asarray(x, dtype=np.float32)


def select_pathway_terms(
    adata: ad.AnnData,
    prior_dir: str | Path,
    max_terms: int = 40,
    min_genes: int = 8,
    max_genes: int = 500,
) -> list[tuple[str, list[str]]]:
    var_genes = {str(gene).upper() for gene in adata.var_names}
    candidates = []
    for path in sorted(Path(prior_dir).glob("*.gmt")):
        if path.stem == "ppi_hub":
            continue
        for term, genes in parse_gmt(path):
            present = sorted(genes & var_genes)
            if len(present) < min_genes or len(present) > max_genes:
                continue
            useful_bonus = int(
                any(
                    word in term.upper()
                    for word in [
                        "INTERFERON",
                        "JAK",
                        "STAT",
                        "MAPK",
                        "TGFB",
                        "APOPTOSIS",
                        "PROLIFERATION",
                        "CELL_CYCLE",
                        "MYC",
                        "E2F",
                        "IMMUNE",
                        "CYTOKINE",
                    ]
                )
            )
            candidates.append(((useful_bonus, len(present)), f"{path.stem}:{term}", present))
    candidates.sort(reverse=True, key=lambda item: item[0])
    selected, seen = [], set()
    for _, term, genes in candidates:
        label = _safe_label(term.split(":", 1)[-1], "pathway")
        if label in seen:
            continue
        selected.append((label, genes))
        seen.add(label)
        if len(selected) >= max_terms:
            break
    return selected


def compute_pathway_scores(adata: ad.AnnData, terms: list[tuple[str, list[str]]]) -> pd.DataFrame:
    var_lookup = {str(gene).upper(): i for i, gene in enumerate(adata.var_names)}
    x = _matrix_to_float(adata.X)
    scores = {}
    for label, genes in terms:
        idx = [var_lookup[gene] for gene in genes if gene in var_lookup]
        if not idx:
            continue
        values = np.asarray(x[:, idx].mean(axis=1)).reshape(-1)
        scores[label] = values
    frame = pd.DataFrame(scores, index=adata.obs_names)
    for col in frame.columns:
        std = frame[col].std()
        if std > 1e-9:
            frame[col] = (frame[col] - frame[col].mean()) / std
        else:
            frame[col] = 0.0
    return frame


def compute_pathway_scores_from_saved_terms(adata: ad.AnnData, terms: list[dict]) -> pd.DataFrame:
    normalized = [(term["state_feature"], [gene.upper() for gene in term["genes"]]) for term in terms]
    return compute_pathway_scores(adata, normalized)


def parse_extra_obsm_specs(value: str | None) -> list[tuple[str, str]]:
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
        if not key or not prefix:
            raise ValueError("Extra modality specs must look like obsm_key:prefix, e.g. protein:protein,atac:atac.")
        specs.append((key, prefix))
    return specs


def _obsm_names(adata: ad.AnnData, obsm_key: str) -> list[str]:
    for key in [f"{obsm_key}_names", "protein_names", "adt_names", "atac_names", "gene_activity_names"]:
        if key in adata.uns:
            return [str(x) for x in adata.uns[key]]
    matrix = adata.obsm[obsm_key]
    return [f"{obsm_key}_{i + 1}" for i in range(matrix.shape[1])]


def _scaled_rank(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    order = np.argsort(np.nan_to_num(values, nan=-np.inf))
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.linspace(0.0, 1.0, len(order), endpoint=True)
    return ranks


def _ko_effect_scores(values: np.ndarray, labels: pd.Series | None) -> np.ndarray | None:
    if labels is None:
        return None
    labels = labels.astype(str).reset_index(drop=True)
    ctrl_mask = control_mask(labels).to_numpy()
    if ctrl_mask.sum() < 3 or (~ctrl_mask).sum() < 3:
        return None
    ctrl = values[ctrl_mask]
    ctrl_mean = np.nanmean(ctrl, axis=0)
    ctrl_var = np.nanvar(ctrl, axis=0)
    scores = []
    for _, idx in labels.loc[~ctrl_mask].groupby(labels.loc[~ctrl_mask], observed=True).groups.items():
        group = values[np.asarray(list(idx), dtype=int)]
        if len(group) < 3:
            continue
        pooled = np.sqrt(ctrl_var + np.nanvar(group, axis=0) + 1e-6)
        scores.append(np.abs(np.nanmean(group, axis=0) - ctrl_mean) / pooled)
    if not scores:
        return None
    return np.nanmax(np.vstack(scores), axis=0)


def select_obsm_feature_indices(
    values: np.ndarray,
    max_features: int | None,
    selection: str = "variance",
    labels: pd.Series | None = None,
) -> np.ndarray:
    keep = np.arange(values.shape[1])
    if max_features is None or max_features <= 0 or values.shape[1] <= max_features:
        return keep
    variance = np.nanvar(values, axis=0)
    if selection == "ko_effect":
        effect = _ko_effect_scores(values, labels)
        score = effect if effect is not None else variance
    elif selection == "hybrid":
        effect = _ko_effect_scores(values, labels)
        if effect is None:
            score = variance
        else:
            score = 0.45 * _scaled_rank(variance) + 0.55 * _scaled_rank(effect)
    elif selection == "atac_peak":
        effect = _ko_effect_scores(values, labels)
        open_fraction = np.nanmean(values > 0, axis=0)
        sparse_window = 1.0 - np.abs(open_fraction - 0.20) / 0.20
        sparse_window = np.clip(sparse_window, 0.0, 1.0)
        if effect is None:
            score = 0.55 * _scaled_rank(variance) + 0.45 * _scaled_rank(sparse_window)
        else:
            score = 0.35 * _scaled_rank(variance) + 0.40 * _scaled_rank(effect) + 0.25 * _scaled_rank(sparse_window)
    else:
        score = variance
    return np.sort(np.argsort(np.nan_to_num(score, nan=-np.inf))[::-1][:max_features])


def append_obsm_scores(
    frame: pd.DataFrame,
    adata: ad.AnnData,
    obsm_key: str,
    prefix: str,
    max_features: int | None = None,
    feature_selection: str = "variance",
    labels: pd.Series | None = None,
) -> pd.DataFrame:
    if obsm_key not in adata.obsm:
        raise ValueError(f"obsm key '{obsm_key}' was not found. Available keys: {list(adata.obsm.keys())}")
    values = np.asarray(adata.obsm[obsm_key], dtype=float)
    names = _obsm_names(adata, obsm_key)
    keep = select_obsm_feature_indices(values, max_features, selection=feature_selection, labels=labels)
    cols = []
    seen = set(frame.columns)
    for i in keep:
        raw = names[i] if i < len(names) else f"{i + 1}"
        col = _safe_label(raw, prefix)
        suffix = 2
        base = col
        while col in seen:
            col = f"{base}_{suffix}"
            suffix += 1
        seen.add(col)
        cols.append(col)
    extra = pd.DataFrame(values[:, keep], columns=cols, index=frame.index)
    return pd.concat([frame, extra], axis=1).copy()


def append_extra_obsm_scores(
    frame: pd.DataFrame,
    adata: ad.AnnData,
    extra_obsm: list[tuple[str, str]],
    max_features_per_obsm: int | None = None,
    feature_selection: str = "variance",
    labels: pd.Series | None = None,
) -> pd.DataFrame:
    seen = set()
    for obsm_key, prefix in extra_obsm:
        if (obsm_key, prefix) in seen:
            continue
        seen.add((obsm_key, prefix))
        frame = append_obsm_scores(
            frame,
            adata,
            obsm_key,
            prefix,
            max_features=max_features_per_obsm,
            feature_selection=feature_selection,
            labels=labels,
        )
    return frame


def extra_obsm_manifest(frame: pd.DataFrame, extra_obsm: list[tuple[str, str]]) -> pd.DataFrame:
    rows = []
    seen = set()
    for obsm_key, prefix in extra_obsm:
        if (obsm_key, prefix) in seen:
            continue
        seen.add((obsm_key, prefix))
        cols = [col for col in frame.columns if col.startswith(f"{prefix}_")]
        rows.extend({"state_feature": col, "source": f"obsm:{obsm_key}", "n_genes_used": np.nan} for col in cols)
    return pd.DataFrame(rows)


def h5ad_to_state_table(
    input_h5ad: str | Path,
    ko_col: str,
    prior_dir: str | Path,
    max_pathways: int = 40,
    protein_obsm: str | None = None,
    protein_prefix: str = "protein",
    extra_obsm: list[tuple[str, str]] | None = None,
    max_extra_features_per_obsm: int | None = None,
    extra_feature_selection: str = "variance",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(input_h5ad)
    if ko_col not in adata.obs:
        raise ValueError(f"KO column '{ko_col}' was not found in adata.obs.")
    terms = select_pathway_terms(adata, prior_dir, max_terms=max_pathways)
    if not terms:
        raise ValueError("No pathway terms overlap the matrix gene names. Check gene symbols and GMT files.")
    frame = compute_pathway_scores(adata, terms)
    frame.insert(0, "ko_target", adata.obs[ko_col].astype(str).to_numpy())
    frame.insert(0, "cell_id", adata.obs_names.astype(str))
    extra = list(extra_obsm or [])
    if protein_obsm:
        extra.append((protein_obsm, protein_prefix))
    if extra:
        frame = append_extra_obsm_scores(
            frame,
            adata,
            extra,
            max_features_per_obsm=max_extra_features_per_obsm,
            feature_selection=extra_feature_selection,
            labels=adata.obs[ko_col].astype(str),
        )
    manifest = pd.DataFrame(
        {
            "state_feature": [label for label, _ in terms],
            "source": "RNA pathway score",
            "n_genes_used": [len(genes) for _, genes in terms],
        }
    )
    if extra:
        manifest = pd.concat([manifest, extra_obsm_manifest(frame, extra)], ignore_index=True)
    return frame, manifest


def h5ad_to_state_scores(
    input_h5ad: str | Path,
    prior_dir: str | Path,
    max_pathways: int = 40,
    protein_obsm: str | None = None,
    protein_prefix: str = "protein",
    extra_obsm: list[tuple[str, str]] | None = None,
    max_extra_features_per_obsm: int | None = None,
    extra_feature_selection: str = "variance",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(input_h5ad)
    terms = select_pathway_terms(adata, prior_dir, max_terms=max_pathways)
    if not terms:
        raise ValueError("No pathway terms overlap the matrix gene names. Check gene symbols and GMT files.")
    frame = compute_pathway_scores(adata, terms)
    frame.insert(0, "cell_id", adata.obs_names.astype(str))
    extra = list(extra_obsm or [])
    if protein_obsm:
        extra.append((protein_obsm, protein_prefix))
    if extra:
        frame = append_extra_obsm_scores(
            frame,
            adata,
            extra,
            max_features_per_obsm=max_extra_features_per_obsm,
            feature_selection=extra_feature_selection,
            labels=None,
        )
    manifest = pd.DataFrame(
        {
            "state_feature": [label for label, _ in terms],
            "source": "RNA pathway score",
            "n_genes_used": [len(genes) for _, genes in terms],
        }
    )
    if extra:
        manifest = pd.concat([manifest, extra_obsm_manifest(frame, extra)], ignore_index=True)
    return frame, manifest


def h5ad_to_state_scores_with_terms(
    input_h5ad: str | Path,
    terms: list[dict],
    protein_obsm: str | None = None,
    protein_prefix: str = "protein",
    extra_obsm: list[tuple[str, str]] | None = None,
    max_extra_features_per_obsm: int | None = None,
    extra_feature_selection: str = "variance",
    obs_columns: list[str] | None = None,
) -> pd.DataFrame:
    adata = ad.read_h5ad(input_h5ad)
    frame = compute_pathway_scores_from_saved_terms(adata, terms)
    frame.insert(0, "cell_id", adata.obs_names.astype(str))
    for col in obs_columns or []:
        if col in adata.obs:
            frame[col] = adata.obs[col].astype(str).to_numpy()
    extra = list(extra_obsm or [])
    if protein_obsm:
        extra.append((protein_obsm, protein_prefix))
    if extra:
        frame = append_extra_obsm_scores(
            frame,
            adata,
            extra,
            max_features_per_obsm=max_extra_features_per_obsm,
            feature_selection=extra_feature_selection,
            labels=None,
        )
    return frame


def csv_matrix_to_state_table(
    input_csv: str | Path,
    ko_col: str,
    prior_dir: str | Path,
    max_pathways: int = 40,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(input_csv)
    if ko_col not in raw.columns:
        raise ValueError(f"KO column '{ko_col}' was not found in the CSV.")
    gene_cols = [
        col
        for col in raw.columns
        if col != ko_col and pd.api.types.is_numeric_dtype(raw[col]) and GENE_RE.match(str(col).upper())
    ]
    if len(gene_cols) < 20:
        raise ValueError("CSV matrix needs numeric gene-symbol columns plus a KO label column.")
    adata = ad.AnnData(raw[gene_cols].to_numpy(dtype=np.float32))
    adata.var_names = [str(col).upper() for col in gene_cols]
    adata.obs_names = raw.index.astype(str)
    terms = select_pathway_terms(adata, prior_dir, max_terms=max_pathways)
    frame = compute_pathway_scores(adata, terms)
    frame.insert(0, "ko_target", raw[ko_col].astype(str).to_numpy())
    frame.insert(0, "cell_id", raw.index.astype(str))
    manifest = pd.DataFrame(
        {
            "state_feature": [label for label, _ in terms],
            "source": "RNA pathway score",
            "n_genes_used": [len(genes) for _, genes in terms],
        }
    )
    return frame, manifest
