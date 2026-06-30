from __future__ import annotations

import pickle
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

from .core import (
    build_delta_table,
    control_mask,
    fit_pls,
    ko_prior_vector,
    select_prior_terms,
    split_ko,
)
from .interaction import fit_interaction_residual_model, predict_interaction_residual_delta
from .preprocess import (
    append_obsm_scores,
    append_extra_obsm_scores,
    compute_pathway_scores,
    h5ad_to_state_scores_with_terms,
    select_pathway_terms,
)
from .visualization import setup_plot, wrap_label


def _serialize_state_terms(terms: list[tuple[str, list[str]]]) -> list[dict]:
    return [{"state_feature": label, "genes": list(genes)} for label, genes in terms]


def _serialize_prior_terms(terms: list[tuple[str, set[str]]]) -> list[dict]:
    rows = []
    for item in terms:
        label, genes = item[0], item[1]
        weight = float(item[2]) if len(item) >= 3 else 1.0
        rows.append({"term": label, "genes": sorted(genes), "weight": weight})
    return rows


def _deserialize_prior_terms(terms: list[dict]) -> list[tuple[str, set[str]]]:
    return [(item["term"], set(item["genes"]), float(item.get("weight", 1.0))) for item in terms]


def _prior_coverage_rows(reference: dict, target_kos: list[str]) -> list[dict]:
    prior_terms = _deserialize_prior_terms(reference.get("prior_terms") or [])
    training_genes = set(reference.get("training_genes") or [])
    rows = []
    for ko in target_kos:
        genes = set(split_ko(ko))
        hits = []
        weighted_hits = 0.0
        for term_name, members, weight in prior_terms:
            overlap = sorted(genes & members)
            if not overlap:
                continue
            hits.append((term_name, overlap, weight))
            weighted_hits += float(weight) * len(overlap) / max(1, len(genes))
        top_terms = "; ".join(f"{name}({'+'.join(overlap)})" for name, overlap, _ in hits[:8])
        rows.append(
            {
                "ko_target": ko,
                "n_genes": len(genes),
                "n_prior_terms_hit": len(hits),
                "weighted_prior_hit_score": weighted_hits,
                "genes_seen_in_reference": "+".join(sorted(genes & training_genes)),
                "genes_unseen_in_reference": "+".join(sorted(genes - training_genes)),
                "top_prior_terms": top_terms,
            }
        )
    return rows


def _load_training_state_from_h5ad(
    input_h5ad: str | Path,
    ko_col: str,
    prior_dir: str | Path,
    max_pathways: int,
    protein_obsm: str | None,
    protein_prefix: str,
    extra_obsm: list[tuple[str, str]] | None = None,
    max_extra_features_per_obsm: int | None = None,
    extra_feature_selection: str = "variance",
    obs_columns: list[str] | None = None,
) -> tuple[pd.DataFrame, list[dict]]:
    adata = ad.read_h5ad(input_h5ad)
    if ko_col not in adata.obs:
        raise ValueError(f"KO column '{ko_col}' was not found in adata.obs.")
    state_terms = select_pathway_terms(adata, prior_dir, max_terms=max_pathways)
    frame = compute_pathway_scores(adata, state_terms)
    frame.insert(0, "ko_target", adata.obs[ko_col].astype(str).to_numpy())
    frame.insert(0, "cell_id", adata.obs_names.astype(str))
    for col in obs_columns or []:
        if col and col in adata.obs and col != ko_col:
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
            labels=adata.obs[ko_col].astype(str),
        )
    return frame, _serialize_state_terms(state_terms)


def _state_features(frame: pd.DataFrame) -> list[str]:
    ignored = {"cell_id", "ko_target", "dataset", "state", "batch", "sample", "donor"}
    return [col for col in frame.columns if col not in ignored and pd.api.types.is_numeric_dtype(frame[col])]


def _apply_batch_control_centering(
    frame: pd.DataFrame,
    features: list[str],
    batch_col: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not batch_col or batch_col not in frame.columns:
        return frame, pd.DataFrame()
    ctrl = frame.loc[control_mask(frame["ko_target"]), features]
    if ctrl.empty:
        return frame, pd.DataFrame()
    adjusted = frame.copy()
    global_control_mean = ctrl.mean().to_numpy(dtype=float)
    rows = []
    for batch, idx in adjusted.groupby(batch_col, observed=True).groups.items():
        batch_idx = list(idx)
        batch_ctrl = adjusted.loc[
            adjusted.index.isin(batch_idx) & control_mask(adjusted["ko_target"]),
            features,
        ]
        if batch_ctrl.empty:
            offset = np.zeros(len(features), dtype=float)
            n_control = 0
            status = "no_batch_controls"
        else:
            offset = batch_ctrl.mean().to_numpy(dtype=float) - global_control_mean
            n_control = len(batch_ctrl)
            status = "centered_to_global_control"
            adjusted.loc[batch_idx, features] = adjusted.loc[batch_idx, features].to_numpy(dtype=float) - offset.reshape(1, -1)
        rows.append(
            {
                "batch": batch,
                "n_cells": len(batch_idx),
                "n_control_cells": n_control,
                "mean_abs_control_offset_removed": float(np.mean(np.abs(offset))),
                "status": status,
            }
        )
    return adjusted, pd.DataFrame(rows)


def _ko_delta_frame(frame: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    ctrl = frame.loc[control_mask(frame["ko_target"]), features]
    if ctrl.empty:
        return pd.DataFrame()
    control_mean = ctrl.mean().to_numpy(dtype=float)
    rows = []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        values = group[features].mean().to_numpy(dtype=float) - control_mean
        row = {
            "ko_genes": str(ko),
            "n_ko_genes": len(split_ko(str(ko))),
            "n_cells": len(group),
        }
        row.update({feature: value for feature, value in zip(features, values)})
        rows.append(row)
    return pd.DataFrame(rows)


def _residual_profile(frame: pd.DataFrame, features: list[str]) -> dict:
    residuals = []
    for _, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        values = group[features].to_numpy(dtype=float)
        if len(values) >= 2:
            residuals.append(values - np.nanmean(values, axis=0, keepdims=True))
    if residuals:
        bank = np.vstack(residuals)
        std = np.nanstd(bank, axis=0)
        radius = np.sqrt(np.nanmean(bank * bank, axis=1))
        return {
            "feature_residual_std": std,
            "median_residual_radius": float(np.nanmedian(radius)),
            "q90_residual_radius": float(np.nanquantile(radius, 0.90)),
        }
    return {
        "feature_residual_std": np.ones(len(features), dtype=float) * 1e-6,
        "median_residual_radius": 0.0,
        "q90_residual_radius": 0.0,
    }


def train_reference_model(
    input_h5ad: str | Path | None,
    state_csv: str | Path | None,
    ko_col: str,
    prior_dir: str | Path,
    output_model: str | Path,
    max_pathways: int = 40,
    protein_obsm: str | None = None,
    protein_prefix: str = "protein",
    extra_obsm: list[tuple[str, str]] | None = None,
    max_extra_features_per_obsm: int | None = None,
    extra_feature_selection: str = "variance",
    dataset_name: str = "reference perturbation dataset",
    batch_col: str | None = None,
    interaction_mode: str = "auto",
) -> dict:
    if input_h5ad:
        frame, state_terms = _load_training_state_from_h5ad(
            input_h5ad,
            ko_col,
            prior_dir,
            max_pathways,
            protein_obsm,
            protein_prefix,
            extra_obsm,
            max_extra_features_per_obsm,
            extra_feature_selection,
            obs_columns=[batch_col] if batch_col else None,
        )
        input_kind = "h5ad"
    elif state_csv:
        frame = pd.read_csv(state_csv).rename(columns={ko_col: "ko_target"})
        state_terms = []
        input_kind = "state_csv"
    else:
        raise ValueError("Provide --input-h5ad or --state-csv.")

    features = _state_features(frame)
    if batch_col and batch_col in features:
        features = [feature for feature in features if feature != batch_col]
    frame, batch_summary = _apply_batch_control_centering(frame, features, batch_col)
    feature_matrix = frame[features].to_numpy(dtype=float)
    feature_mean = np.nanmean(feature_matrix, axis=0)
    feature_std = np.nanstd(feature_matrix, axis=0)
    feature_std = np.where(feature_std < 1e-6, 1.0, feature_std)
    z = (feature_matrix - feature_mean.reshape(1, -1)) / feature_std.reshape(1, -1)
    train_distance = np.sqrt(np.nanmean(z * z, axis=1))
    perturb_genes = {gene for ko in frame["ko_target"].astype(str).unique() for gene in split_ko(ko)}
    prior_terms = select_prior_terms(Path(prior_dir), perturb_genes)
    train_labels, train_delta, control_mean = build_delta_table(frame, "ko_target", features, holdouts=set())
    model = fit_pls(train_labels, train_delta, prior_terms)
    residual_profile = _residual_profile(frame, features)
    interaction_model = None
    interaction_status = "off"
    if interaction_mode != "off":
        try:
            delta_frame = _ko_delta_frame(frame, features)
            if (
                not delta_frame.empty
                and (delta_frame["n_ko_genes"] == 1).sum() >= 3
                and (delta_frame["n_ko_genes"] == 2).sum() >= 3
            ):
                interaction_model = fit_interaction_residual_model(
                    delta_frame,
                    ko_col="ko_genes",
                    n_ko_col="n_ko_genes",
                    target_cols=features,
                    prior_dir=prior_dir,
                )
                interaction_status = "trained"
            else:
                interaction_status = "not_enough_single_and_double_ko_labels"
        except Exception as exc:
            if interaction_mode == "on":
                raise
            interaction_status = f"failed_auto: {exc}"
    reference = {
        "version": 3,
        "model_type": "prior_constrained_residual_pls_reference",
        "supports_batch_ko": True,
        "supports_prediction_only_application": True,
        "supports_cell_type_stratification": True,
        "supports_double_ko": True,
        "supports_interaction_residual": interaction_model is not None,
        "supports_batch_covariate": bool(batch_col),
        "dataset_name": dataset_name,
        "input_kind": input_kind,
        "features": features,
        "state_terms": state_terms,
        "prior_terms": _serialize_prior_terms(prior_terms),
        "model": model,
        "training_ko_labels": train_labels,
        "control_mean": control_mean,
        "protein_obsm": protein_obsm,
        "protein_prefix": protein_prefix,
        "extra_obsm": extra_obsm or [],
        "max_extra_features_per_obsm": max_extra_features_per_obsm,
        "extra_feature_selection": extra_feature_selection,
        "max_pathways": max_pathways,
        "batch_col": batch_col,
        "batch_summary": batch_summary.to_dict(orient="records"),
        "interaction_mode": interaction_mode,
        "interaction_status": interaction_status,
        "interaction_model": interaction_model,
        "uncertainty_model": {
            "model_type": "hard_constraint_residual_band",
            "feature_residual_std": residual_profile["feature_residual_std"],
            "median_residual_radius": residual_profile["median_residual_radius"],
            "q90_residual_radius": residual_profile["q90_residual_radius"],
            "interpretation": "The KO mean delta remains fixed by the residual/PLS baseline; this profile only describes cell-level uncertainty around that hard constraint.",
        },
        "feature_mean": feature_mean,
        "feature_std": feature_std,
        "training_distance_quantiles": {
            "q50": float(np.nanquantile(train_distance, 0.50)),
            "q90": float(np.nanquantile(train_distance, 0.90)),
            "q95": float(np.nanquantile(train_distance, 0.95)),
        },
        "training_genes": sorted(perturb_genes),
        "interpretation": {
            "purpose": "Predict KO-induced state deltas and apply them to ordinary cells.",
            "accuracy_boundary": "Accuracy metrics require real KO labels in the evaluation dataset. Unlabeled 10X/multiome data are prediction-only.",
            "unseen_gene_boundary": "Unseen genes are predicted by pathway/TF/PPI/motif priors and should be treated as lower confidence when prior coverage is weak.",
        },
    }
    output_model = Path(output_model)
    output_model.parent.mkdir(parents=True, exist_ok=True)
    with output_model.open("wb") as handle:
        pickle.dump(reference, handle)
    _write_reference_metadata(reference, output_model)
    return reference


def _write_reference_metadata(reference: dict, output_model: Path) -> None:
    meta = {
        "version": reference.get("version"),
        "model_type": reference.get("model_type"),
        "dataset_name": reference.get("dataset_name"),
        "input_kind": reference.get("input_kind"),
        "n_training_ko_labels": len(reference.get("training_ko_labels", [])),
        "n_training_genes": len(reference.get("training_genes", [])),
        "n_state_features": len(reference.get("features", [])),
        "training_ko_labels": reference.get("training_ko_labels", []),
        "training_genes": reference.get("training_genes", []),
        "extra_obsm": reference.get("extra_obsm", []),
        "max_extra_features_per_obsm": reference.get("max_extra_features_per_obsm"),
        "extra_feature_selection": reference.get("extra_feature_selection"),
        "batch_col": reference.get("batch_col"),
        "batch_summary": reference.get("batch_summary", []),
        "interaction_status": reference.get("interaction_status"),
        "uncertainty_model": {
            key: value
            for key, value in (reference.get("uncertainty_model") or {}).items()
            if key != "feature_residual_std"
        },
        "supports": {
            "batch_ko": bool(reference.get("supports_batch_ko")),
            "double_ko": bool(reference.get("supports_double_ko")),
            "interaction_residual": bool(reference.get("supports_interaction_residual")),
            "batch_covariate": bool(reference.get("supports_batch_covariate")),
            "hard_constrained_uncertainty": bool(reference.get("uncertainty_model")),
            "cell_type_stratification": bool(reference.get("supports_cell_type_stratification")),
            "prediction_only_application": bool(reference.get("supports_prediction_only_application")),
        },
        "interpretation": reference.get("interpretation", {}),
    }
    metadata_path = output_model.with_suffix(output_model.suffix + ".metadata.json")
    metadata_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")


def load_reference_model(path: str | Path) -> dict:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def inspect_reference_model(reference_model: str | Path, out_dir: str | Path | None = None, target_kos: list[str] | None = None) -> dict:
    reference = load_reference_model(reference_model)
    prior_terms = reference.get("prior_terms") or []
    libraries = {}
    for item in prior_terms:
        term = str(item.get("term", "unknown"))
        lib = term.split(":", 1)[0] if ":" in term else "unknown"
        libraries.setdefault(lib, {"n_terms": 0, "mean_weight": []})
        libraries[lib]["n_terms"] += 1
        libraries[lib]["mean_weight"].append(float(item.get("weight", 1.0)))
    prior_summary = pd.DataFrame(
        [
            {
                "prior_library": lib,
                "n_terms": values["n_terms"],
                "mean_weight": float(np.mean(values["mean_weight"])) if values["mean_weight"] else np.nan,
            }
            for lib, values in sorted(libraries.items())
        ]
    )
    feature_summary = pd.DataFrame(
        {
            "state_feature": reference.get("features", []),
            "source_guess": [
                "protein/ADT" if feature.startswith("protein_") else "ATAC/chromVAR/peak" if any(feature.startswith(prefix) for prefix in ["atac_", "tf_", "peak_", "chromvar_"]) else "pathway/program"
                for feature in reference.get("features", [])
            ],
        }
    )
    gene_summary = pd.DataFrame({"training_gene": reference.get("training_genes", [])})
    target_coverage = pd.DataFrame(_prior_coverage_rows(reference, target_kos or []))
    summary = {
        "version": reference.get("version"),
        "model_type": reference.get("model_type"),
        "dataset_name": reference.get("dataset_name"),
        "n_training_ko_labels": len(reference.get("training_ko_labels", [])),
        "n_training_genes": len(reference.get("training_genes", [])),
        "n_state_features": len(reference.get("features", [])),
        "n_prior_terms": len(prior_terms),
    }
    if out_dir:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([summary]).to_csv(out / "reference_summary.csv", index=False)
        prior_summary.to_csv(out / "reference_prior_libraries.csv", index=False)
        feature_summary.to_csv(out / "reference_state_features.csv", index=False)
        gene_summary.to_csv(out / "reference_training_genes.csv", index=False)
        if not target_coverage.empty:
            target_coverage.to_csv(out / "reference_target_prior_coverage.csv", index=False)
        _write_reference_inspection_report(reference, summary, prior_summary, feature_summary, target_coverage, out)
    return {
        "summary": summary,
        "prior_libraries": prior_summary,
        "state_features": feature_summary,
        "training_genes": gene_summary,
        "target_prior_coverage": target_coverage,
    }


def _write_reference_inspection_report(
    reference: dict,
    summary: dict,
    prior_summary: pd.DataFrame,
    feature_summary: pd.DataFrame,
    target_coverage: pd.DataFrame,
    out_dir: Path,
) -> None:
    coverage_text = "No target KO coverage requested."
    if not target_coverage.empty:
        coverage_text = target_coverage[["ko_target", "n_prior_terms_hit", "weighted_prior_hit_score", "genes_unseen_in_reference"]].round(3).to_string(index=False)
    text = f"""# Reference Model Inspection Report

## Summary

- Dataset: {summary.get('dataset_name')}
- Version: {summary.get('version')}
- Model type: {summary.get('model_type')}
- Training KO labels: {summary.get('n_training_ko_labels')}
- Training genes: {summary.get('n_training_genes')}
- State features: {summary.get('n_state_features')}
- Prior terms: {summary.get('n_prior_terms')}

## Prior Libraries

{prior_summary.round(3).to_string(index=False)}

## Target Prior Coverage

{coverage_text}

## Interpretation

This report checks what the reference model has learned and whether requested KO genes are covered by the biological prior network. A target with few prior hits, unseen genes, or weak weighted prior score should be treated as lower confidence.
"""
    (out_dir / "reference_inspection_report.md").write_text(text, encoding="utf-8")


def _align_features(frame: pd.DataFrame, features: list[str], preserve_cols: list[str] | None = None) -> pd.DataFrame:
    aligned = frame.copy()
    for feature in features:
        if feature not in aligned.columns:
            aligned[feature] = 0.0
    preserve = [col for col in preserve_cols or [] if col and col in aligned.columns and col not in {"cell_id", *features}]
    extra_cols = [
        col
        for col in aligned.columns
        if col not in {"cell_id", *features, *preserve} and not pd.api.types.is_numeric_dtype(aligned[col])
    ]
    return aligned[["cell_id", *preserve, *extra_cols, *features]]


def apply_reference_model(
    reference_model: str | Path,
    input_h5ad: str | Path | None,
    state_csv: str | Path | None,
    target_kos: list[str],
    out_dir: str | Path,
    max_cells: int = 800,
    seed: int = 7,
    cell_type_col: str | None = None,
    batch_col: str | None = None,
    uncertainty_method: str = "none",
    uncertainty_scale: float = 0.25,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    reference = load_reference_model(reference_model)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if input_h5ad:
        frame = h5ad_to_state_scores_with_terms(
            input_h5ad,
            reference["state_terms"],
            protein_obsm=reference.get("protein_obsm"),
            protein_prefix=reference.get("protein_prefix", "protein"),
            extra_obsm=reference.get("extra_obsm") or [],
            max_extra_features_per_obsm=reference.get("max_extra_features_per_obsm"),
            extra_feature_selection=reference.get("extra_feature_selection", "variance"),
            obs_columns=[col for col in [cell_type_col, batch_col] if col],
        )
    elif state_csv:
        frame = pd.read_csv(state_csv)
    else:
        raise ValueError("Provide --input-h5ad or --state-csv.")

    features = reference["features"]
    frame = _align_features(frame, features, preserve_cols=[cell_type_col, batch_col])
    x = frame[features].to_numpy(dtype=float)
    prior_terms = _deserialize_prior_terms(reference["prior_terms"])
    rng = np.random.default_rng(seed)
    if len(frame) > max_cells:
        keep = np.sort(rng.choice(len(frame), size=max_cells, replace=False))
        base = frame.iloc[keep].reset_index(drop=True)
        x = x[keep]
    else:
        base = frame.reset_index(drop=True)

    all_cells, delta_rows = [], []
    for ko in target_kos:
        prediction_source = "prior_constrained_pls"
        delta = None
        if len(split_ko(ko)) == 2 and reference.get("interaction_model") is not None:
            interaction_delta = predict_interaction_residual_delta(reference["interaction_model"], [ko])
            if not interaction_delta.empty:
                delta = interaction_delta.iloc[0][features].to_numpy(dtype=float)
                prediction_source = "interaction_residual"
        if delta is None:
            delta = reference["model"].predict(ko_prior_vector(ko, prior_terms).reshape(1, -1)).reshape(-1)
        virtual = x + delta.reshape(1, -1)
        delta_rows.append({"ko_target": ko, "prediction_source": prediction_source, **{f"pred_delta_{feature}": value for feature, value in zip(features, delta)}})
        before_cols = [
            "cell_id",
            *([cell_type_col] if cell_type_col and cell_type_col in base.columns else []),
            *([batch_col] if batch_col and batch_col in base.columns else []),
            *features,
        ]
        before = base[before_cols].copy()
        before["ko_target"] = ko
        before["state"] = "input cells"
        after = pd.DataFrame(virtual, columns=features)
        after.insert(0, "cell_id", base["cell_id"].to_numpy())
        if cell_type_col and cell_type_col in base.columns:
            after[cell_type_col] = base[cell_type_col].to_numpy()
        if batch_col and batch_col in base.columns:
            after[batch_col] = base[batch_col].to_numpy()
        after["ko_target"] = ko
        after["state"] = "virtual KO cells"
        all_cells.extend([before, after])

    cells = pd.concat(all_cells, ignore_index=True)
    deltas = pd.DataFrame(delta_rows)
    cells.to_csv(out_dir / "applied_virtual_cells.csv", index=False)
    deltas.to_csv(out_dir / "predicted_ko_delta.csv", index=False)
    prior_coverage = _write_prior_coverage(reference, target_kos, out_dir)
    confidence = _write_transfer_confidence(reference, base, features, target_kos, out_dir, prior_coverage=prior_coverage)
    _write_target_interpretation(reference, target_kos, out_dir, prior_coverage=prior_coverage)
    if cell_type_col and cell_type_col in cells.columns:
        _write_cell_type_outputs(cells, features, cell_type_col, out_dir)
    if batch_col and batch_col in cells.columns:
        _write_batch_outputs(cells, batch_col, out_dir)
    _plot_apply_delta(deltas, out_dir)
    if uncertainty_method != "none":
        _write_uncertainty_outputs(deltas, reference, features, out_dir, uncertainty_method, uncertainty_scale)
    _plot_apply_pca(cells, features, out_dir)
    _write_apply_report(
        reference,
        target_kos,
        out_dir,
        confidence,
        cell_type_col=cell_type_col,
        batch_col=batch_col,
        uncertainty_method=uncertainty_method,
    )
    return cells, deltas


def _write_target_interpretation(reference: dict, target_kos: list[str], out_dir: Path, prior_coverage: pd.DataFrame | None = None) -> pd.DataFrame:
    training_genes = set(reference.get("training_genes") or [])
    prior_lookup = {}
    if prior_coverage is not None and not prior_coverage.empty:
        prior_lookup = prior_coverage.set_index("ko_target").to_dict(orient="index")
    rows = []
    for ko in target_kos:
        genes = split_ko(ko)
        missing = [gene for gene in genes if gene not in training_genes]
        if len(genes) == 1:
            mode = "single_ko"
            note = "single-gene reference application"
        elif len(genes) == 2:
            if reference.get("interaction_model") is not None:
                mode = "double_ko_interaction_residual"
                note = "double-KO prediction uses the trained interaction residual model on top of single-gene additive effects"
            else:
                mode = "double_ko_prior_delta"
                note = "double-KO prediction uses the learned prior-constrained delta for the pair; train a reference with double-KO labels to enable interaction residuals"
        else:
            mode = "multi_gene_exploratory"
            note = "more than two genes is exploratory and should be treated as low confidence unless specifically validated"
        rows.append(
            {
                "ko_target": ko,
                "n_genes": len(genes),
                "analysis_mode": mode,
                "genes_seen_in_reference": "+".join([gene for gene in genes if gene in training_genes]),
                "genes_unseen_in_reference": "+".join(missing),
                "n_prior_terms_hit": prior_lookup.get(ko, {}).get("n_prior_terms_hit", np.nan),
                "weighted_prior_hit_score": prior_lookup.get(ko, {}).get("weighted_prior_hit_score", np.nan),
                "interpretation": note,
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "target_interpretation.csv", index=False)
    return table


def _write_prior_coverage(reference: dict, target_kos: list[str], out_dir: Path) -> pd.DataFrame:
    table = pd.DataFrame(_prior_coverage_rows(reference, target_kos))
    table.to_csv(out_dir / "prior_coverage.csv", index=False)
    _plot_prior_coverage(table, out_dir)
    return table


def _plot_prior_coverage(table: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    if table.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6.5, 2.5 * len(table)), 4.2), constrained_layout=True)
    bars = ax.bar(table["ko_target"], table["n_prior_terms_hit"], color="#4C78A8", width=0.58)
    ax.set_ylabel("Prior terms hit")
    ax.set_title("KO Prior Coverage")
    ax.set_xlabel("Requested virtual KO")
    for bar, (_, row) in zip(bars, table.iterrows()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            f"score {row['weighted_prior_hit_score']:.1f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )
    fig.savefig(out_dir / "05_prior_coverage.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_cell_type_outputs(cells: pd.DataFrame, features: list[str], cell_type_col: str, out_dir: Path) -> pd.DataFrame:
    rows = []
    grouped = cells.groupby(["ko_target", cell_type_col, "state"], observed=True)
    means = grouped[features].mean().reset_index()
    for (ko, cell_type), sub in means.groupby(["ko_target", cell_type_col], observed=True):
        before = sub.loc[sub["state"] == "input cells", features]
        after = sub.loc[sub["state"] == "virtual KO cells", features]
        if before.empty or after.empty:
            continue
        delta = after.iloc[0].to_numpy(dtype=float) - before.iloc[0].to_numpy(dtype=float)
        for feature, value in zip(features, delta):
            rows.append({"ko_target": ko, cell_type_col: cell_type, "feature": feature, "predicted_delta": value})
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "cell_type_predicted_delta.csv", index=False)
    if not table.empty:
        _plot_cell_type_delta(table, cell_type_col, out_dir)
    return table


def _plot_cell_type_delta(table: pd.DataFrame, cell_type_col: str, out_dir: Path, max_features: int = 12) -> None:
    setup_plot()
    score = table.groupby("feature")["predicted_delta"].apply(lambda x: np.max(np.abs(x))).sort_values(ascending=False)
    keep = list(score.head(max_features).index)
    plot = table.loc[table["feature"].isin(keep)].copy()
    plot["target_cell_type"] = plot["ko_target"].astype(str) + " | " + plot[cell_type_col].astype(str)
    heat = plot.pivot_table(index="feature", columns="target_cell_type", values="predicted_delta", aggfunc="mean")
    heat.index = [wrap_label(idx.replace("pathway_", "").replace("protein_", "").replace("program_", ""), 22) for idx in heat.index]
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(heat.columns) + 4), max(5, 0.45 * len(heat) + 2)), constrained_layout=True)
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=False, cbar_kws={"label": "predicted cell-type delta"}, ax=ax)
    ax.set_title("Prediction-only cell-type stratified virtual KO effects")
    ax.set_xlabel("KO | cell type")
    ax.set_ylabel("")
    fig.savefig(out_dir / "04_cell_type_predicted_delta_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_batch_outputs(cells: pd.DataFrame, batch_col: str, out_dir: Path) -> pd.DataFrame:
    table = (
        cells.groupby(["ko_target", batch_col, "state"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    table.to_csv(out_dir / "batch_composition.csv", index=False)
    setup_plot()
    fig, ax = plt.subplots(figsize=(max(6.5, 2.6 * table["ko_target"].nunique()), 4.6), constrained_layout=True)
    plot = table.loc[table["state"] == "input cells"].copy()
    if not plot.empty:
        sns.barplot(data=plot, x="ko_target", y="n_cells", hue=batch_col, ax=ax)
        ax.set_title("Input Cells by Batch")
        ax.set_xlabel("Requested virtual KO")
        ax.set_ylabel("Cells")
        ax.legend(title=batch_col, bbox_to_anchor=(1.02, 1), loc="upper left")
        fig.savefig(out_dir / "06_batch_composition.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    return table


def _write_transfer_confidence(
    reference: dict,
    base: pd.DataFrame,
    features: list[str],
    target_kos: list[str],
    out_dir: Path,
    prior_coverage: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows = []
    feature_mean = reference.get("feature_mean")
    feature_std = reference.get("feature_std")
    quantiles = reference.get("training_distance_quantiles") or {}
    training_genes = set(reference.get("training_genes") or [])
    if feature_mean is not None and feature_std is not None:
        mean = np.asarray(feature_mean, dtype=float)
        std = np.asarray(feature_std, dtype=float)
        std = np.where(std < 1e-6, 1.0, std)
        x = base[features].to_numpy(dtype=float)
        z = (x - mean.reshape(1, -1)) / std.reshape(1, -1)
        distances = np.sqrt(np.nanmean(z * z, axis=1))
        median_distance = float(np.nanmedian(distances))
        frac_beyond_q95 = float(np.mean(distances > float(quantiles.get("q95", np.inf))))
        reference_q95 = float(quantiles.get("q95", np.nan))
    else:
        median_distance = np.nan
        frac_beyond_q95 = np.nan
        reference_q95 = np.nan
    prior_lookup = {}
    if prior_coverage is not None and not prior_coverage.empty:
        prior_lookup = prior_coverage.set_index("ko_target").to_dict(orient="index")

    for ko in target_kos:
        genes = set(split_ko(ko))
        covered = sorted(genes & training_genes)
        missing = sorted(genes - training_genes)
        prior_hit_count = float(prior_lookup.get(ko, {}).get("n_prior_terms_hit", np.nan))
        prior_hit_score = float(prior_lookup.get(ko, {}).get("weighted_prior_hit_score", np.nan))
        if pd.isna(frac_beyond_q95):
            confidence = "unknown"
            reason = "reference model was trained before confidence metadata existed"
        elif missing:
            confidence = "low"
            reason = "one or more requested KO genes were not present in reference training labels"
        elif pd.notna(prior_hit_count) and prior_hit_count < 2:
            confidence = "low"
            reason = "requested KO has weak pathway/TF/PPI/motif prior coverage"
        elif frac_beyond_q95 > 0.35:
            confidence = "low"
            reason = "many input cells are outside the reference state distribution"
        elif frac_beyond_q95 > 0.15:
            confidence = "medium"
            reason = "some input cells are outside the reference state distribution"
        else:
            confidence = "high"
            reason = "input cells are close to the reference state distribution and KO genes were seen"
        rows.append(
            {
                "ko_target": ko,
                "confidence": confidence,
                "reason": reason,
                "median_reference_distance": median_distance,
                "reference_q95_distance": reference_q95,
                "fraction_cells_beyond_reference_q95": frac_beyond_q95,
                "n_prior_terms_hit": prior_hit_count,
                "weighted_prior_hit_score": prior_hit_score,
                "ko_genes_seen_in_reference": "+".join(covered),
                "ko_genes_missing_from_reference": "+".join(missing),
            }
        )
    confidence = pd.DataFrame(rows)
    confidence.to_csv(out_dir / "transfer_confidence.csv", index=False)
    _plot_transfer_confidence(confidence, out_dir)
    return confidence


def _plot_transfer_confidence(confidence: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    if confidence.empty:
        return
    score_map = {"low": 0.25, "medium": 0.60, "high": 0.90, "unknown": 0.40}
    colors = {"low": "#D55E00", "medium": "#E69F00", "high": "#009E73", "unknown": "#999999"}
    plot = confidence.copy()
    plot["score"] = plot["confidence"].map(score_map).fillna(0.4)
    fig, ax = plt.subplots(figsize=(max(6.5, 2.8 * len(plot)), 4.2), constrained_layout=True)
    bars = ax.bar(plot["ko_target"], plot["score"], width=0.55, color=[colors.get(x, "#999999") for x in plot["confidence"]])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Transfer confidence")
    ax.set_title("Reference transfer confidence")
    ax.set_xlabel("Requested virtual KO")
    ax.set_yticks([0.25, 0.60, 0.90])
    ax.set_yticklabels(["low", "medium", "high"])
    for bar, (_, row) in zip(bars, plot.iterrows()):
        x = bar.get_x() + bar.get_width() / 2
        ax.text(x, bar.get_height() + 0.035, row["confidence"], ha="center", fontsize=11, fontweight="bold")
        beyond = row.get("fraction_cells_beyond_reference_q95", np.nan)
        missing = row.get("ko_genes_missing_from_reference", "")
        note = f"outside ref q95: {beyond:.1%}" if pd.notna(beyond) else "outside ref q95: n/a"
        if isinstance(missing, str) and missing:
            note += f"\nmissing: {missing}"
        else:
            note += "\nKO genes seen"
        ax.text(x, 0.08, note, ha="center", va="bottom", fontsize=8.5, color="white")
    fig.savefig(out_dir / "03_transfer_confidence.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_uncertainty_outputs(
    deltas: pd.DataFrame,
    reference: dict,
    features: list[str],
    out_dir: Path,
    method: str,
    scale: float,
    max_features: int = 12,
) -> pd.DataFrame:
    profile = reference.get("uncertainty_model") or {}
    residual_std = np.asarray(profile.get("feature_residual_std", np.ones(len(features)) * 1e-6), dtype=float)
    if residual_std.shape[0] != len(features):
        residual_std = np.resize(residual_std, len(features))
    rows = []
    for _, row in deltas.iterrows():
        for feature, std in zip(features, residual_std):
            pred = float(row[f"pred_delta_{feature}"])
            width = float(abs(scale) * std)
            rows.append(
                {
                    "ko_target": row["ko_target"],
                    "feature": feature,
                    "uncertainty_method": method,
                    "prediction_source": row.get("prediction_source", "unknown"),
                    "predicted_delta": pred,
                    "lower_delta": pred - width,
                    "upper_delta": pred + width,
                    "interval_half_width": width,
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(out_dir / "uncertainty_intervals.csv", index=False)
    score = table.groupby("feature")["predicted_delta"].apply(lambda x: np.max(np.abs(x))).sort_values(ascending=False)
    keep = set(score.head(max_features).index)
    plot = table.loc[table["feature"].isin(keep)].copy()
    if plot.empty:
        return table
    plot["target_feature"] = plot["ko_target"].astype(str) + " | " + plot["feature"].astype(str).str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    plot = plot.sort_values("predicted_delta")
    setup_plot()
    fig, ax = plt.subplots(figsize=(8.4, max(5, 0.45 * len(plot))), constrained_layout=True)
    ax.errorbar(
        plot["predicted_delta"],
        plot["target_feature"],
        xerr=[plot["predicted_delta"] - plot["lower_delta"], plot["upper_delta"] - plot["predicted_delta"]],
        fmt="o",
        color="#4C78A8",
        ecolor="#9ECAE1",
        elinewidth=2,
        capsize=3,
    )
    ax.axvline(0, color="0.25", linestyle="--", linewidth=1)
    ax.set_title("Hard-Constrained KO Uncertainty Band")
    ax.set_xlabel("Predicted delta with residual uncertainty")
    ax.set_ylabel("")
    fig.savefig(out_dir / "07_uncertainty_intervals.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    return table


def _plot_apply_delta(deltas: pd.DataFrame, out_dir: Path, max_features: int = 10) -> None:
    setup_plot()
    rows = []
    for _, row in deltas.iterrows():
        ko = row["ko_target"]
        for col in [c for c in deltas.columns if c.startswith("pred_delta_")]:
            feature = col.removeprefix("pred_delta_")
            rows.append({"ko_target": ko, "feature": feature, "predicted_delta": float(row[col])})
    plot = pd.DataFrame(rows)
    score = plot.groupby("feature")["predicted_delta"].apply(lambda x: np.max(np.abs(x))).sort_values(ascending=False)
    keep = list(score.head(max_features).index)
    heat = plot.loc[plot["feature"].isin(keep)].pivot(index="feature", columns="ko_target", values="predicted_delta")
    heat.index = [wrap_label(idx.replace("pathway_", "").replace("protein_", "").replace("program_", ""), 20) for idx in heat.index]
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    fig, ax = plt.subplots(figsize=(max(5, 1.3 * len(heat.columns) + 4), max(5, 0.45 * len(heat) + 2)), constrained_layout=True)
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", annot_kws={"fontsize": 8}, cbar_kws={"label": "predicted KO delta"}, ax=ax)
    ax.set_title("Predicted Virtual KO Effects")
    ax.set_xlabel("Target KO")
    ax.set_ylabel("")
    fig.savefig(out_dir / "01_predicted_ko_delta_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_apply_pca(cells: pd.DataFrame, features: list[str], out_dir: Path) -> None:
    setup_plot()
    matrix = StandardScaler().fit_transform(cells[features].to_numpy(dtype=float))
    coords = PCA(n_components=2, random_state=0).fit_transform(matrix)
    plot = cells[["ko_target", "state"]].copy()
    plot["PC1"] = coords[:, 0]
    plot["PC2"] = coords[:, 1]
    kos = list(plot["ko_target"].drop_duplicates())
    fig, axes = plt.subplots(1, len(kos), figsize=(max(6, 5 * len(kos)), 4.6), squeeze=False, constrained_layout=True)
    for ax, ko in zip(axes.flat, kos):
        sub = plot.loc[plot["ko_target"] == ko]
        sns.scatterplot(data=sub, x="PC1", y="PC2", hue="state", palette={"input cells": "#BDBDBD", "virtual KO cells": "#E76F51"}, s=20, alpha=0.65, linewidth=0, ax=ax)
        ax.set_title(f"{ko}: input vs virtual KO")
        ax.legend(title="")
    fig.suptitle("Reference Model Applied to Unlabeled Cells")
    fig.savefig(out_dir / "02_input_vs_virtual_pca.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_apply_report(
    reference: dict,
    target_kos: list[str],
    out_dir: Path,
    confidence: pd.DataFrame | None = None,
    cell_type_col: str | None = None,
    batch_col: str | None = None,
    uncertainty_method: str = "none",
) -> None:
    confidence_text = "not available"
    if confidence is not None and not confidence.empty:
        confidence_text = "\n" + confidence.round(3).to_string(index=False)
    text = f"""# Reference Virtual KO Application Report

This output applies a reference perturbation model to input cells that may not have KO labels.

## Reference Model

- Reference dataset: {reference.get('dataset_name', 'unknown')}
- Training KO labels: {len(reference.get('training_ko_labels', []))}
- State features: {len(reference.get('features', []))}
- Legacy protein/extra modality key: {reference.get('protein_obsm') or 'none'}
- Extra modalities: {reference.get('extra_obsm') or []}
- Batch covariate used during reference training: {reference.get('batch_col') or 'not provided'}
- Double-KO interaction residual: {reference.get('interaction_status') or 'not trained'}
- Hard-constrained uncertainty model: {uncertainty_method}

## Target KO

{', '.join(target_kos)}

## Important Interpretation

These outputs are predictions only. If the input dataset does not contain real KO labels, the software cannot calculate true-vs-virtual accuracy, AUC, or distribution improvement on that dataset.

This is the correct mode for ordinary 10X, DOGMA-seq, TEA-seq or other unlabeled multiome data. The plots show predicted state shifts, not validated accuracy.

## Batch and Cell-Type Output

- Batch KO targets requested: {len(target_kos)}
- Cell-type stratification column: {cell_type_col or 'not provided'}
- Application batch column: {batch_col or 'not provided'}
- Uncertainty output: {uncertainty_method}

## Transfer Confidence

{confidence_text}

Generated files:

- `applied_virtual_cells.csv`
- `predicted_ko_delta.csv`
- `target_interpretation.csv`
- `prior_coverage.csv`
- `transfer_confidence.csv`
- `01_predicted_ko_delta_heatmap.png`
- `02_input_vs_virtual_pca.png`
- `03_transfer_confidence.png`
- `05_prior_coverage.png`
- `06_batch_composition.png` when a batch column is provided.
- `uncertainty_intervals.csv` and `07_uncertainty_intervals.png` when uncertainty is enabled.
- `cell_type_predicted_delta.csv` and `04_cell_type_predicted_delta_heatmap.png` when a cell-type column is provided.
"""
    (out_dir / "apply_report.md").write_text(text, encoding="utf-8")
    prediction = f"""# Prediction-Only 10X / Multiome Application Report

This folder contains virtual knockout predictions for input cells.

## What the Results Mean

- The reference model predicts how pathway/protein/ATAC state features may move after the requested KO.
- The input cells are shifted by the predicted KO delta to create virtual KO cells.
- If the input data do not contain real KO labels, these results are not an internal accuracy benchmark.

## What Is Not Reported

- No true-vs-virtual accuracy.
- No dataset-specific ROC-AUC.
- No dataset-specific R2/MAE.
- No claim that the prediction has been experimentally validated in this input dataset.

## What To Use

- `01_predicted_ko_delta_heatmap.png`: which state features are predicted to increase or decrease.
- `02_input_vs_virtual_pca.png`: visible movement from input cells to virtual KO cells.
- `03_transfer_confidence.png`: whether the input cells and KO targets are close to the reference training regime.
- `04_cell_type_predicted_delta_heatmap.png`: cell-type-specific predicted changes, when available.
- `05_prior_coverage.png`: whether requested KO genes are covered by pathway/TF/PPI/motif priors.
- `06_batch_composition.png`: input-cell batch composition, when a batch column is provided.
- `07_uncertainty_intervals.png`: residual uncertainty band around the hard-constrained KO delta, when enabled.
"""
    (out_dir / "prediction_only_report.md").write_text(prediction, encoding="utf-8")
