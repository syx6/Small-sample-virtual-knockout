from __future__ import annotations

import pickle
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
) -> tuple[pd.DataFrame, list[dict]]:
    adata = ad.read_h5ad(input_h5ad)
    if ko_col not in adata.obs:
        raise ValueError(f"KO column '{ko_col}' was not found in adata.obs.")
    state_terms = select_pathway_terms(adata, prior_dir, max_terms=max_pathways)
    frame = compute_pathway_scores(adata, state_terms)
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
    return frame, _serialize_state_terms(state_terms)


def _state_features(frame: pd.DataFrame) -> list[str]:
    ignored = {"cell_id", "ko_target", "dataset", "state"}
    return [col for col in frame.columns if col not in ignored and pd.api.types.is_numeric_dtype(frame[col])]


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
        )
        input_kind = "h5ad"
    elif state_csv:
        frame = pd.read_csv(state_csv).rename(columns={ko_col: "ko_target"})
        state_terms = []
        input_kind = "state_csv"
    else:
        raise ValueError("Provide --input-h5ad or --state-csv.")

    features = _state_features(frame)
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
    reference = {
        "version": 1,
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
        "feature_mean": feature_mean,
        "feature_std": feature_std,
        "training_distance_quantiles": {
            "q50": float(np.nanquantile(train_distance, 0.50)),
            "q90": float(np.nanquantile(train_distance, 0.90)),
            "q95": float(np.nanquantile(train_distance, 0.95)),
        },
        "training_genes": sorted(perturb_genes),
    }
    output_model = Path(output_model)
    output_model.parent.mkdir(parents=True, exist_ok=True)
    with output_model.open("wb") as handle:
        pickle.dump(reference, handle)
    return reference


def load_reference_model(path: str | Path) -> dict:
    with Path(path).open("rb") as handle:
        return pickle.load(handle)


def _align_features(frame: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    aligned = frame.copy()
    for feature in features:
        if feature not in aligned.columns:
            aligned[feature] = 0.0
    return aligned[["cell_id", *features]]


def apply_reference_model(
    reference_model: str | Path,
    input_h5ad: str | Path | None,
    state_csv: str | Path | None,
    target_kos: list[str],
    out_dir: str | Path,
    max_cells: int = 800,
    seed: int = 7,
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
        )
    elif state_csv:
        frame = pd.read_csv(state_csv)
    else:
        raise ValueError("Provide --input-h5ad or --state-csv.")

    features = reference["features"]
    frame = _align_features(frame, features)
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
        delta = reference["model"].predict(ko_prior_vector(ko, prior_terms).reshape(1, -1)).reshape(-1)
        virtual = x + delta.reshape(1, -1)
        delta_rows.append({"ko_target": ko, **{f"pred_delta_{feature}": value for feature, value in zip(features, delta)}})
        before = base[["cell_id", *features]].copy()
        before["ko_target"] = ko
        before["state"] = "input cells"
        after = pd.DataFrame(virtual, columns=features)
        after.insert(0, "cell_id", base["cell_id"].to_numpy())
        after["ko_target"] = ko
        after["state"] = "virtual KO cells"
        all_cells.extend([before, after])

    cells = pd.concat(all_cells, ignore_index=True)
    deltas = pd.DataFrame(delta_rows)
    cells.to_csv(out_dir / "applied_virtual_cells.csv", index=False)
    deltas.to_csv(out_dir / "predicted_ko_delta.csv", index=False)
    confidence = _write_transfer_confidence(reference, base, features, target_kos, out_dir)
    _plot_apply_delta(deltas, out_dir)
    _plot_apply_pca(cells, features, out_dir)
    _write_apply_report(reference, target_kos, out_dir, confidence)
    return cells, deltas


def _write_transfer_confidence(reference: dict, base: pd.DataFrame, features: list[str], target_kos: list[str], out_dir: Path) -> pd.DataFrame:
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

    for ko in target_kos:
        genes = set(split_ko(ko))
        covered = sorted(genes & training_genes)
        missing = sorted(genes - training_genes)
        if pd.isna(frac_beyond_q95):
            confidence = "unknown"
            reason = "reference model was trained before confidence metadata existed"
        elif missing:
            confidence = "low"
            reason = "one or more requested KO genes were not present in reference training labels"
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


def _write_apply_report(reference: dict, target_kos: list[str], out_dir: Path, confidence: pd.DataFrame | None = None) -> None:
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

## Target KO

{', '.join(target_kos)}

## Important Interpretation

These outputs are predictions only. If the input dataset does not contain real KO labels, the software cannot calculate true-vs-virtual accuracy, AUC, or distribution improvement on that dataset.

## Transfer Confidence

{confidence_text}

Generated files:

- `applied_virtual_cells.csv`
- `predicted_ko_delta.csv`
- `transfer_confidence.csv`
- `01_predicted_ko_delta_heatmap.png`
- `02_input_vs_virtual_pca.png`
- `03_transfer_confidence.png`
"""
    (out_dir / "apply_report.md").write_text(text, encoding="utf-8")
