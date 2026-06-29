from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")


@dataclass
class VirtualKOResult:
    metrics: pd.DataFrame
    summary: pd.DataFrame
    virtual_cells: pd.DataFrame
    delta_table: pd.DataFrame
    auc_points: pd.DataFrame
    calibration: pd.DataFrame


def control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("control|non-target|nontarget|negative|safe|nt")


def clean_gene(label: str) -> str:
    text = str(label).upper()
    text = re.sub(r"^P-SG", "", text)
    text = re.sub(r"^SG", "", text)
    text = re.sub(r"_[0-9]+$", "", text)
    text = re.sub(r"-[0-9]+$", "", text)
    text = text.replace("P_", "").replace("P-", "")
    return text


def split_ko(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    genes = [clean_gene(part.strip()) for part in text.split("_") if part.strip()]
    out = []
    for gene in genes:
        low = gene.lower()
        if low.startswith("nt") or low in {"ctrl", "control", "nan", "none"}:
            continue
        if GENE_RE.match(gene) and not gene.isdigit():
            out.append(gene)
    return out


def parse_gmt(path: Path, include_term_gene: bool = False) -> list[tuple[str, set[str]]]:
    terms = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = {gene.upper() for gene in parts[2:] if GENE_RE.match(gene.upper())}
            if include_term_gene:
                token = term.split()[0].upper() if term.split() else ""
                if GENE_RE.match(token):
                    genes.add(token)
            if genes:
                terms.append((term, genes))
    return terms


def _prior_should_include_term_gene(stem: str) -> bool:
    return stem in {"ppi_hub", "tf_target", "atac_motif_tf", "motif_tf_target"}


def _library_weight(stem: str) -> float:
    weights = {
        "tf_target": 1.45,
        "atac_motif_tf": 1.60,
        "motif_tf_target": 1.60,
        "ppi_hub": 1.25,
        "reactome": 1.05,
        "hallmark": 1.00,
    }
    return weights.get(stem, 1.0)


def _mechanism_weight(term: str) -> float:
    upper = term.upper()
    if any(word in upper for word in ["TGFB", "TGF_BETA", "TGF-BETA", "SMAD"]):
        return 1.35
    if any(word in upper for word in ["MAPK", "ERK", "JNK", "P38", "RAS", "RAF", "MEK"]):
        return 1.30
    if any(word in upper for word in ["MOTIF", "CHROMVAR", "TF_TARGET", "TRANSCRIPTION_FACTOR"]):
        return 1.20
    return 1.0


def _term_gene_token(term: str) -> str | None:
    token = term.split()[0].upper() if term.split() else ""
    return token if GENE_RE.match(token) else None


def _term_parts(term_entry) -> tuple[str, set[str], float]:
    if len(term_entry) >= 3:
        return term_entry[0], term_entry[1], float(term_entry[2])
    return term_entry[0], term_entry[1], 1.0


def select_prior_terms(prior_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 180) -> list[tuple[str, set[str], float]]:
    selected = []
    for path in sorted(prior_dir.glob("*.gmt")):
        scored = []
        include_direct = _prior_should_include_term_gene(path.stem)
        for term, genes in parse_gmt(path, include_term_gene=include_direct):
            overlap = len(genes & perturb_genes)
            if overlap == 0 or len(genes) < 5 or len(genes) > 800:
                continue
            direct_gene = _term_gene_token(term)
            direct_hit = bool(direct_gene and direct_gene in perturb_genes)
            coverage = overlap / max(1, len(genes))
            weight = (
                _library_weight(path.stem)
                * _mechanism_weight(term)
                * (1.75 if direct_hit and include_direct else 1.0)
                * (1.0 + min(1.0, 8.0 * coverage))
            )
            score = (direct_hit, overlap, coverage, -len(genes))
            scored.append((score, f"{path.stem}:{term}", genes, float(weight)))
        scored.sort(reverse=True, key=lambda item: item[0])
        selected.extend((name, genes, weight) for _, name, genes, weight in scored[:max_terms_per_library])
    if not selected:
        selected = [("ko_gene_count", set(), 1.0)]
    return selected


def ko_prior_vector(label: str, terms: list[tuple[str, set[str]]]) -> np.ndarray:
    genes = set(split_ko(label))
    denom = max(1, len(genes))
    values = []
    for term_entry in terms:
        term_name, members, weight = _term_parts(term_entry)
        overlap = len(genes & members)
        overlap_fraction = overlap / denom
        term_coverage = overlap / max(1, len(members))
        direct_gene = _term_gene_token(term_name.split(":", 1)[-1])
        direct_hit = float(bool(direct_gene and direct_gene in genes))
        values.append(weight * (overlap_fraction + 0.35 * direct_hit + 0.15 * term_coverage))
    values.append(float(len(genes)))
    return np.asarray(values, dtype=float)


def fit_pls(train_labels: list[str], train_delta: np.ndarray, terms: list[tuple[str, set[str]]]):
    x = np.vstack([ko_prior_vector(label, terms) for label in train_labels])
    n_components = min(6, x.shape[0] - 1, x.shape[1], train_delta.shape[1])
    return make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True)).fit(x, train_delta)


def loo_training_predictions(train_labels: list[str], train_delta: np.ndarray, terms: list[tuple[str, set[str]]]) -> np.ndarray:
    preds = []
    for i, ko in enumerate(train_labels):
        keep = np.arange(len(train_labels)) != i
        if keep.sum() < 3:
            preds.append(train_delta[i])
            continue
        kept_labels = [label for j, label in enumerate(train_labels) if keep[j]]
        model = fit_pls(kept_labels, train_delta[keep], terms)
        preds.append(model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1))
    return np.vstack(preds)


def calibration_factors(loo_pred: np.ndarray, truth: np.ndarray, mode: str) -> tuple[np.ndarray, str]:
    if mode == "none":
        return np.ones(truth.shape[1]), "none"
    eps = 1e-9
    global_alpha = np.sum(loo_pred * truth) / (np.sum(loo_pred * loo_pred) + eps)
    global_alpha = float(np.clip(global_alpha, 0.15, 2.5))
    feature_alpha = np.sum(loo_pred * truth, axis=0) / (np.sum(loo_pred * loo_pred, axis=0) + eps)
    feature_alpha = np.clip(feature_alpha, 0.15, 2.5)
    candidates = {
        "none": np.ones(truth.shape[1]),
        "global_scale": np.full(truth.shape[1], global_alpha),
        "feature_scale": feature_alpha,
    }
    if mode in candidates and mode != "auto":
        return candidates[mode], mode
    losses = {name: float(np.mean(np.abs(loo_pred * alpha.reshape(1, -1) - truth))) for name, alpha in candidates.items()}
    best = min(losses, key=losses.get)
    return candidates[best], best


def shape_profile(frame: pd.DataFrame, ko_col: str, features: list[str], holdouts: set[str], mode: str) -> dict:
    profile: dict = {"mode": mode}
    if mode == "none":
        profile["shape_alpha"] = np.ones(len(features), dtype=float)
        return profile
    ctrl = frame.loc[control_mask(frame[ko_col]), features].to_numpy(dtype=float)
    ctrl_mean = np.nanmean(ctrl, axis=0)
    ctrl_std = np.nanstd(ctrl, axis=0)
    ctrl_std = np.where(ctrl_std < 1e-6, 1.0, ctrl_std)
    eps = 1e-9
    ctrl_open_fraction = np.nanmean(ctrl > eps, axis=0)
    ctrl_min = np.nanmin(ctrl, axis=0)
    ratios = []
    open_fractions = []
    centered_quantiles = []
    grid = np.asarray([0.00, 0.01, 0.03, 0.05, 0.10, 0.20, 0.35, 0.50, 0.65, 0.80, 0.90, 0.95, 0.97, 0.99, 1.00])
    for ko, group in frame.loc[~control_mask(frame[ko_col])].groupby(ko_col, observed=True):
        ko = str(ko)
        if ko in holdouts or len(group) < 3:
            continue
        values = group[features].to_numpy(dtype=float)
        std = np.nanstd(values, axis=0)
        ratios.append(std / ctrl_std)
        open_fractions.append(np.nanmean(values > eps, axis=0))
        centered = values - np.nanmean(values, axis=0).reshape(1, -1)
        centered_quantiles.append(np.nanquantile(centered, grid, axis=0))
    if not ratios:
        profile["shape_alpha"] = np.ones(len(features), dtype=float)
        return profile
    ratio = np.nanmedian(np.vstack(ratios), axis=0)
    profile["shape_alpha"] = np.clip(np.nan_to_num(ratio, nan=1.0, posinf=1.0, neginf=1.0), 0.45, 2.20)
    if mode == "quantile" and centered_quantiles:
        target_quantiles = np.nanmedian(np.stack(centered_quantiles, axis=0), axis=0)
        target_quantiles = np.nan_to_num(target_quantiles, nan=0.0, posinf=0.0, neginf=0.0)
        # Keep the learned quantile shape anchored to the training KO variance scale but bounded.
        target_quantiles = np.clip(target_quantiles, -6.0, 6.0)
        target_open_fraction = np.nanmedian(np.vstack(open_fractions), axis=0)
        target_open_fraction = np.clip(np.nan_to_num(target_open_fraction, nan=ctrl_open_fraction), 0.0, 1.0)
        sparse_mask = (ctrl_min >= -eps) & ((ctrl_open_fraction < 0.98) | (target_open_fraction < 0.98))
        hard_zero_mask = sparse_mask.copy()
        profile["quantile_grid"] = grid
        profile["target_centered_quantiles"] = target_quantiles
        profile["target_open_fraction"] = target_open_fraction
        profile["control_open_fraction"] = ctrl_open_fraction
        profile["sparse_nonnegative_mask"] = sparse_mask
        profile["hard_zero_mask"] = hard_zero_mask
        profile["control_mean"] = ctrl_mean
    return profile


def apply_shape_calibration(ctrl: np.ndarray, control_mean: np.ndarray, pred_delta: np.ndarray, profile: dict) -> np.ndarray:
    mode = profile.get("mode", "none")
    if mode == "quantile" and "target_centered_quantiles" in profile:
        grid = np.asarray(profile["quantile_grid"], dtype=float)
        target = np.asarray(profile["target_centered_quantiles"], dtype=float)
        shaped = np.zeros_like(ctrl, dtype=float)
        n = ctrl.shape[0]
        if n <= 1:
            return control_mean.reshape(1, -1) + pred_delta.reshape(1, -1)
        target_open = np.asarray(profile.get("target_open_fraction", np.ones(ctrl.shape[1])), dtype=float)
        sparse_mask = np.asarray(profile.get("sparse_nonnegative_mask", np.zeros(ctrl.shape[1], dtype=bool)), dtype=bool)
        hard_zero_mask = np.asarray(profile.get("hard_zero_mask", sparse_mask), dtype=bool)
        for j in range(ctrl.shape[1]):
            order = np.argsort(ctrl[:, j], kind="mergesort")
            pct = np.empty(n, dtype=float)
            pct[order] = np.linspace(0.0, 1.0, n)
            mapped = np.interp(pct, grid, target[:, j])
            if sparse_mask[j]:
                open_fraction = float(np.clip(target_open[j], 0.0, 1.0))
                if open_fraction <= 0.0:
                    mapped = np.full(n, np.interp(0.0, grid, target[:, j]), dtype=float)
                elif open_fraction < 0.98:
                    open_start = 1.0 - open_fraction
                    closed_level = float(np.interp(open_start, grid, target[:, j]))
                    mapped = np.where(pct < open_start, closed_level, mapped)
            shaped[:, j] = mapped
        virtual = control_mean.reshape(1, -1) + shaped + pred_delta.reshape(1, -1)
        for j in np.where(hard_zero_mask)[0]:
            open_fraction = float(np.clip(target_open[j], 0.0, 1.0))
            n_open = int(round(open_fraction * n))
            if n_open <= 0:
                virtual[:, j] = 0.0
                continue
            if n_open >= n:
                virtual[:, j] = np.maximum(virtual[:, j], 0.0)
                continue
            desired_mean = float(control_mean[j] + pred_delta[j])
            order = np.argsort(virtual[:, j], kind="mergesort")
            closed = order[: n - n_open]
            opened = order[n - n_open :]
            virtual[closed, j] = 0.0
            needed = desired_mean * n - float(np.sum(virtual[:, j]))
            virtual[opened, j] = np.maximum(virtual[opened, j] + needed / max(1, len(opened)), 0.0)
        return virtual
    shape_alpha = np.asarray(profile.get("shape_alpha", np.ones(ctrl.shape[1])), dtype=float)
    return control_mean.reshape(1, -1) + (ctrl - control_mean.reshape(1, -1)) * shape_alpha.reshape(1, -1) + pred_delta.reshape(1, -1)


def feature_columns(frame: pd.DataFrame, ko_col: str, explicit: list[str] | None = None) -> list[str]:
    if explicit:
        return explicit
    ignored = {ko_col, "cell_id", "dataset", "state", "calibration_method"}
    return [
        col
        for col in frame.columns
        if col not in ignored and pd.api.types.is_numeric_dtype(frame[col]) and frame[col].notna().any()
    ]


def build_delta_table(frame: pd.DataFrame, ko_col: str, features: list[str], holdouts: set[str]) -> tuple[list[str], np.ndarray, np.ndarray]:
    ctrl = frame.loc[control_mask(frame[ko_col]), features]
    if ctrl.empty:
        raise ValueError("No control cells detected. Use labels like control, ctrl, NT, negative, or non-targeting.")
    control_mean = ctrl.mean().to_numpy(dtype=float)
    labels, deltas = [], []
    for ko, group in frame.loc[~control_mask(frame[ko_col])].groupby(ko_col, observed=True):
        ko = str(ko)
        if ko in holdouts:
            continue
        labels.append(ko)
        deltas.append(group[features].mean().to_numpy(dtype=float) - control_mean)
    if len(labels) < 3:
        raise ValueError("Need at least 3 training KO labels after excluding holdouts.")
    return labels, np.vstack(deltas), control_mean


def make_auc_points(delta_table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for _, row in delta_table.iterrows():
        for feature in features:
            rows.append(
                {
                    "ko_target": row["ko_target"],
                    "feature": feature,
                    "true_abs_delta": abs(float(row[f"true_delta_{feature}"])),
                    "pred_abs_delta": abs(float(row[f"pred_delta_{feature}"])),
                    "true_delta": float(row[f"true_delta_{feature}"]),
                    "pred_delta": float(row[f"pred_delta_{feature}"]),
                }
            )
    return pd.DataFrame(rows)


def run_virtual_ko(
    frame: pd.DataFrame,
    ko_col: str,
    holdouts: list[str],
    prior_dir: str | Path,
    features: list[str] | None = None,
    dataset_name: str = "Virtual KO dataset",
    modality: str = "state score table",
    representation: str = "pathway/program scores",
    calibrate: str = "auto",
    shape_calibrate: str = "none",
    max_cells_per_state: int = 180,
    seed: int = 7,
) -> VirtualKOResult:
    rng = np.random.default_rng(seed)
    features = feature_columns(frame, ko_col, features)
    holdout_set = set(holdouts)
    perturb_genes = {gene for ko in frame[ko_col].astype(str).unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path(prior_dir), perturb_genes)
    train_labels, train_delta, _ = build_delta_table(frame, ko_col, features, holdout_set)
    model = fit_pls(train_labels, train_delta, terms)
    loo_pred = loo_training_predictions(train_labels, train_delta, terms)
    alpha, calibration_method = calibration_factors(loo_pred, train_delta, calibrate)
    shape = shape_profile(frame, ko_col, features, holdout_set, shape_calibrate)
    shape_alpha = np.asarray(shape.get("shape_alpha", np.ones(len(features))), dtype=float)
    target_open_fraction = np.asarray(shape.get("target_open_fraction", np.full(len(features), np.nan)), dtype=float)
    control_open_fraction = np.asarray(shape.get("control_open_fraction", np.full(len(features), np.nan)), dtype=float)
    sparse_nonnegative = np.asarray(shape.get("sparse_nonnegative_mask", np.zeros(len(features), dtype=bool)), dtype=bool)
    hard_zero_shape = np.asarray(shape.get("hard_zero_mask", np.zeros(len(features), dtype=bool)), dtype=bool)

    control = frame.loc[control_mask(frame[ko_col]), features].to_numpy(dtype=float)
    control_mean = np.nanmean(control, axis=0)
    metric_rows, cell_rows, delta_rows = [], [], []
    for ko in holdouts:
        true = frame.loc[frame[ko_col].astype(str) == ko, features].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        ctrl = control[rng.integers(0, len(control), size=len(true))]
        pred_delta = model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1) * alpha
        virtual = apply_shape_calibration(ctrl, control_mean, pred_delta, shape)
        true_delta = true.mean(axis=0) - ctrl.mean(axis=0)
        denom = np.linalg.norm(true_delta) * np.linalg.norm(pred_delta)
        delta_row = {
            "dataset": dataset_name,
            "input_modality": modality,
            "state_representation": representation,
            "ko_target": ko,
            "calibration_method": calibration_method,
            "shape_calibration_method": shape_calibrate,
            "direction_cosine": float(np.dot(true_delta, pred_delta) / denom) if denom > 1e-9 else np.nan,
            "mean_abs_delta_error": float(np.mean(np.abs(pred_delta - true_delta))),
            "true_delta_norm": float(np.linalg.norm(true_delta)),
            "pred_delta_norm": float(np.linalg.norm(pred_delta)),
        }
        for feature, truth, pred in zip(features, true_delta, pred_delta):
            delta_row[f"true_delta_{feature}"] = truth
            delta_row[f"pred_delta_{feature}"] = pred
        delta_rows.append(delta_row)

        for j, feature in enumerate(features):
            w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
            w_virt = wasserstein_distance(true[:, j], virtual[:, j])
            metric_rows.append(
                {
                    "dataset": dataset_name,
                    "input_modality": modality,
                    "state_representation": representation,
                    "ko_target": ko,
                    "feature": feature,
                    "wasserstein_true_vs_virtual": w_virt,
                    "wasserstein_true_vs_control": w_ctrl,
                    "distribution_improvement": 1.0 - w_virt / w_ctrl if w_ctrl > 1e-9 else np.nan,
                }
            )
        for state, matrix in [("control cells", ctrl), ("virtual KO cells", virtual), ("true KO cells", true)]:
            take = min(max_cells_per_state, len(matrix))
            idx = rng.choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=features)
            tmp["dataset"] = dataset_name
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)

    metrics = pd.DataFrame(metric_rows)
    virtual_cells = pd.concat(cell_rows, ignore_index=True) if cell_rows else pd.DataFrame()
    delta_table = pd.DataFrame(delta_rows)
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby(["dataset", "input_modality", "state_representation"], observed=True)
        .agg(
            n_ko=("ko_target", "nunique"),
            n_features=("feature", "nunique"),
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
    )
    if not delta_table.empty:
        summary["mean_direction_cosine"] = float(delta_table["direction_cosine"].mean())
        summary["mean_abs_delta_error"] = float(delta_table["mean_abs_delta_error"].mean())
        summary["calibration_method"] = calibration_method
        summary["shape_calibration_method"] = shape_calibrate
    auc_points = make_auc_points(delta_table, features) if not delta_table.empty else pd.DataFrame()
    calibration = pd.DataFrame(
        {
            "feature": features,
            "alpha": alpha,
            "calibration_method": calibration_method,
            "shape_alpha": shape_alpha,
            "control_open_fraction": control_open_fraction,
            "target_open_fraction": target_open_fraction,
            "sparse_nonnegative_shape_feature": sparse_nonnegative,
            "hard_zero_shape_feature": hard_zero_shape,
            "shape_calibration_method": shape_calibrate,
        }
    )
    return VirtualKOResult(metrics, summary, virtual_cells, delta_table, auc_points, calibration)
