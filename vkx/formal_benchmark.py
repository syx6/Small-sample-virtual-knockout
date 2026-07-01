from __future__ import annotations

from pathlib import Path
import json

import numpy as np
import pandas as pd


EXTERNAL_METHODS = ["scGen", "CPA", "GEARS", "CellOT"]


def run_formal_benchmark(
    state_csv: str | Path,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    out_dir: str | Path,
    features: list[str] | None = None,
    methods: list[str] | None = None,
    external_predictions_csv: str | Path | None = None,
    calibrate: str = "auto",
    shape_calibrate: str = "none",
    seed: int = 7,
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(state_csv)
    features = _feature_columns(frame, ko_col, features)
    methods = methods or ["vkx", "boosted", "ensemble", "pls", "ridge", "additive", "scgen", "cpa", "gears", "cellot"]
    methods = [_canonical_method(method) for method in methods]

    predictions = []
    availability = []
    if "vkx" in methods:
        pred, status = _predict_vkx(frame, ko_col, target_kos, prior_dir, features, calibrate, shape_calibrate, seed)
        predictions.append(pred)
        availability.append(status)
    if "pls" in methods:
        pred, status = _predict_prior_model(frame, ko_col, target_kos, prior_dir, features, model_type="pls")
        predictions.append(pred)
        availability.append(status)
    if "ridge" in methods:
        pred, status = _predict_prior_model(frame, ko_col, target_kos, prior_dir, features, model_type="ridge")
        predictions.append(pred)
        availability.append(status)
    if "ensemble" in methods:
        pred, status = _predict_constrained_ensemble(frame, ko_col, target_kos, prior_dir, features, seed=seed)
        predictions.append(pred)
        availability.append(status)
    if "calibrated" in methods:
        pred, status = _predict_calibrated_ensemble(frame, ko_col, target_kos, prior_dir, features, seed=seed)
        predictions.append(pred)
        availability.append(status)
    if "boosted" in methods:
        pred, status = _predict_response_boosted_anchor(frame, ko_col, target_kos, prior_dir, features, seed=seed)
        predictions.append(pred)
        availability.append(status)
    if "additive" in methods:
        pred, status = _predict_additive(frame, ko_col, target_kos, features)
        predictions.append(pred)
        availability.append(status)

    external = _load_external_predictions(external_predictions_csv, features)
    for method in EXTERNAL_METHODS:
        key = _canonical_method(method)
        if key not in methods:
            continue
        if not external.empty and method in set(external["method"]):
            pred = external.loc[external["method"] == method].copy()
            predictions.append(pred)
            availability.append({"method": method, "status": "provided", "reason": "Loaded from external prediction CSV."})
        else:
            availability.append(
                {
                    "method": method,
                    "status": "not_run",
                    "reason": "No external prediction CSV was provided. Install/run the method separately and pass its predictions for a strict benchmark.",
                }
            )

    prediction_table = pd.concat([p for p in predictions if not p.empty], ignore_index=True) if predictions else pd.DataFrame()
    true_table = _true_delta_table(frame, ko_col, target_kos, features)
    metrics = _score_predictions(prediction_table, true_table, features)
    roc_points = _benchmark_roc_points(prediction_table, true_table, features)
    availability_table = pd.DataFrame(availability)
    prediction_table.to_csv(out / "formal_benchmark_predictions.csv", index=False)
    true_table.to_csv(out / "formal_benchmark_truth.csv", index=False)
    metrics.to_csv(out / "formal_benchmark_metrics.csv", index=False)
    roc_points.to_csv(out / "formal_benchmark_roc_points.csv", index=False)
    _method_metric_comparison(metrics).to_csv(out / "method_metric_comparison.csv", index=False)
    availability_table.to_csv(out / "method_availability.csv", index=False)
    _plot_benchmark_figures(metrics, prediction_table, true_table, availability_table, roc_points, out, features)
    _write_report(metrics, availability_table, out)
    return {
        "predictions": prediction_table,
        "truth": true_table,
        "metrics": metrics,
        "roc_points": roc_points,
        "availability": availability_table,
    }


def _canonical_method(method: str) -> str:
    text = method.strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "vkx": "vkx",
        "pls": "pls",
        "priorpls": "pls",
        "ridge": "ridge",
        "ensemble": "ensemble",
        "constrainedensemble": "ensemble",
        "vkxensemble": "ensemble",
        "calibrated": "calibrated",
        "calibratedensemble": "calibrated",
        "amplitudecalibrated": "calibrated",
        "boosted": "boosted",
        "responseboosted": "boosted",
        "priorboosted": "boosted",
        "additive": "additive",
        "scgen": "scgen",
        "cpa": "cpa",
        "gears": "gears",
        "cellot": "cellot",
    }
    if text not in aliases:
        raise ValueError(f"Unknown benchmark method: {method}")
    return aliases[text]


def _feature_columns(frame: pd.DataFrame, ko_col: str, explicit: list[str] | None) -> list[str]:
    if explicit:
        return explicit
    ignored = {ko_col, "cell_id", "dataset", "state", "batch", "cell_type"}
    return [col for col in frame.columns if col not in ignored and pd.api.types.is_numeric_dtype(frame[col])]


def _predict_vkx(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    calibrate: str,
    shape_calibrate: str,
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    try:
        from .core import run_virtual_ko

        result = run_virtual_ko(
            frame=frame,
            ko_col=ko_col,
            holdouts=target_kos,
            prior_dir=prior_dir,
            features=features,
            dataset_name="formal benchmark",
            modality="state score table",
            representation="state scores",
            calibrate=calibrate,
            shape_calibrate=shape_calibrate,
            seed=seed,
        )
        rows = []
        for _, row in result.delta_table.iterrows():
            out = {"method": "VKX", "ko_target": row["ko_target"], "prediction_status": "ok"}
            for feature in features:
                out[f"pred_delta_{feature}"] = row.get(f"pred_delta_{feature}", np.nan)
            rows.append(out)
        return pd.DataFrame(rows), {"method": "VKX", "status": "run", "reason": "Native VKX prior-constrained residual/PLS baseline."}
    except Exception as exc:
        return pd.DataFrame(), {"method": "VKX", "status": "failed", "reason": str(exc)}


def _predict_prior_model(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    model_type: str,
) -> tuple[pd.DataFrame, dict]:
    try:
        from .core import build_delta_table, fit_pls, ko_prior_vector, select_prior_terms, split_ko
    except Exception as exc:
        return pd.DataFrame(), {"method": model_type.upper(), "status": "failed", "reason": f"Missing core dependency: {exc}"}

    holdout_set = set(target_kos)
    perturb_genes = {gene for ko in frame[ko_col].astype(str).unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path(prior_dir), perturb_genes)
    labels, train_delta, _ = build_delta_table(frame, ko_col, features, holdout_set)
    x_train = np.vstack([ko_prior_vector(label, terms) for label in labels])
    if model_type == "pls":
        model = fit_pls(labels, train_delta, terms)
        method_name = "PLS"
    else:
        model = _fit_numpy_ridge(x_train, train_delta)
        method_name = "Ridge"
    rows = []
    for ko in target_kos:
        x = ko_prior_vector(ko, terms).reshape(1, -1)
        pred = model.predict(x).reshape(-1) if model_type == "pls" else _predict_numpy_ridge(model, x).reshape(-1)
        out = {"method": method_name, "ko_target": ko, "prediction_status": "ok"}
        for feature, value in zip(features, pred):
            out[f"pred_delta_{feature}"] = value
        rows.append(out)
    return pd.DataFrame(rows), {"method": method_name, "status": "run", "reason": f"{method_name} prior-vector delta regression."}


def _predict_constrained_ensemble(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    del seed
    try:
        from .core import build_delta_table, fit_pls, ko_prior_vector, select_prior_terms, split_ko
    except Exception as exc:
        return pd.DataFrame(), {"method": "ConstrainedEnsemble", "status": "failed", "reason": f"Missing core dependency: {exc}"}

    holdout_set = set(target_kos)
    perturb_genes = {gene for ko in frame[ko_col].astype(str).unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path(prior_dir), perturb_genes)
    labels, train_delta, _ = build_delta_table(frame, ko_col, features, holdout_set)
    x_train = np.vstack([ko_prior_vector(label, terms) for label in labels])
    pls_model = fit_pls(labels, train_delta, terms)
    ridge_model = _fit_numpy_ridge(x_train, train_delta)
    pls_train_pred = pls_model.predict(x_train)
    ridge_train_pred = _predict_numpy_ridge(ridge_model, x_train)
    pls_mae = float(np.mean(np.abs(pls_train_pred - train_delta)))
    ridge_mae = float(np.mean(np.abs(ridge_train_pred - train_delta)))
    weights = _inverse_error_weights({"PLS": pls_mae, "Ridge": ridge_mae})
    rows = []
    for ko in target_kos:
        x = ko_prior_vector(ko, terms).reshape(1, -1)
        pls_pred = pls_model.predict(x).reshape(-1)
        ridge_pred = _predict_numpy_ridge(ridge_model, x).reshape(-1)
        pred = weights["PLS"] * pls_pred + weights["Ridge"] * ridge_pred
        out = {
            "method": "ConstrainedEnsemble",
            "ko_target": ko,
            "prediction_status": "ok",
            "ensemble_weight_pls": weights["PLS"],
            "ensemble_weight_ridge": weights["Ridge"],
        }
        for feature, value in zip(features, pred):
            out[f"pred_delta_{feature}"] = value
        rows.append(out)
    reason = f"Training-error weighted PLS/Ridge anchor. weights: PLS={weights['PLS']:.2f}, Ridge={weights['Ridge']:.2f}."
    return pd.DataFrame(rows), {"method": "ConstrainedEnsemble", "status": "run", "reason": reason}


def _inverse_error_weights(errors: dict[str, float]) -> dict[str, float]:
    inv = {name: 1.0 / max(value, 1e-6) for name, value in errors.items() if np.isfinite(value)}
    if not inv:
        return {name: 1.0 / len(errors) for name in errors}
    total = sum(inv.values())
    return {name: inv.get(name, 0.0) / total for name in errors}


def _predict_calibrated_ensemble(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    try:
        from .core import build_delta_table, fit_pls, ko_prior_vector, select_prior_terms, split_ko
    except Exception as exc:
        return pd.DataFrame(), {"method": "CalibratedEnsemble", "status": "failed", "reason": f"Missing core dependency: {exc}"}

    holdout_set = set(target_kos)
    perturb_genes = {gene for ko in frame[ko_col].astype(str).unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path(prior_dir), perturb_genes)
    labels, train_delta, _ = build_delta_table(frame, ko_col, features, holdout_set)
    if len(labels) < 4:
        return pd.DataFrame(), {"method": "CalibratedEnsemble", "status": "failed", "reason": "Need at least four training KOs for calibrated ensemble."}

    x_train = np.vstack([ko_prior_vector(label, terms) for label in labels])
    oof = _kfold_oof_prior_predictions(labels, x_train, train_delta, terms, seed=seed)
    calibration = {}
    errors = {}
    for name, pred in oof.items():
        beta = _global_amplitude_scale(pred, train_delta)
        calibrated_pred = beta * pred
        calibration[name] = beta
        errors[name] = float(np.mean(np.abs(calibrated_pred - train_delta)))
    weights = _inverse_error_weights(errors)

    pls_model = fit_pls(labels, train_delta, terms)
    ridge_model = _fit_numpy_ridge(x_train, train_delta)
    rows = []
    for ko in target_kos:
        x = ko_prior_vector(ko, terms).reshape(1, -1)
        pls_pred = calibration["PLS"] * pls_model.predict(x).reshape(-1)
        ridge_pred = calibration["Ridge"] * _predict_numpy_ridge(ridge_model, x).reshape(-1)
        pred = weights["PLS"] * pls_pred + weights["Ridge"] * ridge_pred
        out = {
            "method": "CalibratedEnsemble",
            "ko_target": ko,
            "prediction_status": "ok",
            "ensemble_weight_pls": weights["PLS"],
            "ensemble_weight_ridge": weights["Ridge"],
            "amplitude_scale_pls": calibration["PLS"],
            "amplitude_scale_ridge": calibration["Ridge"],
        }
        for feature, value in zip(features, pred):
            out[f"pred_delta_{feature}"] = value
        rows.append(out)
    reason = (
        "K-fold out-of-fold weighted PLS/Ridge anchor with global amplitude calibration. "
        f"weights: PLS={weights['PLS']:.2f}, Ridge={weights['Ridge']:.2f}; "
        f"scales: PLS={calibration['PLS']:.2f}, Ridge={calibration['Ridge']:.2f}."
    )
    return pd.DataFrame(rows), {"method": "CalibratedEnsemble", "status": "run", "reason": reason}


def _predict_response_boosted_anchor(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    seed: int,
) -> tuple[pd.DataFrame, dict]:
    from .core import split_ko

    pred, status = _predict_constrained_ensemble(frame, ko_col, target_kos, prior_dir, features, seed=seed)
    if pred.empty:
        return pred, {"method": "ResponseBoosted", "status": status.get("status", "failed"), "reason": status.get("reason", "Base ensemble failed.")}
    pred = pred.copy()
    pred["method"] = "ResponseBoosted"
    boosted_events = []
    for idx, row in pred.iterrows():
        ko_genes = split_ko(str(row["ko_target"]))
        boosted_count = 0
        max_boost = 1.0
        boosted_families = set()
        for feature in features:
            family = _response_boost_family(feature)
            factor = _response_boost_factor(ko_genes, family)
            col = f"pred_delta_{feature}"
            if col in pred.columns and factor != 1.0:
                pred.at[idx, col] = pred.at[idx, col] * factor
                pred.at[idx, f"boost_factor_{feature}"] = factor
                pred.at[idx, f"boost_family_{feature}"] = family
                boosted_count += 1
                max_boost = max(max_boost, factor)
                boosted_families.add(family)
        pred.at[idx, "boosted_feature_count"] = boosted_count
        pred.at[idx, "max_boost_factor"] = max_boost
        pred.at[idx, "boosted_families"] = ";".join(sorted(boosted_families)) if boosted_families else "none"
        boosted_events.append(boosted_count)
    reason = (
        "Constrained ensemble plus adaptive response-strength priors. "
        "Boost factors are selected from KO genes and feature families and are written into prediction metadata. "
        f"mean_boosted_features={np.mean(boosted_events):.1f}."
    )
    return pred, {"method": "ResponseBoosted", "status": "run", "reason": reason}


def _response_boost_family(feature: str) -> str:
    text = str(feature).lower()
    if any(key in text for key in ["interferon", "jak", "stat", "ifn"]):
        return "ifn_jak_stat"
    if any(key in text for key in ["mapk", "tgfb", "tgf_beta", "erk"]):
        return "mapk_tgfb"
    if any(key in text for key in ["myc", "e2f", "g2m", "cell_cycle", "mitotic"]):
        return "cell_cycle_myc_e2f"
    if text.startswith("protein_") and any(key in text for key in ["pdl1", "pdl2", "cd86", "cd366", "pdcd1", "ctla4"]):
        return "protein_checkpoint"
    return "other"


def _response_boost_factor(ko_genes: list[str], family: str) -> float:
    genes = {gene.upper() for gene in ko_genes}
    if family == "ifn_jak_stat" and any(gene.startswith(("STAT", "JAK", "IFN", "IRF")) for gene in genes):
        return 1.5
    if family == "mapk_tgfb" and any(gene.startswith(("MAPK", "RAF", "RAS", "MEK", "ERK", "TGFB", "SMAD")) for gene in genes):
        return 1.35
    if family == "cell_cycle_myc_e2f" and any(gene.startswith(("MYC", "E2F", "CDK", "CCN", "RB")) for gene in genes):
        return 1.25
    if family == "protein_checkpoint" and any(gene.startswith(("STAT", "JAK", "IFN", "IRF", "CD274", "PDCD1LG")) for gene in genes):
        return 1.15
    return 1.0


def _kfold_oof_prior_predictions(
    labels: list[str],
    x_train: np.ndarray,
    train_delta: np.ndarray,
    terms: list[str],
    seed: int,
    max_folds: int = 3,
) -> dict[str, np.ndarray]:
    from .core import fit_pls

    n = len(labels)
    rng = np.random.default_rng(seed)
    order = rng.permutation(n)
    folds = np.array_split(order, min(max_folds, n))
    out = {
        "PLS": np.zeros_like(train_delta, dtype=float),
        "Ridge": np.zeros_like(train_delta, dtype=float),
    }
    for fold in folds:
        keep = np.ones(n, dtype=bool)
        keep[fold] = False
        fold_labels = [labels[i] for i in np.where(keep)[0]]
        pls = fit_pls(fold_labels, train_delta[keep], terms)
        ridge = _fit_numpy_ridge(x_train[keep], train_delta[keep])
        out["PLS"][fold] = pls.predict(x_train[fold])
        out["Ridge"][fold] = _predict_numpy_ridge(ridge, x_train[fold])
    return out


def _global_amplitude_scale(pred: np.ndarray, truth: np.ndarray, min_scale: float = 0.5, max_scale: float = 3.0) -> float:
    mask = np.isfinite(pred) & np.isfinite(truth)
    if mask.sum() == 0:
        return 1.0
    p = pred[mask]
    t = truth[mask]
    denom = float(np.dot(p, p))
    if denom <= 1e-9:
        return 1.0
    scale = float(np.dot(p, t) / denom)
    if not np.isfinite(scale):
        return 1.0
    return float(np.clip(scale, min_scale, max_scale))


def _fit_numpy_ridge(x_train: np.ndarray, y_train: np.ndarray) -> dict:
    x_mean = x_train.mean(axis=0, keepdims=True)
    x_std = x_train.std(axis=0, keepdims=True)
    x_std = np.where(x_std < 1e-8, 1.0, x_std)
    y_mean = y_train.mean(axis=0, keepdims=True)
    x = (x_train - x_mean) / x_std
    y = y_train - y_mean
    alphas = np.logspace(-3, 3, 13)
    best_alpha = alphas[0]
    best_loss = np.inf
    if len(x) >= 4:
        for alpha in alphas:
            losses = []
            for i in range(len(x)):
                keep = np.arange(len(x)) != i
                coef = _ridge_coef(x[keep], y[keep], alpha)
                pred = x[i : i + 1] @ coef + y_mean
                losses.append(np.mean(np.abs(pred - y_train[i : i + 1])))
            loss = float(np.mean(losses))
            if loss < best_loss:
                best_loss = loss
                best_alpha = alpha
    coef = _ridge_coef(x, y, best_alpha)
    return {"x_mean": x_mean, "x_std": x_std, "y_mean": y_mean, "coef": coef, "alpha": best_alpha}


def _ridge_coef(x: np.ndarray, y: np.ndarray, alpha: float) -> np.ndarray:
    reg = alpha * np.eye(x.shape[1])
    return np.linalg.pinv(x.T @ x + reg) @ x.T @ y


def _predict_numpy_ridge(model: dict, x: np.ndarray) -> np.ndarray:
    x_scaled = (x - model["x_mean"]) / model["x_std"]
    return x_scaled @ model["coef"] + model["y_mean"]


def _predict_additive(frame: pd.DataFrame, ko_col: str, target_kos: list[str], features: list[str]) -> tuple[pd.DataFrame, dict]:
    from .core import control_mask, split_ko

    control = frame.loc[control_mask(frame[ko_col]), features].mean().to_numpy(dtype=float)
    single_delta: dict[str, np.ndarray] = {}
    holdout_set = set(target_kos)
    for ko, group in frame.loc[~control_mask(frame[ko_col])].groupby(ko_col, observed=True):
        if str(ko) in holdout_set:
            continue
        genes = split_ko(str(ko))
        if len(genes) == 1:
            single_delta[genes[0]] = group[features].mean().to_numpy(dtype=float) - control
    rows = []
    missing = []
    for ko in target_kos:
        genes = split_ko(ko)
        if genes and all(gene in single_delta for gene in genes):
            pred = np.sum([single_delta[gene] for gene in genes], axis=0)
            status = "ok"
        else:
            pred = np.full(len(features), np.nan)
            status = "missing_single_gene_components"
            missing.append(ko)
        out = {"method": "Additive", "ko_target": ko, "prediction_status": status}
        for feature, value in zip(features, pred):
            out[f"pred_delta_{feature}"] = value
        rows.append(out)
    reason = "Single-gene delta summation for targets whose component single KOs are present."
    if missing:
        reason += f" Missing components for: {', '.join(missing[:8])}."
    return pd.DataFrame(rows), {"method": "Additive", "status": "run_partial" if missing else "run", "reason": reason}


def _true_delta_table(frame: pd.DataFrame, ko_col: str, target_kos: list[str], features: list[str]) -> pd.DataFrame:
    from .core import control_mask

    control = frame.loc[control_mask(frame[ko_col]), features].mean().to_numpy(dtype=float)
    rows = []
    for ko in target_kos:
        group = frame.loc[frame[ko_col].astype(str) == ko, features]
        if group.empty:
            continue
        delta = group.mean().to_numpy(dtype=float) - control
        out = {"ko_target": ko}
        for feature, value in zip(features, delta):
            out[f"true_delta_{feature}"] = value
        rows.append(out)
    return pd.DataFrame(rows)


def _load_external_predictions(path: str | Path | None, features: list[str]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    table = pd.read_csv(path)
    required = {"method", "ko_target"}
    if not required.issubset(table.columns):
        raise ValueError("External prediction CSV must contain method and ko_target columns.")
    for feature in features:
        col = f"pred_delta_{feature}"
        if col not in table.columns:
            table[col] = np.nan
    method_map = {"scgen": "scGen", "cpa": "CPA", "gears": "GEARS", "cellot": "CellOT"}
    table["method"] = table["method"].astype(str).map(lambda x: method_map.get(x.lower(), x))
    table["prediction_status"] = table.get("prediction_status", "external_provided")
    return table


def _score_predictions(pred: pd.DataFrame, truth: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    if pred.empty or truth.empty:
        return pd.DataFrame()
    merged = pred.merge(truth, on="ko_target", how="inner")
    for _, row in merged.iterrows():
        y_true = np.asarray([row.get(f"true_delta_{feature}", np.nan) for feature in features], dtype=float)
        y_pred = np.asarray([row.get(f"pred_delta_{feature}", np.nan) for feature in features], dtype=float)
        mask = np.isfinite(y_true) & np.isfinite(y_pred)
        if mask.sum() == 0:
            continue
        yt = y_true[mask]
        yp = y_pred[mask]
        rows.append(
            {
                "method": row["method"],
                "ko_target": row["ko_target"],
                "prediction_status": row.get("prediction_status", "ok"),
                "n_features": int(mask.sum()),
                "mae": float(np.mean(np.abs(yp - yt))),
                "r2": _r2(yt, yp),
                "direction_cosine": _cosine(yt, yp),
                "roc_auc": _auc(np.abs(yt), np.abs(yp)),
                "feature_hit_rate": _hit_rate(yt, yp),
            }
        )
    return pd.DataFrame(rows)


def _method_metric_comparison(metrics: pd.DataFrame) -> pd.DataFrame:
    if metrics.empty:
        return pd.DataFrame()
    return (
        metrics.groupby("method", observed=True)[["roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate"]]
        .mean(numeric_only=True)
        .reset_index()
        .sort_values(["roc_auc", "direction_cosine", "r2"], ascending=False)
    )


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    denom = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if denom <= 1e-12:
        return np.nan
    return float(1.0 - np.sum((y_true - y_pred) ** 2) / denom)


def _cosine(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    denom = float(np.linalg.norm(y_true) * np.linalg.norm(y_pred))
    return float(np.dot(y_true, y_pred) / denom) if denom > 1e-12 else np.nan


def _auc(abs_true: np.ndarray, abs_pred: np.ndarray) -> float:
    threshold = np.nanquantile(abs_true, 0.65)
    labels = abs_true >= threshold
    if labels.sum() == 0 or (~labels).sum() == 0:
        return np.nan
    order = np.argsort(abs_pred)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(abs_pred) + 1)
    pos_ranks = ranks[labels].sum()
    n_pos = labels.sum()
    n_neg = (~labels).sum()
    return float((pos_ranks - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def _hit_rate(y_true: np.ndarray, y_pred: np.ndarray, top_fraction: float = 0.25) -> float:
    k = max(1, int(round(len(y_true) * top_fraction)))
    true_top = set(np.argsort(np.abs(y_true))[-k:])
    pred_top = set(np.argsort(np.abs(y_pred))[-k:])
    return float(len(true_top & pred_top) / k)


def _benchmark_roc_points(pred: pd.DataFrame, truth: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    if pred.empty or truth.empty:
        return pd.DataFrame()
    merged = pred.merge(truth, on="ko_target", how="inner")
    for _, row in merged.iterrows():
        y_true = np.asarray([row.get(f"true_delta_{feature}", np.nan) for feature in features], dtype=float)
        y_pred = np.asarray([row.get(f"pred_delta_{feature}", np.nan) for feature in features], dtype=float)
        mask = np.isfinite(y_true) & np.isfinite(y_pred)
        if mask.sum() < 3:
            continue
        labels = np.abs(y_true[mask]) >= np.nanquantile(np.abs(y_true[mask]), 0.65)
        scores = np.abs(y_pred[mask])
        if labels.sum() == 0 or (~labels).sum() == 0:
            continue
        order = np.argsort(scores)[::-1]
        tp = 0
        fp = 0
        pos = int(labels.sum())
        neg = int((~labels).sum())
        method = row["method"]
        ko = row["ko_target"]
        rows.append({"method": method, "ko_target": ko, "fpr": 0.0, "tpr": 0.0, "threshold_rank": 0})
        for rank, idx in enumerate(order, start=1):
            if labels[idx]:
                tp += 1
            else:
                fp += 1
            rows.append({"method": method, "ko_target": ko, "fpr": fp / neg, "tpr": tp / pos, "threshold_rank": rank})
        rows.append({"method": method, "ko_target": ko, "fpr": 1.0, "tpr": 1.0, "threshold_rank": len(order) + 1})
    return pd.DataFrame(rows)


def _plot_benchmark_figures(
    metrics: pd.DataFrame,
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    availability: pd.DataFrame,
    roc_points: pd.DataFrame,
    out: Path,
    features: list[str],
) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        _plot_benchmark_fallback(metrics, pred, truth, availability, roc_points, out, features)
        return
    if not metrics.empty:
        summary = _method_metric_comparison(metrics)
        fig, axes = plt.subplots(1, 4, figsize=(17, 4.8), constrained_layout=True)
        for ax, metric, title in zip(
            axes,
            ["roc_auc", "direction_cosine", "r2", "mae"],
            ["AUC", "Direction", "R2", "MAE lower better"],
        ):
            order = summary.sort_values(metric, ascending=(metric == "mae"))["method"]
            sns.barplot(data=summary, x=metric, y="method", order=order, color="#4C78A8", ax=ax)
            ax.set_title(title)
            ax.set_ylabel("")
        fig.suptitle("Formal Virtual KO Benchmark")
        fig.savefig(out / "01_formal_benchmark_metric_panel.png", bbox_inches="tight", dpi=300)
        plt.close(fig)
    _plot_prediction_heatmap(pred, truth, out, features)
    if not availability.empty:
        fig, ax = plt.subplots(figsize=(9, max(3.8, 0.45 * len(availability) + 1.2)), constrained_layout=True)
        colors = availability["status"].map(lambda x: "#2A9D8F" if str(x).startswith("run") or x == "provided" else "#E76F51")
        ax.barh(availability["method"], np.ones(len(availability)), color=colors)
        ax.set_xlim(0, 1)
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.set_title("Method availability in this benchmark")
        for i, (_, row) in enumerate(availability.iterrows()):
            ax.text(0.03, i, f"{row['status']}: {row['reason']}", va="center", color="white" if str(row["status"]).startswith("run") else "black", fontsize=8)
        fig.savefig(out / "03_method_availability.png", bbox_inches="tight", dpi=300)
        plt.close(fig)
    _plot_benchmark_roc_curves(roc_points, metrics, out)


def _plot_benchmark_roc_curves(roc_points: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    if roc_points.empty:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    summary = _method_metric_comparison(metrics)
    auc_map = dict(zip(summary["method"], summary["roc_auc"])) if not summary.empty else {}
    colors = {
        "ConstrainedEnsemble": "#2A9D8F",
        "Ridge": "#4C78A8",
        "PLS": "#F4A261",
        "VKX": "#E76F51",
        "Additive": "#8D99AE",
        "scGen": "#7B2CBF",
        "CPA": "#4361EE",
        "GEARS": "#3A86FF",
        "CellOT": "#6A994E",
    }
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    for method, group in roc_points.groupby("method", observed=True):
        curve = group.groupby("fpr", as_index=False)["tpr"].mean().sort_values("fpr")
        label = f"{method} AUC={auc_map.get(method, np.nan):.2f}" if method in auc_map else str(method)
        ax.plot(curve["fpr"], curve["tpr"], marker="o", linewidth=2, markersize=3.5, color=colors.get(method, None), label=label)
    ax.plot([0, 1], [0, 1], linestyle="--", color="0.55", linewidth=1)
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC curves for strong-response feature ranking")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.savefig(out / "04_formal_benchmark_roc_curves.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_prediction_heatmap(pred: pd.DataFrame, truth: pd.DataFrame, out: Path, features: list[str], max_features: int = 12) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return
    if pred.empty or truth.empty:
        return
    merged = pred.merge(truth, on="ko_target", how="inner")
    if merged.empty:
        return
    feature_scores = []
    for feature in features:
        col = f"true_delta_{feature}"
        if col in merged.columns:
            feature_scores.append((feature, float(np.nanmean(np.abs(merged[col])))))
    chosen = [feature for feature, _ in sorted(feature_scores, key=lambda item: item[1], reverse=True)[:max_features]]
    rows = []
    for _, row in merged.iterrows():
        for feature in chosen:
            rows.append({"method_ko": f"{row['method']} | {row['ko_target']}", "feature": _short(feature), "delta": row.get(f"pred_delta_{feature}", np.nan), "kind": "pred"})
            rows.append({"method_ko": f"TRUE | {row['ko_target']}", "feature": _short(feature), "delta": row.get(f"true_delta_{feature}", np.nan), "kind": "true"})
    table = pd.DataFrame(rows).pivot_table(index="method_ko", columns="feature", values="delta", aggfunc="mean")
    vmax = np.nanmax(np.abs(table.to_numpy()))
    fig, ax = plt.subplots(figsize=(max(9, 0.65 * len(table.columns) + 3), max(5, 0.32 * len(table) + 2)), constrained_layout=True)
    sns.heatmap(table, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, ax=ax, cbar_kws={"label": "delta"})
    ax.set_title("True and Predicted KO Delta by Method")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(out / "02_formal_benchmark_delta_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_benchmark_fallback(
    metrics: pd.DataFrame,
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    availability: pd.DataFrame,
    roc_points: pd.DataFrame,
    out: Path,
    features: list[str],
) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    img = Image.new("RGB", (1450, 760), "white")
    draw = ImageDraw.Draw(img)
    font = _pil_font(16)
    small = _pil_font(13)
    title_font = _pil_font(26)
    draw.text((35, 28), "Formal Virtual KO Benchmark", fill=(20, 20, 20), font=title_font)
    if not metrics.empty:
        summary = _method_metric_comparison(metrics)
        best = summary.iloc[0]
        draw.text(
            (35, 67),
            f"Best ranked method: {best['method']} | AUC {best['roc_auc']:.2f}, direction {best['direction_cosine']:.2f}, R2 {best['r2']:.2f}, MAE {best['mae']:.3f}",
            fill=(70, 70, 70),
            font=font,
        )
        panels = [
            ("roc_auc", "AUC higher is better", 140, 0.0, 1.0),
            ("direction_cosine", "Direction higher is better", 465, 0.0, 1.0),
            ("r2", "R2 higher is better", 790, -1.0, 1.0),
            ("mae", "MAE lower is better", 1115, 0.0, max(0.01, float(summary["mae"].max()))),
        ]
        colors = {
            "ConstrainedEnsemble": (42, 157, 143),
            "VKX": (231, 111, 81),
            "PLS": (76, 120, 168),
            "Ridge": (42, 157, 143),
            "Additive": (158, 158, 158),
        }
        for metric, title, x0, vmin, vmax in panels:
            draw.text((x0, 125), title, fill=(30, 30, 30), font=font)
            for i, row in enumerate(summary.itertuples()):
                y = 170 + i * 66
                method = str(row.method)
                value = float(getattr(row, metric))
                if metric == "mae":
                    width = int(220 * (1.0 - min(max((value - vmin) / (vmax - vmin + 1e-9), 0), 1)))
                    label = f"{value:.3f}"
                else:
                    width = int(220 * min(max((value - vmin) / (vmax - vmin + 1e-9), 0), 1))
                    label = f"{value:.2f}"
                if x0 == 140:
                    draw.text((35, y + 8), _display_method(method), fill=(30, 30, 30), font=font)
                draw.rectangle((x0, y, x0 + 220, y + 28), outline=(220, 220, 220), fill=(246, 246, 246))
                draw.rectangle((x0, y, x0 + width, y + 28), fill=colors.get(method, (90, 90, 90)))
                draw.text((x0 + 232, y + 5), label, fill=(30, 30, 30), font=small)
        draw.text((35, 690), "Deep methods are included only when external prediction CSVs are provided. See 03_method_availability.png for run status.", fill=(80, 80, 80), font=small)
    img.save(out / "01_formal_benchmark_metric_panel.png")
    _plot_delta_heatmap_fallback(pred, truth, out, features)
    _plot_availability_fallback(availability, out)
    _plot_roc_fallback(roc_points, metrics, out)


def _plot_roc_fallback(roc_points: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    if roc_points.empty:
        return
    summary = _method_metric_comparison(metrics)
    auc_map = dict(zip(summary["method"], summary["roc_auc"])) if not summary.empty else {}
    colors = {
        "ConstrainedEnsemble": (42, 157, 143),
        "Ridge": (76, 120, 168),
        "PLS": (244, 162, 97),
        "VKX": (231, 111, 81),
        "Additive": (141, 153, 174),
    }
    img = Image.new("RGB", (1100, 820), "white")
    draw = ImageDraw.Draw(img)
    font = _pil_font(18)
    small = _pil_font(14)
    title_font = _pil_font(24)
    draw.text((45, 28), "ROC curves for strong-response feature ranking", fill=(20, 20, 20), font=title_font)
    left, top, size = 115, 105, 600
    draw.rectangle((left, top, left + size, top + size), outline=(180, 180, 180), width=1)
    draw.line((left, top + size, left + size, top), fill=(160, 160, 160), width=1)
    draw.text((left + 210, top + size + 42), "False positive rate", fill=(30, 30, 30), font=font)
    draw.text((25, top + 265), "True positive rate", fill=(30, 30, 30), font=font)
    for tick in [0.0, 0.5, 1.0]:
        x_tick = left + int(tick * size)
        y_tick = top + size - int(tick * size)
        draw.line((x_tick, top + size, x_tick, top + size + 6), fill=(80, 80, 80))
        draw.text((x_tick - 10, top + size + 10), f"{tick:.1f}", fill=(40, 40, 40), font=small)
        draw.line((left - 6, y_tick, left, y_tick), fill=(80, 80, 80))
        draw.text((left - 42, y_tick - 8), f"{tick:.1f}", fill=(40, 40, 40), font=small)
    legend_y = 115
    for method, group in roc_points.groupby("method", observed=True):
        curve = group.groupby("fpr", as_index=False)["tpr"].mean().sort_values("fpr")
        pts = []
        for row in curve.itertuples():
            x = left + int(float(row.fpr) * size)
            y = top + size - int(float(row.tpr) * size)
            pts.append((x, y))
        color = colors.get(str(method), (80, 80, 80))
        if len(pts) >= 2:
            draw.line(pts, fill=color, width=3)
        for x, y in pts:
            draw.ellipse((x - 3, y - 3, x + 3, y + 3), fill=color)
        draw.rectangle((760, legend_y, 784, legend_y + 16), fill=color)
        auc = auc_map.get(method, np.nan)
        label = f"{method}: AUC {auc:.2f}" if np.isfinite(auc) else str(method)
        draw.text((796, legend_y - 2), label, fill=(30, 30, 30), font=small)
        legend_y += 34
    img.save(out / "04_formal_benchmark_roc_curves.png")


def _plot_delta_heatmap_fallback(pred: pd.DataFrame, truth: pd.DataFrame, out: Path, features: list[str]) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    if pred.empty or truth.empty:
        return
    merged = pred.merge(truth, on="ko_target", how="inner")
    if merged.empty:
        return
    feature_scores = []
    for feature in features:
        values = merged.get(f"true_delta_{feature}")
        if values is not None:
            feature_scores.append((feature, float(np.nanmean(np.abs(values)))))
    chosen = [feature for feature, _ in sorted(feature_scores, key=lambda item: item[1], reverse=True)[:10]]
    rows = []
    for _, row in merged.iterrows():
        rows.append((f"{row['method']} | {row['ko_target']}", [float(row.get(f"pred_delta_{feature}", np.nan)) for feature in chosen]))
    for _, row in truth.iterrows():
        rows.append((f"TRUE | {row['ko_target']}", [float(row.get(f"true_delta_{feature}", np.nan)) for feature in chosen]))
    cell_w = 76
    cell_h = 22
    left = 260
    top = 100
    width = left + cell_w * len(chosen) + 40
    height = top + cell_h * len(rows) + 40
    img = Image.new("RGB", (max(900, width), max(420, height)), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    draw.text((30, 25), "True and Predicted KO Delta by Method", fill=(20, 20, 20), font=font)
    values = np.asarray([vals for _, vals in rows], dtype=float)
    vmax = max(float(np.nanmax(np.abs(values))), 1e-9)
    for j, feature in enumerate(chosen):
        draw.text((left + j * cell_w, 70), _short(feature, 10), fill=(30, 30, 30), font=font)
    for i, (label, vals) in enumerate(rows):
        y = top + i * cell_h
        draw.text((30, y), _short(label, 32), fill=(30, 30, 30), font=font)
        for j, value in enumerate(vals):
            x = left + j * cell_w
            if not np.isfinite(value):
                color = (235, 235, 235)
            elif value >= 0:
                intensity = int(255 - 160 * min(abs(value) / vmax, 1))
                color = (255, intensity, intensity)
            else:
                intensity = int(255 - 160 * min(abs(value) / vmax, 1))
                color = (intensity, intensity, 255)
            draw.rectangle((x, y, x + cell_w - 3, y + cell_h - 3), fill=color, outline=(245, 245, 245))
    img.save(out / "02_formal_benchmark_delta_heatmap.png")


def _plot_availability_fallback(availability: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    if availability.empty:
        return
    row_h = 38
    img = Image.new("RGB", (1300, 80 + row_h * len(availability)), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    draw.text((30, 25), "Method availability", fill=(20, 20, 20), font=font)
    for i, (_, row) in enumerate(availability.iterrows()):
        y = 70 + i * row_h
        status = str(row["status"])
        color = (42, 157, 143) if status.startswith("run") or status == "provided" else (231, 111, 81)
        draw.rectangle((30, y, 190, y + 24), fill=color)
        draw.text((40, y + 6), str(row["method"]), fill=(255, 255, 255), font=font)
        draw.text((210, y + 6), f"{status}: {row['reason']}", fill=(30, 30, 30), font=font)
    img.save(out / "03_method_availability.png")


def _write_report(metrics: pd.DataFrame, availability: pd.DataFrame, out: Path) -> None:
    summary = _method_metric_comparison(metrics)
    payload = {
        "n_methods_scored": int(summary.shape[0]) if not summary.empty else 0,
        "methods": summary["method"].tolist() if not summary.empty else [],
        "external_methods": EXTERNAL_METHODS,
    }
    (out / "formal_benchmark_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    text = [
        "# Formal Virtual KO Benchmark",
        "",
        "This benchmark compares VKX against classical baselines and reserves explicit slots for scGen, CPA, GEARS, and CellOT.",
        "",
        "Deep external methods are marked `not_run` unless their prediction CSV is provided. This prevents accidental over-claiming.",
        "",
        "## Metric Summary",
        "",
        summary.to_string(index=False) if not summary.empty else "No scored methods.",
        "",
        "## Method Availability",
        "",
        availability.to_string(index=False) if not availability.empty else "No method availability table.",
        "",
        "## Figures",
        "",
        "- `01_formal_benchmark_metric_panel.png`",
        "- `02_formal_benchmark_delta_heatmap.png`",
        "- `03_method_availability.png`",
        "- `04_formal_benchmark_roc_curves.png`",
    ]
    (out / "formal_benchmark_report.md").write_text("\n".join(text), encoding="utf-8")


def _short(value: str, max_len: int = 38) -> str:
    value = str(value)
    return value if len(value) <= max_len else value[: max_len - 3] + "..."


def _display_method(value: str) -> str:
    return {"ConstrainedEnsemble": "Ensemble"}.get(str(value), str(value))


def _pil_font(size: int):
    from PIL import ImageFont

    for name in ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()
