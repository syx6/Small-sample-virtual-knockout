from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


DIFFICULT_PROGRAM_KEYWORDS = ("MAPK", "TGFB", "TGF_BETA", "TGF-")
ATAC_KEYWORDS = ("atac", "peak", "chromvar", "motif")


def diagnose_virtual_ko_results(
    delta_csv: str | Path,
    out_dir: str | Path,
    manifest_csv: str | Path | None = None,
    confidence_csv: str | Path | None = None,
    min_true_delta: float = 0.03,
    large_error: float | None = None,
    max_features_per_ko: int = 25,
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    delta = pd.read_csv(delta_csv)
    manifest = pd.read_csv(manifest_csv) if manifest_csv and Path(manifest_csv).exists() else pd.DataFrame()
    confidence = pd.read_csv(confidence_csv) if confidence_csv and Path(confidence_csv).exists() else pd.DataFrame()
    feature_table = _make_feature_diagnosis(delta, manifest, min_true_delta=min_true_delta, large_error=large_error)
    ko_table = _make_ko_diagnosis(feature_table, delta, confidence, max_features_per_ko=max_features_per_ko)
    feature_table.to_csv(out / "feature_failure_diagnosis.csv", index=False)
    ko_table.to_csv(out / "ko_failure_diagnosis.csv", index=False)
    _plot_ko_overview(ko_table, out)
    _plot_feature_error_heatmap(feature_table, out, max_features_per_ko=max_features_per_ko)
    _write_report(ko_table, feature_table, out, Path(delta_csv))
    return {"ko": ko_table, "feature": feature_table}


def _make_feature_diagnosis(
    delta: pd.DataFrame,
    manifest: pd.DataFrame,
    min_true_delta: float,
    large_error: float | None,
) -> pd.DataFrame:
    features = sorted({col.removeprefix("pred_delta_") for col in delta.columns if col.startswith("pred_delta_")})
    meta = _manifest_lookup(manifest)
    rows = []
    for _, row in delta.iterrows():
        ko = str(row.get("ko_target", row.get("target", row.get("KO", "unknown"))))
        for feature in features:
            pred = _num(row.get(f"pred_delta_{feature}", np.nan))
            true = _num(row.get(f"true_delta_{feature}", np.nan))
            has_truth = np.isfinite(true)
            abs_error = abs(pred - true) if has_truth and np.isfinite(pred) else np.nan
            direction_match = np.sign(pred) == np.sign(true) if has_truth and abs(true) >= min_true_delta and np.isfinite(pred) else np.nan
            rows.append(
                {
                    "ko_target": ko,
                    "feature": feature,
                    "modality": _infer_modality(feature, meta.get(feature, {})),
                    "true_delta": true if has_truth else np.nan,
                    "pred_delta": pred,
                    "abs_error": abs_error,
                    "direction_match": direction_match,
                    "true_is_small": bool(has_truth and abs(true) < min_true_delta),
                    "is_difficult_program": _is_difficult_program(feature),
                    "issue_tags": "",
                    "severity": "prediction_only" if not has_truth else "ok",
                    "plain_explanation": "",
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    threshold = large_error
    if threshold is None and table["abs_error"].notna().any():
        threshold = max(float(table["abs_error"].quantile(0.75)), min_true_delta * 2)
    if threshold is None:
        threshold = min_true_delta * 2
    table["large_error_threshold"] = threshold
    diagnoses = [_diagnose_feature(row, min_true_delta=min_true_delta, large_error=threshold) for _, row in table.iterrows()]
    table["issue_tags"] = [item[0] for item in diagnoses]
    table["severity"] = [item[1] for item in diagnoses]
    table["plain_explanation"] = [item[2] for item in diagnoses]
    return table


def _make_ko_diagnosis(
    feature_table: pd.DataFrame,
    delta: pd.DataFrame,
    confidence: pd.DataFrame,
    max_features_per_ko: int,
) -> pd.DataFrame:
    rows = []
    if feature_table.empty:
        return pd.DataFrame()
    for ko, sub in feature_table.groupby("ko_target", observed=True):
        truth = sub["true_delta"].notna().any()
        evaluable = sub.loc[sub["direction_match"].notna()]
        risky = sub.loc[sub["severity"].isin(["warning", "high_risk"])]
        top = (
            sub.assign(rank_value=sub["abs_error"].fillna(sub["pred_delta"].abs()))
            .sort_values("rank_value", ascending=False)
            .head(max_features_per_ko)
        )
        row = {
            "ko_target": ko,
            "analysis_mode": "evaluation" if truth else "prediction_only",
            "n_features": int(len(sub)),
            "n_evaluable_features": int(len(evaluable)),
            "direction_match_fraction": float(evaluable["direction_match"].mean()) if len(evaluable) else np.nan,
            "mean_abs_error": float(sub["abs_error"].mean()) if truth else np.nan,
            "median_abs_error": float(sub["abs_error"].median()) if truth else np.nan,
            "n_high_risk_features": int((sub["severity"] == "high_risk").sum()),
            "n_warning_features": int((sub["severity"] == "warning").sum()),
            "top_issue_features": "; ".join(top.loc[top["severity"] != "ok", "feature"].astype(str).head(8)),
            "top_changed_features": "; ".join(sub.reindex(sub["pred_delta"].abs().sort_values(ascending=False).index)["feature"].astype(str).head(8)),
            "main_issue_tags": _top_tags(risky["issue_tags"]),
            "plain_readout": _plain_ko_readout(ko, sub, truth),
        }
        row.update(_lookup_confidence(confidence, ko))
        rows.append(row)
    table = pd.DataFrame(rows)
    if not table.empty:
        table = table.sort_values(["n_high_risk_features", "mean_abs_error"], ascending=[False, False], na_position="last")
    return table


def _diagnose_feature(row: pd.Series, min_true_delta: float, large_error: float) -> tuple[str, str, str]:
    tags: list[str] = []
    if not np.isfinite(row["true_delta"]):
        tags.append("prediction_only_no_truth")
        if abs(row["pred_delta"]) < min_true_delta:
            tags.append("small_predicted_shift")
        return ";".join(tags), "prediction_only", "No real KO label is available in this input, so this is a virtual application result, not an accuracy estimate."
    if row["true_is_small"]:
        tags.append("near_zero_real_effect")
    if bool(row["direction_match"]) is False:
        tags.append("direction_mismatch")
    if np.isfinite(row["abs_error"]) and row["abs_error"] >= large_error:
        tags.append("large_error")
    if abs(row["true_delta"]) >= large_error and abs(row["pred_delta"]) < min_true_delta:
        tags.append("missed_real_effect")
    if abs(row["pred_delta"]) >= large_error and abs(row["true_delta"]) < min_true_delta:
        tags.append("overcalled_effect")
    modality = str(row.get("modality", "")).lower()
    if any(key in modality for key in ATAC_KEYWORDS) and ("large_error" in tags or "direction_mismatch" in tags):
        tags.append("sparse_atac_shape_or_prior_issue")
    if row["is_difficult_program"] and ("large_error" in tags or "direction_mismatch" in tags):
        tags.append("difficult_mapk_tgfb_program")
    if not tags:
        return "", "ok", "Virtual KO agrees with the real KO direction and has no major error flag."
    severity = "high_risk" if "direction_mismatch" in tags or "missed_real_effect" in tags else "warning"
    explanation = _plain_feature_explanation(tags)
    return ";".join(tags), severity, explanation


def _plain_feature_explanation(tags: list[str]) -> str:
    if "difficult_mapk_tgfb_program" in tags:
        return "This is a MAPK/TGFB-like program where nonlinear signalling and incomplete priors often make prediction harder."
    if "sparse_atac_shape_or_prior_issue" in tags:
        return "This ATAC/peak feature may be affected by sparse accessibility, zero inflation, or incomplete peak-gene/motif priors."
    if "direction_mismatch" in tags:
        return "The predicted change goes in the opposite direction from the observed KO change."
    if "missed_real_effect" in tags:
        return "The real KO effect is strong, but the virtual KO effect is too small."
    if "overcalled_effect" in tags:
        return "The virtual KO predicts a strong change, but the observed KO effect is weak."
    if "near_zero_real_effect" in tags:
        return "The observed KO effect is close to zero, so AUC/direction-style interpretation is unstable for this feature."
    return "The feature has a larger than expected quantitative error."


def _manifest_lookup(manifest: pd.DataFrame) -> dict[str, dict]:
    if manifest.empty:
        return {}
    feature_col = next((col for col in ["state_feature", "feature", "feature_name", "name"] if col in manifest.columns), None)
    if feature_col is None:
        return {}
    return {str(row[feature_col]): row.to_dict() for _, row in manifest.iterrows()}


def _infer_modality(feature: str, metadata: dict) -> str:
    for key in ["modality", "source", "obsm_key", "feature_type"]:
        value = metadata.get(key)
        if value is not None and str(value) and str(value) != "nan":
            return str(value)
    lower = feature.lower()
    if "protein_" in lower or lower.startswith("adt_"):
        return "protein/ADT"
    if any(key in lower for key in ATAC_KEYWORDS):
        return "ATAC/peak"
    if "pathway" in lower or "program" in lower:
        return "pathway/program"
    return "state_score"


def _lookup_confidence(confidence: pd.DataFrame, ko: str) -> dict:
    if confidence.empty or "ko_target" not in confidence.columns:
        return {"transfer_confidence": ""}
    hit = confidence.loc[confidence["ko_target"].astype(str) == ko]
    if hit.empty:
        return {"transfer_confidence": ""}
    row = hit.iloc[0]
    value = row.get("confidence", row.get("transfer_confidence", ""))
    return {"transfer_confidence": value}


def _top_tags(series: pd.Series) -> str:
    counts: dict[str, int] = {}
    for value in series.dropna().astype(str):
        for tag in value.split(";"):
            if tag:
                counts[tag] = counts.get(tag, 0) + 1
    return "; ".join([f"{tag}({count})" for tag, count in sorted(counts.items(), key=lambda item: item[1], reverse=True)[:5]])


def _plain_ko_readout(ko: str, sub: pd.DataFrame, truth: bool) -> str:
    top_changes = sub.reindex(sub["pred_delta"].abs().sort_values(ascending=False).index)["feature"].astype(str).head(3).tolist()
    if not truth:
        return f"{ko}: prediction-only result. Largest predicted changes: {', '.join(top_changes)}."
    high = int((sub["severity"] == "high_risk").sum())
    warn = int((sub["severity"] == "warning").sum())
    evaluable = sub.loc[sub["direction_match"].notna()]
    match = float(evaluable["direction_match"].mean()) if len(evaluable) else np.nan
    if high == 0 and warn <= 2:
        return f"{ko}: generally consistent with real KO; largest predicted changes: {', '.join(top_changes)}."
    match_text = "n/a" if not np.isfinite(match) else f"{match:.2f}"
    return f"{ko}: needs caution. Direction-match fraction is {match_text}, with {high} high-risk and {warn} warning features."


def _plot_ko_overview(ko_table: pd.DataFrame, out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        _plot_ko_overview_pillow(ko_table, out_dir)
        return
    if ko_table.empty:
        return
    plot = ko_table.copy().head(20)
    labels = plot["ko_target"].astype(str)
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4.8, 0.38 * len(plot) + 2.0)), constrained_layout=True)
    axes[0].barh(labels, plot["direction_match_fraction"].fillna(0), color="#2A9D8F")
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Direction match fraction")
    axes[0].set_title("Direction agreement")
    axes[0].invert_yaxis()
    width = np.arange(len(plot))
    axes[1].barh(width, plot["n_warning_features"], color="#E9C46A", label="warning")
    axes[1].barh(width, plot["n_high_risk_features"], left=plot["n_warning_features"], color="#E76F51", label="high risk")
    axes[1].set_yticks(width)
    axes[1].set_yticklabels(labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel("Flagged features")
    axes[1].set_title("Why to be cautious")
    axes[1].legend(frameon=False)
    fig.suptitle("Virtual KO Diagnosis Overview")
    fig.savefig(out_dir / "01_failure_diagnosis_overview.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_feature_error_heatmap(feature_table: pd.DataFrame, out_dir: Path, max_features_per_ko: int) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        _plot_feature_error_heatmap_pillow(feature_table, out_dir, max_features_per_ko=max_features_per_ko)
        return
    if feature_table.empty or not feature_table["abs_error"].notna().any():
        return
    selected = []
    for _, sub in feature_table.groupby("ko_target", observed=True):
        selected.extend(sub.sort_values("abs_error", ascending=False)["feature"].head(max_features_per_ko).tolist())
    selected = list(dict.fromkeys(selected))[:60]
    plot = feature_table.loc[feature_table["feature"].isin(selected)]
    matrix = plot.pivot_table(index="feature", columns="ko_target", values="abs_error", aggfunc="max").fillna(0)
    if matrix.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.75 * matrix.shape[1] + 3), max(5, 0.24 * matrix.shape[0] + 2)), constrained_layout=True)
    vmax = max(float(np.nanpercentile(matrix.values, 95)), 1e-6)
    image = ax.imshow(matrix.values, aspect="auto", cmap="Reds", vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels([_short(label, 45) for label in matrix.index], fontsize=8)
    ax.set_title("Feature-level absolute error heatmap")
    ax.set_xlabel("KO target")
    ax.set_ylabel("State feature")
    fig.colorbar(image, ax=ax, label="absolute error")
    fig.savefig(out_dir / "02_feature_error_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_ko_overview_pillow(ko_table: pd.DataFrame, out_dir: Path) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    if ko_table.empty:
        return
    plot = ko_table.copy().head(20)
    row_h = 34
    width = 1600
    height = 120 + row_h * len(plot)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    draw.text((30, 25), "Virtual KO Diagnosis Overview", fill=(20, 20, 20), font=font)
    draw.text((30, 60), "Direction match", fill=(60, 60, 60), font=font)
    draw.text((710, 60), "Flagged features: warning + high risk", fill=(60, 60, 60), font=font)
    max_flags = int((plot["n_warning_features"].fillna(0) + plot["n_high_risk_features"].fillna(0)).max())
    flag_scale = min(18.0, 360.0 / max(max_flags, 1))
    for i, (_, row) in enumerate(plot.iterrows()):
        y = 95 + i * row_h
        ko = _short(str(row["ko_target"]), 30)
        draw.text((30, y), ko, fill=(40, 40, 40), font=font)
        match = row.get("direction_match_fraction", np.nan)
        match = 0 if not np.isfinite(match) else float(match)
        draw.rectangle((250, y, 650, y + 16), outline=(210, 210, 210), fill=(245, 245, 245))
        draw.rectangle((250, y, 250 + int(400 * min(max(match, 0), 1)), y + 16), fill=(42, 157, 143))
        draw.text((660, y), f"{match:.2f}", fill=(40, 40, 40), font=font)
        warning = int(row.get("n_warning_features", 0))
        high = int(row.get("n_high_risk_features", 0))
        warning_w = int(warning * flag_scale)
        high_w = int(high * flag_scale)
        draw.rectangle((850, y, 850 + warning_w, y + 16), fill=(233, 196, 106))
        draw.rectangle((850 + warning_w, y, 850 + warning_w + high_w, y + 16), fill=(231, 111, 81))
        draw.text((870 + warning_w + high_w, y), f"{warning} warning, {high} high", fill=(40, 40, 40), font=font)
    img.save(out_dir / "01_failure_diagnosis_overview.png")


def _plot_feature_error_heatmap_pillow(feature_table: pd.DataFrame, out_dir: Path, max_features_per_ko: int) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    if feature_table.empty or not feature_table["abs_error"].notna().any():
        return
    selected = []
    for _, sub in feature_table.groupby("ko_target", observed=True):
        selected.extend(sub.sort_values("abs_error", ascending=False)["feature"].head(max_features_per_ko).tolist())
    selected = list(dict.fromkeys(selected))[:50]
    matrix = (
        feature_table.loc[feature_table["feature"].isin(selected)]
        .pivot_table(index="feature", columns="ko_target", values="abs_error", aggfunc="max")
        .fillna(0)
    )
    if matrix.empty:
        return
    cell_w = 90
    cell_h = 18
    left = 360
    top = 90
    width = left + cell_w * matrix.shape[1] + 40
    height = top + cell_h * matrix.shape[0] + 60
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    draw.text((30, 25), "Feature-level absolute error heatmap", fill=(20, 20, 20), font=font)
    vmax = max(float(np.nanpercentile(matrix.values, 95)), 1e-6)
    for j, col in enumerate(matrix.columns):
        draw.text((left + j * cell_w, 60), _short(str(col), 12), fill=(40, 40, 40), font=font)
    for i, feature in enumerate(matrix.index):
        y = top + i * cell_h
        draw.text((30, y), _short(str(feature), 52), fill=(40, 40, 40), font=font)
        for j, col in enumerate(matrix.columns):
            value = float(matrix.loc[feature, col])
            intensity = int(255 - 200 * min(value / vmax, 1))
            color = (255, intensity, intensity)
            x = left + j * cell_w
            draw.rectangle((x, y, x + cell_w - 4, y + cell_h - 2), fill=color, outline=(240, 240, 240))
            if value > 0:
                draw.text((x + 6, y + 2), f"{value:.2f}", fill=(80, 20, 20), font=font)
    img.save(out_dir / "02_feature_error_heatmap.png")


def _write_report(ko_table: pd.DataFrame, feature_table: pd.DataFrame, out_dir: Path, delta_csv: Path) -> None:
    lines = [
        "# Virtual KO Failure Diagnosis Report",
        "",
        "This report explains which virtual KO results are reliable, which need caution, and why.",
        "",
        f"Input delta table: `{delta_csv}`",
        "",
        "## How to read this",
        "",
        "- `ok`: predicted KO direction and magnitude are broadly consistent with the real KO effect.",
        "- `warning`: quantitative error, near-zero true effect, or possible over-calling needs caution.",
        "- `high_risk`: direction mismatch or missed real effect; do not use this feature as strong evidence.",
        "- `prediction_only`: no real KO labels were available, so this is an application report, not an accuracy benchmark.",
        "",
        "## KO-level summary",
        "",
        ko_table.to_string(index=False) if not ko_table.empty else "No KO-level diagnosis was generated.",
        "",
        "## Most important feature-level warnings",
        "",
    ]
    if feature_table.empty:
        lines.append("No feature-level diagnosis was generated.")
    else:
        flagged = feature_table.loc[feature_table["severity"].isin(["warning", "high_risk"])].copy()
        if flagged.empty:
            lines.append("No warning or high-risk features were detected.")
        else:
            cols = ["ko_target", "feature", "modality", "true_delta", "pred_delta", "abs_error", "severity", "issue_tags", "plain_explanation"]
            lines.append(flagged.sort_values(["severity", "abs_error"], ascending=[True, False])[cols].head(60).to_string(index=False))
    lines.extend(
        [
            "",
            "## Generated files",
            "",
            "- `ko_failure_diagnosis.csv`: one row per KO target.",
            "- `feature_failure_diagnosis.csv`: one row per KO-feature pair.",
            "- `01_failure_diagnosis_overview.png`: readable KO-level diagnosis figure.",
            "- `02_feature_error_heatmap.png`: feature-level error heatmap when real KO labels are available.",
        ]
    )
    (out_dir / "failure_diagnosis_report.md").write_text("\n".join(lines), encoding="utf-8")


def _short(value: str, max_len: int) -> str:
    return value if len(value) <= max_len else value[: max_len - 3] + "..."


def _num(value) -> float:
    try:
        return float(value)
    except Exception:
        return np.nan


def _is_difficult_program(feature: str) -> bool:
    upper = feature.upper()
    return any(keyword in upper for keyword in DIFFICULT_PROGRAM_KEYWORDS)
