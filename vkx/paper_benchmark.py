from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np
import pandas as pd


EXTERNAL_METHOD_ORDER = ["scGen", "CPA", "GEARS", "CellOT"]


def make_paper_benchmark_package(
    formal_dir: str | Path,
    result_dirs: list[str | Path],
    out_dir: str | Path,
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    formal_path = Path(formal_dir)
    result_paths = [Path(item) for item in result_dirs if str(item).strip()]

    metrics = _read_csv(formal_path / "method_metric_comparison.csv")
    availability = _read_csv(formal_path / "method_availability.csv")
    roc_points = _read_csv(formal_path / "formal_benchmark_roc_points.csv")
    predictions = _read_csv(formal_path / "formal_benchmark_predictions.csv")
    truth = _read_csv(formal_path / "formal_benchmark_truth.csv")
    result_summary = _collect_result_summaries(result_paths)

    metrics.to_csv(out / "paper_method_metrics.csv", index=False)
    availability.to_csv(out / "paper_method_availability.csv", index=False)
    result_summary.to_csv(out / "paper_result_summary.csv", index=False)

    _plot_method_leaderboard(metrics, availability, out)
    _plot_roc_curves(roc_points, metrics, out)
    _plot_delta_heatmap(predictions, truth, out)
    _plot_result_gallery(result_summary, out)
    _copy_key_figures(formal_path, result_paths, out)
    _write_report(out, formal_path, result_paths, metrics, availability, result_summary)
    return {
        "metrics": metrics,
        "availability": availability,
        "result_summary": result_summary,
    }


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def _collect_result_summaries(result_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    for path in result_dirs:
        summary = _read_csv(path / "summary.csv")
        auc = _read_csv(path / "auc_summary.csv")
        if summary.empty:
            continue
        row = summary.iloc[0].to_dict()
        row["result_id"] = path.name
        row["result_dir"] = str(path)
        row["roc_auc"] = float(auc["roc_auc"].iloc[0]) if not auc.empty and "roc_auc" in auc.columns and pd.notna(auc["roc_auc"].iloc[0]) else np.nan
        row["task_type"] = _guess_task_type(path.name, row)
        rows.append(row)
    cols = [
        "result_id",
        "task_type",
        "dataset",
        "input_modality",
        "state_representation",
        "n_ko",
        "n_features",
        "mean_direction_cosine",
        "mean_abs_delta_error",
        "mean_distribution_improvement",
        "improved_fraction",
        "roc_auc",
        "calibration_method",
        "shape_calibration_method",
        "result_dir",
    ]
    if not rows:
        return pd.DataFrame(columns=cols)
    table = pd.DataFrame(rows)
    for col in cols:
        if col not in table.columns:
            table[col] = np.nan
    return table[cols]


def _guess_task_type(name: str, row: dict) -> str:
    text = name.lower()
    modality = str(row.get("input_modality", "")).lower() + " " + str(row.get("state_representation", "")).lower()
    is_multi = any(key in modality for key in ["adt", "protein", "multi", "gdo"])
    if "double" in text or int(row.get("n_ko", 1) or 1) > 1:
        return "double KO + multimodal" if is_multi else "double KO"
    if "atac" in text or "peak" in text:
        return "ATAC peak"
    return "single KO + multimodal" if is_multi else "single KO"


def _plot_method_leaderboard(metrics: pd.DataFrame, availability: pd.DataFrame, out: Path) -> None:
    if metrics.empty:
        return
    rows = metrics.copy()
    for col in ["roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate"]:
        if col in rows.columns:
            rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows["_score"] = (
        rows.get("roc_auc", np.nan).fillna(-1)
        + 0.25 * rows.get("direction_cosine", np.nan).fillna(-1)
        + 0.15 * rows.get("r2", np.nan).fillna(-1)
        - 0.25 * rows.get("mae", np.nan).fillna(1)
    )
    rows = rows.sort_values("_score", ascending=False)
    missing = []
    if not availability.empty and "method" in availability.columns:
        for method in EXTERNAL_METHOD_ORDER:
            subset = availability.loc[availability["method"].astype(str).str.lower() == method.lower()]
            if not subset.empty:
                missing.append({"method": method, "status": str(subset["status"].iloc[0])})
    rows["status"] = "scored"
    _draw_method_leaderboard(
        rows=rows,
        missing=pd.DataFrame(missing),
        out_path=out / "01_method_leaderboard.png",
        title="VKX formal benchmark against scored baselines",
        subtitle="Only same-dataset predictions are scored; external deep models are tracked as benchmark slots.",
        metric_specs=[
            ("roc_auc", "AUC", 0.0, 1.0, False),
            ("direction_cosine", "Direction", 0.0, 1.0, False),
            ("r2", "R2", -0.5, 1.0, False),
            ("mae", "MAE", 0.0, max(0.01, float(pd.to_numeric(rows["mae"], errors="coerce").max(skipna=True) or 0.01)), True),
        ],
    )


def _plot_roc_curves(roc_points: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    if roc_points.empty:
        return
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    auc_map = {}
    if not metrics.empty and {"method", "roc_auc"}.issubset(metrics.columns):
        auc_map = dict(zip(metrics["method"].astype(str), pd.to_numeric(metrics["roc_auc"], errors="coerce")))
    img = Image.new("RGB", (1180, 860), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(25)
    font = _pil_font(16)
    small = _pil_font(13)
    draw.text((42, 25), "ROC curves: can the method rank strong KO-response features first?", fill=(20, 20, 20), font=title_font)
    draw.text((42, 60), "Higher and more left-up curves are better. This is the curve version of AUC, not a bar-only summary.", fill=(80, 80, 80), font=font)
    left, top, size = 115, 125, 610
    draw.rectangle((left, top, left + size, top + size), outline=(170, 170, 170), width=1)
    for tick in [0.0, 0.25, 0.5, 0.75, 1.0]:
        x = left + int(tick * size)
        y = top + size - int(tick * size)
        draw.line((x, top, x, top + size), fill=(238, 238, 238))
        draw.line((left, y, left + size, y), fill=(238, 238, 238))
        draw.text((x - 12, top + size + 9), f"{tick:.2g}", fill=(50, 50, 50), font=small)
        draw.text((left - 45, y - 8), f"{tick:.2g}", fill=(50, 50, 50), font=small)
    draw.line((left, top + size, left + size, top), fill=(150, 150, 150), width=1)
    colors = _method_colors()
    legend_y = 130
    for method, group in roc_points.groupby("method", observed=True):
        curve = group.groupby("fpr", as_index=False)["tpr"].mean().sort_values("fpr")
        points = [(left + int(float(r.fpr) * size), top + size - int(float(r.tpr) * size)) for r in curve.itertuples()]
        color = colors.get(str(method), (70, 70, 70))
        if len(points) > 1:
            draw.line(points, fill=color, width=4)
        for x, y in points:
            draw.ellipse((x - 3, y - 3, x + 3, y + 3), fill=color)
        draw.rectangle((775, legend_y, 805, legend_y + 18), fill=color)
        auc = auc_map.get(str(method), np.nan)
        label = f"{method}: AUC {auc:.3f}" if np.isfinite(auc) else str(method)
        draw.text((820, legend_y - 2), label, fill=(30, 30, 30), font=font)
        legend_y += 38
    draw.text((left + 205, top + size + 45), "False positive rate", fill=(30, 30, 30), font=font)
    draw.text((18, top + 270), "True positive rate", fill=(30, 30, 30), font=font)
    img.save(out / "02_auc_roc_curves.png")


def _plot_delta_heatmap(predictions: pd.DataFrame, truth: pd.DataFrame, out: Path) -> None:
    if predictions.empty or truth.empty:
        return
    merged = predictions.merge(truth, on="ko_target", how="inner")
    if merged.empty:
        return
    features = [col.removeprefix("true_delta_") for col in truth.columns if col.startswith("true_delta_")]
    scores = []
    for feature in features:
        values = pd.to_numeric(merged.get(f"true_delta_{feature}"), errors="coerce")
        scores.append((feature, float(np.nanmean(np.abs(values)))))
    chosen = [feature for feature, _ in sorted(scores, key=lambda x: x[1], reverse=True)[:8]]
    keep_methods = ["ResponseBoosted", "ConstrainedEnsemble", "VKX"]
    rows: list[tuple[str, list[float]]] = []
    for _, row in truth.iterrows():
        rows.append((f"TRUE | {row['ko_target']}", [float(row.get(f"true_delta_{feature}", np.nan)) for feature in chosen]))
    merged["_method_rank"] = merged["method"].map({name: i for i, name in enumerate(keep_methods)}).fillna(999)
    for _, row in merged.loc[merged["method"].isin(keep_methods)].sort_values(["_method_rank", "ko_target"]).iterrows():
        values = [float(row.get(f"pred_delta_{feature}", np.nan)) for feature in chosen]
        if np.isfinite(values).any():
            rows.append((f"{_method_short(row['method'])} | {row['ko_target']}", values))
    _draw_heatmap(
        rows,
        chosen,
        out / "03_real_vs_virtual_method_heatmap.png",
        "Real vs virtual KO delta heatmap",
        "Limited to the strongest features and key scored methods; full values remain in the CSV tables.",
    )


def _plot_result_gallery(summary: pd.DataFrame, out: Path) -> None:
    if summary.empty:
        return
    rows = summary.copy()
    for col in ["mean_direction_cosine", "mean_abs_delta_error", "mean_distribution_improvement", "improved_fraction", "roc_auc"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows["display_label"] = rows.apply(_display_result_label, axis=1)
    _draw_metric_table(
        rows=rows,
        out_path=out / "04_single_double_multimodal_gallery.png",
        title="Single KO, double KO and multimodal/ATAC result gallery",
        subtitle="This panel asks whether VKX is useful beyond one RNA-only example.",
        label_col="display_label",
        metric_specs=[
            ("roc_auc", "AUC", 0.0, 1.0, False),
            ("mean_direction_cosine", "Direction", 0.0, 1.0, False),
            ("mean_distribution_improvement", "Distrib", -0.5, 1.0, False),
            ("mean_abs_delta_error", "MAE", 0.0, max(0.01, float(rows["mean_abs_delta_error"].max(skipna=True) or 0.01)), True),
        ],
        extra_col="task_type",
    )


def _copy_key_figures(formal_dir: Path, result_dirs: list[Path], out: Path) -> None:
    asset_dir = out / "assets"
    if asset_dir.exists():
        shutil.rmtree(asset_dir)
    asset_dir.mkdir(parents=True, exist_ok=True)
    candidates = [
        formal_dir / "01_formal_benchmark_metric_panel.png",
        formal_dir / "02_formal_benchmark_delta_heatmap.png",
        formal_dir / "03_method_availability.png",
        formal_dir / "04_formal_benchmark_roc_curves.png",
    ]
    for result in result_dirs:
        for name in [
            "02_true_vs_virtual_heatmap.png",
            "03_cell_state_umap.png",
            "04_auc_strong_response_roc.png",
            "05_atac_peak_level_changes.png",
        ]:
            candidates.append(result / name)
    rows = []
    for path in candidates:
        if not path.exists():
            continue
        dest_name = _unique(asset_dir, f"{path.parent.name}_{path.name}")
        shutil.copy2(path, asset_dir / dest_name)
        rows.append({"figure": dest_name, "source": str(path), "why_it_matters": _figure_explanation(path.name)})
    pd.DataFrame(rows).to_csv(out / "paper_figure_index.csv", index=False)


def _write_report(
    out: Path,
    formal_dir: Path,
    result_dirs: list[Path],
    metrics: pd.DataFrame,
    availability: pd.DataFrame,
    result_summary: pd.DataFrame,
) -> None:
    best = "not available"
    if not metrics.empty:
        scored = metrics.copy()
        scored["_score"] = (
            pd.to_numeric(scored.get("roc_auc"), errors="coerce").fillna(-1)
            + 0.25 * pd.to_numeric(scored.get("direction_cosine"), errors="coerce").fillna(-1)
            + 0.15 * pd.to_numeric(scored.get("r2"), errors="coerce").fillna(-1)
            - 0.25 * pd.to_numeric(scored.get("mae"), errors="coerce").fillna(1)
        )
        row = scored.sort_values("_score", ascending=False).iloc[0]
        best = f"{row['method']} (AUC={row.get('roc_auc', np.nan):.3f}, R2={row.get('r2', np.nan):.3f}, MAE={row.get('mae', np.nan):.3f})"
    external_status = "not checked"
    if not availability.empty:
        sub = availability.loc[availability["method"].astype(str).isin(EXTERNAL_METHOD_ORDER)]
        if not sub.empty:
            external_status = sub[["method", "status"]].to_string(index=False)
    text = [
        "# VKX paper benchmark package",
        "",
        "## One-sentence conclusion",
        "",
        "VKX is currently best described as a small-sample, interpretable, prior-constrained virtual KO framework; the current evidence shows improvement over internal classical baselines, while scGen/CPA/GEARS/CellOT still require same-dataset external prediction files before making a strict superiority claim.",
        "",
        "## Main figures",
        "",
        "![Method leaderboard](01_method_leaderboard.png)",
        "",
        "![AUC ROC curves](02_auc_roc_curves.png)",
        "",
        "![Real vs virtual heatmap](03_real_vs_virtual_method_heatmap.png)",
        "",
        "![Result gallery](04_single_double_multimodal_gallery.png)",
        "",
        "## Best current method in the formal benchmark",
        "",
        best,
        "",
        "## External method status",
        "",
        "```text",
        external_status,
        "```",
        "",
        "## What has been done for the three requested items",
        "",
        "1. Formal benchmark: VKX is compared with ridge/PLS/additive and explicit scGen/CPA/GEARS/CellOT slots. External methods are not silently counted unless their prediction CSV is supplied.",
        "2. VKX optimization: boosted/adaptive response anchor, double-KO interaction residual, and ATAC quantile/zero-inflated shape calibration are promoted as the current optimization path.",
        "3. Visualization: the package now contains a method leaderboard, ROC curves, real-vs-virtual heatmap, single/double/multimodal gallery, and copied source figures in `assets/`.",
        "",
        "## Remaining weaknesses",
        "",
        "- This package does not prove VKX is stronger than scGen/CPA/GEARS/CellOT until those methods are run on the same input and imported through `--external-predictions-csv`.",
        "- Double-KO nonlinear effects remain the main scientific difficulty, especially MAPK/TGFB-like programs.",
        "- Full public RNA+ADT+ATAC+perturbation labelled benchmark is still not confirmed.",
        "- ATAC peak-level distribution shape is improved by quantile calibration but remains harder than pathway/protein state prediction.",
        "",
        "## Inputs used",
        "",
        f"- Formal benchmark directory: `{formal_dir}`",
        "- Result directories:",
    ]
    text.extend([f"  - `{path}`" for path in result_dirs])
    text.extend(
        [
            "",
            "## Tables",
            "",
            "- `paper_method_metrics.csv`",
            "- `paper_method_availability.csv`",
            "- `paper_result_summary.csv`",
            "- `paper_figure_index.csv`",
        ]
    )
    if not result_summary.empty:
        text.extend(["", "## Result summary", "", "```text", result_summary.to_string(index=False), "```"])
    (out / "paper_benchmark_report_zh.md").write_text("\n".join(text), encoding="utf-8")


def _draw_method_leaderboard(
    rows: pd.DataFrame,
    missing: pd.DataFrame,
    out_path: Path,
    title: str,
    subtitle: str,
    metric_specs: list[tuple[str, str, float, float, bool]],
) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    rows = rows.head(8).copy()
    row_h = 66
    width = 1680
    height = max(560, 150 + row_h * len(rows) + 80)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(26)
    font = _pil_font(16)
    small = _pil_font(13)
    draw.text((36, 24), title, fill=(18, 18, 18), font=title_font)
    draw.text((36, 62), subtitle, fill=(85, 85, 85), font=font)
    draw.text((36, 104), "Scored methods", fill=(25, 25, 25), font=font)
    x_label = 36
    metric_x = [420, 650, 880, 1110]
    colors = _method_colors()
    for i, (_, row) in enumerate(rows.iterrows()):
        y = 142 + i * row_h
        method = str(row.get("method", ""))
        draw.rectangle((x_label, y + 8, x_label + 10, y + 38), fill=colors.get(method, (42, 157, 143)))
        draw.text((x_label + 18, y + 8), _method_short(method), fill=(25, 25, 25), font=font)
        for x0, spec in zip(metric_x, metric_specs):
            col, name, vmin, vmax, lower_better = spec
            if i == 0:
                draw.text((x0, 104), name, fill=(30, 30, 30), font=small)
            value = pd.to_numeric(pd.Series([row.get(col, np.nan)]), errors="coerce").iloc[0]
            if not np.isfinite(value):
                frac = 0.0
                label = "NA"
            else:
                raw = min(max((float(value) - vmin) / (vmax - vmin + 1e-9), 0.0), 1.0)
                frac = 1.0 - raw if lower_better else raw
                label = f"{float(value):.3f}" if abs(float(value)) < 1 else f"{float(value):.2f}"
            bar_color = (76, 120, 168) if lower_better else (42, 157, 143)
            draw.rectangle((x0, y + 10, x0 + 158, y + 34), fill=(246, 246, 246), outline=(222, 222, 222))
            draw.rectangle((x0, y + 10, x0 + int(158 * frac), y + 34), fill=bar_color)
            draw.text((x0 + 170, y + 12), label, fill=(30, 30, 30), font=small)

    box_x = 1360
    box_y = 104
    draw.rounded_rectangle((box_x, box_y, width - 38, height - 58), radius=8, fill=(248, 248, 248), outline=(222, 222, 222))
    draw.text((box_x + 22, box_y + 22), "External methods", fill=(25, 25, 25), font=font)
    draw.text((box_x + 22, box_y + 50), "Not scored until same-data\npredictions are imported.", fill=(90, 90, 90), font=small)
    if missing.empty:
        draw.text((box_x + 22, box_y + 100), "No external slots found.", fill=(90, 90, 90), font=small)
    else:
        for i, (_, row) in enumerate(missing.iterrows()):
            y = box_y + 105 + i * 46
            method = str(row.get("method", ""))
            status = str(row.get("status", "not_run"))
            draw.rectangle((box_x + 22, y + 6, box_x + 34, y + 26), fill=_method_colors().get(method, (120, 120, 120)))
            draw.text((box_x + 44, y), method, fill=(25, 25, 25), font=font)
            draw.text((box_x + 155, y + 2), status, fill=(120, 45, 45) if status == "not_run" else (45, 110, 75), font=small)
    img.save(out_path)


def _draw_metric_table(
    rows: pd.DataFrame,
    out_path: Path,
    title: str,
    subtitle: str,
    label_col: str,
    metric_specs: list[tuple[str, str, float, float, bool]],
    extra_col: str,
) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    rows = rows.head(16).copy()
    row_h = 58
    width = 1650
    height = max(460, 135 + row_h * len(rows) + 60)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(24)
    font = _pil_font(15)
    small = _pil_font(12)
    draw.text((35, 24), title, fill=(20, 20, 20), font=title_font)
    draw.text((35, 58), subtitle, fill=(80, 80, 80), font=font)
    x_label = 35
    x_extra = 450
    metric_x = [690, 915, 1140, 1365]
    colors = _method_colors()
    for i, (_, row) in enumerate(rows.iterrows()):
        y = 112 + i * row_h
        label = str(row.get(label_col, ""))
        color = colors.get(label, (42, 157, 143))
        draw.rectangle((x_label, y + 7, x_label + 8, y + 34), fill=color)
        draw.text((x_label + 16, y + 8), _short(label, 44), fill=(30, 30, 30), font=font)
        draw.text((x_extra, y + 8), _short(row.get(extra_col, ""), 25), fill=(80, 80, 80), font=font)
        for x0, spec in zip(metric_x, metric_specs):
            col, name, vmin, vmax, lower_better = spec
            if i == 0:
                draw.text((x0, 88), name, fill=(30, 30, 30), font=small)
            value = pd.to_numeric(pd.Series([row.get(col, np.nan)]), errors="coerce").iloc[0]
            if not np.isfinite(value):
                frac = 0.0
                text = "NA"
                bar_color = (190, 190, 190)
            else:
                raw = min(max((float(value) - vmin) / (vmax - vmin + 1e-9), 0.0), 1.0)
                frac = 1.0 - raw if lower_better else raw
                text = f"{float(value):.3f}" if abs(float(value)) < 1 else f"{float(value):.2f}"
                bar_color = (76, 120, 168) if lower_better else (42, 157, 143)
            draw.rectangle((x0, y + 8, x0 + 160, y + 30), fill=(245, 245, 245), outline=(220, 220, 220))
            draw.rectangle((x0, y + 8, x0 + int(160 * frac), y + 30), fill=bar_color)
            draw.text((x0 + 170, y + 10), text, fill=(30, 30, 30), font=small)
    img.save(out_path)


def _draw_heatmap(rows: list[tuple[str, list[float]]], features: list[str], out_path: Path, title: str, subtitle: str) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    if not rows or not features:
        return
    labels = [_feature_label(feature) for feature in features]
    cell_w = 112
    cell_h = 30
    left = 340
    top = 130
    width = max(1200, left + cell_w * len(features) + 70)
    height = max(520, top + cell_h * len(rows) + 80)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(24)
    font = _pil_font(14)
    small = _pil_font(11)
    draw.text((35, 24), title, fill=(20, 20, 20), font=title_font)
    draw.text((35, 58), subtitle, fill=(80, 80, 80), font=font)
    matrix = np.asarray([values for _, values in rows], dtype=float)
    vmax = max(float(np.nanmax(np.abs(matrix))) if matrix.size else 0.0, 1e-9)
    for j, label in enumerate(labels):
        parts = label.split("\n")
        for k, part in enumerate(parts[:2]):
            draw.text((left + j * cell_w, 88 + k * 14), _short(part, 16), fill=(30, 30, 30), font=small)
    for i, (label, values) in enumerate(rows):
        y = top + i * cell_h
        draw.text((35, y + 6), _short(label, 38), fill=(30, 30, 30), font=font)
        for j, value in enumerate(values):
            x = left + j * cell_w
            color = _delta_color(value, vmax)
            draw.rectangle((x, y, x + cell_w - 4, y + cell_h - 4), fill=color, outline=(255, 255, 255))
            if np.isfinite(value) and abs(value) / vmax > 0.45:
                draw.text((x + 18, y + 7), f"{value:.2f}", fill=(255, 255, 255), font=small)
            elif np.isfinite(value):
                draw.text((x + 18, y + 7), f"{value:.2f}", fill=(40, 40, 40), font=small)
    img.save(out_path)


def _delta_color(value: float, vmax: float) -> tuple[int, int, int]:
    if not np.isfinite(value):
        return (235, 235, 235)
    frac = min(abs(float(value)) / vmax, 1.0)
    base = int(255 - 175 * frac)
    if value >= 0:
        return (190, base, base)
    return (base, base, 190)


def _method_colors() -> dict[str, tuple[int, int, int]]:
    return {
        "ResponseBoosted": (231, 111, 81),
        "VKX": (231, 111, 81),
        "ConstrainedEnsemble": (42, 157, 143),
        "CalibratedEnsemble": (42, 157, 143),
        "Ridge": (76, 120, 168),
        "PLS": (244, 162, 97),
        "Additive": (141, 153, 174),
        "scGen": (123, 44, 191),
        "CPA": (67, 97, 238),
        "GEARS": (58, 134, 255),
        "CellOT": (106, 153, 78),
    }


def _method_short(method: object) -> str:
    return {
        "ResponseBoosted": "VKX-Boosted",
        "ConstrainedEnsemble": "Ensemble",
        "CalibratedEnsemble": "Calibrated",
    }.get(str(method), str(method))


def _display_result_label(row: pd.Series) -> str:
    text = str(row.get("result_id", ""))
    dataset = str(row.get("dataset", ""))
    task = str(row.get("task_type", ""))
    if "papalexi" in dataset.lower() or "single_gene" in text:
        return "Papalexi single KO\nRNA+ADT"
    if "hmpcite" in text.lower() or "hmpcite" in dataset.lower():
        return "HMPCITE double KO\nRNA+ADT"
    if "double_gene" in text:
        return "Norman double KO\nRNA programs"
    if "atac" in text.lower() or "ATAC" in task:
        return "scPerturb ATAC\npeak-level"
    return text.replace("_", " ")


def _feature_label(feature: object) -> str:
    text = str(feature)
    for prefix in ["pathway_", "program_", "protein_", "atac_", "tf_", "peak_"]:
        if text.startswith(prefix):
            text = text[len(prefix) :]
    lower = text.lower()
    replacements = [
        ("interferon_gamma_response", "IFN gamma\nresponse"),
        ("interferon_gamma_signaling", "IFN gamma\nsignaling"),
        ("interferon_alpha_beta_signaling", "IFN alpha/beta\nsignaling"),
        ("interferon_alpha_response", "IFN alpha\nresponse"),
        ("interferon_signaling", "IFN\nsignaling"),
        ("innate_immune_system", "Innate\nimmune"),
        ("adaptive_immune_system", "Adaptive\nimmune"),
        ("cytokine_signaling_in_immune_system", "Cytokine\nimmune"),
        ("immune_system", "Immune\nsystem"),
        ("il_6_jak_stat3_signaling", "IL6 JAK\nSTAT3"),
        ("il_2_stat5_signaling", "IL2\nSTAT5"),
        ("myc_targets_v1", "MYC\ntargets"),
        ("myc_encode", "MYC\nENCODE"),
        ("myc_chea", "MYC\nChEA"),
        ("mapk_family_signaling_cascades", "MAPK family\ncascades"),
        ("mapk1_mapk3_signaling", "MAPK1/3\nsignaling"),
    ]
    for key, label in replacements:
        if lower == key:
            return label
    text = text.replace("_Signaling", "").replace("_Response", "").replace("_System", "")
    text = text.replace("_ENCODE", "").replace("_CHEA", "")
    parts = [part for part in text.split("_") if part]
    if len(parts) <= 2:
        return " ".join(parts)
    mid = max(1, len(parts) // 2)
    return " ".join(parts[:mid]) + "\n" + " ".join(parts[mid:])


def _figure_explanation(name: str) -> str:
    lower = name.lower()
    if "heatmap" in lower:
        return "Shows what was knocked out, which features changed, and where virtual KO differs from real KO."
    if "roc" in lower or "auc" in lower:
        return "Curve-based AUC evidence for ranking strong KO-response features."
    if "umap" in lower or "pca" in lower:
        return "Visualizes whether virtual KO cells move from control toward true KO state."
    if "peak" in lower:
        return "ATAC peak-level regulatory view, including locus or peak accessibility changes."
    if "availability" in lower:
        return "Prevents over-claiming by showing which methods were actually run."
    return "Benchmark figure copied from the source result directory."


def _unique(directory: Path, name: str) -> str:
    candidate = name
    stem = Path(name).stem
    suffix = Path(name).suffix
    i = 2
    while (directory / candidate).exists():
        candidate = f"{stem}_{i}{suffix}"
        i += 1
    return candidate


def _short(value: object, max_len: int) -> str:
    text = str(value)
    return text if len(text) <= max_len else text[: max_len - 3] + "..."


def _pil_font(size: int):
    from PIL import ImageFont

    for name in ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()
