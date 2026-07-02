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

    _plot_publication_main_figure(metrics, availability, roc_points, predictions, truth, result_summary, out)
    _plot_method_leaderboard(metrics, availability, out)
    _plot_roc_curves(roc_points, metrics, out)
    _plot_delta_heatmap(predictions, truth, out)
    _plot_result_gallery(result_summary, out)
    _plot_adaptive_improvement(metrics, out)
    _plot_benchmark_completeness(availability, result_summary, out)
    _plot_before_after_umap(result_paths, out)
    _plot_single_double_response_map(result_summary, out)
    _plot_peak_locus_track(result_paths, out)
    _plot_method_radar(metrics, availability, out)
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


def _plot_adaptive_improvement(metrics: pd.DataFrame, out: Path) -> None:
    if metrics.empty:
        return
    rows = metrics.copy()
    for col in ["roc_auc", "direction_cosine", "r2", "mae"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    focus = rows.loc[rows["method"].astype(str).isin(["VKXAdaptive", "ResponseBoosted", "ConstrainedEnsemble", "PLS", "Ridge", "VKX"])].copy()
    if focus.empty:
        return
    focus["_label"] = focus["method"].map(_method_short)
    _draw_metric_table(
        rows=focus.sort_values(["roc_auc", "r2"], ascending=False),
        out_path=out / "05_adaptive_improvement.png",
        title="Adaptive VKX improvement and baseline matching",
        subtitle="VKX-Adaptive selects a stable anchor by cross-validation so the small-sample model can match strong classical baselines.",
        label_col="_label",
        metric_specs=[
            ("roc_auc", "AUC", 0.0, 1.0, False),
            ("direction_cosine", "Direction", 0.0, 1.0, False),
            ("r2", "R2", -0.5, 1.0, False),
            ("mae", "MAE", 0.0, max(0.01, float(focus["mae"].max(skipna=True) or 0.01)), True),
        ],
        extra_col="method",
    )


def _plot_benchmark_completeness(availability: pd.DataFrame, result_summary: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    width, height = 1500, 760
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(26)
    font = _pil_font(17)
    small = _pil_font(13)
    draw.text((36, 28), "Benchmark completeness and remaining gaps", fill=(18, 18, 18), font=title_font)
    draw.text((36, 64), "Completed evidence is separated from external methods or datasets that still need real inputs.", fill=(80, 80, 80), font=font)
    checks = [
        ("Same-data internal benchmark", True, "VKX variants, PLS, Ridge and ensemble are scored together."),
        ("AUC shown as ROC curve", True, "Curve points are generated from held-out KO tasks."),
        ("Single-KO result figure", _has_task(result_summary, "single"), "Papalexi RNA+ADT example."),
        ("Double-KO result figure", _has_task(result_summary, "double"), "Norman/HMPCITE examples."),
        ("ATAC peak-level result figure", _has_task(result_summary, "atac"), "scPerturb ATAC peak-level example."),
        ("scGen/CPA/GEARS/CellOT scored", _external_scored(availability), "Requires external prediction CSV from those methods."),
        ("Full RNA+ADT+ATAC labelled benchmark", False, "Still needs a confirmed public labelled trimodal dataset."),
    ]
    y = 120
    for label, done, note in checks:
        color = (42, 157, 143) if done else (210, 110, 80)
        mark = "DONE" if done else "GAP"
        draw.rounded_rectangle((36, y, 150, y + 36), radius=6, fill=color)
        draw.text((61, y + 9), mark, fill=(255, 255, 255), font=small)
        draw.text((178, y + 3), label, fill=(25, 25, 25), font=font)
        draw.text((585, y + 7), note, fill=(85, 85, 85), font=small)
        y += 70
    img.save(out / "06_benchmark_completeness.png")


def _plot_before_after_umap(result_dirs: list[Path], out: Path) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    entries = []
    for path in result_dirs:
        fig = path / "03_cell_state_umap.png"
        if fig.exists():
            entries.append((path, fig))
    if not entries:
        return
    width, height = 2200, 1500
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(32)
    font = _pil_font(18)
    small = _pil_font(14)
    draw.text((42, 28), "Before/after cell-state movement by UMAP", fill=(18, 18, 18), font=title_font)
    draw.text((42, 68), "Each panel asks whether virtual KO cells move away from control and toward the true KO state.", fill=(80, 80, 80), font=font)
    slots = [(42, 135, 1040, 640), (1120, 135, 1040, 640), (42, 830, 1040, 640), (1120, 830, 1040, 640)]
    for (path, fig), (x, y, w, h) in zip(entries[:4], slots):
        draw.rounded_rectangle((x, y, x + w, y + h), radius=10, outline=(220, 220, 220), fill=(250, 250, 250))
        label = _result_title(path)
        draw.text((x + 22, y + 18), label, fill=(25, 25, 25), font=font)
        try:
            panel = Image.open(fig).convert("RGB")
            panel.thumbnail((w - 48, h - 82))
            px = x + (w - panel.width) // 2
            py = y + 58 + (h - 92 - panel.height) // 2
            img.paste(panel, (px, py))
        except Exception:
            draw.text((x + 22, y + 70), f"Could not read {fig.name}", fill=(120, 70, 70), font=small)
    img.save(out / "07_before_after_umap_panel.png")


def _plot_single_double_response_map(result_summary: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    if result_summary.empty:
        return
    rows = result_summary.copy()
    for col in ["roc_auc", "mean_direction_cosine", "mean_distribution_improvement", "improved_fraction", "mean_abs_delta_error"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows["label"] = rows.apply(_display_result_label, axis=1)
    metrics = [
        ("roc_auc", "AUC", 0.0, 1.0, False),
        ("mean_direction_cosine", "Direction", 0.0, 1.0, False),
        ("improved_fraction", "Hit-rate", 0.0, 1.0, False),
        ("mean_distribution_improvement", "Distribution", -0.25, 0.75, False),
        ("mean_abs_delta_error", "MAE", 0.0, max(0.01, float(rows["mean_abs_delta_error"].max(skipna=True) or 0.01)), True),
    ]
    cell_w, cell_h = 210, 82
    width = 1750
    height = max(520, 170 + cell_h * min(len(rows), 8) + 80)
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(28)
    font = _pil_font(16)
    small = _pil_font(12)
    draw.text((40, 26), "Single-KO vs double-KO response map", fill=(18, 18, 18), font=title_font)
    draw.text((40, 62), "Green means stronger evidence. Blue on MAE means lower error. This summarizes biological task difficulty, not only method rank.", fill=(80, 80, 80), font=font)
    left, top = 360, 145
    for j, (_, name, _, _, _) in enumerate(metrics):
        draw.text((left + j * cell_w + 12, top - 36), name, fill=(35, 35, 35), font=font)
    for i, (_, row) in enumerate(rows.head(8).iterrows()):
        y = top + i * cell_h
        draw.text((40, y + 18), _short(row["label"], 35), fill=(25, 25, 25), font=font)
        draw.text((250, y + 20), _short(row.get("task_type", ""), 16), fill=(90, 90, 90), font=small)
        for j, (col, _, vmin, vmax, lower) in enumerate(metrics):
            x = left + j * cell_w
            value = float(row[col]) if pd.notna(row[col]) else np.nan
            frac = 0.0 if not np.isfinite(value) else min(max((value - vmin) / (vmax - vmin + 1e-9), 0), 1)
            frac = 1 - frac if lower else frac
            color = _blend((245, 245, 245), (76, 120, 168) if lower else (42, 157, 143), frac)
            draw.rounded_rectangle((x, y, x + cell_w - 16, y + cell_h - 12), radius=6, fill=color, outline=(255, 255, 255))
            label = "NA" if not np.isfinite(value) else f"{value:.2f}"
            text_color = (255, 255, 255) if frac > 0.58 else (30, 30, 30)
            draw.text((x + 72, y + 27), label, fill=text_color, font=font)
    img.save(out / "08_single_double_response_map.png")


def _plot_peak_locus_track(result_dirs: list[Path], out: Path) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    peak_rows = _collect_peak_deltas(result_dirs)
    if not peak_rows:
        return
    width, height = 1900, 940
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(28)
    font = _pil_font(15)
    small = _pil_font(12)
    draw.text((38, 28), "Peak locus track: true vs virtual ATAC response", fill=(18, 18, 18), font=title_font)
    draw.text((38, 64), "Selected peak-level features are placed by genomic coordinate when parsable; bar height shows KO accessibility delta.", fill=(80, 80, 80), font=font)
    chrom_groups: dict[str, list[dict]] = {}
    for row in peak_rows[:46]:
        chrom_groups.setdefault(row["chrom"], []).append(row)
    y = 180
    for chrom, rows in list(chrom_groups.items())[:6]:
        rows = sorted(rows, key=lambda r: r["start"])
        min_pos = min(r["start"] for r in rows)
        max_pos = max(r["end"] for r in rows)
        span = max(max_pos - min_pos, 1)
        left, right = 210, 1780
        draw.text((38, y - 8), chrom, fill=(25, 25, 25), font=font)
        draw.line((left, y, right, y), fill=(160, 160, 160), width=2)
        draw.text((left, y + 18), f"{min_pos:,}", fill=(85, 85, 85), font=small)
        draw.text((right - 82, y + 18), f"{max_pos:,}", fill=(85, 85, 85), font=small)
        label_candidates = sorted(rows, key=lambda item: (item["gene"].upper() == "KDM6A", abs(item["true"])), reverse=True)[:2]
        label_ids = {id(item): idx for idx, item in enumerate(label_candidates)}
        for row in rows:
            x = left + int((row["start"] - min_pos) / span * (right - left))
            true_h = int(max(min(abs(row["true"]) * 430, 95), 3))
            pred_h = int(max(min(abs(row["pred"]) * 430, 95), 3))
            true_color = (190, 70, 70) if row["true"] >= 0 else (70, 90, 180)
            pred_color = (230, 150, 110) if row["pred"] >= 0 else (120, 145, 210)
            draw.rectangle((x - 5, y - true_h if row["true"] >= 0 else y, x + 1, y if row["true"] >= 0 else y + true_h), fill=true_color)
            draw.rectangle((x + 2, y - pred_h if row["pred"] >= 0 else y, x + 8, y if row["pred"] >= 0 else y + pred_h), fill=pred_color)
            if id(row) in label_ids:
                offset = [-58, 46][label_ids[id(row)] % 2]
                label_y = y + offset
                if x - left < 130:
                    label_x = x + 28
                elif right - x < 170:
                    label_x = x - 170
                else:
                    label_x = x - 76
                label_x = max(left + 8, min(label_x, right - 170))
                draw.line((x, y - 5 if row["true"] >= 0 else y + 5, label_x + 8, label_y + 14), fill=(155, 155, 155), width=1)
                draw.text((label_x, label_y), _short(row["gene"] + " " + row["type"], 19), fill=(35, 35, 35), font=small)
        y += 118
    draw.rectangle((1285, 805, 1315, 825), fill=(190, 70, 70))
    draw.text((1326, 803), "true KO delta +", fill=(35, 35, 35), font=small)
    draw.rectangle((1285, 835, 1315, 855), fill=(230, 150, 110))
    draw.text((1326, 833), "virtual KO delta +", fill=(35, 35, 35), font=small)
    draw.rectangle((1510, 805, 1540, 825), fill=(70, 90, 180))
    draw.text((1552, 803), "true KO delta -", fill=(35, 35, 35), font=small)
    draw.rectangle((1510, 835, 1540, 855), fill=(120, 145, 210))
    draw.text((1552, 833), "virtual KO delta -", fill=(35, 35, 35), font=small)
    img.save(out / "09_peak_locus_track.png")


def _plot_method_radar(metrics: pd.DataFrame, availability: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    if metrics.empty:
        return
    rows = metrics.copy()
    for col in ["roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows = rows.loc[rows["method"].astype(str).isin(["VKXAdaptive", "ResponseBoosted", "ConstrainedEnsemble", "CalibratedEnsemble", "PLS", "Ridge", "VKX"])].copy()
    if rows.empty:
        return
    max_mae = max(float(rows["mae"].max(skipna=True) or 0.01), 0.01)
    axes = [
        ("AUC", "roc_auc", 0, 1, False),
        ("Direction", "direction_cosine", 0, 1, False),
        ("R2", "r2", -0.5, 1, False),
        ("Hit-rate", "feature_hit_rate", 0, 1, False),
        ("low MAE", "mae", 0, max_mae, True),
    ]
    width, height = 1500, 1050
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(28)
    font = _pil_font(15)
    small = _pil_font(12)
    draw.text((38, 30), "Method comparison radar and leaderboard", fill=(18, 18, 18), font=title_font)
    draw.text((38, 66), "Radar normalizes heterogeneous metrics so method trade-offs are visible; exact values remain in the leaderboard CSV.", fill=(80, 80, 80), font=font)
    cx, cy, radius = 500, 520, 310
    for frac in [0.25, 0.5, 0.75, 1.0]:
        points = [_radar_point(cx, cy, radius * frac, k, len(axes)) for k in range(len(axes))]
        draw.line(points + [points[0]], fill=(225, 225, 225), width=1)
    for k, (name, *_rest) in enumerate(axes):
        x, y = _radar_point(cx, cy, radius + 38, k, len(axes))
        draw.text((x - 36, y - 9), name, fill=(35, 35, 35), font=small)
        draw.line((cx, cy, *_radar_point(cx, cy, radius, k, len(axes))), fill=(230, 230, 230), width=1)
    colors = _method_colors()
    ordered = rows.sort_values(["roc_auc", "r2"], ascending=False).head(6)
    for _, row in ordered.iterrows():
        pts = []
        for k, (_, col, vmin, vmax, lower) in enumerate(axes):
            value = float(row[col]) if pd.notna(row[col]) else np.nan
            frac = 0 if not np.isfinite(value) else min(max((value - vmin) / (vmax - vmin + 1e-9), 0), 1)
            frac = 1 - frac if lower else frac
            pts.append(_radar_point(cx, cy, radius * frac, k, len(axes)))
        color = colors.get(str(row["method"]), (80, 80, 80))
        draw.line(pts + [pts[0]], fill=color, width=3)
        for x, y in pts:
            draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=color)
    legend_x, legend_y = 900, 150
    draw.text((legend_x, legend_y - 44), "Scored methods", fill=(25, 25, 25), font=font)
    for i, (_, row) in enumerate(ordered.iterrows()):
        y = legend_y + i * 52
        color = colors.get(str(row["method"]), (80, 80, 80))
        draw.rectangle((legend_x, y, legend_x + 26, y + 18), fill=color)
        draw.text((legend_x + 40, y - 3), _method_short(row["method"]), fill=(25, 25, 25), font=font)
        draw.text((legend_x + 260, y - 2), f"AUC {_fmt(row['roc_auc'])} | R2 {_fmt(row['r2'])} | MAE {_fmt(row['mae'])}", fill=(75, 75, 75), font=small)
    y = legend_y + 340
    draw.text((legend_x, y), "External method status", fill=(25, 25, 25), font=font)
    y += 38
    for method in EXTERNAL_METHOD_ORDER:
        status = "not_run"
        if not availability.empty and {"method", "status"}.issubset(availability.columns):
            sub = availability.loc[availability["method"].astype(str).str.lower() == method.lower()]
            if not sub.empty:
                status = str(sub["status"].iloc[0])
        color = colors.get(method, (120, 120, 120))
        draw.rectangle((legend_x, y, legend_x + 26, y + 18), fill=color)
        draw.text((legend_x + 40, y - 3), method, fill=(25, 25, 25), font=font)
        draw.text((legend_x + 170, y - 1), status, fill=(130, 70, 70) if status == "not_run" else (45, 110, 75), font=small)
        y += 42
    img.save(out / "10_method_radar_leaderboard.png")


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
        "![Publication main figure](00_publication_main_figure.png)",
        "",
        "![Method leaderboard](01_method_leaderboard.png)",
        "",
        "![AUC ROC curves](02_auc_roc_curves.png)",
        "",
        "![Real vs virtual heatmap](03_real_vs_virtual_method_heatmap.png)",
        "",
        "![Result gallery](04_single_double_multimodal_gallery.png)",
        "",
        "![Adaptive improvement](05_adaptive_improvement.png)",
        "",
        "![Benchmark completeness](06_benchmark_completeness.png)",
        "",
        "![Before/after UMAP](07_before_after_umap_panel.png)",
        "",
        "![Single vs double response map](08_single_double_response_map.png)",
        "",
        "![Peak locus track](09_peak_locus_track.png)",
        "",
        "![Method radar and leaderboard](10_method_radar_leaderboard.png)",
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
        "4. Publication figure extension: the package now also includes before/after UMAP, single-vs-double response map, ATAC peak locus track, and radar/leaderboard panels.",
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


def _plot_publication_main_figure(
    metrics: pd.DataFrame,
    availability: pd.DataFrame,
    roc_points: pd.DataFrame,
    predictions: pd.DataFrame,
    truth: pd.DataFrame,
    result_summary: pd.DataFrame,
    out: Path,
) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    width, height = 2400, 1700
    img = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(img)
    title_font = _pil_font(34)
    panel_font = _pil_font(28)
    font = _pil_font(20)
    small = _pil_font(16)
    tiny = _pil_font(13)
    draw.text((50, 28), "Formal benchmark of VKX virtual knockout", fill=(15, 15, 15), font=title_font)
    draw.text((50, 76), "Same-dataset scoring is separated from external deep-method slots to avoid over-claiming.", fill=(75, 75, 75), font=font)

    panel_a = (50, 135, 1160, 690)
    panel_b = (1240, 135, 2350, 690)
    panel_c = (50, 770, 1530, 1600)
    panel_d = (1610, 770, 2350, 1600)
    _draw_panel_box(draw, panel_a, "A", "Method-level benchmark", panel_font)
    _draw_panel_box(draw, panel_b, "B", "AUC as ROC curve", panel_font)
    _draw_panel_box(draw, panel_c, "C", "Real vs virtual KO response", panel_font)
    _draw_panel_box(draw, panel_d, "D", "Single/double/multimodal coverage", panel_font)

    _draw_main_method_panel(draw, metrics, availability, panel_a, font, small, tiny)
    _draw_main_roc_panel(draw, roc_points, metrics, panel_b, font, small, tiny)
    _draw_main_heatmap_panel(draw, predictions, truth, panel_c, font, small, tiny)
    _draw_main_gallery_panel(draw, result_summary, panel_d, font, small, tiny)
    img.save(out / "00_publication_main_figure.png")


def _draw_panel_box(draw, box: tuple[int, int, int, int], letter: str, title: str, panel_font) -> None:
    x0, y0, x1, y1 = box
    draw.rounded_rectangle((x0, y0, x1, y1), radius=14, outline=(225, 225, 225), width=2, fill=(253, 253, 253))
    draw.text((x0 + 22, y0 + 16), letter, fill=(15, 15, 15), font=panel_font)
    draw.text((x0 + 62, y0 + 18), title, fill=(30, 30, 30), font=panel_font)


def _draw_main_method_panel(draw, metrics: pd.DataFrame, availability: pd.DataFrame, box: tuple[int, int, int, int], font, small, tiny) -> None:
    x0, y0, x1, _ = box
    if metrics.empty:
        draw.text((x0 + 28, y0 + 76), "No benchmark metrics found.", fill=(80, 80, 80), font=font)
        return
    rows = metrics.copy()
    for col in ["roc_auc", "r2", "mae", "direction_cosine"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    rows["_rank"] = rows["roc_auc"].fillna(-1) + 0.2 * rows["r2"].fillna(-1) - 0.2 * rows["mae"].fillna(1)
    rows = rows.sort_values("_rank", ascending=False).head(5)
    draw.text((x0 + 32, y0 + 70), "Scored methods", fill=(80, 80, 80), font=small)
    metric_x = [x0 + 360, x0 + 560, x0 + 760]
    for x, label in zip(metric_x, ["AUC", "R2", "MAE"]):
        draw.text((x, y0 + 105), label, fill=(45, 45, 45), font=tiny)
    colors = _method_colors()
    max_mae = max(float(rows["mae"].max(skipna=True)), 0.01)
    for i, (_, row) in enumerate(rows.iterrows()):
        y = y0 + 138 + i * 66
        method = str(row["method"])
        draw.rectangle((x0 + 32, y + 7, x0 + 46, y + 38), fill=colors.get(method, (80, 80, 80)))
        draw.text((x0 + 58, y + 8), _method_short(method), fill=(25, 25, 25), font=font)
        specs = [("roc_auc", 0, 1, False), ("r2", -0.5, 1, False), ("mae", 0, max_mae, True)]
        for x, (col, vmin, vmax, lower_better) in zip(metric_x, specs):
            val = float(row[col]) if pd.notna(row[col]) else np.nan
            frac = 0 if not np.isfinite(val) else min(max((val - vmin) / (vmax - vmin + 1e-9), 0), 1)
            frac = 1 - frac if lower_better else frac
            bar_color = (76, 120, 168) if lower_better else (42, 157, 143)
            draw.rectangle((x, y + 8, x + 118, y + 32), fill=(242, 242, 242), outline=(220, 220, 220))
            draw.rectangle((x, y + 8, x + int(118 * frac), y + 32), fill=bar_color)
            draw.text((x + 128, y + 10), f"{val:.2f}" if col != "mae" else f"{val:.3f}", fill=(30, 30, 30), font=tiny)
    external = []
    if not availability.empty and "method" in availability.columns:
        for method in EXTERNAL_METHOD_ORDER:
            sub = availability.loc[availability["method"].astype(str).str.lower() == method.lower()]
            if not sub.empty:
                external.append(f"{method}: {sub['status'].iloc[0]}")
    text = "External slots: " + "; ".join(external) if external else "External slots: not provided"
    draw.text((x0 + 32, y0 + 500), text, fill=(105, 65, 65), font=small)


def _draw_main_roc_panel(draw, roc_points: pd.DataFrame, metrics: pd.DataFrame, box: tuple[int, int, int, int], font, small, tiny) -> None:
    x0, y0, x1, _ = box
    if roc_points.empty:
        draw.text((x0 + 28, y0 + 76), "No ROC points found.", fill=(80, 80, 80), font=font)
        return
    auc_map = dict(zip(metrics["method"].astype(str), pd.to_numeric(metrics["roc_auc"], errors="coerce"))) if not metrics.empty else {}
    left, top, size = x0 + 70, y0 + 105, 420
    draw.rectangle((left, top, left + size, top + size), outline=(180, 180, 180))
    for tick in [0, 0.5, 1.0]:
        tx = left + int(tick * size)
        ty = top + size - int(tick * size)
        draw.line((tx, top, tx, top + size), fill=(235, 235, 235))
        draw.line((left, ty, left + size, ty), fill=(235, 235, 235))
        draw.text((tx - 10, top + size + 12), f"{tick:g}", fill=(60, 60, 60), font=tiny)
        draw.text((left - 28, ty - 8), f"{tick:g}", fill=(60, 60, 60), font=tiny)
    draw.line((left, top + size, left + size, top), fill=(165, 165, 165))
    colors = _method_colors()
    order = ["VKXAdaptive", "ResponseBoosted", "ConstrainedEnsemble", "Ridge", "PLS", "VKX"]
    legend_y = top + 15
    for method in order:
        group = roc_points.loc[roc_points["method"].astype(str) == method]
        if group.empty:
            continue
        curve = group.groupby("fpr", as_index=False)["tpr"].mean().sort_values("fpr")
        points = [(left + int(float(r.fpr) * size), top + size - int(float(r.tpr) * size)) for r in curve.itertuples()]
        color = colors.get(method, (70, 70, 70))
        if len(points) > 1:
            draw.line(points, fill=color, width=4)
        draw.rectangle((x0 + 560, legend_y + 4, x0 + 582, legend_y + 22), fill=color)
        auc = auc_map.get(method, np.nan)
        label = f"{_method_short(method)} AUC {auc:.3f}" if np.isfinite(auc) else _method_short(method)
        draw.text((x0 + 596, legend_y), label, fill=(25, 25, 25), font=small)
        legend_y += 38
    draw.text((left + 160, top + size + 48), "False positive rate", fill=(40, 40, 40), font=small)
    draw.text((left - 8, top - 30), "TPR", fill=(40, 40, 40), font=small)


def _draw_main_heatmap_panel(draw, predictions: pd.DataFrame, truth: pd.DataFrame, box: tuple[int, int, int, int], font, small, tiny) -> None:
    if predictions.empty or truth.empty:
        draw.text((box[0] + 28, box[1] + 76), "No prediction/truth table found.", fill=(80, 80, 80), font=font)
        return
    merged = predictions.merge(truth, on="ko_target", how="inner")
    if merged.empty:
        return
    features = [col.removeprefix("true_delta_") for col in truth.columns if col.startswith("true_delta_")]
    scored = []
    for feature in features:
        scored.append((feature, float(np.nanmean(np.abs(pd.to_numeric(truth[f"true_delta_{feature}"], errors="coerce"))))))
    chosen = [feature for feature, _ in sorted(scored, key=lambda x: x[1], reverse=True)[:7]]
    rows = []
    for _, row in truth.iterrows():
        rows.append((f"TRUE {row['ko_target']}", [float(row.get(f"true_delta_{f}", np.nan)) for f in chosen]))
    boosted = merged.loc[merged["method"].astype(str) == "ResponseBoosted"]
    for _, row in boosted.iterrows():
        rows.append((f"VKX-B {row['ko_target']}", [float(row.get(f"pred_delta_{f}", np.nan)) for f in chosen]))
    x0, y0, _, _ = box
    left = x0 + 255
    top = y0 + 105
    cell_w = 130
    cell_h = 48
    matrix = np.asarray([vals for _, vals in rows], dtype=float)
    vmax = max(float(np.nanmax(np.abs(matrix))), 1e-9)
    for j, feature in enumerate(chosen):
        label = _feature_label(feature)
        for k, part in enumerate(label.split("\n")[:2]):
            draw.text((left + j * cell_w, y0 + 68 + 17 * k), _short(part, 16), fill=(35, 35, 35), font=tiny)
    for i, (label, values) in enumerate(rows):
        y = top + i * cell_h
        draw.text((x0 + 32, y + 10), _short(label, 21), fill=(25, 25, 25), font=small)
        for j, value in enumerate(values):
            x = left + j * cell_w
            draw.rectangle((x, y, x + cell_w - 5, y + cell_h - 5), fill=_delta_color(value, vmax), outline=(255, 255, 255))
            text_color = (255, 255, 255) if np.isfinite(value) and abs(value) / vmax > 0.45 else (35, 35, 35)
            draw.text((x + 28, y + 11), f"{value:.2f}", fill=text_color, font=tiny)
    draw.text((x0 + 32, y0 + 690), "Rows compare true KO deltas with VKX-Boosted predictions on the strongest response features.", fill=(80, 80, 80), font=small)


def _draw_main_gallery_panel(draw, result_summary: pd.DataFrame, box: tuple[int, int, int, int], font, small, tiny) -> None:
    x0, y0, _, _ = box
    if result_summary.empty:
        draw.text((x0 + 28, y0 + 76), "No result summaries found.", fill=(80, 80, 80), font=font)
        return
    rows = result_summary.copy()
    rows["display_label"] = rows.apply(_display_result_label, axis=1)
    for col in ["roc_auc", "mean_direction_cosine", "mean_abs_delta_error"]:
        rows[col] = pd.to_numeric(rows[col], errors="coerce")
    draw.text((x0 + 32, y0 + 72), "Dataset / task", fill=(70, 70, 70), font=small)
    draw.text((x0 + 355, y0 + 72), "AUC", fill=(70, 70, 70), font=small)
    draw.text((x0 + 520, y0 + 72), "Direction", fill=(70, 70, 70), font=small)
    max_mae = max(float(rows["mean_abs_delta_error"].max(skipna=True)), 0.01)
    for i, (_, row) in enumerate(rows.head(6).iterrows()):
        y = y0 + 112 + i * 92
        draw.text((x0 + 32, y), str(row["display_label"]), fill=(25, 25, 25), font=small)
        for x, col, vmin, vmax, lower in [
            (x0 + 355, "roc_auc", 0, 1, False),
            (x0 + 520, "mean_direction_cosine", 0, 1, False),
            (x0 + 355, "mean_abs_delta_error", 0, max_mae, True),
        ]:
            yy = y + (36 if col == "mean_abs_delta_error" else 0)
            value = float(row[col]) if pd.notna(row[col]) else np.nan
            frac = 0 if not np.isfinite(value) else min(max((value - vmin) / (vmax - vmin + 1e-9), 0), 1)
            frac = 1 - frac if lower else frac
            color = (76, 120, 168) if lower else (42, 157, 143)
            draw.rectangle((x, yy, x + 108, yy + 22), fill=(244, 244, 244), outline=(220, 220, 220))
            draw.rectangle((x, yy, x + int(108 * frac), yy + 22), fill=color)
            label = "MAE" if col == "mean_abs_delta_error" else ""
            draw.text((x + 116, yy + 2), f"{label} {value:.2f}".strip(), fill=(35, 35, 35), font=tiny)
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
        "VKXAdaptive": (38, 70, 83),
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
        "VKXAdaptive": "VKX-Adaptive",
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


def _has_task(summary: pd.DataFrame, keyword: str) -> bool:
    if summary.empty:
        return False
    fields = []
    for col in ["task_type", "input_modality", "state_representation", "dataset"]:
        if col in summary.columns:
            fields.extend(summary[col].astype(str).str.lower().tolist())
    return keyword.lower() in " ".join(fields)


def _external_scored(availability: pd.DataFrame) -> bool:
    if availability.empty or "method" not in availability.columns or "status" not in availability.columns:
        return False
    sub = availability.loc[availability["method"].astype(str).isin(EXTERNAL_METHOD_ORDER)]
    return bool((sub["status"].astype(str) == "provided").any())


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


def _result_title(path: Path) -> str:
    text = path.name.lower()
    if "single" in text or "papalexi" in text:
        return "Single KO: RNA + ADT"
    if "hmpcite" in text:
        return "Double KO: RNA + ADT"
    if "double" in text or "norman" in text:
        return "Double KO: RNA programs"
    if "atac" in text or "peak" in text:
        return "ATAC peak-level KO"
    return path.name.replace("_", " ")


def _collect_peak_deltas(result_dirs: list[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in result_dirs:
        table_path = path / "delta_table.csv"
        if not table_path.exists():
            continue
        try:
            table = pd.read_csv(table_path, nrows=1)
        except Exception:
            continue
        if table.empty:
            continue
        row = table.iloc[0]
        peaks = []
        for col in table.columns:
            if not col.startswith("true_delta_peak_"):
                continue
            feature = col.removeprefix("true_delta_")
            pred_col = f"pred_delta_{feature}"
            if pred_col not in table.columns:
                continue
            parsed = _parse_peak_feature(feature)
            if parsed is None:
                continue
            true = pd.to_numeric(pd.Series([row.get(col)]), errors="coerce").iloc[0]
            pred = pd.to_numeric(pd.Series([row.get(pred_col)]), errors="coerce").iloc[0]
            if not np.isfinite(true) or not np.isfinite(pred):
                continue
            parsed["true"] = float(true)
            parsed["pred"] = float(pred)
            parsed["score"] = abs(float(true)) + 0.35 * abs(float(pred))
            peaks.append(parsed)
        rows.extend(sorted(peaks, key=lambda item: item["score"], reverse=True)[:48])
    return sorted(rows, key=lambda item: item["score"], reverse=True)


def _parse_peak_feature(feature: str) -> dict | None:
    parts = feature.split("_")
    if len(parts) < 6 or parts[0] != "peak":
        return None
    chrom = parts[1]
    try:
        start = int(parts[2])
        end = int(parts[3])
    except Exception:
        return None
    peak_type = parts[-1]
    gene = "_".join(parts[4:-1]) or "peak"
    return {"chrom": chrom, "start": start, "end": end, "gene": gene, "type": peak_type}


def _radar_point(cx: int, cy: int, radius: float, idx: int, n: int) -> tuple[int, int]:
    angle = -np.pi / 2 + 2 * np.pi * idx / n
    return int(cx + radius * np.cos(angle)), int(cy + radius * np.sin(angle))


def _blend(low: tuple[int, int, int], high: tuple[int, int, int], frac: float) -> tuple[int, int, int]:
    frac = min(max(float(frac), 0.0), 1.0)
    return tuple(int(low[i] + (high[i] - low[i]) * frac) for i in range(3))


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


def _fmt(value: object) -> str:
    try:
        return f"{float(value):.3f}"
    except Exception:
        return "NA"


def _pil_font(size: int):
    from PIL import ImageFont

    for name in ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()
