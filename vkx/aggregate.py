from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def aggregate_benchmarks(result_dirs: list[str | Path], out_dir: str | Path) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    dirs = _expand_result_dirs(result_dirs)
    formal = _collect_formal_metrics(dirs)
    generators = _collect_generator_metrics(dirs)
    best = _best_formal_methods(formal)
    generator_summary = _generator_summary(generators)

    formal.to_csv(out / "formal_method_metrics_aggregate.csv", index=False)
    best.to_csv(out / "formal_best_methods.csv", index=False)
    generators.to_csv(out / "generator_metrics_aggregate.csv", index=False)
    generator_summary.to_csv(out / "generator_summary.csv", index=False)
    _plot_formal_leaderboard(best, out)
    _plot_generator_summary(generator_summary, out)
    _write_report(out, dirs, formal, best, generators, generator_summary)
    return {
        "formal": formal,
        "best": best,
        "generators": generators,
        "generator_summary": generator_summary,
    }


def _expand_result_dirs(result_dirs: list[str | Path]) -> list[Path]:
    expanded: list[Path] = []
    for item in result_dirs:
        path = Path(item)
        if not path.exists():
            continue
        if path.is_file():
            path = path.parent
        if _has_supported_result(path):
            expanded.append(path)
        else:
            for child in sorted(path.iterdir()):
                if child.is_dir() and _has_supported_result(child):
                    expanded.append(child)
    seen = set()
    unique = []
    for path in expanded:
        resolved = path.resolve()
        if resolved not in seen:
            unique.append(path)
            seen.add(resolved)
    return unique


def _has_supported_result(path: Path) -> bool:
    return (path / "method_metric_comparison.csv").exists() or (path / "hard_generator_metrics.csv").exists()


def _collect_formal_metrics(dirs: list[Path]) -> pd.DataFrame:
    rows = []
    for path in dirs:
        csv = path / "method_metric_comparison.csv"
        if not csv.exists():
            continue
        table = pd.read_csv(csv)
        for _, row in table.iterrows():
            out = row.to_dict()
            out["benchmark_id"] = path.name
            out["result_dir"] = str(path)
            out["result_type"] = "formal_benchmark"
            rows.append(out)
    if not rows:
        return pd.DataFrame(columns=["benchmark_id", "method", "roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate", "result_dir", "result_type"])
    table = pd.DataFrame(rows)
    for col in ["roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate"]:
        if col in table.columns:
            table[col] = pd.to_numeric(table[col], errors="coerce")
    return table


def _best_formal_methods(formal: pd.DataFrame) -> pd.DataFrame:
    if formal.empty:
        return pd.DataFrame(columns=["benchmark_id", "best_method", "roc_auc", "direction_cosine", "r2", "mae", "feature_hit_rate", "result_dir"])
    rows = []
    for benchmark_id, group in formal.groupby("benchmark_id", observed=True):
        scored = group.copy()
        scored["_rank_score"] = (
            scored.get("roc_auc", np.nan).fillna(-1)
            + 0.25 * scored.get("direction_cosine", np.nan).fillna(-1)
            + 0.15 * scored.get("r2", np.nan).fillna(-1)
            - 0.25 * scored.get("mae", np.nan).fillna(1)
        )
        row = scored.sort_values("_rank_score", ascending=False).iloc[0]
        rows.append(
            {
                "benchmark_id": benchmark_id,
                "best_method": row.get("method"),
                "roc_auc": row.get("roc_auc", np.nan),
                "direction_cosine": row.get("direction_cosine", np.nan),
                "r2": row.get("r2", np.nan),
                "mae": row.get("mae", np.nan),
                "feature_hit_rate": row.get("feature_hit_rate", np.nan),
                "result_dir": row.get("result_dir"),
            }
        )
    return pd.DataFrame(rows).sort_values(["roc_auc", "r2"], ascending=False)


def _collect_generator_metrics(dirs: list[Path]) -> pd.DataFrame:
    rows = []
    for path in dirs:
        csv = path / "hard_generator_metrics.csv"
        if not csv.exists():
            continue
        table = pd.read_csv(csv)
        anchor = _read_anchor_method(path / "hard_generator_report.md")
        for _, row in table.iterrows():
            out = row.to_dict()
            out["generator_id"] = path.name
            out["anchor_method"] = anchor
            out["result_dir"] = str(path)
            out["result_type"] = "hard_generator"
            rows.append(out)
    if not rows:
        return pd.DataFrame(columns=["generator_id", "ko_target", "mae", "direction_cosine", "r2", "anchor_method", "result_dir", "result_type"])
    table = pd.DataFrame(rows)
    for col in ["mae", "direction_cosine", "r2"]:
        if col in table.columns:
            table[col] = pd.to_numeric(table[col], errors="coerce")
    return table


def _read_anchor_method(report: Path) -> str:
    if not report.exists():
        return "unknown"
    for line in report.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.strip().lower().startswith("- anchor method:"):
            return line.split(":", 1)[1].strip().strip("`")
    return "unknown"


def _generator_summary(generators: pd.DataFrame) -> pd.DataFrame:
    if generators.empty:
        return pd.DataFrame(columns=["generator_id", "anchor_method", "n_ko", "mean_mae", "mean_direction_cosine", "mean_r2", "result_dir"])
    summary = (
        generators.groupby(["generator_id", "anchor_method", "result_dir"], observed=True)
        .agg(n_ko=("ko_target", "nunique"), mean_mae=("mae", "mean"), mean_direction_cosine=("direction_cosine", "mean"), mean_r2=("r2", "mean"))
        .reset_index()
    )
    return summary.sort_values(["mean_r2", "mean_mae"], ascending=[False, True])


def _plot_formal_leaderboard(best: pd.DataFrame, out: Path) -> None:
    if best.empty:
        return
    rows = best.head(18).copy()
    rows["display_id"] = rows["benchmark_id"].map(_display_result_id)
    _plot_table_bars(
        rows=rows,
        out_path=out / "01_formal_benchmark_leaderboard.png",
        title="Formal benchmark leaderboard",
        subtitle="Best method per result directory; higher AUC/R2 and lower MAE are better.",
        label_col="display_id",
        bar_specs=[
            ("roc_auc", "AUC", 0.0, 1.0, False),
            ("r2", "R2", -0.5, 1.0, False),
            ("mae", "MAE", 0.0, max(0.01, float(rows["mae"].max())), True),
        ],
        extra_col="best_method",
    )


def _plot_generator_summary(summary: pd.DataFrame, out: Path) -> None:
    if summary.empty:
        return
    rows = summary.head(18).copy()
    rows["display_id"] = rows["generator_id"].map(_display_result_id)
    _plot_table_bars(
        rows=rows,
        out_path=out / "02_hard_generator_leaderboard.png",
        title="Hard-constrained generator leaderboard",
        subtitle="Mean cell-level metrics across generated KO targets.",
        label_col="display_id",
        bar_specs=[
            ("mean_r2", "R2", -0.5, 1.0, False),
            ("mean_direction_cosine", "Direction", 0.0, 1.0, False),
            ("mean_mae", "MAE", 0.0, max(0.01, float(rows["mean_mae"].max())), True),
        ],
        extra_col="anchor_method",
    )


def _plot_table_bars(
    rows: pd.DataFrame,
    out_path: Path,
    title: str,
    subtitle: str,
    label_col: str,
    bar_specs: list[tuple[str, str, float, float, bool]],
    extra_col: str,
) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    font = _pil_font(15)
    small = _pil_font(12)
    title_font = _pil_font(24)
    row_h = 58
    width = 1500
    height = 130 + row_h * len(rows) + 70
    img = Image.new("RGB", (width, max(420, height)), "white")
    draw = ImageDraw.Draw(img)
    draw.text((35, 25), title, fill=(20, 20, 20), font=title_font)
    draw.text((35, 60), subtitle, fill=(80, 80, 80), font=font)
    x_label = 35
    x_extra = 430
    panel_x = [650, 910, 1170]
    for i, (_, row) in enumerate(rows.iterrows()):
        y = 110 + i * row_h
        draw.text((x_label, y + 8), _short(row.get(label_col, ""), 45), fill=(30, 30, 30), font=font)
        draw.text((x_extra, y + 8), _short(row.get(extra_col, ""), 24), fill=(80, 80, 80), font=font)
        for x0, spec in zip(panel_x, bar_specs):
            col, name, vmin, vmax, lower_better = spec
            value = float(row.get(col, np.nan))
            if i == 0:
                draw.text((x0, 88), name, fill=(30, 30, 30), font=small)
            if not np.isfinite(value):
                frac = 0.0
                label = "NA"
            else:
                raw = min(max((value - vmin) / (vmax - vmin + 1e-9), 0.0), 1.0)
                frac = 1.0 - raw if lower_better else raw
                label = f"{value:.3f}" if abs(value) < 1 else f"{value:.2f}"
            draw.rectangle((x0, y + 8, x0 + 180, y + 30), fill=(245, 245, 245), outline=(220, 220, 220))
            color = (42, 157, 143) if not lower_better else (76, 120, 168)
            draw.rectangle((x0, y + 8, x0 + int(180 * frac), y + 30), fill=color)
            draw.text((x0 + 190, y + 10), label, fill=(30, 30, 30), font=small)
    img.save(out_path)


def _write_report(
    out: Path,
    dirs: list[Path],
    formal: pd.DataFrame,
    best: pd.DataFrame,
    generators: pd.DataFrame,
    generator_summary: pd.DataFrame,
) -> None:
    lines = [
        "# Aggregate Benchmark Report",
        "",
        f"Scanned result directories: {len(dirs)}",
        f"Formal benchmark metric rows: {len(formal)}",
        f"Hard-generator metric rows: {len(generators)}",
        "",
        "## Formal Benchmark Best Methods",
        "",
        "```text",
        best.head(20).to_string(index=False) if not best.empty else "No formal benchmark metrics were found.",
        "```",
        "",
        "## Hard Generator Summary",
        "",
        "```text",
        generator_summary.head(20).to_string(index=False) if not generator_summary.empty else "No hard-generator metrics were found.",
        "```",
        "",
        "## Figures",
        "",
        "![Formal benchmark leaderboard](01_formal_benchmark_leaderboard.png)",
        "",
        "![Hard-constrained generator leaderboard](02_hard_generator_leaderboard.png)",
        "",
        "## Interpretation",
        "",
        "This aggregate report is a bookkeeping layer. It does not replace dataset-specific heatmaps, ROC curves, or cell-state plots.",
        "Use it to decide which result directories deserve deeper inspection and which anchor or method should be promoted as the current default.",
    ]
    (out / "aggregate_benchmark_report.md").write_text("\n".join(lines), encoding="utf-8")


def _short(value: object, max_len: int = 40) -> str:
    text = str(value)
    return text if len(text) <= max_len else text[: max_len - 3] + "..."


def _display_result_id(value: object) -> str:
    text = str(value)
    replacements = [
        ("formal_benchmark_", "formal: "),
        ("hard_generator_", "generator: "),
        ("papalexi_hmpcite_", ""),
        ("papalexi_", ""),
        ("_demo", ""),
        ("_adaptive_boosted", " adaptive boosted"),
        ("_boosted", " boosted"),
        ("_ensemble", " ensemble"),
        ("_calibrated", " calibrated"),
    ]
    for old, new in replacements:
        text = text.replace(old, new)
    return text.replace("_", " ")


def _pil_font(size: int):
    from PIL import ImageFont

    for name in ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf"]:
        try:
            return ImageFont.truetype(name, size=size)
        except Exception:
            continue
    return ImageFont.load_default()
