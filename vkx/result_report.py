from __future__ import annotations

from pathlib import Path

import pandas as pd

from .cards import make_ko_summary_cards
from .diagnostics import diagnose_virtual_ko_results
from .figure_pack import make_figure_package


def summarize_result_directory(result_dir: str | Path, out_dir: str | Path | None = None) -> dict[str, Path | pd.DataFrame]:
    result = Path(result_dir)
    out = Path(out_dir) if out_dir else result / "readable_result_report"
    out.mkdir(parents=True, exist_ok=True)

    delta_csv, mode = _find_delta_csv(result)
    manifest_csv = _first_existing(result, ["derived_state_manifest.csv"])
    auc_csv = _first_existing(result, ["auc_summary.csv"])
    confidence_csv = _first_existing(result, ["transfer_confidence.csv"])
    analysis_mode = _read_analysis_mode(result, fallback=mode)

    ko_cards = pd.DataFrame()
    diagnosis_ko = pd.DataFrame()
    diagnosis_feature = pd.DataFrame()
    if delta_csv is not None:
        ko_cards = make_ko_summary_cards(
            delta_csv=delta_csv,
            out_dir=out / "ko_cards",
            auc_csv=auc_csv,
            confidence_csv=confidence_csv,
        )
        diagnosis = diagnose_virtual_ko_results(
            delta_csv=delta_csv,
            manifest_csv=manifest_csv,
            confidence_csv=confidence_csv,
            out_dir=out / "diagnosis",
        )
        diagnosis_ko = diagnosis["ko"]
        diagnosis_feature = diagnosis["feature"]

    figure_index = make_figure_package(result, out / "figure_package")
    report_path = out / "user_readable_result_report.md"
    _write_user_report(
        report_path=report_path,
        result_dir=result,
        output_dir=out,
        analysis_mode=analysis_mode,
        delta_csv=delta_csv,
        ko_cards=ko_cards,
        diagnosis_ko=diagnosis_ko,
        diagnosis_feature=diagnosis_feature,
        figure_index=figure_index,
    )
    return {
        "report": report_path,
        "figure_index": figure_index,
        "ko_cards": ko_cards,
        "diagnosis_ko": diagnosis_ko,
        "diagnosis_feature": diagnosis_feature,
    }


def _find_delta_csv(result: Path) -> tuple[Path | None, str]:
    candidates = [
        ("delta_table.csv", "evaluation"),
        ("predicted_ko_delta.csv", "prediction_only"),
    ]
    for name, mode in candidates:
        path = result / name
        if path.exists():
            return path, mode
    return None, "unknown"


def _first_existing(result: Path, names: list[str]) -> Path | None:
    for name in names:
        path = result / name
        if path.exists():
            return path
    return None


def _read_analysis_mode(result: Path, fallback: str) -> str:
    path = result / "analysis_mode.md"
    if not path.exists():
        return fallback
    text = path.read_text(encoding="utf-8", errors="ignore")
    for line in text.splitlines():
        if line.strip().startswith("- Mode:"):
            return line.split(":", 1)[1].strip().strip("`")
    return fallback


def _write_user_report(
    report_path: Path,
    result_dir: Path,
    output_dir: Path,
    analysis_mode: str,
    delta_csv: Path | None,
    ko_cards: pd.DataFrame,
    diagnosis_ko: pd.DataFrame,
    diagnosis_feature: pd.DataFrame,
    figure_index: pd.DataFrame,
) -> None:
    lines = [
        "# VKX User-readable Result Report",
        "",
        f"Source result directory: `{result_dir}`",
        f"Report directory: `{output_dir}`",
        "",
        "## 1. What kind of result is this?",
        "",
        _mode_explanation(analysis_mode),
        "",
        "## 2. Bottom-line diagnosis",
        "",
    ]
    if delta_csv is None:
        lines.append("No `delta_table.csv` or `predicted_ko_delta.csv` was found, so this report only packaged existing figures.")
    elif diagnosis_ko.empty:
        lines.append("Diagnosis could not be computed, but KO cards and figures were still collected when possible.")
    else:
        lines.extend(_diagnosis_summary_lines(diagnosis_ko))
    lines.extend(["", "## 3. KO summary cards", ""])
    if ko_cards.empty:
        lines.append("No KO cards were generated.")
    else:
        lines.append("Each KO card lists the most changed state features and, when real KO labels exist, the virtual-vs-real agreement.")
        lines.append("")
        for _, row in ko_cards.head(20).iterrows():
            fig = str(row.get("card_figure", ""))
            csv_path = str(row.get("card_csv", ""))
            lines.append(f"### {row['ko_target']}")
            lines.append("")
            if fig:
                rel = _relative(report_path.parent, Path(fig))
                lines.append(f"![KO card]({rel})")
                lines.append("")
            if csv_path:
                lines.append(f"Card table: `{_relative(report_path.parent, Path(csv_path))}`")
                lines.append("")
    lines.extend(["", "## 4. Reliability diagnosis figures", ""])
    for name in ["01_failure_diagnosis_overview.png", "02_feature_error_heatmap.png"]:
        path = output_dir / "diagnosis" / name
        if path.exists():
            lines.append(f"### {name}")
            lines.append("")
            lines.append(f"![{name}]({_relative(report_path.parent, path)})")
            lines.append("")
    lines.extend(["", "## 5. All collected figures", ""])
    if figure_index.empty:
        lines.append("No PNG figures were found.")
    else:
        lines.append("A complete figure package was generated here:")
        lines.append("")
        lines.append(f"- `{_relative(report_path.parent, output_dir / 'figure_package' / 'figure_package_report.md')}`")
        lines.append(f"- `{_relative(report_path.parent, output_dir / 'figure_package' / 'figure_index.csv')}`")
    lines.extend(["", "## 6. Files generated by this summary command", ""])
    lines.extend(
        [
            "- `ko_cards/`: one readable card per KO target.",
            "- `diagnosis/`: KO-level and feature-level reliability diagnosis.",
            "- `figure_package/`: copied PNG figures with explanations.",
            "- `user_readable_result_report.md`: this integrated report.",
        ]
    )
    if not diagnosis_feature.empty:
        high = diagnosis_feature.loc[diagnosis_feature["severity"] == "high_risk"]
        if not high.empty:
            lines.extend(["", "## 7. Highest-risk feature examples", ""])
            cols = ["ko_target", "feature", "modality", "true_delta", "pred_delta", "abs_error", "issue_tags"]
            lines.append(_markdown_table(high[cols].head(20)))
    report_path.write_text("\n".join(lines), encoding="utf-8")


def _mode_explanation(mode: str) -> str:
    if mode == "evaluation":
        return "This is a labelled KO evaluation. Real KO labels are available, so heatmaps, ROC/AUC, MAE/R2-style metrics, and diagnosis flags can be interpreted as accuracy evidence."
    if mode == "prediction_only":
        return "This is a prediction-only application. The input does not contain matching real KO labels, so the report shows predicted state shifts, transfer confidence, prior coverage, and uncertainty, but not true accuracy."
    if mode == "double_ko_evaluation":
        return "This is a labelled double-KO evaluation. Additive and interaction residual predictions can be compared against real double-KO effects."
    if mode == "state_scoring_only":
        return "This result only converts cells into pathway/program/protein state scores. It does not evaluate or apply a KO model yet."
    return "The analysis mode could not be inferred. Interpret accuracy metrics only if real KO labels are present."


def _diagnosis_summary_lines(diagnosis_ko: pd.DataFrame) -> list[str]:
    lines = []
    for _, row in diagnosis_ko.head(20).iterrows():
        ko = row.get("ko_target", "unknown")
        mode = row.get("analysis_mode", "unknown")
        high = int(row.get("n_high_risk_features", 0))
        warn = int(row.get("n_warning_features", 0))
        match = row.get("direction_match_fraction", pd.NA)
        match_text = "n/a" if pd.isna(match) else f"{float(match):.2f}"
        readout = row.get("plain_readout", "")
        lines.append(f"- `{ko}` ({mode}): direction match {match_text}; {high} high-risk and {warn} warning features. {readout}")
    return lines


def _relative(base: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(base.resolve()).as_posix()
    except Exception:
        try:
            return path.resolve().as_posix()
        except Exception:
            return str(path)


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return ""
    clean = frame.copy()
    for col in clean.columns:
        clean[col] = clean[col].map(_format_cell)
    lines = []
    headers = [str(col) for col in clean.columns]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for _, row in clean.iterrows():
        lines.append("| " + " | ".join(str(row[col]) for col in clean.columns) + " |")
    return "\n".join(lines)


def _format_cell(value) -> str:
    if pd.isna(value):
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    text = str(value).replace("|", "/")
    return text if len(text) <= 80 else text[:77] + "..."
