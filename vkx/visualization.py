from __future__ import annotations

from pathlib import Path
import textwrap
import re

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc, roc_curve
from sklearn.preprocessing import StandardScaler
import umap


PALETTE = {
    "control cells": "#BDBDBD",
    "virtual KO cells": "#E76F51",
    "true KO cells": "#2A9D8F",
}


def setup_plot() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"]:
        if font in available_fonts:
            plt.rcParams["font.family"] = font
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["axes.titlesize"] = 12
    plt.rcParams["axes.labelsize"] = 11
    plt.rcParams["xtick.labelsize"] = 9
    plt.rcParams["ytick.labelsize"] = 9
    plt.rcParams["legend.fontsize"] = 10


def short_feature(feature: str) -> str:
    for prefix in ["pathway_", "protein_", "program_", "atac_", "tf_"]:
        if feature.startswith(prefix):
            feature = feature.replace(prefix, "")
    return feature[:36]


def wrap_label(text: str, width: int = 16) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False, break_on_hyphens=False))


def plot_summary(summary: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    fig = plt.figure(figsize=(14.5, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.45, 0.9])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    row = summary.iloc[0]
    values = pd.DataFrame(
        [
            {"metric": "Distribution\nimprovement", "value": row.get("mean_distribution_improvement", np.nan), "baseline": 0.0},
            {"metric": "Improved\nfeatures", "value": row.get("improved_fraction", np.nan), "baseline": 0.0},
            {"metric": "Direction\ncosine", "value": row.get("mean_direction_cosine", np.nan), "baseline": 0.0},
        ]
    )
    colors = ["#4E79A7" if value >= 0 else "#C65D4B" for value in values["value"]]
    axes[0].bar([wrap_label(x, 13) for x in values["metric"]], values["value"], color=colors, width=0.62)
    axes[0].axhline(0, color="0.25", linewidth=1)
    axes[0].set_title("Performance overview", pad=8)
    axes[0].set_ylabel("score")
    for tick in axes[0].get_xticklabels():
        tick.set_rotation(0)
    for i, value in enumerate(values["value"]):
        axes[0].text(i, value + (0.035 if value >= 0 else -0.035), f"{value:.2f}", ha="center", va="bottom" if value >= 0 else "top", fontsize=10)

    axes[1].axis("off")
    verdict = "Good signal" if row.get("mean_distribution_improvement", 0) > 0 and row.get("mean_direction_cosine", 0) > 0.5 else "Mixed signal"
    text = (
        f"Dataset: {row['dataset']}\n"
        f"State: {row['state_representation']}\n"
        f"KO targets: {int(row['n_ko'])}\n"
        f"Features: {int(row['n_features'])}\n"
        f"Calibration: {row.get('calibration_method', 'none')}\n\n"
        f"Verdict: {verdict}"
    )
    axes[1].text(0.02, 0.92, text, va="top", fontsize=11, linespacing=1.35, transform=axes[1].transAxes)

    mae = row.get("mean_abs_delta_error", np.nan)
    axes[2].bar(["Magnitude\nerror"], [mae], color="#F28E2B", width=0.55)
    axes[2].set_title("Magnitude error", pad=8)
    axes[2].set_ylabel("mean abs. error")
    axes[2].text(0, mae, f"{mae:.3f}", ha="center", va="bottom", fontsize=10)
    fig.suptitle("Virtual KO Result Summary", fontsize=16)
    fig.savefig(out_dir / "01_summary_dashboard.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_true_virtual_heatmap(delta_table: pd.DataFrame, out_dir: Path, max_features: int = 10) -> None:
    setup_plot()
    if delta_table.empty:
        return
    feature_cols = [col.removeprefix("true_delta_") for col in delta_table.columns if col.startswith("true_delta_")]
    max_features = min(max_features, 8)
    if len(feature_cols) > max_features:
        errors = []
        for feature in feature_cols:
            errors.append(
                (
                    feature,
                    float(np.nanmean(np.abs(delta_table[f"true_delta_{feature}"] - delta_table[f"pred_delta_{feature}"]))),
                    float(np.nanmean(np.abs(delta_table[f"true_delta_{feature}"]))),
                )
            )
        errors.sort(key=lambda item: item[2], reverse=True)
        feature_cols = [feature for feature, _, _ in errors[:max_features]]
    true_df = delta_table.set_index("ko_target")[[f"true_delta_{f}" for f in feature_cols]]
    pred_df = delta_table.set_index("ko_target")[[f"pred_delta_{f}" for f in feature_cols]]
    true_df.columns = [wrap_label(short_feature(f), 18) for f in feature_cols]
    pred_df.columns = [wrap_label(short_feature(f), 18) for f in feature_cols]
    true_df = true_df.T
    pred_df = pred_df.T
    err_df = pred_df - true_df
    vmax = np.nanmax(np.abs(pd.concat([true_df, pred_df]).to_numpy()))
    err_vmax = np.nanmax(np.abs(err_df.to_numpy()))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, max(5.8, 0.55 * len(true_df) + 2.8)), sharey=True, constrained_layout=True)
    annot = true_df.size <= 24
    annot_kws = {"fontsize": 8}
    sns.heatmap(true_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=annot, annot_kws=annot_kws, fmt=".2f", cbar=False, ax=axes[0])
    sns.heatmap(pred_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=annot, annot_kws=annot_kws, fmt=".2f", cbar=False, ax=axes[1])
    sns.heatmap(err_df, cmap="vlag", center=0, vmin=-err_vmax, vmax=err_vmax, annot=annot, annot_kws=annot_kws, fmt=".2f", cbar_kws={"label": "virtual - true"}, ax=axes[2])
    axes[0].set_title("Real KO change")
    axes[1].set_title("Virtual KO change")
    axes[2].set_title("Prediction error")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=0)
        ax.tick_params(axis="y", rotation=0)
    fig.suptitle("Real vs Virtual KO Changes", fontsize=16)
    fig.savefig(out_dir / "02_true_vs_virtual_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_umap(cells: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    if cells.empty:
        return
    features = [
        col
        for col in cells.columns
        if pd.api.types.is_numeric_dtype(cells[col]) and col not in {"UMAP1", "UMAP2"}
    ]
    rows = []
    for ko, sub in cells.groupby("ko_target", observed=True):
        matrix = StandardScaler().fit_transform(sub[features].to_numpy(dtype=float))
        reducer = umap.UMAP(n_neighbors=min(30, max(5, len(sub) // 30)), min_dist=0.25, random_state=11)
        coords = reducer.fit_transform(matrix)
        tmp = sub.copy()
        tmp["UMAP1"] = coords[:, 0]
        tmp["UMAP2"] = coords[:, 1]
        rows.append(tmp)
    embedded = pd.concat(rows, ignore_index=True)
    kos = list(embedded["ko_target"].drop_duplicates())[:6]
    fig, axes = plt.subplots(len(kos), 2, figsize=(12, 4.1 * len(kos) + 2.0), squeeze=False, constrained_layout=False)
    for row, ko in enumerate(kos):
        group = embedded.loc[embedded["ko_target"] == ko]
        before = group.loc[group["state"] == "control cells"]
        sns.scatterplot(data=before, x="UMAP1", y="UMAP2", color=PALETTE["control cells"], s=22, alpha=0.78, linewidth=0, ax=axes[row, 0])
        axes[row, 0].set_title(f"{ko}: before KO", pad=6)
        for state, alpha, size in [("control cells", 0.18, 18), ("virtual KO cells", 0.75, 28), ("true KO cells", 0.75, 28)]:
            part = group.loc[group["state"] == state]
            sns.scatterplot(data=part, x="UMAP1", y="UMAP2", color=PALETTE[state], s=size, alpha=alpha, linewidth=0, label=state if row == 0 else None, ax=axes[row, 1])
        centers = group.groupby("state", observed=True)[["UMAP1", "UMAP2"]].mean()
        if {"control cells", "virtual KO cells", "true KO cells"}.issubset(set(centers.index)):
            start = centers.loc["control cells"]
            for state, color, linestyle in [("virtual KO cells", PALETTE["virtual KO cells"], "-"), ("true KO cells", PALETTE["true KO cells"], "--")]:
                end = centers.loc[state]
                axes[row, 1].annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "color": color, "lw": 2.5, "linestyle": linestyle})
        axes[row, 1].set_title(f"{ko}: virtual vs real KO", pad=6)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.035), frameon=True)
        axes[0, 1].legend_.remove()
    fig.suptitle("Single-cell State Movement", fontsize=16, y=0.96)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.18, wspace=0.18)
    fig.savefig(out_dir / "03_cell_state_umap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    embedded.to_csv(out_dir / "umap_cells.csv", index=False)


def plot_auc(auc_points: pd.DataFrame, out_dir: Path, threshold: float | None = None) -> pd.DataFrame:
    setup_plot()
    if auc_points.empty:
        return pd.DataFrame()
    if threshold is None:
        threshold = float(np.nanquantile(auc_points["true_abs_delta"], 0.65))
    labels = auc_points["true_abs_delta"].to_numpy(dtype=float) >= threshold
    scores = auc_points["pred_abs_delta"].to_numpy(dtype=float)
    rows = [{"task": "strong pathway/program response", "threshold": threshold, "n_positive": int(labels.sum()), "n_negative": int((~labels).sum()), "roc_auc": np.nan}]
    fig, ax = plt.subplots(figsize=(6.8, 5.8), constrained_layout=True)
    if labels.sum() > 0 and (~labels).sum() > 0:
        fpr, tpr, _ = roc_curve(labels, scores)
        roc_auc = auc(fpr, tpr)
        rows[0]["roc_auc"] = float(roc_auc)
        ax.plot(fpr, tpr, linewidth=3, color="#4E79A7", label=f"AUC={roc_auc:.2f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="0.45", linewidth=1, label="random")
    ax.set(
        title="ROC Curve for Strong KO Responses",
        xlabel="False positive rate",
        ylabel="True positive rate",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
    )
    ax.legend(loc="lower right", frameon=True)
    fig.savefig(out_dir / "04_auc_strong_response_roc.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    result = pd.DataFrame(rows)
    result.to_csv(out_dir / "auc_summary.csv", index=False)
    return result


PEAK_FEATURE_RE = re.compile(r"peak_(chr[^_]+)_(\d+)_(\d+)_([^_]+)_(.+)")


def parse_peak_feature(feature: str) -> dict | None:
    match = PEAK_FEATURE_RE.match(feature)
    if not match:
        return None
    chrom, start, end, gene, peak_type = match.groups()
    return {
        "feature": feature,
        "chrom": chrom,
        "start": int(start),
        "end": int(end),
        "gene": gene,
        "peak_type": peak_type,
    }


def plot_peak_level_changes(delta_table: pd.DataFrame, out_dir: Path, max_peaks: int = 10) -> None:
    setup_plot()
    if delta_table.empty:
        return
    rows = []
    for col in delta_table.columns:
        if not col.startswith("true_delta_peak_"):
            continue
        feature = col.removeprefix("true_delta_")
        parsed = parse_peak_feature(feature)
        if parsed is None:
            continue
        true_values = delta_table[f"true_delta_{feature}"].to_numpy(dtype=float)
        pred_values = delta_table[f"pred_delta_{feature}"].to_numpy(dtype=float)
        parsed["true_delta"] = float(np.nanmean(true_values))
        parsed["virtual_delta"] = float(np.nanmean(pred_values))
        parsed["error"] = parsed["virtual_delta"] - parsed["true_delta"]
        parsed["score"] = abs(parsed["true_delta"])
        rows.append(parsed)
    if not rows:
        return
    plot = pd.DataFrame(rows).sort_values("score", ascending=False).head(max_peaks).copy()
    plot = plot.sort_values(["chrom", "start"])
    plot["label"] = plot["gene"] + "\n" + plot["peak_type"] + "\n" + plot["start"].astype(str)
    long = plot.melt(
        id_vars=["label", "chrom", "start", "end", "gene", "peak_type"],
        value_vars=["true_delta", "virtual_delta"],
        var_name="change",
        value_name="delta",
    )
    long["change"] = long["change"].map({"true_delta": "Real KO", "virtual_delta": "Virtual KO"})

    fig, axes = plt.subplots(1, 3, figsize=(17.5, max(5.0, 0.42 * len(plot) + 2.2)), constrained_layout=True)
    sns.barplot(data=long, x="label", y="delta", hue="change", palette=["#2A9D8F", "#E76F51"], ax=axes[0])
    axes[0].axhline(0, color="#444", ls="--", lw=1)
    axes[0].set_title("ATAC peak changes: real vs virtual KO")
    axes[0].set_xlabel("Peak")
    axes[0].set_ylabel("Peak accessibility state delta")
    axes[0].tick_params(axis="x", rotation=0)
    axes[0].legend(title="")

    heat = plot.set_index("label")[["true_delta", "virtual_delta", "error"]]
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    sns.heatmap(
        heat,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        annot_kws={"fontsize": 8},
        cbar_kws={"label": "delta"},
        ax=axes[1],
    )
    axes[1].set_title("Peak-level error heatmap")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")

    scatter = plot.copy()
    sns.scatterplot(
        data=scatter,
        x="true_delta",
        y="virtual_delta",
        hue="peak_type",
        size="score",
        sizes=(50, 180),
        palette="tab10",
        ax=axes[2],
    )
    lim = np.nanmax(np.abs(scatter[["true_delta", "virtual_delta"]].to_numpy()))
    if not np.isfinite(lim) or lim < 1e-6:
        lim = 1.0
    axes[2].plot([-lim, lim], [-lim, lim], color="0.35", linestyle="--", linewidth=1)
    axes[2].axhline(0, color="0.75", linewidth=0.8)
    axes[2].axvline(0, color="0.75", linewidth=0.8)
    axes[2].set_xlim(-lim * 1.1, lim * 1.1)
    axes[2].set_ylim(-lim * 1.1, lim * 1.1)
    axes[2].set_title("Peak direction agreement")
    axes[2].set_xlabel("Real KO peak delta")
    axes[2].set_ylabel("Virtual KO peak delta")
    axes[2].legend(title="Peak type", fontsize=8, title_fontsize=9, loc="best")
    fig.savefig(out_dir / "05_atac_peak_level_changes.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def write_report(summary: pd.DataFrame, auc_summary: pd.DataFrame, out_dir: Path) -> None:
    row = summary.iloc[0]
    auc_text = "not available"
    if not auc_summary.empty and pd.notna(auc_summary["roc_auc"].iloc[0]):
        auc_text = f"{auc_summary['roc_auc'].iloc[0]:.3f}"
    text = f"""# Virtual KO software report

## Input

- Dataset: {row['dataset']}
- Input modality: {row['input_modality']}
- State representation: {row['state_representation']}
- KO targets evaluated: {int(row['n_ko'])}
- State features: {int(row['n_features'])}

## Main results

- Mean distribution improvement: {row['mean_distribution_improvement']:.3f}
- Fraction of features improved: {row['improved_fraction']:.3f}
- Mean direction cosine: {row.get('mean_direction_cosine', np.nan):.3f}
- Mean magnitude error: {row.get('mean_abs_delta_error', np.nan):.3f}
- Strong-response ROC-AUC: {auc_text}
- Calibration method: {row.get('calibration_method', 'none')}

## How to read the figures

- `01_summary_dashboard.png`: one-screen overview of whether virtual KO is closer to real KO.
- `02_true_vs_virtual_heatmap.png`: left is real KO change, middle is virtual KO change, right is error.
- `03_cell_state_umap.png`: shows whether virtual KO cells move from control toward real KO cells.
- `04_auc_strong_response_roc.png`: tests whether the model can identify strong KO-response features.
- `05_atac_peak_level_changes.png`: shown when peak-level ATAC features are present.

## Output tables

- `virtual_cells.csv`: generated virtual KO single-cell states with matched control and true KO samples.
- `metrics.csv`: feature-level distribution distances.
- `delta_table.csv`: KO-level true and predicted deltas.
- `summary.csv`: compact result summary.
- `auc_summary.csv`: ROC-AUC summary.
"""
    (out_dir / "report.md").write_text(text, encoding="utf-8")


def make_all_plots(summary: pd.DataFrame, delta_table: pd.DataFrame, cells: pd.DataFrame, auc_points: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_summary(summary, out_dir)
    plot_true_virtual_heatmap(delta_table, out_dir)
    plot_umap(cells, out_dir)
    auc_summary = plot_auc(auc_points, out_dir)
    plot_peak_level_changes(delta_table, out_dir)
    write_report(summary, auc_summary, out_dir)
    return auc_summary
