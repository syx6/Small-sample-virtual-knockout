from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"
PEAK_RUN = ROOT / "results" / "scperturb_atac_regulatory_peak_prior_quantile_kdm6a"


RUNS = [
    ("Gene activity", ROOT / "results" / "scperturb_atac_gene_activity_kdm6a"),
    ("+ all chromVAR", ROOT / "results" / "scperturb_atac_gene_activity_chromvar_kdm6a"),
    ("+ top100 chromVAR", ROOT / "results" / "scperturb_atac_gene_activity_chromvar_top100_kdm6a"),
    ("+ regulatory peaks + quantile", PEAK_RUN),
]


PEAK_RE = re.compile(r"peak_(chr[^_]+)_(\d+)_(\d+)_([^_]+)_(.+)")


def parse_peak_feature(feature: str) -> dict | None:
    match = PEAK_RE.match(feature)
    if not match:
        return None
    chrom, start, end, gene, peak_type = match.groups()
    return {
        "feature": feature,
        "chrom": chrom,
        "start": int(start),
        "end": int(end),
        "mid": (int(start) + int(end)) / 2,
        "gene": gene,
        "peak_type": peak_type,
    }


def load_metric_rows() -> pd.DataFrame:
    rows = []
    for label, path in RUNS:
        summary = pd.read_csv(path / "summary.csv").iloc[0]
        auc = pd.read_csv(path / "auc_summary.csv").iloc[0]
        rows.extend(
            [
                {"model": label, "metric": "AUC", "value": float(auc["roc_auc"])},
                {"model": label, "metric": "Direction", "value": float(summary["mean_direction_cosine"])},
                {"model": label, "metric": "MAE\n(lower better)", "value": float(summary["mean_abs_delta_error"])},
                {"model": label, "metric": "Feature\nhit-rate", "value": float(summary["improved_fraction"])},
            ]
        )
    return pd.DataFrame(rows)


def peak_delta_table() -> pd.DataFrame:
    delta = pd.read_csv(PEAK_RUN / "delta_table.csv").iloc[0]
    rows = []
    for col in delta.index:
        if not col.startswith("true_delta_peak_"):
            continue
        feature = col.removeprefix("true_delta_")
        parsed = parse_peak_feature(feature)
        if parsed is None:
            continue
        parsed["true_delta"] = float(delta[f"true_delta_{feature}"])
        parsed["virtual_delta"] = float(delta[f"pred_delta_{feature}"])
        parsed["error"] = parsed["virtual_delta"] - parsed["true_delta"]
        rows.append(parsed)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")
    peak_delta = peak_delta_table()
    kdm6a = peak_delta.loc[peak_delta["gene"].str.upper().eq("KDM6A")].sort_values("start")
    if kdm6a.empty:
        kdm6a = peak_delta.sort_values("true_delta", key=lambda x: np.abs(x), ascending=False).head(8)

    fig = plt.figure(figsize=(15.5, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.05], hspace=0.34, wspace=0.28)
    fig.suptitle("ATAC peak-level virtual KO view: KDM6A example", fontsize=17, fontweight="bold")

    ax = fig.add_subplot(gs[0, :2])
    metric_rows = load_metric_rows()
    sns.barplot(data=metric_rows, x="metric", y="value", hue="model", palette=["#6b9ac4", "#984ea3", "#2a9d8f", "#d95f02"], ax=ax)
    ax.set_title("ATAC model variants after adding real peak features", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Metric value")
    ax.axhline(0, color="#444", ls="--", lw=1)
    ax.legend(title="", ncol=2, fontsize=8.8, loc="upper right")
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=7.5, padding=2)

    ax = fig.add_subplot(gs[0, 2])
    ax.set_axis_off()
    text = (
        "What changed?\n\n"
        "- Peak selection now combines target locus, marker peaks, KO effect, accessibility and TF/motif prior.\n"
        "- The regulatory peak model reaches AUC 0.67 and direction 0.77.\n"
        "- Quantile/zero-inflated calibration improves feature-level distribution hit-rate.\n"
        "- AUC stays stable because shape calibration changes single-cell distributions, not the KO mean direction.\n\n"
        "This is a real peak-level regulatory-prior view, not only gene activity or motif proxy."
    )
    ax.text(0.02, 0.96, text, va="top", fontsize=10.5, linespacing=1.3)

    ax = fig.add_subplot(gs[1, 0])
    plot = kdm6a.copy()
    plot["label"] = plot["peak_type"] + "\n" + plot["start"].astype(str)
    long = plot.melt(id_vars=["label", "peak_type", "start"], value_vars=["true_delta", "virtual_delta"], var_name="change", value_name="delta")
    long["change"] = long["change"].map({"true_delta": "Real KO", "virtual_delta": "Virtual KO"})
    sns.barplot(data=long, x="label", y="delta", hue="change", palette=["#2A9D8F", "#E76F51"], ax=ax)
    ax.axhline(0, color="#444", ls="--", lw=1)
    ax.set_title("KDM6A locus peaks: real vs virtual delta", fontsize=12, fontweight="bold")
    ax.set_xlabel("Peak type and genomic start")
    ax.set_ylabel("Peak accessibility state delta")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(title="")

    ax = fig.add_subplot(gs[1, 1])
    heat = plot.set_index("label")[["true_delta", "virtual_delta", "error"]]
    vmax = np.nanmax(np.abs(heat.to_numpy()))
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", annot_kws={"fontsize": 8}, cbar_kws={"label": "delta"}, ax=ax)
    ax.set_title("Peak delta heatmap", fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")

    ax = fig.add_subplot(gs[1, 2])
    cells = pd.read_csv(PEAK_RUN / "virtual_cells.csv")
    chosen_row = kdm6a.loc[kdm6a["peak_type"].str.contains("Promoter", case=False, na=False)].head(1)
    if chosen_row.empty:
        chosen_row = kdm6a.head(1)
    chosen_info = chosen_row.iloc[0]
    chosen = chosen_info["feature"]
    short_name = chosen_info["peak_type"] + f" peak\n{chosen_info['chrom']}:{int(chosen_info['start'])}-{int(chosen_info['end'])}"
    order = ["control cells", "virtual KO cells", "true KO cells"]
    sns.violinplot(
        data=cells,
        x="state",
        y=chosen,
        hue="state",
        order=order,
        hue_order=order,
        palette={"control cells": "#BDBDBD", "virtual KO cells": "#E76F51", "true KO cells": "#2A9D8F"},
        inner="quartile",
        cut=0,
        legend=False,
        ax=ax,
    )
    low, high = np.nanquantile(cells[chosen].to_numpy(dtype=float), [0.01, 0.99])
    if np.isfinite(low) and np.isfinite(high) and high > low:
        pad = (high - low) * 0.08
        ax.set_ylim(low - pad, high + pad)
    ax.set_title("Single-cell peak accessibility distribution", fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel(short_name)
    ax.tick_params(axis="x", rotation=18)

    fig.savefig(OUT / "19_atac_peak_level_visualization.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT / '19_atac_peak_level_visualization.png'}")


if __name__ == "__main__":
    main()
