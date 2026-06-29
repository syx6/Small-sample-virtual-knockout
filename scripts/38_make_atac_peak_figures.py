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
    sparse_rows = []
    for state in order:
        values = cells.loc[cells["state"].eq(state), chosen].to_numpy(dtype=float)
        open_like = values > 0
        sparse_rows.append(
            {
                "state": state,
                "open_fraction": float(open_like.mean()),
                "open_count": int(open_like.sum()),
                "n_cells": int(len(values)),
                "median": float(np.nanmedian(values)),
                "p99": float(np.nanquantile(values, 0.99)),
            }
        )
    sparse = pd.DataFrame(sparse_rows)
    palette = {"control cells": "#BDBDBD", "virtual KO cells": "#E76F51", "true KO cells": "#2A9D8F"}
    sns.barplot(data=sparse, x="state", y="open_fraction", hue="state", order=order, hue_order=order, palette=palette, legend=False, ax=ax)
    ax.set_ylim(0, max(0.06, float(sparse["open_fraction"].max()) * 1.45))
    for i, row in sparse.reset_index(drop=True).iterrows():
        ax.text(
            i,
            row["open_fraction"] + ax.get_ylim()[1] * 0.035,
            f"{row['open_fraction']*100:.1f}%\n({int(row['open_count'])}/{int(row['n_cells'])})",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    ax.set_title("Sparse peak open-like cell fraction", fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Open-like cells\n(score > 0)")
    ax.tick_params(axis="x", rotation=18)
    note = (
        short_name.replace("\n", " ") + "\n"
        "Sparse peak: true KO is flat when all sampled cells are closed."
    )
    ax.text(0.02, 0.97, note, transform=ax.transAxes, va="top", fontsize=8.0, color="#555")

    fig.savefig(OUT / "19_atac_peak_level_visualization.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT / '19_atac_peak_level_visualization.png'}")


if __name__ == "__main__":
    main()
