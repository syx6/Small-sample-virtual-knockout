from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"


def card(ax, title: str, value: str, note: str, color: str) -> None:
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.10, edgecolor=color, lw=1.6))
    ax.text(0.06, 0.76, title, fontsize=11, fontweight="bold", color="#222")
    ax.text(0.06, 0.42, value, fontsize=20, fontweight="bold", color=color)
    ax.text(0.06, 0.15, note, fontsize=9, color="#555", va="bottom")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")
    summary = pd.read_csv(ROOT / "results" / "hmpcite_multimodal_doubleko_cebp_med12" / "summary.csv").iloc[0]
    auc = pd.read_csv(ROOT / "results" / "hmpcite_multimodal_doubleko_cebp_med12" / "auc_summary.csv").iloc[0]
    counts = pd.read_csv(ROOT / "data" / "hmpcite_gse243244" / "hmpcite_ko_counts_threshold10.csv")
    metrics = pd.read_csv(ROOT / "results" / "hmpcite_multimodal_doubleko_interaction" / "double_interaction_metrics.csv")
    agg = metrics.loc[metrics["subset"].eq("all_combos")].groupby("model")[["mae", "r2", "roc_auc_abs_gt_0.15"]].mean()

    fig = plt.figure(figsize=(14, 7.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.9, 1.1, 1.3], hspace=0.42, wspace=0.28)
    fig.suptitle("HMPCITE-seq multimodal double-KO dataset integration", fontsize=18, fontweight="bold", y=0.98)

    n_cells = int(counts["n_cells"].sum())
    n_double = int((counts["n_ko_genes"] == 2).sum())
    card(fig.add_subplot(gs[0, 0]), "Dataset", f"{n_cells:,}", "cells with RNA + ADT + GDO labels", "#377eb8")
    card(fig.add_subplot(gs[0, 1]), "Double KOs", f"{n_double}", "real double-gene combinations", "#4daf4a")
    card(fig.add_subplot(gs[0, 2]), "Cebpb+Med12 AUC", f"{auc['roc_auc']:.2f}", "strong-response ROC curve", "#984ea3")
    card(fig.add_subplot(gs[0, 3]), "Direction match", f"{summary['mean_direction_cosine']:.2f}", "virtual vs real KO direction", "#e41a1c")

    ax = fig.add_subplot(gs[1, :2])
    dist = counts.groupby("n_ko_genes")["n_cells"].sum().rename(index={0: "control", 1: "single KO", 2: "double KO"})
    ax.bar(dist.index.astype(str), dist.values, color=["#bdbdbd", "#377eb8", "#4daf4a"])
    ax.set_title("Usable perturbation labels", fontsize=13, fontweight="bold")
    ax.set_ylabel("Cells")
    for i, v in enumerate(dist.values):
        ax.text(i, v + 250, f"{int(v):,}", ha="center", fontsize=10)

    ax = fig.add_subplot(gs[1, 2:])
    top_double = counts.loc[counts["n_ko_genes"].eq(2)].head(10)
    sns.barplot(data=top_double, x="n_cells", y="ko_target", color="#4daf4a", ax=ax)
    ax.set_title("Top double-KO combinations", fontsize=13, fontweight="bold")
    ax.set_xlabel("Cells")
    ax.set_ylabel("")

    ax = fig.add_subplot(gs[2, :])
    plot = agg.rename(index={"single_gene_additive": "Additive baseline", "interaction_residual": "Interaction residual"}).reset_index(names="model")
    long = plot.melt(id_vars="model", value_vars=["mae", "r2", "roc_auc_abs_gt_0.15"], var_name="metric", value_name="value")
    long["metric"] = long["metric"].replace({"mae": "MAE\nlower better", "r2": "R2\nhigher better", "roc_auc_abs_gt_0.15": "ROC-AUC\nhigher better"})
    sns.barplot(
        data=long,
        x="metric",
        y="value",
        hue="model",
        hue_order=["Additive baseline", "Interaction residual"],
        palette={"Additive baseline": "#9e9e9e", "Interaction residual": "#377eb8"},
        ax=ax,
    )
    ax.axhline(0, color="#444", lw=1, ls="--")
    ax.set_title("Across 55 real multimodal double KOs: interaction model improves magnitude prediction", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean metric")
    ax.legend(title="", loc="upper left")
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=9, padding=2)

    fig.text(
        0.5,
        0.02,
        "Source: GSE243244 HMPCITE-seq perturbation sample. RNA is scored as pathways; ADT proteins are included as extra state features.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.savefig(OUT / "14_hmpcite_multimodal_doubleko_summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
