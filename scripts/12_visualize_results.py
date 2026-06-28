from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


FIG_DIR = Path("results/figures")


def setup() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300


def savefig(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIG_DIR / name, bbox_inches="tight")
    plt.close()


def plot_papalexi_auc() -> None:
    auc = pd.read_csv("results/papalexi_multimodal_auc.csv")
    keep = auc[
        (
            (auc["target"] == "delta_pathway_IFNG_JAK_STAT")
            & (auc["direction"].isin(["absolute", "decrease"]))
        )
        | ((auc["target"] == "delta_protein_PDL1") & (auc["direction"].isin(["absolute", "decrease"])))
        | ((auc["target"] == "delta_protein_CD86") & (auc["direction"] == "absolute"))
        | ((auc["target"] == "delta_protein_CD366") & (auc["direction"] == "absolute"))
    ].copy()
    keep = keep.dropna(subset=["roc_auc"])
    label_map = {
        "delta_pathway_IFNG_JAK_STAT": "IFNG-JAK-STAT",
        "delta_protein_PDL1": "PDL1 protein",
        "delta_protein_CD86": "CD86 protein",
        "delta_protein_CD366": "CD366 protein",
    }
    keep["label"] = keep["target"].map(label_map) + " / " + keep["direction"]
    keep["model"] = keep["model"].map({"pls_pred": "PLS", "ridge_pred": "Ridge"})
    plt.figure(figsize=(9, 5.5))
    ax = sns.barplot(data=keep, x="roc_auc", y="label", hue="model", palette="Set2")
    ax.axvline(0.5, color="0.35", linestyle="--", linewidth=1)
    ax.set(xlim=(0, 1.05), xlabel="ROC-AUC", ylabel="", title="Papalexi: strong-response detection")
    ax.legend(title="")
    savefig("papalexi_auc_strong_response.png")


def plot_papalexi_corr() -> None:
    corr = pd.read_csv("results/papalexi_pathway_protein_correlation.csv")
    mat = corr.pivot(index="pathway", columns="protein", values="pearson_r")
    mat.index = mat.index.str.replace("delta_pathway_", "", regex=False)
    mat.columns = mat.columns.str.replace("delta_protein_", "", regex=False)
    plt.figure(figsize=(8, 6))
    ax = sns.heatmap(mat, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=0.5)
    ax.set(title="Papalexi: pathway-protein coupling", xlabel="Protein delta", ylabel="Pathway delta")
    savefig("papalexi_pathway_protein_correlation.png")


def plot_norman_r2() -> None:
    metrics = pd.read_csv("results/norman_system_prior_metrics.csv")
    metrics = metrics[metrics["subset"].isin(["all_combos", "has_unseen_gene", "all_genes_seen"])].copy()
    metrics["target"] = metrics["target"].str.replace("delta_program_", "", regex=False)
    metrics["subset"] = metrics["subset"].map(
        {
            "all_combos": "All combos",
            "has_unseen_gene": "Has unseen gene",
            "all_genes_seen": "All genes seen",
        }
    )
    metrics["model"] = metrics["model"].map(
        {
            "single_gene_additive": "Additive",
            "system_prior_ridge": "System prior",
        }
    )
    g = sns.catplot(
        data=metrics,
        x="r2",
        y="target",
        hue="model",
        col="subset",
        kind="bar",
        palette="Set2",
        height=5,
        aspect=1.1,
        sharex=True,
        sharey=True,
    )
    for ax in g.axes.flat:
        ax.axvline(0, color="0.25", linewidth=1)
        ax.set_ylabel("")
    g.set_axis_labels("R2 on double-gene combos", "")
    g.set_titles("{col_name}")
    g.fig.suptitle("Norman: single-gene training -> double-gene prediction", y=1.05)
    g.savefig(FIG_DIR / "norman_system_prior_r2.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_norman_scatter() -> None:
    pred = pd.read_csv("results/norman_system_prior_predictions.csv")
    targets = [
        "delta_program_ERYTHROID",
        "delta_program_GRANULOCYTE_APOPTOSIS",
        "delta_program_MAPK_TGFB",
        "delta_program_PRO_GROWTH",
    ]
    rows = []
    for target in targets:
        for _, row in pred.iterrows():
            rows.append(
                {
                    "target": target.replace("delta_program_", ""),
                    "true": row[f"true_{target}"],
                    "system_prior": row[f"system_pred_{target}"],
                    "additive": row[f"additive_pred_{target}"],
                    "all_genes_seen": row["all_genes_seen_in_single"],
                }
            )
    long = pd.DataFrame(rows)
    g = sns.FacetGrid(long, col="target", col_wrap=2, height=4.2, sharex=False, sharey=False)
    g.map_dataframe(
        sns.scatterplot,
        x="true",
        y="system_prior",
        hue="all_genes_seen",
        palette={True: "#2A9D8F", False: "#E76F51"},
        alpha=0.85,
    )
    for ax in g.axes.flat:
        lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
        hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([lo, hi], [lo, hi], color="0.25", linestyle="--", linewidth=1)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
    g.add_legend(title="All genes seen")
    g.set_axis_labels("True delta", "System-prior predicted delta")
    g.set_titles("{col_name}")
    g.fig.suptitle("Norman: system-prior predictions", y=1.03)
    g.savefig(FIG_DIR / "norman_system_prior_scatter.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_prior_hits() -> None:
    hits = pd.read_csv("results/norman_system_prior_term_hits.csv")
    top = hits.sort_values("n_prior_terms_hit", ascending=False).head(15)
    plt.figure(figsize=(9, 7))
    ax = sns.barplot(data=top, x="n_prior_terms_hit", y="ko_genes", color="#457B9D")
    ax.set(xlabel="Number of matched prior terms", ylabel="", title="Norman: prior coverage for combo perturbations")
    savefig("norman_prior_term_hits.png")


def main() -> None:
    setup()
    plot_papalexi_auc()
    plot_papalexi_corr()
    plot_norman_r2()
    plot_norman_scatter()
    plot_prior_hits()
    print(f"Saved figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
