from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"


def read_metric(result_dir: Path) -> tuple[pd.Series, pd.Series]:
    return pd.read_csv(result_dir / "summary.csv").iloc[0], pd.read_csv(result_dir / "auc_summary.csv").iloc[0]


def card(ax, title: str, value: str, note: str, color: str) -> None:
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.10, edgecolor=color, lw=1.6))
    ax.text(0.06, 0.76, title, fontsize=10.5, fontweight="bold", color="#222")
    ax.text(0.06, 0.42, value, fontsize=20, fontweight="bold", color=color)
    ax.text(0.06, 0.14, note, fontsize=8.8, color="#555", va="bottom")


def image_panel(ax, path: Path, title: str) -> None:
    ax.imshow(plt.imread(path))
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_axis_off()


def make_extension_summary() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    hmp_dir = ROOT / "results" / "hmpcite_multimodal_doubleko_extra_obsm_demo"
    atac_dir = ROOT / "results" / "scperturb_atac_gene_activity_kdm6a"
    hmp, hmp_auc = read_metric(hmp_dir)
    atac, atac_auc = read_metric(atac_dir)

    fig = plt.figure(figsize=(15.2, 7.6))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.9, 1.0, 1.25], hspace=0.38, wspace=0.26)
    fig.suptitle("Multimodal extension: RNA+ADT double-KO and ATAC perturbation", fontsize=17, fontweight="bold", y=0.98)

    card(fig.add_subplot(gs[0, 0]), "RNA+ADT double KO", f"{hmp_auc['roc_auc']:.2f}", "Cebpb+Med12 ROC-AUC", "#377eb8")
    card(fig.add_subplot(gs[0, 1]), "Direction match", f"{hmp['mean_direction_cosine']:.2f}", "RNA pathway + ADT state", "#4daf4a")
    card(fig.add_subplot(gs[0, 2]), "ATAC perturbation", f"{atac_auc['roc_auc']:.2f}", "KDM6A gene activity ROC-AUC", "#984ea3")
    card(fig.add_subplot(gs[0, 3]), "ATAC direction", f"{atac['mean_direction_cosine']:.2f}", "regulatory layer is harder", "#e41a1c")

    rows = pd.DataFrame(
        [
            {
                "dataset": "HMPCITE\nRNA+ADT",
                "AUC": hmp_auc["roc_auc"],
                "Direction": hmp["mean_direction_cosine"],
                "Feature\nhit-rate": hmp["improved_fraction"],
                "State\nimprovement": hmp["mean_distribution_improvement"],
            },
            {
                "dataset": "scPerturb\nATAC",
                "AUC": atac_auc["roc_auc"],
                "Direction": atac["mean_direction_cosine"],
                "Feature\nhit-rate": atac["improved_fraction"],
                "State\nimprovement": atac["mean_distribution_improvement"],
            },
        ]
    )
    metric_long = rows.melt(id_vars="dataset", var_name="metric", value_name="value")
    ax = fig.add_subplot(gs[1:, :2])
    sns.barplot(data=metric_long, x="metric", y="value", hue="dataset", palette=["#377eb8", "#984ea3"], ax=ax)
    ax.axhline(0, color="#444", lw=1, ls="--")
    ax.set_title("What improves, and what remains hard?", fontsize=13, fontweight="bold", pad=8)
    ax.set_xlabel("")
    ax.set_ylabel("Metric value")
    ax.tick_params(axis="x", rotation=0, labelsize=9)
    ax.legend(title="", loc="upper right")
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=8, padding=2)

    ax = fig.add_subplot(gs[1, 2:])
    ax.set_axis_off()
    text = (
        "Interpretation\n\n"
        "1. RNA+ADT double-KO is now a real labeled evaluation.\n"
        "2. ATAC gene activity is also supported and evaluated.\n"
        "3. ATAC is harder: weaker direction match and near-zero distribution gain.\n"
        "4. A future true RNA+ADT+ATAC perturbation dataset can use the same interface."
    )
    ax.text(0.02, 0.92, text, va="top", fontsize=11, linespacing=1.35)

    ax = fig.add_subplot(gs[2, 2:])
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0.02, 0.25), 0.28, 0.42, facecolor="#e8f1fb", edgecolor="#377eb8", lw=1.5))
    ax.add_patch(plt.Rectangle((0.36, 0.25), 0.28, 0.42, facecolor="#eaf5ea", edgecolor="#4daf4a", lw=1.5))
    ax.add_patch(plt.Rectangle((0.70, 0.25), 0.28, 0.42, facecolor="#f2e8f5", edgecolor="#984ea3", lw=1.5))
    ax.text(0.16, 0.53, "RNA", ha="center", fontsize=16, fontweight="bold", color="#377eb8")
    ax.text(0.50, 0.53, "ADT", ha="center", fontsize=16, fontweight="bold", color="#4daf4a")
    ax.text(0.84, 0.53, "ATAC", ha="center", fontsize=16, fontweight="bold", color="#984ea3")
    ax.text(0.16, 0.37, "pathway scores", ha="center", fontsize=9)
    ax.text(0.50, 0.37, "protein scores", ha="center", fontsize=9)
    ax.text(0.84, 0.37, "gene activity /\nchromVAR", ha="center", fontsize=9)
    ax.text(0.5, 0.08, "--extra-obsm protein:protein,atac:atac", ha="center", fontsize=11, fontweight="bold", color="#333")

    fig.text(
        0.5,
        0.015,
        "Use labeled perturbation data for accuracy; use unlabeled tri-modal data only for state conversion or reference-model application.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.savefig(OUT / "15_rna_adt_atac_extension_summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_result_gallery() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hmp_dir = ROOT / "results" / "hmpcite_multimodal_doubleko_extra_obsm_demo"
    atac_dir = ROOT / "results" / "scperturb_atac_gene_activity_kdm6a"
    fig, axes = plt.subplots(2, 3, figsize=(16.4, 9.0))
    image_panel(axes[0, 0], hmp_dir / "02_true_vs_virtual_heatmap.png", "RNA+ADT double KO: heatmap")
    image_panel(axes[0, 1], hmp_dir / "03_cell_state_umap.png", "RNA+ADT double KO: UMAP")
    image_panel(axes[0, 2], hmp_dir / "04_auc_strong_response_roc.png", "RNA+ADT double KO: ROC")
    image_panel(axes[1, 0], atac_dir / "02_true_vs_virtual_heatmap.png", "ATAC gene activity: heatmap")
    image_panel(axes[1, 1], atac_dir / "03_cell_state_umap.png", "ATAC gene activity: UMAP")
    image_panel(axes[1, 2], atac_dir / "04_auc_strong_response_roc.png", "ATAC gene activity: ROC")
    fig.suptitle("Modality extension result gallery: what the user can actually inspect", fontsize=17, fontweight="bold")
    fig.text(
        0.5,
        0.02,
        "Top row: real RNA+ADT double-KO evaluation. Bottom row: ATAC/gene-activity perturbation evaluation.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.88, bottom=0.08, wspace=0.08, hspace=0.22)
    fig.savefig(OUT / "16_modality_extension_result_gallery.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_input_mode_matrix() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [
        ("RNA only + KO labels", "Evaluate", "AUC / heatmap / UMAP"),
        ("RNA + ADT + KO labels", "Evaluate", "multimodal state accuracy"),
        ("ATAC gene activity + KO labels", "Evaluate", "regulatory state accuracy"),
        ("RNA + ADT + ATAC + KO labels", "Evaluate", "full tri-modal accuracy"),
        ("RNA + ADT + ATAC without KO labels", "Apply only", "predicted shift, no internal AUC"),
    ]
    fig, ax = plt.subplots(figsize=(12, 4.7))
    ax.set_axis_off()
    ax.text(0.02, 0.92, "Which visualizations are valid for each input?", fontsize=17, fontweight="bold")
    headers = ["Input data", "Mode", "Valid visual output"]
    xs = [0.02, 0.47, 0.67]
    widths = [0.40, 0.16, 0.29]
    y = 0.76
    for x, w, h in zip(xs, widths, headers):
        ax.add_patch(plt.Rectangle((x, y), w, 0.10, color="#333"))
        ax.text(x + 0.012, y + 0.052, h, va="center", fontsize=10, color="white", fontweight="bold")
    y -= 0.12
    for i, row in enumerate(rows):
        bg = "#f5f7fb" if i % 2 == 0 else "white"
        for x, w in zip(xs, widths):
            ax.add_patch(plt.Rectangle((x, y), w, 0.105, facecolor=bg, edgecolor="#dddddd", lw=0.8))
        for x, text in zip(xs, row):
            ax.text(x + 0.012, y + 0.054, text, va="center", fontsize=10, color="#222")
        y -= 0.105
    ax.text(0.02, 0.06, "Rule: AUC/R2/MAE need real perturbation labels. UMAP without labels is only a visual prediction.", fontsize=10, color="#555")
    fig.savefig(OUT / "17_multimodal_input_visualization_matrix.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    make_extension_summary()
    make_result_gallery()
    make_input_mode_matrix()
    print(f"Wrote modality extension figures to {OUT}")


if __name__ == "__main__":
    main()
