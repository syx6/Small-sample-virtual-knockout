from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"


RUNS = [
    ("Gene activity only", ROOT / "results" / "scperturb_atac_gene_activity_kdm6a"),
    ("Gene activity + all chromVAR", ROOT / "results" / "scperturb_atac_gene_activity_chromvar_kdm6a"),
    ("Gene activity + top100 chromVAR", ROOT / "results" / "scperturb_atac_gene_activity_chromvar_top100_kdm6a"),
]


def load_rows() -> pd.DataFrame:
    rows = []
    for label, path in RUNS:
        summary = pd.read_csv(path / "summary.csv").iloc[0]
        auc = pd.read_csv(path / "auc_summary.csv").iloc[0]
        rows.extend(
            [
                {"model": label, "metric": "AUC", "value": float(auc["roc_auc"])},
                {"model": label, "metric": "Direction", "value": float(summary["mean_direction_cosine"])},
                {"model": label, "metric": "Feature\nhit-rate", "value": float(summary["improved_fraction"])},
                {"model": label, "metric": "State\nimprovement", "value": float(summary["mean_distribution_improvement"])},
            ]
        )
    return pd.DataFrame(rows)


def image_panel(ax, path: Path, title: str) -> None:
    ax.imshow(plt.imread(path))
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_axis_off()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")
    data = load_rows()
    fig = plt.figure(figsize=(15, 8.3))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0], hspace=0.34, wspace=0.22)
    ax = fig.add_subplot(gs[0, :2])
    palette = ["#6b9ac4", "#984ea3", "#2a9d8f"]
    sns.barplot(data=data, x="metric", y="value", hue="model", palette=palette, ax=ax)
    ax.axhline(0, color="#444", ls="--", lw=1)
    ax.set_ylim(min(-0.05, data["value"].min() - 0.05), 0.75)
    ax.set_title("ATAC prior and chromVAR ablation for KDM6A virtual KO", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Metric value")
    ax.get_legend().remove()
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=8, padding=2)

    ax = fig.add_subplot(gs[0, 2])
    ax.set_axis_off()
    labels = ["Gene activity only", "Gene activity + all chromVAR", "Gene activity + top100 chromVAR"]
    for i, (label, color) in enumerate(zip(labels, palette)):
        y = 0.98 - i * 0.08
        ax.add_patch(plt.Rectangle((0.02, y - 0.025), 0.045, 0.025, color=color, transform=ax.transAxes, clip_on=False))
        ax.text(0.08, y - 0.012, label, transform=ax.transAxes, va="center", fontsize=9.2)
    text = (
        "How to read this\n\n"
        "1. Gene activity alone keeps the best AUC.\n"
        "2. All chromVAR motifs add noise and lower AUC.\n"
        "3. Selecting top-variable motifs improves direction and MAE,\n"
        "   but strong-response ranking is still weaker.\n\n"
        "Conclusion: motif/TF information should be filtered or weighted;\n"
        "it should not be added blindly."
    )
    ax.text(0.02, 0.68, text, va="top", fontsize=10.5, linespacing=1.28)

    image_panel(
        fig.add_subplot(gs[1, 0]),
        ROOT / "results" / "scperturb_atac_gene_activity_kdm6a" / "04_auc_strong_response_roc.png",
        "Gene activity only: ROC",
    )
    image_panel(
        fig.add_subplot(gs[1, 1]),
        ROOT / "results" / "scperturb_atac_gene_activity_chromvar_top100_kdm6a" / "02_true_vs_virtual_heatmap.png",
        "Top100 chromVAR: heatmap",
    )
    image_panel(
        fig.add_subplot(gs[1, 2]),
        ROOT / "results" / "scperturb_atac_gene_activity_chromvar_top100_kdm6a" / "04_auc_strong_response_roc.png",
        "Top100 chromVAR: ROC",
    )
    fig.suptitle("ATAC-specific optimization: TF/motif priors and chromVAR input", fontsize=17, fontweight="bold")
    fig.savefig(OUT / "18_atac_chromvar_prior_ablation.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT / '18_atac_chromvar_prior_ablation.png'}")


if __name__ == "__main__":
    main()
