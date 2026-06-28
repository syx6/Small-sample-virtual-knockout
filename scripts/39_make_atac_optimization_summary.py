from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"


RUNS = [
    ("Gene activity\nbaseline", ROOT / "results" / "scperturb_atac_gene_activity_kdm6a"),
    ("Top100 chromVAR\nvariance", ROOT / "results" / "scperturb_atac_gene_activity_chromvar_top100_kdm6a"),
    ("Weighted + hybrid\nchromVAR", ROOT / "results" / "scperturb_atac_weighted_hybrid_chromvar_kdm6a"),
    ("Weighted + hybrid\n+ locus peaks", ROOT / "results" / "scperturb_atac_weighted_hybrid_peak_kdm6a"),
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
                {"model": label, "metric": "MAE\n(lower better)", "value": float(summary["mean_abs_delta_error"])},
                {"model": label, "metric": "Distribution\nimprovement", "value": float(summary["mean_distribution_improvement"])},
                {"model": label, "metric": "Feature\nhit-rate", "value": float(summary["improved_fraction"])},
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
    fig = plt.figure(figsize=(16, 9.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.05], hspace=0.34, wspace=0.26)
    fig.suptitle("ATAC optimization summary: weighted priors, hybrid selection, and locus-aware peaks", fontsize=17, fontweight="bold")

    ax = fig.add_subplot(gs[0, :2])
    palette = ["#6b9ac4", "#2a9d8f", "#4c78a8", "#d95f02"]
    sns.barplot(data=data, x="metric", y="value", hue="model", palette=palette, ax=ax)
    ax.axhline(0, color="#444", ls="--", lw=1)
    ax.set_title("KDM6A ATAC virtual KO ablation", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Metric value")
    ax.get_legend().remove()
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", fontsize=7, padding=2)

    ax = fig.add_subplot(gs[0, 2])
    ax.set_axis_off()
    for i, ((label, _), color) in enumerate(zip(RUNS, palette)):
        y = 0.98 - i * 0.075
        ax.add_patch(plt.Rectangle((0.02, y - 0.023), 0.04, 0.022, color=color, transform=ax.transAxes, clip_on=False))
        ax.text(0.08, y - 0.012, label.replace("\n", " "), transform=ax.transAxes, va="center", fontsize=8.6)
    text = (
        "Interpretation\n\n"
        "1. Weighted TF/motif priors + hybrid feature selection improve the overall ATAC signal.\n"
        "2. Hybrid chromVAR gives the best direction, MAE, and distribution improvement.\n"
        "3. Adding locus-aware peaks gives the best AUC, so it helps rank strong-response features.\n"
        "4. Peak-level cell distributions remain difficult; use peak plots as evidence, not as a solved generator."
    )
    ax.text(0.02, 0.64, text, va="top", fontsize=10.5, linespacing=1.3)

    image_panel(
        fig.add_subplot(gs[1, 0]),
        ROOT / "results" / "scperturb_atac_weighted_hybrid_chromvar_kdm6a" / "04_auc_strong_response_roc.png",
        "Weighted + hybrid chromVAR: ROC",
    )
    image_panel(
        fig.add_subplot(gs[1, 1]),
        ROOT / "results" / "scperturb_atac_weighted_hybrid_peak_kdm6a" / "05_atac_peak_level_changes.png",
        "Weighted + peaks: peak-level changes",
    )
    image_panel(
        fig.add_subplot(gs[1, 2]),
        ROOT / "results" / "scperturb_atac_weighted_hybrid_peak_kdm6a" / "04_auc_strong_response_roc.png",
        "Weighted + peaks: ROC",
    )

    fig.savefig(OUT / "20_atac_weighted_prior_feature_selection_summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT / '20_atac_weighted_prior_feature_selection_summary.png'}")


if __name__ == "__main__":
    main()
