from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc, r2_score, roc_curve


FIG_DIR = Path("results/figures")


def setup() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"]:
        if font in available_fonts:
            plt.rcParams["font.family"] = font
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300


def signed_score(values: pd.Series, direction: str) -> np.ndarray:
    arr = values.to_numpy(dtype=float)
    if direction == "decrease":
        return -arr
    if direction == "increase":
        return arr
    if direction == "absolute":
        return np.abs(arr)
    raise ValueError(direction)


def add_roc_curve(ax, df: pd.DataFrame, true_col: str, pred_col: str, direction: str, threshold: float, label: str) -> None:
    y_score_true = signed_score(df[true_col], direction)
    y_score_pred = signed_score(df[pred_col], direction)
    y_true = y_score_true >= threshold
    if y_true.sum() == 0 or (~y_true).sum() == 0:
        return
    fpr, tpr, _ = roc_curve(y_true, y_score_pred)
    roc_auc = auc(fpr, tpr)
    ax.plot(fpr, tpr, linewidth=2.5, label=f"{label} (AUC={roc_auc:.2f})")


def plot_papalexi_roc_curves() -> None:
    df = pd.read_csv("results/papalexi_multimodal_joint_predictions.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharex=True, sharey=True)
    tasks = [
        {
            "title": "IFNG pathway decrease",
            "true": "true_delta_pathway_IFNG_JAK_STAT",
            "pls": "pls_pred_delta_pathway_IFNG_JAK_STAT",
            "ridge": "ridge_pred_delta_pathway_IFNG_JAK_STAT",
            "direction": "decrease",
            "threshold": 0.15,
        },
        {
            "title": "PDL1 protein decrease",
            "true": "true_delta_protein_PDL1",
            "pls": "pls_pred_delta_protein_PDL1",
            "ridge": "ridge_pred_delta_protein_PDL1",
            "direction": "decrease",
            "threshold": 0.15,
        },
        {
            "title": "CD86 protein strong change",
            "true": "true_delta_protein_CD86",
            "pls": "pls_pred_delta_protein_CD86",
            "ridge": "ridge_pred_delta_protein_CD86",
            "direction": "absolute",
            "threshold": 0.15,
        },
    ]
    for ax, task in zip(axes, tasks):
        add_roc_curve(ax, df, task["true"], task["pls"], task["direction"], task["threshold"], "PLS")
        add_roc_curve(ax, df, task["true"], task["ridge"], task["direction"], task["threshold"], "Ridge")
        ax.plot([0, 1], [0, 1], linestyle="--", color="0.45", linewidth=1)
        ax.set_title(task["title"], fontsize=14)
        ax.set_xlabel("False positive rate")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.legend(loc="lower right", fontsize=10)
    axes[0].set_ylabel("True positive rate")
    fig.suptitle("Papalexi: ROC curves for strong-response detection", fontsize=18, y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "papalexi_roc_curves.png", bbox_inches="tight")
    plt.close(fig)


def plot_norman_true_vs_pred() -> None:
    df = pd.read_csv("results/norman_system_prior_predictions.csv")
    tasks = [
        ("delta_program_GRANULOCYTE_APOPTOSIS", "Granulocyte/apoptosis"),
        ("delta_program_MAPK_TGFB", "MAPK/TGFB"),
        ("delta_program_PRO_GROWTH", "Pro-growth"),
        ("delta_program_ERYTHROID", "Erythroid"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, (target, title) in zip(axes.flat, tasks):
        true = df[f"true_{target}"].to_numpy(dtype=float)
        pred = df[f"system_pred_{target}"].to_numpy(dtype=float)
        seen = df["all_genes_seen_in_single"].astype(bool)
        sns.scatterplot(
            x=true,
            y=pred,
            hue=seen.map({True: "seen genes", False: "unseen gene"}),
            palette={"seen genes": "#2A9D8F", "unseen gene": "#E76F51"},
            s=70,
            alpha=0.85,
            ax=ax,
        )
        lo = min(true.min(), pred.min())
        hi = max(true.max(), pred.max())
        pad = 0.05 * (hi - lo + 1e-9)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], linestyle="--", color="0.35", linewidth=1)
        r2_all = r2_score(true, pred)
        unseen_mask = ~seen.to_numpy()
        r2_unseen = r2_score(true[unseen_mask], pred[unseen_mask]) if unseen_mask.sum() >= 3 else np.nan
        ax.text(
            0.04,
            0.94,
            f"R2 all={r2_all:.2f}\nR2 unseen={r2_unseen:.2f}",
            transform=ax.transAxes,
            va="top",
            fontsize=11,
            bbox=dict(facecolor="white", edgecolor="0.8", boxstyle="round,pad=0.25"),
        )
        ax.set_title(title)
        ax.set_xlabel("True delta")
        ax.set_ylabel("Predicted delta")
        ax.legend(title="", fontsize=9, loc="lower right")
    fig.suptitle("Norman: true vs predicted double-gene KO effects", fontsize=18, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "norman_true_vs_pred_system_prior.png", bbox_inches="tight")
    plt.close(fig)


def plot_norman_r2_delta() -> None:
    metrics = pd.read_csv("results/norman_system_prior_metrics.csv")
    subset = metrics[metrics["subset"] == "has_unseen_gene"].copy()
    wide = subset.pivot(index="target", columns="model", values="r2").reset_index()
    wide["delta_r2"] = wide["system_prior_ridge"] - wide["single_gene_additive"]
    wide["target"] = wide["target"].str.replace("delta_program_", "", regex=False)
    wide = wide.sort_values("delta_r2", ascending=True)
    plt.figure(figsize=(8, 5.5))
    ax = sns.barplot(data=wide, x="delta_r2", y="target", color="#457B9D")
    ax.axvline(0, color="0.25", linewidth=1)
    ax.set(
        xlabel="R2 improvement from system prior\n(has-unseen-gene double KO)",
        ylabel="",
        title="Does system prior improve combo extrapolation?",
    )
    for container in ax.containers:
        ax.bar_label(container, fmt="%.2f", padding=4, fontsize=10)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "norman_r2_improvement_system_prior.png", bbox_inches="tight")
    plt.close()


def main() -> None:
    setup()
    plot_papalexi_roc_curves()
    plot_norman_true_vs_pred()
    plot_norman_r2_delta()
    print("Saved standard evaluation plots.")


if __name__ == "__main__":
    main()
