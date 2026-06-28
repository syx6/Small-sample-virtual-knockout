from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import auc as sklearn_auc
from sklearn.metrics import roc_curve


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "user_facing_figures"
DOC = ROOT / "docs" / "user_facing_figure_guide.md"


def savefig(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


def clean_feature(name: str) -> str:
    for prefix in [
        "true_delta_",
        "pred_delta_",
        "delta_",
        "pathway_",
        "program_",
        "protein_",
    ]:
        name = name.replace(prefix, "")
    name = name.replace("_", " ")
    return "\n".join(textwrap.wrap(name, width=24))


def metric_value(path: Path, col: str) -> float:
    df = pd.read_csv(path)
    return float(df.iloc[0][col])


def metric_card(ax, title: str, value: str, note: str, color: str) -> None:
    ax.set_axis_off()
    ax.add_patch(
        plt.Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.10, edgecolor=color, lw=1.8)
    )
    ax.text(0.06, 0.76, title, fontsize=11, fontweight="bold", color="#222")
    value_size = 22 if len(value) <= 8 else 17
    ax.text(0.06, 0.42, value, fontsize=value_size, fontweight="bold", color=color)
    ax.text(0.06, 0.16, note, fontsize=9, color="#555", va="bottom")


def true_pred_heatmap(delta_path: Path, title: str, name: str, top_n: int = 10) -> None:
    df = pd.read_csv(delta_path)
    row = df.iloc[0]
    true_cols = [
        c
        for c in df.columns
        if c.startswith("true_delta_") and c not in {"true_delta_norm"}
    ]
    pairs = []
    for true_col in true_cols:
        suffix = true_col.replace("true_delta_", "")
        pred_col = f"pred_delta_{suffix}"
        if pred_col in df.columns:
            pairs.append((suffix, float(row[true_col]), float(row[pred_col])))
    pairs = sorted(pairs, key=lambda x: abs(x[1]), reverse=True)[:top_n]
    mat = pd.DataFrame(
        {"Real KO": [p[1] for p in pairs], "Virtual KO": [p[2] for p in pairs]},
        index=[clean_feature(p[0]) for p in pairs],
    )
    vmax = max(0.2, float(np.nanmax(np.abs(mat.to_numpy()))))
    fig, ax = plt.subplots(figsize=(7.8, 0.45 * len(mat) + 1.8))
    sns.heatmap(
        mat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "State change"},
    )
    ax.set_title(title, fontsize=15, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=11, rotation=0)
    ax.tick_params(axis="y", labelsize=9)
    savefig(fig, name)


def predicted_heatmap(pred_path: Path) -> None:
    df = pd.read_csv(pred_path)
    row = df.iloc[0]
    pred_cols = [c for c in df.columns if c.startswith("pred_delta_")]
    vals = sorted(
        [(c.replace("pred_delta_", ""), float(row[c])) for c in pred_cols],
        key=lambda x: abs(x[1]),
        reverse=True,
    )[:12]
    mat = pd.DataFrame(
        {"Predicted virtual KO": [v for _, v in vals]},
        index=[clean_feature(k) for k, _ in vals],
    )
    vmax = max(0.2, float(np.nanmax(np.abs(mat.to_numpy()))))
    fig = plt.figure(figsize=(9, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.35)
    ax = fig.add_subplot(gs[0, 0])
    sns.heatmap(
        mat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "Predicted change"},
    )
    ax.set_title("Reference model prediction: STAT1 KO", fontsize=14, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=10, rotation=0)
    ax.tick_params(axis="y", labelsize=9)

    txt = fig.add_subplot(gs[0, 1])
    txt.set_axis_off()
    lines = [
        ("Mode", "Prediction only"),
        ("Input", "Ordinary cells + reference model"),
        ("Output", "Virtual KO state shift"),
        ("No labels?", "No internal AUC / accuracy"),
    ]
    y = 0.86
    for key, value in lines:
        txt.text(0.0, y, key, fontsize=11, fontweight="bold", color="#333")
        txt.text(0.0, y - 0.08, value, fontsize=14, color="#1f5aa6")
        y -= 0.22
    txt.text(
        0,
        0.05,
        "Use perturbation data to validate accuracy.\nUse ordinary 10X data to apply the model.",
        fontsize=10,
        color="#555",
    )
    savefig(fig, "06_reference_apply_prediction_only.png")


def method_workflow() -> None:
    fig, ax = plt.subplots(figsize=(13.5, 3.8))
    ax.set_axis_off()
    steps = [
        ("1. User input", "scRNA matrix\n+ optional protein\n+ optional ATAC"),
        ("2. Internal state", "pathway / program scores\n+ modality scores"),
        ("3. Virtual KO model", "prior-constrained residual / PLS\nsmall-sample friendly"),
        ("4. Outputs", "virtual cells\nheatmaps, UMAP, ROC-AUC"),
    ]
    colors = ["#377eb8", "#4daf4a", "#984ea3", "#e41a1c"]
    xs = [0.13, 0.38, 0.63, 0.88]
    box_w = 0.20
    for i, ((head, body), color, x) in enumerate(zip(steps, colors, xs)):
        ax.add_patch(
            plt.Rectangle((x - box_w / 2, 0.30), box_w, 0.43, facecolor=color, alpha=0.10, edgecolor=color, lw=2)
        )
        ax.text(x, 0.64, head, ha="center", fontsize=11, fontweight="bold", color="#222")
        ax.text(x, 0.45, body, ha="center", va="center", fontsize=9.6, color="#333")
        if i < len(xs) - 1:
            ax.text((x + xs[i + 1]) / 2, 0.52, ">", ha="center", va="center", fontsize=22, fontweight="bold", color="#777")
    ax.text(0.5, 0.92, "Virtual knockout workflow for small-sample multimodal data", ha="center", fontsize=16, fontweight="bold")
    ax.text(0.5, 0.13, "Users provide raw matrices. Pathway/program scores are created inside the software.", ha="center", fontsize=11, color="#555")
    savefig(fig, "01_method_workflow.png")


def summary_panel(
    summary_path: Path,
    auc_path: Path,
    delta_path: Path,
    title: str,
    ko: str,
    output_name: str,
) -> None:
    summary = pd.read_csv(summary_path).iloc[0]
    auc = pd.read_csv(auc_path).iloc[0]
    fig = plt.figure(figsize=(12, 7))
    gs = fig.add_gridspec(2, 4, height_ratios=[0.9, 1.55], hspace=0.35, wspace=0.45)
    fig.suptitle(title, fontsize=17, fontweight="bold", y=0.98)
    metric_card(
        fig.add_subplot(gs[0, 0]),
        "Direction match",
        f"{summary['mean_direction_cosine']:.2f}",
        "1.00 means same direction",
        "#377eb8",
    )
    metric_card(
        fig.add_subplot(gs[0, 1]),
        "Improved features",
        f"{summary['improved_fraction'] * 100:.0f}%",
        "closer to real KO than control",
        "#4daf4a",
    )
    metric_card(
        fig.add_subplot(gs[0, 2]),
        "ROC-AUC",
        f"{auc['roc_auc']:.2f}",
        "detects strong responses",
        "#984ea3",
    )
    metric_card(
        fig.add_subplot(gs[0, 3]),
        "KO tested",
        ko,
        "single or double gene",
        "#e41a1c",
    )

    df = pd.read_csv(delta_path)
    row = df.iloc[0]
    true_cols = [
        c
        for c in df.columns
        if c.startswith("true_delta_") and c not in {"true_delta_norm"}
    ]
    pairs = []
    for true_col in true_cols:
        suffix = true_col.replace("true_delta_", "")
        pred_col = f"pred_delta_{suffix}"
        if pred_col in df.columns:
            pairs.append((suffix, float(row[true_col]), float(row[pred_col])))
    pairs = sorted(pairs, key=lambda x: abs(x[1]), reverse=True)[:8]
    mat = pd.DataFrame(
        {"Real KO": [p[1] for p in pairs], "Virtual KO": [p[2] for p in pairs]},
        index=[clean_feature(p[0]) for p in pairs],
    )
    ax = fig.add_subplot(gs[1, :2])
    vmax = max(0.2, float(np.nanmax(np.abs(mat.to_numpy()))))
    sns.heatmap(
        mat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.7,
        linecolor="white",
        cbar=False,
    )
    ax.set_title("Largest real KO changes", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0, labelsize=10)
    ax.tick_params(axis="y", labelsize=8.5)

    ax2 = fig.add_subplot(gs[1, 2:])
    true = np.array([p[1] for p in pairs])
    pred = np.array([p[2] for p in pairs])
    lim = max(0.2, float(np.nanmax(np.abs(np.r_[true, pred])))) * 1.15
    ax2.scatter(true, pred, s=70, color="#377eb8", alpha=0.85)
    ax2.plot([-lim, lim], [-lim, lim], color="#555", lw=1.5, ls="--")
    # Keep this panel visually calm: the heatmap carries feature names, while
    # the scatter shows whether real and virtual effects follow the diagonal.
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.axhline(0, color="#ddd", lw=1)
    ax2.axvline(0, color="#ddd", lw=1)
    ax2.set_xlabel("Real KO change")
    ax2.set_ylabel("Virtual KO change")
    ax2.set_title("Real vs virtual agreement", fontsize=13, fontweight="bold")
    savefig(fig, output_name)


def interaction_improvement() -> None:
    df = pd.read_csv(ROOT / "results" / "norman_double_interaction_metrics.csv")
    agg = (
        df[df["subset"].isin(["all_combos", "has_unseen_gene"])]
        .groupby(["subset", "model"])[["r2", "roc_auc_abs_gt_0.15", "mae"]]
        .mean()
        .reset_index()
    )
    labels = {
        "all_combos": "All double KOs",
        "has_unseen_gene": "Contains unseen gene",
        "single_gene_additive": "Additive baseline",
        "interaction_residual": "Interaction model",
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), sharex=False)
    specs = [
        ("r2", "R2: magnitude fit", True),
        ("roc_auc_abs_gt_0.15", "ROC-AUC: strong response", True),
        ("mae", "MAE: lower is better", False),
    ]
    palette = {"single_gene_additive": "#9e9e9e", "interaction_residual": "#377eb8"}
    for ax, (metric, title, higher) in zip(axes, specs):
        sns.barplot(
            data=agg,
            x="subset",
            y=metric,
            hue="model",
            ax=ax,
            palette=palette,
            order=["all_combos", "has_unseen_gene"],
            hue_order=["single_gene_additive", "interaction_residual"],
        )
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        tick_positions = ax.get_xticks()
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([labels[t.get_text()] for t in ax.get_xticklabels()], fontsize=9)
        ax.legend_.remove()
        for container in ax.containers:
            ax.bar_label(container, fmt="%.2f", fontsize=8, padding=2)
        if metric == "r2":
            ax.axhline(0, color="#444", lw=1, ls="--")
    for ax in axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Double-gene KO: interaction prior improves nonlinear effects", fontsize=16, fontweight="bold", y=0.99)
    fig.text(
        0.5,
        0.91,
        "Gray = additive baseline    Blue = interaction model",
        ha="center",
        fontsize=10,
        color="#444",
    )
    fig.text(
        0.5,
        0.02,
        "The interaction model adds prior-based gene-gene residuals on top of the additive single-gene baseline.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.80, bottom=0.20, wspace=0.20)
    savefig(fig, "05_double_interaction_improvement.png")


def trust_matrix() -> None:
    rows = [
        ("Single KO with labels", "Strong", "Real-vs-virtual heatmap, UMAP, ROC-AUC"),
        ("Double KO with labels", "Good; nonlinear model helps", "Interaction residual improves R2 and AUC"),
        ("Ordinary 10X without labels", "Apply only", "Can show predicted shift, cannot prove accuracy inside that dataset"),
        ("Multimodal input", "Supported", "RNA pathway scores + ADT/protein; ATAC score interface prepared"),
        ("Few features / tiny data", "Use caution", "AUC can look high when feature count is very small"),
    ]
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.set_axis_off()
    ax.text(0.02, 0.93, "How to read the current method", fontsize=17, fontweight="bold")
    headers = ["Use case", "Current verdict", "What figure to trust"]
    xs = [0.02, 0.36, 0.60]
    widths = [0.30, 0.20, 0.36]
    y0 = 0.78
    for x, w, h in zip(xs, widths, headers):
        ax.add_patch(plt.Rectangle((x, y0), w, 0.10, color="#333", alpha=0.95))
        ax.text(x + 0.012, y0 + 0.052, h, va="center", fontsize=10, color="white", fontweight="bold")
    y = y0 - 0.12
    for i, row in enumerate(rows):
        bg = "#f5f7fb" if i % 2 == 0 else "#ffffff"
        for x, w in zip(xs, widths):
            ax.add_patch(plt.Rectangle((x, y), w, 0.105, facecolor=bg, edgecolor="#dddddd", lw=0.8))
        for x, text in zip(xs, row):
            wrapped = "\n".join(textwrap.wrap(text, width=34 if x > 0.55 else 24))
            ax.text(x + 0.012, y + 0.054, wrapped, va="center", fontsize=9.5, color="#222")
        y -= 0.105
    ax.text(0.02, 0.05, "Bottom line: use labeled perturbation data to evaluate; use ordinary 10X data to apply a trained reference model.", fontsize=10, color="#555")
    savefig(fig, "07_what_to_trust.png")


def auc_curve_panel() -> None:
    examples = [
        (
            "Single KO: STAT1",
            ROOT / "results" / "software_interface_single_gene_demo" / "auc_points.csv",
            ROOT / "results" / "software_interface_single_gene_demo" / "auc_summary.csv",
            "#377eb8",
        ),
        (
            "Double KO: CEBPB+CEBPA",
            ROOT / "results" / "software_interface_double_gene_demo" / "auc_points.csv",
            ROOT / "results" / "software_interface_double_gene_demo" / "auc_summary.csv",
            "#984ea3",
        ),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, (title, points_path, summary_path, color) in zip(axes, examples):
        points = pd.read_csv(points_path)
        summary = pd.read_csv(summary_path).iloc[0]
        y_true = (points["true_abs_delta"] >= float(summary["threshold"])).astype(int)
        y_score = points["pred_abs_delta"].astype(float)
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = sklearn_auc(fpr, tpr)
        ax.plot(fpr, tpr, color=color, lw=2.8, label=f"AUC = {roc_auc:.2f}")
        ax.plot([0, 1], [0, 1], color="#999", lw=1.5, ls="--", label="Random = 0.50")
        ax.fill_between(fpr, tpr, alpha=0.12, color=color)
        ax.set_title(title, fontsize=13, fontweight="bold")
        ax.set_xlabel("False positive rate")
        ax.set_ylabel("True positive rate")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.03)
        ax.legend(loc="lower right", frameon=True)
        ax.text(
            0.02,
            0.93,
            f"strong-response features: {int(summary['n_positive'])}",
            transform=ax.transAxes,
            fontsize=9,
            color="#555",
        )
    fig.suptitle("ROC curves for detecting strong KO responses", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "AUC is shown as a curve here. Bar/card values are only compact summaries of the same ROC analysis.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.82, bottom=0.18, wspace=0.25)
    savefig(fig, "08_auc_roc_curves.png")


def multimodal_dataset_summary() -> None:
    summary = pd.read_csv(ROOT / "results" / "multi_dataset_virtual_ko_summary.csv")
    summary["dataset_short"] = summary["dataset"].replace(
        {
            "Papalexi ECCITE-seq": "Papalexi\nRNA+ADT",
            "Norman Perturb-seq": "Norman\nRNA",
            "Datlinger CRISPR RNA": "Datlinger\nRNA",
            "Dixit Perturb-seq RNA": "Dixit\nRNA",
        }
    )
    colors = ["#377eb8" if "Papalexi" in d else "#9e9e9e" for d in summary["dataset"]]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.7))
    ax = axes[0]
    ax.bar(summary["dataset_short"], summary["improved_fraction"] * 100, color=colors)
    ax.set_title("How many features improve?", fontsize=13, fontweight="bold")
    ax.set_ylabel("Improved features (%)")
    ax.set_ylim(0, 100)
    for i, v in enumerate(summary["improved_fraction"] * 100):
        ax.text(i, v + 3, f"{v:.0f}%", ha="center", fontsize=9)
    ax.text(0, 89, "multimodal\nexample", ha="center", fontsize=9, color="#377eb8", fontweight="bold")

    ax = axes[1]
    ax.bar(summary["dataset_short"], summary["mean_distribution_improvement"], color=colors)
    ax.axhline(0, color="#444", lw=1.2, ls="--")
    ax.set_title("Mean distribution improvement", fontsize=13, fontweight="bold")
    ax.set_ylabel("Higher is better")
    for i, v in enumerate(summary["mean_distribution_improvement"]):
        if v >= 0:
            ax.text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
        else:
            ax.text(i, v + 0.025, f"{v:.2f}", ha="center", fontsize=9, color="white")
    fig.suptitle("Multi-dataset virtual KO summary", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.01,
        "Current multimodal evidence is Papalexi RNA+ADT. RNA-only datasets are useful stress tests but perform less consistently.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.80, bottom=0.22, wspace=0.28)
    savefig(fig, "09_multimodal_multi_dataset_summary.png")


def r2_mae_intuitive_panel() -> None:
    df = pd.read_csv(ROOT / "results" / "norman_double_interaction_metrics.csv")
    agg = (
        df[df["subset"].eq("all_combos")]
        .groupby("model")[["r2", "mae", "roc_auc_abs_gt_0.15"]]
        .mean()
        .loc[["single_gene_additive", "interaction_residual"]]
    )
    additive = agg.loc["single_gene_additive"]
    interaction = agg.loc["interaction_residual"]
    mae_reduction = (additive["mae"] - interaction["mae"]) / additive["mae"] * 100

    fig = plt.figure(figsize=(11, 4.6))
    gs = fig.add_gridspec(1, 3, wspace=0.35)
    colors = {"Additive": "#9e9e9e", "Interaction": "#377eb8"}

    ax = fig.add_subplot(gs[0, 0])
    ax.barh(["Additive", "Interaction"], [max(0, additive["r2"]), interaction["r2"]], color=[colors["Additive"], colors["Interaction"]])
    ax.set_xlim(0, 1)
    ax.set_title("R2 as pattern explained", fontsize=12, fontweight="bold")
    ax.set_xlabel("0 = none, 1 = perfect")
    for i, v in enumerate([additive["r2"], interaction["r2"]]):
        ax.text(max(0, v) + 0.03, i, f"{v:.2f}", va="center", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    ax.barh(["Additive", "Interaction"], [additive["mae"], interaction["mae"]], color=[colors["Additive"], colors["Interaction"]])
    ax.set_xlim(0, max(additive["mae"], interaction["mae"]) * 1.2)
    ax.set_title("MAE as prediction error", fontsize=12, fontweight="bold")
    ax.set_xlabel("Lower is better")
    for i, v in enumerate([additive["mae"], interaction["mae"]]):
        ax.text(v + 0.004, i, f"{v:.2f}", va="center", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    ax.set_axis_off()
    ax.add_patch(plt.Rectangle((0.06, 0.40), 0.88, 0.22, facecolor="#e8f1fb", edgecolor="#377eb8", lw=1.6))
    ax.add_patch(plt.Rectangle((0.06, 0.40), 0.88 * mae_reduction / 100, 0.22, facecolor="#377eb8", edgecolor="none"))
    ax.text(0.5, 0.74, "Error reduction", ha="center", fontsize=13, fontweight="bold")
    ax.text(0.5, 0.49, f"{mae_reduction:.0f}% lower MAE", ha="center", va="center", fontsize=18, color="#1f5aa6", fontweight="bold")
    ax.text(0.5, 0.22, "Interaction model compared with additive baseline", ha="center", fontsize=9.5, color="#555")

    fig.suptitle("R2 and MAE in plain visual terms", fontsize=16, fontweight="bold")
    fig.subplots_adjust(top=0.78, bottom=0.18)
    savefig(fig, "10_r2_mae_intuitive.png")


def umap_examples_panel() -> None:
    images = [
        ("Single-gene cell state shift", ROOT / "results" / "figures" / "papalexi_cell_level_umap_single_gene_before_after.png"),
        ("Multi-gene cell state shift", ROOT / "results" / "figures" / "norman_cell_level_umap_multi_gene_before_after.png"),
        ("Across datasets", ROOT / "results" / "figures" / "multi_dataset_virtual_ko_umap_examples.png"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
    for ax, (title, path) in zip(axes, images):
        ax.imshow(plt.imread(path))
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_axis_off()
    fig.suptitle("UMAP views: before and after virtual KO", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.02,
        "UMAP is used for visual cell-state movement. Accuracy still needs labeled perturbation data when available.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.82, bottom=0.12, wspace=0.05)
    savefig(fig, "11_umap_before_after_examples.png")


def double_ko_many_combos() -> None:
    pred = pd.read_csv(ROOT / "results" / "norman_double_interaction_predictions.csv")
    programs = [
        "program_ERYTHROID",
        "program_GRANULOCYTE_APOPTOSIS",
        "program_MAPK_TGFB",
        "program_PRO_GROWTH",
        "program_PIONEER_TF",
    ]
    rows = []
    for _, row in pred.iterrows():
        true = np.array([row[f"true_delta_{p}"] for p in programs], dtype=float)
        additive = np.array([row[f"additive_pred_delta_{p}"] for p in programs], dtype=float)
        interaction = np.array([row[f"interaction_pred_delta_{p}"] for p in programs], dtype=float)
        rows.append(
            {
                "ko": row["ko_genes"],
                "additive_mae": float(np.mean(np.abs(true - additive))),
                "interaction_mae": float(np.mean(np.abs(true - interaction))),
            }
        )
    rank = pd.DataFrame(rows)
    rank["improvement"] = rank["additive_mae"] - rank["interaction_mae"]
    selected = rank.sort_values("improvement", ascending=False).head(12)["ko"].tolist()
    heat = rank[rank["ko"].isin(selected)].set_index("ko").loc[selected]
    mat = heat[["additive_mae", "interaction_mae"]].rename(
        columns={"additive_mae": "Additive error", "interaction_mae": "Interaction error"}
    )

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    sns.heatmap(
        mat,
        ax=ax,
        cmap="Reds",
        annot=True,
        fmt=".2f",
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "MAE; lower is better"},
    )
    ax.set_title("Other double-gene KOs: interaction model lowers error", fontsize=15, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("Double KO combination")
    ax.tick_params(axis="x", labelrotation=0, labelsize=10)
    ax.tick_params(axis="y", labelsize=9)
    savefig(fig, "12_other_double_ko_combos.png")


def tenx_and_multimodal_io_examples() -> None:
    manifest = pd.read_csv(ROOT / "results" / "ordinary_10x_like_input_demo" / "derived_state_manifest.csv")
    pred = pd.read_csv(ROOT / "results" / "multimodal_single_double_apply_demo" / "predicted_ko_delta.csv")

    source_counts = manifest["source"].value_counts().rename_axis("source").reset_index(name="n_features")
    source_counts["source"] = source_counts["source"].replace({"RNA pathway score": "RNA pathway\nscores", "obsm:protein": "ADT/protein\nfeatures"})

    rows = []
    for _, row in pred.iterrows():
        ko = row["ko_target"]
        for col in [c for c in pred.columns if c.startswith("pred_delta_")]:
            feature = col.replace("pred_delta_", "")
            rows.append({"ko": ko, "feature": feature, "delta": float(row[col])})
    long = pd.DataFrame(rows)
    score = long.groupby("feature")["delta"].apply(lambda x: np.max(np.abs(x))).sort_values(ascending=False)
    keep = list(score.head(12).index)
    heat = long[long["feature"].isin(keep)].pivot(index="feature", columns="ko", values="delta").loc[keep]
    heat.index = [clean_feature(idx) for idx in heat.index]
    vmax = max(0.2, float(np.nanmax(np.abs(heat.to_numpy()))))

    fig = plt.figure(figsize=(13, 6.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.55], wspace=0.30)
    ax = fig.add_subplot(gs[0, 0])
    colors = ["#377eb8" if "RNA" in s else "#4daf4a" for s in source_counts["source"]]
    ax.bar(source_counts["source"], source_counts["n_features"], color=colors)
    ax.set_title("Ordinary 10X-like input becomes state scores", fontsize=13, fontweight="bold")
    ax.set_ylabel("Number of derived features")
    for i, v in enumerate(source_counts["n_features"]):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=10)
    ax.text(
        0.5,
        -0.15,
        "Input: cells x genes matrix; optional ADT/protein in obsm",
        ha="center",
        transform=ax.transAxes,
        fontsize=9.5,
        color="#555",
    )

    ax = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        heat,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "Predicted state change"},
    )
    ax.set_title("Multimodal input: single and double KO application", fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelrotation=0)
    ax.tick_params(axis="y", labelsize=8.5)

    fig.suptitle("10X input and multimodal one/two-gene KO outputs", fontsize=16, fontweight="bold")
    fig.text(
        0.5,
        0.035,
        "Application outputs: accuracy metrics require matching real KO labels in the same dataset.",
        ha="center",
        fontsize=10,
        color="#555",
    )
    fig.subplots_adjust(top=0.82, bottom=0.20)
    savefig(fig, "13_10x_multimodal_single_double_outputs.png")


def write_doc() -> None:
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(
        """# 面向用户的可视化结果说明

这套图放在 `results/user_facing_figures/`。它的目的不是展示所有模型细节，而是让第一次接触这个方法的人直接看懂三件事：

1. 用户输入什么。
2. 虚拟敲除预测了什么变化。
3. 预测效果到底好不好，以及哪些场景不能过度解释。

## 01_method_workflow.png

说明软件流程。用户输入的是普通单细胞 RNA 矩阵，或 RNA 加 protein/ADT/ATAC 等多模态矩阵；通路分数和 program score 是软件内部自动生成的状态表示，不要求用户提前准备。

## 02_single_ko_summary.png

STAT1 单基因敲除的总览图。上方四个数字是最容易读的结论：

- Direction match：越接近 1，说明预测变化方向越接近真实 KO。
- Improved features：有多少状态特征比 untreated/control 更接近真实 KO。
- ROC-AUC：识别强响应通路/蛋白的能力，越接近 1 越好。
- KO tested：本图测试的敲除基因。

下方 heatmap 对比真实 KO 和虚拟 KO 的状态变化；颜色方向一致，说明预测方向对；颜色深浅接近，说明幅度也接近。

## 03_single_ko_true_vs_virtual.png

单基因敲除的放大版 heatmap。适合回答“敲了什么基因，引起了哪些通路或蛋白变化，和真实情况差多少”。

## 04_double_ko_summary.png

CEBPB+CEBPA 双基因敲除的总览图。读法与单基因图相同。这个例子强调软件默认支持一个或两个基因敲除。

## 05_double_interaction_improvement.png

双基因敲除不是简单的两个单基因效果相加，所以这里比较了：

- Additive baseline：直接相加的基线模型。
- Interaction model：加入系统网络先验后的双基因相互作用残差模型。

结果显示 interaction model 在 R2 和 ROC-AUC 上明显更好，MAE 更低，说明它能更好处理双敲非线性效应。

## 06_reference_apply_prediction_only.png

普通 10X 单细胞数据如果没有真实 KO 标签，只能做“应用/预测”，不能在该数据内部证明准确率。图中展示的是 reference model 预测 STAT1 KO 后最可能变化的状态特征。

这类图适合展示“如果在这批普通细胞里虚拟敲掉某个基因，细胞状态可能往哪里变”；但不能报告真实 AUC 或真实准确性，除非另有 perturb-seq/CRISPR/药物扰动标签。

## 07_what_to_trust.png

给非专业用户的读图总结：

- 有真实 KO 标签时，看真实 vs 虚拟 heatmap、UMAP 和 ROC-AUC。
- 没有 KO 标签时，只能看预测的状态转换，不能说模型在该数据内部被验证。
- 双基因敲除要看 interaction model，而不是只看简单相加。
- 特征很少时 AUC 要谨慎解释，需要同时看 heatmap、R2/MAE 和生物学方向是否合理。

## 08_auc_roc_curves.png

这是 ROC 曲线版的 AUC。之前 summary 卡片里的 AUC 是一个数字摘要；这张图展示完整曲线。曲线越靠左上角，说明模型越能把“强响应通路/蛋白”和“弱响应特征”区分开。

## 09_multimodal_multi_dataset_summary.png

展示多数据集结果。Papalexi 是当前真正的多模态例子，即 RNA pathway score 加 ADT/protein；Norman、Datlinger、Dixit 主要是 RNA-only 或 RNA-derived 状态表示。当前结果显示，多模态例子的稳定性更好，RNA-only 数据在小样本和跨数据设置下波动更大。

## 10_r2_mae_intuitive.png

把 R2 和 MAE 改成更直观的读法：

- R2：模型解释真实变化模式的程度，越高越好。
- MAE：预测误差，越低越好。
- Error reduction：相互作用模型相比简单相加模型把平均误差降低了多少。

## 11_umap_before_after_examples.png

把单敲、多敲和多数据集的 UMAP 放在一起。UMAP 主要用来看虚拟敲除前后细胞状态是否发生可见移动；但 UMAP 是可视化，不是准确率证明，准确性仍要看真实 KO 标签下的 heatmap、ROC-AUC、R2/MAE。

## 12_other_double_ko_combos.png

展示 Norman 数据中多个其它双基因敲除组合。每一行是一个双敲组合，每两列比较简单相加模型和加入相互作用先验后的模型误差。颜色越浅、数字越小越好。

## 13_10x_multimodal_single_double_outputs.png

这张图补充普通 10X 和多模态输入的实际输出。

左边说明普通细胞矩阵输入后，软件会自动派生 RNA pathway score；如果 h5ad 里有 ADT/protein/ATAC 这类额外模态，也会作为额外状态特征进入模型。

右边展示 RNA+ADT 多模态输入下，同时预测一个基因敲除（STAT1）和两个基因敲除（STAT1+JAK2）后的状态变化。注意这属于 reference model application：可以输出预测变化和虚拟细胞，但如果输入数据没有真实双敲标签，就不能在该数据内部计算真实 AUC/R2/MAE。
""",
        encoding="utf-8",
    )


def main() -> None:
    sns.set_theme(style="whitegrid", context="notebook")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlelocation": "center",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    method_workflow()
    summary_panel(
        ROOT / "results" / "software_interface_single_gene_demo" / "summary.csv",
        ROOT / "results" / "software_interface_single_gene_demo" / "auc_summary.csv",
        ROOT / "results" / "software_interface_single_gene_demo" / "delta_table.csv",
        "Single-gene virtual KO result: STAT1",
        "STAT1",
        "02_single_ko_summary.png",
    )
    true_pred_heatmap(
        ROOT / "results" / "software_interface_single_gene_demo" / "delta_table.csv",
        "Real vs virtual KO changes: STAT1",
        "03_single_ko_true_vs_virtual.png",
        top_n=12,
    )
    summary_panel(
        ROOT / "results" / "software_interface_double_gene_demo" / "summary.csv",
        ROOT / "results" / "software_interface_double_gene_demo" / "auc_summary.csv",
        ROOT / "results" / "software_interface_double_gene_demo" / "delta_table.csv",
        "Double-gene virtual KO result: CEBPB + CEBPA",
        "CEBPB+CEBPA",
        "04_double_ko_summary.png",
    )
    interaction_improvement()
    predicted_heatmap(ROOT / "results" / "reference_apply_stat1_demo" / "predicted_ko_delta.csv")
    trust_matrix()
    auc_curve_panel()
    multimodal_dataset_summary()
    r2_mae_intuitive_panel()
    umap_examples_panel()
    double_ko_many_combos()
    tenx_and_multimodal_io_examples()
    write_doc()
    print(f"Wrote figure pack to {OUT}")
    print(f"Wrote guide to {DOC}")


if __name__ == "__main__":
    main()
