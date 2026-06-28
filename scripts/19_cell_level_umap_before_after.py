from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import StandardScaler
import umap


FIG_DIR = Path("results/figures")


SINGLE_GENE_KOS = ["STAT1", "JAK2", "IFNGR2", "IRF1"]
MULTI_GENE_KOS = ["AHR+KLF1", "CEBPB+CEBPA", "MAPK1+TGFBR2", "CBL+UBASH3B"]


PALETTE = {
    "control cells": "#BDBDBD",
    "virtual KO cells": "#E76F51",
    "true KO cells": "#2A9D8F",
}


def setup_plot() -> None:
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


def control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("nt|control|non|safe|neg")


def state_columns(frame: pd.DataFrame) -> list[str]:
    return [
        col
        for col in frame.columns
        if col.startswith("pathway_") or col.startswith("protein_") or col.startswith("program_")
    ]


def fit_umap(cells: pd.DataFrame, features: list[str], seed: int) -> pd.DataFrame:
    matrix = cells[features].to_numpy(dtype=float)
    scaled = StandardScaler().fit_transform(matrix)
    reducer = umap.UMAP(
        n_neighbors=min(30, max(5, len(cells) // 25)),
        min_dist=0.25,
        metric="euclidean",
        random_state=seed,
    )
    coords = reducer.fit_transform(scaled)
    out = cells.copy()
    out["UMAP1"] = coords[:, 0]
    out["UMAP2"] = coords[:, 1]
    return out


def add_centroid_arrows(ax, plot: pd.DataFrame) -> None:
    centers = plot.groupby("state", observed=True)[["UMAP1", "UMAP2"]].mean()
    if "control cells" not in centers.index:
        return
    start = centers.loc["control cells"]
    ax.scatter(
        start["UMAP1"],
        start["UMAP2"],
        s=110,
        marker="X",
        color="#6E6E6E",
        edgecolor="white",
        linewidth=0.8,
        zorder=8,
    )
    for state, color, linestyle in [
        ("virtual KO cells", PALETTE["virtual KO cells"], "-"),
        ("true KO cells", PALETTE["true KO cells"], "--"),
    ]:
        if state not in centers.index:
            continue
        end = centers.loc[state]
        ax.scatter(
            end["UMAP1"],
            end["UMAP2"],
            s=120,
            marker="X",
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=9,
        )
        ax.annotate(
            "",
            xy=(end["UMAP1"], end["UMAP2"]),
            xytext=(start["UMAP1"], start["UMAP2"]),
            arrowprops={
                "arrowstyle": "->",
                "color": color,
                "lw": 3.2,
                "linestyle": linestyle,
                "shrinkA": 2,
                "shrinkB": 2,
                "mutation_scale": 17,
            },
            zorder=10,
        )


def centroid_metrics(cells: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rows = []
    for ko, group in cells.groupby("ko_target", observed=True):
        centers = group.groupby("state", observed=True)[["UMAP1", "UMAP2"]].mean()
        if not {"control cells", "virtual KO cells", "true KO cells"}.issubset(set(centers.index)):
            continue
        control_to_true = float(np.linalg.norm(centers.loc["control cells"] - centers.loc["true KO cells"]))
        virtual_to_true = float(np.linalg.norm(centers.loc["virtual KO cells"] - centers.loc["true KO cells"]))
        rows.append(
            {
                "dataset": dataset,
                "ko_target": ko,
                "control_to_true_umap_distance": control_to_true,
                "virtual_to_true_umap_distance": virtual_to_true,
                "umap_centroid_improvement": 1.0 - virtual_to_true / control_to_true if control_to_true > 1e-9 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_before_after_umap(cells: pd.DataFrame, title: str, out_path: Path) -> None:
    kos = list(cells["ko_target"].drop_duplicates())
    fig, axes = plt.subplots(len(kos), 2, figsize=(11.8, 3.2 * len(kos)), sharex=True, sharey=True)
    if len(kos) == 1:
        axes = np.asarray([axes])

    x_pad = (cells["UMAP1"].max() - cells["UMAP1"].min()) * 0.08
    y_pad = (cells["UMAP2"].max() - cells["UMAP2"].min()) * 0.08
    xlim = (cells["UMAP1"].min() - x_pad, cells["UMAP1"].max() + x_pad)
    ylim = (cells["UMAP2"].min() - y_pad, cells["UMAP2"].max() + y_pad)

    for row, ko in enumerate(kos):
        group = cells.loc[cells["ko_target"] == ko]
        before = group.loc[group["state"] == "control cells"]
        after = group.loc[group["state"].isin(["control cells", "virtual KO cells", "true KO cells"])]

        ax = axes[row, 0]
        sns.scatterplot(
            data=before,
            x="UMAP1",
            y="UMAP2",
            color=PALETTE["control cells"],
            s=28,
            alpha=0.75,
            linewidth=0,
            ax=ax,
        )
        ax.set_title(f"{ko}: before KO\ncontrol cells")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")

        ax = axes[row, 1]
        for state, alpha, size in [
            ("control cells", 0.20, 20),
            ("virtual KO cells", 0.72, 30),
            ("true KO cells", 0.72, 30),
        ]:
            sub = after.loc[after["state"] == state]
            sns.scatterplot(
                data=sub,
                x="UMAP1",
                y="UMAP2",
                color=PALETTE[state],
                label=state if row == 0 else None,
                s=size,
                alpha=alpha,
                linewidth=0,
                ax=ax,
            )
        add_centroid_arrows(ax, after)
        ax.set_title(f"{ko}: after KO\nvirtual vs real")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("")

    handles, labels = axes[0, 1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, title="", loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995))
        axes[0, 1].legend_.remove()
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_improvement_summary(metrics: pd.DataFrame) -> None:
    plot = metrics.copy()
    plot["kind"] = np.where(plot["dataset"].str.contains("multi", case=False), "multi", "single")
    plot["label"] = plot["ko_target"] + " (" + plot["kind"] + ")"
    plot = plot.sort_values("umap_centroid_improvement")
    plt.figure(figsize=(9.2, 5.4))
    colors = np.where(plot["umap_centroid_improvement"] >= 0, "#457B9D", "#C65D4B")
    ax = plt.barh(plot["label"], plot["umap_centroid_improvement"], color=colors)
    plt.axvline(0, color="0.25", linewidth=1)
    plt.xlabel("UMAP centroid improvement\n(>0 means virtual KO centroid is closer to true KO)")
    plt.ylabel("")
    plt.title("Before/after UMAP: which KOs move in the right direction?")
    plt.tick_params(axis="y", labelsize=11)
    for bar, value in zip(ax, plot["umap_centroid_improvement"]):
        x = bar.get_width()
        plt.text(
            x + (0.02 if x >= 0 else -0.02),
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}",
            va="center",
            ha="left" if x >= 0 else "right",
            fontsize=10,
        )
    plt.tight_layout()
    plt.savefig(FIG_DIR / "cell_level_umap_before_after_improvement_summary.png", bbox_inches="tight", dpi=300)
    plt.close()


def load_single_gene_cells() -> tuple[pd.DataFrame, list[str]]:
    cells = pd.read_csv("results/papalexi_cell_level_residual_generator_cells.csv")
    cells = cells.loc[cells["ko_target"].isin(SINGLE_GENE_KOS)].copy()
    features = state_columns(cells)
    keep = ["ko_target", "state"] + features
    cells = cells[keep].dropna()
    cells["dataset"] = "Papalexi single-gene KO"
    return cells, features


def load_norman_frame() -> tuple[pd.DataFrame, list[str]]:
    adata = ad.read_h5ad("data/norman_small_program.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_genes"].astype(str).values
    for col in [col for col in obs.columns if col.startswith("program_")]:
        frame[col] = obs[col].astype(float).values
    features = [col for col in frame.columns if col.startswith("program_")]
    return frame, features


def build_multi_gene_cells(seed: int = 17) -> tuple[pd.DataFrame, list[str]]:
    rng = np.random.default_rng(seed)
    frame, features = load_norman_frame()
    pred = pd.read_csv("results/norman_system_prior_predictions.csv")
    control = frame.loc[control_mask(frame["ko_target"]), features].to_numpy(dtype=float)
    rows = []
    for ko in MULTI_GENE_KOS:
        true = frame.loc[frame["ko_target"].astype(str) == ko, features].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        pred_row = pred.loc[pred["ko_genes"].astype(str) == ko]
        if pred_row.empty:
            continue
        n = len(true)
        ctrl = control[rng.integers(0, len(control), size=n)]
        delta = np.asarray([float(pred_row[f"system_pred_delta_{feature}"].iloc[0]) for feature in features])
        virtual = ctrl + delta.reshape(1, -1)
        for state, matrix in [
            ("control cells", ctrl),
            ("virtual KO cells", virtual),
            ("true KO cells", true),
        ]:
            tmp = pd.DataFrame(matrix, columns=features)
            tmp["ko_target"] = ko
            tmp["state"] = state
            tmp["dataset"] = "Norman multi-gene KO"
            rows.append(tmp)
    return pd.concat(rows, ignore_index=True), features


def write_doc(single_metrics: pd.DataFrame, multi_metrics: pd.DataFrame) -> None:
    all_metrics = pd.concat([single_metrics, multi_metrics], ignore_index=True)
    single_mean = single_metrics["umap_centroid_improvement"].mean()
    multi_mean = multi_metrics["umap_centroid_improvement"].mean()
    rounded = all_metrics.round(3)
    table_lines = [
        "| dataset | KO | control->true | virtual->true | improvement |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, row in rounded.iterrows():
        table_lines.append(
            "| "
            f"{row['dataset']} | {row['ko_target']} | "
            f"{row['control_to_true_umap_distance']} | "
            f"{row['virtual_to_true_umap_distance']} | "
            f"{row['umap_centroid_improvement']} |"
        )
    table = "\n".join(table_lines)
    text = f"""# Cell-level UMAP：敲除前后单细胞状态移动

这组图专门回答一个直观问题：

```text
敲除前的 control cells 在哪里？
虚拟敲除后的 cells 移动到哪里？
真实实验 KO cells 在哪里？
```

## 图怎么看

每个 KO 都有两列：

- 左列：敲除前，只显示 control cells。
- 右列：敲除后，显示 virtual KO cells 和 true KO cells，同时淡灰色保留 control cells 作为参照。

箭头含义：

- 橙色实线箭头：control 质心 -> virtual KO 质心。
- 绿色虚线箭头：control 质心 -> true KO 质心。

如果橙色箭头和绿色箭头方向接近，并且橙色点云靠近绿色点云，说明虚拟敲除捕捉到了单细胞状态变化。

## 单基因敲除

图：`results/figures/papalexi_cell_level_umap_single_gene_before_after.png`

使用 Papalexi 多模态数据：pathway score + protein score。

平均 UMAP 质心改进：`{single_mean:.3f}`

## 多基因组合敲除

图：`results/figures/norman_cell_level_umap_multi_gene_before_after.png`

使用 Norman 组合扰动数据：gene program score。

平均 UMAP 质心改进：`{multi_mean:.3f}`

## 辅助解释表

注意：UMAP 质心距离只用于解释可视化，不作为正式性能指标。正式评价仍然看前面的 Wasserstein、ROC-AUC、R2。

{table}

## 当前结论

这两张 UMAP 图能让人直接看到：模型不是只输出一个分数，而是在 cell-level 状态空间里把 control cells 移向 KO-like 状态。

但也要诚实：如果某些 KO 的橙色云团没有靠近绿色云团，说明当前稳定基线还不够，需要在这个基础上继续接 VAE / flow matching / diffusion，而不是直接从零训练复杂模型。
"""
    Path("docs/cell_level_umap_before_after.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()

    single_cells, single_features = load_single_gene_cells()
    single_umap = fit_umap(single_cells, single_features, seed=41)
    single_umap.to_csv("results/papalexi_cell_level_umap_single_gene_cells.csv", index=False)
    plot_before_after_umap(
        single_umap,
        "Single-gene virtual KO: before and after in cell-level UMAP",
        FIG_DIR / "papalexi_cell_level_umap_single_gene_before_after.png",
    )
    single_metrics = centroid_metrics(single_umap, "Papalexi single-gene KO")

    multi_cells, multi_features = build_multi_gene_cells()
    multi_umap = fit_umap(multi_cells, multi_features, seed=43)
    multi_umap.to_csv("results/norman_cell_level_umap_multi_gene_cells.csv", index=False)
    plot_before_after_umap(
        multi_umap,
        "Multi-gene virtual KO: before and after in cell-level UMAP",
        FIG_DIR / "norman_cell_level_umap_multi_gene_before_after.png",
    )
    multi_metrics = centroid_metrics(multi_umap, "Norman multi-gene KO")

    metrics = pd.concat([single_metrics, multi_metrics], ignore_index=True)
    metrics.to_csv("results/cell_level_umap_before_after_metrics.csv", index=False)
    plot_improvement_summary(metrics)
    write_doc(single_metrics, multi_metrics)

    print(
        "single-gene mean UMAP centroid improvement="
        f"{single_metrics['umap_centroid_improvement'].mean():.3f}"
    )
    print(
        "multi-gene mean UMAP centroid improvement="
        f"{multi_metrics['umap_centroid_improvement'].mean():.3f}"
    )
    print("Saved cell-level before/after UMAP figures, tables, and docs.")


if __name__ == "__main__":
    main()
