from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns


FIG_DIR = Path("results/figures")


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


def feature_label(feature: str) -> str:
    for prefix in ["pathway_", "protein_", "program_", "svd_state_"]:
        if feature.startswith(prefix):
            return feature.replace(prefix, "")
    return feature


def dataset_short(dataset: str) -> str:
    return {
        "Papalexi ECCITE-seq": "Papalexi",
        "Norman Perturb-seq": "Norman",
        "Datlinger CRISPR RNA": "Datlinger",
        "Dixit Perturb-seq RNA": "Dixit",
    }.get(dataset, dataset)


def feature_subset(dataset: str, features: list[str]) -> list[str]:
    if dataset == "Papalexi ECCITE-seq":
        preferred = [
            "pathway_IFNG_JAK_STAT",
            "pathway_IMMUNE_CHECKPOINT",
            "pathway_NRF2_STRESS",
            "pathway_CELL_CYCLE_G2M",
            "protein_PDL1",
            "protein_CD86",
            "protein_PDL2",
            "protein_CD366",
        ]
        return [feature for feature in preferred if feature in features]
    if dataset == "Norman Perturb-seq":
        preferred = [
            "program_ERYTHROID",
            "program_GRANULOCYTE_APOPTOSIS",
            "program_MAPK_TGFB",
            "program_PRO_GROWTH",
            "program_PIONEER_TF",
        ]
        return [feature for feature in preferred if feature in features]
    return features[:6]


def mean_state_tables(cells: pd.DataFrame, dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sub = cells.loc[cells["dataset"] == dataset].copy()
    feature_cols = [
        col
        for col in sub.columns
        if col.startswith("pathway_") or col.startswith("protein_") or col.startswith("program_") or col.startswith("svd_state_")
    ]
    feature_cols = [col for col in feature_cols if sub[col].notna().any()]
    feature_cols = feature_subset(dataset, feature_cols)
    means = sub.groupby(["ko_target", "state"], observed=True)[feature_cols].mean()
    true_delta = []
    virtual_delta = []
    for ko in sub["ko_target"].drop_duplicates():
        if (ko, "control cells") not in means.index:
            continue
        control = means.loc[(ko, "control cells")]
        true = means.loc[(ko, "true KO cells")] - control
        virtual = means.loc[(ko, "virtual KO cells")] - control
        true_delta.append(pd.Series(true.to_numpy(), index=feature_cols, name=ko))
        virtual_delta.append(pd.Series(virtual.to_numpy(), index=feature_cols, name=ko))
    true_df = pd.DataFrame(true_delta)
    virtual_df = pd.DataFrame(virtual_delta)
    error_df = virtual_df - true_df
    true_df.columns = [feature_label(col) for col in true_df.columns]
    virtual_df.columns = [feature_label(col) for col in virtual_df.columns]
    error_df.columns = [feature_label(col) for col in error_df.columns]
    return true_df, virtual_df, error_df


def plot_dataset_heatmap(cells: pd.DataFrame, dataset: str) -> None:
    true_df, virtual_df, error_df = mean_state_tables(cells, dataset)
    if true_df.empty:
        return
    vmax = np.nanmax(np.abs(pd.concat([true_df, virtual_df]).to_numpy()))
    err_vmax = np.nanmax(np.abs(error_df.to_numpy()))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, max(3.2, 0.55 * len(true_df) + 2.0)), sharey=True)
    sns.heatmap(true_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[0], cbar=False)
    sns.heatmap(virtual_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[1], cbar=False)
    sns.heatmap(error_df, cmap="vlag", center=0, vmin=-err_vmax, vmax=err_vmax, annot=True, fmt=".2f", ax=axes[2], cbar_kws={"label": "virtual - true"})
    axes[0].set_title("Real KO change")
    axes[1].set_title("Virtual KO change")
    axes[2].set_title("Prediction error")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=35)
    fig.suptitle(f"{dataset}: real vs virtual knockout effects", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"multi_dataset_heatmap_{dataset_short(dataset).lower()}_true_virtual_error.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_combined_direction_heatmap(cells: pd.DataFrame) -> None:
    rows = []
    for dataset in cells["dataset"].drop_duplicates():
        true_df, virtual_df, error_df = mean_state_tables(cells, dataset)
        for ko in true_df.index:
            true_vec = true_df.loc[ko].to_numpy(dtype=float)
            virt_vec = virtual_df.loc[ko].to_numpy(dtype=float)
            denom = (np.linalg.norm(true_vec) * np.linalg.norm(virt_vec)) + 1e-9
            cosine = float(np.dot(true_vec, virt_vec) / denom)
            mae = float(np.mean(np.abs(virt_vec - true_vec)))
            rows.append(
                {
                    "dataset": dataset,
                    "ko_target": ko,
                    "label": f"{ko} ({dataset_short(dataset)})",
                    "direction_cosine": cosine,
                    "mean_abs_delta_error": mae,
                }
            )
    summary = pd.DataFrame(rows)
    summary.to_csv("results/multi_dataset_true_vs_virtual_heatmap_summary.csv", index=False)
    heat = summary.set_index("label")[["direction_cosine", "mean_abs_delta_error"]]
    fig, axes = plt.subplots(1, 2, figsize=(9, max(4.5, 0.42 * len(heat) + 1.8)), sharey=True)
    sns.heatmap(heat[["direction_cosine"]], cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", ax=axes[0], cbar_kws={"label": "cosine"})
    vmax = np.nanmax(heat[["mean_abs_delta_error"]].to_numpy())
    sns.heatmap(heat[["mean_abs_delta_error"]], cmap="rocket_r", vmin=0, vmax=vmax, annot=True, fmt=".2f", ax=axes[1], cbar_kws={"label": "MAE"})
    axes[0].set_title("Direction agreement")
    axes[1].set_title("Magnitude error")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Real vs virtual KO effect agreement across datasets", y=1.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "multi_dataset_true_vs_virtual_agreement_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def write_doc(datasets: list[str]) -> None:
    lines = [f"- `results/figures/multi_dataset_heatmap_{dataset_short(dataset).lower()}_true_virtual_error.png`" for dataset in datasets]
    text = f"""# True KO vs virtual KO heatmap visualization

这组图专门展示：

```text
真实敲除后状态变化
vs
虚拟敲除预测状态变化
vs
预测误差
```

每个数据集都有一张三栏 heatmap：

- Real KO change：真实 KO 相对 control 的变化。
- Virtual KO change：虚拟 KO 相对 control 的变化。
- Prediction error：virtual - true。

另外有一张总览图：

- `results/figures/multi_dataset_true_vs_virtual_agreement_heatmap.png`

单数据集 heatmap：

{chr(10).join(lines)}

读图规则：

- 左右两张热图颜色方向一致，说明模型预测到了正确变化方向。
- 误差热图颜色越浅，说明虚拟 KO 越接近真实 KO。
- Direction agreement 越接近 1，说明整体变化方向越一致。
"""
    Path("docs/multi_dataset_true_vs_virtual_heatmaps.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    cells = pd.read_csv("results/multi_dataset_virtual_ko_cells.csv")
    datasets = list(cells["dataset"].drop_duplicates())
    for dataset in datasets:
        plot_dataset_heatmap(cells, dataset)
    plot_combined_direction_heatmap(cells)
    write_doc(datasets)
    print("Saved multi-dataset true-vs-virtual KO heatmaps.")


if __name__ == "__main__":
    main()
