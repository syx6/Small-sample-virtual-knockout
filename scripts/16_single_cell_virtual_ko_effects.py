from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


FIG_DIR = Path("results/figures")


PAPALEXI_FEATURES = {
    "pathway_IFNG_JAK_STAT": "IFNG-JAK-STAT pathway",
    "pathway_IMMUNE_CHECKPOINT": "Immune checkpoint pathway",
    "protein_PDL1": "PDL1 protein",
    "protein_CD86": "CD86 protein",
}


NORMAN_FEATURES = {
    "program_ERYTHROID": "Erythroid program",
    "program_GRANULOCYTE_APOPTOSIS": "Granulocyte/apoptosis program",
    "program_MAPK_TGFB": "MAPK/TGFB program",
    "program_PRO_GROWTH": "Pro-growth program",
}


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


def control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("nt|control|non|safe|neg")


def papalexi_feature_frame() -> pd.DataFrame:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_target"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("pathway_")]:
        frame[col] = obs[col].astype(float).values
    if "protein" in adata.obsm:
        protein = np.asarray(adata.obsm["protein"])
        names = list(adata.uns["protein_names"])
        for i, name in enumerate(names):
            frame[f"protein_{name}"] = protein[:, i]
    return frame


def norman_feature_frame() -> pd.DataFrame:
    adata = ad.read_h5ad("data/norman_small_program.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_genes"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("program_")]:
        frame[col] = obs[col].astype(float).values
    return frame


def papalexi_predicted_delta(ko: str, feature: str, model: str = "ridge") -> float:
    pred = pd.read_csv("results/papalexi_multimodal_joint_predictions.csv")
    row = pred.loc[pred["ko_target"].astype(str) == ko]
    if row.empty:
        raise KeyError(f"KO not found in Papalexi predictions: {ko}")
    prefix = "ridge_pred" if model == "ridge" else "pls_pred"
    if feature.startswith("pathway_"):
        col = f"{prefix}_delta_{feature}"
    elif feature.startswith("protein_"):
        protein = feature.removeprefix("protein_")
        col = f"{prefix}_delta_protein_{protein}"
    else:
        raise ValueError(feature)
    return float(row[col].iloc[0])


def norman_predicted_delta(ko: str, feature: str) -> float:
    pred = pd.read_csv("results/norman_system_prior_predictions.csv")
    row = pred.loc[pred["ko_genes"].astype(str) == ko]
    if row.empty:
        raise KeyError(f"KO not found in Norman predictions: {ko}")
    col = f"system_pred_delta_{feature}"
    return float(row[col].iloc[0])


def build_virtual_distribution(
    frame: pd.DataFrame,
    ko_col: str,
    ko: str,
    feature: str,
    predicted_delta: float,
    max_cells: int = 500,
) -> pd.DataFrame:
    ctrl = frame.loc[control_mask(frame[ko_col]), feature].dropna().astype(float)
    actual = frame.loc[frame[ko_col].astype(str) == ko, feature].dropna().astype(float)
    if len(ctrl) > max_cells:
        ctrl = ctrl.sample(max_cells, random_state=1)
    if len(actual) > max_cells:
        actual = actual.sample(max_cells, random_state=2)
    virtual = ctrl + predicted_delta
    return pd.DataFrame(
        {
            "value": pd.concat([ctrl, actual, virtual], ignore_index=True),
            "state": ["control"] * len(ctrl) + ["true KO cells"] * len(actual) + ["virtual KO cells"] * len(virtual),
            "ko": ko,
            "feature": feature,
        }
    )


def plot_papalexi_single_cell_distributions() -> None:
    frame = papalexi_feature_frame()
    selected = [
        ("STAT1", "pathway_IFNG_JAK_STAT", "pls"),
        ("JAK2", "pathway_IFNG_JAK_STAT", "pls"),
        ("IFNGR2", "protein_PDL1", "ridge"),
        ("JAK2", "protein_PDL1", "ridge"),
    ]
    rows = []
    for ko, feature, model in selected:
        rows.append(
            build_virtual_distribution(
                frame,
                "ko_target",
                ko,
                feature,
                papalexi_predicted_delta(ko, feature, model=model),
            )
        )
    long = pd.concat(rows, ignore_index=True)
    long["feature_label"] = long["feature"].map(PAPALEXI_FEATURES)
    long["panel"] = long["ko"] + "\n" + long["feature_label"]

    g = sns.catplot(
        data=long,
        x="state",
        y="value",
        col="panel",
        col_wrap=2,
        kind="violin",
        inner="quartile",
        cut=0,
        height=4,
        aspect=1.15,
        palette={"control": "#BDBDBD", "true KO cells": "#2A9D8F", "virtual KO cells": "#E76F51"},
        sharey=False,
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=20)
        ax.set_xlabel("")
        ax.set_ylabel("Single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Papalexi: single-cell virtual KO distributions", y=1.03)
    g.savefig(FIG_DIR / "papalexi_single_cell_virtual_ko_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_papalexi_single_cell_pca() -> None:
    frame = papalexi_feature_frame()
    features = [feature for feature in PAPALEXI_FEATURES if feature in frame.columns]
    selected_kos = ["NT", "STAT1", "JAK2", "IFNGR2", "IRF1"]
    mask = frame["ko_target"].isin(selected_kos) | control_mask(frame["ko_target"])
    plot = frame.loc[mask, ["ko_target"] + features].dropna().copy()
    plot["state"] = np.where(control_mask(plot["ko_target"]), "control", plot["ko_target"])
    if len(plot) > 1200:
        plot = plot.sample(1200, random_state=3)
    coords = PCA(n_components=2, random_state=4).fit_transform(plot[features].to_numpy(dtype=float))
    plot["PC1"] = coords[:, 0]
    plot["PC2"] = coords[:, 1]
    plt.figure(figsize=(8, 6))
    ax = sns.scatterplot(data=plot, x="PC1", y="PC2", hue="state", s=30, alpha=0.75)
    ax.set_title("Papalexi: cells in pathway/protein state space")
    ax.set_xlabel("PC1 from pathway/protein scores")
    ax.set_ylabel("PC2 from pathway/protein scores")
    ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_single_cell_state_space_pca.png", bbox_inches="tight")
    plt.close()


def plot_norman_single_cell_combo_distributions() -> None:
    frame = norman_feature_frame()
    selected = [
        ("AHR+KLF1", "program_ERYTHROID"),
        ("CEBPB+CEBPA", "program_GRANULOCYTE_APOPTOSIS"),
        ("MAPK1+TGFBR2", "program_MAPK_TGFB"),
        ("CBL+UBASH3B", "program_ERYTHROID"),
    ]
    rows = []
    for ko, feature in selected:
        rows.append(
            build_virtual_distribution(
                frame,
                "ko_target",
                ko,
                feature,
                norman_predicted_delta(ko, feature),
            )
        )
    long = pd.concat(rows, ignore_index=True)
    long["feature_label"] = long["feature"].map(NORMAN_FEATURES)
    long["panel"] = long["ko"] + "\n" + long["feature_label"]

    g = sns.catplot(
        data=long,
        x="state",
        y="value",
        col="panel",
        col_wrap=2,
        kind="violin",
        inner="quartile",
        cut=0,
        height=4,
        aspect=1.15,
        palette={"control": "#BDBDBD", "true KO cells": "#2A9D8F", "virtual KO cells": "#E76F51"},
        sharey=False,
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=20)
        ax.set_xlabel("")
        ax.set_ylabel("Single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Norman: single-cell distributions for double-gene virtual KO", y=1.03)
    g.savefig(FIG_DIR / "norman_single_cell_virtual_combo_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_single_cell_doc() -> None:
    text = """# 单细胞层面的虚拟敲除效果

之前的多数结果是 perturbation-level，也就是把同一个 KO 下的细胞聚合后看平均变化。为了展示单细胞效果，现在补充了单细胞分布图。

## 图怎么看

文件：

- `results/figures/papalexi_single_cell_virtual_ko_distributions.png`
- `results/figures/papalexi_single_cell_state_space_pca.png`
- `results/figures/norman_single_cell_virtual_combo_distributions.png`

每个小图比较三类细胞状态：

- `control`：真实 control 细胞的单细胞分布。
- `true KO cells`：真实 KO 细胞的单细胞分布。
- `virtual KO cells`：把模型预测的 KO delta 加到 control 细胞上，得到的虚拟 KO 单细胞分布。

如果 `virtual KO cells` 和 `true KO cells` 的分布接近，说明模型在单细胞状态空间里预测得较好。

## 当前结论

Papalexi 中，STAT1/JAK2/IFNGR2 相关 KO 对 IFNG-JAK-STAT pathway 和 PDL1 protein 的单细胞分布有明显移动，虚拟 KO 能捕捉主要方向，但幅度有时偏低。

Norman 中，部分双基因组合如 `AHR+KLF1`、`CEBPB+CEBPA` 的主要程序变化方向能被捕捉；但 `MAPK1+TGFBR2` 等组合仍预测不足。

## 重要限制

当前 virtual KO cells 是一个简化版本：

```text
virtual single-cell score = control single-cell score + predicted perturbation delta
```

也就是说，它展示的是“模型预测的 KO 后单细胞状态分布”，但还不是最终的 cell-level 条件生成模型。下一步需要训练真正的 cell-state-conditioned model，让每个 control 细胞都有自己的条件化 KO 响应。
"""
    Path("docs/single_cell_virtual_ko_effects.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup()
    plot_papalexi_single_cell_distributions()
    plot_papalexi_single_cell_pca()
    plot_norman_single_cell_combo_distributions()
    write_single_cell_doc()
    print("Saved single-cell virtual KO effect plots.")


if __name__ == "__main__":
    main()
