from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns


FIG_DIR = Path("results/figures")


PROGRAM_LABELS = {
    "delta_program_ERYTHROID": "Erythroid",
    "delta_program_GRANULOCYTE_APOPTOSIS": "Granulocyte/apoptosis",
    "delta_program_MAPK_TGFB": "MAPK/TGFB",
    "delta_program_PRO_GROWTH": "Pro-growth",
    "delta_program_PIONEER_TF": "Pioneer TF",
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


def classify_combo(label: str, seen_single_genes: set[str]) -> str:
    genes = [gene.strip().upper() for gene in label.split("+") if gene.strip()]
    unseen = [gene for gene in genes if gene not in seen_single_genes]
    if not unseen:
        return "seen combo: both genes have single-gene training data"
    return "unseen combo: at least one gene lacks single-gene training data"


def build_explanation_table() -> pd.DataFrame:
    pred = pd.read_csv("results/norman_system_prior_predictions.csv")
    delta = pd.read_csv("results/norman_program_delta.csv")
    hits = pd.read_csv("results/norman_system_prior_term_hits.csv")

    single = delta[delta["n_ko_genes"] == 1]
    seen_single_genes = {
        str(label).replace("+ctrl", "").replace("ctrl+", "").upper()
        for label in single["ko_genes"]
    }

    rows = []
    for _, row in pred.iterrows():
        ko = row["ko_genes"]
        combo_type = classify_combo(ko, seen_single_genes)
        program_errors = {}
        true_values = {}
        pred_values = {}
        for target, label in PROGRAM_LABELS.items():
            true = row[f"true_{target}"]
            system_pred = row[f"system_pred_{target}"]
            true_values[label] = true
            pred_values[label] = system_pred
            program_errors[label] = abs(system_pred - true)
        max_program = max(true_values, key=lambda key: abs(true_values[key]))
        worst_program = max(program_errors, key=program_errors.get)
        rows.append(
            {
                "ko_genes": ko,
                "combo_type": combo_type,
                "n_cells": row["n_cells"],
                "strongest_true_program": max_program,
                "strongest_true_delta": true_values[max_program],
                "worst_predicted_program": worst_program,
                "worst_abs_error": program_errors[worst_program],
                **{f"true_{k}": v for k, v in true_values.items()},
                **{f"pred_{k}": v for k, v in pred_values.items()},
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(hits, on="ko_genes", how="left")
    out.to_csv("results/norman_combo_explanation_table.csv", index=False)
    return out


def plot_top_combo_heatmap(explain: pd.DataFrame) -> None:
    top = explain.reindex(explain["strongest_true_delta"].abs().sort_values(ascending=False).index).head(16)
    true_cols = [f"true_{label}" for label in PROGRAM_LABELS.values()]
    pred_cols = [f"pred_{label}" for label in PROGRAM_LABELS.values()]
    true_mat = top.set_index("ko_genes")[true_cols]
    pred_mat = top.set_index("ko_genes")[pred_cols]
    true_mat.columns = [c.replace("true_", "") for c in true_mat.columns]
    pred_mat.columns = [c.replace("pred_", "") for c in pred_mat.columns]

    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharey=True)
    vmax = max(float(true_mat.abs().max().max()), float(pred_mat.abs().max().max()))
    sns.heatmap(true_mat, cmap="vlag", center=0, vmax=vmax, vmin=-vmax, ax=axes[0], cbar=False)
    sns.heatmap(pred_mat, cmap="vlag", center=0, vmax=vmax, vmin=-vmax, ax=axes[1], cbar=True)
    axes[0].set_title("真实通路/程序变化")
    axes[1].set_title("模型预测变化")
    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[0].set_ylabel("双基因 KO")
    axes[1].set_ylabel("")
    fig.suptitle("Norman: 具体双基因 KO 引起的程序变化", fontsize=18, y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "norman_combo_true_vs_pred_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def plot_representative_combos(explain: pd.DataFrame) -> None:
    requested = [
        "CEBPB+CEBPA",
        "AHR+KLF1",
        "MAPK1+TGFBR2",
        "CBL+UBASH3B",
        "CDKN1B+CDKN1A",
        "KLF1+BAK1",
    ]
    selected = explain[explain["ko_genes"].isin(requested)].copy()
    if selected.empty:
        selected = explain.reindex(explain["strongest_true_delta"].abs().sort_values(ascending=False).index).head(6)
    rows = []
    for _, row in selected.iterrows():
        for label in PROGRAM_LABELS.values():
            rows.append({"ko_genes": row["ko_genes"], "program": label, "type": "true", "delta": row[f"true_{label}"]})
            rows.append({"ko_genes": row["ko_genes"], "program": label, "type": "predicted", "delta": row[f"pred_{label}"]})
    long = pd.DataFrame(rows)
    g = sns.catplot(
        data=long,
        x="delta",
        y="program",
        hue="type",
        col="ko_genes",
        col_wrap=2,
        kind="bar",
        height=3.8,
        aspect=1.25,
        palette={"true": "#2A9D8F", "predicted": "#E76F51"},
        sharex=True,
    )
    for ax in g.axes.flat:
        ax.axvline(0, color="0.25", linewidth=1)
        ax.set_ylabel("")
    g.set_axis_labels("Program delta", "")
    g.set_titles("{col_name}")
    g.fig.suptitle("Norman: 具体 KO 组合的真实变化 vs 预测变化", y=1.03)
    g.savefig(FIG_DIR / "norman_representative_combo_bars.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_seen_unseen_doc() -> None:
    text = """# seen gene / unseen gene 是什么意思

在 Norman 2019 的多基因组合实验里，我们故意做了一个外推测试：

```text
训练集：单基因扰动
测试集：双基因组合扰动
```

因此：

- **seen gene**：这个基因在训练集中出现过单基因扰动。例如训练里有 `KLF1` 单基因扰动，那么在 `AHR+KLF1` 里，`KLF1` 是 seen gene。
- **unseen gene**：这个基因没有单基因训练样本，只在双基因组合里出现。例如训练里没有 `AHR` 单基因扰动，那么 `AHR+KLF1` 这个组合属于含 unseen gene 的组合。
- **seen combo**：双基因组合里的两个基因都在单基因训练集中出现过。
- **unseen combo**：双基因组合里至少一个基因没有单基因训练样本。

这不是说生物学上“没见过这个基因”，而是说模型训练时没有见过这个基因的单独扰动效果。

## 为什么要这样分

因为多基因虚拟敲除最难的不是“见过的基因怎么相加”，而是：

1. 没有某个基因的单基因 KO 数据时，能不能靠 pathway/TF/PPI 先验外推？
2. 两个基因组合时，效应是否仍然接近加和？
3. 哪些组合出现非线性，模型预测会失败？

## 具体结果怎么看

推荐看两个文件：

- `results/norman_combo_explanation_table.csv`
- `results/figures/norman_combo_true_vs_pred_heatmap.png`
- `results/figures/norman_representative_combo_bars.png`

这些文件会告诉你：

```text
敲了什么基因组合
真实引起了哪个通路/程序变化
模型预测了什么变化
预测和真实差多少
这个组合命中了哪些系统先验 term
```
"""
    Path("docs/seen_unseen_gene_explanation.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup()
    explain = build_explanation_table()
    plot_top_combo_heatmap(explain)
    plot_representative_combos(explain)
    write_seen_unseen_doc()
    print("Saved combo-level explanations.")


if __name__ == "__main__":
    main()
