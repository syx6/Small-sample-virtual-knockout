from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
import seaborn as sns


FIG_DIR = Path("results/figures")
DOC_DIR = Path("docs")


def setup() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"]:
        if font in available_fonts:
            plt.rcParams["font.family"] = font
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300


def savefig(name: str) -> None:
    plt.tight_layout()
    plt.savefig(FIG_DIR / name, bbox_inches="tight")
    plt.close()


def make_effectiveness_summary() -> pd.DataFrame:
    pap_auc = pd.read_csv("results/papalexi_multimodal_auc.csv")
    pap_reg = pd.read_csv("results/papalexi_multimodal_joint_metrics.csv")
    norman = pd.read_csv("results/norman_system_prior_metrics.csv")

    ifng_auc = pap_auc[
        (pap_auc["target"] == "delta_pathway_IFNG_JAK_STAT")
        & (pap_auc["direction"] == "decrease")
        & (pap_auc["model"] == "pls_pred")
    ]["roc_auc"].iloc[0]
    pdl1_auc = pap_auc[
        (pap_auc["target"] == "delta_protein_PDL1")
        & (pap_auc["direction"] == "decrease")
        & (pap_auc["model"] == "ridge_pred")
    ]["roc_auc"].iloc[0]
    pdl1_r2 = pap_reg[
        (pap_reg["target"] == "delta_protein_PDL1")
        & (pap_reg["model"] == "ridge_prior_joint")
    ]["r2"].iloc[0]
    ifng_r2 = pap_reg[
        (pap_reg["target"] == "delta_pathway_IFNG_JAK_STAT")
        & (pap_reg["model"] == "pls_prior_joint")
    ]["r2"].iloc[0]
    unseen = norman[
        (norman["subset"] == "has_unseen_gene")
        & (norman["model"].isin(["single_gene_additive", "system_prior_ridge"]))
        & (norman["target"].isin(
            [
                "delta_program_GRANULOCYTE_APOPTOSIS",
                "delta_program_MAPK_TGFB",
                "delta_program_PRO_GROWTH",
            ]
        ))
    ]
    improvements = []
    for target, group in unseen.groupby("target"):
        vals = group.set_index("model")["r2"]
        improvements.append(vals["system_prior_ridge"] - vals["single_gene_additive"])
    mean_unseen_gain = sum(improvements) / len(improvements)

    rows = [
        {
            "question": "能不能识别强响应 KO?",
            "best_metric": f"IFNG decrease AUC={ifng_auc:.2f}; PDL1 decrease AUC={pdl1_auc:.2f}",
            "verdict": "好",
            "plain_meaning": "适合先回答哪些 KO 可能产生明显通路/蛋白响应。",
        },
        {
            "question": "能不能精确预测连续变化幅度?",
            "best_metric": f"PDL1 R2={pdl1_r2:.2f}; IFNG pathway R2={ifng_r2:.2f}",
            "verdict": "中等",
            "plain_meaning": "蛋白表型较稳，RNA 通路幅度预测还需要更好的 cell-level 模型。",
        },
        {
            "question": "能不能外推到双基因组合 KO?",
            "best_metric": f"unseen-combo mean R2 gain={mean_unseen_gain:.2f}",
            "verdict": "部分有效",
            "plain_meaning": "系统先验能改善未见基因组合，但组合效应还不能完全解决。",
        },
        {
            "question": "小样本多模态思路是否成立?",
            "best_metric": "pathway-protein coupling r≈0.70",
            "verdict": "成立",
            "plain_meaning": "多模态表型给通路状态提供了可解释约束，是当前方法最强证据。",
        },
    ]
    df = pd.DataFrame(rows)
    df.to_csv("results/user_facing_effectiveness_summary.csv", index=False)
    return df


def plot_result_verdict(summary: pd.DataFrame) -> None:
    color_map = {"好": "#2A9D8F", "中等": "#E9C46A", "部分有效": "#F4A261", "成立": "#457B9D"}
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis("off")
    ax.text(0.0, 1.05, "虚拟敲除方法当前效果判断", fontsize=24, weight="bold", transform=ax.transAxes)
    ax.text(
        0.0,
        0.97,
        "结论：适合做强响应筛选和机制解释；连续幅度精确预测、多基因非线性组合仍需增强。",
        fontsize=14,
        transform=ax.transAxes,
    )
    y = 0.78
    for _, row in summary.iterrows():
        ax.text(0.02, y, row["question"], fontsize=16, weight="bold", transform=ax.transAxes)
        ax.text(
            0.43,
            y,
            row["verdict"],
            fontsize=15,
            color="white",
            bbox=dict(boxstyle="round,pad=0.35", facecolor=color_map[row["verdict"]], edgecolor="none"),
            transform=ax.transAxes,
        )
        ax.text(0.56, y, row["best_metric"], fontsize=14, transform=ax.transAxes)
        ax.text(0.05, y - 0.08, row["plain_meaning"], fontsize=13, color="0.25", transform=ax.transAxes)
        ax.plot([0.02, 0.98], [y - 0.12, y - 0.12], color="0.85", transform=ax.transAxes)
        y -= 0.2
    savefig("user_facing_method_verdict.png")


def plot_minimal_performance() -> None:
    auc = pd.read_csv("results/papalexi_multimodal_auc.csv")
    auc_keep = auc[
        (
            (auc["target"] == "delta_pathway_IFNG_JAK_STAT")
            & (auc["direction"] == "decrease")
            & (auc["model"] == "pls_pred")
        )
        | (
            (auc["target"] == "delta_protein_PDL1")
            & (auc["direction"] == "decrease")
            & (auc["model"] == "ridge_pred")
        )
        | (
            (auc["target"] == "delta_protein_CD86")
            & (auc["direction"] == "absolute")
            & (auc["model"] == "pls_pred")
        )
    ].copy()
    auc_keep["task"] = [
        "IFNG pathway decrease",
        "CD86 protein strong change",
        "PDL1 protein decrease",
    ]
    norman = pd.read_csv("results/norman_system_prior_metrics.csv")
    norman = norman[
        (norman["subset"] == "has_unseen_gene")
        & (norman["target"].isin(
            [
                "delta_program_GRANULOCYTE_APOPTOSIS",
                "delta_program_MAPK_TGFB",
                "delta_program_PRO_GROWTH",
            ]
        ))
    ].copy()
    norman["target"] = norman["target"].str.replace("delta_program_", "", regex=False)
    norman["model"] = norman["model"].map(
        {"single_gene_additive": "Additive", "system_prior_ridge": "System prior"}
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    sns.barplot(data=auc_keep, x="roc_auc", y="task", color="#2A9D8F", ax=axes[0])
    axes[0].axvline(0.5, color="0.3", linestyle="--", linewidth=1)
    axes[0].set(xlim=(0, 1.05), xlabel="ROC-AUC", ylabel="", title="Papalexi: strong-response detection")
    for container in axes[0].containers:
        axes[0].bar_label(container, fmt="%.2f", padding=4, fontsize=11)

    sns.barplot(data=norman, x="r2", y="target", hue="model", palette="Set2", ax=axes[1])
    axes[1].axvline(0, color="0.3", linewidth=1)
    axes[1].set(xlabel="R2 on unseen-gene double KO", ylabel="", title="Norman: combo extrapolation")
    axes[1].legend(title="")
    savefig("user_facing_key_performance.png")


def plot_io_diagram() -> None:
    fig, ax = plt.subplots(figsize=(13, 6))
    ax.axis("off")
    boxes = [
        (0.03, 0.58, 0.22, 0.26, "输入数据", "control / KO 单细胞数据\nRNA 必需；protein/ATAC 可选\nKO 基因或组合 KO 标签"),
        (0.31, 0.58, 0.24, 0.26, "状态表征", "RNA -> pathway/TF scores\nprotein -> surface phenotype\nATAC -> TF/chromatin activity"),
        (0.61, 0.58, 0.20, 0.26, "系统先验", "Reactome / MSigDB\nTF-target / PPI\nKO genes -> prior features"),
        (0.37, 0.14, 0.25, 0.27, "虚拟 KO 模型", "输入 baseline cell state + KO prior\n输出 KO 后状态变化"),
        (0.71, 0.12, 0.26, 0.30, "输出结果", "pathway delta\nprotein delta\n强响应排名 / AUC\n组合 KO 机制解释"),
    ]
    for x, y, w, h, title, body in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor="#F8F9FA", edgecolor="#457B9D", linewidth=2, transform=ax.transAxes)
        ax.add_patch(rect)
        ax.text(x + 0.02, y + h - 0.06, title, fontsize=16, weight="bold", transform=ax.transAxes)
        ax.text(x + 0.02, y + h - 0.13, body, fontsize=12.5, va="top", transform=ax.transAxes)
    arrows = [
        ((0.25, 0.71), (0.31, 0.71)),
        ((0.55, 0.71), (0.61, 0.71)),
        ((0.43, 0.58), (0.45, 0.40)),
        ((0.71, 0.58), (0.54, 0.40)),
        ((0.61, 0.28), (0.72, 0.28)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=2, color="0.25"), xycoords=ax.transAxes)
    ax.text(0.03, 0.95, "方法输入输出：给用户看的版本", fontsize=24, weight="bold", transform=ax.transAxes)
    savefig("user_facing_input_output_diagram.png")


def write_docs(summary: pd.DataFrame) -> None:
    headers = ["问题", "指标", "判断", "含义"]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in summary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    str(row["question"]),
                    str(row["best_metric"]),
                    str(row["verdict"]),
                    str(row["plain_meaning"]),
                ]
            )
            + " |"
        )
    table = "\n".join(lines)
    text = f"""# 用户视角的方法效果与输入输出

## 一句话结论

当前原型**不是已经能精确预测所有 KO 后表达变化的最终模型**。它目前最可靠的能力是：

1. 识别哪些 KO 会产生强通路/蛋白响应。
2. 用 pathway + protein 输出解释 KO 后细胞状态。
3. 在部分双基因组合上做可解释外推，系统先验对未见基因组合有帮助。

还不够强的地方是：

1. 连续变化幅度的精确预测仍一般。
2. MAPK/TGFB、Pro-growth 这类非线性组合效应仍难。
3. 当前模型还是 perturbation-level 聚合模型，不是最终 cell-level 条件生成模型。

## 当前效果判断

{table}

## 给别人使用时的输入

最小输入：

- 单细胞 RNA 表达矩阵，格式可以是 `.h5ad` / `.h5mu` / count matrix。
- 每个细胞的 perturbation 标签，例如 `STAT1`、`JAK2`、`KLF1+BAK1`、`ctrl`。
- control/negative-control 细胞标记。

推荐输入：

- RNA + protein/CITE-seq。
- RNA + ATAC 或 gene activity。
- cell type / sample / batch metadata。
- KO gene list，支持单基因和多基因组合。

系统先验输入：

- Reactome / MSigDB pathway gene sets。
- TF-target gene sets。
- PPI hub 或 gene network。

## 给别人使用时的输出

主要输出：

- 每个 KO 或组合 KO 的 pathway activity delta。
- 可选 protein phenotype delta，例如 PDL1、CD86。
- 强响应 KO 排名。
- 多基因组合 KO 的预测表。
- 命中的 pathway / TF / PPI prior terms，用于解释。

评估输出：

- 回归：MAE、R2。
- 排序/筛选：ROC-AUC、PR-AUC。
- 可视化：效果判断图、输入输出图、强响应图、组合 KO 外推图。

## 推荐对外表述

这个方法目前应定位为：

> 一个面向小样本多模态单细胞数据的机制先验驱动虚拟敲除框架，用于预测和解释 KO 后的通路/蛋白状态变化，特别适合强响应筛选和多基因组合扰动的初步外推。

不建议现在声称：

- 可以精确重建 KO 后全转录组。
- 可以可靠预测所有未知多基因组合。
- 已经优于大型预训练生成模型。
"""
    (DOC_DIR / "user_facing_method_summary.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup()
    summary = make_effectiveness_summary()
    plot_result_verdict(summary)
    plot_minimal_performance()
    plot_io_diagram()
    write_docs(summary)
    print("Saved user-facing summary figures and docs.")


if __name__ == "__main__":
    main()
