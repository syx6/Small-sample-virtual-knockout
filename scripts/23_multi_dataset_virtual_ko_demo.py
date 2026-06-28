from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.stats import wasserstein_distance
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import umap


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")
FIG_DIR = Path("results/figures")


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


def clean_gene(label: str) -> str:
    text = str(label).upper()
    text = re.sub(r"^P-SG", "", text)
    text = re.sub(r"^SG", "", text)
    text = re.sub(r"_[0-9]+$", "", text)
    text = re.sub(r"-[0-9]+$", "", text)
    text = text.replace("P_", "").replace("P-", "")
    return text


def split_ko(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    genes = [clean_gene(part.strip()) for part in text.split("_") if part.strip()]
    out = []
    for gene in genes:
        if gene.lower().startswith("nt") or gene.lower() in {"ctrl", "control", "nan", "none"}:
            continue
        if GENE_RE.match(gene) and not gene.isdigit():
            out.append(gene)
    return out


def parse_gmt(path: Path, include_term_gene: bool = False) -> list[tuple[str, set[str]]]:
    terms = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = {gene.upper() for gene in parts[2:] if GENE_RE.match(gene.upper())}
            if include_term_gene:
                token = term.split()[0].upper() if term.split() else ""
                if GENE_RE.match(token):
                    genes.add(token)
            if genes:
                terms.append((term, genes))
    return terms


def select_prior_terms(priors_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 180) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(priors_dir.glob("*.gmt")):
        scored = []
        for term, genes in parse_gmt(path, include_term_gene=path.stem == "ppi_hub"):
            overlap = len(genes & perturb_genes)
            if overlap == 0 or len(genes) < 5 or len(genes) > 800:
                continue
            scored.append(((overlap, -len(genes)), f"{path.stem}:{term}", genes))
        scored.sort(reverse=True, key=lambda item: item[0])
        selected.extend((name, genes) for _, name, genes in scored[:max_terms_per_library])
    return selected


def ko_prior_vector(label: str, terms: list[tuple[str, set[str]]]) -> np.ndarray:
    genes = set(split_ko(label))
    denom = max(1, len(genes))
    values = [len(genes & members) / denom for _, members in terms]
    values.append(float(len(genes)))
    return np.asarray(values, dtype=float)


def fit_pls(train_labels: list[str], train_delta: np.ndarray, terms: list[tuple[str, set[str]]]):
    x = np.vstack([ko_prior_vector(label, terms) for label in train_labels])
    n_components = min(6, x.shape[0] - 1, x.shape[1], train_delta.shape[1])
    return make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True)).fit(x, train_delta)


def evaluate_virtual_cells(
    dataset: str,
    modality: str,
    state_representation: str,
    frame: pd.DataFrame,
    state_cols: list[str],
    holdouts: list[str],
    predicted_delta_by_ko: dict[str, np.ndarray],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    metric_rows, cell_rows, pred_rows = [], [], []
    for ko in holdouts:
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0 or ko not in predicted_delta_by_ko:
            continue
        ctrl = control[rng.integers(0, len(control), size=len(true))]
        delta = predicted_delta_by_ko[ko]
        virtual = ctrl + delta.reshape(1, -1)
        pred_rows.append(
            {
                "dataset": dataset,
                "input_modality": modality,
                "state_representation": state_representation,
                "ko_target": ko,
                "n_true_cells": len(true),
                "n_virtual_cells": len(virtual),
                "output": "virtual KO single-cell state table",
            }
        )
        for j, feature in enumerate(state_cols):
            w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
            w_virt = wasserstein_distance(true[:, j], virtual[:, j])
            metric_rows.append(
                {
                    "dataset": dataset,
                    "input_modality": modality,
                    "state_representation": state_representation,
                    "ko_target": ko,
                    "feature": feature,
                    "wasserstein_true_vs_virtual": w_virt,
                    "wasserstein_true_vs_control": w_ctrl,
                    "distribution_improvement": 1.0 - w_virt / w_ctrl if w_ctrl > 1e-9 else np.nan,
                }
            )
        for state, matrix in [("control cells", ctrl), ("virtual KO cells", virtual), ("true KO cells", true)]:
            take = min(180, len(matrix))
            idx = rng.choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=state_cols)
            tmp["dataset"] = dataset
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True), pd.DataFrame(pred_rows)


def load_papalexi() -> tuple[pd.DataFrame, list[str], list[str], dict[str, np.ndarray], str, str]:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_target"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("pathway_")]:
        frame[col] = obs[col].astype(float).values
    protein = np.asarray(adata.obsm["protein"])
    for i, name in enumerate(adata.uns["protein_names"]):
        frame[f"protein_{name}"] = protein[:, i]
    state_cols = [c for c in frame.columns if c.startswith("pathway_") or c.startswith("protein_")]
    holdouts = ["STAT1", "JAK2", "IFNGR2", "IRF1"]
    genes = {gene for ko in frame["ko_target"].unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), genes)
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    train_labels, train_delta = [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdouts:
            continue
        train_labels.append(ko)
        train_delta.append(group[state_cols].mean().to_numpy(dtype=float) - control_mean)
    model = fit_pls(train_labels, np.vstack(train_delta), terms)
    pred = {ko: model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1) for ko in holdouts}
    return frame, state_cols, holdouts, pred, "RNA pathway + ADT protein", "pathway/protein scores"


def load_norman() -> tuple[pd.DataFrame, list[str], list[str], dict[str, np.ndarray], str, str]:
    adata = ad.read_h5ad("data/norman_small_program.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_genes"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("program_")]:
        frame[col] = obs[col].astype(float).values
    state_cols = [c for c in frame.columns if c.startswith("program_")]
    holdouts = ["AHR+KLF1", "CEBPB+CEBPA", "MAPK1+TGFBR2", "CBL+UBASH3B"]
    pred_table = pd.read_csv("results/norman_system_prior_predictions.csv")
    pred = {}
    for ko in holdouts:
        row = pred_table.loc[pred_table["ko_genes"].astype(str) == ko]
        if row.empty:
            continue
        pred[ko] = np.asarray([float(row[f"system_pred_delta_{feature}"].iloc[0]) for feature in state_cols])
    return frame, state_cols, holdouts, pred, "single-cell RNA perturb-seq", "gene program scores"


def matrix_subset_to_float(x):
    if sparse.issparse(x):
        return x.astype(np.float32)
    return np.asarray(x, dtype=np.float32)


def external_labels(adata: ad.AnnData, dataset_name: str) -> tuple[pd.Series, pd.Series]:
    obs = adata.obs.copy()
    if "DatlingerBock2021" in dataset_name:
        label = obs["perturbation"].astype(str).map(clean_gene)
        is_control = obs["perturbation"].astype(str).str.lower().eq("control")
    elif "DixitRegev2016" in dataset_name:
        label = obs["target"].astype(str).map(clean_gene)
        is_control = obs["perturbation"].astype(str).str.lower().eq("control")
        label = label.mask(label.str.lower().isin(["nan", "none"]), "control")
    else:
        label = obs["perturbation"].astype(str).map(clean_gene)
        is_control = control_mask(label)
    return label, is_control


def pathway_label(term: str) -> str:
    label = term.split(":", 1)[-1]
    label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
    return label[:46]


def select_state_terms(
    adata: ad.AnnData,
    perturbations: list[str],
    max_terms: int = 14,
) -> list[tuple[str, list[str]]]:
    var_genes = {str(gene).upper() for gene in adata.var_names}
    perturb_genes = {gene for ko in perturbations for gene in split_ko(ko)}
    candidates = []
    for path in [Path("data/priors/hallmark.gmt"), Path("data/priors/reactome.gmt"), Path("data/priors/tf_target.gmt")]:
        for term, genes in parse_gmt(path):
            present = sorted(genes & var_genes)
            if len(present) < 8 or len(present) > 350:
                continue
            overlap = len(genes & perturb_genes)
            immune_bonus = int(any(word in term.upper() for word in ["T_CELL", "INTERFERON", "CYTOKINE", "MAPK", "JAK", "STAT", "APOPTOSIS", "PROLIFERATION", "E2F", "MYC"]))
            if overlap == 0 and immune_bonus == 0:
                continue
            candidates.append(((overlap, immune_bonus, -len(present)), f"{path.stem}:{term}", present))
    candidates.sort(reverse=True, key=lambda item: item[0])
    selected = []
    seen = set()
    for _, name, genes in candidates:
        label = f"pathway_{pathway_label(name)}"
        if label in seen:
            continue
        selected.append((label, genes))
        seen.add(label)
        if len(selected) >= max_terms:
            break
    return selected


def compute_pathway_scores(adata: ad.AnnData, keep_idx: np.ndarray, terms: list[tuple[str, list[str]]]) -> pd.DataFrame:
    var_lookup = {str(gene).upper(): i for i, gene in enumerate(adata.var_names)}
    x = matrix_subset_to_float(adata.X[keep_idx])
    scores = {}
    for label, genes in terms:
        idx = [var_lookup[gene] for gene in genes if gene in var_lookup]
        if not idx:
            continue
        sub = x[:, idx]
        values = np.asarray(sub.mean(axis=1)).reshape(-1)
        scores[label] = values
    frame = pd.DataFrame(scores)
    for col in frame.columns:
        std = frame[col].std()
        if std > 1e-9:
            frame[col] = (frame[col] - frame[col].mean()) / std
        else:
            frame[col] = 0.0
    return frame


def load_external_pathway(path: Path, holdouts: list[str], seed: int) -> tuple[pd.DataFrame, list[str], list[str], dict[str, np.ndarray], str, str]:
    rng = np.random.default_rng(seed)
    adata = ad.read_h5ad(path)
    label, is_control = external_labels(adata, path.name)
    counts = label.loc[~is_control & ~label.str.upper().str.startswith("INTERGENIC")].value_counts()
    selected = counts.loc[counts >= 180].head(14).index.tolist()
    holdouts = [ko for ko in holdouts if ko in selected]
    selected = sorted(set(selected) | set(holdouts))
    ctrl_idx = np.flatnonzero(is_control.to_numpy())
    keep_idx = list(rng.choice(ctrl_idx, size=min(1400, len(ctrl_idx)), replace=False))
    for ko in selected:
        idx = np.flatnonzero((label == ko).to_numpy())
        if len(idx) > 0:
            keep_idx.extend(rng.choice(idx, size=min(420, len(idx)), replace=False).tolist())
    keep_idx = np.asarray(sorted(set(keep_idx)))
    state_terms = select_state_terms(adata, selected, max_terms=14)
    frame = compute_pathway_scores(adata, keep_idx, state_terms)
    frame["ko_target"] = label.iloc[keep_idx].to_numpy()
    frame.loc[is_control.iloc[keep_idx].to_numpy(), "ko_target"] = "control"
    state_cols = [c for c in frame.columns if c.startswith("pathway_")]
    genes = {gene for ko in selected for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), genes)
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    train_labels, train_delta = [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdouts:
            continue
        train_labels.append(ko)
        train_delta.append(group[state_cols].mean().to_numpy(dtype=float) - control_mean)
    model = fit_pls(train_labels, np.vstack(train_delta), terms)
    pred = {ko: model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1) for ko in holdouts}
    return frame, state_cols, holdouts, pred, "single-cell RNA CRISPR screen", "pathway/program scores"


def fit_umap(cells: pd.DataFrame, features: list[str], seed: int) -> pd.DataFrame:
    matrix = cells[features].to_numpy(dtype=float)
    matrix = StandardScaler().fit_transform(matrix)
    reducer = umap.UMAP(n_neighbors=min(30, max(5, len(cells) // 30)), min_dist=0.25, random_state=seed)
    coords = reducer.fit_transform(matrix)
    out = cells.copy()
    out["UMAP1"] = coords[:, 0]
    out["UMAP2"] = coords[:, 1]
    return out


def add_centroid_arrows(ax, plot: pd.DataFrame) -> None:
    centers = plot.groupby("state", observed=True)[["UMAP1", "UMAP2"]].mean()
    if not {"control cells", "virtual KO cells", "true KO cells"}.issubset(set(centers.index)):
        return
    start = centers.loc["control cells"]
    for state, color, linestyle in [
        ("virtual KO cells", PALETTE["virtual KO cells"], "-"),
        ("true KO cells", PALETTE["true KO cells"], "--"),
    ]:
        end = centers.loc[state]
        ax.scatter(end["UMAP1"], end["UMAP2"], s=90, marker="X", color=color, edgecolor="white", zorder=5)
        ax.annotate(
            "",
            xy=(end["UMAP1"], end["UMAP2"]),
            xytext=(start["UMAP1"], start["UMAP2"]),
            arrowprops={"arrowstyle": "->", "color": color, "lw": 2.8, "linestyle": linestyle, "mutation_scale": 14},
            zorder=6,
        )


def plot_multi_dataset_umap(cells: pd.DataFrame, example_kos: dict[str, str], feature_cols_by_dataset: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    metric_rows = []
    for dataset, ko in example_kos.items():
        sub = cells.loc[(cells["dataset"] == dataset) & (cells["ko_target"] == ko)].copy()
        if sub.empty:
            continue
        features = feature_cols_by_dataset[dataset]
        embedded = fit_umap(sub, features, seed=101 + len(rows))
        rows.append(embedded)
        centers = embedded.groupby("state", observed=True)[["UMAP1", "UMAP2"]].mean()
        control_to_true = float(np.linalg.norm(centers.loc["control cells"] - centers.loc["true KO cells"]))
        virtual_to_true = float(np.linalg.norm(centers.loc["virtual KO cells"] - centers.loc["true KO cells"]))
        metric_rows.append(
            {
                "dataset": dataset,
                "example_ko": ko,
                "control_to_true_umap_distance": control_to_true,
                "virtual_to_true_umap_distance": virtual_to_true,
                "umap_centroid_improvement": 1.0 - virtual_to_true / control_to_true if control_to_true > 1e-9 else np.nan,
            }
        )
    plot_cells = pd.concat(rows, ignore_index=True)
    metric = pd.DataFrame(metric_rows)

    datasets = list(example_kos.keys())
    fig, axes = plt.subplots(len(datasets), 2, figsize=(11.2, 3.0 * len(datasets)), sharex=False, sharey=False)
    for row, dataset in enumerate(datasets):
        ko = example_kos[dataset]
        group = plot_cells.loc[(plot_cells["dataset"] == dataset) & (plot_cells["ko_target"] == ko)]
        if group.empty:
            continue
        for col, mode in enumerate(["before", "after"]):
            ax = axes[row, col]
            if mode == "before":
                sub = group.loc[group["state"] == "control cells"]
                sns.scatterplot(data=sub, x="UMAP1", y="UMAP2", color=PALETTE["control cells"], s=24, alpha=0.78, linewidth=0, ax=ax)
                ax.set_title(f"{dataset}\n{ko}: before KO")
            else:
                for state, alpha, size in [("control cells", 0.18, 18), ("virtual KO cells", 0.72, 28), ("true KO cells", 0.72, 28)]:
                    sub = group.loc[group["state"] == state]
                    sns.scatterplot(data=sub, x="UMAP1", y="UMAP2", color=PALETTE[state], s=size, alpha=alpha, linewidth=0, label=state if row == 0 else None, ax=ax)
                add_centroid_arrows(ax, group)
                ax.set_title(f"{dataset}\n{ko}: virtual vs real")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
    handles, labels = axes[0, 1].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.0))
        axes[0, 1].legend_.remove()
    fig.suptitle("Virtual KO cell-state changes across modalities and datasets", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "multi_dataset_virtual_ko_umap_examples.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    plot_cells.to_csv("results/multi_dataset_virtual_ko_umap_cells.csv", index=False)
    return metric


def plot_effect_summary(summary: pd.DataFrame) -> None:
    plot = summary.sort_values("mean_distribution_improvement", ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sns.barplot(data=plot, x="dataset", y="mean_distribution_improvement", ax=axes[0], color="#4E79A7")
    sns.barplot(data=plot, x="dataset", y="improved_fraction", ax=axes[1], color="#4E79A7")
    axes[0].axhline(0, color="0.25", linewidth=1)
    axes[0].set_title("Mean distribution improvement")
    axes[1].set_title("Fraction improved")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=25)
    fig.suptitle("Virtual KO performance across input modalities", y=1.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "multi_dataset_virtual_ko_effect_summary.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_dataset_ko_heatmap(metrics: pd.DataFrame) -> None:
    plot = (
        metrics.groupby(["dataset", "ko_target"], observed=True)["distribution_improvement"]
        .mean()
        .reset_index()
    )
    short_dataset = {
        "Papalexi ECCITE-seq": "Papalexi",
        "Norman Perturb-seq": "Norman",
        "Datlinger CRISPR RNA": "Datlinger",
        "Dixit Perturb-seq RNA": "Dixit",
    }
    plot["dataset_short"] = plot["dataset"].map(short_dataset).fillna(plot["dataset"])
    plot["label"] = plot["ko_target"] + " (" + plot["dataset_short"] + ")"
    plot = plot.sort_values(["dataset", "distribution_improvement"], ascending=[True, False])
    plt.figure(figsize=(8.8, 6.8))
    colors = np.where(plot["distribution_improvement"] >= 0, "#4E79A7", "#C65D4B")
    bars = plt.barh(plot["label"], plot["distribution_improvement"], color=colors)
    plt.axvline(0, color="0.25", linewidth=1)
    plt.xlabel("Mean distribution improvement\n(>0 means virtual KO is closer to true KO)")
    plt.ylabel("")
    plt.title("Virtual KO effect by dataset and target")
    plt.tick_params(axis="y", labelsize=10)
    for bar, value in zip(bars, plot["distribution_improvement"]):
        x = bar.get_width()
        plt.text(
            x + (0.015 if x >= 0 else -0.015),
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}",
            ha="left" if x >= 0 else "right",
            va="center",
            fontsize=9,
        )
    plt.tight_layout()
    plt.savefig(FIG_DIR / "multi_dataset_virtual_ko_by_target.png", bbox_inches="tight", dpi=300)
    plt.close()


def write_doc(input_output: pd.DataFrame, summary: pd.DataFrame, umap_metrics: pd.DataFrame) -> None:
    text = f"""# 多数据集虚拟敲除实例：输入、输出和效果

这一步把同一套 hard-constrained residual/PLS 虚拟敲除方法放到多种数据上测试。

## 方法原则

```text
虚拟 KO cell = control cell + PLS/residual baseline 预测的 KO 方向
```

方向由系统先验和训练 KO 学到，生成模型不允许自由改变方向。对于小样本数据，当前更适合输出虚拟 KO 细胞状态和不确定性范围，而不是自由生成复杂分布。

## 实例输入和输出

输入表：

{input_output.to_string(index=False)}

输出包括：

- 每个目标 KO 的虚拟单细胞状态表
- 每个数据集的效果指标
- before/after UMAP 图
- UMAP 质心移动指标

## 效果摘要

{summary.round(3).to_string(index=False)}

## UMAP 质心移动

{umap_metrics.round(3).to_string(index=False)}

## 图

- `results/figures/multi_dataset_virtual_ko_effect_summary.png`
- `results/figures/multi_dataset_virtual_ko_by_target.png`
- `results/figures/multi_dataset_virtual_ko_umap_examples.png`

## 怎么解释

如果 mean distribution improvement 大于 0，说明虚拟 KO 细胞比原始 control cells 更接近真实 KO cells。

如果 UMAP centroid improvement 大于 0，说明在单细胞状态空间里，虚拟 KO 的细胞云团质心向真实 KO 云团移动。

当前结果用于展示方法可以接收多种输入模态：RNA+protein、多基因 RNA perturb-seq、普通 CRISPR RNA 数据。不同数据集效果不同，这正是需要用实例输出给别人看的原因。
"""
    Path("docs/multi_dataset_virtual_ko_demo.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    datasets = []
    datasets.append(("Papalexi ECCITE-seq", *load_papalexi()))
    datasets.append(("Norman Perturb-seq", *load_norman()))
    datasets.append(("Datlinger CRISPR RNA", *load_external_pathway(Path("data/scperturb_extra/DatlingerBock2021.h5ad"), ["LAT", "LCK", "JUND", "FOS"], seed=31)))
    datasets.append(("Dixit Perturb-seq RNA", *load_external_pathway(Path("data/scperturb_extra/DixitRegev2016_K562_TFs_13_days.h5ad"), ["ELF1", "CREB1", "ELK1", "GABPA"], seed=37)))

    all_metrics, all_cells, all_io = [], [], []
    feature_cols_by_dataset = {}
    for i, (dataset, frame, state_cols, holdouts, pred, modality, representation) in enumerate(datasets):
        metrics, cells, io = evaluate_virtual_cells(dataset, modality, representation, frame, state_cols, holdouts, pred, seed=131 + i)
        all_metrics.append(metrics)
        all_cells.append(cells)
        all_io.append(io)
        feature_cols_by_dataset[dataset] = state_cols

    metrics = pd.concat(all_metrics, ignore_index=True)
    cells = pd.concat(all_cells, ignore_index=True)
    input_output = pd.concat(all_io, ignore_index=True)
    metrics.to_csv("results/multi_dataset_virtual_ko_metrics.csv", index=False)
    cells.to_csv("results/multi_dataset_virtual_ko_cells.csv", index=False)
    input_output.to_csv("results/multi_dataset_virtual_ko_input_output_examples.csv", index=False)

    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby(["dataset", "input_modality", "state_representation"], observed=True)
        .agg(
            n_ko=("ko_target", "nunique"),
            n_features=("feature", "nunique"),
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
    )
    summary.to_csv("results/multi_dataset_virtual_ko_summary.csv", index=False)
    plot_effect_summary(summary)
    plot_dataset_ko_heatmap(metrics)

    example_kos = {
        "Papalexi ECCITE-seq": "STAT1",
        "Norman Perturb-seq": "CEBPB+CEBPA",
        "Datlinger CRISPR RNA": "LAT",
        "Dixit Perturb-seq RNA": "ELF1",
    }
    umap_metrics = plot_multi_dataset_umap(cells, example_kos, feature_cols_by_dataset)
    umap_metrics.to_csv("results/multi_dataset_virtual_ko_umap_metrics.csv", index=False)

    write_doc(input_output, summary, umap_metrics)
    print(summary.round(3).to_string(index=False))
    print("\nUMAP examples")
    print(umap_metrics.round(3).to_string(index=False))
    print("Saved multi-dataset virtual KO demo results.")


if __name__ == "__main__":
    main()
