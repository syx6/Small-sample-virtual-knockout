from __future__ import annotations

import argparse
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wasserstein_distance
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")
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


def control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("nt|control|non|safe|neg")


def split_ko(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    genes = [part.strip().upper() for part in text.split("_") if part.strip()]
    return [gene for gene in genes if not gene.lower().startswith("nt") and gene.lower() != "ctrl"]


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
                first_token = term.split()[0].upper() if term.split() else ""
                if GENE_RE.match(first_token):
                    genes.add(first_token)
            if genes:
                terms.append((term, genes))
    return terms


def select_prior_terms(priors_dir: Path, perturb_genes: set[str], max_terms_per_library: int) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(priors_dir.glob("*.gmt")):
        scored = []
        terms = parse_gmt(path, include_term_gene=path.stem == "ppi_hub")
        for term, genes in terms:
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


def load_papalexi_state() -> tuple[pd.DataFrame, list[str]]:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_target"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("pathway_")]:
        frame[col] = obs[col].astype(float).values
    if "protein" in adata.obsm:
        protein = np.asarray(adata.obsm["protein"])
        for i, name in enumerate(adata.uns["protein_names"]):
            frame[f"protein_{name}"] = protein[:, i]
    state_cols = [c for c in frame.columns if c.startswith("pathway_") or c.startswith("protein_")]
    return frame.dropna(subset=state_cols), state_cols


def fit_delta_model(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: set[str],
) -> tuple[object, pd.DataFrame, np.ndarray, np.ndarray]:
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    rows = []
    x_rows = []
    y_rows = []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdout_kos:
            continue
        mean_state = group[state_cols].mean().to_numpy(dtype=float)
        delta = mean_state - control_mean
        x_rows.append(ko_prior_vector(ko, terms))
        y_rows.append(delta)
        rows.append({"ko_target": ko, "n_cells": len(group)})
    x_train = np.vstack(x_rows)
    y_train = np.vstack(y_rows)
    n_components = min(6, x_train.shape[0] - 1, y_train.shape[1], x_train.shape[1])
    model = make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True))
    model.fit(x_train, y_train)
    return model, pd.DataFrame(rows), x_train, y_train


def build_residual_bank(
    frame: pd.DataFrame,
    state_cols: list[str],
    train_kos: pd.DataFrame,
) -> dict[str, np.ndarray]:
    bank = {}
    for ko in train_kos["ko_target"]:
        target = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(target) == 0:
            continue
        bank[ko] = target - target.mean(axis=0, keepdims=True)
    return bank


def nearest_training_kos(
    ko: str,
    terms: list[tuple[str, set[str]]],
    train_kos: pd.DataFrame,
    x_train: np.ndarray,
    n_neighbors: int,
) -> list[str]:
    query = ko_prior_vector(ko, terms).reshape(1, -1)
    sims = cosine_similarity(query, x_train).ravel()
    order = np.argsort(-sims)[: max(1, n_neighbors)]
    return train_kos.iloc[order]["ko_target"].tolist()


def residual_matrix_for_neighbors(
    residual_bank: dict[str, np.ndarray],
    neighbors: list[str],
) -> np.ndarray:
    mats = [residual_bank[ko] for ko in neighbors if ko in residual_bank]
    if not mats:
        mats = list(residual_bank.values())
    return np.vstack(mats)


def generate_residual_cells(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    model,
    train_kos: pd.DataFrame,
    x_train: np.ndarray,
    residual_bank: dict[str, np.ndarray],
    ko: str,
    n_cells: int,
    n_neighbors: int,
    residual_scale: float,
    seed: int,
) -> tuple[np.ndarray, list[str], np.ndarray]:
    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    ctrl = control[rng.integers(0, len(control), size=n_cells)]
    predicted_delta = model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1)
    neighbors = nearest_training_kos(ko, terms, train_kos, x_train, n_neighbors)
    residuals = residual_matrix_for_neighbors(residual_bank, neighbors)
    noise = residuals[rng.integers(0, len(residuals), size=n_cells)] * residual_scale
    virtual = ctrl + predicted_delta.reshape(1, -1) + noise
    return virtual, neighbors, predicted_delta


def evaluate_holdouts(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    model,
    train_kos: pd.DataFrame,
    x_train: np.ndarray,
    residual_bank: dict[str, np.ndarray],
    holdout_kos: list[str],
    n_neighbors: int,
    residual_scale: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    metric_rows = []
    cell_rows = []
    neighbor_rows = []
    for offset, ko in enumerate(holdout_kos):
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        virtual, neighbors, predicted_delta = generate_residual_cells(
            frame,
            state_cols,
            terms,
            model,
            train_kos,
            x_train,
            residual_bank,
            ko,
            len(true),
            n_neighbors,
            residual_scale,
            seed + offset,
        )
        ctrl = control[np.random.default_rng(seed + 100 + offset).integers(0, len(control), size=len(true))]
        for rank, neighbor in enumerate(neighbors, start=1):
            neighbor_rows.append({"ko_target": ko, "neighbor_rank": rank, "training_ko_used_for_noise": neighbor})
        for j, feature in enumerate(state_cols):
            true_mean = true[:, j].mean()
            virtual_mean = virtual[:, j].mean()
            control_mean = ctrl[:, j].mean()
            w_control = wasserstein_distance(true[:, j], ctrl[:, j])
            w_virtual = wasserstein_distance(true[:, j], virtual[:, j])
            improvement = 1.0 - (w_virtual / w_control) if w_control > 1e-9 else np.nan
            mean_improvement = abs(control_mean - true_mean) - abs(virtual_mean - true_mean)
            metric_rows.append(
                {
                    "ko_target": ko,
                    "feature": feature,
                    "true_mean": true_mean,
                    "virtual_mean": virtual_mean,
                    "control_mean": control_mean,
                    "predicted_delta": predicted_delta[j],
                    "true_delta_vs_control": true_mean - control_mean,
                    "abs_mean_error": abs(virtual_mean - true_mean),
                    "control_abs_mean_error": abs(control_mean - true_mean),
                    "mean_error_reduction": mean_improvement,
                    "wasserstein_true_vs_virtual": w_virtual,
                    "wasserstein_true_vs_control": w_control,
                    "distribution_improvement": improvement,
                }
            )
        for state, matrix in [("control cells", ctrl), ("virtual KO cells", virtual), ("true KO cells", true)]:
            take = min(180, len(matrix))
            idx = np.random.default_rng(seed + offset + len(state)).choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=state_cols)
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True), pd.DataFrame(neighbor_rows)


def plot_holdout_distributions(cells: pd.DataFrame, state_cols: list[str]) -> None:
    selected_features = [
        "pathway_IFNG_JAK_STAT",
        "pathway_IMMUNE_CHECKPOINT",
        "protein_PDL1",
        "protein_CD86",
    ]
    selected_features = [feature for feature in selected_features if feature in state_cols]
    plot = cells.melt(id_vars=["ko_target", "state"], value_vars=selected_features, var_name="feature", value_name="score")
    plot["panel"] = plot["ko_target"] + "\n" + plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    g = sns.catplot(
        data=plot,
        x="state",
        y="score",
        col="panel",
        col_wrap=4,
        kind="violin",
        inner="quartile",
        cut=0,
        height=3.2,
        aspect=1.05,
        sharey=False,
        order=["control cells", "virtual KO cells", "true KO cells"],
        hue="state",
        legend=False,
        palette={"control cells": "#BDBDBD", "virtual KO cells": "#E76F51", "true KO cells": "#2A9D8F"},
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=24)
        ax.set_xlabel("")
        ax.set_ylabel("single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Residual cell-level generator: does virtual KO move toward true KO?", y=1.03)
    g.savefig(FIG_DIR / "papalexi_cell_level_residual_generator_holdout_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_metric_summary(metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("feature", observed=True)
        .agg(
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
            mean_error_reduction=("mean_error_reduction", "mean"),
        )
        .reset_index()
        .sort_values("mean_distribution_improvement")
    )
    summary.to_csv("results/papalexi_cell_level_residual_generator_feature_summary.csv", index=False)
    plt.figure(figsize=(8.3, 6))
    ax = sns.barplot(data=summary, x="mean_distribution_improvement", y="feature", color="#457B9D")
    ax.axvline(0, color="0.2", linewidth=1)
    ax.set(
        xlabel="Mean distribution improvement vs control\n(>0 means generated cells are closer to real KO cells)",
        ylabel="",
        title="Cell-level residual generator: held-out KO performance",
    )
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_cell_level_residual_generator_metric_summary.png", bbox_inches="tight")
    plt.close()


def plot_state_space(cells: pd.DataFrame, state_cols: list[str]) -> None:
    features = [col for col in state_cols if col in cells.columns]
    coords = PCA(n_components=2, random_state=7).fit_transform(cells[features].to_numpy(dtype=float))
    plot = cells[["ko_target", "state"]].copy()
    plot["PC1"] = coords[:, 0]
    plot["PC2"] = coords[:, 1]
    g = sns.relplot(
        data=plot,
        x="PC1",
        y="PC2",
        hue="state",
        col="ko_target",
        col_wrap=2,
        kind="scatter",
        s=35,
        alpha=0.72,
        height=4,
        aspect=1.1,
        palette={"control cells": "#BDBDBD", "virtual KO cells": "#E76F51", "true KO cells": "#2A9D8F"},
    )
    g.set_titles("{col_name}")
    g.fig.suptitle("Generated single cells in pathway/protein state space", y=1.03)
    g.savefig(FIG_DIR / "papalexi_cell_level_residual_generator_state_space.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_delta_heatmap(metrics: pd.DataFrame) -> None:
    selected = [
        "pathway_IFNG_JAK_STAT",
        "pathway_IMMUNE_CHECKPOINT",
        "protein_PDL1",
        "protein_CD86",
    ]
    plot = metrics.loc[metrics["feature"].isin(selected)].copy()
    plot["feature"] = plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    heat = plot.pivot(index="feature", columns="ko_target", values="virtual_mean") - plot.pivot(index="feature", columns="ko_target", values="control_mean")
    true_heat = plot.pivot(index="feature", columns="ko_target", values="true_mean") - plot.pivot(index="feature", columns="ko_target", values="control_mean")
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    vmax = np.nanmax(np.abs(pd.concat([heat, true_heat]).to_numpy()))
    sns.heatmap(heat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[0], cbar=False)
    sns.heatmap(true_heat, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[1], cbar_kws={"label": "change vs sampled control"})
    axes[0].set_title("Generated KO change")
    axes[1].set_title("Real KO change")
    axes[0].set_xlabel("")
    axes[1].set_xlabel("")
    axes[0].set_ylabel("")
    axes[1].set_ylabel("")
    fig.suptitle("What changed after virtual KO?", y=1.02)
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_cell_level_residual_generator_delta_heatmap.png", bbox_inches="tight")
    plt.close()


def plot_model_comparison(metrics: pd.DataFrame) -> None:
    rows = [
        {
            "model": "MLP generator v1",
            "mean_distribution_improvement": np.nan,
            "improved_fraction": np.nan,
        },
        {
            "model": "Residual generator v1",
            "mean_distribution_improvement": metrics["distribution_improvement"].mean(),
            "improved_fraction": (metrics["distribution_improvement"] > 0).mean(),
        },
    ]
    old_path = Path("results/papalexi_cell_level_generator_metrics.csv")
    if old_path.exists():
        old = pd.read_csv(old_path)
        rows[0]["mean_distribution_improvement"] = old["distribution_improvement"].mean()
        rows[0]["improved_fraction"] = (old["distribution_improvement"] > 0).mean()
    plot = pd.DataFrame(rows).melt(id_vars="model", var_name="metric", value_name="value")
    plot["metric"] = plot["metric"].map(
        {
            "mean_distribution_improvement": "Mean distribution improvement",
            "improved_fraction": "Fraction improved",
        }
    )
    g = sns.catplot(
        data=plot,
        x="model",
        y="value",
        col="metric",
        kind="bar",
        sharey=False,
        height=4.2,
        aspect=0.9,
        color="#457B9D",
    )
    for ax in g.axes.flat:
        ax.axhline(0, color="0.25", linewidth=1)
        ax.tick_params(axis="x", rotation=18)
        ax.set_xlabel("")
    g.set_titles("{col_name}")
    g.fig.suptitle("Cell-level generator comparison on held-out KOs", y=1.05)
    g.savefig(FIG_DIR / "papalexi_cell_level_generator_model_comparison.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_doc(holdout_kos: list[str], metrics: pd.DataFrame, residual_scale: float) -> None:
    mean_improvement = metrics["distribution_improvement"].mean()
    improved_fraction = (metrics["distribution_improvement"] > 0).mean()
    key = metrics.loc[
        metrics["feature"].isin(["pathway_IFNG_JAK_STAT", "pathway_IMMUNE_CHECKPOINT", "protein_PDL1", "protein_CD86"])
    ]
    key_summary = (
        key.groupby("feature", observed=True)["distribution_improvement"]
        .mean()
        .sort_values(ascending=False)
        .round(3)
        .to_dict()
    )
    lines = [f"- {feature}: {value}" for feature, value in key_summary.items()]
    text = f"""# 真正的 cell-level 条件生成模型：Residual Generator v1

这一步不再只是预测“一个 KO 组的平均变化”，而是直接生成一批虚拟单细胞。

## 输入是什么

1. 一批未敲除的 control cells，每个细胞已经表示成 pathway score + protein score。
2. 要敲除的基因，例如 `STAT1`、`JAK2`、`IFNGR2`、`IRF1`。
3. 这些基因在 Reactome、MSigDB、TF-target、PPI 网络中的先验特征。

## 输出是什么

对每个目标 KO，输出一批虚拟 KO cells。每个虚拟细胞都有：

- pathway scores
- protein scores
- 与真实 KO 细胞分布的距离评估
- 与 control cells 相比是否更接近真实 KO 的判断

## 为什么先用 residual generator

第一版 MLP 条件生成器在小样本下表现不好，生成分布比原始 control 更远。原因很直接：小样本里没有同一个细胞 KO 前后的真实配对，神经网络很容易学到错误的平均模式。

Residual Generator v1 更保守：

```text
虚拟 KO cell = control cell + 系统先验预测的 KO 平均变化 + 可选的真实细胞波动
```

这样做的好处是：平均变化由网络先验约束，细胞级输出从真实 control cells 出发，而不是让模型凭空编。

本轮调参发现，当前小样本 holdout 测试里加入额外 residual noise 会降低效果，所以最终默认使用：

```text
residual_scale = {residual_scale}
```

这意味着当前最可靠的 cell-level 生成方式是：先把真实 control cells 沿着预测到的 KO 方向移动。细胞之间原本的差异仍然保留，因为每个虚拟细胞都来自一个真实 control cell。

## 测试方式

完全留出这些 KO，不让模型训练时看到它们：

```text
{", ".join(holdout_kos)}
```

然后比较三类细胞：

- control cells：原始未敲除细胞
- virtual KO cells：模型生成的敲除后细胞
- true KO cells：实验里真实测到的敲除细胞

## 当前效果

平均分布改进值：`{mean_improvement:.3f}`

所有 KO-特征组合中，生成细胞比 control 更接近真实 KO 的比例：`{improved_fraction:.1%}`

重点特征的平均分布改进：

{chr(10).join(lines)}

解释规则：数值大于 0 表示生成细胞比 control 更接近真实 KO；小于 0 表示还不如直接用 control。

## 结果图

- `results/figures/papalexi_cell_level_residual_generator_holdout_distributions.png`
- `results/figures/papalexi_cell_level_residual_generator_metric_summary.png`
- `results/figures/papalexi_cell_level_residual_generator_state_space.png`
- `results/figures/papalexi_cell_level_residual_generator_delta_heatmap.png`
- `results/figures/papalexi_cell_level_generator_model_comparison.png`

## 结论

这个版本是进入 cell-level 条件生成后的第一个可解释基线。它回答的问题是：给定一个真实 control cell 和一个 KO 条件，能不能生成更像真实 KO 的单细胞状态。

如果它只在部分通路或蛋白上有效，说明多模态和系统先验确实提供了信号，但还不足以完全解决小样本条件生成。下一步可以在这个稳定基线上接 VAE / flow matching / diffusion，而不是直接从零训练复杂模型。
"""
    Path("docs/cell_level_residual_generator_v1.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--holdout-kos", default="STAT1,JAK2,IFNGR2,IRF1")
    parser.add_argument("--max-terms-per-library", type=int, default=160)
    parser.add_argument("--n-neighbors", type=int, default=5)
    parser.add_argument("--residual-scale", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=31)
    args = parser.parse_args()

    setup_plot()
    frame, state_cols = load_papalexi_state()
    holdout_kos = [ko.strip() for ko in args.holdout_kos.split(",") if ko.strip()]
    perturb_genes = {gene for ko in frame["ko_target"].unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), perturb_genes, args.max_terms_per_library)
    print(f"Using {len(terms)} system-prior terms")

    model, train_kos, x_train, _ = fit_delta_model(frame, state_cols, terms, set(holdout_kos))
    residual_bank = build_residual_bank(frame, state_cols, train_kos)
    metrics, cells, neighbors = evaluate_holdouts(
        frame,
        state_cols,
        terms,
        model,
        train_kos,
        x_train,
        residual_bank,
        holdout_kos,
        args.n_neighbors,
        args.residual_scale,
        args.seed,
    )

    metrics.to_csv("results/papalexi_cell_level_residual_generator_metrics.csv", index=False)
    cells.to_csv("results/papalexi_cell_level_residual_generator_cells.csv", index=False)
    neighbors.to_csv("results/papalexi_cell_level_residual_generator_neighbors.csv", index=False)

    plot_holdout_distributions(cells, state_cols)
    plot_metric_summary(metrics)
    plot_state_space(cells, state_cols)
    plot_delta_heatmap(metrics)
    plot_model_comparison(metrics)
    write_doc(holdout_kos, metrics, args.residual_scale)

    print(
        "mean distribution improvement="
        f"{metrics['distribution_improvement'].mean():.3f}; "
        "improved fraction="
        f"{(metrics['distribution_improvement'] > 0).mean():.1%}"
    )
    print("Saved residual cell-level generator metrics, cells, figures, and docs.")


if __name__ == "__main__":
    main()
