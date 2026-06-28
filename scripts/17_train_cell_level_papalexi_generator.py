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
from sklearn.decomposition import PCA
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.neural_network import MLPRegressor
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
                tokens = term.split()
                if tokens and GENE_RE.match(tokens[0].upper()):
                    genes.add(tokens[0].upper())
            if genes:
                terms.append((term, genes))
    return terms


def select_prior_terms(priors_dir: Path, perturb_genes: set[str], max_terms_per_library: int) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(priors_dir.glob("*.gmt")):
        terms = parse_gmt(path, include_term_gene=path.stem == "ppi_hub")
        scored = []
        for term, genes in terms:
            overlap = len(genes & perturb_genes)
            if overlap == 0:
                continue
            if len(genes) < 5 or len(genes) > 800:
                continue
            scored.append(((overlap, -len(genes)), f"{path.stem}:{term}", genes))
        scored.sort(reverse=True, key=lambda x: x[0])
        selected.extend((name, genes) for _, name, genes in scored[:max_terms_per_library])
    return selected


def split_ko(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    genes = [part.strip().upper() for part in text.split("_") if part.strip()]
    return [gene for gene in genes if not gene.lower().startswith("nt") and gene.lower() != "ctrl"]


def ko_prior_vector(label: str, terms: list[tuple[str, set[str]]]) -> np.ndarray:
    genes = set(split_ko(label))
    denom = max(1, len(genes))
    values = []
    for _, members in terms:
        values.append(len(genes & members) / denom)
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


def make_training_pairs(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: set[str],
    pairs_per_ko: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    xs = []
    ys = []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdout_kos:
            continue
        target = group[state_cols].to_numpy(dtype=float)
        n = min(pairs_per_ko, max(len(target), 1) * 3)
        ctrl_idx = rng.integers(0, len(control), size=n)
        target_idx = rng.integers(0, len(target), size=n)
        prior = ko_prior_vector(ko, terms)
        xs.append(np.hstack([control[ctrl_idx], np.tile(prior, (n, 1))]))
        ys.append(target[target_idx])
    return np.vstack(xs), np.vstack(ys)


def generate_virtual_cells(
    model,
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    ko: str,
    n_cells: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    ctrl = control[rng.integers(0, len(control), size=n_cells)]
    prior = ko_prior_vector(ko, terms)
    x = np.hstack([ctrl, np.tile(prior, (n_cells, 1))])
    return model.predict(x)


def evaluate_holdouts(
    model,
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: list[str],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    metric_rows = []
    cell_rows = []
    for offset, ko in enumerate(holdout_kos):
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        virtual = generate_virtual_cells(model, frame, state_cols, terms, ko, len(true), seed + offset)
        ctrl = control[np.random.default_rng(seed + 100 + offset).integers(0, len(control), size=len(true))]
        for j, feature in enumerate(state_cols):
            true_mean = true[:, j].mean()
            virt_mean = virtual[:, j].mean()
            ctrl_mean = ctrl[:, j].mean()
            w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
            w_virtual = wasserstein_distance(true[:, j], virtual[:, j])
            improvement = 1.0 - (w_virtual / w_ctrl) if w_ctrl > 1e-9 else np.nan
            metric_rows.append(
                {
                    "ko_target": ko,
                    "feature": feature,
                    "true_mean": true_mean,
                    "virtual_mean": virt_mean,
                    "control_mean": ctrl_mean,
                    "abs_mean_error": abs(virt_mean - true_mean),
                    "control_abs_mean_error": abs(ctrl_mean - true_mean),
                    "wasserstein_true_vs_virtual": w_virtual,
                    "wasserstein_true_vs_control": w_ctrl,
                    "distribution_improvement": improvement,
                }
            )
        for state, matrix in [("true KO cells", true), ("virtual KO cells", virtual), ("control cells", ctrl)]:
            take = min(160, len(matrix))
            idx = np.random.default_rng(seed + offset + len(state)).choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=state_cols)
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True)


def plot_holdout_distributions(cells: pd.DataFrame, state_cols: list[str]) -> None:
    selected_features = [
        "pathway_IFNG_JAK_STAT",
        "pathway_IMMUNE_CHECKPOINT",
        "protein_PDL1",
        "protein_CD86",
    ]
    selected_features = [f for f in selected_features if f in state_cols]
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
        height=3.3,
        aspect=1.05,
        sharey=False,
        palette={"control cells": "#BDBDBD", "true KO cells": "#2A9D8F", "virtual KO cells": "#E76F51"},
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=25)
        ax.set_xlabel("")
        ax.set_ylabel("single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Cell-level conditional generator: held-out KO distributions", y=1.03)
    g.savefig(FIG_DIR / "papalexi_cell_level_generator_holdout_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_state_space(cells: pd.DataFrame, state_cols: list[str]) -> None:
    features = [c for c in state_cols if c in cells.columns]
    x = cells[features].to_numpy(dtype=float)
    coords = PCA(n_components=2, random_state=7).fit_transform(x)
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
        alpha=0.75,
        height=4,
        aspect=1.1,
        palette={"control cells": "#BDBDBD", "true KO cells": "#2A9D8F", "virtual KO cells": "#E76F51"},
    )
    g.set_titles("{col_name}")
    g.fig.suptitle("Held-out KO cells in generated pathway/protein state space", y=1.03)
    g.savefig(FIG_DIR / "papalexi_cell_level_generator_state_space.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_metric_summary(metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("feature", observed=True)
        .agg(
            mean_distribution_improvement=("distribution_improvement", "mean"),
            mean_abs_error=("abs_mean_error", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
        .sort_values("mean_distribution_improvement")
    )
    summary.to_csv("results/papalexi_cell_level_generator_feature_summary.csv", index=False)
    plt.figure(figsize=(8, 6))
    ax = sns.barplot(data=summary, x="mean_distribution_improvement", y="feature", color="#457B9D")
    ax.axvline(0, color="0.25", linewidth=1)
    ax.set(
        xlabel="Mean distribution improvement vs control\n(>0 means virtual KO is closer to true KO)",
        ylabel="",
        title="Cell-level generator: distribution-level improvement",
    )
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_cell_level_generator_metric_summary.png", bbox_inches="tight")
    plt.close()


def write_doc(holdout_kos: list[str]) -> None:
    text = f"""# Cell-level 条件生成模型第一版

这一版开始从 KO 组平均预测，升级到单细胞条件生成。

## 模型输入

```text
单个 control 细胞的 pathway/protein state
+ KO gene 的 Reactome/MSigDB/TF-target/PPI 先验特征
```

## 模型输出

```text
该细胞在 KO 条件下的 pathway/protein state
```

## 训练方式

Papalexi 数据中没有同一个细胞 KO 前后的真实配对，因此训练时使用随机配对：

```text
control cell state + KO condition -> sampled true KO cell state
```

这使模型学习条件分布，而不是单个细胞的精确一一对应。

## 测试方式

完全留出这些 KO，不参与训练：

```text
{", ".join(holdout_kos)}
```

然后用 control 细胞生成这些 KO 的虚拟单细胞状态，并与真实 KO 细胞分布比较。

## 输出图

- `results/figures/papalexi_cell_level_generator_holdout_distributions.png`
- `results/figures/papalexi_cell_level_generator_state_space.png`
- `results/figures/papalexi_cell_level_generator_metric_summary.png`

## 如何解释

如果 virtual KO cells 的分布比 control cells 更接近 true KO cells，说明模型学到了 KO 条件下的单细胞状态移动。

当前模型是第一版轻量 MLP 条件生成器，还不是最终的 VAE / flow matching / diffusion 模型。
"""
    Path("docs/cell_level_generator_v1.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--holdout-kos", default="STAT1,JAK2,IFNGR2,IRF1")
    parser.add_argument("--pairs-per-ko", type=int, default=360)
    parser.add_argument("--max-terms-per-library", type=int, default=160)
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args()

    setup_plot()
    frame, state_cols = load_papalexi_state()
    holdout_kos = [ko.strip() for ko in args.holdout_kos.split(",") if ko.strip()]
    perturb_genes = {gene for ko in frame["ko_target"].unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), perturb_genes, args.max_terms_per_library)
    print(f"Using {len(terms)} system-prior terms")

    x_train, y_train = make_training_pairs(
        frame,
        state_cols,
        terms,
        set(holdout_kos),
        args.pairs_per_ko,
        args.seed,
    )
    model = make_pipeline(
        StandardScaler(),
        MLPRegressor(
            hidden_layer_sizes=(96, 48),
            activation="relu",
            alpha=1e-3,
            learning_rate_init=1e-3,
            max_iter=1200,
            early_stopping=True,
            validation_fraction=0.15,
            random_state=args.seed,
        ),
    )
    model.fit(x_train, y_train)
    train_pred = model.predict(x_train)
    print(f"train MAE={mean_absolute_error(y_train, train_pred):.3f}; train R2={r2_score(y_train, train_pred):.3f}")

    metrics, cells = evaluate_holdouts(model, frame, state_cols, terms, holdout_kos, args.seed)
    metrics.to_csv("results/papalexi_cell_level_generator_metrics.csv", index=False)
    cells.to_csv("results/papalexi_cell_level_generator_cells.csv", index=False)

    plot_holdout_distributions(cells, state_cols)
    plot_state_space(cells, state_cols)
    plot_metric_summary(metrics)
    write_doc(holdout_kos)
    print("Saved cell-level generator metrics, cells, figures, and docs.")


if __name__ == "__main__":
    import argparse

    main()
