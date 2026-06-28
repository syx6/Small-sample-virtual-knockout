from __future__ import annotations

import itertools
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
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")
FIG_DIR = Path("results/figures")
HOLDOUT_KOS = ["STAT1", "JAK2", "IFNGR2", "IRF1"]


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


def papalexi_delta_table(frame: pd.DataFrame, state_cols: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    rows = []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        rows.append(
            {
                "ko_target": ko,
                "n_cells": len(group),
                **{
                    f"delta_{feature}": value
                    for feature, value in zip(state_cols, group[state_cols].mean().to_numpy(dtype=float) - control_mean)
                },
            }
        )
    return pd.DataFrame(rows), control_mean


def norman_delta_table() -> pd.DataFrame:
    table = pd.read_csv("results/norman_program_delta.csv")
    return table.rename(columns={"ko_genes": "ko_target"})


def fit_pls(x: np.ndarray, y: np.ndarray, max_components: int = 6) -> object:
    n_components = min(max_components, x.shape[0] - 1, x.shape[1], y.shape[1])
    return make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True)).fit(x, y)


def fit_cross_response_model(norman: pd.DataFrame, terms: list[tuple[str, set[str]]]) -> object:
    x = np.vstack([ko_prior_vector(ko, terms) for ko in norman["ko_target"]])
    y_cols = [col for col in norman.columns if col.startswith("delta_program_")]
    y = norman[y_cols].to_numpy(dtype=float)
    return fit_pls(x, y, max_components=5)


def append_cross_embedding(x: np.ndarray, cross_model: object) -> np.ndarray:
    embedding = cross_model.predict(x)
    return np.hstack([x, embedding])


def weighted_knn_predict(x_train: np.ndarray, y_train: np.ndarray, x_query: np.ndarray, k: int) -> np.ndarray:
    sims = cosine_similarity(x_query.reshape(1, -1), x_train).ravel()
    order = np.argsort(-sims)[: min(k, len(sims))]
    selected = sims[order]
    weights = np.maximum(selected, 0)
    if weights.sum() < 1e-9:
        weights = np.ones_like(weights)
    weights = weights / weights.sum()
    return weights @ y_train[order]


def candidate_predictions(
    train_x: np.ndarray,
    train_y: np.ndarray,
    query_x: np.ndarray,
    cross_model: object,
) -> dict[str, np.ndarray]:
    base = fit_pls(train_x, train_y)
    cross = fit_pls(append_cross_embedding(train_x, cross_model), train_y)
    return {
        "pls": base.predict(query_x.reshape(1, -1)).reshape(-1),
        "cross_pls": cross.predict(append_cross_embedding(query_x.reshape(1, -1), cross_model)).reshape(-1),
        "knn3": weighted_knn_predict(train_x, train_y, query_x, k=3),
        "knn7": weighted_knn_predict(train_x, train_y, query_x, k=7),
        "zero": np.zeros(train_y.shape[1]),
    }


def tune_ensemble_weights(x: np.ndarray, y: np.ndarray, cross_model: object) -> tuple[dict[str, float], pd.DataFrame]:
    names = ["pls", "cross_pls", "knn3", "knn7", "zero"]
    loo_rows = []
    pred_by_candidate = {name: [] for name in names}
    truth = []
    for i in range(len(x)):
        keep = np.arange(len(x)) != i
        preds = candidate_predictions(x[keep], y[keep], x[i], cross_model)
        for name in names:
            pred_by_candidate[name].append(preds[name])
        truth.append(y[i])
        loo_rows.append({"heldout_training_ko": i})
    pred_by_candidate = {name: np.vstack(values) for name, values in pred_by_candidate.items()}
    truth = np.vstack(truth)

    best_weights = None
    best_loss = np.inf
    grid = np.linspace(0, 1, 5)
    for raw in itertools.product(grid, repeat=len(names)):
        total = sum(raw)
        if total <= 0:
            continue
        weights = np.asarray(raw) / total
        pred = sum(weights[j] * pred_by_candidate[name] for j, name in enumerate(names))
        loss = np.mean(np.abs(pred - truth))
        if loss < best_loss:
            best_loss = loss
            best_weights = weights
    weights = {name: float(best_weights[j]) for j, name in enumerate(names)}
    rows = []
    for name in names:
        rows.append({"candidate": name, "loo_mae": float(np.mean(np.abs(pred_by_candidate[name] - truth))), "weight": weights[name]})
    pred = sum(weights[name] * pred_by_candidate[name] for name in names)
    rows.append({"candidate": "ensemble", "loo_mae": float(np.mean(np.abs(pred - truth))), "weight": 1.0})
    return weights, pd.DataFrame(rows)


def predict_delta(
    train_x: np.ndarray,
    train_y: np.ndarray,
    query_x: np.ndarray,
    cross_model: object,
    weights: dict[str, float],
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    preds = candidate_predictions(train_x, train_y, query_x, cross_model)
    ensemble = sum(weights[name] * preds[name] for name in weights)
    return ensemble, preds


def evaluate_holdouts(
    frame: pd.DataFrame,
    state_cols: list[str],
    delta_table: pd.DataFrame,
    terms: list[tuple[str, set[str]]],
    cross_model: object,
    weights: dict[str, float],
    seed: int = 71,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    train = delta_table.loc[~delta_table["ko_target"].isin(HOLDOUT_KOS)].copy()
    train_x = np.vstack([ko_prior_vector(ko, terms) for ko in train["ko_target"]])
    y_cols = [f"delta_{feature}" for feature in state_cols]
    train_y = train[y_cols].to_numpy(dtype=float)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)

    metric_rows, cell_rows, pred_rows = [], [], []
    for ko in HOLDOUT_KOS:
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        ctrl = control[rng.integers(0, len(control), size=len(true))]
        query_x = ko_prior_vector(ko, terms)
        ensemble_delta, candidate_delta = predict_delta(train_x, train_y, query_x, cross_model, weights)
        candidate_delta["constrained_ensemble"] = ensemble_delta

        for model, delta in candidate_delta.items():
            virtual = ctrl + delta.reshape(1, -1)
            pred_rows.append({"ko_target": ko, "model": model, **{f"pred_delta_{feature}": value for feature, value in zip(state_cols, delta)}})
            for j, feature in enumerate(state_cols):
                w_control = wasserstein_distance(true[:, j], ctrl[:, j])
                w_virtual = wasserstein_distance(true[:, j], virtual[:, j])
                metric_rows.append(
                    {
                        "model": model,
                        "ko_target": ko,
                        "feature": feature,
                        "true_mean": true[:, j].mean(),
                        "virtual_mean": virtual[:, j].mean(),
                        "control_mean": ctrl[:, j].mean(),
                        "wasserstein_true_vs_virtual": w_virtual,
                        "wasserstein_true_vs_control": w_control,
                        "distribution_improvement": 1.0 - w_virtual / w_control if w_control > 1e-9 else np.nan,
                    }
                )
            if model in ["pls", "cross_pls", "knn3", "constrained_ensemble"]:
                for state, matrix in [("control cells", ctrl), ("true KO cells", true), (f"{model} virtual cells", virtual)]:
                    take = min(130, len(matrix))
                    idx = rng.choice(len(matrix), size=take, replace=False)
                    tmp = pd.DataFrame(matrix[idx], columns=state_cols)
                    tmp["ko_target"] = ko
                    tmp["state"] = state
                    cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True), pd.DataFrame(pred_rows)


def plot_summary(metrics: pd.DataFrame) -> None:
    label_map = {
        "pls": "PLS",
        "cross_pls": "Cross-data PLS",
        "knn3": "KNN-3",
        "knn7": "KNN-7",
        "zero": "Zero shrink",
        "constrained_ensemble": "Constrained ensemble",
    }
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(mean_distribution_improvement=("distribution_improvement", "mean"), improved_fraction=("improved", "mean"))
        .reset_index()
        .sort_values("mean_distribution_improvement", ascending=False)
    )
    summary["label"] = summary["model"].map(label_map)
    summary.to_csv("results/papalexi_constrained_direction_prior_summary.csv", index=False)
    plot = summary.melt(id_vars=["model", "label"], var_name="metric", value_name="value")
    plot["metric"] = plot["metric"].map(
        {"mean_distribution_improvement": "Mean distribution improvement", "improved_fraction": "Fraction improved"}
    )
    order = summary["label"].tolist()
    g = sns.catplot(data=plot, x="label", y="value", col="metric", kind="bar", order=order, sharey=False, height=4.1, aspect=1.05)
    for ax in g.axes.flat:
        ax.axhline(0, color="0.25", linewidth=1)
        ax.tick_params(axis="x", rotation=25)
        ax.set_xlabel("")
    g.set_titles("{col_name}")
    g.fig.suptitle("Stronger KO direction priors on held-out single-gene KOs", y=1.04)
    g.savefig(FIG_DIR / "papalexi_constrained_direction_prior_summary.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_by_ko_heatmap(metrics: pd.DataFrame) -> None:
    label_map = {
        "pls": "PLS",
        "cross_pls": "Cross PLS",
        "knn3": "KNN-3",
        "knn7": "KNN-7",
        "zero": "Zero",
        "constrained_ensemble": "Ensemble",
    }
    table = metrics.groupby(["model", "ko_target"], observed=True)["distribution_improvement"].mean().unstack("ko_target")
    order = metrics.groupby("model", observed=True)["distribution_improvement"].mean().sort_values(ascending=False).index
    table = table.loc[order]
    table.index = [label_map.get(idx, idx) for idx in table.index]
    vmax = np.nanmax(np.abs(table.to_numpy()))
    plt.figure(figsize=(7.3, 4.4))
    ax = sns.heatmap(table, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", cbar_kws={"label": "mean distribution improvement"})
    ax.set_title("Which constrained prior works for which held-out KO?")
    ax.set_xlabel("held-out KO")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_constrained_direction_prior_by_ko_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close()


def plot_distributions(cells: pd.DataFrame) -> None:
    selected = ["pathway_IFNG_JAK_STAT", "protein_PDL1", "protein_CD86"]
    selected = [feature for feature in selected if feature in cells.columns]
    plot = cells.melt(id_vars=["ko_target", "state"], value_vars=selected, var_name="feature", value_name="score")
    plot["state"] = plot["state"].replace(
        {
            "pls virtual cells": "PLS",
            "cross_pls virtual cells": "Cross PLS",
            "knn3 virtual cells": "KNN-3",
            "constrained_ensemble virtual cells": "Ensemble",
        }
    )
    plot["panel"] = plot["ko_target"] + "\n" + plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    g = sns.catplot(data=plot, x="state", y="score", col="panel", col_wrap=3, kind="box", showfliers=False, height=3.25, aspect=1.05, sharey=False, color="#8DA0CB")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Cell distributions from stronger KO direction priors", y=1.02)
    g.savefig(FIG_DIR / "papalexi_constrained_direction_prior_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_doc(weights: dict[str, float], loo: pd.DataFrame, metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(mean_distribution_improvement=("distribution_improvement", "mean"), improved_fraction=("improved", "mean"))
        .reset_index()
        .sort_values("mean_distribution_improvement", ascending=False)
    )
    lines = [f"- {row.model}: 平均分布改进 {row.mean_distribution_improvement:.3f}, 改进比例 {row.improved_fraction:.1%}" for row in summary.itertuples(index=False)]
    weight_lines = [f"- {name}: {value:.2f}" for name, value in weights.items() if value > 1e-9]
    loo_lines = [f"- {row.candidate}: LOO MAE {row.loo_mae:.3f}, ensemble weight {row.weight:.2f}" for row in loo.itertuples(index=False)]
    text = f"""# 增强 KO 方向先验：constrained direction prior

这一步不再加深生成模型，而是先把“KO 后平均状态往哪个方向移动”预测得更稳。

## 做了什么

对每个 held-out KO，模型融合了五类方向预测：

- `pls`：只用 Papalexi 训练 KO 的系统先验 PLS。
- `cross_pls`：先用 Norman 数据预训练一个 perturbation response embedding，再用于 Papalexi。
- `knn3`：找系统先验最相似的 3 个训练 KO，做加权平均。
- `knn7`：找系统先验最相似的 7 个训练 KO，做加权平均。
- `zero`：收缩到 0，防止方向预测过猛。

融合权重不是手调的，而是在 Papalexi 训练 KO 上做 leave-one-KO-out 自动选择。

## 自动选择的权重

{chr(10).join(weight_lines)}

## 训练 KO 留一验证

{chr(10).join(loo_lines)}

## Held-out KO 测试结果

{chr(10).join(lines)}

## 图

- `results/figures/papalexi_constrained_direction_prior_summary.png`
- `results/figures/papalexi_constrained_direction_prior_by_ko_heatmap.png`
- `results/figures/papalexi_constrained_direction_prior_distributions.png`

## 当前结论

这个实验检验了三件事：更强 KO 方向先验、更多训练 KO 的留一调权、以及 Norman -> Papalexi 的跨数据预训练 embedding。

如果 constrained ensemble 超过普通 PLS，说明方向先验增强有效；如果没有超过，说明当前限制主要来自跨数据状态不一致或训练 KO 数量仍然不足。无论哪种结果，都比继续盲目加深 VAE/flow/diffusion 更有信息量。
"""
    Path("docs/constrained_direction_prior.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    frame, state_cols = load_papalexi_state()
    pap_delta, _ = papalexi_delta_table(frame, state_cols)
    norman = norman_delta_table()
    all_genes = {gene for ko in pd.concat([pap_delta["ko_target"], norman["ko_target"]]) for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), all_genes)
    cross_model = fit_cross_response_model(norman, terms)

    train = pap_delta.loc[~pap_delta["ko_target"].isin(HOLDOUT_KOS)].copy()
    x_train = np.vstack([ko_prior_vector(ko, terms) for ko in train["ko_target"]])
    y_cols = [f"delta_{feature}" for feature in state_cols]
    y_train = train[y_cols].to_numpy(dtype=float)
    weights, loo = tune_ensemble_weights(x_train, y_train, cross_model)
    loo.to_csv("results/papalexi_constrained_direction_prior_loo.csv", index=False)

    metrics, cells, predictions = evaluate_holdouts(frame, state_cols, pap_delta, terms, cross_model, weights)
    metrics.to_csv("results/papalexi_constrained_direction_prior_metrics.csv", index=False)
    cells.to_csv("results/papalexi_constrained_direction_prior_cells.csv", index=False)
    predictions.to_csv("results/papalexi_constrained_direction_prior_predictions.csv", index=False)

    plot_summary(metrics)
    plot_by_ko_heatmap(metrics)
    plot_distributions(cells)
    write_doc(weights, loo, metrics)

    summary = metrics.groupby("model", observed=True)["distribution_improvement"].mean().sort_values(ascending=False)
    print("Selected ensemble weights:", {k: round(v, 3) for k, v in weights.items() if v > 1e-9})
    print(summary.round(3).to_string())
    print("Saved constrained direction prior metrics, figures, and docs.")


if __name__ == "__main__":
    main()
