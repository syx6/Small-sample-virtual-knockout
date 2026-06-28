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
from sklearn.decomposition import TruncatedSVD
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


def select_prior_terms(priors_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 160) -> list[tuple[str, set[str]]]:
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


def fit_pls_delta(frame: pd.DataFrame, state_cols: list[str], terms: list[tuple[str, set[str]]], holdouts: set[str]):
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    x_rows, y_rows, kos = [], [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdouts:
            continue
        x_rows.append(ko_prior_vector(ko, terms))
        y_rows.append(group[state_cols].mean().to_numpy(dtype=float) - control_mean)
        kos.append(ko)
    x = np.vstack(x_rows)
    y = np.vstack(y_rows)
    n_components = min(6, x.shape[0] - 1, x.shape[1], y.shape[1])
    model = make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True)).fit(x, y)
    return model, pd.DataFrame({"ko_target": kos}), x, y


def papalexi_residual_bank(frame: pd.DataFrame, state_cols: list[str], train_kos: pd.DataFrame) -> np.ndarray:
    residuals = []
    for ko in train_kos["ko_target"]:
        target = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(target) == 0:
            continue
        residuals.append(target - target.mean(axis=0, keepdims=True))
    return np.vstack(residuals)


def external_dataset_labels(adata: ad.AnnData, name: str) -> tuple[pd.Series, pd.Series]:
    obs = adata.obs.copy()
    if "DatlingerBock2021" in name:
        label = obs["perturbation"].astype(str)
        is_control = label.str.lower().eq("control")
    elif "DixitRegev2016" in name:
        label = obs["target"].astype(str)
        is_control = obs["perturbation"].astype(str).str.lower().eq("control")
        label = label.mask(label.str.lower().isin(["nan", "none"]), "control")
    else:
        label = obs["perturbation"].astype(str)
        is_control = control_mask(label)
    return label, is_control


def matrix_subset_to_float(x):
    if sparse.issparse(x):
        return x.astype(np.float32)
    return np.asarray(x, dtype=np.float32)


def estimate_external_uncertainty(seed: int = 83) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for path in sorted(Path("data/scperturb_extra").glob("*.h5ad")):
        if "DatlingerBock2021" not in path.name and "DixitRegev2016" not in path.name:
            continue
        adata = ad.read_h5ad(path)
        label, is_control = external_dataset_labels(adata, path.name)
        perturb_counts = label.loc[~is_control & ~label.str.upper().str.startswith("INTERGENIC")].value_counts()
        perturbations = perturb_counts.loc[perturb_counts >= 180].head(12).index.tolist()
        ctrl_idx = np.flatnonzero(is_control.to_numpy())
        keep_idx = list(rng.choice(ctrl_idx, size=min(1200, len(ctrl_idx)), replace=False))
        for pert in perturbations:
            idx = np.flatnonzero((label == pert).to_numpy())
            keep_idx.extend(rng.choice(idx, size=min(350, len(idx)), replace=False).tolist())
        keep_idx = np.asarray(sorted(set(keep_idx)))
        x = matrix_subset_to_float(adata.X[keep_idx])
        svd = TruncatedSVD(n_components=9, random_state=seed)
        emb = svd.fit_transform(x)
        kept_label = label.iloc[keep_idx].reset_index(drop=True)
        kept_control = is_control.iloc[keep_idx].reset_index(drop=True)
        control_emb = emb[kept_control.to_numpy()]
        control_mean = control_emb.mean(axis=0)
        for pert in perturbations:
            target = emb[(kept_label == pert).to_numpy()]
            if len(target) < 80:
                continue
            delta = target.mean(axis=0) - control_mean
            effect_norm = np.linalg.norm(delta)
            residual_radius = np.median(np.linalg.norm(target - target.mean(axis=0), axis=1))
            control_radius = np.median(np.linalg.norm(control_emb - control_mean, axis=1))
            rows.append(
                {
                    "dataset": path.name,
                    "perturbation": pert,
                    "n_cells": len(target),
                    "effect_norm": effect_norm,
                    "target_residual_radius": residual_radius,
                    "control_radius": control_radius,
                    "residual_to_effect_ratio": residual_radius / (effect_norm + 1e-9),
                    "target_to_control_spread_ratio": residual_radius / (control_radius + 1e-9),
                }
            )
    return pd.DataFrame(rows)


def tune_noise_scale(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    residual_bank: np.ndarray,
    max_scale: float,
    seed: int = 89,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    non_control = sorted([ko for ko in frame["ko_target"].unique() if not control_mask(pd.Series([ko])).iloc[0]])
    rows = []
    scales = np.round(np.linspace(0, max_scale, 7), 3)
    for scale in scales:
        improvements = []
        for ko in non_control:
            train_holdouts = set(HOLDOUT_KOS) | {ko}
            if ko in HOLDOUT_KOS:
                continue
            model, _, _, _ = fit_pls_delta(frame, state_cols, terms, train_holdouts)
            true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
            control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
            ctrl = control[rng.integers(0, len(control), size=len(true))]
            delta = model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1)
            noise = residual_bank[rng.integers(0, len(residual_bank), size=len(true))] * scale
            noise = hard_bound_noise(noise, delta, max_fraction=max_scale)
            virtual = ctrl + delta.reshape(1, -1) + noise
            for j in range(len(state_cols)):
                w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
                w_virt = wasserstein_distance(true[:, j], virtual[:, j])
                if w_ctrl > 1e-9:
                    improvements.append(1 - w_virt / w_ctrl)
        rows.append({"noise_scale": scale, "loo_mean_distribution_improvement": float(np.nanmean(improvements))})
    return pd.DataFrame(rows)


def hard_bound_noise(noise: np.ndarray, delta: np.ndarray, max_fraction: float) -> np.ndarray:
    delta_norm = np.linalg.norm(delta) + 1e-9
    noise_norm = np.linalg.norm(noise, axis=1, keepdims=True) + 1e-9
    max_norm = max_fraction * delta_norm
    noise = noise * np.minimum(1.0, max_norm / noise_norm)
    direction = delta.reshape(1, -1) / delta_norm
    projection = noise @ direction.T
    noise = noise - np.minimum(projection, 0) * direction
    return noise


def evaluate_hard_constraint(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    residual_bank: np.ndarray,
    noise_scale: float,
    interval_scale: float,
    seed: int = 97,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    model, train_kos, _, _ = fit_pls_delta(frame, state_cols, terms, set(HOLDOUT_KOS))
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    residual_std = residual_bank.std(axis=0)
    metric_rows, cell_rows, interval_rows = [], [], []
    for ko in HOLDOUT_KOS:
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0:
            continue
        ctrl = control[rng.integers(0, len(control), size=len(true))]
        delta = model.predict(ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1)
        baseline = ctrl + delta.reshape(1, -1)
        noise = residual_bank[rng.integers(0, len(residual_bank), size=len(true))] * noise_scale
        noise = hard_bound_noise(noise, delta, max_fraction=max(noise_scale, 1e-9))
        virtual = baseline + noise
        lower = baseline.mean(axis=0) - interval_scale * residual_std
        upper = baseline.mean(axis=0) + interval_scale * residual_std
        for j, feature in enumerate(state_cols):
            w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
            w_base = wasserstein_distance(true[:, j], baseline[:, j])
            w_virt = wasserstein_distance(true[:, j], virtual[:, j])
            true_mean = true[:, j].mean()
            interval_rows.append(
                {
                    "ko_target": ko,
                    "feature": feature,
                    "true_mean": true_mean,
                    "predicted_mean": baseline[:, j].mean(),
                    "interval_lower": lower[j],
                    "interval_upper": upper[j],
                    "covered": lower[j] <= true_mean <= upper[j],
                    "interval_width": upper[j] - lower[j],
                }
            )
            for model_name, dist in [("hard_mean", baseline), ("hard_uncertainty_samples", virtual)]:
                w_model = w_base if model_name == "hard_mean" else w_virt
                metric_rows.append(
                    {
                        "model": model_name,
                        "ko_target": ko,
                        "feature": feature,
                        "wasserstein_true_vs_virtual": w_model,
                        "wasserstein_true_vs_control": w_ctrl,
                        "distribution_improvement": 1.0 - w_model / w_ctrl if w_ctrl > 1e-9 else np.nan,
                    }
                )
        for state, matrix in [("control cells", ctrl), ("true KO cells", true), ("hard mean cells", baseline), ("uncertainty samples", virtual)]:
            take = min(150, len(matrix))
            idx = rng.choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=state_cols)
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True), pd.DataFrame(interval_rows)


def plot_summary(metrics: pd.DataFrame, intervals: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(mean_distribution_improvement=("distribution_improvement", "mean"), improved_fraction=("improved", "mean"))
        .reset_index()
    )
    coverage = intervals["covered"].mean()
    summary.to_csv("results/papalexi_hard_constraint_uncertainty_summary.csv", index=False)
    summary["label"] = summary["model"].map({"hard_mean": "Hard mean", "hard_uncertainty_samples": "Uncertainty samples"})
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2))
    sns.barplot(data=summary, x="label", y="mean_distribution_improvement", ax=axes[0], color="#4E79A7")
    sns.barplot(data=summary, x="label", y="improved_fraction", ax=axes[1], color="#4E79A7")
    axes[0].axhline(0, color="0.25", linewidth=1)
    axes[0].set_title("Mean distribution improvement")
    axes[1].set_title("Fraction improved")
    for ax in axes:
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", rotation=12)
    fig.suptitle(f"Hard-constrained generator; interval coverage={coverage:.1%}", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "papalexi_hard_constraint_uncertainty_summary.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_intervals(intervals: pd.DataFrame) -> None:
    selected = ["pathway_IFNG_JAK_STAT", "pathway_IMMUNE_CHECKPOINT", "protein_PDL1", "protein_CD86"]
    plot = intervals.loc[intervals["feature"].isin(selected)].copy()
    plot["feature_label"] = plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    plot["panel"] = plot["ko_target"] + "\n" + plot["feature_label"]
    n = len(plot)
    fig, axes = plt.subplots(4, 4, figsize=(12, 9), sharex=False)
    for ax, (_, row) in zip(axes.flat, plot.iterrows()):
        ax.hlines(0, row["interval_lower"], row["interval_upper"], color="#4E79A7", linewidth=5, alpha=0.75)
        ax.scatter(row["predicted_mean"], 0, color="#E76F51", s=55, label="predicted mean")
        ax.scatter(row["true_mean"], 0, color="#2A9D8F", s=55, label="true mean")
        ax.set_yticks([])
        ax.set_title(row["panel"], fontsize=10)
    for ax in axes.flat[n:]:
        ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.01))
    fig.suptitle("Hard constraint uncertainty intervals: true mean should fall inside the band", y=1.05)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "papalexi_hard_constraint_uncertainty_intervals.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_distributions(cells: pd.DataFrame) -> None:
    selected = ["pathway_IFNG_JAK_STAT", "protein_PDL1", "protein_CD86"]
    selected = [feature for feature in selected if feature in cells.columns]
    plot = cells.melt(id_vars=["ko_target", "state"], value_vars=selected, var_name="feature", value_name="score")
    plot["panel"] = plot["ko_target"] + "\n" + plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    g = sns.catplot(data=plot, x="state", y="score", col="panel", col_wrap=3, kind="box", showfliers=False, height=3.2, aspect=1.05, sharey=False, color="#8DA0CB")
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=35, labelsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Hard-constrained uncertainty samples stay near the fixed KO direction", y=1.02)
    g.savefig(FIG_DIR / "papalexi_hard_constraint_uncertainty_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_doc(external: pd.DataFrame, tuning: pd.DataFrame, metrics: pd.DataFrame, intervals: pd.DataFrame, max_scale: float, chosen_scale: float, interval_scale: float) -> None:
    external_summary = external.groupby("dataset", observed=True).agg(
        n_perturbations=("perturbation", "nunique"),
        median_residual_to_effect=("residual_to_effect_ratio", "median"),
        median_target_to_control_spread=("target_to_control_spread_ratio", "median"),
    )
    model_summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(mean_distribution_improvement=("distribution_improvement", "mean"), improved_fraction=("improved", "mean"))
    )
    text = f"""# Hard constraint uncertainty generator

这一步把 residual/PLS baseline 作为 hard constraint：KO 方向由 PLS 决定，生成模型不能把方向推翻，只能在这个方向附近给出不确定性范围。

## 引入的外部 perturbation 数据

外部数据来自 scPerturb/Zenodo 的标准化 h5ad 集合：

- DatlingerBock2021
- DixitRegev2016_K562_TFs_13_days

这些数据不用于学习 Papalexi 的 KO 方向，只用于估计单细胞扰动后“细胞云团围绕平均方向的波动比例”。

外部不确定性摘要：

{external_summary.round(3).to_string()}

根据外部数据和保守上限，本轮允许的最大噪声比例为：`{max_scale:.3f}`。

## 训练 KO 留一调参

候选噪声强度中，最佳值为：`{chosen_scale:.3f}`。

{tuning.round(3).to_string(index=False)}

## Held-out KO 测试结果

{model_summary.round(3).to_string()}

不确定性区间覆盖率：`{intervals['covered'].mean():.1%}`

区间宽度系数：`{interval_scale:.2f}`

## 图

- `results/figures/papalexi_hard_constraint_uncertainty_summary.png`
- `results/figures/papalexi_hard_constraint_uncertainty_intervals.png`
- `results/figures/papalexi_hard_constraint_uncertainty_distributions.png`

## 当前结论

当前数据下，最优噪声强度仍然非常保守。这说明方向预测比随机生成更重要：虚拟 KO 的均值应由 residual/PLS baseline 固定，生成模型现阶段更适合作为 uncertainty band，而不是自由改变细胞状态。
"""
    Path("docs/hard_constraint_uncertainty_generator.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    frame, state_cols = load_papalexi_state()
    perturb_genes = {gene for ko in frame["ko_target"].unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), perturb_genes)
    _, train_kos, _, _ = fit_pls_delta(frame, state_cols, terms, set(HOLDOUT_KOS))
    residual_bank = papalexi_residual_bank(frame, state_cols, train_kos)

    external = estimate_external_uncertainty()
    external.to_csv("results/external_perturbation_uncertainty_ratios.csv", index=False)
    median_ratio = float(external["target_to_control_spread_ratio"].median()) if len(external) else 1.0
    max_scale = float(np.clip(0.12 * median_ratio, 0.03, 0.20))
    tuning = tune_noise_scale(frame, state_cols, terms, residual_bank, max_scale=max_scale)
    tuning.to_csv("results/papalexi_hard_constraint_uncertainty_tuning.csv", index=False)
    chosen_scale = float(tuning.sort_values("loo_mean_distribution_improvement", ascending=False)["noise_scale"].iloc[0])
    interval_scale = 1.25

    metrics, cells, intervals = evaluate_hard_constraint(frame, state_cols, terms, residual_bank, chosen_scale, interval_scale)
    metrics.to_csv("results/papalexi_hard_constraint_uncertainty_metrics.csv", index=False)
    cells.to_csv("results/papalexi_hard_constraint_uncertainty_cells.csv", index=False)
    intervals.to_csv("results/papalexi_hard_constraint_uncertainty_intervals.csv", index=False)

    plot_summary(metrics, intervals)
    plot_intervals(intervals)
    plot_distributions(cells)
    write_doc(external, tuning, metrics, intervals, max_scale, chosen_scale, interval_scale)

    print(f"external median target/control spread ratio={median_ratio:.3f}; max_scale={max_scale:.3f}; chosen_scale={chosen_scale:.3f}")
    print(metrics.groupby("model", observed=True)["distribution_improvement"].mean().round(3).to_string())
    print(f"interval coverage={intervals['covered'].mean():.1%}")
    print("Saved hard-constraint uncertainty generator results.")


if __name__ == "__main__":
    main()
