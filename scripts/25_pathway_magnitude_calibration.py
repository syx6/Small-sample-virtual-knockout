from __future__ import annotations

import importlib.util
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wasserstein_distance


FIG_DIR = Path("results/figures")


def load_demo_module():
    path = Path("scripts/23_multi_dataset_virtual_ko_demo.py")
    spec = importlib.util.spec_from_file_location("multi_dataset_demo", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


m23 = load_demo_module()


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


def dataset_short(dataset: str) -> str:
    return {
        "Papalexi ECCITE-seq": "Papalexi",
        "Norman Perturb-seq": "Norman",
        "Datlinger CRISPR RNA": "Datlinger",
        "Dixit Perturb-seq RNA": "Dixit",
    }.get(dataset, dataset)


def state_delta_table(
    frame: pd.DataFrame,
    state_cols: list[str],
    holdouts: list[str],
) -> tuple[np.ndarray, list[str], np.ndarray, list[str], np.ndarray]:
    control_mean = frame.loc[m23.control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    labels, delta = [], []
    for ko, group in frame.loc[~m23.control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdouts:
            continue
        labels.append(str(ko))
        delta.append(group[state_cols].mean().to_numpy(dtype=float) - control_mean)
    return control_mean, labels, np.vstack(delta), holdouts, control_mean


def fit_predict_delta(
    train_labels: list[str],
    train_delta: np.ndarray,
    query_labels: list[str],
    terms: list[tuple[str, set[str]]],
) -> dict[str, np.ndarray]:
    model = m23.fit_pls(train_labels, train_delta, terms)
    return {
        ko: model.predict(m23.ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1)
        for ko in query_labels
    }


def loo_training_predictions(
    train_labels: list[str],
    train_delta: np.ndarray,
    terms: list[tuple[str, set[str]]],
) -> np.ndarray:
    preds = []
    for i, ko in enumerate(train_labels):
        keep = np.arange(len(train_labels)) != i
        if keep.sum() < 3:
            preds.append(train_delta[i])
            continue
        kept_labels = [label for j, label in enumerate(train_labels) if keep[j]]
        model = m23.fit_pls(kept_labels, train_delta[keep], terms)
        preds.append(model.predict(m23.ko_prior_vector(ko, terms).reshape(1, -1)).reshape(-1))
    return np.vstack(preds)


def calibration_factors(loo_pred: np.ndarray, truth: np.ndarray) -> dict[str, np.ndarray]:
    eps = 1e-9
    global_alpha = np.sum(loo_pred * truth) / (np.sum(loo_pred * loo_pred) + eps)
    # Keep a non-zero KO effect. A fitted alpha of exactly zero can reduce error
    # by predicting "no change", but that is not useful for virtual knockout.
    global_alpha = float(np.clip(global_alpha, 0.15, 2.5))
    denom = np.sum(loo_pred * loo_pred, axis=0) + eps
    feature_alpha = np.sum(loo_pred * truth, axis=0) / denom
    feature_alpha = np.clip(feature_alpha, 0.15, 2.5)
    return {
        "uncalibrated": np.ones(truth.shape[1]),
        "global_scale": np.full(truth.shape[1], global_alpha),
        "feature_scale": feature_alpha,
    }


def choose_method(loo_pred: np.ndarray, truth: np.ndarray, factors: dict[str, np.ndarray]) -> tuple[str, pd.DataFrame]:
    rows = []
    for method, alpha in factors.items():
        pred = loo_pred * alpha.reshape(1, -1)
        rows.append(
            {
                "calibration_method": method,
                "training_delta_mae": float(np.mean(np.abs(pred - truth))),
                "training_direction_cosine": mean_cosine(pred, truth),
            }
        )
    table = pd.DataFrame(rows).sort_values(["training_delta_mae", "calibration_method"])
    return str(table.iloc[0]["calibration_method"]), table


def mean_cosine(pred: np.ndarray, truth: np.ndarray) -> float:
    scores = []
    for p, t in zip(pred, truth):
        denom = np.linalg.norm(p) * np.linalg.norm(t)
        if denom > 1e-9:
            scores.append(float(np.dot(p, t) / denom))
    return float(np.mean(scores)) if scores else np.nan


def load_papalexi_calibration():
    frame, state_cols, holdouts, _, modality, representation = m23.load_papalexi()
    genes = {gene for ko in frame["ko_target"].unique() for gene in m23.split_ko(ko)}
    terms = m23.select_prior_terms(Path("data/priors"), genes)
    _, train_labels, train_delta, _, _ = state_delta_table(frame, state_cols, holdouts)
    pred = fit_predict_delta(train_labels, train_delta, holdouts, terms)
    return frame, state_cols, holdouts, pred, train_labels, train_delta, terms, modality, representation


def load_norman_calibration():
    adata = ad.read_h5ad("data/norman_small_program.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_genes"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("program_")]:
        frame[col] = obs[col].astype(float).values
    state_cols = [c for c in frame.columns if c.startswith("program_")]
    holdouts = ["AHR+KLF1", "CEBPB+CEBPA", "MAPK1+TGFBR2", "CBL+UBASH3B"]
    genes = {gene for ko in frame["ko_target"].unique() for gene in m23.split_ko(ko)}
    terms = m23.select_prior_terms(Path("data/priors"), genes)
    _, train_labels, train_delta, _, _ = state_delta_table(frame, state_cols, holdouts)
    pred = fit_predict_delta(train_labels, train_delta, holdouts, terms)
    return frame, state_cols, holdouts, pred, train_labels, train_delta, terms, "single-cell RNA perturb-seq", "gene program scores"


def load_external_calibration(path: Path, holdouts: list[str], seed: int):
    rng = np.random.default_rng(seed)
    adata = ad.read_h5ad(path)
    label, is_control = m23.external_labels(adata, path.name)
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
    state_terms = m23.select_state_terms(adata, selected, max_terms=14)
    frame = m23.compute_pathway_scores(adata, keep_idx, state_terms)
    frame["ko_target"] = label.iloc[keep_idx].to_numpy()
    frame.loc[is_control.iloc[keep_idx].to_numpy(), "ko_target"] = "control"
    state_cols = [c for c in frame.columns if c.startswith("pathway_")]
    genes = {gene for ko in selected for gene in m23.split_ko(ko)}
    terms = m23.select_prior_terms(Path("data/priors"), genes)
    _, train_labels, train_delta, _, _ = state_delta_table(frame, state_cols, holdouts)
    pred = fit_predict_delta(train_labels, train_delta, holdouts, terms)
    return frame, state_cols, holdouts, pred, train_labels, train_delta, terms, "single-cell RNA CRISPR screen", "pathway/program scores"


def evaluate_methods(
    dataset: str,
    modality: str,
    representation: str,
    frame: pd.DataFrame,
    state_cols: list[str],
    holdouts: list[str],
    pred: dict[str, np.ndarray],
    factors: dict[str, np.ndarray],
    best_method: str,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    control = frame.loc[m23.control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=float)
    metric_rows, cell_rows, delta_rows = [], [], []
    for ko in holdouts:
        true = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=float)
        if len(true) == 0 or ko not in pred:
            continue
        ctrl = control[rng.integers(0, len(control), size=len(true))]
        true_delta = true.mean(axis=0) - ctrl.mean(axis=0)
        for method, alpha in factors.items():
            delta = pred[ko] * alpha
            virtual = ctrl + delta.reshape(1, -1)
            virt_delta = virtual.mean(axis=0) - ctrl.mean(axis=0)
            denom = np.linalg.norm(true_delta) * np.linalg.norm(virt_delta)
            delta_rows.append(
                {
                    "dataset": dataset,
                    "ko_target": ko,
                    "calibration_method": method,
                    "is_selected_method": method == best_method,
                    "direction_cosine": float(np.dot(true_delta, virt_delta) / denom) if denom > 1e-9 else np.nan,
                    "mean_abs_delta_error": float(np.mean(np.abs(virt_delta - true_delta))),
                    "predicted_delta_norm": float(np.linalg.norm(virt_delta)),
                    "true_delta_norm": float(np.linalg.norm(true_delta)),
                }
            )
            for j, feature in enumerate(state_cols):
                w_ctrl = wasserstein_distance(true[:, j], ctrl[:, j])
                w_virt = wasserstein_distance(true[:, j], virtual[:, j])
                metric_rows.append(
                    {
                        "dataset": dataset,
                        "input_modality": modality,
                        "state_representation": representation,
                        "ko_target": ko,
                        "feature": feature,
                        "calibration_method": method,
                        "is_selected_method": method == best_method,
                        "wasserstein_true_vs_virtual": w_virt,
                        "wasserstein_true_vs_control": w_ctrl,
                        "distribution_improvement": 1.0 - w_virt / w_ctrl if w_ctrl > 1e-9 else np.nan,
                    }
                )
            if method in {"uncalibrated", best_method}:
                for state, matrix in [("control cells", ctrl), ("virtual KO cells", virtual), ("true KO cells", true)]:
                    take = min(180, len(matrix))
                    idx = rng.choice(len(matrix), size=take, replace=False)
                    tmp = pd.DataFrame(matrix[idx], columns=state_cols)
                    tmp["dataset"] = dataset
                    tmp["ko_target"] = ko
                    tmp["state"] = state
                    tmp["calibration_method"] = method
                    cell_rows.append(tmp)
    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True), pd.DataFrame(delta_rows)


def plot_summary(summary: pd.DataFrame) -> None:
    selected = summary.loc[summary["is_selected_method"]].copy()
    uncal = summary.loc[summary["calibration_method"] == "uncalibrated"].copy()
    plot = pd.concat(
        [
            uncal.assign(display_method="before calibration"),
            selected.assign(display_method="after calibration"),
        ],
        ignore_index=True,
    )
    plot["dataset_short"] = plot["dataset"].map(dataset_short)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    sns.barplot(data=plot, x="dataset_short", y="mean_distribution_improvement", hue="display_method", ax=axes[0])
    sns.barplot(data=plot, x="dataset_short", y="mean_abs_delta_error", hue="display_method", ax=axes[1])
    axes[0].axhline(0, color="0.25", linewidth=1)
    axes[0].set_title("Distribution closeness")
    axes[1].set_title("Pathway/program magnitude error")
    axes[0].set_ylabel("mean improvement")
    axes[1].set_ylabel("mean absolute delta error")
    for ax in axes:
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)
    axes[0].legend_.remove()
    axes[1].legend(title="")
    fig.suptitle("Effect of non-negative magnitude calibration", y=1.03)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pathway_magnitude_calibration_summary.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_delta_heatmap(delta_metrics: pd.DataFrame) -> None:
    selected = delta_metrics.loc[
        delta_metrics["is_selected_method"] | (delta_metrics["calibration_method"] == "uncalibrated")
    ].copy()
    selected["label"] = selected["ko_target"] + " (" + selected["dataset"].map(dataset_short) + ")"
    pivot = selected.pivot_table(
        index="label",
        columns="calibration_method",
        values="mean_abs_delta_error",
        aggfunc="mean",
    )
    ordered_cols = [col for col in ["uncalibrated", "global_scale", "feature_scale"] if col in pivot.columns]
    pivot = pivot[ordered_cols]
    fig, ax = plt.subplots(figsize=(7.2, max(5, 0.42 * len(pivot) + 1.5)))
    vmax = float(np.nanmax(pivot.to_numpy()))
    sns.heatmap(pivot, cmap="rocket_r", vmin=0, vmax=vmax, annot=True, fmt=".2f", cbar_kws={"label": "MAE"}, ax=ax)
    ax.set_title("Magnitude error before and after calibration")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pathway_magnitude_calibration_delta_error_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def feature_label(feature: str) -> str:
    for prefix in ["pathway_", "protein_", "program_"]:
        if feature.startswith(prefix):
            return feature.replace(prefix, "")
    return feature


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


def mean_state_tables(cells: pd.DataFrame, dataset: str, method: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sub = cells.loc[(cells["dataset"] == dataset) & (cells["calibration_method"] == method)].copy()
    feature_cols = [
        col
        for col in sub.columns
        if col.startswith("pathway_") or col.startswith("protein_") or col.startswith("program_")
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


def plot_selected_calibrated_heatmaps(cells: pd.DataFrame, factor_table: pd.DataFrame) -> None:
    selected = factor_table.loc[factor_table["selected_for_dataset"]].set_index("dataset")["calibration_method"].to_dict()
    for dataset, method in selected.items():
        true_df, virtual_df, error_df = mean_state_tables(cells, dataset, method)
        if true_df.empty:
            continue
        vmax = np.nanmax(np.abs(pd.concat([true_df, virtual_df]).to_numpy()))
        err_vmax = np.nanmax(np.abs(error_df.to_numpy()))
        fig, axes = plt.subplots(1, 3, figsize=(15.5, max(3.2, 0.55 * len(true_df) + 2.0)), sharey=True)
        sns.heatmap(true_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[0], cbar=False)
        sns.heatmap(virtual_df, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", ax=axes[1], cbar=False)
        sns.heatmap(error_df, cmap="vlag", center=0, vmin=-err_vmax, vmax=err_vmax, annot=True, fmt=".2f", ax=axes[2], cbar_kws={"label": "virtual - true"})
        axes[0].set_title("Real KO change")
        axes[1].set_title("Calibrated virtual KO")
        axes[2].set_title("Prediction error")
        for ax in axes:
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="x", rotation=35)
        fig.suptitle(f"{dataset}: calibrated real vs virtual KO ({method})", y=1.04)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"pathway_magnitude_calibration_heatmap_{dataset_short(dataset).lower()}_{method}.png", bbox_inches="tight", dpi=300)
        plt.close(fig)


def write_doc(summary: pd.DataFrame, factors: pd.DataFrame, loo: pd.DataFrame) -> None:
    selected = summary.loc[summary["is_selected_method"]].copy()
    before = summary.loc[summary["calibration_method"] == "uncalibrated"].copy()
    compare = before.merge(
        selected,
        on=["dataset", "input_modality", "state_representation"],
        suffixes=("_before", "_after"),
    )
    compare = compare[
        [
            "dataset",
            "input_modality",
            "state_representation",
            "mean_distribution_improvement_before",
            "mean_distribution_improvement_after",
            "mean_abs_delta_error_before",
            "mean_abs_delta_error_after",
            "calibration_method_after",
        ]
    ]
    text = f"""# Pathway/program score magnitude calibration

这一步没有改变 KO 方向模型，只做一件事：给 PLS/residual 预测出来的 pathway/program 变化幅度加一个非负倍率。

## 为什么要做

前一版 RNA-only 结果的主要问题不是完全预测错方向，而是预测变化幅度不稳。也就是说，模型知道细胞状态大概应该往哪里移动，但移动多远还需要校准。

## 校准规则

```text
virtual KO state = control state + alpha * predicted KO delta
```

这里的 `alpha` 限制在 `0.15 到 2.5` 之间，所以校准不会把 KO 效应方向反过来，也不会把虚拟敲除缩成“几乎没有敲除”。我们比较了三种版本：

- `uncalibrated`: 不校准，alpha = 1。
- `global_scale`: 每个数据集一个统一倍率。
- `feature_scale`: 每个 pathway/program 一个倍率。

倍率只用训练 KO 的 leave-one-KO-out 误差学习，holdout KO 的真实结果只用于最后评估。

## 结果摘要

{compare.round(3).to_string(index=False)}

## 训练 KO 上选择的校准方式

{loo.round(3).to_string(index=False)}

## 校准倍率

{factors.round(3).to_string(index=False)}

## 图

- `results/figures/pathway_magnitude_calibration_summary.png`
- `results/figures/pathway_magnitude_calibration_delta_error_heatmap.png`
- `results/figures/pathway_magnitude_calibration_heatmap_papalexi_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_norman_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_datlinger_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_dixit_global_scale.png`

## 现在怎么看效果

如果校准后 `mean_abs_delta_error` 下降，说明 pathway/program 变化幅度更接近真实 KO。
如果 `mean_distribution_improvement` 上升，说明虚拟 KO 单细胞分布也更接近真实 KO。
如果某个数据集方向本来就错，幅度校准救不了方向错误；它只负责把已经大致正确的方向调到更合适的大小。
"""
    Path("docs/pathway_magnitude_calibration.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    datasets = [
        ("Papalexi ECCITE-seq", *load_papalexi_calibration()),
        ("Norman Perturb-seq", *load_norman_calibration()),
        ("Datlinger CRISPR RNA", *load_external_calibration(Path("data/scperturb_extra/DatlingerBock2021.h5ad"), ["LAT", "LCK", "JUND", "FOS"], seed=31)),
        ("Dixit Perturb-seq RNA", *load_external_calibration(Path("data/scperturb_extra/DixitRegev2016_K562_TFs_13_days.h5ad"), ["ELF1", "CREB1", "ELK1", "GABPA"], seed=37)),
    ]

    all_metrics, all_cells, all_delta, factor_rows, loo_rows = [], [], [], [], []
    for i, (dataset, frame, state_cols, holdouts, pred, train_labels, train_delta, terms, modality, representation) in enumerate(datasets):
        loo_pred = loo_training_predictions(train_labels, train_delta, terms)
        factors = calibration_factors(loo_pred, train_delta)
        best_method, loo_table = choose_method(loo_pred, train_delta, factors)
        loo_table["dataset"] = dataset
        loo_rows.append(loo_table)
        for method, alpha in factors.items():
            factor_rows.append(
                {
                    "dataset": dataset,
                    "calibration_method": method,
                    "selected_for_dataset": method == best_method,
                    "mean_alpha": float(np.mean(alpha)),
                    "min_alpha": float(np.min(alpha)),
                    "max_alpha": float(np.max(alpha)),
                }
            )
        metrics, cells, delta_metrics = evaluate_methods(
            dataset,
            modality,
            representation,
            frame,
            state_cols,
            holdouts,
            pred,
            factors,
            best_method,
            seed=251 + i,
        )
        all_metrics.append(metrics)
        all_cells.append(cells)
        all_delta.append(delta_metrics)

    metrics = pd.concat(all_metrics, ignore_index=True)
    cells = pd.concat(all_cells, ignore_index=True)
    delta_metrics = pd.concat(all_delta, ignore_index=True)
    factor_table = pd.DataFrame(factor_rows)
    loo_table = pd.concat(loo_rows, ignore_index=True)

    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .merge(delta_metrics[["dataset", "ko_target", "calibration_method", "mean_abs_delta_error"]], on=["dataset", "ko_target", "calibration_method"], how="left")
        .groupby(["dataset", "input_modality", "state_representation", "calibration_method", "is_selected_method"], observed=True)
        .agg(
            n_ko=("ko_target", "nunique"),
            n_features=("feature", "nunique"),
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
            mean_abs_delta_error=("mean_abs_delta_error", "mean"),
        )
        .reset_index()
    )

    metrics.to_csv("results/pathway_magnitude_calibration_metrics.csv", index=False)
    cells.to_csv("results/pathway_magnitude_calibration_cells.csv", index=False)
    delta_metrics.to_csv("results/pathway_magnitude_calibration_delta_metrics.csv", index=False)
    factor_table.to_csv("results/pathway_magnitude_calibration_factors.csv", index=False)
    loo_table.to_csv("results/pathway_magnitude_calibration_loo.csv", index=False)
    summary.to_csv("results/pathway_magnitude_calibration_summary.csv", index=False)

    plot_summary(summary)
    plot_delta_heatmap(delta_metrics)
    plot_selected_calibrated_heatmaps(cells, factor_table)
    write_doc(summary, factor_table, loo_table)

    print(summary.round(3).to_string(index=False))
    print("\nSelected calibration")
    print(factor_table.loc[factor_table["selected_for_dataset"]].round(3).to_string(index=False))
    print("Saved pathway/program magnitude calibration results.")


if __name__ == "__main__":
    main()
