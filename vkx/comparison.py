from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


METHOD_REGISTRY = [
    {
        "method": "VKX hard-constrained multimodal baseline",
        "category": "small-sample interpretable baseline",
        "input_modalities": "RNA, ADT/protein, ATAC gene activity, chromVAR, peak features",
        "multi_gene": "single and double KO; interaction residual for double KO",
        "small_sample_fit": "strong",
        "interpretability": "high",
        "needs_large_training_set": "no",
        "main_strength": "Works as a conservative, explainable small-sample baseline with explicit benchmark/application boundaries.",
        "main_limitation": "Not a free-form high-capacity generative model.",
    },
    {
        "method": "Linear / ridge / PLS baseline",
        "category": "classical baseline",
        "input_modalities": "usually RNA or state scores",
        "multi_gene": "mostly additive unless extended",
        "small_sample_fit": "strong",
        "interpretability": "medium-high",
        "needs_large_training_set": "no",
        "main_strength": "Stable and easy to audit.",
        "main_limitation": "Weak nonlinear and distribution-shape modeling.",
    },
    {
        "method": "scGen",
        "category": "VAE latent shift",
        "input_modalities": "primarily RNA",
        "multi_gene": "not primary focus",
        "small_sample_fit": "medium",
        "interpretability": "medium-low",
        "needs_large_training_set": "usually yes",
        "main_strength": "Learns perturbation shifts in latent space.",
        "main_limitation": "Less tailored to small-sample multimodal and peak-level regulatory priors.",
    },
    {
        "method": "CPA / compositional perturbation autoencoder",
        "category": "deep compositional model",
        "input_modalities": "primarily RNA; covariates supported",
        "multi_gene": "supports compositional perturbations",
        "small_sample_fit": "medium-low",
        "interpretability": "medium",
        "needs_large_training_set": "yes",
        "main_strength": "Models perturbation and covariate composition.",
        "main_limitation": "Training data and tuning requirements are higher than conservative baselines.",
    },
    {
        "method": "GEARS",
        "category": "graph neural perturbation model",
        "input_modalities": "primarily RNA plus gene graph",
        "multi_gene": "strong focus on gene combinations",
        "small_sample_fit": "medium",
        "interpretability": "medium",
        "needs_large_training_set": "usually yes",
        "main_strength": "Strong conceptual fit for combinatorial perturbation prediction.",
        "main_limitation": "Not designed as a lightweight multimodal ATAC/ADT reporting workflow.",
    },
    {
        "method": "CellOT / optimal transport",
        "category": "distribution transport",
        "input_modalities": "usually RNA/state embeddings",
        "multi_gene": "not primary focus",
        "small_sample_fit": "medium",
        "interpretability": "medium-low",
        "needs_large_training_set": "depends",
        "main_strength": "Models distributional transport from control to perturbed states.",
        "main_limitation": "Biological prior, multi-gene interaction, and multimodal reporting need extra engineering.",
    },
    {
        "method": "Free diffusion / flow matching",
        "category": "high-capacity generative model",
        "input_modalities": "varies",
        "multi_gene": "possible with enough data",
        "small_sample_fit": "weak unless strongly constrained",
        "interpretability": "low-medium",
        "needs_large_training_set": "yes",
        "main_strength": "Can model complex cell-level distributions.",
        "main_limitation": "Easy to overfit or drift in small-sample settings; needs careful hard constraints.",
    },
]


METRIC_ALIASES = {
    "roc_auc": ["roc_auc", "auc", "mean_auc", "roc_auc_abs_gt_0.15"],
    "direction": ["direction", "direction_cosine", "mean_direction_cosine"],
    "mae": ["mae", "mean_abs_delta_error", "mean_mae"],
    "r2": ["r2", "mean_r2"],
    "distribution_improvement": ["distribution_improvement", "mean_distribution_improvement"],
    "improved_fraction": ["improved_fraction", "feature_hit_rate"],
}


def method_registry_table() -> pd.DataFrame:
    return pd.DataFrame(METHOD_REGISTRY)


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    for col in candidates:
        if col in columns:
            return col
    return None


def load_metric_tables(metric_csvs: list[str | Path] | None) -> pd.DataFrame:
    if not metric_csvs:
        return pd.DataFrame()
    rows = []
    for path in metric_csvs:
        table = pd.read_csv(path)
        source = Path(path).stem
        method_col = _first_existing(list(table.columns), ["method", "model", "model_name", "variant"])
        dataset_col = _first_existing(list(table.columns), ["dataset", "input_modality", "source_dataset"])
        for _, row in table.iterrows():
            out = {
                "source_file": str(path),
                "method": str(row[method_col]) if method_col else source,
                "dataset": str(row[dataset_col]) if dataset_col else source,
            }
            for metric, aliases in METRIC_ALIASES.items():
                col = _first_existing(list(table.columns), aliases)
                out[metric] = float(row[col]) if col and pd.notna(row[col]) else np.nan
            rows.append(out)
    return pd.DataFrame(rows)


def make_method_comparison(metric_csvs: list[str | Path] | None, out_dir: str | Path) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    registry = method_registry_table()
    metrics = load_metric_tables(metric_csvs)
    registry.to_csv(out / "method_registry.csv", index=False)
    if not metrics.empty:
        metrics.to_csv(out / "method_metric_comparison.csv", index=False)
        _plot_metric_comparison(metrics, out)
    _plot_method_positioning(registry, out)
    _write_comparison_report(registry, metrics, out)
    return {"registry": registry, "metrics": metrics}


def _plot_method_positioning(registry: pd.DataFrame, out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return
    score_map = {"strong": 3, "medium": 2, "medium-low": 1.4, "weak unless strongly constrained": 1.0}
    interp_map = {"high": 3, "medium-high": 2.6, "medium": 2, "medium-low": 1.4, "low-medium": 1.2}
    plot = registry.copy()
    plot["small_sample_score"] = plot["small_sample_fit"].map(score_map).fillna(1.5)
    plot["interpretability_score"] = plot["interpretability"].map(interp_map).fillna(1.5)
    colors = ["#E76F51" if "VKX" in method else "#4C78A8" for method in plot["method"]]
    fig, ax = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
    ax.scatter(plot["small_sample_score"], plot["interpretability_score"], s=120, c=colors, alpha=0.85)
    for _, row in plot.iterrows():
        ax.text(row["small_sample_score"] + 0.035, row["interpretability_score"] + 0.035, row["method"].split(" / ")[0], fontsize=8)
    ax.set_xlim(0.8, 3.35)
    ax.set_ylim(0.8, 3.35)
    ax.set_xlabel("Small-sample suitability")
    ax.set_ylabel("Interpretability / reportability")
    ax.set_title("Method Positioning for Small-sample Multimodal Virtual KO")
    fig.savefig(out_dir / "01_method_positioning.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _plot_metric_comparison(metrics: pd.DataFrame, out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return
    available = [metric for metric in ["roc_auc", "direction", "r2", "mae", "distribution_improvement"] if metrics[metric].notna().any()]
    if not available:
        return
    plot = metrics.melt(id_vars=["method", "dataset"], value_vars=available, var_name="metric", value_name="value").dropna()
    fig, axes = plt.subplots(1, len(available), figsize=(max(5 * len(available), 8), 4.8), squeeze=False, constrained_layout=True)
    for ax, metric in zip(axes.flat, available):
        sub = plot.loc[plot["metric"] == metric]
        order = sub.groupby("method")["value"].mean().sort_values(ascending=(metric == "mae")).index
        sns.barplot(data=sub, x="value", y="method", order=order, color="#4C78A8", ax=ax)
        ax.set_title(metric.replace("_", " ").upper())
        ax.set_xlabel("lower is better" if metric == "mae" else "higher is better")
        ax.set_ylabel("")
    fig.suptitle("Empirical Method / Variant Comparison")
    fig.savefig(out_dir / "02_metric_comparison.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def _write_comparison_report(registry: pd.DataFrame, metrics: pd.DataFrame, out_dir: Path) -> None:
    metric_text = "No metric CSV was provided. This report contains conceptual positioning only."
    if not metrics.empty:
        summary = (
            metrics.groupby("method", observed=True)
            [["roc_auc", "direction", "r2", "mae", "distribution_improvement", "improved_fraction"]]
            .mean(numeric_only=True)
            .round(4)
            .reset_index()
        )
        metric_text = summary.to_string(index=False)
    text = f"""# Virtual KO Method Comparison Report

## Why Compare Methods?

The goal is not to claim that one method dominates every setting. The key question is which method is appropriate for small-sample, multimodal, interpretable virtual knockout.

## Conceptual Registry

{registry.to_string(index=False)}

## Empirical Metrics

{metric_text}

## Recommended Claim Boundary

VKX should be positioned as a small-sample, multimodal, prior-constrained baseline with explicit uncertainty and benchmark/application boundaries. Deep models such as VAE/CPA/GEARS/flow/diffusion are important comparison points, but they usually require more labelled perturbations and careful tuning.

## Generated Figures

- `01_method_positioning.png`
- `02_metric_comparison.png`, when metric CSVs are provided
"""
    (out_dir / "method_comparison_report.md").write_text(text, encoding="utf-8")
