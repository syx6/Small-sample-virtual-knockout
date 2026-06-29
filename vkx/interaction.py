from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from sklearn.linear_model import RidgeCV
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler

from .core import parse_gmt, split_ko
from .visualization import setup_plot


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")


def _mechanism_weight(term: str, library: str) -> float:
    upper = term.upper()
    weight = 1.0
    if library in {"tf_target", "atac_motif_tf", "motif_tf_target"}:
        weight *= 1.45
    if any(word in upper for word in ["TGFB", "TGF_BETA", "TGF-BETA", "SMAD"]):
        weight *= 1.35
    if any(word in upper for word in ["MAPK", "ERK", "JNK", "P38", "RAS", "RAF", "MEK"]):
        weight *= 1.30
    return weight


def _select_terms(prior_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 180) -> list[tuple[str, set[str], float]]:
    selected = []
    for path in sorted(prior_dir.glob("*.gmt")):
        scored = []
        for term, genes in parse_gmt(path, include_term_gene=path.stem == "ppi_hub"):
            overlap = len(genes & perturb_genes)
            if overlap == 0 or len(genes) < 5 or len(genes) > 800:
                continue
            weight = _mechanism_weight(term, path.stem) * (1.0 + min(1.0, 8.0 * overlap / max(1, len(genes))))
            scored.append(((overlap, weight, -len(genes)), f"{path.stem}:{term}", genes, weight))
        scored.sort(reverse=True, key=lambda item: item[0])
        selected.extend((name, genes, weight) for _, name, genes, weight in scored[:max_terms_per_library])
    return selected


def _term_features(labels: pd.Series, terms: list[tuple[str, set[str]]]) -> sparse.csr_matrix:
    rows, cols, data = [], [], []
    for i, label in enumerate(labels.astype(str)):
        genes = set(split_ko(label))
        denom = max(1, len(genes))
        for j, term_entry in enumerate(terms):
            members = term_entry[1]
            weight = float(term_entry[2]) if len(term_entry) >= 3 else 1.0
            overlap = len(genes & members)
            if overlap:
                rows.append(i)
                cols.append(j)
                data.append(weight * overlap / denom)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(terms)))


def _pair_interaction_features(labels: pd.Series, terms: list[tuple[str, set[str]]]) -> sparse.csr_matrix:
    rows, cols, data = [], [], []
    for i, label in enumerate(labels.astype(str)):
        genes = split_ko(label)
        if len(genes) < 2:
            continue
        g1, g2 = genes[:2]
        for j, term_entry in enumerate(terms):
            members = term_entry[1]
            weight = float(term_entry[2]) if len(term_entry) >= 3 else 1.0
            if g1 in members and g2 in members:
                rows.append(i)
                cols.append(j)
                data.append(weight)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(terms)))


def _combo_features(combo: pd.DataFrame, ko_col: str, terms: list[tuple[str, set[str]]], gene_mlb: MultiLabelBinarizer | None = None):
    genes = [split_ko(label) for label in combo[ko_col].astype(str)]
    if gene_mlb is None:
        gene_mlb = MultiLabelBinarizer()
        gene = sparse.csr_matrix(gene_mlb.fit_transform(genes), dtype=float)
    else:
        gene = sparse.csr_matrix(gene_mlb.transform(genes), dtype=float)
    term = _term_features(combo[ko_col], terms)
    pair = _pair_interaction_features(combo[ko_col], terms)
    n_gene = sparse.csr_matrix(np.asarray([[len(g)] for g in genes], dtype=float))
    return sparse.hstack([gene, term, pair, n_gene], format="csr"), gene_mlb


def _additive_predictions(single: pd.DataFrame, combo: pd.DataFrame, ko_col: str, target_cols: list[str]) -> np.ndarray:
    gene_to_delta = {}
    for _, row in single.iterrows():
        genes = split_ko(row[ko_col])
        if len(genes) == 1:
            gene_to_delta[genes[0]] = row[target_cols].to_numpy(dtype=float)
    mean_single = single[target_cols].mean().to_numpy(dtype=float)
    pred = []
    for label in combo[ko_col]:
        values = [gene_to_delta.get(gene, mean_single) for gene in split_ko(label)]
        pred.append(np.sum(values, axis=0))
    return np.asarray(pred)


def _loo_interaction_residual(
    combo: pd.DataFrame,
    additive: np.ndarray,
    ko_col: str,
    target_cols: list[str],
    terms: list[tuple[str, set[str]]],
) -> np.ndarray:
    truth = combo[target_cols].to_numpy(dtype=float)
    residual = truth - additive
    preds = []
    alphas = np.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])
    for i in range(len(combo)):
        train_idx = np.arange(len(combo)) != i
        test_idx = np.array([i])
        x_train, mlb = _combo_features(combo.iloc[train_idx].reset_index(drop=True), ko_col, terms)
        x_test, _ = _combo_features(combo.iloc[test_idx].reset_index(drop=True), ko_col, terms, mlb)
        x_scaler = StandardScaler(with_mean=False)
        y_scaler = StandardScaler()
        xt = x_scaler.fit_transform(x_train)
        xv = x_scaler.transform(x_test)
        yt = y_scaler.fit_transform(residual[train_idx])
        model = RidgeCV(alphas=alphas)
        model.fit(xt, yt)
        pred_residual = y_scaler.inverse_transform(model.predict(xv)).reshape(-1)
        preds.append(additive[i] + pred_residual)
    return np.vstack(preds)


def _metric_rows(combo: pd.DataFrame, pred: np.ndarray, ko_col: str, target_cols: list[str], model: str) -> list[dict]:
    rows = []
    truth = combo[target_cols].to_numpy(dtype=float)
    seen_col = "all_genes_seen_in_single"
    subsets = {"all_combos": np.ones(len(combo), dtype=bool)}
    if seen_col in combo.columns:
        subsets["all_genes_seen"] = combo[seen_col].to_numpy(dtype=bool)
        subsets["has_unseen_gene"] = ~combo[seen_col].to_numpy(dtype=bool)
    for subset, mask in subsets.items():
        if mask.sum() < 3:
            continue
        y = truth[mask]
        p = pred[mask]
        for i, col in enumerate(target_cols):
            label = np.abs(y[:, i]) >= 0.15
            roc = np.nan
            pr = np.nan
            if label.sum() > 0 and (~label).sum() > 0:
                roc = roc_auc_score(label, np.abs(p[:, i]))
                pr = average_precision_score(label, np.abs(p[:, i]))
            rows.append(
                {
                    "model": model,
                    "subset": subset,
                    "n_combos": int(mask.sum()),
                    "target": col,
                    "mae": mean_absolute_error(y[:, i], p[:, i]),
                    "r2": r2_score(y[:, i], p[:, i]),
                    "roc_auc_abs_gt_0.15": roc,
                    "pr_auc_abs_gt_0.15": pr,
                }
            )
    return rows


def _plot_double_interaction(metrics: pd.DataFrame, out_dir: Path) -> None:
    setup_plot()
    subset = metrics.loc[metrics["subset"] == "all_combos"].copy()
    subset["program"] = (
        subset["target"]
        .str.replace("delta_", "", regex=False)
        .str.replace("program_", "", regex=False)
        .str.replace("pathway_", "", regex=False)
    )
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    palette = {"single_gene_additive": "#9e9e9e", "interaction_residual": "#377eb8"}
    sns.barplot(data=subset, x="r2", y="program", hue="model", palette=palette, ax=axes[0])
    axes[0].axvline(0, color="0.25", linewidth=1, linestyle="--")
    axes[0].set_title("R2: effect magnitude")
    axes[0].set_xlabel("Higher is better")
    axes[0].set_ylabel("")
    sns.barplot(data=subset, x="mae", y="program", hue="model", palette=palette, ax=axes[1])
    axes[1].set_title("MAE: prediction error")
    axes[1].set_xlabel("Lower is better")
    axes[1].set_ylabel("")
    sns.barplot(data=subset, x="roc_auc_abs_gt_0.15", y="program", hue="model", palette=palette, ax=axes[2])
    axes[2].axvline(0.5, color="0.25", linewidth=1, linestyle="--")
    axes[2].set_title("ROC-AUC: strong response")
    axes[2].set_xlabel("Higher is better")
    axes[2].set_ylabel("")
    axes[0].legend(title="", loc="lower right")
    axes[1].legend_.remove()
    axes[2].legend_.remove()
    fig.suptitle("Double-KO interaction residual model", fontsize=16)
    fig.savefig(out_dir / "double_interaction_metrics.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def run_double_interaction(
    delta_csv: str | Path,
    ko_col: str,
    n_ko_col: str,
    target_prefix: str,
    prior_dir: str | Path,
    out_dir: str | Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    delta = pd.read_csv(delta_csv)
    target_cols = [col for col in delta.columns if col.startswith(target_prefix)]
    if not target_cols:
        raise ValueError(f"No target columns found with prefix '{target_prefix}'.")
    if ko_col not in delta.columns or n_ko_col not in delta.columns:
        raise ValueError(f"Input must contain '{ko_col}' and '{n_ko_col}'.")

    single = delta.loc[delta[n_ko_col] == 1].copy()
    combo = delta.loc[delta[n_ko_col] == 2].copy().reset_index(drop=True)
    if len(single) < 3 or len(combo) < 3:
        raise ValueError("Need at least 3 single-gene rows and 3 double-gene rows.")
    seen = {split_ko(label)[0] for label in single[ko_col] if len(split_ko(label)) == 1}
    combo["all_genes_seen_in_single"] = combo[ko_col].map(lambda label: all(gene in seen for gene in split_ko(label)))
    perturb_genes = {gene for label in delta[ko_col] for gene in split_ko(label)}
    terms = _select_terms(Path(prior_dir), perturb_genes)

    additive = _additive_predictions(single, combo, ko_col, target_cols)
    interaction = _loo_interaction_residual(combo, additive, ko_col, target_cols, terms)
    metrics = pd.DataFrame(
        _metric_rows(combo, additive, ko_col, target_cols, "single_gene_additive")
        + _metric_rows(combo, interaction, ko_col, target_cols, "interaction_residual")
    )
    out = combo[[ko_col, "all_genes_seen_in_single"]].copy()
    for optional in ["n_cells", n_ko_col]:
        if optional in combo.columns and optional not in out.columns:
            out[optional] = combo[optional]
    truth = combo[target_cols].to_numpy(dtype=float)
    for i, col in enumerate(target_cols):
        out[f"true_{col}"] = truth[:, i]
        out[f"additive_pred_{col}"] = additive[:, i]
        out[f"interaction_pred_{col}"] = interaction[:, i]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(out_dir / "double_interaction_metrics.csv", index=False)
    out.to_csv(out_dir / "double_interaction_predictions.csv", index=False)
    _plot_double_interaction(metrics, out_dir)
    _write_double_interaction_report(metrics, out_dir)
    return metrics, out


def _write_double_interaction_report(metrics: pd.DataFrame, out_dir: Path) -> None:
    summary = (
        metrics.groupby(["model", "subset"], observed=True)
        .agg(
            mean_mae=("mae", "mean"),
            mean_r2=("r2", "mean"),
            mean_roc_auc=("roc_auc_abs_gt_0.15", "mean"),
            mean_pr_auc=("pr_auc_abs_gt_0.15", "mean"),
        )
        .reset_index()
    )
    text = f"""# Double-KO interaction residual report

This command evaluates double-gene virtual knockout effects with two models:

- `single_gene_additive`: sum of the available single-gene effects.
- `interaction_residual`: additive baseline plus a prior-based gene-gene interaction residual.

Summary:

{summary.round(3).to_string(index=False)}

Generated files:

- `double_interaction_metrics.csv`
- `double_interaction_predictions.csv`
- `double_interaction_metrics.png`

Interpretation:

Use the interaction model when it lowers MAE and improves R2/ROC-AUC over the additive baseline. It is most useful for double knockouts whose effects are not simple sums of two single-gene knockouts.
"""
    (out_dir / "double_interaction_report.md").write_text(text, encoding="utf-8")
