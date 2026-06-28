from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from sklearn.linear_model import RidgeCV
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")
FIG_DIR = Path("results/figures")


def setup_plot() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"]:
        if font in available_fonts:
            plt.rcParams["font.family"] = font
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300


def split_genes(label: str) -> list[str]:
    return [part.strip().upper() for part in str(label).split("+") if part.strip() and part.lower() != "ctrl"]


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


def select_terms(prior_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 180) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(prior_dir.glob("*.gmt")):
        scored = []
        for term, genes in parse_gmt(path, include_term_gene=path.stem == "ppi_hub"):
            overlap = len(genes & perturb_genes)
            if overlap == 0 or len(genes) < 5 or len(genes) > 800:
                continue
            scored.append(((overlap, -len(genes)), f"{path.stem}:{term}", genes))
        scored.sort(reverse=True, key=lambda item: item[0])
        selected.extend((name, genes) for _, name, genes in scored[:max_terms_per_library])
    return selected


def build_term_features(labels: pd.Series, terms: list[tuple[str, set[str]]]) -> sparse.csr_matrix:
    rows, cols, data = [], [], []
    for i, label in enumerate(labels):
        genes = set(split_genes(label))
        denom = max(1, len(genes))
        for j, (_, members) in enumerate(terms):
            overlap = len(genes & members)
            if overlap:
                rows.append(i)
                cols.append(j)
                data.append(overlap / denom)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(terms)))


def build_pair_interaction_features(labels: pd.Series, terms: list[tuple[str, set[str]]]) -> sparse.csr_matrix:
    rows, cols, data = [], [], []
    for i, label in enumerate(labels):
        genes = split_genes(label)
        if len(genes) < 2:
            continue
        g1, g2 = genes[:2]
        for j, (_, members) in enumerate(terms):
            hit1 = g1 in members
            hit2 = g2 in members
            if hit1 and hit2:
                rows.append(i)
                cols.append(j)
                data.append(1.0)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(terms)))


def additive_predictions(single: pd.DataFrame, combo: pd.DataFrame, target_cols: list[str]) -> np.ndarray:
    gene_to_delta = {}
    for _, row in single.iterrows():
        genes = split_genes(row["ko_genes"])
        if len(genes) == 1:
            gene_to_delta[genes[0]] = row[target_cols].to_numpy(dtype=float)
    mean_single = single[target_cols].mean().to_numpy(dtype=float)
    pred = []
    for label in combo["ko_genes"]:
        values = [gene_to_delta.get(gene, mean_single) for gene in split_genes(label)]
        pred.append(np.sum(values, axis=0))
    return np.asarray(pred)


def build_combo_features(combo: pd.DataFrame, terms: list[tuple[str, set[str]]], gene_mlb: MultiLabelBinarizer | None = None):
    genes = [split_genes(label) for label in combo["ko_genes"]]
    if gene_mlb is None:
        gene_mlb = MultiLabelBinarizer()
        gene = sparse.csr_matrix(gene_mlb.fit_transform(genes), dtype=float)
    else:
        gene = sparse.csr_matrix(gene_mlb.transform(genes), dtype=float)
    term = build_term_features(combo["ko_genes"], terms)
    pair = build_pair_interaction_features(combo["ko_genes"], terms)
    n_gene = sparse.csr_matrix(np.asarray([[len(g)] for g in genes], dtype=float))
    return sparse.hstack([gene, term, pair, n_gene], format="csr"), gene_mlb


def loo_interaction_residual(combo: pd.DataFrame, additive: np.ndarray, target_cols: list[str], terms: list[tuple[str, set[str]]]) -> np.ndarray:
    truth = combo[target_cols].to_numpy(dtype=float)
    residual = truth - additive
    preds = []
    alphas = np.array([0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0])
    for i in range(len(combo)):
        train_idx = np.arange(len(combo)) != i
        test_idx = np.array([i])
        x_train, mlb = build_combo_features(combo.iloc[train_idx].reset_index(drop=True), terms)
        x_test, _ = build_combo_features(combo.iloc[test_idx].reset_index(drop=True), terms, mlb)
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


def metric_rows(combo: pd.DataFrame, pred: np.ndarray, target_cols: list[str], model: str) -> list[dict]:
    rows = []
    truth = combo[target_cols].to_numpy(dtype=float)
    subsets = {
        "all_combos": np.ones(len(combo), dtype=bool),
        "all_genes_seen": combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
        "has_unseen_gene": ~combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
    }
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


def plot_metrics(metrics: pd.DataFrame) -> None:
    subset = metrics.loc[metrics["subset"] == "all_combos"].copy()
    subset["program"] = subset["target"].str.replace("delta_program_", "", regex=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    sns.barplot(data=subset, x="r2", y="program", hue="model", ax=axes[0])
    axes[0].axvline(0, color="0.25", linewidth=1)
    axes[0].set_title("Double-KO continuous effect prediction")
    axes[0].set_xlabel("R2")
    axes[0].set_ylabel("")
    sns.barplot(data=subset, x="roc_auc_abs_gt_0.15", y="program", hue="model", ax=axes[1])
    axes[1].axvline(0.5, color="0.25", linestyle="--", linewidth=1)
    axes[1].set_title("Strong-response ranking")
    axes[1].set_xlabel("ROC-AUC")
    axes[1].set_ylabel("")
    axes[0].legend(title="", loc="lower right")
    axes[1].legend_.remove()
    fig.savefig(FIG_DIR / "norman_double_interaction_model_comparison.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def write_doc(metrics: pd.DataFrame) -> None:
    all_summary = (
        metrics.groupby(["model", "subset"], observed=True)
        .agg(
            mean_mae=("mae", "mean"),
            mean_r2=("r2", "mean"),
            mean_roc_auc=("roc_auc_abs_gt_0.15", "mean"),
            mean_pr_auc=("pr_auc_abs_gt_0.15", "mean"),
        )
        .reset_index()
    )
    text = f"""# Norman double-KO interaction model optimization

This experiment adds a constrained interaction correction for double-gene knockouts.

Method:

```text
single-gene additive prediction
+ prior-based gene-gene interaction residual
= optimized double-KO prediction
```

The correction is trained by leave-one-combo-out evaluation across 52 Norman double-gene KO combinations.

Summary:

{all_summary.round(3).to_string(index=False)}

Key output files:

- `results/norman_double_interaction_metrics.csv`
- `results/norman_double_interaction_predictions.csv`
- `results/figures/norman_double_interaction_model_comparison.png`

Interpretation:

- If `interaction_residual` improves R2 or ROC-AUC over `single_gene_additive`, then explicit gene-gene interaction features help.
- If it improves some programs but worsens others, the interaction layer should be used selectively rather than as a universal replacement.
"""
    Path("docs/norman_double_interaction_model.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    delta = pd.read_csv("results/norman_program_delta.csv")
    target_cols = [col for col in delta.columns if col.startswith("delta_program_")]
    single = delta.loc[delta["n_ko_genes"] == 1].copy()
    combo = delta.loc[delta["n_ko_genes"] == 2].copy().reset_index(drop=True)
    seen = {split_genes(label)[0] for label in single["ko_genes"] if len(split_genes(label)) == 1}
    combo["all_genes_seen_in_single"] = combo["ko_genes"].map(lambda label: all(gene in seen for gene in split_genes(label)))
    perturb_genes = {gene for label in delta["ko_genes"] for gene in split_genes(label)}
    terms = select_terms(Path("data/priors"), perturb_genes)

    additive = additive_predictions(single, combo, target_cols)
    interaction = loo_interaction_residual(combo, additive, target_cols, terms)

    metric = pd.DataFrame(
        metric_rows(combo, additive, target_cols, "single_gene_additive")
        + metric_rows(combo, interaction, target_cols, "interaction_residual")
    )
    metric.to_csv("results/norman_double_interaction_metrics.csv", index=False)
    out = combo[["ko_genes", "n_cells", "all_genes_seen_in_single"]].copy()
    truth = combo[target_cols].to_numpy(dtype=float)
    for i, col in enumerate(target_cols):
        out[f"true_{col}"] = truth[:, i]
        out[f"additive_pred_{col}"] = additive[:, i]
        out[f"interaction_pred_{col}"] = interaction[:, i]
    out.to_csv("results/norman_double_interaction_predictions.csv", index=False)
    plot_metrics(metric)
    write_doc(metric)
    print(metric.groupby(["model", "subset"], observed=True)[["mae", "r2", "roc_auc_abs_gt_0.15", "pr_auc_abs_gt_0.15"]].mean().round(3).to_string())
    print("Saved Norman double-KO interaction model results.")


if __name__ == "__main__":
    main()
