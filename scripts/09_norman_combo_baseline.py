from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler


GENE_MODULES = {
    "erythroid_tf": {"KLF1", "GATA1"},
    "granulocyte_tf": {"CEBPE", "CEBPB", "SPI1"},
    "mapk_tgfb": {"MAPK1", "TGFBR2", "SMAD4"},
    "pioneer_tf": {"FOXA3", "AHR", "CEBPA", "CEBPB", "SPI1", "KLF1"},
    "growth_cell_cycle": {"MYC", "CDK1", "CDKN1A", "CDKN1B"},
    "apoptosis": {"BAK1", "BAX", "CASP3", "CASP8"},
    "ubash": {"UBASH3A", "UBASH3B"},
}


def split_genes(label: str) -> list[str]:
    if str(label).lower() == "ctrl":
        return []
    return [part.strip().upper() for part in str(label).split("+") if part.strip() and part.lower() != "ctrl"]


def infer_control_mask(values: pd.Series) -> pd.Series:
    return values.astype(str).str.lower().eq("ctrl")


def make_delta_table(adata: ad.AnnData) -> pd.DataFrame:
    program_cols = list(adata.uns.get("program_columns", []))
    control = infer_control_mask(adata.obs["ko_genes"])
    base = adata.obs.loc[control, program_cols].mean()
    rows = []
    for target, frame in adata.obs.groupby("ko_genes", observed=True):
        genes = split_genes(target)
        mean_scores = frame[program_cols].mean()
        row = {"ko_genes": target, "n_cells": len(frame), "n_ko_genes": len(genes), "genes": "+".join(genes) or "ctrl"}
        row.update({f"delta_{col}": mean_scores[col] - base[col] for col in program_cols})
        rows.append(row)
    return pd.DataFrame(rows)


def module_features(gene_sets: list[list[str]]) -> tuple[np.ndarray, list[str]]:
    names = [f"module_{name}" for name in GENE_MODULES] + ["n_ko_genes", "same_module_pair"]
    rows = []
    for genes in gene_sets:
        gset = set(genes)
        row = [float(bool(gset & members)) for members in GENE_MODULES.values()]
        row.append(float(len(gset)))
        same = 0.0
        if len(gset) >= 2:
            for members in GENE_MODULES.values():
                if len(gset & members) >= 2:
                    same = 1.0
                    break
        row.append(same)
        rows.append(row)
    return np.asarray(rows, dtype=float), names


def train_additive(single: pd.DataFrame, combo: pd.DataFrame, target_cols: list[str]) -> np.ndarray:
    gene_to_delta = {}
    for _, row in single.iterrows():
        genes = split_genes(row["ko_genes"])
        if len(genes) == 1:
            gene_to_delta[genes[0]] = row[target_cols].to_numpy(dtype=float)
    mean_single = single[target_cols].mean().to_numpy(dtype=float)
    preds = []
    for label in combo["ko_genes"]:
        parts = split_genes(label)
        vals = [gene_to_delta.get(gene, mean_single) for gene in parts]
        preds.append(np.sum(vals, axis=0))
    return np.asarray(preds, dtype=float)


def train_ridge_prior(train: pd.DataFrame, test: pd.DataFrame, target_cols: list[str]) -> np.ndarray:
    train_genes = [split_genes(label) for label in train["ko_genes"]]
    test_genes = [split_genes(label) for label in test["ko_genes"]]
    mlb = MultiLabelBinarizer()
    x_gene_train = mlb.fit_transform(train_genes)
    x_gene_test = mlb.transform(test_genes)
    x_mod_train, _ = module_features(train_genes)
    x_mod_test, _ = module_features(test_genes)
    x_train = np.hstack([x_gene_train, x_mod_train])
    x_test = np.hstack([x_gene_test, x_mod_test])
    y_train = train[target_cols].to_numpy(dtype=float)

    x_scaler = StandardScaler(with_mean=False)
    y_scaler = StandardScaler()
    x_train = x_scaler.fit_transform(x_train)
    y_train = y_scaler.fit_transform(y_train)
    model = Ridge(alpha=1.0)
    model.fit(x_train, y_train)
    return y_scaler.inverse_transform(model.predict(x_scaler.transform(x_test)))


def metric_rows(y: np.ndarray, pred: np.ndarray, cols: list[str], model: str) -> list[dict]:
    rows = []
    for i, col in enumerate(cols):
        row = {
            "model": model,
            "target": col,
            "mae": mean_absolute_error(y[:, i], pred[:, i]),
            "r2": r2_score(y[:, i], pred[:, i]),
            "roc_auc_abs_gt_0.15": np.nan,
            "pr_auc_abs_gt_0.15": np.nan,
        }
        labels = np.abs(y[:, i]) >= 0.15
        if labels.sum() > 0 and (~labels).sum() > 0:
            score = np.abs(pred[:, i])
            row["roc_auc_abs_gt_0.15"] = roc_auc_score(labels, score)
            row["pr_auc_abs_gt_0.15"] = average_precision_score(labels, score)
        rows.append(row)
    return rows


def metric_rows_by_subset(combo: pd.DataFrame, pred: np.ndarray, cols: list[str], model: str) -> list[dict]:
    rows = []
    y_all = combo[cols].to_numpy(dtype=float)
    for subset_name, mask in {
        "all_combos": np.ones(len(combo), dtype=bool),
        "all_genes_seen": combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
        "has_unseen_gene": ~combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
    }.items():
        if mask.sum() < 3:
            continue
        for row in metric_rows(y_all[mask], pred[mask], cols, model):
            row["subset"] = subset_name
            row["n_combos"] = int(mask.sum())
            rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/norman_small_program.h5ad")
    parser.add_argument("--delta", default="results/norman_program_delta.csv")
    parser.add_argument("--metrics", default="results/norman_combo_metrics.csv")
    parser.add_argument("--predictions", default="results/norman_combo_predictions.csv")
    args = parser.parse_args()

    adata = ad.read_h5ad(args.input)
    delta = make_delta_table(adata)
    Path(args.delta).parent.mkdir(parents=True, exist_ok=True)
    delta.to_csv(args.delta, index=False)

    target_cols = [col for col in delta.columns if col.startswith("delta_program_")]
    single = delta.loc[delta["n_ko_genes"] == 1].copy()
    combo = delta.loc[delta["n_ko_genes"] == 2].copy()
    seen_single_genes = {split_genes(label)[0] for label in single["ko_genes"] if len(split_genes(label)) == 1}
    combo["all_genes_seen_in_single"] = combo["ko_genes"].map(
        lambda label: all(gene in seen_single_genes for gene in split_genes(label))
    )
    y = combo[target_cols].to_numpy(dtype=float)

    additive_pred = train_additive(single, combo, target_cols)
    ridge_pred = train_ridge_prior(pd.concat([single, delta.loc[delta["n_ko_genes"] == 0]], ignore_index=True), combo, target_cols)

    metrics = pd.DataFrame(
        metric_rows_by_subset(combo, additive_pred, target_cols, "single_gene_additive")
        + metric_rows_by_subset(combo, ridge_pred, target_cols, "single_gene_ridge_prior")
    )
    metrics.to_csv(args.metrics, index=False)

    out = combo[["ko_genes", "n_cells", "all_genes_seen_in_single"]].copy()
    for i, col in enumerate(target_cols):
        out[f"true_{col}"] = y[:, i]
        out[f"additive_pred_{col}"] = additive_pred[:, i]
        out[f"ridge_pred_{col}"] = ridge_pred[:, i]
    out.to_csv(args.predictions, index=False)
    print(f"Saved delta table to {args.delta}")
    print(f"Saved combo metrics to {args.metrics}")
    print(f"Saved combo predictions to {args.predictions}")


if __name__ == "__main__":
    main()
