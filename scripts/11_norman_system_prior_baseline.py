from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.metrics import average_precision_score, mean_absolute_error, r2_score, roc_auc_score
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")


def split_genes(label: str) -> list[str]:
    if str(label).lower() == "ctrl":
        return []
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
                tokens = term.split()
                if tokens:
                    hub = tokens[0].upper()
                    if GENE_RE.match(hub):
                        genes.add(hub)
            if genes:
                terms.append((term, genes))
    return terms


def select_terms(
    priors_dir: Path,
    perturb_genes: set[str],
    max_terms_per_library: int,
) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(priors_dir.glob("*.gmt")):
        include_term_gene = path.stem == "ppi_hub"
        terms = parse_gmt(path, include_term_gene=include_term_gene)
        scored = []
        for term, genes in terms:
            overlap = len(genes & perturb_genes)
            if overlap == 0:
                continue
            size = len(genes)
            if size < 5 or size > 800:
                continue
            score = (overlap, -size)
            scored.append((score, f"{path.stem}:{term}", genes))
        scored.sort(reverse=True, key=lambda x: x[0])
        for _, name, genes in scored[:max_terms_per_library]:
            selected.append((name, genes))
    return selected


def build_prior_matrix(labels: pd.Series, terms: list[tuple[str, set[str]]]) -> sparse.csr_matrix:
    rows = []
    cols = []
    data = []
    for i, label in enumerate(labels):
        genes = set(split_genes(label))
        if not genes:
            continue
        denom = max(1, len(genes))
        for j, (_, members) in enumerate(terms):
            overlap = len(genes & members)
            if overlap:
                rows.append(i)
                cols.append(j)
                data.append(overlap / denom)
    return sparse.csr_matrix((data, (rows, cols)), shape=(len(labels), len(terms)))


def build_feature_matrices(train: pd.DataFrame, test: pd.DataFrame, terms: list[tuple[str, set[str]]]):
    train_genes = [split_genes(label) for label in train["ko_genes"]]
    test_genes = [split_genes(label) for label in test["ko_genes"]]
    mlb = MultiLabelBinarizer()
    x_gene_train = sparse.csr_matrix(mlb.fit_transform(train_genes), dtype=float)
    x_gene_test = sparse.csr_matrix(mlb.transform(test_genes), dtype=float)
    x_prior_train = build_prior_matrix(train["ko_genes"], terms)
    x_prior_test = build_prior_matrix(test["ko_genes"], terms)
    train_n = sparse.csr_matrix(np.asarray([[len(g)] for g in train_genes], dtype=float))
    test_n = sparse.csr_matrix(np.asarray([[len(g)] for g in test_genes], dtype=float))
    x_train = sparse.hstack([x_gene_train, x_prior_train, train_n], format="csr")
    x_test = sparse.hstack([x_gene_test, x_prior_test, test_n], format="csr")
    names = [f"gene:{name}" for name in mlb.classes_] + [f"prior:{name}" for name, _ in terms] + ["n_ko_genes"]
    return x_train, x_test, names


def predict_ridge(train: pd.DataFrame, test: pd.DataFrame, target_cols: list[str], terms: list[tuple[str, set[str]]]) -> np.ndarray:
    x_train, x_test, _ = build_feature_matrices(train, test, terms)
    y_train = train[target_cols].to_numpy(dtype=float)
    x_scaler = StandardScaler(with_mean=False)
    y_scaler = StandardScaler()
    x_train = x_scaler.fit_transform(x_train)
    x_test = x_scaler.transform(x_test)
    y_train = y_scaler.fit_transform(y_train)
    model = Ridge(alpha=10.0)
    model.fit(x_train, y_train)
    return y_scaler.inverse_transform(model.predict(x_test))


def train_additive(single: pd.DataFrame, combo: pd.DataFrame, target_cols: list[str]) -> np.ndarray:
    gene_to_delta = {}
    for _, row in single.iterrows():
        genes = split_genes(row["ko_genes"])
        if len(genes) == 1:
            gene_to_delta[genes[0]] = row[target_cols].to_numpy(dtype=float)
    mean_single = single[target_cols].mean().to_numpy(dtype=float)
    preds = []
    for label in combo["ko_genes"]:
        vals = [gene_to_delta.get(gene, mean_single) for gene in split_genes(label)]
        preds.append(np.sum(vals, axis=0))
    return np.asarray(preds, dtype=float)


def metric_rows(combo: pd.DataFrame, pred: np.ndarray, cols: list[str], model: str) -> list[dict]:
    rows = []
    y_all = combo[cols].to_numpy(dtype=float)
    subsets = {
        "all_combos": np.ones(len(combo), dtype=bool),
        "all_genes_seen": combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
        "has_unseen_gene": ~combo["all_genes_seen_in_single"].to_numpy(dtype=bool),
    }
    for subset, mask in subsets.items():
        if mask.sum() < 3:
            continue
        y = y_all[mask]
        p = pred[mask]
        for i, col in enumerate(cols):
            labels = np.abs(y[:, i]) >= 0.15
            roc = np.nan
            pr = np.nan
            if labels.sum() > 0 and (~labels).sum() > 0:
                roc = roc_auc_score(labels, np.abs(p[:, i]))
                pr = average_precision_score(labels, np.abs(p[:, i]))
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


def summarize_term_hits(labels: pd.Series, terms: list[tuple[str, set[str]]], output: Path) -> None:
    hits = []
    for label in labels:
        genes = set(split_genes(label))
        if not genes:
            continue
        matched = []
        for name, members in terms:
            if genes & members:
                matched.append(name)
        hits.append({"ko_genes": label, "n_prior_terms_hit": len(matched), "top_terms": "; ".join(matched[:8])})
    pd.DataFrame(hits).to_csv(output, index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delta", default="results/norman_program_delta.csv")
    parser.add_argument("--priors-dir", default="data/priors")
    parser.add_argument("--max-terms-per-library", type=int, default=200)
    parser.add_argument("--metrics", default="results/norman_system_prior_metrics.csv")
    parser.add_argument("--predictions", default="results/norman_system_prior_predictions.csv")
    parser.add_argument("--term-hits", default="results/norman_system_prior_term_hits.csv")
    args = parser.parse_args()

    delta = pd.read_csv(args.delta)
    target_cols = [col for col in delta.columns if col.startswith("delta_program_")]
    single = delta.loc[delta["n_ko_genes"] == 1].copy()
    combo = delta.loc[delta["n_ko_genes"] == 2].copy()
    seen_single = {split_genes(label)[0] for label in single["ko_genes"] if len(split_genes(label)) == 1}
    combo["all_genes_seen_in_single"] = combo["ko_genes"].map(lambda label: all(g in seen_single for g in split_genes(label)))

    perturb_genes = {gene for label in delta["ko_genes"] for gene in split_genes(label)}
    terms = select_terms(Path(args.priors_dir), perturb_genes, args.max_terms_per_library)
    print(f"Selected {len(terms)} prior terms")

    train = pd.concat([delta.loc[delta["n_ko_genes"] == 0], single], ignore_index=True)
    additive_pred = train_additive(single, combo, target_cols)
    system_pred = predict_ridge(train, combo, target_cols, terms)

    metrics = pd.DataFrame(
        metric_rows(combo, additive_pred, target_cols, "single_gene_additive")
        + metric_rows(combo, system_pred, target_cols, "system_prior_ridge")
    )
    Path(args.metrics).parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.metrics, index=False)

    out = combo[["ko_genes", "n_cells", "all_genes_seen_in_single"]].copy()
    y = combo[target_cols].to_numpy(dtype=float)
    for i, col in enumerate(target_cols):
        out[f"true_{col}"] = y[:, i]
        out[f"additive_pred_{col}"] = additive_pred[:, i]
        out[f"system_pred_{col}"] = system_pred[:, i]
    out.to_csv(args.predictions, index=False)
    summarize_term_hits(combo["ko_genes"], terms, Path(args.term_hits))
    print(f"Saved metrics to {args.metrics}")
    print(f"Saved predictions to {args.predictions}")
    print(f"Saved term hits to {args.term_hits}")


if __name__ == "__main__":
    main()
