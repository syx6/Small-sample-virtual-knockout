from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


GENESETS = {
    "ERYTHROID": [
        "KLF1",
        "GATA1",
        "HBB",
        "HBA1",
        "HBA2",
        "ALAS2",
        "GYPA",
        "EPB42",
        "SLC4A1",
    ],
    "GRANULOCYTE_APOPTOSIS": [
        "CEBPE",
        "CEBPB",
        "SPI1",
        "ELANE",
        "MPO",
        "BAX",
        "BAK1",
        "CASP3",
        "CASP8",
    ],
    "MAPK_TGFB": [
        "MAPK1",
        "TGFBR2",
        "SMAD4",
        "DUSP1",
        "JUN",
        "FOS",
        "EGR1",
    ],
    "PRO_GROWTH": [
        "MYC",
        "CDK1",
        "MKI67",
        "TOP2A",
        "PCNA",
        "CCNB1",
        "CCND1",
    ],
    "MEGAKARYOCYTE": [
        "GATA2",
        "TAL1",
        "RUNX1",
        "ITGA2B",
        "PF4",
        "PPBP",
    ],
    "PIONEER_TF": [
        "FOXA1",
        "FOXA3",
        "AHR",
        "CEBPA",
        "CEBPB",
        "SPI1",
        "KLF1",
    ],
}


def split_perturbation(label: str) -> list[str]:
    parts = [part.strip() for part in str(label).split("+") if part.strip()]
    return [part for part in parts if part.lower() != "ctrl"]


def normalize_log1p(matrix):
    counts = matrix.astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    if sparse.issparse(counts):
        normalized = counts.multiply((1e4 / totals)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        return normalized
    return np.log1p(counts * (1e4 / totals)[:, None])


def score_genes(matrix, var_names: pd.Index, genes: list[str]) -> np.ndarray | None:
    present = [gene for gene in genes if gene in var_names]
    if len(present) < 2:
        return None
    idx = var_names.get_indexer(present)
    subset = matrix[:, idx]
    if sparse.issparse(subset):
        subset = subset.toarray()
    subset = np.asarray(subset, dtype=np.float64)
    mean = subset.mean(axis=0)
    std = subset.std(axis=0)
    std[std == 0] = 1.0
    return ((subset - mean) / std).mean(axis=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/norman_2019.h5ad")
    parser.add_argument("--output", default="data/norman_small_program.h5ad")
    parser.add_argument("--cells-per-perturbation", type=int, default=120)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    adata = ad.read_h5ad(args.input)
    adata.obs["ko_target"] = adata.obs["guide_merged"].astype(str)
    adata.obs["ko_genes"] = adata.obs["ko_target"].map(lambda x: "+".join(split_perturbation(x)) or "ctrl")
    adata.obs["n_ko_genes"] = adata.obs["ko_target"].map(lambda x: len(split_perturbation(x)))

    chosen = []
    for _, index in adata.obs.groupby("ko_genes", observed=True).indices.items():
        index = np.asarray(index)
        n = min(args.cells_per_perturbation, len(index))
        chosen.extend(rng.choice(index, size=n, replace=False))
    adata = adata[np.asarray(chosen)].copy()

    if "gene_name" in adata.var.columns:
        adata.var_names = adata.var["gene_name"].astype(str)
        adata.var_names_make_unique()

    adata.X = normalize_log1p(adata.X)
    for name, genes in GENESETS.items():
        score = score_genes(adata.X, adata.var_names, genes)
        if score is not None:
            adata.obs[f"program_{name}"] = score

    program_cols = [col for col in adata.obs.columns if col.startswith("program_")]
    adata.uns["program_columns"] = program_cols
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)

    summary = (
        adata.obs.groupby(["ko_genes", "n_ko_genes"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["n_ko_genes", "n_cells", "ko_genes"], ascending=[True, False, True])
    )
    summary.to_csv("results/norman_small_summary.csv", index=False)
    print(f"Saved {adata.n_obs} cells x {adata.n_vars} genes to {args.output}")
    print(f"Program columns: {program_cols}")


if __name__ == "__main__":
    main()
