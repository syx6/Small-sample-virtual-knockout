from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import requests
from mudata import read_h5mu
from scipy import sparse


PAPALEXI_URL = "https://figshare.com/ndownloader/files/36509460"


PATHWAY_GENESETS = {
    "IFNG_JAK_STAT": [
        "IFNGR1",
        "IFNGR2",
        "STAT1",
        "STAT2",
        "STAT3",
        "JAK1",
        "JAK2",
        "IRF1",
        "IRF7",
        "CXCL10",
        "CXCL11",
        "GBP1",
        "GBP2",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
    ],
    "ANTIGEN_PRESENTATION": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "PSMB8",
        "PSMB9",
        "NLRC5",
    ],
    "NRF2_STRESS": [
        "NFE2L2",
        "KEAP1",
        "NQO1",
        "HMOX1",
        "GCLC",
        "GCLM",
        "TXNRD1",
        "SLC7A11",
    ],
    "APOPTOSIS": [
        "BAX",
        "BAK1",
        "BCL2",
        "BCL2L1",
        "CASP3",
        "CASP7",
        "CASP8",
        "FAS",
        "TNFRSF10B",
    ],
    "CELL_CYCLE_G2M": [
        "MKI67",
        "TOP2A",
        "CDK1",
        "CCNB1",
        "CCNB2",
        "AURKA",
        "AURKB",
        "UBE2C",
    ],
    "MYELOID_INFLAMMATION": [
        "IL1B",
        "TNF",
        "NFKBIA",
        "NFKB1",
        "CXCL8",
        "CCL2",
        "ICAM1",
        "CD83",
    ],
    "IMMUNE_CHECKPOINT": [
        "CD274",
        "PDCD1LG2",
        "CD86",
        "CD80",
        "LGALS9",
        "HAVCR2",
    ],
}


def pick_obs_column(obs: pd.DataFrame, candidates: list[str]) -> str:
    for name in candidates:
        if name in obs.columns:
            return name
    raise KeyError(f"None of these columns are present: {candidates}")


def download_file(url: str, output: Path) -> None:
    if output.exists() and output.stat().st_size > 0:
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with output.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)


def normalize_log1p(matrix):
    counts = matrix.astype(np.float64)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    totals[totals == 0] = 1.0
    if sparse.issparse(counts):
        normalized = counts.multiply((1e4 / totals)[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        return normalized
    return np.log1p(counts * (1e4 / totals)[:, None])


def score_genes_z(matrix, var_names: pd.Index, genes: list[str]) -> np.ndarray | None:
    present = [gene for gene in genes if gene in var_names]
    if len(present) < 3:
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
    parser.add_argument("--cells-per-perturbation", type=int, default=120)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--raw", default="data/papalexi_2021.h5mu")
    parser.add_argument("--output", default="data/papalexi_small_pathway.h5ad")
    args = parser.parse_args()

    np.random.seed(args.seed)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    download_file(PAPALEXI_URL, Path(args.raw))
    mdata = read_h5mu(args.raw)
    mdata.pull_obs()
    mdata.pull_var()

    rna = mdata["rna"].copy()
    obs = mdata.obs.copy()
    perturb_col = pick_obs_column(obs, ["gene_target", "perturbation", "guide_identity", "gene", "rna:gene"])

    shared_cells = rna.obs_names.intersection(obs.index)
    rna = rna[shared_cells].copy()
    rna.obs = rna.obs.join(obs.loc[shared_cells], how="left", rsuffix="_mdata")
    rna.obs["ko_target"] = rna.obs[perturb_col].astype(str)

    chosen = []
    for target, index in rna.obs.groupby("ko_target", observed=True).indices.items():
        index = np.array(index)
        n = min(args.cells_per_perturbation, len(index))
        chosen.extend(np.random.choice(index, size=n, replace=False))
    rna = rna[np.array(chosen)].copy()

    protein_mod = "adt" if "adt" in mdata.mod else "protein" if "protein" in mdata.mod else None
    if protein_mod is not None:
        adt = mdata[protein_mod]
        common = rna.obs_names.intersection(adt.obs_names)
        rna = rna[common].copy()
        adt = adt[common].copy()
        protein = pd.DataFrame(
            adt.X.toarray() if hasattr(adt.X, "toarray") else adt.X,
            index=adt.obs_names,
            columns=adt.var_names,
        )
        rna.obsm["protein"] = protein.loc[rna.obs_names].to_numpy()
        rna.uns["protein_names"] = protein.columns.astype(str).tolist()

    rna.X = normalize_log1p(rna.X)

    for pathway, genes in PATHWAY_GENESETS.items():
        score = score_genes_z(rna.X, rna.var_names, genes)
        if score is not None:
            rna.obs[f"pathway_{pathway}"] = score

    pathway_cols = [col for col in rna.obs.columns if col.startswith("pathway_")]
    rna.uns["pathway_columns"] = pathway_cols
    rna = ad.AnnData(X=rna.X, obs=rna.obs.copy(), var=rna.var.copy(), obsm=rna.obsm.copy(), uns=rna.uns.copy())
    rna.write_h5ad(args.output)

    summary = rna.obs["ko_target"].value_counts().rename_axis("ko_target").reset_index(name="n_cells")
    summary.to_csv("results/papalexi_small_summary.csv", index=False)
    print(f"Saved {rna.n_obs} cells x {rna.n_vars} genes to {args.output}")
    print(f"Detected pathway columns: {pathway_cols}")


if __name__ == "__main__":
    main()
