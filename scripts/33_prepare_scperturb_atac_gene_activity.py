from __future__ import annotations

from pathlib import Path
import gzip

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "scperturb_atac" / "Liscovitch-BrauerSanjana2021_K562_1" / "gene_scores"
OUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity.h5ad"


def main() -> None:
    obs = pd.read_csv(BASE / "obs.csv")
    var = pd.read_csv(BASE / "var.csv", index_col=0)
    with gzip.open(BASE / "counts.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    x = matrix.astype(np.float32)
    obs = obs.copy()
    obs["ko_target"] = obs["perturbation"].fillna("unassigned").astype(str)
    keep = obs["ko_target"].ne("unassigned").to_numpy()
    obs = obs.loc[keep].reset_index(drop=True)
    x = x[keep]
    adata = ad.AnnData(X=x, obs=obs)
    adata.obs_names = obs["cell_barcode"].astype(str).to_numpy()
    adata.var_names = [str(x).upper() for x in var.index]
    adata.var_names_make_unique()
    adata.uns["source"] = "scPerturb ATAC Liscovitch-Brauer/Sanjana 2021 K562 gene activity"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")
    counts = obs["ko_target"].value_counts().rename_axis("ko_target").reset_index(name="n_cells")
    counts.to_csv(OUT.with_name("liscovitch_k562_gene_activity_ko_counts.csv"), index=False)
    print(adata)
    print(f"Saved {OUT}")
    print(counts.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
