from __future__ import annotations

from pathlib import Path
import gzip

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "scperturb_atac" / "Liscovitch-BrauerSanjana2021_K562_1"
GENE_ACTIVITY = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity.h5ad"
OUT = ROOT / "data" / "scperturb_atac" / "liscovitch_k562_gene_activity_chromvar.h5ad"


def clean_motif_name(name: str) -> str:
    text = str(name)
    if "#" in text:
        text = text.split("#", 1)[0]
    return text.strip()


def main() -> None:
    adata = ad.read_h5ad(GENE_ACTIVITY)
    chrom_obs = pd.read_csv(BASE / "ChromVar" / "obs.csv")
    chrom_var = pd.read_csv(BASE / "ChromVar" / "var.csv", index_col=0)
    with gzip.open(BASE / "ChromVar" / "counts.mtx.gz", "rb") as handle:
        chrom = mmread(handle)
    if sparse.issparse(chrom):
        chrom = chrom.tocsr().astype(np.float32)
    else:
        chrom = np.asarray(chrom, dtype=np.float32)

    chrom_obs_names = chrom_obs["cell_barcode"].astype(str).to_numpy()
    pos = pd.Series(np.arange(len(chrom_obs_names)), index=chrom_obs_names)
    keep = pos.reindex(adata.obs_names.astype(str)).to_numpy()
    if np.isnan(keep).any():
        missing = int(np.isnan(keep).sum())
        raise ValueError(f"{missing} gene-activity cells were not found in ChromVar obs.csv.")
    chrom = chrom[keep.astype(int)]
    if sparse.issparse(chrom):
        chrom = chrom.toarray()

    names = [clean_motif_name(idx) for idx in chrom_var.index]
    adata.obsm["chromvar"] = np.asarray(chrom, dtype=np.float32)
    adata.uns["chromvar_names"] = names
    adata.uns["chromvar_source"] = "ChromVar motif activity from scPerturb ATAC Liscovitch-Brauer/Sanjana 2021"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")
    print(adata)
    print(f"Saved {OUT}")
    print(f"ChromVar motif features: {len(names)}")


if __name__ == "__main__":
    main()
