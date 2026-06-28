from __future__ import annotations

import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "hmpcite_gse243244"
OUT = BASE / "hmpcite_perturbation_rna_adt_doubleko.h5ad"


def read_lines(path: Path) -> list[str]:
    return [line.rstrip("\n").split("\t")[0] for line in path.open(encoding="utf-8")]


def guide_gene(feature: str) -> str | None:
    gene = re.split(r"-", str(feature).strip())[0]
    if gene.lower().startswith("random"):
        return None
    if not re.match(r"^[A-Za-z0-9.-]+$", gene):
        return None
    return gene[:1].upper() + gene[1:].lower()


def assign_ko_labels(gdo: sparse.spmatrix, guide_features: list[str], threshold: int = 10) -> pd.DataFrame:
    gdo = gdo.tocsc()
    labels = []
    for cell_idx in range(gdo.shape[1]):
        col = gdo[:, cell_idx]
        genes = set()
        for row_idx, value in zip(col.indices, col.data):
            if value >= threshold:
                gene = guide_gene(guide_features[row_idx])
                if gene:
                    genes.add(gene)
        genes = sorted(genes)
        labels.append(
            {
                "ko_target": "+".join(genes) if genes else "control",
                "n_ko_genes": len(genes),
                "ko_genes": "+".join(genes),
            }
        )
    return pd.DataFrame(labels)


def normalize_log1p_counts(x: sparse.spmatrix) -> sparse.csr_matrix:
    x = x.tocsr().astype(np.float32)
    totals = np.asarray(x.sum(axis=1)).reshape(-1)
    totals[totals <= 0] = 1.0
    scale = 1e4 / totals
    x = sparse.diags(scale).dot(x).tocsr()
    x.data = np.log1p(x.data)
    return x


def read_mtx(path: Path) -> sparse.spmatrix:
    with path.open("rb") as handle:
        return mmread(handle)


def main() -> None:
    rna_dir = BASE / "Pertubation_cDNA"
    adt_dir = BASE / "Pertubation_ADT"
    gdo_dir = BASE / "Pertubation_GDO"

    barcodes = read_lines(rna_dir / "barcodes.tsv")
    rna_features = read_lines(rna_dir / "features.tsv")
    adt_features = read_lines(adt_dir / "features.tsv")
    guide_features = read_lines(gdo_dir / "features.tsv")

    rna = read_mtx(rna_dir / "matrix.mtx").tocsr()
    adt = read_mtx(adt_dir / "matrix.mtx").tocsr()
    gdo = read_mtx(gdo_dir / "matrix.mtx").tocsr()
    if rna.shape[1] != len(barcodes) or adt.shape[1] != len(barcodes) or gdo.shape[1] != len(barcodes):
        raise ValueError("RNA, ADT and GDO barcode counts do not match.")

    obs = assign_ko_labels(gdo, guide_features, threshold=10)
    keep = obs["n_ko_genes"].isin([0, 1, 2]).to_numpy()
    obs = obs.loc[keep].reset_index(drop=True)
    kept_barcodes = np.asarray(barcodes, dtype=object)[keep]

    x = normalize_log1p_counts(rna[:, keep].T)
    protein = np.asarray(adt[:, keep].T.toarray(), dtype=np.float32)
    protein = np.log1p(protein)

    adata = ad.AnnData(X=x, obs=obs)
    adata.obs_names = [str(x) for x in kept_barcodes]
    adata.var_names = [str(gene).upper() for gene in rna_features]
    adata.var_names_make_unique()
    adata.obsm["protein"] = protein
    adata.uns["protein_names"] = [str(x) for x in adt_features]
    adata.uns["source"] = "GSE243244 HMPCITE-seq perturbation sample; RNA + ADT + GDO-derived KO labels"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")

    counts = obs["ko_target"].value_counts().rename_axis("ko_target").reset_index(name="n_cells")
    counts["n_ko_genes"] = counts["ko_target"].map(lambda x: 0 if x == "control" else str(x).count("+") + 1)
    counts.to_csv(BASE / "hmpcite_ko_counts_threshold10.csv", index=False)
    print(adata)
    print(f"Saved {OUT}")
    print(counts.groupby("n_ko_genes")["n_cells"].agg(["count", "sum"]).to_string())
    print(counts.loc[counts["n_ko_genes"].eq(2)].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
