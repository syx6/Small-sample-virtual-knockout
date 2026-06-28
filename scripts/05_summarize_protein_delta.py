from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main() -> None:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    if "protein" not in adata.obsm:
        raise ValueError("No protein matrix found in adata.obsm['protein'].")

    protein = np.asarray(adata.obsm["protein"])
    protein_names = list(adata.uns["protein_names"])
    obs = adata.obs.copy()
    text = obs["ko_target"].astype(str).str.lower()
    control = (
        text.str.contains("nt")
        | text.str.contains("control")
        | text.str.contains("non")
        | text.str.contains("safe")
        | text.str.contains("neg")
    )
    base = protein[control].mean(axis=0)

    rows = []
    for target, idx in obs.groupby("ko_target", observed=True).indices.items():
        vals = protein[list(idx)].mean(axis=0) - base
        row = {
            "ko_target": target,
            "n_cells": len(idx),
            "is_control": bool(str(target).lower().find("nt") >= 0),
        }
        row.update({f"delta_protein_{name}": vals[i] for i, name in enumerate(protein_names)})
        rows.append(row)

    out = pd.DataFrame(rows).sort_values("ko_target")
    Path("results").mkdir(exist_ok=True)
    out.to_csv("results/papalexi_protein_delta.csv", index=False)
    print("Saved protein deltas to results/papalexi_protein_delta.csv")


if __name__ == "__main__":
    main()
