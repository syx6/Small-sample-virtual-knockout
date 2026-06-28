from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main() -> None:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["cell_id"] = adata.obs_names
    frame["ko_target"] = obs["ko_target"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("pathway_")]:
        frame[col] = obs[col].astype(float).values
    protein = np.asarray(adata.obsm["protein"])
    for i, name in enumerate(adata.uns["protein_names"]):
        frame[f"protein_{name}"] = protein[:, i]
    out = Path("data/examples/papalexi_pathway_protein_state.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    print(f"Saved {out} with {len(frame)} cells and {len(frame.columns) - 2} state features.")


if __name__ == "__main__":
    main()
