from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


def main() -> None:
    adata = ad.read_h5ad("data/norman_small_program.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["cell_id"] = adata.obs_names
    frame["ko_target"] = obs["ko_genes"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("program_")]:
        frame[col] = obs[col].astype(float).values
    out = Path("data/examples/norman_gene_program_state.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    print(f"Saved {out} with {len(frame)} cells and {len(frame.columns) - 2} program features.")


if __name__ == "__main__":
    main()
