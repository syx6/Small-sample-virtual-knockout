from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd


def infer_control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return (
        text.str.contains("nt")
        | text.str.contains("non")
        | text.str.contains("control")
        | text.str.contains("safe")
        | text.str.contains("neg")
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/papalexi_small_pathway.h5ad")
    parser.add_argument("--output", default="results/papalexi_pathway_delta.csv")
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input)
    pathway_cols = list(adata.uns.get("pathway_columns", []))
    if not pathway_cols:
        pathway_cols = [col for col in adata.obs.columns if col.startswith("pathway_")]
    if not pathway_cols:
        raise ValueError("No pathway score columns found.")

    control_mask = infer_control_mask(adata.obs["ko_target"])
    if control_mask.sum() == 0:
        raise ValueError("Could not infer negative-control cells from ko_target.")

    control_mean = adata.obs.loc[control_mask, pathway_cols].mean()
    rows = []
    for target, frame in adata.obs.groupby("ko_target", observed=True):
        mean_scores = frame[pathway_cols].mean()
        delta = mean_scores - control_mean
        row = {"ko_target": target, "n_cells": len(frame), "is_control": bool(infer_control_mask(pd.Series([target])).iloc[0])}
        row.update({col: mean_scores[col] for col in pathway_cols})
        row.update({f"delta_{col}": delta[col] for col in pathway_cols})
        rows.append(row)

    result = pd.DataFrame(rows).sort_values(["is_control", "ko_target"])
    result.to_csv(args.output, index=False)
    print(f"Saved pathway deltas to {args.output}")


if __name__ == "__main__":
    main()
