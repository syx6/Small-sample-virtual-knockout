from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vkx.core import control_mask


ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "results" / "hmpcite_multimodal_doubleko_cebp_med12" / "derived_state_scores.csv"
OUT = ROOT / "results" / "hmpcite_multimodal_doubleko_state_delta.csv"


def n_ko(label: str) -> int:
    if str(label).lower() == "control":
        return 0
    return str(label).count("+") + 1


def main() -> None:
    frame = pd.read_csv(STATE)
    features = [
        col
        for col in frame.columns
        if col not in {"cell_id", "ko_target"} and pd.api.types.is_numeric_dtype(frame[col])
    ]
    ctrl = frame.loc[control_mask(frame["ko_target"]), features]
    if ctrl.empty:
        raise ValueError("No control cells found.")
    control_mean = ctrl.mean()
    rows = []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if len(group) < 20:
            continue
        delta = group[features].mean() - control_mean
        row = {
            "ko_genes": ko,
            "ko_target": ko,
            "n_cells": len(group),
            "n_ko_genes": n_ko(ko),
        }
        for feature in features:
            row[f"delta_{feature}"] = float(delta[feature])
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"Saved {OUT}")
    print(out["n_ko_genes"].value_counts().sort_index().to_string())
    print(out.loc[out["n_ko_genes"].eq(2), ["ko_genes", "n_cells"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
