from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler


PATHWAY_GENESETS = {
    "IFNG_JAK_STAT": {
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
    },
    "ANTIGEN_PRESENTATION": {
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "B2M",
        "TAP1",
        "TAP2",
        "PSMB8",
        "PSMB9",
        "NLRC5",
    },
    "NRF2_STRESS": {
        "NFE2L2",
        "KEAP1",
        "NQO1",
        "HMOX1",
        "GCLC",
        "GCLM",
        "TXNRD1",
        "SLC7A11",
        "CUL3",
    },
    "CELL_CYCLE_G2M": {
        "MKI67",
        "TOP2A",
        "CDK1",
        "CCNB1",
        "CCNB2",
        "AURKA",
        "AURKB",
        "UBE2C",
        "MYC",
    },
    "MYELOID_INFLAMMATION": {
        "IL1B",
        "TNF",
        "NFKBIA",
        "NFKB1",
        "CXCL8",
        "CCL2",
        "ICAM1",
        "CD83",
        "SPI1",
    },
    "IMMUNE_CHECKPOINT": {
        "CD274",
        "PDCD1LG2",
        "CD86",
        "CD80",
        "LGALS9",
        "HAVCR2",
        "CMTM6",
        "MARCH8",
    },
}


def build_features(ko_targets: pd.Series) -> tuple[np.ndarray, list[str]]:
    names = list(PATHWAY_GENESETS)
    rows = []
    for target in ko_targets.astype(str):
        genes = {part.strip().upper() for part in target.replace("+", "_").replace("|", "_").split("_")}
        row = []
        for pathway in names:
            row.append(float(bool(genes & PATHWAY_GENESETS[pathway])))
        row.append(float(len(genes)))
        rows.append(row)
    return np.asarray(rows, dtype=float), [f"prior_{name}" for name in names] + ["n_ko_genes"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/papalexi_pathway_delta.csv")
    parser.add_argument("--predictions", default="results/papalexi_gene_prior_predictions.csv")
    parser.add_argument("--metrics", default="results/papalexi_gene_prior_metrics.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df = df.loc[~df["is_control"].astype(bool)].copy()
    target_cols = [col for col in df.columns if col.startswith("delta_pathway_")]
    x, feature_names = build_features(df["ko_target"])
    y = df[target_cols].to_numpy(dtype=float)

    loo = LeaveOneOut()
    pred = np.zeros_like(y)
    for train_idx, test_idx in loo.split(x):
        x_scaler = StandardScaler()
        y_scaler = StandardScaler()
        x_train = x_scaler.fit_transform(x[train_idx])
        y_train = y_scaler.fit_transform(y[train_idx])
        model = Ridge(alpha=1.0)
        model.fit(x_train, y_train)
        pred[test_idx] = y_scaler.inverse_transform(model.predict(x_scaler.transform(x[test_idx])))

    pred_df = df[["ko_target", "n_cells"]].copy()
    for i, col in enumerate(target_cols):
        pred_df[f"true_{col}"] = y[:, i]
        pred_df[f"pred_{col}"] = pred[:, i]
    for i, name in enumerate(feature_names):
        pred_df[name] = x[:, i]
    Path(args.predictions).parent.mkdir(parents=True, exist_ok=True)
    pred_df.to_csv(args.predictions, index=False)

    metrics = []
    for i, col in enumerate(target_cols):
        metrics.append(
            {
                "target": col,
                "mae": mean_absolute_error(y[:, i], pred[:, i]),
                "r2": r2_score(y[:, i], pred[:, i]),
            }
        )
    pd.DataFrame(metrics).to_csv(args.metrics, index=False)
    print(f"Saved predictions to {args.predictions}")
    print(f"Saved metrics to {args.metrics}")


if __name__ == "__main__":
    main()
