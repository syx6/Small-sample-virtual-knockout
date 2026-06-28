from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
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


GENE_MODULES = {
    "receptor_ifn": {"IFNGR1", "IFNGR2"},
    "jak_stat_core": {"JAK1", "JAK2", "STAT1", "STAT2", "STAT3", "STAT5A"},
    "irf_tf": {"IRF1", "IRF7", "ETV7"},
    "checkpoint_ligand": {"CD274", "PDCD1LG2", "CD86"},
    "checkpoint_surface_regulator": {"CMTM6", "MARCH8"},
    "ubiquitin_proteostasis": {"UBE2L6", "CUL3"},
    "chromatin_transcription": {"BRD4", "MYC", "POU2F2", "ATF2"},
    "tgfb_smad": {"SMAD4"},
    "myeloid_tf": {"SPI1"},
    "membrane_context": {"CAV1", "TNFRSF14"},
}


PROTEIN_TARGETS = {
    "CD86": {"CD86"},
    "PDL1": {"CD274", "CMTM6", "JAK2", "STAT1", "IFNGR1", "IFNGR2", "IRF1"},
    "PDL2": {"PDCD1LG2"},
    "CD366": {"HAVCR2"},
}


def split_targets(label: str) -> set[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    return {part.strip().upper() for part in text.split("_") if part.strip()}


def build_features(ko_targets: pd.Series) -> tuple[np.ndarray, list[str]]:
    names: list[str] = []
    rows = []
    pathway_names = list(PATHWAY_GENESETS)
    module_names = list(GENE_MODULES)
    protein_names = list(PROTEIN_TARGETS)

    names.extend([f"pathway_prior_{name}" for name in pathway_names])
    names.extend([f"module_{name}" for name in module_names])
    names.extend([f"protein_prior_{name}" for name in protein_names])
    names.extend(["n_ko_genes", "is_known_ifn_axis"])

    ifn_axis = PATHWAY_GENESETS["IFNG_JAK_STAT"] | GENE_MODULES["receptor_ifn"] | GENE_MODULES["jak_stat_core"]
    for label in ko_targets:
        genes = split_targets(label)
        row = []
        row.extend(float(bool(genes & PATHWAY_GENESETS[name])) for name in pathway_names)
        row.extend(float(bool(genes & GENE_MODULES[name])) for name in module_names)
        row.extend(float(bool(genes & PROTEIN_TARGETS[name])) for name in protein_names)
        row.append(float(len(genes)))
        row.append(float(bool(genes & ifn_axis)))
        rows.append(row)
    return np.asarray(rows, dtype=float), names


def fit_predict_loo(x: np.ndarray, y: np.ndarray, model_name: str) -> np.ndarray:
    loo = LeaveOneOut()
    pred = np.zeros_like(y, dtype=float)
    for train_idx, test_idx in loo.split(x):
        x_scaler = StandardScaler()
        y_scaler = StandardScaler()
        x_train = x_scaler.fit_transform(x[train_idx])
        y_train = y_scaler.fit_transform(y[train_idx])

        if model_name == "ridge":
            model = Ridge(alpha=1.0)
        elif model_name == "pls":
            n_components = min(4, x_train.shape[0] - 1, x_train.shape[1], y_train.shape[1])
            model = PLSRegression(n_components=max(1, n_components), scale=False)
        else:
            raise ValueError(f"Unknown model: {model_name}")

        model.fit(x_train, y_train)
        pred[test_idx] = y_scaler.inverse_transform(model.predict(x_scaler.transform(x[test_idx])))
    return pred


def metric_table(y: np.ndarray, pred: np.ndarray, target_cols: list[str], model: str) -> pd.DataFrame:
    rows = []
    for i, col in enumerate(target_cols):
        rows.append(
            {
                "model": model,
                "target": col,
                "mae": mean_absolute_error(y[:, i], pred[:, i]),
                "r2": r2_score(y[:, i], pred[:, i]),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pathway", default="results/papalexi_pathway_delta.csv")
    parser.add_argument("--protein", default="results/papalexi_protein_delta.csv")
    parser.add_argument("--metrics", default="results/papalexi_multimodal_joint_metrics.csv")
    parser.add_argument("--predictions", default="results/papalexi_multimodal_joint_predictions.csv")
    parser.add_argument("--correlation", default="results/papalexi_pathway_protein_correlation.csv")
    args = parser.parse_args()

    pathway = pd.read_csv(args.pathway)
    protein = pd.read_csv(args.protein)
    df = pathway.merge(protein, on=["ko_target", "n_cells", "is_control"], how="inner")
    df = df.loc[~df["is_control"].astype(bool)].copy()

    pathway_cols = [col for col in df.columns if col.startswith("delta_pathway_")]
    protein_cols = [col for col in df.columns if col.startswith("delta_protein_")]
    target_cols = pathway_cols + protein_cols

    x, feature_names = build_features(df["ko_target"])
    y = df[target_cols].to_numpy(dtype=float)

    ridge_pred = fit_predict_loo(x, y, "ridge")
    pls_pred = fit_predict_loo(x, y, "pls")

    metrics = pd.concat(
        [
            metric_table(y, ridge_pred, target_cols, "ridge_prior_joint"),
            metric_table(y, pls_pred, target_cols, "pls_prior_joint"),
        ],
        ignore_index=True,
    )
    Path(args.metrics).parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.metrics, index=False)

    pred_df = df[["ko_target", "n_cells"]].copy()
    for i, col in enumerate(target_cols):
        pred_df[f"true_{col}"] = y[:, i]
        pred_df[f"ridge_pred_{col}"] = ridge_pred[:, i]
        pred_df[f"pls_pred_{col}"] = pls_pred[:, i]
    for i, name in enumerate(feature_names):
        pred_df[name] = x[:, i]
    pred_df.to_csv(args.predictions, index=False)

    corr = df[pathway_cols + protein_cols].corr().loc[pathway_cols, protein_cols]
    corr.reset_index(names="pathway").melt(
        id_vars="pathway", var_name="protein", value_name="pearson_r"
    ).to_csv(args.correlation, index=False)

    print(f"Saved metrics to {args.metrics}")
    print(f"Saved predictions to {args.predictions}")
    print(f"Saved pathway-protein correlations to {args.correlation}")


if __name__ == "__main__":
    main()
