from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import MultiLabelBinarizer, StandardScaler


def split_targets(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    parts = [part.strip() for part in text.split("_") if part.strip()]
    return parts or [str(label)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/papalexi_pathway_delta.csv")
    parser.add_argument("--predictions", default="results/papalexi_baseline_predictions.csv")
    parser.add_argument("--metrics", default="results/papalexi_baseline_metrics.csv")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    df = df.loc[~df["is_control"].astype(bool)].copy()
    target_cols = [col for col in df.columns if col.startswith("delta_pathway_")]
    if len(df) < 3:
        raise ValueError("Need at least 3 non-control perturbations for leave-one-out baseline.")

    labels = [split_targets(label) for label in df["ko_target"]]
    mlb = MultiLabelBinarizer()
    x = mlb.fit_transform(labels)
    y = df[target_cols].to_numpy(dtype=float)

    loo = LeaveOneOut()
    pred = np.zeros_like(y)
    for train_idx, test_idx in loo.split(x):
        x_scaler = StandardScaler(with_mean=False)
        y_scaler = StandardScaler()
        x_train = x_scaler.fit_transform(x[train_idx])
        y_train = y_scaler.fit_transform(y[train_idx])
        model = Ridge(alpha=1.0)
        model.fit(x_train, y_train)
        y_hat = model.predict(x_scaler.transform(x[test_idx]))
        pred[test_idx] = y_scaler.inverse_transform(y_hat)

    pred_df = df[["ko_target", "n_cells"]].copy()
    for i, col in enumerate(target_cols):
        pred_df[f"true_{col}"] = y[:, i]
        pred_df[f"pred_{col}"] = pred[:, i]
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
