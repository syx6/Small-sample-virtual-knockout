from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def signed_score(values: pd.Series, direction: str) -> np.ndarray:
    arr = values.to_numpy(dtype=float)
    if direction == "decrease":
        return -arr
    if direction == "increase":
        return arr
    if direction == "absolute":
        return np.abs(arr)
    raise ValueError(f"Unknown direction: {direction}")


def eval_target(df: pd.DataFrame, target: str, pred_prefix: str, threshold: float, direction: str) -> dict:
    true_col = f"true_{target}"
    pred_col = f"{pred_prefix}_{target}"
    y_true_score = signed_score(df[true_col], direction)
    y_pred_score = signed_score(df[pred_col], direction)
    label = y_true_score >= threshold
    n_pos = int(label.sum())
    n_neg = int((~label).sum())
    row = {
        "model": pred_prefix,
        "target": target,
        "direction": direction,
        "threshold": threshold,
        "n_positive": n_pos,
        "n_negative": n_neg,
        "roc_auc": np.nan,
        "pr_auc": np.nan,
    }
    if n_pos > 0 and n_neg > 0:
        row["roc_auc"] = roc_auc_score(label, y_pred_score)
        row["pr_auc"] = average_precision_score(label, y_pred_score)
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/papalexi_multimodal_joint_predictions.csv")
    parser.add_argument("--output", default="results/papalexi_multimodal_auc.csv")
    parser.add_argument("--threshold", type=float, default=0.15)
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    true_targets = [col.removeprefix("true_") for col in df.columns if col.startswith("true_delta_")]
    pred_prefixes = sorted(
        {
            col.split("_delta_", 1)[0]
            for col in df.columns
            if col.endswith(tuple(true_targets)) and not col.startswith("true_")
        }
    )
    if not pred_prefixes:
        pred_prefixes = ["ridge_pred", "pls_pred"]

    rows = []
    for target in true_targets:
        directions = ["absolute"]
        if target in {"delta_pathway_IFNG_JAK_STAT", "delta_protein_PDL1"}:
            directions.extend(["decrease", "increase"])
        for pred_prefix in pred_prefixes:
            if f"{pred_prefix}_{target}" not in df.columns:
                continue
            for direction in directions:
                rows.append(eval_target(df, target, pred_prefix, args.threshold, direction))

    result = pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Saved AUC evaluation to {args.output}")


if __name__ == "__main__":
    main()
