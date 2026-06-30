from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def make_ko_summary_cards(
    delta_csv: str | Path,
    out_dir: str | Path,
    auc_csv: str | Path | None = None,
    confidence_csv: str | Path | None = None,
    max_features: int = 8,
) -> pd.DataFrame:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    delta = pd.read_csv(delta_csv)
    auc = pd.read_csv(auc_csv) if auc_csv and Path(auc_csv).exists() else pd.DataFrame()
    confidence = pd.read_csv(confidence_csv) if confidence_csv and Path(confidence_csv).exists() else pd.DataFrame()
    rows = []
    for _, row in delta.iterrows():
        ko = str(row["ko_target"])
        card = _card_table(row, max_features=max_features)
        card_csv = out / f"ko_card_{_safe_name(ko)}.csv"
        card.to_csv(card_csv, index=False)
        fig_path = _plot_card(ko, card, row, auc, confidence, out)
        rows.append(
            {
                "ko_target": ko,
                "card_csv": str(card_csv),
                "card_figure": str(fig_path) if fig_path else "",
                "top_predicted_increase": card.sort_values("pred_delta", ascending=False)["feature"].iloc[0] if not card.empty else "",
                "top_predicted_decrease": card.sort_values("pred_delta", ascending=True)["feature"].iloc[0] if not card.empty else "",
            }
        )
    index = pd.DataFrame(rows)
    index.to_csv(out / "ko_summary_cards_index.csv", index=False)
    _write_cards_report(index, out)
    return index


def _safe_name(value: str) -> str:
    return "".join(ch if ch.isalnum() else "_" for ch in value)[:80]


def _card_table(row: pd.Series, max_features: int) -> pd.DataFrame:
    features = [col.removeprefix("pred_delta_") for col in row.index if col.startswith("pred_delta_")]
    rows = []
    for feature in features:
        pred = float(row.get(f"pred_delta_{feature}", np.nan))
        true = row.get(f"true_delta_{feature}", np.nan)
        true = float(true) if pd.notna(true) else np.nan
        rows.append(
            {
                "feature": feature,
                "true_delta": true,
                "pred_delta": pred,
                "error": pred - true if pd.notna(true) else np.nan,
                "abs_pred_delta": abs(pred),
                "direction_match": np.sign(pred) == np.sign(true) if pd.notna(true) and abs(true) > 1e-9 else np.nan,
            }
        )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    return table.sort_values("abs_pred_delta", ascending=False).head(max_features).drop(columns=["abs_pred_delta"]).reset_index(drop=True)


def _plot_card(
    ko: str,
    card: pd.DataFrame,
    row: pd.Series,
    auc: pd.DataFrame,
    confidence: pd.DataFrame,
    out_dir: Path,
) -> Path | None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return None
    if card.empty:
        return None
    plot = card.copy()
    plot["feature_short"] = plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False).str.slice(0, 28)
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4.8, 0.45 * len(plot) + 2.2)), gridspec_kw={"width_ratios": [1.6, 1]}, constrained_layout=True)
    melted = plot.melt(id_vars=["feature_short"], value_vars=[col for col in ["true_delta", "pred_delta"] if plot[col].notna().any()], var_name="state", value_name="delta")
    sns.barplot(data=melted, x="delta", y="feature_short", hue="state", palette={"true_delta": "#2A9D8F", "pred_delta": "#E76F51"}, ax=axes[0])
    axes[0].axvline(0, color="0.2", linewidth=1)
    axes[0].set_title("Top predicted state changes")
    axes[0].set_xlabel("KO delta")
    axes[0].set_ylabel("")
    axes[1].axis("off")
    conf_text = "not available"
    if not confidence.empty and "ko_target" in confidence.columns:
        hit = confidence.loc[confidence["ko_target"].astype(str) == ko]
        if not hit.empty:
            conf_text = str(hit.iloc[0].get("confidence", "not available"))
    auc_text = "not available"
    if not auc.empty and "roc_auc" in auc.columns and auc["roc_auc"].notna().any():
        auc_text = f"{float(auc['roc_auc'].dropna().iloc[0]):.3f}"
    summary = [
        f"KO target: {ko}",
        f"Prediction source: {row.get('prediction_source', 'evaluation_delta')}",
        f"Direction cosine: {_fmt(row.get('direction_cosine', np.nan))}",
        f"Mean abs error: {_fmt(row.get('mean_abs_delta_error', np.nan))}",
        f"AUC: {auc_text}",
        f"Confidence: {conf_text}",
    ]
    axes[1].text(0.02, 0.95, "\n".join(summary), va="top", fontsize=12, linespacing=1.55)
    fig.suptitle(f"Virtual KO Summary Card: {ko}", fontsize=15)
    path = out_dir / f"ko_card_{_safe_name(ko)}.png"
    fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return path


def _fmt(value) -> str:
    return "n/a" if pd.isna(value) else f"{float(value):.3f}"


def _write_cards_report(index: pd.DataFrame, out_dir: Path) -> None:
    text = f"""# KO Summary Cards

This folder contains one concise card per KO target.

Each card answers:

- What gene or gene pair was virtually knocked out?
- Which pathway/protein/ATAC features changed the most?
- If true KO labels are available, how close was virtual KO to real KO?
- If transfer confidence is available, how reliable is the reference application?

## Card Index

{index.to_string(index=False) if not index.empty else 'No cards generated.'}
"""
    (out_dir / "ko_summary_cards_report.md").write_text(text, encoding="utf-8")
