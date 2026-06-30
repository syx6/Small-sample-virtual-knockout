from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


PEAK_RE = re.compile(r"^(?P<chrom>chr[^:_-]+|[0-9XYM]+)[:_-](?P<start>[0-9]+)[-_](?P<end>[0-9]+)$", re.IGNORECASE)


def _obsm_names(adata, obsm_key: str) -> list[str]:
    for key in [f"{obsm_key}_names", "peak_names", "atac_names"]:
        if key in adata.uns:
            return [str(x) for x in adata.uns[key]]
    if obsm_key not in adata.obsm:
        raise ValueError(f"obsm key '{obsm_key}' was not found. Available keys: {list(adata.obsm.keys())}")
    return [f"{obsm_key}_{i + 1}" for i in range(adata.obsm[obsm_key].shape[1])]


def _load_peak_names(input_h5ad: str | Path | None, obsm_key: str, feature_names_csv: str | Path | None) -> list[str]:
    if feature_names_csv:
        table = pd.read_csv(feature_names_csv)
        for col in ["feature_name", "peak", "name", "raw_feature_name"]:
            if col in table.columns:
                return table[col].astype(str).tolist()
        raise ValueError("Feature-name CSV must contain feature_name, peak, name, or raw_feature_name.")
    if not input_h5ad:
        raise ValueError("Provide --input-h5ad or --feature-names-csv.")
    import anndata as ad

    adata = ad.read_h5ad(input_h5ad, backed="r")
    return _obsm_names(adata, obsm_key)


def _parse_peak(name: str) -> tuple[str | None, float, float, float]:
    text = str(name).replace(",", "")
    match = PEAK_RE.match(text)
    if not match:
        return None, np.nan, np.nan, np.nan
    chrom = match.group("chrom")
    if not chrom.lower().startswith("chr"):
        chrom = f"chr{chrom}"
    start = float(match.group("start"))
    end = float(match.group("end"))
    return chrom, start, end, (start + end) / 2.0


def _standardize_gene_tss(path: str | Path | None) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    table = pd.read_csv(path)
    lower = {col.lower(): col for col in table.columns}
    gene_col = lower.get("gene") or lower.get("gene_name") or lower.get("symbol") or lower.get("target_gene")
    chrom_col = lower.get("chrom") or lower.get("chr") or lower.get("chromosome")
    tss_col = lower.get("tss") or lower.get("tss_position")
    start_col = lower.get("start")
    end_col = lower.get("end")
    if not gene_col or not chrom_col or (not tss_col and not start_col):
        raise ValueError("Gene TSS CSV needs gene/gene_name, chrom, and tss or start/end columns.")
    out = pd.DataFrame(
        {
            "gene": table[gene_col].astype(str).str.upper(),
            "chrom": table[chrom_col].astype(str).map(lambda x: x if x.lower().startswith("chr") else f"chr{x}"),
        }
    )
    if tss_col:
        out["tss"] = pd.to_numeric(table[tss_col], errors="coerce")
    else:
        start = pd.to_numeric(table[start_col], errors="coerce")
        end = pd.to_numeric(table[end_col], errors="coerce") if end_col else start
        out["tss"] = (start + end) / 2.0
    return out.dropna(subset=["tss"]).reset_index(drop=True)


def _target_set(value: str | None) -> set[str]:
    if not value:
        return set()
    return {part.strip().upper() for part in re.split(r"[,;+|]", value) if part.strip()}


def _nearest_gene_scores(peaks: pd.DataFrame, gene_tss: pd.DataFrame, target_genes: set[str], max_distance: int) -> pd.DataFrame:
    if gene_tss.empty:
        peaks["target_gene"] = ""
        peaks["nearest_gene"] = ""
        peaks["distance_to_tss"] = np.nan
        peaks["peak_gene_link_score"] = 0.0
        peaks["locus_score"] = 0.0
        return peaks
    rows = []
    for _, peak in peaks.iterrows():
        chrom = peak["chrom"]
        center = peak["center"]
        if not isinstance(chrom, str) or pd.isna(center):
            rows.append(("", "", np.nan, 0.0, 0.0))
            continue
        genes = gene_tss.loc[gene_tss["chrom"] == chrom].copy()
        if genes.empty:
            rows.append(("", "", np.nan, 0.0, 0.0))
            continue
        genes["distance"] = (genes["tss"] - center).abs()
        nearest = genes.sort_values("distance").iloc[0]
        distance = float(nearest["distance"])
        nearest_gene = str(nearest["gene"])
        target_gene = nearest_gene if not target_genes or nearest_gene in target_genes else ""
        distance_score = float(np.exp(-distance / max(1.0, max_distance)))
        promoter_bonus = 1.0 if distance <= 2000 else 0.65 if distance <= 20000 else 0.35
        target_bonus = 1.15 if nearest_gene in target_genes else 1.0
        link_score = min(1.0, distance_score * promoter_bonus * target_bonus)
        locus_score = min(1.0, promoter_bonus * target_bonus)
        rows.append((target_gene, nearest_gene, distance, link_score, locus_score))
    peaks[["target_gene", "nearest_gene", "distance_to_tss", "peak_gene_link_score", "locus_score"]] = pd.DataFrame(rows, index=peaks.index)
    return peaks


def _score_by_peak_name(path: str | Path | None, names: list[str], target_genes: set[str], default_score_col: str = "score") -> np.ndarray:
    if not path:
        return np.zeros(len(names), dtype=float)
    table = pd.read_csv(path)
    name_col = next((col for col in ["feature_name", "peak", "name", "raw_feature_name"] if col in table.columns), None)
    if not name_col:
        raise ValueError(f"{path} must contain feature_name, peak, name, or raw_feature_name.")
    score_col = next((col for col in [default_score_col, "weight", "prior_score", "regulatory_prior_score"] if col in table.columns), None)
    tf_col = next((col for col in ["tf", "motif", "gene", "target_gene"] if col in table.columns), None)
    table = table.copy()
    if score_col:
        table["_score"] = pd.to_numeric(table[score_col], errors="coerce").fillna(1.0)
    else:
        table["_score"] = 1.0
    if target_genes and tf_col:
        table["_score"] = np.where(table[tf_col].astype(str).str.upper().isin(target_genes), table["_score"] * 1.25, table["_score"])
    best = table.groupby(table[name_col].astype(str))["_score"].max().to_dict()
    raw = np.asarray([float(best.get(name, 0.0)) for name in names], dtype=float)
    if np.nanmax(raw) > 0:
        raw = raw / np.nanmax(raw)
    return raw


def annotate_peaks(
    input_h5ad: str | Path | None,
    obsm_key: str,
    out_csv: str | Path,
    feature_names_csv: str | Path | None = None,
    gene_tss_csv: str | Path | None = None,
    motif_hits_csv: str | Path | None = None,
    marker_peaks_csv: str | Path | None = None,
    target_genes: str | None = None,
    max_distance: int = 250_000,
) -> pd.DataFrame:
    names = _load_peak_names(input_h5ad, obsm_key, feature_names_csv)
    target = _target_set(target_genes)
    parsed = [_parse_peak(name) for name in names]
    peaks = pd.DataFrame(
        {
            "obsm_key": obsm_key,
            "feature_index": np.arange(len(names), dtype=int),
            "feature_name": names,
            "peak": names,
            "chrom": [item[0] for item in parsed],
            "start": [item[1] for item in parsed],
            "end": [item[2] for item in parsed],
            "center": [item[3] for item in parsed],
        }
    )
    gene_tss = _standardize_gene_tss(gene_tss_csv)
    peaks = _nearest_gene_scores(peaks, gene_tss, target, max_distance=max_distance)
    peaks["motif_to_peak_score"] = _score_by_peak_name(motif_hits_csv, names, target, default_score_col="motif_score")
    peaks["marker_score"] = _score_by_peak_name(marker_peaks_csv, names, target, default_score_col="marker_score")
    peaks["regulatory_prior_score"] = (
        0.40 * peaks["peak_gene_link_score"].astype(float)
        + 0.25 * peaks["motif_to_peak_score"].astype(float)
        + 0.20 * peaks["marker_score"].astype(float)
        + 0.15 * peaks["locus_score"].astype(float)
    ).clip(0.0, 1.0)
    peaks = peaks.drop(columns=["center"])
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    peaks.to_csv(out_csv, index=False)
    _write_report(peaks, out_csv)
    return peaks


def _write_report(peaks: pd.DataFrame, out_csv: Path) -> None:
    top = peaks.sort_values("regulatory_prior_score", ascending=False).head(12)
    text = f"""# Peak Annotation Report

Generated feature annotation file:

```text
{out_csv}
```

## Summary

- Peaks/features annotated: {len(peaks)}
- Peaks with parsed genomic coordinates: {int(peaks['chrom'].notna().sum())}
- Peaks with nonzero peak-gene score: {int((peaks['peak_gene_link_score'] > 0).sum())}
- Peaks with nonzero motif score: {int((peaks['motif_to_peak_score'] > 0).sum())}
- Peaks with nonzero marker score: {int((peaks['marker_score'] > 0).sum())}

## Top Regulatory Peaks

{top[['feature_name', 'target_gene', 'nearest_gene', 'distance_to_tss', 'regulatory_prior_score']].round(4).to_string(index=False)}

## How To Use

Pass this CSV to `run`, `score`, or `train-reference`:

```bash
python -m vkx.cli train-reference \\
  --input-h5ad perturb_multiome.h5ad \\
  --extra-obsm peak:peak,chromvar:tf \\
  --extra-feature-selection atac_peak \\
  --extra-feature-metadata-csv {out_csv}
```
"""
    out_csv.with_suffix(".report.md").write_text(text, encoding="utf-8")
