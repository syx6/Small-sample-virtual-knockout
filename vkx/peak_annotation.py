from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd


PEAK_RE = re.compile(r"^(?P<chrom>chr[^:_-]+|[0-9XYM]+)[:_-](?P<start>[0-9]+)[-_](?P<end>[0-9]+)$", re.IGNORECASE)
GTF_ATTR_RE = re.compile(r'(?P<key>[A-Za-z0-9_.-]+)\s+(?P<quote>"?)(?P<value>[^";]+)(?P=quote)')


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


def _parse_gtf_attributes(value: str) -> dict[str, str]:
    attrs = {}
    for match in GTF_ATTR_RE.finditer(str(value)):
        attrs[match.group("key")] = match.group("value").strip()
    return attrs


def make_gene_tss_from_gtf(
    gtf: str | Path,
    out_csv: str | Path,
    feature_type: str = "gene",
    gene_name_attr: str = "gene_name",
    gene_id_attr: str = "gene_id",
) -> pd.DataFrame:
    rows = []
    gtf = Path(gtf)
    with gtf.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != feature_type:
                continue
            chrom, _, _, start, end, _, strand, _, attrs_text = parts
            chrom = chrom.lstrip("\ufeff")
            attrs = _parse_gtf_attributes(attrs_text)
            gene = attrs.get(gene_name_attr) or attrs.get(gene_id_attr)
            if not gene:
                continue
            start_i = int(start)
            end_i = int(end)
            tss = start_i if strand != "-" else end_i
            if not chrom.lower().startswith("chr"):
                chrom = f"chr{chrom}"
            rows.append(
                {
                    "gene": str(gene).upper(),
                    "gene_id": attrs.get(gene_id_attr, ""),
                    "chrom": chrom,
                    "tss": tss,
                    "strand": strand,
                    "start": start_i,
                    "end": end_i,
                }
            )
    table = pd.DataFrame(rows).drop_duplicates(subset=["gene", "chrom", "tss"]).reset_index(drop=True)
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(out_csv, index=False)
    _write_tss_report(table, out_csv, gtf)
    return table


def _write_tss_report(table: pd.DataFrame, out_csv: Path, source_gtf: Path) -> None:
    chroms = ", ".join(table["chrom"].dropna().astype(str).value_counts().head(8).index.tolist())
    text = f"""# Gene TSS Report

Generated gene TSS file:

```text
{out_csv}
```

Source annotation:

```text
{source_gtf}
```

## Summary

- Genes/TSS rows: {len(table)}
- Unique genes: {table['gene'].nunique() if not table.empty else 0}
- Top chromosomes: {chroms or 'not available'}

## Columns

- `gene`: upper-case gene symbol used by VKX.
- `gene_id`: original gene id when available.
- `chrom`: chromosome with `chr` prefix.
- `tss`: transcription start site.
- `strand`, `start`, `end`: original gene model coordinates.
"""
    out_csv.with_suffix(".report.md").write_text(text, encoding="utf-8")


def standardize_peak_score_table(
    input_csv: str | Path,
    out_csv: str | Path,
    table_type: str = "motif",
    peak_col: str | None = None,
    score_col: str | None = None,
    tf_col: str | None = None,
) -> pd.DataFrame:
    table = pd.read_csv(input_csv)
    peak_col = peak_col or next((col for col in ["peak", "feature_name", "name", "raw_feature_name", "region"] if col in table.columns), None)
    if not peak_col:
        raise ValueError("Input score table needs a peak/feature_name/name/raw_feature_name/region column.")
    score_col = score_col or next((col for col in ["score", "motif_score", "marker_score", "weight", "qvalue", "pvalue"] if col in table.columns), None)
    out = pd.DataFrame({"peak": table[peak_col].astype(str), "feature_name": table[peak_col].astype(str)})
    if tf_col and tf_col in table.columns:
        out["tf"] = table[tf_col].astype(str).str.upper()
    else:
        inferred_tf = next((col for col in ["tf", "motif", "gene", "target_gene"] if col in table.columns), None)
        if inferred_tf:
            out["tf"] = table[inferred_tf].astype(str).str.upper()
    if score_col:
        raw = pd.to_numeric(table[score_col], errors="coerce")
        if score_col.lower() in {"pvalue", "p_value", "qvalue", "q_value"}:
            score = -np.log10(raw.clip(lower=1e-300))
        else:
            score = raw
        score = score.fillna(0.0)
        max_score = float(score.max()) if len(score) else 0.0
        score = score / max_score if max_score > 0 else score
    else:
        score = pd.Series(np.ones(len(table), dtype=float))
    if table_type == "marker":
        out["marker_score"] = score
    else:
        out["motif_score"] = score
    out_csv = Path(out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_csv, index=False)
    return out


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


def build_peak_annotation_pipeline(
    out_csv: str | Path,
    input_h5ad: str | Path | None = None,
    obsm_key: str = "peak",
    feature_names_csv: str | Path | None = None,
    gtf: str | Path | None = None,
    gene_tss_csv: str | Path | None = None,
    raw_motif_hits_csv: str | Path | None = None,
    motif_hits_csv: str | Path | None = None,
    raw_marker_peaks_csv: str | Path | None = None,
    marker_peaks_csv: str | Path | None = None,
    target_genes: str | None = None,
    max_distance: int = 250_000,
) -> pd.DataFrame:
    out_csv = Path(out_csv)
    work_dir = out_csv.parent
    work_dir.mkdir(parents=True, exist_ok=True)
    if gene_tss_csv is None and gtf is not None:
        gene_tss_csv = work_dir / "gene_tss.csv"
        make_gene_tss_from_gtf(gtf, gene_tss_csv)
    if motif_hits_csv is None and raw_motif_hits_csv is not None:
        motif_hits_csv = work_dir / "motif_to_peak.standardized.csv"
        standardize_peak_score_table(raw_motif_hits_csv, motif_hits_csv, table_type="motif")
    if marker_peaks_csv is None and raw_marker_peaks_csv is not None:
        marker_peaks_csv = work_dir / "marker_peaks.standardized.csv"
        standardize_peak_score_table(raw_marker_peaks_csv, marker_peaks_csv, table_type="marker")
    return annotate_peaks(
        input_h5ad=input_h5ad,
        obsm_key=obsm_key,
        out_csv=out_csv,
        feature_names_csv=feature_names_csv,
        gene_tss_csv=gene_tss_csv,
        motif_hits_csv=motif_hits_csv,
        marker_peaks_csv=marker_peaks_csv,
        target_genes=target_genes,
        max_distance=max_distance,
    )


def _write_report(peaks: pd.DataFrame, out_csv: Path) -> None:
    top = peaks.sort_values("regulatory_prior_score", ascending=False).head(12)
    figure_path = _plot_peak_annotation_summary(peaks, out_csv)
    figure_text = f"\n![Peak annotation summary]({figure_path.name})\n" if figure_path is not None else "\nPlot was not generated because matplotlib is not available in this environment.\n"
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

## QC Figure

{figure_text}

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


def _plot_peak_annotation_summary(peaks: pd.DataFrame, out_csv: Path) -> Path | None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        return None
    if peaks.empty:
        return None
    plot = peaks.sort_values("regulatory_prior_score", ascending=False).head(14).copy()
    plot["label"] = plot["feature_name"].astype(str).str.replace("chr", "", regex=False)
    plot["label"] = plot["label"].where(plot["label"].str.len() <= 28, plot["label"].str.slice(0, 25) + "...")
    components = ["peak_gene_link_score", "motif_to_peak_score", "marker_score", "locus_score"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, max(5.2, 0.38 * len(plot) + 2.5)), constrained_layout=True)
    sns.barplot(data=plot, x="regulatory_prior_score", y="label", color="#4C78A8", ax=axes[0])
    axes[0].set_title("Top regulatory-prior peaks")
    axes[0].set_xlabel("Regulatory prior score")
    axes[0].set_ylabel("")
    heat = plot.set_index("label")[components]
    sns.heatmap(heat, cmap="YlGnBu", vmin=0, vmax=1, annot=True, fmt=".2f", cbar_kws={"label": "component score"}, ax=axes[1])
    axes[1].set_title("Why these peaks were selected")
    axes[1].set_xlabel("")
    axes[1].set_ylabel("")
    fig.suptitle("ATAC Peak Annotation QC", fontsize=15)
    figure_path = out_csv.with_suffix(".summary.png")
    fig.savefig(figure_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return figure_path
