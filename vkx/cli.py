from __future__ import annotations

import argparse
from pathlib import Path


def parse_features(value: str | None) -> list[str] | None:
    if not value:
        return None
    return [part.strip() for part in value.split(",") if part.strip()]


def parse_holdouts(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def target_labels(args: argparse.Namespace) -> list[str]:
    value = getattr(args, "target_kos", None) or getattr(args, "holdouts", None)
    if not value:
        raise ValueError("Provide --target-kos. Use one gene like STAT1 or one pair like STAT1+JAK2.")
    return parse_holdouts(value)


def _auto_make_cards(out_dir: Path, delta_name: str, auc_name: str | None = None, confidence_name: str | None = None) -> None:
    from .cards import make_ko_summary_cards

    delta_csv = out_dir / delta_name
    if not delta_csv.exists():
        return
    auc_csv = out_dir / auc_name if auc_name and (out_dir / auc_name).exists() else None
    confidence_csv = out_dir / confidence_name if confidence_name and (out_dir / confidence_name).exists() else None
    try:
        make_ko_summary_cards(
            delta_csv=delta_csv,
            out_dir=out_dir / "ko_cards",
            auc_csv=auc_csv,
            confidence_csv=confidence_csv,
        )
    except Exception as exc:
        print(f"  note: KO summary cards were skipped: {exc}")


def _auto_diagnose(
    out_dir: Path,
    delta_name: str,
    manifest_name: str | None = None,
    confidence_name: str | None = None,
) -> None:
    from .diagnostics import diagnose_virtual_ko_results

    delta_csv = out_dir / delta_name
    if not delta_csv.exists():
        return
    manifest_csv = out_dir / manifest_name if manifest_name and (out_dir / manifest_name).exists() else None
    confidence_csv = out_dir / confidence_name if confidence_name and (out_dir / confidence_name).exists() else None
    try:
        diagnose_virtual_ko_results(
            delta_csv=delta_csv,
            manifest_csv=manifest_csv,
            confidence_csv=confidence_csv,
            out_dir=out_dir / "diagnosis",
        )
    except Exception as exc:
        print(f"  note: failure diagnosis was skipped: {exc}")


def _auto_summarize(out_dir: Path) -> None:
    from .result_report import summarize_result_directory

    try:
        summarize_result_directory(out_dir, out_dir / "readable_result_report")
    except Exception as exc:
        print(f"  note: readable result report was skipped: {exc}")


def run_fit(args: argparse.Namespace) -> None:
    import pandas as pd

    from .core import run_virtual_ko
    from .visualization import make_all_plots

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.state_csv)
    result = run_virtual_ko(
        frame=frame,
        ko_col=args.ko_col,
        holdouts=target_labels(args),
        prior_dir=args.prior_dir,
        features=parse_features(args.features),
        dataset_name=args.dataset_name,
        modality=args.modality,
        representation=args.representation,
        calibrate=args.calibrate,
        shape_calibrate=args.shape_calibrate,
        max_cells_per_state=args.max_cells_per_state,
        seed=args.seed,
    )
    result.metrics.to_csv(out_dir / "metrics.csv", index=False)
    result.summary.to_csv(out_dir / "summary.csv", index=False)
    result.virtual_cells.to_csv(out_dir / "virtual_cells.csv", index=False)
    result.delta_table.to_csv(out_dir / "delta_table.csv", index=False)
    result.auc_points.to_csv(out_dir / "auc_points.csv", index=False)
    result.calibration.to_csv(out_dir / "calibration.csv", index=False)
    auc_summary = make_all_plots(result.summary, result.delta_table, result.virtual_cells, result.auc_points, out_dir)
    _auto_make_cards(out_dir, "delta_table.csv", auc_name="auc_summary.csv")
    _auto_diagnose(out_dir, "delta_table.csv")
    _auto_summarize(out_dir)
    print("Saved reusable virtual KO outputs:")
    print(f"  report: {out_dir / 'report.md'}")
    print(f"  summary: {out_dir / 'summary.csv'}")
    if not auc_summary.empty:
        print(f"  AUC: {auc_summary['roc_auc'].iloc[0]:.3f}")
    write_analysis_mode(out_dir, "evaluation", "KO labels were provided; true-vs-virtual metrics and AUC/R2-style summaries are valid for these held-out KO targets.")


def run_import_data(args: argparse.Namespace) -> None:
    from .importers import import_single_cell_data, write_import_outputs

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_h5ad = Path(args.output_h5ad) if args.output_h5ad else out_dir / "imported_data.h5ad"
    adata = import_single_cell_data(
        input_path=args.input,
        input_format=args.format,
        out_h5ad=out_h5ad,
        ko_metadata_csv=args.metadata_csv,
        cell_id_col=args.cell_id_col,
    )
    write_import_outputs(adata, out_dir, args.format)
    print("Imported single-cell data:")
    print(f"  h5ad: {out_h5ad}")
    print(f"  overview figure: {out_dir}\\input_overview.png")
    print(f"  report: {out_dir}\\import_report.md")
    print(f"  cells: {adata.n_obs}")
    print(f"  RNA features: {adata.n_vars}")


def run_validate_benchmark(args: argparse.Namespace) -> None:
    from .benchmark import validate_multimodal_benchmark

    result = validate_multimodal_benchmark(
        input_h5ad=args.input_h5ad,
        ko_col=args.ko_col,
        extra_obsm=args.extra_obsm,
        out_dir=args.out_dir,
    )
    summary = result["summary"]
    verdict = summary.loc[summary["check"] == "benchmark_mode", "status"].iloc[0]
    print("Validated multimodal perturbation benchmark input:")
    print(f"  verdict: {verdict}")
    print(f"  overview figure: {args.out_dir}\\benchmark_overview.png")
    print(f"  modality figure: {args.out_dir}\\benchmark_modalities.png")
    print(f"  report: {args.out_dir}\\benchmark_readiness_report.md")


def run_benchmark_registry(args: argparse.Namespace) -> None:
    from .benchmark import public_benchmark_registry

    table = public_benchmark_registry(args.out_dir)
    print("Saved public benchmark registry:")
    print(f"  registry: {args.out_dir}\\public_benchmark_registry.csv")
    print(f"  report: {args.out_dir}\\public_benchmark_registry_report.md")
    print(table[["dataset", "status", "recommended_use"]].to_string(index=False))


def run_workflow_template(args: argparse.Namespace) -> None:
    from .workflow import write_workflow_template

    path = write_workflow_template(args.mode, args.out_dir)
    print("Saved workflow template:")
    print(f"  template: {path}")


def run_method_comparison(args: argparse.Namespace) -> None:
    from .comparison import make_method_comparison

    metrics = [item.strip() for item in (args.metric_csv or "").split(",") if item.strip()]
    result = make_method_comparison(metrics, args.out_dir)
    print("Saved method comparison outputs:")
    print(f"  registry: {args.out_dir}\\method_registry.csv")
    print(f"  report: {args.out_dir}\\method_comparison_report.md")
    print(f"  method positioning figure: {args.out_dir}\\01_method_positioning.png")
    if not result["metrics"].empty:
        print(f"  metric comparison: {args.out_dir}\\method_metric_comparison.csv")
        print(f"  metric figure: {args.out_dir}\\02_metric_comparison.png")


def run_ko_cards(args: argparse.Namespace) -> None:
    from .cards import make_ko_summary_cards

    index = make_ko_summary_cards(
        delta_csv=args.delta_csv,
        out_dir=args.out_dir,
        auc_csv=args.auc_csv,
        confidence_csv=args.confidence_csv,
        max_features=args.max_features,
    )
    print("Saved KO summary cards:")
    print(f"  index: {args.out_dir}\\ko_summary_cards_index.csv")
    print(f"  report: {args.out_dir}\\ko_summary_cards_report.md")
    print(f"  cards: {len(index)}")


def run_figure_package(args: argparse.Namespace) -> None:
    from .figure_pack import make_figure_package

    index = make_figure_package(args.result_dir, args.out_dir)
    out = Path(args.out_dir) if args.out_dir else Path(args.result_dir) / "figure_package"
    print("Saved figure package:")
    print(f"  report: {out}\\figure_package_report.md")
    print(f"  index: {out}\\figure_index.csv")
    print(f"  figures: {len(index)}")


def run_diagnose_results(args: argparse.Namespace) -> None:
    from .diagnostics import diagnose_virtual_ko_results

    result = diagnose_virtual_ko_results(
        delta_csv=args.delta_csv,
        manifest_csv=args.manifest_csv,
        confidence_csv=args.confidence_csv,
        out_dir=args.out_dir,
        min_true_delta=args.min_true_delta,
        large_error=args.large_error,
        max_features_per_ko=args.max_features_per_ko,
    )
    print("Saved virtual KO diagnosis outputs:")
    print(f"  KO diagnosis: {args.out_dir}\\ko_failure_diagnosis.csv")
    print(f"  feature diagnosis: {args.out_dir}\\feature_failure_diagnosis.csv")
    print(f"  report: {args.out_dir}\\failure_diagnosis_report.md")
    print(f"  KO rows: {len(result['ko'])}")
    print(f"  feature rows: {len(result['feature'])}")


def run_summarize_result(args: argparse.Namespace) -> None:
    from .result_report import summarize_result_directory

    result = summarize_result_directory(args.result_dir, args.out_dir)
    out = Path(args.out_dir) if args.out_dir else Path(args.result_dir) / "readable_result_report"
    print("Saved user-readable result report:")
    print(f"  report: {result['report']}")
    print(f"  figure package: {out}\\figure_package\\figure_package_report.md")
    print(f"  KO cards: {len(result['ko_cards'])}")
    print(f"  KO diagnosis rows: {len(result['diagnosis_ko'])}")
    print(f"  feature diagnosis rows: {len(result['diagnosis_feature'])}")


def run_formal_benchmark_command(args: argparse.Namespace) -> None:
    from .formal_benchmark import run_formal_benchmark
    from .figure_pack import make_figure_package

    methods = [part.strip() for part in (args.methods or "").split(",") if part.strip()] or None
    result = run_formal_benchmark(
        state_csv=args.state_csv,
        ko_col=args.ko_col,
        target_kos=target_labels(args),
        prior_dir=args.prior_dir,
        out_dir=args.out_dir,
        features=parse_features(args.features),
        methods=methods,
        external_predictions_csv=args.external_predictions_csv,
        calibrate=args.calibrate,
        shape_calibrate=args.shape_calibrate,
        seed=args.seed,
    )
    make_figure_package(args.out_dir, Path(args.out_dir) / "figure_package")
    print("Saved formal benchmark outputs:")
    print(f"  metrics: {args.out_dir}\\formal_benchmark_metrics.csv")
    print(f"  method comparison: {args.out_dir}\\method_metric_comparison.csv")
    print(f"  report: {args.out_dir}\\formal_benchmark_report.md")
    print(f"  figure package: {args.out_dir}\\figure_package\\figure_package_report.md")
    print(f"  scored rows: {len(result['metrics'])}")


def run_train_hard_generator_command(args: argparse.Namespace) -> None:
    from .hard_generator import train_hard_constrained_generator
    from .figure_pack import make_figure_package

    state_csvs = [part.strip() for part in args.state_csvs.split(",") if part.strip()]
    result = train_hard_constrained_generator(
        state_csvs=state_csvs,
        ko_col=args.ko_col,
        target_kos=target_labels(args),
        prior_dir=args.prior_dir,
        out_dir=args.out_dir,
        features=parse_features(args.features),
        samples_per_ko=args.samples_per_ko,
        max_residual_fraction=args.max_residual_fraction,
        epochs=args.epochs,
        seed=args.seed,
    )
    make_figure_package(args.out_dir, Path(args.out_dir) / "figure_package")
    print("Saved hard-constrained generator outputs:")
    print(f"  virtual cells: {args.out_dir}\\hard_generator_virtual_cells.csv")
    print(f"  intervals: {args.out_dir}\\hard_generator_intervals.csv")
    print(f"  metrics: {args.out_dir}\\hard_generator_metrics.csv")
    print(f"  report: {args.out_dir}\\hard_generator_report.md")
    print(f"  backend: {result['model_status']}")


def run_assemble_multiome(args: argparse.Namespace) -> None:
    from .multiome import assemble_multiome_benchmark

    adata = assemble_multiome_benchmark(
        rna_input=args.rna_input,
        rna_format=args.rna_format,
        atac_input=args.atac_input,
        atac_format=args.atac_format,
        metadata_csv=args.metadata_csv,
        output_h5ad=args.output_h5ad,
        out_dir=args.out_dir,
        cell_id_col=args.cell_id_col,
        ko_col=args.ko_col,
        max_atac_features=args.max_atac_features,
    )
    print("Assembled RNA+ATAC perturbation benchmark:")
    print(f"  h5ad: {args.output_h5ad}")
    print(f"  overview figure: {args.out_dir}\\multiome_assembly_overview.png")
    print(f"  report: {args.out_dir}\\multiome_assembly_report.md")
    print(f"  shared cells: {adata.n_obs}")
    print(f"  RNA features: {adata.n_vars}")
    print(f"  ATAC features: {adata.obsm['atac'].shape[1] if 'atac' in adata.obsm else 0}")


def run_annotate_peaks(args: argparse.Namespace) -> None:
    from .peak_annotation import annotate_peaks

    table = annotate_peaks(
        input_h5ad=args.input_h5ad,
        obsm_key=args.obsm_key,
        out_csv=args.out_csv,
        feature_names_csv=args.feature_names_csv,
        gene_tss_csv=args.gene_tss_csv,
        motif_hits_csv=args.motif_hits_csv,
        marker_peaks_csv=args.marker_peaks_csv,
        target_genes=args.target_genes,
        max_distance=args.max_distance,
    )
    print("Annotated ATAC peak features:")
    print(f"  annotation CSV: {args.out_csv}")
    print(f"  report: {Path(args.out_csv).with_suffix('.report.md')}")
    print(f"  features: {len(table)}")
    print(f"  nonzero regulatory prior: {(table['regulatory_prior_score'] > 0).sum()}")


def run_make_gene_tss(args: argparse.Namespace) -> None:
    from .peak_annotation import make_gene_tss_from_gtf

    table = make_gene_tss_from_gtf(
        gtf=args.gtf,
        out_csv=args.out_csv,
        feature_type=args.feature_type,
        gene_name_attr=args.gene_name_attr,
        gene_id_attr=args.gene_id_attr,
    )
    print("Created gene TSS table:")
    print(f"  gene TSS CSV: {args.out_csv}")
    print(f"  report: {Path(args.out_csv).with_suffix('.report.md')}")
    print(f"  rows: {len(table)}")


def run_standardize_peak_scores(args: argparse.Namespace) -> None:
    from .peak_annotation import standardize_peak_score_table

    table = standardize_peak_score_table(
        input_csv=args.input_csv,
        out_csv=args.out_csv,
        table_type=args.table_type,
        peak_col=args.peak_col,
        score_col=args.score_col,
        tf_col=args.tf_col,
    )
    print("Standardized peak score table:")
    print(f"  output CSV: {args.out_csv}")
    print(f"  rows: {len(table)}")


def run_build_peak_annotation(args: argparse.Namespace) -> None:
    from .peak_annotation import build_peak_annotation_pipeline

    table = build_peak_annotation_pipeline(
        out_csv=args.out_csv,
        input_h5ad=args.input_h5ad,
        obsm_key=args.obsm_key,
        feature_names_csv=args.feature_names_csv,
        gtf=args.gtf,
        gene_tss_csv=args.gene_tss_csv,
        raw_motif_hits_csv=args.raw_motif_hits_csv,
        motif_hits_csv=args.motif_hits_csv,
        raw_marker_peaks_csv=args.raw_marker_peaks_csv,
        marker_peaks_csv=args.marker_peaks_csv,
        target_genes=args.target_genes,
        max_distance=args.max_distance,
    )
    print("Built peak annotation workflow outputs:")
    print(f"  annotation CSV: {args.out_csv}")
    print(f"  report: {Path(args.out_csv).with_suffix('.report.md')}")
    print(f"  summary figure: {Path(args.out_csv).with_suffix('.summary.png')}")
    print(f"  features: {len(table)}")


def run_from_raw(args: argparse.Namespace) -> None:
    from .core import run_virtual_ko
    from .preprocess import csv_matrix_to_state_table, h5ad_to_state_table, parse_extra_obsm_specs
    from .visualization import make_all_plots

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.input_h5ad:
        frame, manifest = h5ad_to_state_table(
            input_h5ad=args.input_h5ad,
            ko_col=args.ko_col,
            prior_dir=args.prior_dir,
            max_pathways=args.max_pathways,
            protein_obsm=args.protein_obsm,
            protein_prefix=args.protein_prefix,
            extra_obsm=parse_extra_obsm_specs(args.extra_obsm),
            max_extra_features_per_obsm=args.max_extra_features_per_obsm,
            extra_feature_selection=args.extra_feature_selection,
            extra_feature_metadata_csv=args.extra_feature_metadata_csv,
        )
    elif args.input_csv:
        frame, manifest = csv_matrix_to_state_table(
            input_csv=args.input_csv,
            ko_col=args.ko_col,
            prior_dir=args.prior_dir,
            max_pathways=args.max_pathways,
        )
    else:
        raise ValueError("Provide either --input-h5ad or --input-csv.")

    state_csv = out_dir / "derived_state_scores.csv"
    frame.to_csv(state_csv, index=False)
    manifest.to_csv(out_dir / "derived_state_manifest.csv", index=False)
    result = run_virtual_ko(
        frame=frame,
        ko_col="ko_target",
        holdouts=target_labels(args),
        prior_dir=args.prior_dir,
        features=None,
        dataset_name=args.dataset_name,
        modality=args.modality,
        representation=args.representation,
        calibrate=args.calibrate,
        shape_calibrate=args.shape_calibrate,
        max_cells_per_state=args.max_cells_per_state,
        seed=args.seed,
    )
    result.metrics.to_csv(out_dir / "metrics.csv", index=False)
    result.summary.to_csv(out_dir / "summary.csv", index=False)
    result.virtual_cells.to_csv(out_dir / "virtual_cells.csv", index=False)
    result.delta_table.to_csv(out_dir / "delta_table.csv", index=False)
    result.auc_points.to_csv(out_dir / "auc_points.csv", index=False)
    result.calibration.to_csv(out_dir / "calibration.csv", index=False)
    auc_summary = make_all_plots(result.summary, result.delta_table, result.virtual_cells, result.auc_points, out_dir)
    _auto_make_cards(out_dir, "delta_table.csv", auc_name="auc_summary.csv")
    _auto_diagnose(out_dir, "delta_table.csv", manifest_name="derived_state_manifest.csv")
    _auto_summarize(out_dir)
    print("Saved raw-matrix virtual KO outputs:")
    print(f"  derived state scores: {state_csv}")
    print(f"  report: {out_dir / 'report.md'}")
    print(f"  summary: {out_dir / 'summary.csv'}")
    if not auc_summary.empty:
        print(f"  AUC: {auc_summary['roc_auc'].iloc[0]:.3f}")
    write_analysis_mode(out_dir, "evaluation", "Raw matrix included KO labels; software derived state scores internally and evaluated held-out KO targets.")


def run_score(args: argparse.Namespace) -> None:
    from .preprocess import h5ad_to_state_scores, parse_extra_obsm_specs

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frame, manifest = h5ad_to_state_scores(
        input_h5ad=args.input_h5ad,
        prior_dir=args.prior_dir,
        max_pathways=args.max_pathways,
        protein_obsm=args.protein_obsm,
        protein_prefix=args.protein_prefix,
        extra_obsm=parse_extra_obsm_specs(args.extra_obsm),
        max_extra_features_per_obsm=args.max_extra_features_per_obsm,
        extra_feature_selection=args.extra_feature_selection,
        extra_feature_metadata_csv=args.extra_feature_metadata_csv,
    )
    frame.to_csv(out_dir / "derived_state_scores.csv", index=False)
    manifest.to_csv(out_dir / "derived_state_manifest.csv", index=False)
    print("Saved derived state scores for ordinary single-cell data:")
    print(f"  state scores: {out_dir / 'derived_state_scores.csv'}")
    print(f"  manifest: {out_dir / 'derived_state_manifest.csv'}")
    print("Note: without KO/perturbation labels, this step prepares the cells for virtual KO application but cannot evaluate against real KO.")
    write_analysis_mode(out_dir, "state_scoring_only", "No KO labels are used here. The output prepares ordinary cells for reference-model virtual KO, but no accuracy metrics can be computed.")


def run_train_reference(args: argparse.Namespace) -> None:
    from .preprocess import parse_extra_obsm_specs
    from .reference import train_reference_model

    reference = train_reference_model(
        input_h5ad=args.input_h5ad,
        state_csv=args.state_csv,
        ko_col=args.ko_col,
        prior_dir=args.prior_dir,
        output_model=args.output_model,
        max_pathways=args.max_pathways,
        protein_obsm=args.protein_obsm,
        protein_prefix=args.protein_prefix,
        extra_obsm=parse_extra_obsm_specs(args.extra_obsm),
        max_extra_features_per_obsm=args.max_extra_features_per_obsm,
        extra_feature_selection=args.extra_feature_selection,
        dataset_name=args.dataset_name,
        batch_col=args.batch_col,
        interaction_mode=args.interaction_mode,
        extra_feature_metadata_csv=args.extra_feature_metadata_csv,
    )
    print("Saved reference virtual KO model:")
    print(f"  model: {args.output_model}")
    print(f"  metadata: {args.output_model}.metadata.json")
    print(f"  training KO labels: {len(reference['training_ko_labels'])}")
    print(f"  state features: {len(reference['features'])}")
    print(f"  interaction residual: {reference.get('interaction_status')}")
    print(f"  batch covariate: {reference.get('batch_col') or 'not provided'}")


def run_apply_reference(args: argparse.Namespace) -> None:
    from .reference import apply_reference_model

    cells, deltas = apply_reference_model(
        reference_model=args.reference_model,
        input_h5ad=args.input_h5ad,
        state_csv=args.state_csv,
        target_kos=target_labels(args),
        out_dir=args.out_dir,
        max_cells=args.max_cells,
        seed=args.seed,
        cell_type_col=args.cell_type_col,
        batch_col=args.batch_col,
        uncertainty_method=args.uncertainty_method,
        uncertainty_scale=args.uncertainty_scale,
        uncertainty_samples_per_ko=args.uncertainty_samples_per_ko,
        extra_feature_metadata_csv=args.extra_feature_metadata_csv,
    )
    print("Applied reference virtual KO model:")
    print(f"  virtual cells: {args.out_dir}\\applied_virtual_cells.csv")
    print(f"  predicted delta: {args.out_dir}\\predicted_ko_delta.csv")
    print(f"  target interpretation: {args.out_dir}\\target_interpretation.csv")
    print(f"  prediction-only report: {args.out_dir}\\prediction_only_report.md")
    print(f"  n cells rows: {len(cells)}")
    print(f"  n target KO: {len(deltas)}")
    write_analysis_mode(Path(args.out_dir), "prediction_only", "A saved perturbation reference model was applied to input cells. Without matching real KO labels in this input dataset, AUC/R2/MAE accuracy metrics are not computed.")
    _auto_make_cards(Path(args.out_dir), "predicted_ko_delta.csv", confidence_name="transfer_confidence.csv")
    _auto_diagnose(Path(args.out_dir), "predicted_ko_delta.csv", confidence_name="transfer_confidence.csv")
    _auto_summarize(Path(args.out_dir))


def run_inspect_reference(args: argparse.Namespace) -> None:
    from .reference import inspect_reference_model

    result = inspect_reference_model(
        reference_model=args.reference_model,
        out_dir=args.out_dir,
        target_kos=target_labels(args) if args.target_kos or args.holdouts else None,
    )
    summary = result["summary"]
    print("Inspected reference virtual KO model:")
    print(f"  dataset: {summary.get('dataset_name')}")
    print(f"  version: {summary.get('version')}")
    print(f"  training KO labels: {summary.get('n_training_ko_labels')}")
    print(f"  training genes: {summary.get('n_training_genes')}")
    print(f"  state features: {summary.get('n_state_features')}")
    print(f"  prior terms: {summary.get('n_prior_terms')}")
    if args.out_dir:
        print(f"  report: {args.out_dir}\\reference_inspection_report.md")


def run_double_interaction_command(args: argparse.Namespace) -> None:
    from .interaction import run_double_interaction

    metrics, predictions = run_double_interaction(
        delta_csv=args.delta_csv,
        ko_col=args.ko_col,
        n_ko_col=args.n_ko_col,
        target_prefix=args.target_prefix,
        prior_dir=args.prior_dir,
        out_dir=args.out_dir,
    )
    summary = metrics.groupby(["model", "subset"], observed=True)[["mae", "r2", "roc_auc_abs_gt_0.15"]].mean().round(3)
    print("Saved double-KO interaction outputs:")
    print(f"  metrics: {args.out_dir}\\double_interaction_metrics.csv")
    print(f"  predictions: {args.out_dir}\\double_interaction_predictions.csv")
    print(f"  figure: {args.out_dir}\\double_interaction_metrics.png")
    print(summary.to_string())
    write_analysis_mode(Path(args.out_dir), "double_ko_evaluation", "Double-KO labels are present; this command compares additive single-gene prediction with a prior-based interaction residual model.")


def write_analysis_mode(out_dir: Path, mode: str, interpretation: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    text = f"""# Analysis mode

- Mode: `{mode}`
- Interpretation: {interpretation}

Use this file to decide which plots and metrics are appropriate:

- `evaluation`: real KO labels are available, so heatmaps, UMAP, AUC and error metrics can be interpreted as accuracy evidence.
- `prediction_only`: virtual KO is applied to ordinary or unlabeled cells, so outputs show predicted state shifts only.
- `state_scoring_only`: the input matrix is converted into pathway/program/protein state scores for later use.
- `double_ko_evaluation`: true double-KO labels are available, so additive and interaction models can be compared.
"""
    (out_dir / "analysis_mode.md").write_text(text, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m vkx.cli",
        description="Prior-constrained pathway/program-level virtual knockout interface.",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    importer = sub.add_parser("import-data", help="Import h5ad, 10x mtx, 10x h5, or basic h5Seurat into workflow-ready h5ad with an overview figure.")
    importer.add_argument("--input", required=True, help="Input file or folder.")
    importer.add_argument("--format", required=True, choices=["h5ad", "10x_mtx", "10x_h5", "h5seurat"], help="Input format.")
    importer.add_argument("--out-dir", default="results/import_data_demo")
    importer.add_argument("--output-h5ad", default=None, help="Optional output h5ad path. Default: out-dir/imported_data.h5ad.")
    importer.add_argument("--metadata-csv", default=None, help="Optional metadata CSV to join by cell id, e.g. KO labels or cell types.")
    importer.add_argument("--cell-id-col", default="cell_id", help="Cell id column in metadata CSV.")
    importer.set_defaults(func=run_import_data)
    bench = sub.add_parser("validate-benchmark", help="Check whether an h5ad is ready for labelled multimodal perturbation benchmarking and make overview plots.")
    bench.add_argument("--input-h5ad", required=True)
    bench.add_argument("--ko-col", default="ko_target")
    bench.add_argument("--extra-obsm", default=None, help="Comma-separated modality specs to check, e.g. protein:protein,atac:atac,chromvar:tf,peak:peak.")
    bench.add_argument("--out-dir", default="results/benchmark_readiness")
    bench.set_defaults(func=run_validate_benchmark)
    registry = sub.add_parser("benchmark-registry", help="Write the public multimodal perturbation benchmark registry and interpretation boundaries.")
    registry.add_argument("--out-dir", default="results/public_benchmark_registry")
    registry.set_defaults(func=run_benchmark_registry)
    workflow = sub.add_parser("workflow-template", help="Write a copy-paste workflow template for common labelled, reference, prediction-only, and ATAC peak use cases.")
    workflow.add_argument("--mode", choices=["all", "labelled-benchmark", "reference", "prediction-only", "atac-peak"], default="all")
    workflow.add_argument("--out-dir", default="results/workflow_template")
    workflow.set_defaults(func=run_workflow_template)
    method_cmp = sub.add_parser("method-comparison", help="Create conceptual and empirical comparison figures against existing virtual perturbation methods.")
    method_cmp.add_argument("--metric-csv", default=None, help="Optional comma-separated metric CSV files with method/model and AUC/R2/MAE/direction columns.")
    method_cmp.add_argument("--out-dir", default="results/method_comparison")
    method_cmp.set_defaults(func=run_method_comparison)
    cards = sub.add_parser("ko-cards", help="Create one readable summary card per KO from delta_table.csv or predicted_ko_delta.csv.")
    cards.add_argument("--delta-csv", required=True, help="delta_table.csv from evaluation or predicted_ko_delta.csv from apply-reference.")
    cards.add_argument("--auc-csv", default=None, help="Optional auc_summary.csv.")
    cards.add_argument("--confidence-csv", default=None, help="Optional transfer_confidence.csv.")
    cards.add_argument("--max-features", type=int, default=8)
    cards.add_argument("--out-dir", default="results/ko_summary_cards")
    cards.set_defaults(func=run_ko_cards)
    figures = sub.add_parser("figure-package", help="Collect generated PNG figures into one readable figure package with explanations.")
    figures.add_argument("--result-dir", required=True, help="Existing result directory to scan.")
    figures.add_argument("--out-dir", default=None, help="Output directory. Default: result-dir/figure_package.")
    figures.set_defaults(func=run_figure_package)
    diagnose = sub.add_parser("diagnose-results", help="Explain which virtual KO results are reliable or risky and why.")
    diagnose.add_argument("--delta-csv", required=True, help="delta_table.csv from evaluation or predicted_ko_delta.csv from apply-reference.")
    diagnose.add_argument("--manifest-csv", default=None, help="Optional derived_state_manifest.csv for modality-aware explanations.")
    diagnose.add_argument("--confidence-csv", default=None, help="Optional transfer_confidence.csv from apply-reference.")
    diagnose.add_argument("--min-true-delta", type=float, default=0.03, help="True-delta magnitude below which direction/AUC interpretation is unstable.")
    diagnose.add_argument("--large-error", type=float, default=None, help="Optional absolute-error threshold. Default: data-adaptive 75th percentile.")
    diagnose.add_argument("--max-features-per-ko", type=int, default=25)
    diagnose.add_argument("--out-dir", default="results/failure_diagnosis")
    diagnose.set_defaults(func=run_diagnose_results)
    summarize = sub.add_parser("summarize-result", help="Build one user-readable report from an existing VKX result directory.")
    summarize.add_argument("--result-dir", required=True, help="Existing VKX result directory containing delta_table.csv or predicted_ko_delta.csv.")
    summarize.add_argument("--out-dir", default=None, help="Output directory. Default: result-dir/readable_result_report.")
    summarize.set_defaults(func=run_summarize_result)
    formal = sub.add_parser("formal-benchmark", help="Run a formal method benchmark: VKX vs PLS/Ridge/Additive plus external scGen/CPA/GEARS/CellOT predictions when provided.")
    formal.add_argument("--state-csv", required=True, help="State score CSV with one row per cell, KO labels, and numeric state features.")
    formal.add_argument("--ko-col", default="ko_target")
    formal.add_argument("--target-kos", required=True, help="Comma-separated held-out KO targets, e.g. STAT1,JAK2,STAT1+JAK2.")
    formal.add_argument("--prior-dir", default="data/priors")
    formal.add_argument("--methods", default="vkx,pls,ridge,additive,scgen,cpa,gears,cellot")
    formal.add_argument("--external-predictions-csv", default=None, help="Optional predictions from scGen/CPA/GEARS/CellOT with method, ko_target, and pred_delta_* columns.")
    formal.add_argument("--features", default=None)
    formal.add_argument("--calibrate", choices=["auto", "none", "global_scale", "feature_scale"], default="auto")
    formal.add_argument("--shape-calibrate", choices=["none", "variance", "quantile"], default="none")
    formal.add_argument("--seed", type=int, default=7)
    formal.add_argument("--out-dir", default="results/formal_method_benchmark")
    formal.set_defaults(func=run_formal_benchmark_command)
    hard_gen = sub.add_parser("train-hard-generator", help="Train a lightweight hard-constrained residual generator from one or more perturbation state CSVs.")
    hard_gen.add_argument("--state-csvs", required=True, help="Comma-separated state CSV files. They must share state feature names.")
    hard_gen.add_argument("--ko-col", default="ko_target")
    hard_gen.add_argument("--target-kos", required=True, help="Comma-separated KO targets to generate.")
    hard_gen.add_argument("--prior-dir", default="data/priors")
    hard_gen.add_argument("--features", default=None)
    hard_gen.add_argument("--samples-per-ko", type=int, default=300)
    hard_gen.add_argument("--max-residual-fraction", type=float, default=0.35, help="Hard bound on residual norm relative to baseline KO delta norm.")
    hard_gen.add_argument("--epochs", type=int, default=80)
    hard_gen.add_argument("--seed", type=int, default=11)
    hard_gen.add_argument("--out-dir", default="results/hard_constrained_generator")
    hard_gen.set_defaults(func=run_train_hard_generator_command)
    multiome = sub.add_parser("assemble-multiome", help="Assemble RNA and ATAC matrices plus KO metadata into one benchmark-ready h5ad.")
    multiome.add_argument("--rna-input", required=True)
    multiome.add_argument("--rna-format", required=True, choices=["h5ad", "10x_mtx", "10x_h5"])
    multiome.add_argument("--atac-input", required=True)
    multiome.add_argument("--atac-format", required=True, choices=["h5ad", "10x_mtx", "10x_h5"])
    multiome.add_argument("--metadata-csv", required=True)
    multiome.add_argument("--cell-id-col", default="cell_id")
    multiome.add_argument("--ko-col", default="ko_target")
    multiome.add_argument("--max-atac-features", type=int, default=500)
    multiome.add_argument("--output-h5ad", required=True)
    multiome.add_argument("--out-dir", default="results/multiome_assembly")
    multiome.set_defaults(func=run_assemble_multiome)
    annotate = sub.add_parser("annotate-peaks", help="Build a peak annotation CSV for ATAC peak-gene linkage and motif-to-peak prior weighting.")
    annotate.add_argument("--input-h5ad", default=None, help="h5ad containing peak-level features in obsm.")
    annotate.add_argument("--obsm-key", default="peak", help="obsm key for peak-level ATAC features.")
    annotate.add_argument("--feature-names-csv", default=None, help="Optional CSV with feature_name/peak/name if peak names are not stored in h5ad.uns.")
    annotate.add_argument("--gene-tss-csv", default=None, help="Optional gene annotation CSV with gene, chrom, and tss or start/end columns.")
    annotate.add_argument("--motif-hits-csv", default=None, help="Optional motif-to-peak CSV with peak/feature_name and score/weight columns.")
    annotate.add_argument("--marker-peaks-csv", default=None, help="Optional marker peak CSV with peak/feature_name and marker_score/score columns.")
    annotate.add_argument("--target-genes", default=None, help="Optional comma-separated target genes used to upweight nearby genes or matching TF motifs.")
    annotate.add_argument("--max-distance", type=int, default=250000, help="Maximum distance scale for peak-gene linkage scoring.")
    annotate.add_argument("--out-csv", required=True, help="Output peak annotation CSV for --extra-feature-metadata-csv.")
    annotate.set_defaults(func=run_annotate_peaks)
    tss = sub.add_parser("make-gene-tss", help="Convert a GTF/GFF-style gene annotation file into the gene_tss.csv used by annotate-peaks.")
    tss.add_argument("--gtf", required=True, help="GENCODE/Ensembl-style GTF file.")
    tss.add_argument("--out-csv", required=True, help="Output gene_tss.csv.")
    tss.add_argument("--feature-type", default="gene", help="Feature type to read from the GTF, usually gene.")
    tss.add_argument("--gene-name-attr", default="gene_name", help="GTF attribute used as gene symbol.")
    tss.add_argument("--gene-id-attr", default="gene_id", help="GTF attribute used as stable gene id.")
    tss.set_defaults(func=run_make_gene_tss)
    peak_scores = sub.add_parser("standardize-peak-scores", help="Convert motif/marker peak score tables into the columns expected by annotate-peaks.")
    peak_scores.add_argument("--input-csv", required=True)
    peak_scores.add_argument("--out-csv", required=True)
    peak_scores.add_argument("--table-type", choices=["motif", "marker"], default="motif")
    peak_scores.add_argument("--peak-col", default=None, help="Column containing peak names. Auto-detected when omitted.")
    peak_scores.add_argument("--score-col", default=None, help="Column containing motif/marker score. Auto-detected when omitted.")
    peak_scores.add_argument("--tf-col", default=None, help="Optional TF/motif/gene column.")
    peak_scores.set_defaults(func=run_standardize_peak_scores)
    peak_pipeline = sub.add_parser("build-peak-annotation", help="Run the full helper workflow to create peak_annotation.csv from GTF, motif hits, marker peaks, and peak names.")
    peak_pipeline.add_argument("--input-h5ad", default=None, help="h5ad containing peak-level features in obsm.")
    peak_pipeline.add_argument("--obsm-key", default="peak")
    peak_pipeline.add_argument("--feature-names-csv", default=None, help="Optional CSV with feature_name/peak/name.")
    peak_pipeline.add_argument("--gtf", default=None, help="Optional GTF used to build gene_tss.csv.")
    peak_pipeline.add_argument("--gene-tss-csv", default=None, help="Optional precomputed gene_tss.csv.")
    peak_pipeline.add_argument("--raw-motif-hits-csv", default=None, help="Optional raw motif hit table to standardize.")
    peak_pipeline.add_argument("--motif-hits-csv", default=None, help="Optional already standardized motif-to-peak CSV.")
    peak_pipeline.add_argument("--raw-marker-peaks-csv", default=None, help="Optional raw marker peak table to standardize.")
    peak_pipeline.add_argument("--marker-peaks-csv", default=None, help="Optional already standardized marker peak CSV.")
    peak_pipeline.add_argument("--target-genes", default=None)
    peak_pipeline.add_argument("--max-distance", type=int, default=250000)
    peak_pipeline.add_argument("--out-csv", required=True)
    peak_pipeline.set_defaults(func=run_build_peak_annotation)
    fit = sub.add_parser("fit", help="Fit/evaluate virtual KO from a state-score CSV table.")
    fit.add_argument("--state-csv", required=True, help="CSV with one row per cell, one KO label column, and numeric state-score columns.")
    fit.add_argument("--ko-col", default="ko_target", help="Column containing control/KO labels.")
    fit.add_argument("--target-kos", default=None, help="KO target(s) to predict/evaluate. Usually one gene (STAT1) or one pair (STAT1+JAK2).")
    fit.add_argument("--holdouts", default=None, help="Alias for --target-kos, kept for evaluation scripts.")
    fit.add_argument("--prior-dir", default="data/priors", help="Directory containing Reactome/MSigDB/TF/PPI GMT files.")
    fit.add_argument("--out-dir", default="results/software_interface_demo", help="Output directory.")
    fit.add_argument("--features", default=None, help="Optional comma-separated state-score columns. Default: all numeric columns.")
    fit.add_argument("--dataset-name", default="Virtual KO dataset")
    fit.add_argument("--modality", default="state score table")
    fit.add_argument("--representation", default="pathway/program scores")
    fit.add_argument("--calibrate", choices=["auto", "none", "global_scale", "feature_scale"], default="auto")
    fit.add_argument("--shape-calibrate", choices=["none", "variance", "quantile"], default="none", help="Optional distribution-shape calibration. Use variance or quantile for sparse ATAC/peak features.")
    fit.add_argument("--max-cells-per-state", type=int, default=180)
    fit.add_argument("--seed", type=int, default=7)
    fit.set_defaults(func=run_fit)
    raw = sub.add_parser("run", help="Run virtual KO directly from a raw single-cell matrix.")
    raw.add_argument("--input-h5ad", default=None, help="Raw AnnData h5ad with cells x genes in X and KO labels in obs.")
    raw.add_argument("--input-csv", default=None, help="Raw CSV with one row per cell, gene columns, and a KO label column.")
    raw.add_argument("--ko-col", default="ko_target", help="KO label column in obs or CSV.")
    raw.add_argument("--target-kos", default=None, help="KO target(s) to predict/evaluate. Usually one gene (STAT1) or one pair (STAT1+JAK2).")
    raw.add_argument("--holdouts", default=None, help="Alias for --target-kos, kept for evaluation scripts.")
    raw.add_argument("--prior-dir", default="data/priors", help="Directory containing pathway/network GMT files.")
    raw.add_argument("--out-dir", default="results/software_interface_raw_demo")
    raw.add_argument("--dataset-name", default="Virtual KO dataset")
    raw.add_argument("--modality", default="raw RNA matrix")
    raw.add_argument("--representation", default="auto-derived pathway/program scores")
    raw.add_argument("--max-pathways", type=int, default=40)
    raw.add_argument("--protein-obsm", default=None, help="Optional obsm key for protein/ADT or other cell-by-feature modality.")
    raw.add_argument("--protein-prefix", default="protein")
    raw.add_argument("--extra-obsm", default=None, help="Comma-separated extra obsm specs, e.g. protein:protein,atac:atac. Supersedes --protein-obsm for new workflows.")
    raw.add_argument("--max-extra-features-per-obsm", type=int, default=None, help="Optional cap for each extra obsm modality, useful for large chromVAR/motif/peak matrices.")
    raw.add_argument("--extra-feature-selection", choices=["variance", "ko_effect", "hybrid", "atac_peak"], default="variance", help="How to choose capped extra obsm features. Use atac_peak for sparse raw peak/count matrices.")
    raw.add_argument("--extra-feature-metadata-csv", default=None, help="Optional feature annotation CSV for extra obsm features, e.g. peak-gene linkage or motif-to-peak scores.")
    raw.add_argument("--calibrate", choices=["auto", "none", "global_scale", "feature_scale"], default="auto")
    raw.add_argument("--shape-calibrate", choices=["none", "variance", "quantile"], default="none", help="Optional distribution-shape calibration. Use variance or quantile for sparse ATAC/peak features.")
    raw.add_argument("--max-cells-per-state", type=int, default=180)
    raw.add_argument("--seed", type=int, default=7)
    raw.set_defaults(func=run_from_raw)
    score = sub.add_parser("score", help="Convert ordinary h5ad single-cell data into auto-derived pathway/protein state scores.")
    score.add_argument("--input-h5ad", required=True, help="Ordinary 10X/scRNA h5ad with cells x genes in X.")
    score.add_argument("--prior-dir", default="data/priors")
    score.add_argument("--out-dir", default="results/software_interface_score_only")
    score.add_argument("--max-pathways", type=int, default=40)
    score.add_argument("--protein-obsm", default=None, help="Optional obsm key for protein/ADT/ATAC/gene activity.")
    score.add_argument("--protein-prefix", default="protein")
    score.add_argument("--extra-obsm", default=None, help="Comma-separated extra obsm specs, e.g. protein:protein,atac:atac.")
    score.add_argument("--max-extra-features-per-obsm", type=int, default=None, help="Optional cap for each extra obsm modality.")
    score.add_argument("--extra-feature-selection", choices=["variance", "ko_effect", "hybrid", "atac_peak"], default="variance", help="Feature selection for extra obsm. Without KO labels, KO-effect methods fall back toward unsupervised scores.")
    score.add_argument("--extra-feature-metadata-csv", default=None, help="Optional feature annotation CSV for extra obsm features.")
    score.set_defaults(func=run_score)
    train_ref = sub.add_parser("train-reference", help="Train and save a reference KO-delta model from perturbation data.")
    train_ref.add_argument("--input-h5ad", default=None, help="Perturbation h5ad with RNA matrix and KO labels.")
    train_ref.add_argument("--state-csv", default=None, help="Optional precomputed state score CSV with KO labels.")
    train_ref.add_argument("--ko-col", default="ko_target")
    train_ref.add_argument("--prior-dir", default="data/priors")
    train_ref.add_argument("--output-model", required=True)
    train_ref.add_argument("--dataset-name", default="reference perturbation dataset")
    train_ref.add_argument("--max-pathways", type=int, default=40)
    train_ref.add_argument("--protein-obsm", default=None)
    train_ref.add_argument("--protein-prefix", default="protein")
    train_ref.add_argument("--extra-obsm", default=None, help="Comma-separated extra obsm specs, e.g. protein:protein,atac:atac.")
    train_ref.add_argument("--max-extra-features-per-obsm", type=int, default=None, help="Optional cap for each extra obsm modality during reference training.")
    train_ref.add_argument("--extra-feature-selection", choices=["variance", "ko_effect", "hybrid", "atac_peak"], default="variance", help="How to choose capped extra obsm features during reference training.")
    train_ref.add_argument("--extra-feature-metadata-csv", default=None, help="Optional feature annotation CSV for peak-gene linkage, motif-to-peak, marker, or regulatory prior weights.")
    train_ref.add_argument("--batch-col", default=None, help="Optional obs/state CSV batch/sample/donor column. Control cells are used to remove batch-specific baseline offsets.")
    train_ref.add_argument("--interaction-mode", choices=["auto", "on", "off"], default="auto", help="Train a double-KO interaction residual model when single- and double-KO labels are available.")
    train_ref.set_defaults(func=run_train_reference)

    apply_ref = sub.add_parser("apply-reference", help="Apply a saved reference KO model to ordinary unlabeled cells.")
    apply_ref.add_argument("--reference-model", required=True)
    apply_ref.add_argument("--input-h5ad", default=None, help="Ordinary h5ad to receive virtual KO.")
    apply_ref.add_argument("--state-csv", default=None, help="Optional precomputed state score CSV to receive virtual KO.")
    apply_ref.add_argument("--target-kos", default=None, help="Comma-separated targets. Each target can be one gene or one pair, e.g. STAT1,JAK2,STAT1+JAK2.")
    apply_ref.add_argument("--holdouts", default=None, help="Alias for --target-kos.")
    apply_ref.add_argument("--out-dir", default="results/reference_apply_demo")
    apply_ref.add_argument("--max-cells", type=int, default=800)
    apply_ref.add_argument("--cell-type-col", default=None, help="Optional obs/state CSV column for cell-type stratified prediction-only outputs.")
    apply_ref.add_argument("--batch-col", default=None, help="Optional obs/state CSV batch/sample/donor column for batch-composition output.")
    apply_ref.add_argument("--extra-feature-metadata-csv", default=None, help="Optional feature annotation CSV for extra obsm features; defaults to the path stored in the reference model when available.")
    apply_ref.add_argument("--uncertainty-method", choices=["none", "hard-residual", "vae", "flow", "diffusion"], default="none", help="Optional hard-constrained uncertainty band. VAE/flow/diffusion currently share the same hard residual anchor unless a custom generator is plugged in.")
    apply_ref.add_argument("--uncertainty-scale", type=float, default=0.25, help="Width multiplier for the hard-constrained residual uncertainty band.")
    apply_ref.add_argument("--uncertainty-samples-per-ko", type=int, default=250, help="Number of hard-constrained uncertainty virtual cells to sample per requested KO. Use 0 for interval-only output.")
    apply_ref.add_argument("--seed", type=int, default=7)
    apply_ref.set_defaults(func=run_apply_reference)
    inspect_ref = sub.add_parser("inspect-reference", help="Inspect a saved reference KO model before applying it.")
    inspect_ref.add_argument("--reference-model", required=True)
    inspect_ref.add_argument("--target-kos", default=None, help="Optional comma-separated KO targets to check prior coverage.")
    inspect_ref.add_argument("--holdouts", default=None, help="Alias for --target-kos.")
    inspect_ref.add_argument("--out-dir", default="results/reference_inspection")
    inspect_ref.set_defaults(func=run_inspect_reference)
    double_interaction = sub.add_parser("double-interaction", help="Evaluate double-KO effects with additive and prior-based interaction residual models.")
    double_interaction.add_argument("--delta-csv", required=True, help="KO-level delta table containing single and double KO rows.")
    double_interaction.add_argument("--ko-col", default="ko_genes", help="Column containing KO gene labels such as CEBPB+CEBPA.")
    double_interaction.add_argument("--n-ko-col", default="n_ko_genes", help="Column containing the number of KO genes.")
    double_interaction.add_argument("--target-prefix", default="delta_program_", help="Prefix of numeric effect columns to predict.")
    double_interaction.add_argument("--prior-dir", default="data/priors")
    double_interaction.add_argument("--out-dir", default="results/double_interaction_demo")
    double_interaction.set_defaults(func=run_double_interaction_command)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
