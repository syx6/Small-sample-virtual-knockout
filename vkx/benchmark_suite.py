from __future__ import annotations

from dataclasses import dataclass
import html
from pathlib import Path
import shutil

import pandas as pd

from .aggregate import aggregate_benchmarks
from .formal_benchmark import EXTERNAL_METHODS, _feature_columns, run_formal_benchmark
from .paper_benchmark import make_paper_benchmark_package


DEFAULT_METHODS = [
    "adaptive",
    "boosted",
    "ensemble",
    "calibrated",
    "vkx",
    "pls",
    "ridge",
    "additive",
    "scgen",
    "cpa",
    "gears",
    "cellot",
]

DEFAULT_RESULT_DIRS = [
    "results/software_interface_single_gene_demo",
    "results/software_interface_double_gene_demo",
    "results/hmpcite_multimodal_doubleko_cebp_med12",
    "results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a",
]


@dataclass
class BenchmarkJob:
    dataset_id: str
    state_csv: Path
    ko_col: str
    target_kos: list[str]
    features: list[str] | None = None
    prior_dir: Path | None = None


def run_benchmark_suite(
    out_dir: str | Path,
    suite_csv: str | Path | None = None,
    prior_dir: str | Path = "data/priors",
    external_predictions_csv: str | Path | None = None,
    include_default_examples: bool = True,
    include_long_examples: bool = False,
    methods: list[str] | None = None,
    seed: int = 7,
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    default_prior = Path(prior_dir)
    selected_methods = methods or DEFAULT_METHODS

    jobs = _load_jobs(suite_csv, default_prior) if suite_csv else []
    if include_default_examples:
        jobs.extend(_default_jobs(default_prior, include_long_examples=include_long_examples))
    jobs = _deduplicate_jobs(jobs)
    if not jobs:
        _write_config_template(out)
        _write_empty_report(out)
        return {
            "jobs": pd.DataFrame(),
            "formal": pd.DataFrame(),
            "best": pd.DataFrame(),
            "paper_metrics": pd.DataFrame(),
            "paper_summary": pd.DataFrame(),
        }

    job_table = _job_table(jobs)
    job_table.to_csv(out / "benchmark_suite_jobs.csv", index=False)
    _write_config_template(out)

    formal_dirs: list[Path] = []
    template_paths: list[Path] = []
    for job in jobs:
        formal_dir = out / f"formal_{_safe_id(job.dataset_id)}"
        run_formal_benchmark(
            state_csv=job.state_csv,
            ko_col=job.ko_col,
            target_kos=job.target_kos,
            prior_dir=job.prior_dir or default_prior,
            out_dir=formal_dir,
            features=job.features,
            methods=selected_methods,
            external_predictions_csv=external_predictions_csv,
            calibrate="auto",
            shape_calibrate="none",
            seed=seed,
        )
        formal_dirs.append(formal_dir)
        template_paths.append(_write_external_prediction_template(job, formal_dir))

    aggregate = aggregate_benchmarks(formal_dirs, out / "aggregate")
    paper_result_dirs = _existing_result_dirs(DEFAULT_RESULT_DIRS)
    paper = make_paper_benchmark_package(
        formal_dir=formal_dirs[0],
        result_dirs=[str(path) for path in paper_result_dirs],
        out_dir=out / "paper_figures",
    )
    figure_index = _write_figure_index(out, formal_dirs, paper_result_dirs)
    _write_suite_report(out, jobs, formal_dirs, template_paths, aggregate, paper, figure_index, external_predictions_csv)
    _write_suite_html_report(out, jobs, formal_dirs, template_paths, aggregate, paper, figure_index, external_predictions_csv)
    return {
        "jobs": job_table,
        "formal": aggregate["formal"],
        "best": aggregate["best"],
        "paper_metrics": paper["metrics"],
        "paper_summary": paper["result_summary"],
    }


def _load_jobs(suite_csv: str | Path, default_prior: Path) -> list[BenchmarkJob]:
    table = pd.read_csv(suite_csv)
    required = {"dataset_id", "state_csv", "target_kos"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(f"suite CSV missing columns: {', '.join(sorted(missing))}")
    jobs = []
    for _, row in table.iterrows():
        targets = _parse_list(row.get("target_kos", ""))
        if not targets:
            continue
        jobs.append(
            BenchmarkJob(
                dataset_id=str(row["dataset_id"]),
                state_csv=Path(str(row["state_csv"])),
                ko_col=str(row.get("ko_col", "ko_target") or "ko_target"),
                target_kos=targets,
                features=_parse_list(row.get("features", "")) or None,
                prior_dir=Path(str(row.get("prior_dir", "") or default_prior)),
            )
        )
    return jobs


def _default_jobs(default_prior: Path, include_long_examples: bool = False) -> list[BenchmarkJob]:
    jobs = []
    papalexi = Path("data/examples/papalexi_pathway_protein_state.csv")
    if papalexi.exists():
        jobs.append(
            BenchmarkJob(
                dataset_id="papalexi_rna_protein_single_ko",
                state_csv=papalexi,
                ko_col="ko_target",
                target_kos=["STAT1", "JAK2", "IFNGR2", "IRF1"],
                prior_dir=default_prior,
            )
        )
    norman = Path("data/examples/norman_gene_program_state.csv")
    if include_long_examples and norman.exists():
        jobs.append(
            BenchmarkJob(
                dataset_id="norman_rna_program_double_ko",
                state_csv=norman,
                ko_col="ko_target",
                target_kos=["AHR+KLF1", "CEBPB+MAPK1", "CEBPE+CEBPA", "DUSP9+SNAI1"],
                prior_dir=default_prior,
            )
        )
    return jobs


def _deduplicate_jobs(jobs: list[BenchmarkJob]) -> list[BenchmarkJob]:
    seen = set()
    unique = []
    for job in jobs:
        key = (_safe_id(job.dataset_id), str(job.state_csv), ",".join(job.target_kos))
        if key in seen:
            continue
        if not job.state_csv.exists():
            continue
        unique.append(job)
        seen.add(key)
    return unique


def _job_table(jobs: list[BenchmarkJob]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "dataset_id": job.dataset_id,
                "state_csv": str(job.state_csv),
                "ko_col": job.ko_col,
                "target_kos": ",".join(job.target_kos),
                "features": ",".join(job.features or []),
                "prior_dir": str(job.prior_dir or ""),
            }
            for job in jobs
        ]
    )


def _write_external_prediction_template(job: BenchmarkJob, out: Path) -> Path:
    frame = pd.read_csv(job.state_csv, nrows=200)
    features = _feature_columns(frame, job.ko_col, job.features)
    rows = []
    for method in EXTERNAL_METHODS:
        for ko in job.target_kos:
            row = {
                "method": method,
                "ko_target": ko,
                "prediction_status": "fill_external_method_prediction_here",
            }
            for feature in features:
                row[f"pred_delta_{feature}"] = ""
            rows.append(row)
    path = out / "external_prediction_template.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _write_config_template(out: Path) -> None:
    template = pd.DataFrame(
        [
            {
                "dataset_id": "your_labelled_perturbation_dataset",
                "state_csv": "path/to/state_score_table.csv",
                "ko_col": "ko_target",
                "target_kos": "GENE1,GENE2,GENE1+GENE2",
                "features": "",
                "prior_dir": "data/priors",
            }
        ]
    )
    template.to_csv(out / "benchmark_suite_config_template.csv", index=False)


def _existing_result_dirs(paths: list[str]) -> list[Path]:
    return [Path(path) for path in paths if Path(path).exists()]


def _write_figure_index(out: Path, formal_dirs: list[Path], result_dirs: list[Path]) -> pd.DataFrame:
    rows = []
    for fig in sorted((out / "paper_figures").glob("*.png")):
        rows.append({"section": "paper_figures", "figure": fig.name, "path": str(fig)})
    for formal_dir in formal_dirs:
        for fig in sorted(formal_dir.glob("*.png")):
            rows.append({"section": formal_dir.name, "figure": fig.name, "path": str(fig)})
    for result_dir in result_dirs:
        for fig in sorted(result_dir.glob("*.png")):
            rows.append({"section": result_dir.name, "figure": fig.name, "path": str(fig)})
    table = pd.DataFrame(rows)
    table.to_csv(out / "benchmark_suite_figure_index.csv", index=False)
    _copy_top_figures(out, table)
    return table


def _copy_top_figures(out: Path, figure_index: pd.DataFrame) -> None:
    gallery = out / "top_figures"
    gallery.mkdir(parents=True, exist_ok=True)
    preferred = [
        "00_publication_main_figure.png",
        "01_method_leaderboard.png",
        "02_auc_roc_curves.png",
        "03_real_vs_virtual_method_heatmap.png",
        "04_single_double_multimodal_gallery.png",
        "05_adaptive_improvement.png",
        "06_benchmark_completeness.png",
    ]
    for name in preferred:
        rows = figure_index.loc[figure_index["figure"] == name] if not figure_index.empty else pd.DataFrame()
        if rows.empty:
            continue
        src = Path(str(rows.iloc[0]["path"]))
        if src.exists():
            shutil.copy2(src, gallery / name)


def _write_empty_report(out: Path) -> None:
    text = """# VKX benchmark suite report

No runnable benchmark jobs were found.

Use `benchmark_suite_config_template.csv` to provide labelled perturbation state tables. A labelled table needs one row per cell, a KO label column, and numeric state features such as pathway/program scores, protein scores, ATAC gene activity, chromVAR, or peak features.
"""
    (out / "benchmark_suite_report_zh.md").write_text(text, encoding="utf-8")
    _write_empty_html_report(out)


def _write_suite_report(
    out: Path,
    jobs: list[BenchmarkJob],
    formal_dirs: list[Path],
    template_paths: list[Path],
    aggregate: dict[str, pd.DataFrame],
    paper: dict[str, pd.DataFrame],
    figure_index: pd.DataFrame,
    external_predictions_csv: str | Path | None,
) -> None:
    best = aggregate.get("best", pd.DataFrame())
    metrics = aggregate.get("formal", pd.DataFrame())
    lines = [
        "# VKX 补强总包 benchmark 报告",
        "",
        "## 这次补齐了什么",
        "",
        "- 一条命令运行多数据集正式 benchmark。",
        "- 同时比较 VKXAdaptive、ResponseBoosted、ConstrainedEnsemble、Calibrated、VKX、PLS、Ridge 和 Additive。",
        "- 为 scGen、CPA、GEARS、CellOT 生成严格横向比较模板；只有用户提供同一 held-out KO 的外部预测后才计入分数。",
        "- 自动生成聚合 leaderboard、AUC ROC 曲线、真实/虚拟 KO heatmap 和 7 张论文级主图。",
        "- 保留 batch covariate、double-KO interaction residual、ATAC peak annotation/quantile calibration 等已集成能力的入口边界。",
        "",
        "## 本次实际运行的数据",
        "",
    ]
    for job, formal_dir in zip(jobs, formal_dirs):
        lines.append(f"- `{job.dataset_id}`: `{job.state_csv}`; held-out KO = `{', '.join(job.target_kos)}`; 结果目录 `{formal_dir}`。")
    lines.extend(["", "## 当前结论", ""])
    if not best.empty:
        for _, row in best.iterrows():
            lines.append(
                f"- `{row['benchmark_id']}` 最优方法为 `{row['best_method']}`："
                f"AUC={_fmt(row.get('roc_auc'))}, R2={_fmt(row.get('r2'))}, MAE={_fmt(row.get('mae'))}。"
            )
    else:
        lines.append("- 暂无可聚合分数。")
    lines.extend(["", "## 外部方法 benchmark 怎么补齐", ""])
    if external_predictions_csv:
        lines.append(f"- 本次已读取外部预测表：`{external_predictions_csv}`。")
    else:
        lines.append("- 本次没有外部预测表，所以 scGen/CPA/GEARS/CellOT 标记为 `not_run`，不会被假装计分。")
    for path in template_paths:
        lines.append(f"- 外部方法预测模板：`{path}`。把对应方法预测的 `pred_delta_*` 填入后重新运行 suite 即可。")
    lines.extend(
        [
            "",
            "## 图在哪里",
            "",
            "- 7 张主图在 `paper_figures/` 和 `top_figures/`。",
            "- 完整图片索引在 `benchmark_suite_figure_index.csv`。",
            f"- 本次共索引 `{len(figure_index)}` 张 PNG 图。",
            "",
            "## 仍然不能伪造的部分",
            "",
            "- 真正公开 RNA+ADT+ATAC 且带 perturbation 标签的数据集如果没有下载并确认，不能当 full trimodal labelled benchmark。",
            "- 无 KO 标签的普通 10X、DOGMA/TEA-seq 只能做 reference application 或 prediction-only report，不能在该数据内部报告真实准确率。",
            "- 深度方法 scGen/CPA/GEARS/CellOT 需要按各自官方流程独立训练/推理，然后把预测表交给本 suite 做同口径评分。",
            "",
            "## 输出文件",
            "",
            "- `benchmark_suite_jobs.csv`: 本次运行任务表。",
            "- `aggregate/formal_method_metrics_aggregate.csv`: 所有正式 benchmark 分数。",
            "- `aggregate/formal_best_methods.csv`: 每个 benchmark 的最优方法。",
            "- `paper_figures/paper_benchmark_report_zh.md`: 论文图包解释。",
            "- `benchmark_suite_config_template.csv`: 用户自有 labelled perturbation 数据配置模板。",
        ]
    )
    if not metrics.empty:
        scored = metrics.loc[metrics["method"].astype(str).str.contains("VKX|PLS|Ridge|Additive|Ensemble|Boosted|Calibrated", regex=True, na=False)]
        lines.extend(["", f"本次 scored baseline 行数：`{len(scored)}`。"])
    if not paper.get("metrics", pd.DataFrame()).empty:
        lines.append(f"论文图包方法行数：`{len(paper['metrics'])}`。")
    (out / "benchmark_suite_report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_empty_html_report(out: Path) -> None:
    html_text = """<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <title>VKX benchmark suite report</title>
  <style>
    body { font-family: Arial, "Microsoft YaHei", sans-serif; margin: 40px; color: #222; line-height: 1.65; }
    main { max-width: 1080px; margin: 0 auto; }
    .note { border-left: 5px solid #2f9e8f; background: #f3fbf9; padding: 14px 18px; }
    code { background: #f5f5f5; padding: 2px 5px; border-radius: 4px; }
  </style>
</head>
<body><main>
  <h1>VKX benchmark suite report</h1>
  <p class="note">No runnable benchmark jobs were found.</p>
  <p>Use <code>benchmark_suite_config_template.csv</code> to provide labelled perturbation state tables.</p>
</main></body></html>
"""
    (out / "benchmark_suite_report_zh.html").write_text(html_text, encoding="utf-8")


def _write_suite_html_report(
    out: Path,
    jobs: list[BenchmarkJob],
    formal_dirs: list[Path],
    template_paths: list[Path],
    aggregate: dict[str, pd.DataFrame],
    paper: dict[str, pd.DataFrame],
    figure_index: pd.DataFrame,
    external_predictions_csv: str | Path | None,
) -> None:
    best = aggregate.get("best", pd.DataFrame())
    metrics = aggregate.get("formal", pd.DataFrame())
    paper_metrics = paper.get("metrics", pd.DataFrame())
    main_figures = [
        ("00_publication_main_figure.png", "Figure 1. 总览主图：方法、横向比较、真实/虚拟 KO 和多场景覆盖。"),
        ("01_method_leaderboard.png", "Figure 2. 方法排行榜：AUC/R2/MAE 放在同一张图里。"),
        ("02_auc_roc_curves.png", "Figure 3. AUC 曲线：展示强响应特征排序能力。"),
        ("03_real_vs_virtual_method_heatmap.png", "Figure 4. 真实 KO 与虚拟 KO 的 delta heatmap。"),
        ("04_single_double_multimodal_gallery.png", "Figure 5. 单敲、双敲、多模态和 ATAC 结果汇总。"),
        ("05_adaptive_improvement.png", "Figure 6. VKXAdaptive 相对原始 VKX 和 classical baseline 的变化。"),
        ("06_benchmark_completeness.png", "Figure 7. 哪些方法和数据已经评分，哪些仍是待补齐槽位。"),
    ]
    cards = []
    if not best.empty:
        for _, row in best.iterrows():
            cards.append(
                f"<li><b>{_h(row['benchmark_id'])}</b>: best = <b>{_h(row['best_method'])}</b>, "
                f"AUC={_fmt(row.get('roc_auc'))}, R2={_fmt(row.get('r2'))}, MAE={_fmt(row.get('mae'))}</li>"
            )
    else:
        cards.append("<li>暂无可聚合分数。</li>")
    jobs_html = "\n".join(
        f"<tr><td>{_h(job.dataset_id)}</td><td><code>{_h(job.state_csv)}</code></td>"
        f"<td>{_h(', '.join(job.target_kos))}</td><td><code>{_h(formal_dir)}</code></td></tr>"
        for job, formal_dir in zip(jobs, formal_dirs)
    )
    templates_html = "\n".join(f"<li><code>{_h(path)}</code></li>" for path in template_paths)
    figure_html = "\n".join(
        _figure_block(out / "top_figures" / filename, caption)
        for filename, caption in main_figures
        if (out / "top_figures" / filename).exists()
    )
    scored_rows = 0
    if not metrics.empty and "method" in metrics.columns:
        scored_rows = int(
            metrics["method"]
            .astype(str)
            .str.contains("VKX|PLS|Ridge|Additive|Ensemble|Boosted|Calibrated", regex=True, na=False)
            .sum()
        )
    external_status = "已提供外部预测表。" if external_predictions_csv else "未提供外部预测表，scGen/CPA/GEARS/CellOT 不会被假装计分。"
    html_text = f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <title>VKX 补强总包 benchmark 报告</title>
  <style>
    :root {{ --ink: #202124; --muted: #5f6368; --line: #dadce0; --brand: #2f9e8f; --soft: #f3fbf9; }}
    body {{ font-family: Arial, "Microsoft YaHei", "PingFang SC", sans-serif; margin: 0; color: var(--ink); background: #fff; line-height: 1.65; }}
    main {{ max-width: 1220px; margin: 0 auto; padding: 42px 34px 70px; }}
    h1 {{ font-size: 34px; margin: 0 0 8px; }}
    h2 {{ font-size: 23px; margin-top: 34px; border-bottom: 1px solid var(--line); padding-bottom: 8px; }}
    p, li, td, th {{ font-size: 15px; }}
    .subtitle {{ color: var(--muted); margin-top: 0; }}
    .summary {{ display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin: 24px 0; }}
    .card {{ border: 1px solid var(--line); border-radius: 8px; padding: 14px 16px; background: #fff; }}
    .card b {{ font-size: 24px; display: block; color: var(--brand); }}
    .note {{ border-left: 5px solid var(--brand); background: var(--soft); padding: 14px 18px; }}
    table {{ width: 100%; border-collapse: collapse; margin: 14px 0 20px; }}
    th, td {{ border-bottom: 1px solid var(--line); text-align: left; padding: 10px 8px; vertical-align: top; }}
    th {{ background: #f8f9fa; }}
    code {{ background: #f5f5f5; padding: 2px 5px; border-radius: 4px; }}
    figure {{ margin: 28px 0 36px; border: 1px solid var(--line); border-radius: 8px; padding: 14px; }}
    figure img {{ display: block; max-width: 100%; height: auto; margin: 0 auto; }}
    figcaption {{ color: var(--muted); font-size: 14px; margin-top: 10px; }}
    @media (max-width: 850px) {{ .summary {{ grid-template-columns: 1fr 1fr; }} main {{ padding: 26px 18px; }} }}
  </style>
</head>
<body><main>
  <h1>VKX 补强总包 benchmark 报告</h1>
  <p class="subtitle">一条命令生成正式横向 benchmark、外部方法模板、聚合表和论文级 7 张主图。</p>
  <section class="summary">
    <div class="card"><b>{len(jobs)}</b>运行数据集</div>
    <div class="card"><b>{scored_rows}</b>已评分 baseline 行</div>
    <div class="card"><b>{len(figure_index)}</b>索引 PNG 图</div>
    <div class="card"><b>{len(paper_metrics)}</b>论文图包方法行</div>
  </section>
  <section class="note">
    <b>当前边界：</b>{_h(external_status)} 没有 KO 标签的普通 10X、DOGMA/TEA-seq 只能做 prediction-only 或 reference application，不能在该数据内部报告真实准确率。
  </section>

  <h2>本次实际运行的数据</h2>
  <table><thead><tr><th>dataset_id</th><th>输入 state table</th><th>held-out KO</th><th>结果目录</th></tr></thead><tbody>{jobs_html}</tbody></table>

  <h2>当前结论</h2>
  <ul>{''.join(cards)}</ul>

  <h2>7 张主图</h2>
  {figure_html}

  <h2>外部方法怎么补齐</h2>
  <p>scGen、CPA、GEARS、CellOT 需要先按各自官方流程训练/推理，再把同一 held-out KO 的 <code>pred_delta_*</code> 填入模板。suite 只会评分真实提供的外部预测，不会把未运行方法当成结果。</p>
  <ul>{templates_html}</ul>

  <h2>输出文件</h2>
  <ul>
    <li><code>benchmark_suite_report_zh.md</code>: Markdown 报告。</li>
    <li><code>benchmark_suite_report_zh.html</code>: 本 HTML 图文报告。</li>
    <li><code>aggregate/formal_method_metrics_aggregate.csv</code>: 所有正式 benchmark 分数。</li>
    <li><code>aggregate/formal_best_methods.csv</code>: 每个 benchmark 的最优方法。</li>
    <li><code>benchmark_suite_figure_index.csv</code>: 全部图片索引。</li>
  </ul>
</main></body></html>
"""
    (out / "benchmark_suite_report_zh.html").write_text(html_text, encoding="utf-8")


def _figure_block(path: Path, caption: str) -> str:
    rel = path.parent.name + "/" + path.name
    return f'<figure><img src="{_h(rel)}" alt="{_h(caption)}"><figcaption>{_h(caption)}</figcaption></figure>'


def _h(value: object) -> str:
    return html.escape(str(value), quote=True)


def _parse_list(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).replace(";", ",")
    return [part.strip() for part in text.split(",") if part.strip()]


def _safe_id(value: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(value).strip())
    return safe or "benchmark"


def _fmt(value: object) -> str:
    try:
        return f"{float(value):.3f}"
    except Exception:
        return "NA"
