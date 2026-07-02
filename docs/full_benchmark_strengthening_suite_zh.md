# VKX 补强总包：正式横向 benchmark 与论文级图包

本文档说明当前版本新增的 `benchmark-suite` 流程。它的目标是把分散的模型优化、横向比较、外部方法接入模板和可视化结果整理成一个可复用的软件接口。

## 1. 为什么需要这个总包

前面的实验已经说明：VKX 在小样本 perturbation 数据里可以稳定预测 KO 方向，但如果要说服别人，必须回答三个问题：

1. 和 ridge、PLS、additive、scGen、CPA、GEARS、CellOT 等方法相比怎么样？
2. 单敲、双敲、多模态、ATAC peak 这些场景的结果能不能放到同一套图里看？
3. 哪些方法真的跑了，哪些只是因为缺少外部预测或缺少 labelled trimodal 数据而暂时不能评估？

`benchmark-suite` 就是为这三个问题设计的统一入口。

## 2. 一条命令运行

默认会运行仓库内已有的快速示例数据：

```bash
python -m vkx.cli benchmark-suite --out-dir results/full_benchmark_suite
```

输出包括：

- `benchmark_suite_report_zh.md`: 中文总报告。
- `benchmark_suite_jobs.csv`: 本次运行了哪些数据、哪些 KO。
- `aggregate/formal_method_metrics_aggregate.csv`: 所有正式 benchmark 的分数。
- `aggregate/formal_best_methods.csv`: 每个 benchmark 的最佳方法。
- `paper_figures/`: 论文级主图包。
- `top_figures/`: 复制出的核心 11 张主图，方便直接查看。
- `formal_*/external_prediction_template.csv`: scGen、CPA、GEARS、CellOT 的预测填表模板。

## 3. 默认运行哪些数据

当前默认包含两个 labelled perturbation 示例：

| dataset_id | 数据类型 | 任务 |
|---|---|---|
| `papalexi_rna_protein_single_ko` | RNA pathway + protein/ADT state | 单基因 KO benchmark |
| `norman_rna_program_double_ko` | RNA program score | 双基因 KO benchmark，属于较慢示例，需要加 `--include-long-examples` |

运行较慢的内置双敲示例：

```bash
python -m vkx.cli benchmark-suite \
  --include-long-examples \
  --out-dir results/full_benchmark_suite_with_double_ko
```

如果用户有自己的 labelled perturbation 数据，可以提供配置表：

```csv
dataset_id,state_csv,ko_col,target_kos,features,prior_dir
your_dataset,path/to/state_score_table.csv,ko_target,"GENE1,GENE2,GENE1+GENE2",,data/priors
```

然后运行：

```bash
python -m vkx.cli benchmark-suite \
  --suite-csv your_suite.csv \
  --out-dir results/your_full_benchmark
```

## 4. 输入文件到底是什么

`benchmark-suite` 使用的是“已经转换好的 state score table”，每一行是一个细胞：

```csv
cell_id,ko_target,pathway_IFNG_JAK_STAT,protein_CD86,protein_PDL1,chromvar_IRF1,peak_chrX_123_456
cell_001,control,0.12,0.44,0.31,-0.05,0
cell_002,STAT1,-0.35,0.20,0.18,-0.42,1
```

这里的通路分数不是用户必须手工准备的生信结果。普通用户可以先用 `from-raw`、`score` 或 reference workflow 从原始 10X/Seurat/AnnData 矩阵自动转换：

- RNA matrix -> pathway/program score。
- ADT/protein -> protein state feature。
- ATAC gene activity -> regulatory gene activity score。
- chromVAR -> motif/TF activity score。
- raw peak/count -> peak-level accessibility feature。

也就是说，普通用户的原始输入仍然可以是单细胞矩阵；benchmark 阶段用的是软件内部统一后的可解释状态表示。

## 5. 横向比较方法

当前 suite 会直接评分：

- `VKXAdaptive`: 自适应选择 PLS/Ridge/mean anchor，保证小样本下不被固定模型拖累。
- `ResponseBoosted`: 对强 KO 响应方向进行增强的 hard-constrained baseline。
- `ConstrainedEnsemble`: PLS/Ridge 误差加权 ensemble。
- `CalibratedEnsemble`: 进一步做幅度校准。
- `VKX`: 原始 prior-constrained residual baseline。
- `PLS`: prior state 到 KO delta 的 PLS baseline。
- `Ridge`: ridge regression baseline。
- `Additive`: 双敲时的单基因效应相加 baseline。

外部深度方法：

- `scGen`
- `CPA`
- `GEARS`
- `CellOT`

这些方法不会被假装运行。suite 会生成 `external_prediction_template.csv`，用户或后续脚本把外部方法对同一 held-out KO 的 `pred_delta_*` 填进去后，再重新运行 suite，才会纳入同口径评分。

## 6. 核心 11 张主图

`paper_figures/` 和 `top_figures/` 中会自动生成：

1. `00_publication_main_figure.png`: 总览主图，展示方法、benchmark 和主要结论。
2. `01_method_leaderboard.png`: 方法排行榜，直接看谁更好。
3. `02_auc_roc_curves.png`: AUC 的曲线形式，不再只给柱状图。
4. `03_real_vs_virtual_method_heatmap.png`: 真实 KO delta 与虚拟 KO delta 的 heatmap。
5. `04_single_double_multimodal_gallery.png`: 单敲、双敲、多模态结果汇总。
6. `05_adaptive_improvement.png`: VKXAdaptive 相比原始 VKX/PLS/Ridge 的改进。
7. `06_benchmark_completeness.png`: 哪些 benchmark 已完成，哪些因为缺少外部预测或数据而不能评分。
8. `07_before_after_umap_panel.png`: control、virtual KO、true KO 的单细胞状态移动图。
9. `08_single_double_response_map.png`: single KO、double KO、多模态和 ATAC 场景的响应强弱图。
10. `09_peak_locus_track.png`: ATAC peak locus track，展示真实和虚拟 peak accessibility delta。
11. `10_method_radar_leaderboard.png`: 多指标 radar/leaderboard，展示方法权衡和外部方法状态。

## 7. 仍然不能跳过的现实边界

当前版本支持多模态输入，也能在有多模态 labelled perturbation 数据时把多模态信息用于更可靠的 KO 预测。但是：

- 没有 KO 标签的普通 10X、DOGMA/TEA-seq 只能做 prediction-only 或 reference application，不能在该数据内部报告真实准确率。
- 公开 RNA+ADT+ATAC 且带 perturbation 标签的数据集仍需要确认和下载；在确认前不能声称完成 full trimodal labelled benchmark。
- scGen、CPA、GEARS、CellOT 必须按各自官方流程训练/推理，不能用占位分数比较。
- ATAC peak 的单细胞分布形状仍然困难，quantile/zero-inflated calibration 可以改善分布，但不会保证每个 peak 的真实开放比例都完美恢复。

## 8. 当前方法意义

VKX 的定位不是替代所有大规模深度 perturbation 模型，而是在小样本、多模态、算力有限的场景里提供一个可解释、可校准、可评估的稳定 baseline。它的核心优势是：

- 输入门槛低：可以从普通单细胞矩阵转换到 pathway/program state。
- 多模态友好：RNA、ADT、ATAC、chromVAR、peak 都能进入统一状态空间。
- 小样本稳定：使用 PLS/residual/prior hard constraint，而不是从零训练自由生成模型。
- 结果可解释：输出 KO 方向、真实/预测 delta、通路变化、AUC 曲线、heatmap、UMAP/figure package。
- 边界清楚：没有标签的数据不假装 benchmark；没有外部预测的方法不假装比较。

这个版本使 VKX 更像一个可以交给别人用的软件流程，而不是只在本地脚本里跑通的实验。
