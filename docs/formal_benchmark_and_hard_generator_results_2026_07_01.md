# Formal Benchmark 与 Hard-constrained Generator 初版结果

更新时间：2026-07-01

本次更新新增两个正式接口：

- `formal-benchmark`
- `train-hard-generator`

并在已有 Papalexi/HMPCITE state-score 结果上做了 smoke-test 运行。

## 1. Papalexi 正式横向 benchmark

命令：

```bash
python -m vkx.cli formal-benchmark \
  --state-csv results/software_interface_raw_papalexi/derived_state_scores.csv \
  --ko-col ko_target \
  --target-kos STAT1,JAK2,IFNGR2,IRF1 \
  --prior-dir data/priors \
  --methods vkx,pls,ridge,additive,scgen,cpa,gears,cellot \
  --out-dir results/formal_benchmark_papalexi
```

平均结果：

| method | ROC-AUC | direction cosine | R2 | MAE | feature hit-rate |
|---|---:|---:|---:|---:|---:|
| Ridge | 0.822 | 0.751 | 0.251 | 0.309 | 0.679 |
| PLS | 0.807 | 0.751 | 0.220 | 0.312 | 0.643 |
| VKX | 0.704 | 0.462 | -0.080 | 0.386 | 0.357 |

解释：

- 在这个 Papalexi 单基因 holdout benchmark 上，Ridge/PLS 当前优于 VKX 默认校准版本。
- 这说明正式 benchmark 必须保留 classical baseline，不能只展示 VKX 自己的结果。
- Additive baseline 在单基因 holdout 下没有合法组件，因此被标记为 `run_partial`，不参与平均指标。
- scGen、CPA、GEARS、CellOT 当前被标记为 `not_run`，因为没有外部预测 CSV。后续需要真实运行这些方法后再导入。

输出图：

- `results/formal_benchmark_papalexi/01_formal_benchmark_metric_panel.png`
- `results/formal_benchmark_papalexi/02_formal_benchmark_delta_heatmap.png`
- `results/formal_benchmark_papalexi/03_method_availability.png`
- `results/formal_benchmark_papalexi/figure_package/figure_package_report.md`

## 2. Hard-constrained generator demo

命令：

```bash
python -m vkx.cli train-hard-generator \
  --state-csvs results/software_interface_raw_papalexi/derived_state_scores.csv,results/hmpcite_multimodal_doubleko_cebp_med12/derived_state_scores.csv \
  --ko-col ko_target \
  --target-kos STAT1,JAK2 \
  --prior-dir data/priors \
  --samples-per-ko 120 \
  --epochs 10 \
  --out-dir results/hard_generator_papalexi_hmpcite_demo
```

结果：

| KO | direction cosine | R2 | MAE |
|---|---:|---:|---:|
| STAT1 | 0.813 | 0.106 | 0.556 |
| JAK2 | 0.870 | -0.046 | 0.523 |

解释：

- 这个 generator 保持 VKX baseline KO delta 为 hard constraint，只学习 bounded residual。
- 当前本机没有可用 PyTorch 后端，因此自动退回 `pca_residual_fallback_no_torch`。
- 方向 cosine 仍较高，说明 hard constraint 有效；但 MAE 偏大，说明它不是幅度修正器，也不能替代 baseline。
- 后续应在可用 PyTorch 环境和更多同类型 perturbation 数据上训练 residual VAE，并与 PCA fallback 比较。

输出图：

- `results/hard_generator_papalexi_hmpcite_demo/01_hard_generator_metric_panel.png`
- `results/hard_generator_papalexi_hmpcite_demo/02_hard_generator_intervals.png`
- `results/hard_generator_papalexi_hmpcite_demo/figure_package/figure_package_report.md`

## 3. 当前结论

这次更新的意义不是宣布 VKX 全面优于已有方法，而是建立了一个更诚实、更可复用的比较框架：

1. VKX、Ridge、PLS、Additive 可以在同一输入和指标下比较。
2. scGen、CPA、GEARS、CellOT 有明确外部预测导入口，不再停留在概念比较。
3. 所有 benchmark 图会进入 figure package，方便组织论文级主图。
4. hard-constrained generator 已经可以从多数据 residual bank 采样，但仍保持 KO 方向约束。

下一步重点：

- 为 scGen/CPA/GEARS/CellOT 生成真实 external prediction CSV。
- 在 Norman/HMPCITE/ATAC 数据上批量运行 `formal-benchmark`。
- 在可用 PyTorch 环境下训练 residual VAE，并比较 PCA fallback。
