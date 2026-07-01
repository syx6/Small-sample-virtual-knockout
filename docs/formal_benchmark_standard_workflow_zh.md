# 正式横向 benchmark 标准流程

本项目现在把正式 benchmark 和论文主图拆成两个连续步骤：

## 1. 运行正式横向 benchmark

```bash
python -m vkx.cli formal-benchmark \
  --state-csv data/examples/papalexi_pathway_protein_state.csv \
  --ko-col ko_target \
  --target-kos STAT1,JAK2,IFNGR2,IRF1 \
  --prior-dir data/priors \
  --out-dir results/formal_benchmark_papalexi_standard
```

默认比较方法已经改为：

```text
adaptive, boosted, ensemble, calibrated, vkx, pls, ridge, additive, scgen, cpa, gears, cellot
```

其中：

- `adaptive` 是 VKX-Adaptive：在训练 KO 上做交叉验证，自动选择 PLS/Ridge/mean anchor；
- `boosted` 是当前 VKX-Boosted；
- `ensemble` 是 constrained PLS/Ridge anchor；
- `calibrated` 是带幅度校准的 ensemble；
- `pls` 和 `ridge` 是必须保留的经典 baseline；
- `scgen`、`cpa`、`gears`、`cellot` 是外部方法槽位，只有提供同数据预测文件时才会评分。

## 2. 生成论文级主图

```bash
python -m vkx.cli paper-benchmark \
  --formal-dir results/formal_benchmark_papalexi_standard \
  --result-dirs results/software_interface_single_gene_demo,results/software_interface_double_gene_demo,results/hmpcite_multimodal_doubleko_cebp_med12,results/scperturb_atac_regulatory_peak_prior_quantile_kdm6a \
  --out-dir results/paper_benchmark_package_standard
```

最重要的输出是：

- `00_publication_main_figure.png`：论文级 A-D 主图；
- `01_method_leaderboard.png`：方法指标比较；
- `02_auc_roc_curves.png`：ROC 曲线；
- `03_real_vs_virtual_method_heatmap.png`：真实 KO vs 虚拟 KO heatmap；
- `04_single_double_multimodal_gallery.png`：单敲、双敲、多模态和 ATAC 总览。

## 当前标准 benchmark 结果

在 `data/examples/papalexi_pathway_protein_state.csv` 这个小型 Papalexi benchmark 上：

| method | AUC | direction | R2 | MAE |
|---|---:|---:|---:|---:|
| PLS | 0.736 | 0.690 | 0.234 | 0.140 |
| VKX-Adaptive | 0.736 | 0.690 | 0.234 | 0.140 |
| CalibratedEnsemble | 0.667 | 0.659 | 0.065 | 0.162 |
| VKX-Boosted | 0.625 | 0.608 | 0.184 | 0.139 |
| Ensemble | 0.583 | 0.611 | 0.172 | 0.145 |
| Ridge | 0.583 | 0.611 | 0.172 | 0.145 |
| VKX | 0.389 | 0.494 | -0.013 | 0.171 |

解释：

- PLS 在这个小型 benchmark 上 AUC 最高，说明传统 baseline 很强；
- VKX-Adaptive 自动选择 PLS anchor，因此在这个 benchmark 上与 PLS 持平；
- VKX-Boosted 不是 AUC 第一，但 MAE 最低，说明它在变化幅度误差上有优势；
- 原始 VKX 明显较弱，说明 boosted/prior/calibrated anchor 是必要优化；
- scGen/CPA/GEARS/CellOT 目前没有同数据预测文件，因此不能宣称 VKX 强过它们。

## 外部方法导入要求

要让 scGen/CPA/GEARS/CellOT 进入正式评分，需要提供：

```text
method,ko_target,pred_delta_feature1,pred_delta_feature2,...
```

其中 `feature1/feature2` 必须和 state CSV 中的 feature 名一致。导入后，`formal-benchmark` 会自动把这些方法纳入同一套 AUC/R2/MAE/direction 指标和论文主图。
