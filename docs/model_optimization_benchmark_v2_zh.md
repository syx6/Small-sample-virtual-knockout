# 模型深度优化与 7 张大图 benchmark v2

本次优化目标是：让 VKX 在小样本 benchmark 中不要弱于强 classical baseline，并把 benchmark 结果整理成 7 张读者可直接理解的大图。

## 1. 模型优化：VKX-Adaptive

新增 `VKX-Adaptive`。它不是固定使用某一个 anchor，而是在训练 KO 上做交叉验证，比较：

- PLS anchor；
- Ridge anchor；
- PLS/Ridge mean anchor。

然后选择交叉验证综合分数最高的 anchor 去预测 held-out KO。

在当前 Papalexi 小型 benchmark 上，VKX-Adaptive 自动选择了 PLS：

```text
selected=PLS
MeanPLSRidge=0.441
PLS=0.446
Ridge=0.441
```

## 2. 当前 benchmark 结果

输出目录：

```text
results/formal_benchmark_papalexi_adaptive_v2
results/paper_benchmark_package_adaptive_v2
```

方法比较：

| method | AUC | direction | R2 | MAE | hit-rate |
|---|---:|---:|---:|---:|---:|
| PLS | 0.736 | 0.690 | 0.234 | 0.140 | 0.625 |
| VKX-Adaptive | 0.736 | 0.690 | 0.234 | 0.140 | 0.625 |
| CalibratedEnsemble | 0.667 | 0.659 | 0.065 | 0.162 | 0.625 |
| VKX-Boosted | 0.625 | 0.608 | 0.184 | 0.139 | 0.625 |
| Ensemble | 0.583 | 0.611 | 0.172 | 0.145 | 0.500 |
| Ridge | 0.583 | 0.611 | 0.172 | 0.145 | 0.500 |
| VKX | 0.389 | 0.494 | -0.013 | 0.171 | 0.000 |

结论：

- VKX-Adaptive 已经达到 PLS 水平；
- VKX-Boosted 的 AUC 不是第一，但 MAE 最低；
- 原始 VKX 明显不足，必须使用 adaptive/boosted/calibrated anchor；
- scGen/CPA/GEARS/CellOT 仍然需要同数据 prediction CSV 才能正式评分。

## 3. 7 张大图

1. `00_publication_main_figure.png`  
   论文级 A-D 主图。

2. `01_method_leaderboard.png`  
   正式方法指标比较。

3. `02_auc_roc_curves.png`  
   AUC 的 ROC 曲线，不再用柱状图冒充曲线。

4. `03_real_vs_virtual_method_heatmap.png`  
   真实 KO 与虚拟 KO 的 heatmap。

5. `04_single_double_multimodal_gallery.png`  
   单敲、双敲、多模态和 ATAC peak 总览。

6. `05_adaptive_improvement.png`  
   VKX-Adaptive 是否追上强 baseline。

7. `06_benchmark_completeness.png`  
   哪些 benchmark 已完成，哪些仍是缺口。

## 4. 仍未补齐的真实缺口

- scGen、CPA、GEARS、CellOT 需要真实运行并提供 external prediction CSV；
- 真正公开 RNA+ADT+ATAC 且带 perturbation 标签的数据集仍需确认；
- 多数据集正式 benchmark 需要继续扩展到 Norman、HMPCITE、ATAC 和更多 Perturb-seq 数据；
- 双敲非线性和 MAPK/TGFB program 仍需要专项增强。
