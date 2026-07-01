# Adaptive Response-Boosted Anchor

## 为什么需要这个模块

在 Papalexi RNA+ADT benchmark 中，Ridge/PLS/ensemble 已经能预测较好的 KO 方向，但 STAT1、JAK2、IFNGR2 这类 KO 的 interferon/JAK/STAT 相关通路变化幅度仍然偏小。

这说明当前主要问题不是“方向完全错”，而是“小样本下强响应 program 被保守估计”。如果直接让 neural generator 自由放大，会破坏 hard constraint 的可解释性；因此我们加入 adaptive response-strength prior，只在有生物学依据的 feature family 上增强响应幅度。

## 方法思想

基础预测仍然来自 constrained ensemble：

```text
base_delta = weighted(PLS_delta, Ridge_delta)
```

然后根据 KO gene 和 state feature family 计算 boost factor：

```text
final_delta[j] = base_delta[j] * boost_factor(KO genes, feature_j)
```

当前内置 family：

| family | 触发 feature | 触发 KO gene | boost |
|---|---|---|---|
| `ifn_jak_stat` | interferon / JAK / STAT / IFN | STAT*, JAK*, IFN*, IRF* | 1.50 |
| `mapk_tgfb` | MAPK / TGFB / ERK | MAPK*, RAF*, RAS*, MEK*, ERK*, TGFB*, SMAD* | 1.35 |
| `cell_cycle_myc_e2f` | MYC / E2F / G2M / cell cycle | MYC*, E2F*, CDK*, CCN*, RB* | 1.25 |
| `protein_checkpoint` | PDL1 / PDL2 / CD86 / CD366 / checkpoint protein | STAT*, JAK*, IFN*, IRF*, CD274*, PDCD1LG* | 1.15 |

没有匹配到 family 的 feature 不会被增强。

## 输入

输入和 formal benchmark / hard generator 相同：

```powershell
python -m vkx.cli formal-benchmark `
  --state-csv results\software_interface_raw_papalexi\derived_state_scores.csv `
  --ko-col ko_target `
  --target-kos STAT1,JAK2,IFNGR2,IRF1 `
  --prior-dir data\priors `
  --methods vkx,boosted,ensemble,pls,ridge,additive,scgen,cpa,gears,cellot `
  --out-dir results\formal_benchmark_papalexi_adaptive_boosted
```

Cell-level generator:

```powershell
python -m vkx.cli train-hard-generator `
  --state-csvs results\software_interface_raw_papalexi\derived_state_scores.csv,results\hmpcite_multimodal_doubleko_cebp_med12\derived_state_scores.csv `
  --ko-col ko_target `
  --target-kos STAT1,JAK2 `
  --prior-dir data\priors `
  --anchor-method boosted `
  --samples-per-ko 120 `
  --out-dir results\hard_generator_papalexi_hmpcite_adaptive_boosted_demo
```

## 输出

`formal_benchmark_predictions.csv` 会额外包含：

| column | 含义 |
|---|---|
| `boosted_feature_count` | 这个 KO 有多少 state feature 被增强 |
| `max_boost_factor` | 最大增强倍数 |
| `boosted_families` | 被增强的 feature family |
| `boost_factor_<feature>` | 某个 feature 的增强倍数 |
| `boost_family_<feature>` | 某个 feature 属于哪个增强 family |

这些列让用户可以检查模型为什么增强某些通路，而不是把 boosted 结果当作黑盒。

## 当前 Papalexi 结果

| method | AUC | direction | R2 | MAE |
|---|---:|---:|---:|---:|
| ResponseBoosted | 0.840 | 0.738 | 0.286 | 0.283 |
| ConstrainedEnsemble | 0.822 | 0.751 | 0.251 | 0.309 |
| Ridge | 0.822 | 0.751 | 0.251 | 0.309 |
| PLS | 0.807 | 0.751 | 0.220 | 0.312 |
| VKX | 0.704 | 0.462 | -0.080 | 0.386 |

Cell-level hard generator also improved:

| KO | old ensemble MAE | boosted MAE | old ensemble R2 | boosted R2 |
|---|---:|---:|---:|---:|
| JAK2 | 0.478 | 0.410 | 0.160 | 0.422 |
| STAT1 | 0.536 | 0.480 | 0.218 | 0.398 |

## 适用边界

Adaptive boost 适合“方向对但强响应幅度偏小”的场景。它不能替代更多 perturbation 数据，也不能保证所有通路都改善。

如果 boosted 提高 MAE/R2 但降低 direction cosine，说明它更重视强响应幅度，用户应结合 heatmap、ROC 曲线和 cell-state plot 一起判断。

下一步更可靠的提升方向是：接入更多同类 perturbation 数据，让 response-strength prior 从数据中学习，而不是依赖内置 family 规则。
