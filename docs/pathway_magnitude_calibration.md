# Pathway/program score magnitude calibration

这一步没有改变 KO 方向模型，只做一件事：给 PLS/residual 预测出来的 pathway/program 变化幅度加一个非负倍率。

## 为什么要做

前一版 RNA-only 结果的主要问题不是完全预测错方向，而是预测变化幅度不稳。也就是说，模型知道细胞状态大概应该往哪里移动，但移动多远还需要校准。

## 校准规则

```text
virtual KO state = control state + alpha * predicted KO delta
```

这里的 `alpha` 限制在 `0.15 到 2.5` 之间，所以校准不会把 KO 效应方向反过来，也不会把虚拟敲除缩成“几乎没有敲除”。我们比较了三种版本：

- `uncalibrated`: 不校准，alpha = 1。
- `global_scale`: 每个数据集一个统一倍率。
- `feature_scale`: 每个 pathway/program 一个倍率。

倍率只用训练 KO 的 leave-one-KO-out 误差学习，holdout KO 的真实结果只用于最后评估。

## 结果摘要

              dataset                input_modality   state_representation  mean_distribution_improvement_before  mean_distribution_improvement_after  mean_abs_delta_error_before  mean_abs_delta_error_after calibration_method_after
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores                                -0.170                               -0.059                        0.059                       0.052            feature_scale
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores                                -0.172                               -0.024                        0.044                       0.035             global_scale
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores                                -0.766                               -0.643                        0.180                       0.171            feature_scale
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores                                 0.004                               -0.046                        0.145                       0.175            feature_scale

## 训练 KO 上选择的校准方式

calibration_method  training_delta_mae  training_direction_cosine               dataset
     feature_scale               0.056                      0.330   Papalexi ECCITE-seq
      global_scale               0.058                      0.261   Papalexi ECCITE-seq
      uncalibrated               0.062                      0.261   Papalexi ECCITE-seq
     feature_scale               0.117                      0.822    Norman Perturb-seq
      global_scale               0.118                      0.819    Norman Perturb-seq
      uncalibrated               0.119                      0.819    Norman Perturb-seq
     feature_scale               0.049                      0.437  Datlinger CRISPR RNA
      global_scale               0.050                      0.402  Datlinger CRISPR RNA
      uncalibrated               0.057                      0.402  Datlinger CRISPR RNA
      global_scale               0.058                      0.080 Dixit Perturb-seq RNA
     feature_scale               0.059                      0.085 Dixit Perturb-seq RNA
      uncalibrated               0.064                      0.080 Dixit Perturb-seq RNA

## 校准倍率

              dataset calibration_method  selected_for_dataset  mean_alpha  min_alpha  max_alpha
  Papalexi ECCITE-seq       uncalibrated                 False       1.000      1.000      1.000
  Papalexi ECCITE-seq       global_scale                 False       0.415      0.415      0.415
  Papalexi ECCITE-seq      feature_scale                  True       0.471      0.150      0.906
   Norman Perturb-seq       uncalibrated                 False       1.000      1.000      1.000
   Norman Perturb-seq       global_scale                 False       0.910      0.910      0.910
   Norman Perturb-seq      feature_scale                  True       0.874      0.813      0.967
 Datlinger CRISPR RNA       uncalibrated                 False       1.000      1.000      1.000
 Datlinger CRISPR RNA       global_scale                 False       0.340      0.340      0.340
 Datlinger CRISPR RNA      feature_scale                  True       0.356      0.150      0.658
Dixit Perturb-seq RNA       uncalibrated                 False       1.000      1.000      1.000
Dixit Perturb-seq RNA       global_scale                  True       0.150      0.150      0.150
Dixit Perturb-seq RNA      feature_scale                 False       0.253      0.150      0.902

## 图

- `results/figures/pathway_magnitude_calibration_summary.png`
- `results/figures/pathway_magnitude_calibration_delta_error_heatmap.png`
- `results/figures/pathway_magnitude_calibration_heatmap_papalexi_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_norman_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_datlinger_feature_scale.png`
- `results/figures/pathway_magnitude_calibration_heatmap_dixit_global_scale.png`

## 现在怎么看效果

如果校准后 `mean_abs_delta_error` 下降，说明 pathway/program 变化幅度更接近真实 KO。
如果 `mean_distribution_improvement` 上升，说明虚拟 KO 单细胞分布也更接近真实 KO。
如果某个数据集方向本来就错，幅度校准救不了方向错误；它只负责把已经大致正确的方向调到更合适的大小。
