# Hard constraint uncertainty generator

这一步把 residual/PLS baseline 作为 hard constraint：KO 方向由 PLS 决定，生成模型不能把方向推翻，只能在这个方向附近给出不确定性范围。

## 引入的外部 perturbation 数据

外部数据来自 scPerturb/Zenodo 的标准化 h5ad 集合：

- DatlingerBock2021
- DixitRegev2016_K562_TFs_13_days

这些数据不用于学习 Papalexi 的 KO 方向，只用于估计单细胞扰动后“细胞云团围绕平均方向的波动比例”。

外部不确定性摘要：

                                      n_perturbations  median_residual_to_effect  median_target_to_control_spread
dataset                                                                                                          
DatlingerBock2021.h5ad                             12                      5.990                            0.940
DixitRegev2016_K562_TFs_13_days.h5ad               11                      9.708                            0.994

根据外部数据和保守上限，本轮允许的最大噪声比例为：`0.118`。

## 训练 KO 留一调参

候选噪声强度中，最佳值为：`0.000`。

 noise_scale  loo_mean_distribution_improvement
       0.000                             -0.116
       0.020                             -0.180
       0.039                             -0.159
       0.059                             -0.135
       0.079                             -0.118
       0.098                             -0.146
       0.118                             -0.150

## Held-out KO 测试结果

                          mean_distribution_improvement  improved_fraction
model                                                                     
hard_mean                                         0.061              0.694
hard_uncertainty_samples                          0.061              0.694

不确定性区间覆盖率：`97.2%`

区间宽度系数：`1.25`

## 图

- `results/figures/papalexi_hard_constraint_uncertainty_summary.png`
- `results/figures/papalexi_hard_constraint_uncertainty_intervals.png`
- `results/figures/papalexi_hard_constraint_uncertainty_distributions.png`

## 当前结论

当前数据下，最优噪声强度仍然非常保守。这说明方向预测比随机生成更重要：虚拟 KO 的均值应由 residual/PLS baseline 固定，生成模型现阶段更适合作为 uncertainty band，而不是自由改变细胞状态。
