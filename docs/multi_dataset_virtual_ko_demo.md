# 多数据集虚拟敲除实例：输入、输出和效果

这一步把同一套 hard-constrained residual/PLS 虚拟敲除方法放到多种数据上测试。

## 方法原则

```text
虚拟 KO cell = control cell + PLS/residual baseline 预测的 KO 方向
```

方向由系统先验和训练 KO 学到，生成模型不允许自由改变方向。对于小样本数据，当前更适合输出虚拟 KO 细胞状态和不确定性范围，而不是自由生成复杂分布。

## 实例输入和输出

输入表：

              dataset                input_modality   state_representation    ko_target  n_true_cells  n_virtual_cells                             output
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores        STAT1           120              120 virtual KO single-cell state table
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores         JAK2           120              120 virtual KO single-cell state table
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores       IFNGR2           120              120 virtual KO single-cell state table
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores         IRF1           120              120 virtual KO single-cell state table
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores     AHR+KLF1           120              120 virtual KO single-cell state table
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores  CEBPB+CEBPA            50               50 virtual KO single-cell state table
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores MAPK1+TGFBR2           120              120 virtual KO single-cell state table
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores  CBL+UBASH3B           120              120 virtual KO single-cell state table
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores          LAT           420              420 virtual KO single-cell state table
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores          LCK           420              420 virtual KO single-cell state table
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores         JUND           420              420 virtual KO single-cell state table
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores          FOS           420              420 virtual KO single-cell state table
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores         ELF1           420              420 virtual KO single-cell state table
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores        CREB1           420              420 virtual KO single-cell state table
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores         ELK1           420              420 virtual KO single-cell state table
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores        GABPA           420              420 virtual KO single-cell state table

输出包括：

- 每个目标 KO 的虚拟单细胞状态表
- 每个数据集的效果指标
- before/after UMAP 图
- UMAP 质心移动指标

## 效果摘要

              dataset                input_modality   state_representation  n_ko  n_features  mean_distribution_improvement  improved_fraction
 Datlinger CRISPR RNA single-cell RNA CRISPR screen pathway/program scores     4          14                         -0.196              0.179
Dixit Perturb-seq RNA single-cell RNA CRISPR screen pathway/program scores     4          14                         -0.114              0.232
   Norman Perturb-seq   single-cell RNA perturb-seq    gene program scores     4           5                         -0.233              0.500
  Papalexi ECCITE-seq     RNA pathway + ADT protein pathway/protein scores     4           9                          0.137              0.750

## UMAP 质心移动

              dataset  example_ko  control_to_true_umap_distance  virtual_to_true_umap_distance  umap_centroid_improvement
  Papalexi ECCITE-seq       STAT1                          2.898                          2.844                      0.019
   Norman Perturb-seq CEBPB+CEBPA                          4.458                          2.807                      0.370
 Datlinger CRISPR RNA         LAT                          0.501                          0.578                     -0.155
Dixit Perturb-seq RNA        ELF1                          0.652                          0.133                      0.795

## 图

- `results/figures/multi_dataset_virtual_ko_effect_summary.png`
- `results/figures/multi_dataset_virtual_ko_by_target.png`
- `results/figures/multi_dataset_virtual_ko_umap_examples.png`

## 怎么解释

如果 mean distribution improvement 大于 0，说明虚拟 KO 细胞比原始 control cells 更接近真实 KO cells。

如果 UMAP centroid improvement 大于 0，说明在单细胞状态空间里，虚拟 KO 的细胞云团质心向真实 KO 云团移动。

当前结果用于展示方法可以接收多种输入模态：RNA+protein、多基因 RNA perturb-seq、普通 CRISPR RNA 数据。不同数据集效果不同，这正是需要用实例输出给别人看的原因。
