# 面向用户的可视化结果说明

这套图放在 `results/user_facing_figures/`。它的目的不是展示所有模型细节，而是让第一次接触这个方法的人直接看懂三件事：

1. 用户输入什么。
2. 虚拟敲除预测了什么变化。
3. 预测效果到底好不好，以及哪些场景不能过度解释。

## 01_method_workflow.png

说明软件流程。用户输入的是普通单细胞 RNA 矩阵，或 RNA 加 protein/ADT/ATAC 等多模态矩阵；通路分数和 program score 是软件内部自动生成的状态表示，不要求用户提前准备。

## 02_single_ko_summary.png

STAT1 单基因敲除的总览图。上方四个数字是最容易读的结论：

- Direction match：越接近 1，说明预测变化方向越接近真实 KO。
- Improved features：有多少状态特征比 untreated/control 更接近真实 KO。
- ROC-AUC：识别强响应通路/蛋白的能力，越接近 1 越好。
- KO tested：本图测试的敲除基因。

下方 heatmap 对比真实 KO 和虚拟 KO 的状态变化；颜色方向一致，说明预测方向对；颜色深浅接近，说明幅度也接近。

## 03_single_ko_true_vs_virtual.png

单基因敲除的放大版 heatmap。适合回答“敲了什么基因，引起了哪些通路或蛋白变化，和真实情况差多少”。

## 04_double_ko_summary.png

CEBPB+CEBPA 双基因敲除的总览图。读法与单基因图相同。这个例子强调软件默认支持一个或两个基因敲除。

## 05_double_interaction_improvement.png

双基因敲除不是简单的两个单基因效果相加，所以这里比较了：

- Additive baseline：直接相加的基线模型。
- Interaction model：加入系统网络先验后的双基因相互作用残差模型。

结果显示 interaction model 在 R2 和 ROC-AUC 上明显更好，MAE 更低，说明它能更好处理双敲非线性效应。

## 06_reference_apply_prediction_only.png

普通 10X 单细胞数据如果没有真实 KO 标签，只能做“应用/预测”，不能在该数据内部证明准确率。图中展示的是 reference model 预测 STAT1 KO 后最可能变化的状态特征。

这类图适合展示“如果在这批普通细胞里虚拟敲掉某个基因，细胞状态可能往哪里变”；但不能报告真实 AUC 或真实准确性，除非另有 perturb-seq/CRISPR/药物扰动标签。

## 07_what_to_trust.png

给非专业用户的读图总结：

- 有真实 KO 标签时，看真实 vs 虚拟 heatmap、UMAP 和 ROC-AUC。
- 没有 KO 标签时，只能看预测的状态转换，不能说模型在该数据内部被验证。
- 双基因敲除要看 interaction model，而不是只看简单相加。
- 特征很少时 AUC 要谨慎解释，需要同时看 heatmap、R2/MAE 和生物学方向是否合理。

## 08_auc_roc_curves.png

这是 ROC 曲线版的 AUC。之前 summary 卡片里的 AUC 是一个数字摘要；这张图展示完整曲线。曲线越靠左上角，说明模型越能把“强响应通路/蛋白”和“弱响应特征”区分开。

## 09_multimodal_multi_dataset_summary.png

展示多数据集结果。Papalexi 是当前真正的多模态例子，即 RNA pathway score 加 ADT/protein；Norman、Datlinger、Dixit 主要是 RNA-only 或 RNA-derived 状态表示。当前结果显示，多模态例子的稳定性更好，RNA-only 数据在小样本和跨数据设置下波动更大。

## 10_r2_mae_intuitive.png

把 R2 和 MAE 改成更直观的读法：

- R2：模型解释真实变化模式的程度，越高越好。
- MAE：预测误差，越低越好。
- Error reduction：相互作用模型相比简单相加模型把平均误差降低了多少。

## 11_umap_before_after_examples.png

把单敲、多敲和多数据集的 UMAP 放在一起。UMAP 主要用来看虚拟敲除前后细胞状态是否发生可见移动；但 UMAP 是可视化，不是准确率证明，准确性仍要看真实 KO 标签下的 heatmap、ROC-AUC、R2/MAE。

## 12_other_double_ko_combos.png

展示 Norman 数据中多个其它双基因敲除组合。每一行是一个双敲组合，每两列比较简单相加模型和加入相互作用先验后的模型误差。颜色越浅、数字越小越好。

## 13_10x_multimodal_single_double_outputs.png

这张图补充普通 10X 和多模态输入的实际输出。

左边说明普通细胞矩阵输入后，软件会自动派生 RNA pathway score；如果 h5ad 里有 ADT/protein/ATAC 这类额外模态，也会作为额外状态特征进入模型。

右边展示 RNA+ADT 多模态输入下，同时预测一个基因敲除（STAT1）和两个基因敲除（STAT1+JAK2）后的状态变化。注意这属于 reference model application：可以输出预测变化和虚拟细胞，但如果输入数据没有真实双敲标签，就不能在该数据内部计算真实 AUC/R2/MAE。

## 14_hmpcite_multimodal_doubleko_summary.png

真实多模态 double-KO 结果总览。这里接入的是 HMPCITE-seq GSE243244，输入包含 RNA、ADT 和 guide-derived KO 标签。

这张图的重点是说明：我们已经不只是 RNA-only，也不只是单敲。Cebpb+Med12 双敲在 RNA pathway + ADT 状态空间里可以做真实评估，heatmap、UMAP 和 ROC 曲线都能看。

## 15_rna_adt_atac_extension_summary.png

RNA+ADT 和 ATAC 扩展结果的总览图。

- RNA+ADT double-KO：ROC-AUC 0.98，方向一致性 0.98，说明真实多模态双敲效果较好。
- ATAC gene activity：ROC-AUC 0.64，方向一致性 0.51，说明 ATAC 层已经支持，但调控层预测更难，需要继续加入 motif/TF-target 先验。
- 图底部说明三模态输入接口：`--extra-obsm protein:protein,atac:atac`。

## 16_modality_extension_result_gallery.png

多模态扩展结果画廊。每个数据都展示三类用户最关心的图：

- heatmap：敲了什么基因，哪些通路/状态变化，虚拟 KO 和真实 KO 差多少。
- UMAP：敲除前后单细胞状态是否发生可见移动。
- ROC 曲线：AUC 对应的完整曲线，而不是只有一个柱状图或数字。

上排是 RNA+ADT double-KO，下排是 ATAC gene activity。用户可以直接比较哪一类输入更稳定、哪一类仍然困难。

## 17_multimodal_input_visualization_matrix.png

输入和输出边界表。它回答“用户给什么文件，软件能生成什么可视化”。

核心规则是：AUC/R2/MAE 需要真实 KO 标签；没有 KO 标签的普通 10X 或 Multiome 数据可以做虚拟状态变化展示，但不能在该数据内部证明真实准确性。

## 18_atac_chromvar_prior_ablation.png

ATAC/chromVAR 优化消融图。它比较了三种 KDM6A 虚拟 KO 输入：

- ATAC gene activity only。
- ATAC gene activity + 全量 chromVAR motif activity。
- ATAC gene activity + 方差最高的 100 个 chromVAR motif。

这张图的重点结论是：motif/TF 信息不是越多越好。全量 chromVAR 会引入噪声并降低 AUC；筛选 top100 motif 后方向一致性和 MAE 改善，但 AUC 仍低于 gene activity only。也就是说 ATAC 调控层已经能接入，但需要更强的 motif/TF-target 筛选和加权。

## 19_atac_peak_level_visualization.png

真正的 ATAC peak-level 可视化。之前只展示 gene activity 和 chromVAR，所以看不到“峰图”；现在已经把 `peak_bc` peak matrix 接入，并把 KDM6A 附近 peaks 和高可变 peaks 作为 `obsm["peak"]` 状态特征进入模型。

这张图展示四件事：

- 加入 selected peaks 后，模型 AUC 为 0.62，方向一致性为 0.62，MAE 约 0.05。
- KDM6A locus 的 promoter、intronic、distal peaks 可以直接看真实 KO delta 和虚拟 KO delta。
- peak delta heatmap 显示哪些峰方向对、哪些峰误差大。
- 单细胞 peak accessibility distribution 显示 control、virtual KO、true KO 在某个 promoter peak 上的分布差异。

结论是：现在 ATAC 不只是 gene activity proxy，已经能显示真实 peak-level 结果；但 peak-level 单细胞分布仍然比 RNA/ADT 更难预测。

## 20_atac_weighted_prior_feature_selection_summary.png

ATAC 优化总览图。它比较四个版本：

- gene activity baseline。
- top100 chromVAR variance selection。
- weighted TF/motif prior + hybrid chromVAR selection。
- weighted TF/motif prior + hybrid selection + locus-aware peaks。

这张图的核心结论是：加权 prior 和 hybrid feature selection 明显改善 ATAC/chromVAR 的整体信号；hybrid chromVAR 版本的方向一致性、MAE 和分布改进最好；加入 locus-aware peaks 后 AUC 最高，说明 peak-level 特征更有助于强响应排序，但峰层面的单细胞分布仍然难。
