# SCI 论文草稿：小样本多模态虚拟敲除

## 写作定位

论文类型：算法/方法论文。

一句话论点：

在小样本单细胞扰动研究中，我们提出一种网络先验约束的多模态状态空间虚拟敲除方法，通过 pathway/program/protein/ATAC 状态分数和 hard-constrained residual/PLS baseline，在无需大规模预训练的情况下预测单基因和双基因敲除后的细胞状态变化；该方法在 RNA+ADT double-KO 数据上表现较好，在 ATAC/chromVAR 调控层上可运行但仍存在明显挑战。

English one-sentence argument:

In small-sample single-cell perturbation studies, we present a network-prior-constrained multimodal state-space virtual knockout framework that predicts single- and double-gene knockout-induced cellular state shifts using pathway, protein and regulatory state scores with a hard-constrained residual/PLS baseline, achieving strong performance in RNA+ADT double-knockout data while revealing the remaining difficulty of ATAC/chromVAR regulatory-state prediction.

## 题目备选

中文：

1. 面向小样本多模态单细胞扰动数据的网络先验约束虚拟敲除方法
2. 基于可解释状态空间的小样本单细胞虚拟敲除建模
3. 多模态状态分数与网络先验约束的单细胞虚拟敲除框架

English:

1. Network-prior-constrained virtual knockout in interpretable multimodal single-cell state spaces
2. Small-sample virtual knockout modelling with pathway, protein and regulatory state constraints
3. A hard-constrained residual baseline for multimodal single-cell virtual knockout prediction

## 摘要（中文）

虚拟敲除旨在预测基因敲除后细胞状态的变化，可用于在真实 CRISPR 或 perturb-seq 实验前筛选候选基因和基因组合。然而，现有生成式或大规模预训练方法通常依赖大量扰动样本和较高计算资源，在小样本、多模态数据场景中难以直接应用。本文提出一种面向小样本单细胞扰动数据的多模态虚拟敲除框架。该框架将原始 RNA、ADT/protein、ATAC gene activity 和 chromVAR motif activity 自动转换为可解释的 pathway、protein 和 regulatory state scores，并使用 Reactome、MSigDB、TF-target、PPI 和 motif/TF 先验构建 KO 方向约束。模型采用 residual/PLS baseline 作为 hard constraint，用于预测单基因和双基因敲除后的状态变化；后续 cell-level 生成模型可在该方向附近学习不确定性，而不是从零学习自由分布。我们在 HMPCITE-seq RNA+ADT double-KO 数据上验证了方法，Cebpb+Med12 双敲达到 ROC-AUC 0.978、方向一致性 0.976，且 87.5% 状态特征较 control 更接近真实 KO。双敲 interaction residual 模型相比简单相加模型将 MAE 从 0.195 降至 0.113，并将 R2 从 -0.334 提高到 0.507。在 scPerturb ATAC K562 数据中，KDM6A KO 的 gene activity 层达到 ROC-AUC 0.641；加入全量 chromVAR motif 后 AUC 降至 0.552，而筛选 top100 motif 后方向一致性提高至 0.618，但 AUC 仍为 0.586，提示调控层虚拟敲除仍需更强的 motif/TF-target 筛选和加权。本文进一步明确了无 perturbation 标签的 DOGMA/TEA-seq 等三模态数据只能用于状态转换和 reference model application，不能用于真实准确性验证。总体而言，该方法为小样本多模态扰动研究提供了一种可解释、可复用、边界明确的虚拟敲除基线。

## Abstract (English)

Virtual knockout aims to predict cellular state changes after gene perturbation and can prioritize candidate genes or gene combinations before experimental CRISPR or perturb-seq assays. Existing generative or large-scale pretrained approaches often require many perturbation samples and substantial computational resources, limiting their use in small-sample multimodal studies. Here we present a multimodal virtual knockout framework for small-sample single-cell perturbation data. The framework converts raw RNA, ADT/protein, ATAC gene activity and chromVAR motif activity into interpretable pathway, protein and regulatory state scores, and constrains knockout directions using Reactome, MSigDB, TF-target, PPI and motif/TF priors. Instead of training an unconstrained generator from scratch, the method uses a residual/PLS baseline as a hard constraint to predict single- and double-knockout state shifts; downstream cell-level generative models can then learn uncertainty around this constrained direction. In HMPCITE-seq RNA+ADT double-knockout data, prediction of Cebpb+Med12 achieved a ROC-AUC of 0.978, a direction cosine of 0.976 and moved 87.5% of state features closer to the real knockout than control cells. A double-knockout interaction residual model reduced MAE from 0.195 to 0.113 and increased R2 from -0.334 to 0.507 compared with an additive baseline. In scPerturb ATAC K562 data, KDM6A knockout prediction using gene activity reached a ROC-AUC of 0.641. Adding all chromVAR motif features reduced AUC to 0.552, whereas selecting the top 100 variable motifs improved direction cosine to 0.618 but retained a lower AUC of 0.586, indicating that regulatory-state virtual knockout remains challenging and requires stronger motif/TF-target filtering and weighting. We further distinguish labelled perturbation datasets, which support accuracy evaluation, from unlabeled DOGMA-seq or TEA-seq datasets, which support only state conversion and reference-model application. Overall, this work provides an interpretable, reusable and explicitly bounded virtual knockout baseline for small-sample multimodal perturbation studies.

## Introduction（中文草稿）

单细胞 CRISPR 扰动技术可以直接观测基因敲除或激活后的细胞状态变化，但真实实验通常受到样本量、组合空间、测序成本和实验条件的限制。尤其是在多基因敲除中，可能的组合数量随候选基因数快速增加，完全依赖实验筛选并不现实。因此，虚拟敲除成为一个有实际价值的问题：给定已有扰动数据或普通单细胞数据，预测敲除一个或两个目标基因后细胞状态可能如何变化。

现有虚拟扰动方法大体可以分为三类。第一类是基于线性或低秩模型的扰动效应预测，优点是稳定、可解释、适合小样本，但对复杂非线性分布的刻画能力有限。第二类是基于深度生成模型的方法，如 VAE、normalizing flow、diffusion 或条件生成模型，优点是能够模拟 cell-level 分布，但通常需要较多训练扰动和更高算力。第三类是大规模预训练或跨数据迁移模型，它们能够利用大量已有扰动资源，但对数据规模、批次对齐和预训练资源要求较高。对于许多实际实验室场景，尤其是小样本、多模态、候选基因有限的项目，直接训练复杂生成模型并不一定是最稳妥的选择。

本文的出发点是：小样本虚拟敲除不应首先追求自由生成完整高维表达矩阵，而应先在可解释状态空间中学习可靠的 KO 方向和幅度。多模态数据在这里具有特殊价值。RNA 表达可以被压缩为 pathway 或 program score，ADT/protein 可以直接提供表面蛋白状态，ATAC gene activity 和 chromVAR motif activity 可以提供调控层信息。与原始高维基因表达相比，这些状态分数更低维、更可解释，也更适合在小样本下加入系统先验。

基于这一思想，我们提出一个网络先验约束的多模态虚拟敲除框架。该框架接受普通单细胞 RNA 矩阵、RNA+ADT、ATAC gene activity、chromVAR motif activity 或未来 RNA+ADT+ATAC 三模态 perturbation 数据作为输入，自动构建 pathway/protein/regulatory state scores，并以 Reactome、MSigDB、TF-target、PPI 和 motif/TF prior 约束 KO 方向。模型核心不是自由生成器，而是 hard-constrained residual/PLS baseline。该 baseline 可以单独作为小样本方法使用，也可以作为后续 VAE、flow matching 或 diffusion 的约束骨架。

## Introduction (English Draft)

Single-cell CRISPR perturbation assays directly measure cellular state changes after gene knockout or activation, but experimental screens are constrained by sample size, combinatorial space, sequencing cost and assay conditions. These constraints become particularly severe for multi-gene knockouts, where the number of possible combinations grows rapidly with the number of candidate genes. Virtual knockout is therefore a practical modelling problem: given existing perturbation data or ordinary single-cell profiles, can we predict how the cellular state would change after knocking out one or two target genes?

Existing virtual perturbation methods can be broadly grouped into three categories. Linear or low-rank perturbation models are stable, interpretable and suitable for small samples, but have limited capacity to represent complex nonlinear cell-state distributions. Deep generative approaches, including VAE, flow-based and diffusion models, can model cell-level distributions but typically require many perturbation labels and substantial compute. Large-scale pretrained or cross-dataset transfer models can exploit external perturbation resources, but depend on large training corpora, batch alignment and pretrained infrastructure. For many practical laboratory settings, especially small-sample multimodal studies with a limited number of candidate genes, directly training an unconstrained generator may not be the most reliable starting point.

Our premise is that small-sample virtual knockout should first learn a reliable perturbation direction and magnitude in an interpretable state space, rather than immediately generating full high-dimensional expression profiles. Multimodal data are particularly useful in this setting. RNA expression can be summarized into pathway or program scores, ADT/protein measurements directly capture surface protein states, and ATAC gene activity or chromVAR motif activity provide regulatory-state information. Compared with raw gene expression, these state scores are lower-dimensional, more interpretable and more amenable to biological priors under small-sample conditions.

We therefore developed a network-prior-constrained multimodal virtual knockout framework. The framework accepts ordinary single-cell RNA matrices, RNA+ADT data, ATAC gene activity, chromVAR motif activity or future RNA+ADT+ATAC perturbation datasets; derives pathway, protein and regulatory state scores; and constrains knockout directions using Reactome, MSigDB, TF-target, PPI and motif/TF priors. The core model is not a free-form generator but a hard-constrained residual/PLS baseline. This baseline can be used directly as a small-sample method and can also serve as the constrained backbone for subsequent VAE, flow-matching or diffusion-based cell-level generators.

## Methods（中文草稿）

### 输入和状态表示

输入数据为 AnnData h5ad 或 CSV。对于 h5ad，`adata.X` 存放 RNA 表达矩阵或 ATAC gene activity 矩阵，`adata.obs["ko_target"]` 存放 control、single-KO 或 double-KO 标签。额外模态通过 `adata.obsm` 输入，例如 `adata.obsm["protein"]` 表示 ADT/CITE-seq protein，`adata.obsm["atac"]` 表示 ATAC-derived state features，`adata.obsm["chromvar"]` 表示 chromVAR motif activity。软件通过 `--extra-obsm protein:protein,atac:atac,chromvar:tf` 接口统一读取这些模态。

RNA 或 gene activity 主矩阵首先被转换为 pathway/program scores。软件从 Reactome、MSigDB Hallmark、TF-target 和 PPI 先验中选择与输入基因集重叠的 term，并计算每个细胞在对应 gene set 上的平均表达或 gene activity，随后进行标准化。ADT、ATAC 和 chromVAR 等额外模态作为额外状态特征拼接到状态表中。对于高维 chromVAR 或 ATAC feature matrix，用户可以使用 `--max-extra-features-per-obsm` 进行方差筛选，避免将大量弱噪声 motif 直接加入模型。

### KO 先验向量

每个 KO 标签被解析为一个或两个基因。模型根据 KO 基因与 Reactome、MSigDB、TF-target、PPI 和 motif/TF prior 的重叠构建 prior vector。为了增强 ATAC 调控层先验，TF-target 与 motif/TF 类 term 不仅使用 target genes，也将 term 名称中的 TF/gene 纳入匹配。这使模型能够更直接识别 TF 或 motif 对 KO 方向的贡献。

### Hard-constrained residual/PLS baseline

对每个训练 KO，模型先计算其相对于 control cells 的平均状态变化。随后以 KO prior vector 为输入、KO delta 为输出训练 PLSRegression。对于 held-out KO，模型预测状态 delta，并将该 delta 加到 sampled control cells 上，生成 virtual KO cells。该过程构成 hard constraint：后续生成模型只能在预测方向附近学习不确定性，而不是完全自由生成。

### 双敲 interaction residual

对于 double-KO，简单 additive baseline 将两个 single-KO delta 相加。interaction residual 模型进一步学习 double-KO truth 与 additive prediction 之间的残差，并使用网络先验描述基因对之间的关系。该模型用于捕捉双敲非线性效应。

### 评估指标

在有真实 KO 标签时，模型报告：

- ROC-AUC：识别强响应状态特征的能力。
- Direction cosine：预测 KO delta 与真实 KO delta 的方向一致性。
- MAE：预测 delta 与真实 delta 的平均绝对误差。
- Distribution improvement：virtual KO cells 相比 control cells 是否更接近真实 KO cells。
- Heatmap 和 UMAP：分别展示状态特征变化与单细胞状态移动。

在无 KO 标签数据中，软件只输出状态转换和 reference model application 结果，不报告真实准确性指标。

## Methods (English Draft)

### Input and State Representation

The framework accepts AnnData h5ad or CSV files. For h5ad input, `adata.X` stores the RNA expression matrix or ATAC gene activity matrix, and `adata.obs["ko_target"]` stores control, single-knockout or double-knockout labels. Additional modalities are provided through `adata.obsm`, for example `adata.obsm["protein"]` for ADT/CITE-seq proteins, `adata.obsm["atac"]` for ATAC-derived state features and `adata.obsm["chromvar"]` for chromVAR motif activity. These modalities are read through a unified interface, `--extra-obsm protein:protein,atac:atac,chromvar:tf`.

The RNA or gene-activity matrix is first converted into pathway or program scores. Terms overlapping the input gene set are selected from Reactome, MSigDB Hallmark, TF-target and PPI priors. For each cell, the mean expression or gene activity over each gene set is computed and standardized. ADT, ATAC and chromVAR modalities are appended as additional state features. For high-dimensional chromVAR or ATAC feature matrices, the optional argument `--max-extra-features-per-obsm` selects the most variable features and reduces the effect of weak motif-level noise.

### Knockout Prior Vector

Each knockout label is parsed into one or two genes. The model constructs a prior vector from the overlap between knockout genes and Reactome, MSigDB, TF-target, PPI and motif/TF terms. To strengthen regulatory priors for ATAC data, TF-target and motif/TF terms use both their target genes and the TF or gene encoded in the term name. This allows the model to more directly recognize TF- or motif-related knockout effects.

### Hard-constrained Residual/PLS Baseline

For each training knockout, the model computes its mean state change relative to control cells. PLSRegression is trained to map the knockout prior vector to the observed knockout delta. For a held-out knockout, the model predicts a state delta and adds this delta to sampled control cells to generate virtual knockout cells. This procedure acts as a hard constraint: downstream generative models can learn uncertainty around the predicted direction rather than generating unconstrained states.

### Double-Knockout Interaction Residual

For double knockouts, the additive baseline sums the two single-knockout deltas. The interaction residual model learns the residual between the true double-knockout effect and this additive prediction using network-prior features describing the gene pair. This module captures nonlinear double-knockout effects.

### Evaluation

When true knockout labels are available, the model reports ROC-AUC for strong-response feature ranking, direction cosine for delta orientation, MAE for magnitude error, distribution improvement for cell-level distribution shift, and heatmap/UMAP visualizations. For unlabeled datasets, including ordinary 10X, DOGMA-seq or TEA-seq without perturbation labels, the software reports only predicted state shifts or reference-model application results and does not compute accuracy metrics.

## Results（中文草稿）

### RNA+ADT double-KO 数据验证

我们首先在 HMPCITE-seq GSE243244 上评估模型。该数据包含 RNA、ADT 和 guide-derived perturbation labels，因此可以进行真实准确性验证。对 Cebpb+Med12 双敲，模型达到 ROC-AUC 0.978、方向一致性 0.976、平均分布改进 0.548，并使 87.5% 状态特征相比 control 更接近真实 KO。heatmap 显示 interferon、TGFB、STAT3 和 CD86 等状态特征的变化方向与真实 KO 基本一致，UMAP 显示 virtual KO cells 从 control state 向真实 KO state 移动。

### 双敲非线性效应

为了处理双敲非线性，我们比较了 additive baseline 和 interaction residual model。在 HMPCITE double-KO 评估中，additive baseline 的 MAE 为 0.195、R2 为 -0.334、ROC-AUC 为 0.763；interaction residual model 将 MAE 降至 0.113、R2 提高至 0.507、ROC-AUC 为 0.768。这表明 interaction residual 主要改善了幅度和模式拟合，而对强响应排序的改善较小。

### ATAC/gene activity 层评估

在 scPerturb ATAC K562 数据中，我们使用 gene activity 主矩阵预测 KDM6A KO。模型达到 ROC-AUC 0.641、方向一致性 0.513、MAE 0.053，但分布改进接近 0。这说明 ATAC 层可以接入和评估，但调控层响应比 RNA+ADT 更弱、更稀疏，也更难预测。

### chromVAR motif activity 与 peak-level 输入消融

我们进一步将 chromVAR motif activity 作为 `--extra-obsm chromvar:tf` 输入。全量 2174 个 motif feature 使 ROC-AUC 降至 0.552、方向一致性降至 0.432。使用 `--max-extra-features-per-obsm 100` 筛选 top-variable motifs 后，方向一致性升至 0.618，MAE 降至 0.050，improved feature fraction 升至 61.5%，但 ROC-AUC 仍为 0.586。该结果说明 motif 层不是简单加入越多越好；需要 motif/TF-target 筛选、加权或更强先验建模。

进一步地，我们将 TF-target/motif prior 从简单 overlap 升级为加权 prior，并加入 direct TF hit 与 term coverage 权重；同时将 chromVAR/ATAC feature selection 从方差筛选扩展为 `variance`、`ko_effect` 和 `hybrid`。在不加入 peak 的情况下，weighted prior + hybrid chromVAR 使 KDM6A KO 的 ROC-AUC 达到 0.645、方向一致性达到 0.728、MAE 为 0.048，整体分布改进为 0.113。随后，我们从 `peak_bc` 接入细胞×peak matrix，将 KDM6A locus peaks、`markerpeak_target` peaks 和全局稳定 peaks 作为 `obsm["peak"]` 输入模型。加入 locus-aware peaks 后，ROC-AUC 进一步提高至 0.658，方向一致性为 0.684。该结果说明 peak-level 输入更有利于强响应特征排序，但分布改进仍为负值，提示 peak-level 单细胞分布预测仍然困难。该实验使 ATAC 结果不再只依赖 gene activity 或 motif proxy，而能直接展示真实 peak-level delta 和峰图。

### 无标签三模态数据的应用边界

DOGMA-seq 和 TEA-seq 可以提供 RNA、protein/epitope 和 ATAC 三模态信息，但如果没有 perturbation labels，只能用于状态分数生成和 reference model application。此时可视化可以展示 virtual KO 后预测状态移动，但不能报告真实 AUC、R2 或 MAE。

## Results (English Draft)

### Validation in RNA+ADT Double-Knockout Data

We first evaluated the framework on HMPCITE-seq GSE243244, which contains RNA, ADT and guide-derived perturbation labels and therefore supports direct accuracy evaluation. For Cebpb+Med12 double knockout, the model achieved a ROC-AUC of 0.978, a direction cosine of 0.976, a mean distribution improvement of 0.548 and moved 87.5% of state features closer to the real knockout than control cells. Heatmaps showed concordant changes in interferon, TGFB, STAT3 and CD86-related state features, and UMAP visualization showed virtual knockout cells moving from the control state toward the real knockout state.

### Double-Knockout Nonlinearity

To model nonlinear double-knockout effects, we compared an additive baseline with an interaction residual model. In HMPCITE double-knockout evaluation, the additive baseline obtained MAE 0.195, R2 -0.334 and ROC-AUC 0.763, whereas the interaction residual model reduced MAE to 0.113, increased R2 to 0.507 and achieved ROC-AUC 0.768. Thus, the interaction residual model mainly improved magnitude and pattern fitting, with a smaller effect on strong-response ranking.

### ATAC/Gene-Activity Evaluation

In scPerturb ATAC K562 data, we used the gene-activity matrix to predict KDM6A knockout. The model reached a ROC-AUC of 0.641, a direction cosine of 0.513 and MAE 0.053, but distribution improvement was close to zero. This indicates that ATAC-layer virtual knockout is feasible but more difficult than RNA+ADT prediction, likely because regulatory responses are weaker, sparser and more indirect.

### chromVAR Motif-Activity and Peak-Level Ablation

We then added chromVAR motif activity through `--extra-obsm chromvar:tf`. Adding all 2174 motif features reduced ROC-AUC to 0.552 and direction cosine to 0.432. Selecting the top 100 variable motifs with `--max-extra-features-per-obsm 100` improved direction cosine to 0.618, reduced MAE to 0.050 and increased the improved feature fraction to 61.5%, but ROC-AUC remained 0.586. These results show that motif-level information should be filtered or weighted rather than added blindly.

We further upgraded the TF-target/motif prior from simple overlap to a weighted prior with library weights, direct TF-hit bonuses and term-coverage weights, and extended chromVAR/ATAC feature selection from variance filtering to `variance`, `ko_effect` and `hybrid` modes. Without peak features, weighted priors with hybrid chromVAR selection achieved ROC-AUC 0.645, direction cosine 0.728, MAE 0.048 and mean distribution improvement 0.113 for KDM6A knockout. We then incorporated the cell-by-peak matrix from `peak_bc` and used KDM6A locus peaks, `markerpeak_target` peaks and globally stable peaks as `obsm["peak"]` model inputs. Adding locus-aware peak features increased ROC-AUC to 0.658 with direction cosine 0.684. Thus, peak-level input improved strong-response ranking but did not solve peak-level cell-distribution modelling. This experiment moves the ATAC analysis beyond gene-activity or motif proxies and enables direct visualization of true and virtual peak-level deltas.

### Boundary for Unlabeled Trimodal Data

DOGMA-seq and TEA-seq can provide RNA, protein/epitope and ATAC measurements, but if perturbation labels are absent they support only state scoring and reference-model application. In this setting, visualizations can show predicted virtual knockout state shifts, but true AUC, R2 and MAE cannot be reported.

## Discussion（中文草稿）

本研究支持一个较稳妥的观点：小样本虚拟敲除不应直接从自由生成模型开始，而应首先建立一个可解释、可校准、可验证的状态空间 baseline。RNA+ADT double-KO 结果表明，当数据同时具有真实 perturbation 标签和互补模态时，网络先验约束的 residual/PLS baseline 能够较好恢复 KO 方向和主要状态变化。双敲 interaction residual 的改善进一步说明，多基因敲除需要显式处理非线性相互作用。

与此同时，ATAC/chromVAR 结果也指出了方法边界。调控层信号更稀疏，motif activity 更容易受到噪声、motif redundancy 和间接调控影响。全量 chromVAR 降低 AUC，说明多模态不是简单“越多越好”；模态信息需要被选择、加权和约束。top100 chromVAR 提高方向一致性但没有提高 AUC，提示 motif 层可能更有助于整体方向和幅度，而不一定改善强响应特征排序。

本文另一个重要边界是标签依赖。没有 perturbation 标签的普通 10X、DOGMA-seq 或 TEA-seq 数据可以用于 reference application，但不能在数据内部证明模型准确性。因此，软件需要在输出中明确 analysis mode，避免用户把预测展示误读为真实验证。

未来工作应沿三条路线推进。第一，构建 ATAC-specific motif/TF-target 加权先验，避免全量 motif 噪声。第二，引入更多带标签 perturbation 数据进行跨数据 reference training。第三，在 hard-constrained residual baseline 上接入 VAE、flow matching 或 diffusion，使生成模型只学习预测方向附近的 cell-level 不确定性。

## Discussion (English Draft)

This study supports a conservative strategy for small-sample virtual knockout: instead of starting from an unconstrained generator, one should first establish an interpretable, calibratable and evaluable state-space baseline. The RNA+ADT double-knockout results show that, when real perturbation labels and complementary modalities are available, a network-prior-constrained residual/PLS baseline can recover knockout direction and major state changes. The improvement from the interaction residual model further indicates that multi-gene knockouts require explicit modelling of nonlinear interactions.

The ATAC/chromVAR results also define the boundary of the current approach. Regulatory-layer signals are sparser and motif activity is affected by noise, motif redundancy and indirect regulation. Adding all chromVAR features reduced AUC, showing that multimodality is not automatically beneficial; modality-specific information must be selected, weighted and constrained. Selecting the top 100 chromVAR motifs improved direction and magnitude metrics but did not improve AUC, suggesting that motif-level information may better support overall direction than strong-response feature ranking.

A second boundary is label dependence. Ordinary 10X, DOGMA-seq or TEA-seq datasets without perturbation labels can be used for reference-model application, but they cannot internally validate model accuracy. The software therefore records the analysis mode in each output directory to prevent prediction-only visualizations from being interpreted as evidence of accuracy.

Future work should proceed along three directions. First, ATAC-specific motif/TF-target priors should be weighted and filtered to reduce motif-level noise. Second, more labelled perturbation datasets should be incorporated for cross-dataset reference training. Third, VAE, flow-matching or diffusion modules should be built on top of the hard-constrained residual baseline so that cell-level generators learn uncertainty around a predicted knockout direction rather than unconstrained distributions.

## Figure Plan

Figure 1. Method workflow: raw single-cell/multimodal input, state scoring, prior-constrained KO delta, virtual KO cells, evaluation/application mode.

Figure 2. RNA+ADT double-KO validation: HMPCITE Cebpb+Med12 summary, heatmap, UMAP and ROC curve.

Figure 3. Double-KO interaction residual: additive baseline vs interaction model, MAE/R2/AUC comparison.

Figure 4. Multimodal input boundary: RNA-only, RNA+ADT, ATAC, RNA+ADT+ATAC with labels, and unlabeled DOGMA/TEA-seq.

Figure 5. ATAC/chromVAR ablation: gene activity only, all chromVAR, top100 chromVAR; metrics and representative ROC/heatmap.

## Limitations

1. Current validation is strongest for RNA+ADT perturbation data and weaker for ATAC-only regulatory-state prediction.
2. True public RNA+ADT+ATAC perturbation datasets with clear KO labels remain harder to identify than RNA+ADT or RNA+ATAC perturbation datasets; the current software interface is ready for such data but the present report does not claim a completed accuracy benchmark on a true labelled trimodal dataset.
3. chromVAR motif activity requires feature selection or weighting; full motif matrices can reduce AUC.
4. Unlabeled multiome or trimodal data cannot provide internal accuracy metrics.
5. The current model is a hard-constrained baseline rather than a free cell-level diffusion/VAE generator.

## Key Evidence Table

| Dataset / setting | Modality | KO task | Main result | Interpretation |
|---|---|---|---|---|
| HMPCITE GSE243244 | RNA+ADT+GDO labels | Cebpb+Med12 | AUC 0.978, direction 0.976 | Strong real labelled multimodal double-KO validation |
| HMPCITE double-KO interaction | RNA+ADT state deltas | 55 double-KO combinations | MAE 0.113 vs 0.195 additive; R2 0.507 vs -0.334 | Interaction residual improves double-KO magnitude/pattern |
| scPerturb ATAC | gene activity | KDM6A | AUC 0.641, direction 0.513 | ATAC layer is supported but difficult |
| scPerturb ATAC + all chromVAR | gene activity + 2174 motifs | KDM6A | AUC 0.552, direction 0.432 | Full motif matrix adds noise |
| scPerturb ATAC + top100 chromVAR | gene activity + selected motifs | KDM6A | AUC 0.586, direction 0.618 | Motif selection helps direction but not AUC |
| scPerturb ATAC + weighted hybrid chromVAR | gene activity + selected motifs | KDM6A | AUC 0.645, direction 0.728, MAE 0.048 | Weighted prior and hybrid feature selection improve overall ATAC prediction |
| scPerturb ATAC + weighted hybrid + locus-aware peaks | gene activity + motifs + peak accessibility | KDM6A | AUC 0.658, direction 0.684, MAE 0.057 | Peak-level input gives best AUC and enables peak plots, but distribution remains hard |
| DOGMA/TEA-seq without labels | RNA+protein+ATAC | prediction only | no internal AUC/R2/MAE | Compatible input, not accuracy validation |
