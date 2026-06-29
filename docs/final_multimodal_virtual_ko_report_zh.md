# 小样本多模态虚拟敲除方法优化报告

## 1. 一句话结论

我们开发并迭代了一个面向小样本单细胞扰动数据的虚拟敲除框架：用户输入普通单细胞矩阵或多模态 h5ad，软件自动转换为可解释的 pathway/program/protein/ATAC 状态分数，并用网络先验约束的 residual/PLS baseline 预测单基因或双基因敲除后的细胞状态变化。

当前结果表明：

- 在真实带标签 RNA+ADT double-KO 数据上，方法表现较好。
- 在 ATAC/gene activity 数据上，方法可以运行和评估，但调控层预测明显更难。
- chromVAR motif activity 已经可以通过 `--extra-obsm chromvar:tf` 接入；全量 motif 会引入噪声，筛选后的 top motif 对方向和误差有帮助，但 AUC 仍未超过 gene activity only。
- 对无 perturbation 标签的 DOGMA/TEA-seq 这类三模态数据，只能做三模态输入兼容和 reference application，不能报告真实 AUC/R2/MAE。

## 2. 为什么做虚拟敲除

真实 CRISPR/perturb-seq 实验成本高、周期长，而且很多组合基因敲除不可能穷尽实验。虚拟敲除的目标不是替代实验，而是在实验前给出一个可解释的候选筛选结果：

1. 敲除某个基因后，细胞状态可能向哪里移动。
2. 哪些 pathway、program、蛋白或调控 motif 变化最大。
3. 单敲和双敲是否可能产生不同于简单相加的非线性效应。
4. 在小样本或资源有限场景下，优先测试哪些基因组合。

## 3. 当前方法选择

我们没有直接从零训练大模型式 diffusion/VAE，而是采用更稳的“小样本 hard-constrained baseline”：

```text
输入矩阵
  -> 自动生成状态分数
  -> 网络先验约束 KO 方向
  -> residual/PLS baseline 预测 KO delta
  -> 用真实 KO 标签评估，或对无标签数据做 reference application
```

这样设计的原因是：

- 小样本下，自由生成模型容易学到噪声。
- pathway/program/protein/ATAC 状态比原始高维基因表达更容易解释。
- 网络先验可以把 KO gene 与 Reactome/MSigDB/TF-target/PPI/motif 关系连接起来。
- residual baseline 可以作为后续 VAE、flow matching、diffusion 的 hard constraint，而不是 soft guidance。

### 3.1 当前模型学过多少 KO

当前原型接入过多个真实 perturbation 数据集，但这不等于一个模型一次性学过全基因组 KO。更准确的说法是：

```text
本地实验已经接入百级别扰动基因覆盖；
每一次 reference model 训练只在具体数据集和具体状态空间中学习 KO direction。
```

当前已接入的数据覆盖：

| 数据集 | 模态 | 真实扰动标签 |
|---|---|---:|
| Papalexi ECCITE-seq | RNA + ADT protein | 25 个单基因 KO |
| Norman Perturb-seq | RNA | 25 个单基因 KO + 52 个双基因 KO |
| HMPCITE-seq | RNA + ADT + guide-derived labels | 11 个基因组成的 66 个单敲/双敲条件 |
| Liscovitch ATAC / scPerturb | ATAC gene activity + chromVAR + peak features | 21 个单基因 KO |
| Datlinger CRISPR | RNA | 约 20 个基因，40 个 guide/扰动标签 |
| Dixit TF perturbation | RNA | 10 个 TF KO |

合并看，当前接入过大约 129 个基因名。由于不同数据集存在物种、大小写、命名和状态空间差异，这个数字应理解为“数据覆盖范围”，不是一个统一全基因组预训练模型的训练基因数。

### 3.2 没学过的基因怎么预测

对训练集中没出现过的基因，模型不会用 one-hot 记忆，而是构建系统先验向量：

```text
target gene
-> Reactome/MSigDB pathway membership
-> TF-target relationship
-> PPI neighborhood
-> motif / chromVAR / ATAC regulatory relationship
-> KO prior vector
```

模型从已见 KO 中学习：

```text
KO prior vector -> pathway / protein / ATAC state delta
```

所以未见基因的虚拟敲除过程是：

```text
unseen gene 的网络位置
-> predicted KO delta
-> control cell + predicted KO delta
-> virtual KO cell
```

这属于基于系统先验的 zero-shot / few-shot 外推。它适合候选筛选和方向判断，但必须标记置信度。

### 3.3 未见基因预测的置信度边界

相对可信：

- 基因出现在 Reactome/MSigDB/TF-target/PPI/motif 先验里。
- 与训练过的 KO 基因处于相似 pathway 或 network module。
- 当前细胞类型中该基因或相关通路有表达/活性。
- 多模态信号一致，例如 RNA、protein、ATAC 支持同一方向。

需要谨慎：

- 基因没有真实 KO 训练样本。
- 先验网络覆盖很弱。
- 当前细胞类型中表达很低或调控机制未知。
- 双敲组合中两个基因都没见过。
- KO 效应主要来自现有先验没有覆盖的机制。

因此，对外描述时不能说“模型已经学会所有未见基因 KO”。更准确的说法是：

```text
模型支持对未见基因做网络先验驱动的虚拟敲除外推；
结果适合做候选优先级排序；
缺少先验覆盖或相似训练 KO 时，应输出低置信度。
```

## 4. 输入和输出到底是什么

### 4.1 用户输入

用户不需要提前准备通路分数。推荐输入是 h5ad：

```text
adata.X                 RNA matrix 或 ATAC gene activity matrix
adata.var_names         gene symbols
adata.obs["ko_target"]  control / single KO / double KO 标签
adata.obsm["protein"]   可选，ADT/CITE-seq protein matrix
adata.obsm["atac"]      可选，ATAC state / gene activity / LSI
adata.obsm["chromvar"]  可选，chromVAR motif activity
```

普通 10X 单细胞 RNA 数据也支持，但如果没有 KO 标签，就只能做预测应用，不能做真实准确性评估。

### 4.2 关键命令

RNA+ADT：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\hmpcite_gse243244\hmpcite_perturbation_rna_adt_doubleko.h5ad `
  --ko-col ko_target `
  --target-kos Cebpb+Med12 `
  --prior-dir data\priors `
  --out-dir results\hmpcite_multimodal_doubleko_extra_obsm_demo `
  --extra-obsm protein:protein
```

ATAC gene activity + chromVAR motif activity：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\scperturb_atac\liscovitch_k562_gene_activity_chromvar.h5ad `
  --ko-col ko_target `
  --target-kos KDM6A `
  --prior-dir data\priors `
  --out-dir results\scperturb_atac_gene_activity_chromvar_top100_kdm6a `
  --extra-obsm chromvar:tf `
  --max-extra-features-per-obsm 100
```

真正 RNA+ADT+ATAC 且带 perturbation 标签的数据：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad your_rna_adt_atac_perturbation.h5ad `
  --ko-col ko_target `
  --target-kos GENE1+GENE2 `
  --prior-dir data\priors `
  --out-dir results\your_trimodal_doubleko `
  --extra-obsm protein:protein,atac:atac,chromvar:tf `
  --max-extra-features-per-obsm 100
```

无 perturbation 标签的 DOGMA/TEA-seq：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli score `
  --input-h5ad your_unlabeled_dogma_or_teaseq.h5ad `
  --prior-dir data\priors `
  --out-dir results\trimodal_state_scores `
  --extra-obsm protein:protein,atac:atac,chromvar:tf `
  --max-extra-features-per-obsm 100
```

这一步只生成状态分数，不报告 AUC/R2/MAE。

### 4.3 输出文件

每次评估会输出：

- `summary.csv`：总体指标。
- `delta_table.csv`：真实 KO delta 与虚拟 KO delta。
- `metrics.csv`：每个状态特征的分布改进。
- `auc_summary.csv` 和 `auc_points.csv`：AUC 曲线数据。
- `virtual_cells.csv`：控制细胞、虚拟 KO 细胞、真实 KO 细胞的状态矩阵。
- `01_summary_dashboard.png`：总览图。
- `02_true_vs_virtual_heatmap.png`：真实 KO vs 虚拟 KO heatmap。
- `03_cell_state_umap.png`：单细胞状态变化 UMAP。
- `04_auc_strong_response_roc.png`：AUC/ROC 曲线。
- `analysis_mode.md`：说明该结果是评估模式、预测模式还是仅状态打分。

## 5. 已完成的优化

### 5.1 TF-target/motif prior

已修改 `vkx/core.py` 的先验选择逻辑。现在 `tf_target`、`atac_motif_tf`、`motif_tf_target` 与 `ppi_hub` 一样，会把 term 名称中的 TF/gene 也纳入匹配，而不是只看 target gene overlap。

这对 ATAC/chromVAR 很重要，因为 motif 或 TF term 本身就携带调控因子信息。

### 5.2 chromVAR motif activity 输入

新增脚本：

```text
scripts/35_prepare_scperturb_atac_with_chromvar.py
```

生成文件：

```text
data/scperturb_atac/liscovitch_k562_gene_activity_chromvar.h5ad
```

该文件包含：

```text
7108 cells
23127 gene activity features
2174 chromVAR motif activity features in adata.obsm["chromvar"]
```

### 5.3 大型 obsm 特征筛选

新增 CLI 参数：

```text
--max-extra-features-per-obsm 100
```

用途：当 chromVAR、ATAC LSI 或 motif 矩阵特征很多时，可以按方差保留最有信息的前 N 个特征，减少噪声和运行负担。

### 5.4 真三模态接口

软件现在支持：

```text
--extra-obsm protein:protein,atac:atac,chromvar:tf
```

因此拿到真正 RNA+ADT+ATAC 且带 perturbation 标签的数据后，不需要再改模型代码，可以直接运行评估。

### 5.5 DOGMA/TEA-seq 边界

DOGMA-seq 和 TEA-seq 是三模态测量技术，可以同时获得 RNA、蛋白/表位和 ATAC 信息。但如果数据没有 perturbation/KO 标签，它们不能用于证明虚拟 KO 准确性，只能用于：

- 三模态状态分数生成。
- reference model application。
- 展示虚拟 KO 后预测的状态移动。

## 6. 主要实验结果

### 6.1 HMPCITE RNA+ADT double-KO

数据：GSE243244 HMPCITE-seq，包含 RNA、ADT、GDO-derived perturbation labels。

任务：预测并评估 Cebpb+Med12 双敲。

结果：

```text
ROC-AUC: 0.978
direction cosine: 0.976
mean distribution improvement: 0.548
improved features: 87.5%
mean abs delta error: 0.114
```

解释：

- 真实 KO 和虚拟 KO 的方向高度一致。
- 大多数状态特征从 control 移向真实 KO。
- RNA pathway + ADT 状态空间下，双敲效果是目前最清楚、最强的一组验证。

关键图：

```text
results/user_facing_figures/14_hmpcite_multimodal_doubleko_summary.png
results/user_facing_figures/15_rna_adt_atac_extension_summary.png
results/user_facing_figures/16_modality_extension_result_gallery.png
```

### 6.2 HMPCITE double-KO interaction residual

任务：比较简单相加模型与 interaction residual 模型。

结果：

```text
additive:     MAE 0.195, R2 -0.334, ROC-AUC 0.763
interaction:  MAE 0.113, R2  0.507, ROC-AUC 0.768
```

解释：

- 双敲不是两个单敲简单相加。
- interaction residual 明显改善 MAE 和 R2。
- AUC 只轻微改善，说明 interaction 主要改善幅度和模式拟合，不一定显著改善强响应排序。

### 6.3 scPerturb ATAC gene activity

数据：scPerturb ATAC K562，Liscovitch-Brauer/Sanjana 2021 gene activity。

任务：KDM6A KO。

结果：

```text
ROC-AUC: 0.641
direction cosine: 0.513
mean distribution improvement: -0.000
improved features: 46.7%
mean abs delta error: 0.053
```

解释：

- ATAC/gene activity 层可以评估，但预测难度明显高于 RNA+ADT。
- 分布改进接近 0，说明虚拟细胞整体分布没有明显比 control 更接近真实 KO。

### 6.4 scPerturb ATAC + 全量 chromVAR

输入：gene activity 主矩阵 + 2174 个 chromVAR motif activity。

结果：

```text
ROC-AUC: 0.552
direction cosine: 0.432
mean distribution improvement: 0.019
improved features: 59.4%
mean abs delta error: 0.058
```

解释：

- 全量 motif activity 不是自动增强，反而降低 AUC 和方向一致性。
- 原因可能是 motif 特征很多、弱噪声多，而且 KDM6A 不是典型 sequence-specific TF。

### 6.5 scPerturb ATAC + top100 chromVAR

输入：gene activity 主矩阵 + 方差最高 100 个 chromVAR motif。

结果：

```text
ROC-AUC: 0.586
direction cosine: 0.618
mean distribution improvement: 0.043
improved features: 61.5%
mean abs delta error: 0.050
```

解释：

- top100 chromVAR 比全量 chromVAR 更稳。
- 方向一致性和 MAE 优于 gene activity only。
- 但 AUC 仍低于 gene activity only，说明强响应排序还不稳定。

关键图：

```text
results/user_facing_figures/18_atac_chromvar_prior_ablation.png
```

### 6.6 scPerturb ATAC + weighted prior + hybrid chromVAR + locus-aware peaks

用户指出加入 ATAC 后还应该看到 peak-level 峰图，这个问题是正确的。之前的 ATAC 结果主要使用 gene activity 和 chromVAR，它们是 peak-derived summary 或 motif summary，不是真正的 peak-level view。我们随后从 `peak_bc` 中接入细胞×peak matrix，选择 KDM6A 附近 peaks 和高可变 peaks，并作为 `obsm["peak"]` 输入模型。

在最新版本中，这一部分进一步优化为：

- TF-target/motif prior 不再只是 overlap，而是带 library weight、direct TF hit bonus 和 term coverage。
- chromVAR/ATAC 额外模态不再只按方差筛选，而支持 `variance`、`ko_effect`、`hybrid` 三种选择；带 KO 标签的 ATAC 数据推荐 `hybrid`。
- peak-level ATAC 不再只按全局方差或可及细胞数筛选，而使用 KDM6A locus peaks、`markerpeak_target` 和全局稳定 peaks 的组合。

输入文件：

```text
data/scperturb_atac/liscovitch_k562_gene_activity_chromvar_peaks.h5ad
```

运行方式：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\scperturb_atac\liscovitch_k562_gene_activity_chromvar_peaks.h5ad `
  --ko-col ko_target `
  --target-kos KDM6A `
  --prior-dir data\priors `
  --out-dir results\scperturb_atac_weighted_hybrid_peak_kdm6a `
  --extra-obsm chromvar:tf,peak:peak `
  --max-extra-features-per-obsm 100 `
  --extra-feature-selection hybrid
```

结果：

```text
ROC-AUC: 0.658
direction cosine: 0.684
mean abs delta error: 0.057
improved features: 37.8%
mean distribution improvement: -0.159
```

解释：

- 加权 prior + hybrid chromVAR 不加 peaks 时，AUC 为 0.645，方向一致性为 0.728，MAE 为 0.048，整体分布改进为 0.113。
- 加入 locus-aware peaks 后，AUC 进一步提高到 0.658，说明 peak-level 输入更有利于强响应特征排序。
- 但 distribution improvement 为负，说明 peak-level 单细胞分布仍难预测。
- 这证明 peak-level 输入和可视化已经接入，但也提示 ATAC 峰层面的生成需要更强的 peak selection 和局部调控先验。

新增可视化：

```text
results/scperturb_atac_gene_activity_chromvar_peak_kdm6a/05_atac_peak_level_changes.png
results/scperturb_atac_weighted_hybrid_peak_kdm6a/05_atac_peak_level_changes.png
results/user_facing_figures/19_atac_peak_level_visualization.png
results/user_facing_figures/20_atac_weighted_prior_feature_selection_summary.png
```

## 7. 如何判断方法好不好

不能只看一个指标。

### 7.1 AUC

AUC 衡量模型能否把强响应状态特征排在前面。适合看“哪些 pathway/protein/motif 最可能真的响应 KO”。

注意：

- 特征数很少时 AUC 不稳定。
- 加入大量 motif 后，AUC 可能下降，因为排序任务变难。
- AUC 高不代表幅度完全正确。

### 7.2 Direction cosine

衡量虚拟 KO delta 和真实 KO delta 的方向是否一致。越接近 1 越好。

这对虚拟敲除很重要，因为用户通常关心“细胞状态往哪边变”。

### 7.3 MAE

衡量预测变化幅度与真实变化幅度的平均误差。越低越好。

### 7.4 Heatmap

最适合非专业用户阅读：

- 左列：真实 KO 变化。
- 中列：虚拟 KO 变化。
- 右列：预测误差。

颜色方向一致说明方向对；颜色深浅接近说明幅度也接近。

### 7.5 UMAP

UMAP 显示敲除前后单细胞状态是否发生可见移动。

注意：UMAP 是可视化，不是准确性证明。准确性仍要看真实 KO 标签下的 heatmap、AUC、R2/MAE。

## 8. 方法适用性

### 适合

- 小样本 perturb-seq / CRISPR screen。
- RNA-only 或 RNA+ADT perturbation 数据。
- 有单敲和部分双敲标签的数据。
- 需要解释“敲除什么基因，引起哪些通路/蛋白/调控状态变化”的场景。
- 需要在普通 10X 数据上应用 reference model 做虚拟状态转换的场景。

### 谨慎

- ATAC-only 或 gene activity-only 数据。
- chromVAR motif 特征很多但缺少筛选时。
- KO 基因在 reference training 中从未出现时。
- 只有普通无标签 10X/Multiome/DOGMA/TEA-seq 数据时。

### 不适合直接声称准确性

- 没有真实 perturbation 标签的数据。
- 只有预测结果、没有真实 KO 对照的数据。
- 特征很少但只报告 AUC 的情况。

## 9. 与生成模型的关系

当前模型不是自由 diffusion/VAE，而是 hard-constrained residual/PLS baseline。

这不是弱点，而是小样本场景下的稳妥选择：

- 先学 KO 方向和平均 delta。
- 再在这个方向附近学习不确定性。
- 未来 VAE/flow matching/diffusion 只负责生成“方向附近的细胞状态分布”，不能自由偏离 baseline。

推荐下一步：

```text
residual/PLS hard constraint
  + calibrated pathway/protein/ATAC delta
  + uncertainty generator
  + cell-level conditional sampling
```

## 10. 已生成的核心图

```text
results/user_facing_figures/01_method_workflow.png
results/user_facing_figures/08_auc_roc_curves.png
results/user_facing_figures/11_umap_before_after_examples.png
results/user_facing_figures/14_hmpcite_multimodal_doubleko_summary.png
results/user_facing_figures/15_rna_adt_atac_extension_summary.png
results/user_facing_figures/16_modality_extension_result_gallery.png
results/user_facing_figures/17_multimodal_input_visualization_matrix.png
results/user_facing_figures/18_atac_chromvar_prior_ablation.png
results/user_facing_figures/19_atac_peak_level_visualization.png
```

## 11. 数据来源

- HMPCITE-seq GSE243244：GEO 页面说明该数据包含 Perturbation cDNA、ADT 和 GDO 文件，并用于研究 Cebpb、Med12 等组合扰动。
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243244
- scPerturb ATAC files：Figshare 页面说明该资源包含 ChromVar、gene_scores 等 ATAC feature matrices，并来自单细胞扰动数据库。
  https://plus.figshare.com/articles/dataset/scPerturb_Single-Cell_Perturbation_Data_ATAC_files/24160968
- DOGMA-seq：文献介绍可同时测量 chromatin accessibility、gene expression 和 protein。
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8763625/
- TEA-seq：文献介绍可同时测量 transcriptomics、epitopes 和 chromatin accessibility。
  https://elifesciences.org/articles/63632

## 12. 最终判断

这个方法的核心思路是成立的，但适用边界要写清楚。

最有说服力的结果是 RNA+ADT double-KO，因为它有真实标签、有多模态状态、有清晰 heatmap/UMAP/ROC。ATAC/chromVAR 是必要扩展，但目前还属于困难场景，需要 motif/TF-target 筛选、加权和更好的调控先验。DOGMA/TEA-seq 这类无标签三模态数据可以展示“模型能吃进去多模态并输出虚拟状态变化”，但不能当作准确性验证。

因此，当前版本适合定位为：

```text
一种面向小样本单细胞扰动数据的、网络先验约束的、多模态状态空间虚拟敲除方法。
```
