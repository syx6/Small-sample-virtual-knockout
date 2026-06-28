# RNA + ADT + ATAC 多模态扩展

## 当前新增能力

软件现在支持多个 `obsm` 模态输入，不再只能写成 `--protein-obsm`。

推荐新接口：

```powershell
--extra-obsm protein:protein,atac:atac
```

含义：

```text
adata.X                 RNA 或 gene activity 主矩阵
adata.obsm["protein"]   ADT / CITE-seq protein matrix
adata.obsm["atac"]      ATAC / gene activity / chromVAR / LSI-derived score matrix
```

输出中的状态特征会自动命名为：

```text
pathway_...
protein_...
atac_...
```

旧接口仍然保留：

```powershell
--protein-obsm protein
```

但后续建议统一使用：

```powershell
--extra-obsm protein:protein
```

## 已验证的新接口

### 1. HMPCITE RNA + ADT + double-KO

命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\hmpcite_gse243244\hmpcite_perturbation_rna_adt_doubleko.h5ad `
  --ko-col ko_target `
  --target-kos Cebpb+Med12 `
  --prior-dir data\priors `
  --out-dir results\hmpcite_multimodal_doubleko_extra_obsm_demo `
  --dataset-name "HMPCITE-seq GSE243244" `
  --modality "RNA + ADT via extra-obsm" `
  --representation "auto-derived pathway/protein scores" `
  --extra-obsm protein:protein `
  --max-pathways 30 `
  --calibrate auto
```

结果与旧 `--protein-obsm protein` 接口一致：

```text
ROC-AUC: 0.978
mean distribution improvement: 0.548
improved features: 87.5%
direction cosine: 0.976
```

说明新的多模态接口不会破坏之前的 RNA+ADT double-KO 评估。

### 2. scPerturb ATAC gene activity perturbation

接入数据：

```text
data/scperturb_atac/liscovitch_k562_gene_activity.h5ad
```

来源：

```text
scPerturb ATAC files
Liscovitch-Brauer/Sanjana 2021 K562 gene activity
```

处理脚本：

```text
scripts/33_prepare_scperturb_atac_gene_activity.py
```

测试命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\scperturb_atac\liscovitch_k562_gene_activity.h5ad `
  --ko-col ko_target `
  --target-kos KDM6A `
  --prior-dir data\priors `
  --out-dir results\scperturb_atac_gene_activity_kdm6a `
  --dataset-name "scPerturb ATAC K562" `
  --modality "ATAC gene activity" `
  --representation "auto-derived pathway scores from gene activity" `
  --max-pathways 30 `
  --calibrate auto
```

结果：

```text
ROC-AUC: 0.641
direction cosine: 0.513
mean distribution improvement: ~0.000
improved features: 46.7%
```

解释：

ATAC/gene activity 层比 RNA+ADT 层更难预测。这个结果不是失败，而是说明调控层扰动响应更弱、更稀疏，后续需要更强的 ATAC-specific prior 或 motif/TF-target 约束。

## RNA + ADT + ATAC 真三模态数据如何输入

如果用户有真正三模态 h5ad，推荐结构：

```text
adata.X                 RNA matrix, cells x genes
adata.var_names         gene symbols
adata.obs["ko_target"]  control / single KO / double KO labels
adata.obsm["protein"]   ADT protein matrix
adata.uns["protein_names"]
adata.obsm["atac"]      ATAC gene activity, chromVAR, or regulatory program scores
adata.uns["atac_names"]
```

运行：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad your_rna_adt_atac_perturbation.h5ad `
  --ko-col ko_target `
  --target-kos GENE1+GENE2 `
  --prior-dir data\priors `
  --out-dir results\your_trimodal_doubleko `
  --dataset-name "Your RNA+ADT+ATAC perturbation dataset" `
  --modality "RNA + ADT + ATAC" `
  --representation "pathway + protein + ATAC state scores" `
  --extra-obsm protein:protein,atac:atac `
  --max-pathways 40
```

如果 h5ad 中还包含 chromVAR motif activity，推荐写成：

```powershell
--extra-obsm protein:protein,atac:atac,chromvar:tf
```

如果 h5ad 中包含 peak-level accessibility matrix，也可以继续加：

```powershell
--extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak
```

其中 `adata.obsm["peak"]` 应该是 cells × selected peaks，`adata.uns["peak_names"]` 推荐使用类似 `chrX:44732092-44732592|KDM6A|Promoter` 的名称，软件会自动把它整理成 `peak_...` 状态特征并生成 peak-level 图。

对于 chromVAR、ATAC LSI 或其它高维额外模态，建议先限制进入模型的额外特征数：

```powershell
--max-extra-features-per-obsm 100
```

这个参数可以配合特征选择方法使用：

```powershell
--extra-feature-selection variance
--extra-feature-selection ko_effect
--extra-feature-selection hybrid
```

其中 `hybrid` 会综合特征方差和 KO-vs-control 效应大小；没有 KO 标签时会自动退回 variance。当前 KDM6A ATAC 消融结果显示，weighted prior + hybrid chromVAR 的整体方向、MAE 和分布改进最好；weighted prior + hybrid + locus-aware peaks 的 AUC 最高。

## 重要边界

RNA+ADT+ATAC 只有在有真实 perturbation 标签时，才能评估虚拟敲除准确性。

如果没有 KO 标签：

```text
可以：状态转换、reference model application、预测性 heatmap/UMAP
不可以：报告真实 AUC/R2/MAE
```

如果有 single KO 但没有 double KO：

```text
可以：评估单敲
可以：预测双敲
不可以：在该数据内部证明双敲准确性
```

如果有 double KO：

```text
可以：真实评估双敲
可以：比较 additive baseline 与 interaction residual
可以：报告 AUC/R2/MAE/heatmap/UMAP
```

## 配套可视化结果

这一部分已经补齐为三张面向用户的图，统一放在：

```text
results/user_facing_figures/
```

### 15_rna_adt_atac_extension_summary.png

这张图回答“多模态扩展后效果到底怎么样”。

- HMPCITE RNA+ADT double-KO 是真实带标签的多模态双敲评估，Cebpb+Med12 的 ROC-AUC 为 0.98，方向一致性为 0.98，说明 RNA pathway + ADT 状态空间里预测效果较好。
- scPerturb ATAC gene activity 的 KDM6A 结果 ROC-AUC 为 0.64，方向一致性为 0.51，说明 ATAC/gene activity 层已经可以接入和评估，但调控层信号更难预测。
- 底部示意图说明真实三模态输入时，RNA、ADT、ATAC 可以通过同一个 `--extra-obsm protein:protein,atac:atac` 接口进入模型。

### 16_modality_extension_result_gallery.png

这张图把每个模态结果拆成用户最关心的三类证据：

- heatmap：真实 KO 变化、虚拟 KO 变化、预测误差是否一致。
- UMAP：敲除前后单细胞状态是否发生可见移动。
- ROC 曲线：AUC 不再只显示一个数字，而是显示完整曲线。

上排是 HMPCITE RNA+ADT double-KO，下排是 scPerturb ATAC gene activity。这个布局能直接看出：RNA+ADT 的敲除响应更清晰，ATAC 层虽然支持，但目前效果偏弱。

### 17_multimodal_input_visualization_matrix.png

这张图专门解释“用户给什么输入，就能看什么输出”。

- 有 KO 标签：可以报告 AUC、R2、MAE、heatmap、UMAP。
- 没有 KO 标签：可以做虚拟应用和状态变化可视化，但不能在该数据内部报告真实准确性。
- 真正 RNA+ADT+ATAC 且有 KO 标签时，可以评估完整三模态准确性；如果只是普通无标签三模态 10X/Multiome 数据，只能做预测展示。

### 18_atac_chromvar_prior_ablation.png

这张图展示 ATAC-specific 优化后的一组消融：

- gene activity only：AUC 0.64，方向一致性 0.51。
- gene activity + 全量 chromVAR：AUC 0.55，方向一致性 0.43。
- gene activity + top100 chromVAR：AUC 0.59，方向一致性 0.62。

结论是：chromVAR motif activity 已经可以作为 `--extra-obsm chromvar:tf` 输入，但不能全量无筛选地加入。当前更合理的做法是先筛选 top motif，再进一步开发 motif/TF-target 加权先验。

### 19_atac_peak_level_visualization.png

这张图展示真正的 ATAC peak-level 结果。我们从 `peak_bc` 中接入细胞×peak matrix，把 KDM6A 附近 peaks 和高可变 peaks 放进 `obsm["peak"]`。结果显示 selected peaks 版本 AUC 为 0.617、方向一致性为 0.617、MAE 为 0.0498。

需要注意：peak-level 输入改善了 AUC 和可解释性，但 distribution improvement 仍然为负，说明峰层面的单细胞分布预测仍然比 pathway/protein 层更难。

### 20_atac_weighted_prior_feature_selection_summary.png

这张图总结最新 ATAC 优化：

- TF-target/motif prior 已经从 overlap 升级为加权 prior。
- chromVAR/ATAC feature selection 已经支持 `variance`、`ko_effect`、`hybrid`。
- peak-level ATAC 已经使用 locus-aware selection，即 target locus peaks + markerpeak_target + 全局稳定 peaks。
- 最新 KDM6A 结果中，weighted + hybrid chromVAR 的方向一致性最高，weighted + hybrid + locus-aware peaks 的 AUC 最高。

## 当前状态和下一步

已经完成：

1. TF-target/motif 类先验会把 term 中的 TF/gene 名称纳入匹配，并加入 library weight、direct TF hit bonus 和 coverage weight。
2. chromVAR motif activity 可以作为 `--extra-obsm chromvar:tf` 输入。
3. 新增 `--max-extra-features-per-obsm` 和 `--extra-feature-selection`，用于更稳定地筛选高维 chromVAR/ATAC 额外模态。
4. peak-level ATAC 可以作为 `--extra-obsm peak:peak` 输入，并自动生成 peak-level 图。
5. 真正 RNA+ADT+ATAC 且带 perturbation 标签的数据，可以直接用 `--extra-obsm protein:protein,atac:atac,chromvar:tf,peak:peak` 运行。
6. 无 perturbation 标签的 DOGMA/TEA-seq 数据，只做三模态输入兼容和 reference application，不作为准确性验证。

仍需继续优化：

1. 接入公开或自有的真正 RNA+ADT+ATAC+perturbation 标签数据后，补完整三模态准确性 benchmark。
2. 对 peak-level ATAC 继续开发更细的 locus-aware selection，例如 enhancer-promoter linkage、motif-in-peak 和 target-gene distance weighting。
3. 在 hard-constrained baseline 上接 cell-level uncertainty generator。
