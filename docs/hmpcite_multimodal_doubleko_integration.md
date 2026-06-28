# HMPCITE-seq 多模态 double-KO 数据接入记录

## 数据来源

接入数据：GSE243244 HMPCITE-seq perturbation sample。

下载地址：

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243244
```

原始补充包：

```text
data/hmpcite_gse243244/GSE243244_RAW.tar
```

本次使用的矩阵：

- `GSM7781666_Pertubation_cDNA.tar.gz`: RNA/cDNA matrix
- `GSM7781668_Pertubation_ADT.tar.gz`: ADT/protein matrix
- `GSM7781669_Pertubation_GDO.tar.gz`: guide-derived oligo matrix

## 为什么这个数据合适

这个数据同时具备：

- RNA 表达矩阵
- ADT 蛋白矩阵
- GDO guide 标签
- 大量单基因和双基因 perturbation 组合

因此它比普通 Perturb-seq 更适合测试“多模态 + 双基因虚拟敲除”。

## 处理方式

脚本：

```text
scripts/30_prepare_hmpcite_double_ko_multimodal.py
```

处理规则：

1. 从 GDO matrix 中解析 guide 对应基因。
2. GDO count 阈值设为 10。
3. 保留 control、single KO、double KO 细胞。
4. 去掉三基因及以上组合，避免复杂多扰动影响当前双敲评估。
5. RNA 做 log-normalization。
6. ADT 蛋白放入 `obsm["protein"]`。
7. KO 标签写入 `obs["ko_target"]`。

生成文件：

```text
data/hmpcite_gse243244/hmpcite_perturbation_rna_adt_doubleko.h5ad
data/hmpcite_gse243244/hmpcite_ko_counts_threshold10.csv
```

数据规模：

```text
19,714 cells
32,286 RNA genes
2 ADT proteins
11 single-gene perturbations
55 real double-gene perturbation combinations
```

## 单个真实多模态双敲评估

测试目标：

```text
Cebpb+Med12
```

运行命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\hmpcite_gse243244\hmpcite_perturbation_rna_adt_doubleko.h5ad `
  --ko-col ko_target `
  --target-kos Cebpb+Med12 `
  --prior-dir data\priors `
  --out-dir results\hmpcite_multimodal_doubleko_cebp_med12 `
  --dataset-name "HMPCITE-seq GSE243244" `
  --modality "RNA + ADT + GDO-derived double KO" `
  --representation "auto-derived pathway/protein scores" `
  --protein-obsm protein `
  --max-pathways 30 `
  --calibrate auto
```

结果：

```text
mean distribution improvement: 0.548
improved features: 87.5%
direction cosine: 0.976
mean absolute delta error: 0.114
ROC-AUC: 0.978
```

解释：

这个结果说明，在 `Cebpb+Med12` 真实双敲上，虚拟 KO 对主要 pathway/protein 状态变化的方向预测较好，强响应特征排序也较好。

## 55 个真实多模态双敲组合的 interaction residual 评估

先导出所有 KO 的状态变化：

```powershell
.\.venv\Scripts\python.exe scripts\31_export_hmpcite_state_delta.py
```

再运行 double-interaction：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli double-interaction `
  --delta-csv results\hmpcite_multimodal_doubleko_state_delta.csv `
  --ko-col ko_genes `
  --n-ko-col n_ko_genes `
  --target-prefix delta_ `
  --prior-dir data\priors `
  --out-dir results\hmpcite_multimodal_doubleko_interaction
```

总体结果：

```text
single_gene_additive: MAE 0.195, R2 -0.334, ROC-AUC 0.763
interaction_residual: MAE 0.113, R2  0.507, ROC-AUC 0.768
```

解释：

- interaction residual 明显降低 MAE。
- interaction residual 明显提升 R2。
- ROC-AUC 变化不大，说明这个改进主要来自“变化幅度预测更准”，而不是强响应排序显著改变。

## 用户展示图

用户版总结图：

```text
results/user_facing_figures/14_hmpcite_multimodal_doubleko_summary.png
```

关键结果图：

```text
results/hmpcite_multimodal_doubleko_cebp_med12/02_true_vs_virtual_heatmap.png
results/hmpcite_multimodal_doubleko_cebp_med12/03_cell_state_umap.png
results/hmpcite_multimodal_doubleko_cebp_med12/04_auc_strong_response_roc.png
results/hmpcite_multimodal_doubleko_interaction/double_interaction_metrics.png
```

## 当前结论

现在已经真正接入了一个可评估的多模态 double-KO 数据集。

这一步解决了之前的主要短板：

- 不再只是 RNA-only 双敲。
- 不再只是多模态单敲。
- 不再只是普通 10X prediction-only。

现在可以在同一个数据中同时评估：

```text
RNA pathway state
ADT protein state
真实 double KO 标签
虚拟 double KO 预测
interaction residual 是否优于 additive baseline
```

## 注意事项

GDO 阈值会影响 KO 标签分配。当前阈值为 10，是为了在标签可靠性和细胞数量之间折中。后续可以做阈值敏感性分析，例如 threshold 5、10、20 下结果是否稳定。
