# 单基因敲除和双基因敲除结果整理

## 1. 总体结论

当前方法已经整理成一个可复用的虚拟敲除流程：

```text
用户输入原始单细胞/多组学矩阵
-> 软件自动转换为可解释状态表示
-> 用系统先验 + residual/PLS 预测 KO 方向
-> 生成虚拟 KO 单细胞状态
-> 输出 heatmap、UMAP、AUC/ROC 和指标表
```

单基因敲除目前表现更稳。以 Papalexi 的 `STAT1` 为例，虚拟 KO 方向和真实 KO 方向高度一致，AUC 也比较好。

双基因敲除可以支持，但难度更高。以 Norman 的 `CEBPB+CEBPA` 为例，效果很好；但在所有双基因组合上，MAPK/TGFB 和 Pro-growth 这类非线性组合效应仍然是主要难点。

需要特别说明：

- 当前真正跑通的多模态示例主要是 Papalexi 的 RNA + ADT protein。
- Norman、Datlinger、Dixit 结果主要是 RNA 或 RNA-derived program/pathway。
- 普通 10X 无 KO 标签数据目前支持状态转换，但不能在数据内部评估真实准确性。
- 普通 10X 完整虚拟 KO 应用需要外部 reference perturbation model，这一步还需要继续做成正式接口。
- 双敲非线性效应还没有完全解决，尤其是 MAPK/TGFB 这类 program。

更完整的范围说明见：

```text
docs/current_scope_and_limitations.md
```

## 2. 方法说明

### 2.1 用户真正输入什么

普通用户不需要输入 pathway score。

推荐输入是：

| 输入 | 说明 |
|---|---|
| `.h5ad` | 单细胞 RNA 矩阵，`adata.X = cells x genes` |
| `adata.var_names` | gene symbols |
| `adata.obs[ko_col]` | 每个细胞的 KO / perturbation 标签 |
| `adata.obsm["protein"]` | 可选，CITE-seq / ADT protein |
| `adata.obsm["atac"]` 或 gene activity | 可选，ATAC / regulatory score |

如果是 CSV，也可以是：

```text
ko_target,STAT1,JAK2,IRF1,...
control,0.1,1.3,0.4,...
STAT1,0.0,1.2,0.2,...
```

### 2.2 软件内部自动生成什么

软件会自动把原始 RNA 矩阵转成 pathway/program score。

```text
RNA gene expression
-> Reactome / MSigDB / TF-target gene sets
-> pathway/program score
```

如果有 protein/ADT/ATAC，也会合并成状态表示：

```text
state = RNA pathway/program score + protein score + ATAC/regulatory score
```

这些 state score 是模型内部表示，不是要求用户自己准备。

### 2.3 模型如何做虚拟敲除

核心公式：

```text
virtual KO cell = control cell + predicted KO delta
```

其中 `predicted KO delta` 来自：

1. 训练 KO 中观测到的真实状态变化。
2. KO gene 对应的系统先验特征。
3. Reactome / MSigDB / TF-target / PPI 等网络先验。
4. residual / PLS 模型预测出的 KO 方向。

这个方法不是让深度生成模型自由生成细胞，而是把 KO 方向强约束住。这样更适合小样本。

### 2.4 单敲和双敲分别怎么输入

单基因敲除：

```powershell
--target-kos STAT1
```

双基因敲除：

```powershell
--target-kos CEBPB+CEBPA
```

批量评估多个 KO 时才写成逗号分隔：

```powershell
--target-kos STAT1,JAK2,IFNGR2,IRF1
```

## 3. 单基因敲除结果

### 3.1 示例：Papalexi `STAT1` 单敲

运行命令：

```powershell
.\.venv\Scripts\python.exe -m vkx.cli run `
  --input-h5ad data\papalexi_small_pathway.h5ad `
  --ko-col ko_target `
  --target-kos STAT1 `
  --prior-dir data\priors `
  --out-dir results\software_interface_single_gene_demo `
  --dataset-name "Papalexi ECCITE-seq" `
  --modality "raw RNA matrix + ADT protein obsm" `
  --representation "auto-derived pathway/protein scores" `
  --protein-obsm protein `
  --calibrate none `
  --max-pathways 24
```

输入：

| 项目 | 内容 |
|---|---|
| 数据集 | Papalexi ECCITE-seq |
| 原始输入 | `.h5ad` RNA 矩阵 |
| KO 标签 | `adata.obs["ko_target"]` |
| 额外模态 | `adata.obsm["protein"]` |
| 目标 KO | `STAT1` |
| 内部状态 | 自动生成的 pathway/protein score |

输出目录：

```text
results/software_interface_single_gene_demo
```

主要指标：

| 指标 | 值 | 含义 |
|---|---:|---|
| n_ko | 1 | 评估 1 个单基因 KO |
| n_features | 28 | 24 个 RNA pathway + 4 个 protein |
| mean distribution improvement | 0.144 | 虚拟 KO 比 control 更接近真实 KO |
| improved fraction | 0.821 | 82.1% 的状态特征得到改善 |
| mean direction cosine | 0.879 | 虚拟 KO 方向和真实 KO 方向高度一致 |
| mean abs delta error | 0.424 | KO 变化幅度仍有误差 |
| strong-response ROC-AUC | 0.878 | 能较好识别强响应状态特征 |

解释：

- `STAT1` 单敲整体效果较好。
- 方向一致性高，说明模型知道 KO 后细胞状态应该往哪里变。
- AUC 高，说明强响应 pathway/protein 能被排在前面。
- 主要不足是变化幅度偏保守，真实 KO 的强变化比虚拟 KO 更大。

### 3.2 单敲图怎么看

主要图：

| 图 | 路径 | 作用 |
|---|---|---|
| Summary dashboard | `results/software_interface_single_gene_demo/01_summary_dashboard.png` | 一眼看总体效果 |
| True vs virtual heatmap | `results/software_interface_single_gene_demo/02_true_vs_virtual_heatmap.png` | 看真实 KO 和虚拟 KO 的 pathway/protein 变化 |
| UMAP | `results/software_interface_single_gene_demo/03_cell_state_umap.png` | 看虚拟 KO 细胞是否向真实 KO 状态移动 |
| ROC curve | `results/software_interface_single_gene_demo/04_auc_strong_response_roc.png` | 看强响应识别能力 |

heatmap 的读法：

```text
左：真实 STAT1 KO 相对 control 的变化
中：虚拟 STAT1 KO 相对 control 的变化
右：虚拟 - 真实 的误差
```

如果左图和中图颜色方向一致，说明方向预测对了。
如果右图颜色浅，说明幅度也接近。

### 3.3 多个单基因 KO 的评估

除了 `STAT1` 单独示例，我们还评估了 Papalexi 中多个单基因 KO：

| KO | direction cosine | mean abs delta error | UMAP centroid improvement |
|---|---:|---:|---:|
| STAT1 | 0.931 | 0.196 | 0.297 |
| JAK2 | 0.600 | 0.185 | 0.076 |
| IFNGR2 | 0.966 | 0.102 | -0.396 |
| IRF1 | 0.189 | 0.101 | 0.676 |

解释：

- `STAT1` 和 `IFNGR2` 的方向一致性很好。
- `JAK2` 中等。
- `IRF1` 在 heatmap 方向一致性低，但 UMAP 质心移动较好，说明不同指标看到的是不同层面。
- `IFNGR2` 方向一致但 UMAP 变差，说明方向对不代表单细胞云团一定移动到真实 KO 区域。

所以单基因 KO 不能只看一个指标，至少要同时看：

1. direction cosine
2. mean abs delta error
3. UMAP centroid improvement
4. heatmap
5. AUC/ROC

### 3.4 其他 RNA-only 单敲数据

在 Datlinger 和 Dixit RNA-only 数据上，方向通常还可以，但幅度需要校准。

代表结果：

| 数据集 | KO | direction cosine | mean abs delta error |
|---|---|---:|---:|
| Datlinger | LAT | 0.877 | 0.052 |
| Datlinger | JUND | 0.863 | 0.044 |
| Datlinger | FOS | 0.944 | 0.048 |
| Datlinger | LCK | 0.290 | 0.212 |
| Dixit | ELF1 | 0.994 | 0.015 |
| Dixit | ELK1 | 0.969 | 0.016 |
| Dixit | GABPA | 0.936 | 0.024 |
| Dixit | CREB1 | -0.413 | 0.068 |

解释：

- RNA-only 并不是不能做，很多 KO 的方向是对的。
- 失败例子包括 `CREB1` 和 `LCK`。
- RNA-only 的主要问题是变化幅度和分布形状，不是所有 KO 都能靠 pathway score 解决。
- 因此 RNA-only 默认建议 `--calibrate auto`。

## 4. 双基因敲除结果

### 4.1 示例：Norman `CEBPB+CEBPA` 双敲

运行命令：

```powershell
.\.venv\Scripts\python.exe scripts\27_export_norman_program_state_csv.py

.\.venv\Scripts\python.exe -m vkx.cli fit `
  --state-csv data\examples\norman_gene_program_state.csv `
  --ko-col ko_target `
  --target-kos CEBPB+CEBPA `
  --prior-dir data\priors `
  --out-dir results\software_interface_double_gene_demo `
  --dataset-name "Norman Perturb-seq" `
  --modality "single-cell RNA perturb-seq" `
  --representation "gene program scores" `
  --calibrate none
```

这里用的是 `fit` 而不是 `run`，因为当前 Norman 示例已经被整理成 gene program score。普通用户的原始 `.h5ad` 仍然推荐用 `run`。

输入：

| 项目 | 内容 |
|---|---|
| 数据集 | Norman Perturb-seq |
| 状态表示 | gene program scores |
| 目标 KO | `CEBPB+CEBPA` |
| 特征数 | 5 个 gene program |

输出目录：

```text
results/software_interface_double_gene_demo
```

主要指标：

| 指标 | 值 | 含义 |
|---|---:|---|
| n_ko | 1 | 评估 1 个双基因 KO |
| n_features | 5 | 5 个 gene program |
| mean distribution improvement | 0.241 | 虚拟 KO 明显比 control 更接近真实 KO |
| improved fraction | 0.600 | 60% 的 program 得到改善 |
| mean direction cosine | 0.960 | 双敲方向预测非常好 |
| mean abs delta error | 0.127 | 幅度误差中等 |
| strong-response ROC-AUC | 1.000 | 在 5 个 program 中强响应排序完全正确 |

注意：

`AUC=1.000` 不能过度解读，因为这里只有 5 个 program 特征。它说明在这个双敲例子中，模型把强响应 program 排对了，但不是大规模分类任务。

### 4.2 `CEBPB+CEBPA` 具体 program 变化

| Program | 真实 KO delta | 虚拟 KO delta | 解释 |
|---|---:|---:|---|
| ERYTHROID | -0.300 | -0.468 | 方向正确，幅度略大 |
| GRANULOCYTE_APOPTOSIS | 0.503 | 0.720 | 方向正确，幅度略大 |
| MAPK_TGFB | 0.248 | 0.090 | 方向正确，但低估 |
| PRO_GROWTH | -0.211 | -0.295 | 方向正确 |
| PIONEER_TF | 0.235 | 0.229 | 很接近 |

解释：

- 5 个 program 的方向全部基本正确。
- `GRANULOCYTE_APOPTOSIS` 是真实最强正响应，模型也预测为强正响应。
- `MAPK_TGFB` 被明显低估，是这个组合里主要误差之一。

### 4.3 双敲图怎么看

主要图：

| 图 | 路径 | 作用 |
|---|---|---|
| Summary dashboard | `results/software_interface_double_gene_demo/01_summary_dashboard.png` | 双敲总体效果 |
| True vs virtual heatmap | `results/software_interface_double_gene_demo/02_true_vs_virtual_heatmap.png` | 5 个 program 的真实/虚拟变化 |
| UMAP | `results/software_interface_double_gene_demo/03_cell_state_umap.png` | 虚拟双敲细胞是否靠近真实双敲细胞 |
| ROC curve | `results/software_interface_double_gene_demo/04_auc_strong_response_roc.png` | 强响应 program 排序 |

## 5. 52 个 Norman 双基因组合的整体评估

为了看双敲是否只是在一个例子上好，我们还评估了 Norman 中 52 个双基因组合。

系统先验模型 `system_prior_ridge` 在全部 52 个组合上的平均表现：

| 指标 | 值 |
|---|---:|
| n_combos | 52 |
| mean MAE | 0.131 |
| mean R2 | 0.274 |
| mean ROC-AUC | 0.682 |
| mean PR-AUC | 0.704 |

分 program 看：

| Program | MAE | R2 | ROC-AUC | PR-AUC |
|---|---:|---:|---:|---:|
| ERYTHROID | 0.239 | 0.391 | 0.631 | 0.872 |
| GRANULOCYTE_APOPTOSIS | 0.132 | 0.745 | 0.675 | 0.770 |
| MAPK_TGFB | 0.084 | -0.139 | 0.536 | 0.325 |
| PRO_GROWTH | 0.128 | -0.077 | 0.880 | 0.920 |
| PIONEER_TF | 0.074 | 0.451 | 0.688 | 0.632 |

解释：

- `GRANULOCYTE_APOPTOSIS` 和 `PIONEER_TF` 的连续变化预测较好。
- `PRO_GROWTH` 虽然 R2 不好，但强响应识别 AUC/PR-AUC 很好，说明排序能力比幅度预测更好。
- `MAPK_TGFB` 是主要难点，R2 和 AUC 都不稳定。

### 5.1 seen combo vs unseen combo

双敲组合分两类：

| 类型 | 含义 | n |
|---|---|---:|
| all_genes_seen | 双敲里的两个基因都有单基因训练信息 | 12 |
| has_unseen_gene | 至少一个基因缺少单基因训练信息 | 40 |

`all_genes_seen` 平均表现：

| 指标 | 值 |
|---|---:|
| mean MAE | 0.092 |
| mean R2 | 0.382 |
| mean ROC-AUC | 0.945 |
| mean PR-AUC | 0.897 |

`has_unseen_gene` 平均表现：

| 指标 | 值 |
|---|---:|
| mean MAE | 0.143 |
| mean R2 | 0.156 |
| mean ROC-AUC | 0.617 |
| mean PR-AUC | 0.621 |

解释：

- 双敲里两个基因都在训练集中见过时，效果明显更好。
- 有 unseen gene 时，模型仍可用，但泛化难度更大。
- 这说明多基因 KO 的主要挑战是组合泛化，而不是单纯把两个单基因相加。

### 5.2 系统先验对双敲有没有帮助

在 `has_unseen_gene` 的 40 个组合上，system prior 相比 single-gene additive 的 R2 改变量：

| Program | delta R2 |
|---|---:|
| GRANULOCYTE_APOPTOSIS | +0.325 |
| MAPK_TGFB | +0.398 |
| PRO_GROWTH | +0.843 |
| ERYTHROID | -0.190 |
| PIONEER_TF | -0.091 |

解释：

- 系统先验对 `GRANULOCYTE_APOPTOSIS`、`MAPK_TGFB`、`PRO_GROWTH` 有明显帮助。
- 对 `ERYTHROID` 和 `PIONEER_TF` 有时会变差。
- 所以系统先验不是万能增强器，而是对部分程序有帮助。

## 6. 双敲优化：interaction residual 模型

针对双敲非线性问题，我们新增了一个交互修正模型：

```text
single-gene additive prediction
+ prior-based gene-gene interaction residual
= optimized double-KO prediction
```

直观理解：

- 先用单基因效应相加得到基础预测。
- 再根据两个基因是否共同命中 Reactome/MSigDB/TF-target/PPI 先验，学习一个“组合交互修正”。
- 评估方式是 leave-one-combo-out：每次拿出一个双敲组合做测试，其余 51 个组合训练交互修正。

结果文件：

```text
results/norman_double_interaction_metrics.csv
results/norman_double_interaction_predictions.csv
results/figures/norman_double_interaction_model_comparison.png
docs/norman_double_interaction_model.md
scripts/28_norman_double_interaction_model.py
```

整体效果：

| Model | subset | mean MAE | mean R2 | mean ROC-AUC | mean PR-AUC |
|---|---|---:|---:|---:|---:|
| single_gene_additive | all_combos | 0.150 | 0.008 | 0.707 | 0.682 |
| interaction_residual | all_combos | 0.076 | 0.617 | 0.894 | 0.845 |
| single_gene_additive | has_unseen_gene | 0.159 | -0.101 | 0.638 | 0.586 |
| interaction_residual | has_unseen_gene | 0.067 | 0.697 | 0.895 | 0.862 |

这说明显式建模双基因交互以后，双敲效果明显提升。

分 program 看 `all_combos`：

| Program | additive R2 | interaction R2 | additive AUC | interaction AUC |
|---|---:|---:|---:|---:|
| ERYTHROID | 0.508 | 0.841 | 0.592 | 0.858 |
| GRANULOCYTE_APOPTOSIS | 0.549 | 0.909 | 0.815 | 0.976 |
| MAPK_TGFB | -0.636 | 0.042 | 0.617 | 0.836 |
| PRO_GROWTH | -0.888 | 0.536 | 0.756 | 0.898 |
| PIONEER_TF | 0.504 | 0.756 | 0.756 | 0.902 |

重点：

- `MAPK_TGFB` 不再是严重负 R2，虽然仍不是最强 program，但已经明显改善。
- `PRO_GROWTH` 从明显失败变成可用。
- `has_unseen_gene` 子集改善最大，说明 interaction residual 对泛化组合特别有帮助。

当前推荐：

```text
单敲：prior-constrained residual/PLS
双敲：如果有双敲训练数据，优先 additive + interaction residual
双敲：如果没有双敲训练数据，退回 residual/PLS 或 system prior baseline
```

## 7. 单敲 vs 双敲：适用性比较

| 维度 | 单敲 | 双敲 |
|---|---|---|
| 默认用户输入 | `--target-kos STAT1` | `--target-kos CEBPB+CEBPA` |
| 当前稳定性 | 更稳 | 交互模型后明显改善 |
| 最好示例 | Papalexi STAT1 | Norman CEBPB+CEBPA |
| 方向预测 | 通常较好 | 好坏取决于组合 |
| 幅度预测 | 仍偏保守 | interaction residual 后明显改善 |
| UMAP 状态移动 | 多数可见 | 好组合很明显，难组合较弱 |
| 推荐用途 | 单基因机制筛选、强响应排序 | 双基因组合初筛、组合机制解释 |

## 8. 当前推荐使用方式

### 7.1 单基因 KO

普通用户首选：

```powershell
python -m vkx.cli run --input-h5ad your_data.h5ad --ko-col ko_target --target-kos STAT1 ...
```

如果是 RNA-only：

```powershell
--calibrate auto
```

如果是 RNA + protein / multiome：

```powershell
--calibrate none
```

### 7.2 双基因 KO

输入：

```powershell
--target-kos GENE1+GENE2
```

例如：

```powershell
--target-kos CEBPB+CEBPA
```

双敲结果必须看 heatmap，不能只看 AUC。因为双敲常见问题不是“能不能识别强响应”，而是组合效应幅度和方向是否准确。

## 9. 当前方法的不足

1. 双敲非线性效应已有 interaction residual 改善，但还需要进一步集成到主 CLI 和 reference model。
2. `MAPK_TGFB` 已明显改善，但仍是相对较难的 program。
3. 普通 10X 没有 KO 标签时，只能做状态转换和虚拟应用，不能在该数据内部评估真实准确性。
4. 当前生成模型不是自由 diffusion/VAE，而是 hard-constrained residual/PLS，更适合小样本，但复杂分布形状模拟能力有限。
5. AUC 在特征数很少时需要谨慎解释。

## 10. 关键文件索引

### 单敲结果

```text
results/software_interface_single_gene_demo/summary.csv
results/software_interface_single_gene_demo/delta_table.csv
results/software_interface_single_gene_demo/auc_summary.csv
results/software_interface_single_gene_demo/01_summary_dashboard.png
results/software_interface_single_gene_demo/02_true_vs_virtual_heatmap.png
results/software_interface_single_gene_demo/03_cell_state_umap.png
results/software_interface_single_gene_demo/04_auc_strong_response_roc.png
```

### 双敲结果

```text
results/software_interface_double_gene_demo/summary.csv
results/software_interface_double_gene_demo/delta_table.csv
results/software_interface_double_gene_demo/auc_summary.csv
results/software_interface_double_gene_demo/01_summary_dashboard.png
results/software_interface_double_gene_demo/02_true_vs_virtual_heatmap.png
results/software_interface_double_gene_demo/03_cell_state_umap.png
results/software_interface_double_gene_demo/04_auc_strong_response_roc.png
```

### 双敲整体评估

```text
results/norman_system_prior_metrics.csv
results/norman_combo_explanation_table.csv
results/norman_combo_metrics.csv
results/figures/norman_true_vs_pred_system_prior.png
results/figures/norman_r2_improvement_system_prior.png
results/figures/norman_combo_true_vs_pred_heatmap.png
results/figures/norman_cell_level_umap_multi_gene_before_after.png
```

### 双敲交互优化

```text
results/norman_double_interaction_metrics.csv
results/norman_double_interaction_predictions.csv
results/figures/norman_double_interaction_model_comparison.png
docs/norman_double_interaction_model.md
scripts/28_norman_double_interaction_model.py
```

### 方法代码

```text
vkx/preprocess.py
vkx/core.py
vkx/visualization.py
vkx/cli.py
```

## 11. 最终判断

单敲：当前已经可以作为比较清楚的可复用 demo。`STAT1` 示例方向一致性高，AUC 好，图也能直观看到真实/虚拟 KO 的差别。

双敲：已经可以支持，并且 `CEBPB+CEBPA` 这种组合效果很好。新增 interaction residual 后，52 个 Norman 双敲组合整体明显提升，说明“强约束 + 交互修正”是当前比自由生成模型更合适的优化方向。
