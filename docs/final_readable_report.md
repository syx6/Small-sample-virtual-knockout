# 小样本多模态虚拟敲除方法：阶段性结果说明

这份文档面向不熟悉生信建模的人，解释当前原型做了什么、输入输出是什么、结果好在哪里、还差在哪里。

## 1. 一句话结论

当前原型已经证明：

> 用小样本多模态数据，加上通路/转录因子/PPI 等系统先验，可以比较好地识别哪些基因敲除会造成明显的通路或蛋白状态变化。

但是它还没有达到：

> 精确生成每一个细胞在任意多基因敲除后的完整状态。

换句话说，现在这个方法更像一个 **“机制解释型虚拟敲除筛选器”**，而不是最终版的 **“单细胞级虚拟细胞生成器”**。

当前最可靠的能力：

- 找出强响应 KO。
- 解释 KO 后哪些通路/蛋白变化。
- 对部分双基因组合 KO 做初步外推。
- 用多模态结果证明通路变化和蛋白表型确实有关联。

当前还需要增强的能力：

- 精确预测连续变化幅度。
- 处理复杂的非线性多基因组合效应。
- 真正做到 cell-level 条件生成。

---

## 2. 方法到底输入什么、输出什么？

### 输入

最小输入：

- 单细胞 RNA 表达矩阵。
- 每个细胞对应的扰动标签，例如 `STAT1`、`JAK2`、`KLF1+BAK1`、`ctrl`。
- control 或 negative-control 细胞标记。

推荐输入：

- RNA + protein/CITE-seq。
- RNA + ATAC 或 gene activity。
- 细胞类型、样本来源、批次等 metadata。
- 单基因或多基因 KO 标签。

系统先验：

- Reactome pathway。
- MSigDB Hallmark。
- TF-target。
- PPI hub 或 gene network。

### 输出

主要输出：

- 每个 KO 或组合 KO 后的通路活性变化。
- 可选的蛋白表型变化，例如 PDL1、CD86。
- 强响应 KO 排名。
- 多基因组合 KO 的预测结果。
- 该 KO 命中的通路、转录因子或 PPI 先验解释。

### 方法流程图

![方法输入输出](../results/figures/user_facing_input_output_diagram.png)

这张图可以这样理解：

1. 输入单细胞数据和 KO 标签。
2. 把细胞状态压缩成更容易解释的通路、TF、蛋白状态。
3. 把 KO 基因映射到 Reactome/MSigDB/TF-target/PPI 等先验。
4. 模型预测 KO 后细胞状态会怎样变化。
5. 输出通路变化、蛋白变化、强响应排名和机制解释。

---

## 3. 当前实验用了哪些数据？

### Papalexi 2021 ECCITE-seq

特点：

- 有 RNA。
- 有蛋白/CITE-seq。
- 有 CRISPR 扰动标签。
- 适合测试“多模态虚拟敲除”。

我们主要用它回答：

> 多模态信息是否真的能帮助解释 KO 后状态变化？

### Norman 2019 Perturb-seq

特点：

- 有真实双基因组合扰动。
- 主要是 RNA。
- 适合测试“多基因组合 KO 外推”。

我们主要用它回答：

> 只看单基因扰动训练，能不能预测双基因组合扰动？

---

## 4. 总体效果判断

![方法效果判断](../results/figures/user_facing_method_verdict.png)

这张图是最推荐先看的总览。

结论：

- 强响应 KO 识别：好。
- 连续变化幅度预测：中等。
- 双基因组合外推：部分有效。
- 小样本多模态思路：成立。

简单说：

> 这个方法现在能较好回答“哪些 KO 值得关注、影响哪些通路/蛋白”，但还不能完全精确回答“每个细胞的所有状态数值会变成多少”。

---

## 5. AUC 曲线：能不能识别强响应 KO？

![ROC 曲线](../results/figures/papalexi_roc_curves.png)

这张图是标准 ROC 曲线，不是简单柱状图。

它回答的问题是：

> 如果有很多 KO，模型能不能把真正强响应的 KO 排到前面？

图里虚线表示随机猜测。曲线越靠近左上角，说明模型越好。

当前结果：

- IFNG-JAK-STAT 通路下降：AUC 约 0.95。
- PDL1 蛋白下降：AUC 约 0.96-1.00。
- CD86 蛋白强变化：AUC 约 0.95-0.98。

解释：

> 模型在“筛选强响应 KO”这个任务上表现较好。

这对实际使用很重要，因为很多情况下用户首先想知道的是：

> 哪些 KO 最可能产生明显生物学效应？

---

## 6. 多模态是否真的有用？

![通路蛋白相关](../results/figures/papalexi_pathway_protein_correlation.png)

这张图展示 RNA 通路变化和蛋白变化之间的关系。

最关键的结果：

- IFNG-JAK-STAT 通路变化和 PDL1 蛋白变化相关性约 0.70。
- Immune checkpoint 通路变化和 PDL1 蛋白变化相关性约 0.68。

解释：

> RNA 通路变化不是孤立的，它和蛋白表型变化能对上。

这支持我们的核心想法：

> 小样本之所以有希望，是因为多模态数据提供了更高信息密度和更强的生物学约束。

---

## 7. 连续变化幅度预测准不准？

![Norman 真实预测散点](../results/figures/norman_true_vs_pred_system_prior.png)

这张图比单纯 R2 柱状图更重要。

怎么看：

- 横轴是真实双基因 KO 后的变化。
- 纵轴是模型预测变化。
- 虚线是完美预测线。
- 点越靠近虚线，说明预测越准。

当前结果：

- Granulocyte/apoptosis：预测较好。
- Erythroid：有一定趋势。
- MAPK/TGFB：预测较弱。
- Pro-growth：预测较弱。

解释：

> 模型能抓住一部分程序变化方向，但连续幅度预测仍不稳定。

这说明下一步需要更强的 cell-level 条件模型，而不是只在 KO 组平均水平上建模。

---

## 8. 系统先验有没有帮助？

![系统先验提升](../results/figures/norman_r2_improvement_system_prior.png)

这张图回答：

> Reactome/MSigDB/TF-target/PPI 这些系统先验有没有比简单加和更好？

横轴是 R2 改善：

- 大于 0：系统先验比 additive baseline 好。
- 小于 0：系统先验反而不如简单加和。

当前结果：

- Pro-growth 改善最大。
- MAPK/TGFB 有改善。
- Granulocyte/apoptosis 有改善。
- Erythroid 和 Pioneer TF 没有改善。

解释：

> 系统先验对部分未见基因组合有帮助，但不是所有任务都提升。

这很合理，因为有些组合接近简单加和，有些组合需要更复杂的非线性模型。

---

## 9. seen gene / unseen gene 是什么意思？

在 Norman 的多基因组合实验中，我们故意做了这样一个测试：

```text
训练：单基因扰动
测试：双基因组合扰动
```

因此：

- seen gene：这个基因在训练集中有单基因扰动样本。
- unseen gene：这个基因没有单基因训练样本，只在双基因组合中出现。
- seen combo：双基因组合里的两个基因都见过单基因扰动。
- unseen combo：双基因组合里至少一个基因没见过单基因扰动。

这不是说生物学上不知道这个基因，而是说：

> 模型训练时没有看过这个基因单独敲除会怎样。

为什么要这样做？

因为真正有价值的虚拟 KO 方法，应该能在数据不完整时借助系统先验进行外推。

---

## 10. 具体敲了什么基因？真实变化和预测有什么区别？

### 组合 KO 热图

![组合 KO 热图](../results/figures/norman_combo_true_vs_pred_heatmap.png)

这张图每一行是一个双基因 KO，例如：

- `AHR+KLF1`
- `CBL+UBASH3B`
- `CEBPB+CEBPA`
- `MAPK1+TGFBR2`

左边是真实变化，右边是模型预测变化。

红色表示该程序增强，蓝色表示该程序降低。

可以看到：

- `CEBPB+CEBPA` 真实增强 Granulocyte/apoptosis，模型也预测到这个方向。
- `AHR+KLF1` 真实增强 Erythroid，模型也预测增强，但幅度偏低。
- `MAPK1+TGFBR2` 真实 MAPK/TGFB 增强明显，但模型预测不足。
- `CBL+UBASH3B` 真实 Erythroid 增强很强，模型明显低估。

### 代表组合条形图

![代表组合真实 vs 预测](../results/figures/norman_representative_combo_bars.png)

这张图更适合逐个看例子：

- 绿色是真实双基因 KO 变化。
- 橙色是模型预测变化。
- 正值表示程序增强。
- 负值表示程序降低。

这张图说明：

> 模型能捕捉一些组合 KO 的主要方向，但对幅度大的组合容易低估。

---

## 11. 单细胞层面的虚拟敲除效果

之前很多图是 KO 组平均结果。下面这些图展示单细胞层面的状态分布。

### Papalexi 单细胞虚拟 KO

![Papalexi 单细胞虚拟 KO](../results/figures/papalexi_single_cell_virtual_ko_distributions.png)

三种颜色分别是：

- 灰色：真实 control 细胞。
- 绿色：真实 KO 细胞。
- 橙色：虚拟 KO 细胞。

虚拟 KO 细胞目前这样构造：

```text
virtual KO cell = control cell score + predicted KO delta
```

当前结果：

- `STAT1` KO 后 IFNG-JAK-STAT 通路明显下降，虚拟 KO 也预测下降，但幅度偏小。
- `JAK2` KO 后 IFNG-JAK-STAT 通路下降，虚拟 KO 抓到了方向。
- `IFNGR2` / `JAK2` KO 后 PDL1 蛋白下降，虚拟 KO 也能预测到下降趋势。

解释：

> 模型能捕捉部分单细胞分布移动方向，但还不是完整的单细胞生成模型。

### Norman 双基因组合单细胞分布

![Norman 单细胞组合 KO](../results/figures/norman_single_cell_virtual_combo_distributions.png)

当前结果：

- `AHR+KLF1`：真实 Erythroid 上升，虚拟 KO 也上升，但幅度偏低。
- `CEBPB+CEBPA`：真实 Granulocyte/apoptosis 上升，虚拟 KO 方向基本对。
- `MAPK1+TGFBR2`：真实 MAPK/TGFB 上升明显，虚拟 KO 预测不足。
- `CBL+UBASH3B`：真实 Erythroid 上升很强，虚拟 KO 明显低估。

---

## 12. 当前方法到底好不好？

### 好的地方

1. 强响应筛选能力好。  
   ROC 曲线显示 IFNG pathway、PDL1 protein、CD86 protein 的强响应识别效果较好。

2. 多模态思路成立。  
   RNA 通路变化和蛋白表型变化能对上。

3. 系统先验确实有帮助。  
   对部分未见基因组合，Reactome/MSigDB/TF-target/PPI 先验能改善外推。

4. 结果可解释。  
   输出不是一堆基因表达，而是通路变化、蛋白变化、命中的系统先验。

### 不足的地方

1. 连续幅度预测还不够准。  
   很多组合预测方向对，但幅度低估。

2. 复杂多基因组合仍然困难。  
   例如 `MAPK1+TGFBR2`、`CBL+UBASH3B` 等组合仍预测不足。

3. 单细胞虚拟 KO 还是简化版。  
   当前是 `control cell + predicted delta`，还不是最终的 cell-level 条件生成模型。

4. 使用的数据是轻量预处理版。  
   Papalexi 和 Norman 都只有约 2000 个基因，限制了更完整的 pathway/TF activity 计算。

---

## 13. 对外应该怎么介绍这个方法？

推荐表述：

> 这是一个面向小样本多模态单细胞数据的机制先验驱动虚拟敲除框架。它利用 RNA、蛋白、扰动标签和通路/网络先验，预测 KO 后的通路和蛋白状态变化，尤其适合强响应筛选、多基因组合扰动的初步外推和机制解释。

不建议现在这样说：

- 可以精确重建 KO 后全转录组。
- 可以可靠预测任意多基因组合。
- 已经是完整的单细胞虚拟细胞生成模型。

---

## 14. 下一步应该做什么？

下一步才进入真正的 cell-level 条件生成模型。

目标应该从：

```text
KO 组平均状态 + KO prior -> 平均 delta
```

升级为：

```text
单个 control cell state + KO genes + system prior
-> 该细胞 KO 后 pathway/protein state
```

推荐模型方向：

- conditional VAE。
- flow matching。
- lightweight conditional diffusion。
- graph-prior encoder + cell-state decoder。

输出仍然建议保持在 pathway/protein/TF activity 空间，而不是直接生成全基因表达。

这样更适合小样本，也更容易解释。

