# Papalexi 小样本多模态实验记录

## 当前设置

- 数据：Papalexi 2021 ECCITE-seq，RNA + protein + sgRNA。
- 小样本抽样：每个 KO 最多 120 个细胞。
- 规模：3037 cells x 2000 genes，4 个蛋白读数。
- 输出：
  - RNA pathway delta。
  - protein delta。

## 关键观察

### 1. 通路信号成立

IFNG-JAK-STAT 通路变化有明确生物学方向。下降最强的 KO 包括 `STAT1`、`JAK2`、`STAT2`、`IFNGR2` 和 `IFNGR1`。

### 2. 蛋白读数提供强约束

PDL1 蛋白变化与 IFNG-JAK-STAT 通路变化高度相关：

- `delta_pathway_IFNG_JAK_STAT` vs `delta_protein_PDL1`: Pearson r ≈ 0.70。
- `delta_pathway_IMMUNE_CHECKPOINT` vs `delta_protein_PDL1`: Pearson r ≈ 0.68。

这支持“用多模态表型约束虚拟 KO 状态”的思路。

### 3. 简单 one-hot KO 不能泛化

普通 one-hot 留一法在未见 KO 上基本只能预测训练均值，各 pathway 的 R2 均为负。

### 4. 加入机制先验后有改善

加入 KO gene 的 pathway/module/protein-prior 特征后：

- IFNG-JAK-STAT 的 PLS joint baseline R2 ≈ 0.23。
- PDL1 protein delta 的 ridge joint baseline R2 ≈ 0.68。
- CD86 protein delta 的 ridge joint baseline R2 ≈ 0.32。
- CD366 protein delta 的 ridge joint baseline R2 ≈ 0.60。

说明蛋白表型比部分 RNA pathway 更稳定，也更适合小样本监督信号。

## 当前瓶颈

1. 当前 Papalexi 文件只保留 2000 个预处理基因，很多 pathway gene 缺失。
2. KO 先验仍是手工规则，缺少系统性网络特征。
3. 目前是 perturbation-level 聚合预测，还没有建 cell-state-conditioned model。
4. Papalexi 主要是单基因 KO，不能充分验证多基因组合泛化。

## 下一步建议

1. 用 full gene/raw 数据重新计算更完整的 pathway/TF activity。
2. 引入 Reactome、MSigDB、DoRothEA、STRING/PPI 或 TF-target 网络特征。
3. 建立 cell-level 条件模型：

   ```text
   baseline RNA pathway + protein state + KO prior -> KO pathway delta + protein delta
   ```

4. 在 Norman 2019 上做单基因训练、双基因组合 KO 测试。
5. 最后再把 ridge/PLS baseline 替换成 conditional VAE 或 flow matching。

## AUC 评估

AUC 不适合作为连续值回归的唯一指标，但适合回答“模型能不能识别强响应 KO”。当前用 `abs(delta) >= 0.15` 作为强响应阈值。

Papalexi 多模态联合模型中：

- IFNG-JAK-STAT 下降响应：PLS ROC-AUC ≈ 0.95。
- PDL1 蛋白下降响应：ridge ROC-AUC = 1.00。
- CD86 蛋白强响应：PLS ROC-AUC ≈ 0.98。

这说明模型对“哪些 KO 会造成强方向性响应”的排序能力比连续值精确回归更稳定。

## Norman 2019 组合扰动验证

已下载轻量预处理版 Norman 2019，规模为 27,658 cells x 2,000 genes。小样本抽样后为 7,969 cells。

设计：

- 训练：单基因扰动。
- 测试：双基因组合扰动。
- 目标：预测细胞程序分数变化。

全体 52 个双基因组合中，只有 12 个组合的两个基因都在单基因训练集中出现过。这个拆分很重要，因为未见基因组合需要网络/模块先验才能外推。

在 `all_genes_seen` 的 12 个组合上：

- Erythroid program：R2 ≈ 0.76。
- Granulocyte/apoptosis program：R2 ≈ 0.78。
- Pioneer TF program：R2 ≈ 0.55。

这支持“单基因扰动学习到的程序变化可以部分外推到双基因组合”。但 MAPK/TGFB 和 Pro-growth 仍然较差，提示组合效应并非简单加和。
