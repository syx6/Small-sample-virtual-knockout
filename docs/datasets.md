# 推荐实践数据集

## 1. Papalexi 2021 ECCITE-seq

最推荐作为第一版原型。

- 数据类型：RNA + ADT 蛋白 + sgRNA。
- 扰动类型：CRISPR KO / sgRNA perturbation。
- 场景：受刺激 THP-1 细胞，研究免疫检查点调控。
- 优点：天然多模态，扰动数量适中，适合小样本方法开发。
- 局限：主要是单基因 KO，不是真实多基因组合 KO。
- 获取方式：`pertpy.data.papalexi_2021()`。

适合任务：

- 输入 RNA/蛋白状态 + KO gene，预测 KO 后通路分数。
- 比较 RNA-only、protein-only、RNA+protein 三种条件。
- 把单基因 KO 作为多基因 KO 模型的预训练/校准数据。

## 2. Norman 2019 Perturb-seq

推荐作为多基因组合 KO 的验证数据。

- 数据类型：scRNA-seq + guide assignment。
- 扰动类型：单基因和双基因组合 perturbation。
- 优点：有真实组合扰动，适合测试组合泛化。
- 局限：不是多模态。
- 获取方式：`pertpy.data.norman_2019()` 或 `pertpy.data.norman_2019_raw()`。

适合任务：

- 训练单基因 KO，预测双基因 KO。
- 比较 additive model 与 nonlinear model。
- 检查通路分数是否比基因表达更稳。

## 3. Frangieh 2021 Perturb-CITE-seq

适合作为第二阶段扩展。

- 数据类型：RNA + 20 个表面蛋白 + perturbation。
- 扰动规模：约 218,000 cells，约 750 perturbations。
- 优点：大规模多模态扰动数据。
- 局限：对初始实验偏大，建议只抽小样本子集。
- 获取方式：`pertpy.data.frangieh_2021()`。

适合任务：

- 小样本抽样实验。
- 跨扰动泛化。
- 与 Papalexi 数据做外部验证。

## 4. Perturb-ATAC / CAT-ATAC / Multiome Perturb-seq

适合作为 RNA+ATAC 方向的后续扩展。

- 数据类型：ATAC 或 RNA+ATAC + CRISPR guide。
- 优点：更贴近“多模态调控状态”。
- 局限：下载、预处理、peak-to-gene link、gene activity 计算更重。

适合任务：

- KO 后 chromatin accessibility 或 gene activity 变化。
- 用 ATAC 推断 TF activity，再预测通路状态。
- 建立 RNA pathway + ATAC TF activity 的双输出模型。
