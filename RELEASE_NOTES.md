# 当前版本说明

本版本是“小样本多模态虚拟敲除”方法的可复用原型包，重点面向：

- 单基因虚拟敲除
- 双基因虚拟敲除与 interaction residual 分析
- RNA-only 输入自动转换为 pathway/program score
- RNA + ADT/protein、RNA + ATAC/gene activity、RNA + chromVAR/motif activity 等多模态输入
- 小样本 perturbation 数据下的 hard-constrained residual/PLS baseline
- ROC/AUC、heatmap、UMAP、summary dashboard、ATAC peak-level 可视化

## 仓库包含内容

- `vkx/`: 可复用 Python 软件接口与命令行入口
- `scripts/`: 复现实验、生成报告和整理图件的脚本
- `docs/`: 中文方法说明、实验总结、适用范围、SCI 论文草稿
- `data/priors/`: Reactome/MSigDB-like pathway、TF-target、PPI 等小型网络先验文件
- `results/user_facing_figures/`: 面向用户阅读的核心结果图
- `README.md`: 用户输入、输出和快速运行说明

## 仓库不包含内容

以下内容被 `.gitignore` 排除，不随 GitHub 版本上传：

- 原始单细胞数据、h5ad/h5mu/mtx 文件
- 下载得到的公共数据压缩包
- 大型中间结果表
- Python 虚拟环境
- 训练缓存和临时模型文件

这样做是为了让 GitHub 仓库保持轻量、可读、可维护。用户需要复现实验时，可根据 `README.md` 和 `docs/` 中的命令重新下载或放入自己的数据。

## 当前方法定位

这个版本不是从零训练的自由 diffusion/VAE 生成模型，而是把 residual/PLS 方向模型作为 hard constraint。生成层只在这个方向附近表达不确定性，因此更适合小样本、多模态、有系统先验的数据场景。

主要优势：

- 对小样本更稳
- 输入更接近普通用户的真实数据格式
- 可解释性强，输出以 pathway/program/protein/ATAC 状态变化为核心
- 可视化完整，便于判断真实 KO 与虚拟 KO 的差异

主要限制：

- 没有 perturbation 标签的普通 10X 数据只能做 reference application，不能在该数据内部验证准确性
- 复杂双敲非线性效应仍需要更多真实 double-KO 数据
- 真正 RNA + ADT + ATAC + perturbation 标签的三模态 benchmark 仍是后续最重要的验证方向
