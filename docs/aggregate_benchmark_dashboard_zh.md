# VKX 横向 benchmark 汇总图说明

这个文档说明 `aggregate-benchmarks` 命令的作用。它解决的问题很简单：前面每次实验都会生成一组结果目录，如果只逐个打开，用户很难判断当前 VKX 到底有没有变好、哪一个版本更适合作为默认方法、还差在哪里。因此我们新增了一个总览入口，把多个正式 benchmark 和 hard-constrained generator 结果放到同一张图里。

## 输入是什么

用户不需要重新整理数据，只需要给一个结果目录，例如：

```bash
python -m vkx.cli aggregate-benchmarks --result-dirs results --out-dir results/aggregate_benchmark_report
```

`--result-dirs` 可以是一个总目录，也可以是多个结果目录，用逗号隔开。程序会自动寻找这些文件：

- `method_metric_comparison.csv`：正式横向 benchmark 的方法比较结果。
- `hard_generator_metrics.csv`：hard constraint 细胞级生成模型的 KO 预测结果。
- `hard_generator_report.md`：用于读取 generator 的 anchor method。

## 输出是什么

输出目录会包含：

- `01_formal_benchmark_leaderboard.png`：正式 benchmark 总览图。
- `02_hard_generator_leaderboard.png`：hard-constrained generator 总览图。
- `formal_method_metrics_aggregate.csv`：所有正式 benchmark 的原始汇总表。
- `formal_best_methods.csv`：每个 benchmark 目录中综合表现最好的方法。
- `generator_metrics_aggregate.csv`：所有 generator 的逐 KO 指标。
- `generator_summary.csv`：每个 generator 目录的平均表现。
- `aggregate_benchmark_report.md`：带图片引用和表格摘要的结果报告。

## 两张图怎么看

第一张图回答：在正式横向 benchmark 中，VKX 当前哪个版本最好？

- `AUC` 越高越好，表示真实 KO 改变较大的 feature 是否能被模型排到前面。
- `R2` 越高越好，表示预测变化幅度和真实变化幅度是否一致。
- `MAE` 越低越好，表示平均误差更小。
- 当前 adaptive/boosted 版本相对原始 constrained ensemble 有提升，说明 KO-response prior 和 hard anchor 对小样本 perturbation 场景是有帮助的。

第二张图回答：在 cell-level 条件生成阶段，hard constraint 附近的不确定性生成有没有改善？

- `R2` 看生成细胞的平均 KO response 是否接近真实 response。
- `Direction` 看预测方向是否一致。
- `MAE` 看数值误差。
- boosted anchor 的表现优于 ensemble anchor，说明轻量生成器应该继续围绕更强的 baseline 约束学习，而不是自由生成。

## 当前结论

当前 VKX 的强项是：小样本、可解释、多模态输入兼容、单敲/双敲、reference model、普通 10X prediction-only 应用，以及围绕 hard constraint 的细胞级虚拟 KO 生成。

当前仍需要补强的是：

- 与 scGen、CPA、GEARS、CellOT 等方法做更完整的同数据集横向比较。
- 接入更多有 perturbation 标签的数据，尤其是真正 RNA+ADT+ATAC 且带 KO 标签的数据。
- 继续增强双敲非线性交互，尤其是 MAPK/TGFB 这类难预测 program。
- 让 ATAC peak-level 图继续从“全局指标好看”变成“位点、motif、target gene 解释更清楚”。
- 在更多 perturbation 数据上训练轻量 neural generator，但仍保持 residual/PLS baseline 作为 hard constraint。

## 为什么这个入口重要

用户最关心的不是某个技术模块有没有实现，而是“这个方法到底有没有比之前更好”。`aggregate-benchmarks` 把分散实验收敛成两张主图：一张看正式 benchmark，一张看细胞级 generator。它会成为后续论文级 figure package 的主索引。
