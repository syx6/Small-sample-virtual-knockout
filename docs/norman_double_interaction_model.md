# Norman double-KO interaction model optimization

This experiment adds a constrained interaction correction for double-gene knockouts.

Method:

```text
single-gene additive prediction
+ prior-based gene-gene interaction residual
= optimized double-KO prediction
```

The correction is trained by leave-one-combo-out evaluation across 52 Norman double-gene KO combinations.

Summary:

               model          subset  mean_mae  mean_r2  mean_roc_auc  mean_pr_auc
interaction_residual      all_combos     0.076    0.617         0.894        0.845
interaction_residual  all_genes_seen     0.105    0.026         0.872        0.830
interaction_residual has_unseen_gene     0.067    0.697         0.895        0.862
single_gene_additive      all_combos     0.150    0.008         0.707        0.682
single_gene_additive  all_genes_seen     0.121   -0.083         0.982        0.900
single_gene_additive has_unseen_gene     0.159   -0.101         0.638        0.586

Key output files:

- `results/norman_double_interaction_metrics.csv`
- `results/norman_double_interaction_predictions.csv`
- `results/figures/norman_double_interaction_model_comparison.png`

Interpretation:

- If `interaction_residual` improves R2 or ROC-AUC over `single_gene_additive`, then explicit gene-gene interaction features help.
- If it improves some programs but worsens others, the interaction layer should be used selectively rather than as a universal replacement.
