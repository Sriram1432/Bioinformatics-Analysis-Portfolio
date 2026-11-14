# Single-cell RNA-seq of human lung disease

This workflow analyses the GSE136831 human lung atlas with clinical samples as the biological replication unit. It covers data preparation, cellular composition, pseudobulk differential expression, pathway enrichment, macrophage subgroups, ATII heterogeneity, robustness checks and figure generation.

Run the complete workflow from this directory:

```text
Rscript run_pipeline.R
```

The runner executes the eight numbered scripts in order and stops at the first error. Source data and generated results are not included in this repository; the original inputs and the R packages listed in the repository dependency files are required for numerical reproduction.
