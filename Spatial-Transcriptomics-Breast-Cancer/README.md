# Spatial transcriptomics of breast cancer

This workflow combines tissue morphology and spatial gene expression to characterise breast tumour architecture, cellular programmes, spatially variable genes and pathway enrichment. Results are interpreted at the tissue-section level and do not imply patient-level replication where it is unavailable.

Inspect the execution order without running an analysis:

```text
python run_pipeline.py --dry-run
```

Run all eight scripts with `python run_pipeline.py`. The runner uses the active Python interpreter, executes scripts from the project directory and stops at the first failure. Original spatial datasets and generated outputs are excluded from the repository and are required for numerical reproduction.
