# Shotgun metagenomics of colorectal cancer

This multi-cohort workflow examines colorectal-cancer-associated microbial composition, diversity, taxonomic biomarkers, predicted functional potential, adenoma comparisons and external validation. Cohort identity is retained in the analyses to avoid treating study-specific measurements as interchangeable.

Run the workflow from this directory:

```text
Rscript run_pipeline.R
```

The runner discovers exactly eight numbered scripts under `scripts/`, executes them in order and stops if the expected set is incomplete or a script fails. The harmonised abundance tables and generated results are not included; numerical reproduction requires the original data and documented R dependencies.
