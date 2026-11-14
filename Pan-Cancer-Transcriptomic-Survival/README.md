# Pan-cancer transcriptomic survival analysis

This workflow evaluates cancer-specific and pan-cancer associations between transcript abundance and survival. It includes cohort preparation, expression checks, Cox modelling, clinical adjustment, proportional-hazards diagnostics, internal validation and pathway analysis.

Run the workflow from this directory:

```text
Rscript run_pipeline.R
```

The runner executes the eight numbered scripts in dependency order and stops if a script fails. TCGA source data, intermediate objects and generated figures are not distributed here. Reproducing the numerical results requires those inputs and the package versions documented at repository level.
