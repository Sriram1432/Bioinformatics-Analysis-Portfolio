#!/usr/bin/env Rscript
root <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), winslash = "/", mustWork = TRUE)
Sys.setenv(PROJECT3_ROOT = root)
scripts <- file.path(root, sprintf("%02d_%s.R", 1:8, c("prepare_cohorts", "expression_integrity", "cancer_specific_cox", "pan_cancer_analysis", "adjusted_models_diagnostics", "internal_validation", "pathway_analysis", "figures")))
missing <- scripts[!file.exists(scripts)]
if (length(missing)) stop("Missing pipeline scripts: ", paste(missing, collapse = ", "))
for (script in scripts) {
    status <- system2(file.path(R.home("bin"), "Rscript"), script)
    if (!identical(status, 0L)) stop("Pipeline stopped after failure in ", basename(script))
}
