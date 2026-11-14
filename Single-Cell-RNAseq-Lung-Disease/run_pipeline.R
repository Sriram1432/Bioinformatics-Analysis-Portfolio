#!/usr/bin/env Rscript
command_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(command_file)) stop("Run this pipeline with Rscript run_pipeline.R", call. = FALSE)
root <- normalizePath(dirname(sub("^--file=", "", command_file[[1]])), winslash = "/", mustWork = TRUE)
Sys.setenv(PROJECT1_ROOT = root)
scripts <- file.path(root, sprintf("%02d_%s.R", 1:8, c("prepare_data", "cellular_composition", "primary_pseudobulk", "pathway_enrichment", "macrophage_subgroups", "atii_heterogeneity", "robustness_diagnostics", "figures")))
missing <- scripts[!file.exists(scripts)]
if (length(missing)) stop("Missing pipeline scripts: ", paste(missing, collapse = ", "))
for (index in seq_along(scripts)) {
    script <- scripts[[index]]
    message(sprintf("[%d/%d] %s", index, length(scripts), basename(script)))
    status <- system2(file.path(R.home("bin"), "Rscript"), script)
    if (!identical(status, 0L)) stop("Pipeline stopped after failure in ", basename(script))
}
