#!/usr/bin/env Rscript
configured_root <- Sys.getenv("PROJECT4_ROOT", unset = "")
command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_root <- if (length(command_file)) {
    dirname(sub("^--file=", "", command_file[[1]]))
} else {
    getwd()
}
root <- if (nzchar(configured_root)) configured_root else script_root
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
Sys.setenv(PROJECT4_ROOT = root)
scripts <- sort(list.files(file.path(root, "scripts"), pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
if (length(scripts) != 8L) {
    stop("Expected 8 numbered analysis scripts but found ", length(scripts), " in ", file.path(root, "scripts"), call. = FALSE)
}
for (index in seq_along(scripts)) {
    script <- scripts[[index]]
    message(sprintf("[%d/%d] %s", index, length(scripts), basename(script)))
    local_env <- new.env(parent = globalenv())
    sys.source(script, envir = local_env)
}
