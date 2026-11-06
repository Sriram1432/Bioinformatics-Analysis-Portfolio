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
stopifnot(length(scripts) == 8L)
for (script in scripts) {
    local_env <- new.env(parent = globalenv())
    sys.source(script, envir = local_env)
}
