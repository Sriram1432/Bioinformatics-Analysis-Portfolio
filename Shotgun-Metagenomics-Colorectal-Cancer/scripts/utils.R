project_root <- function() {
    root <- Sys.getenv("PROJECT4_ROOT", unset = getwd())
    normalizePath(root, mustWork = TRUE)
}

require_packages <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing))
        stop("Install required R packages: ", paste(missing, collapse = ", "))
}

read_results <- function(...) {
    path <- file.path(project_root(), "results", ...)
    if (!file.exists(path))
        stop("Required result not found: ", path)
    data.table::fread(path, showProgress = FALSE)
}
