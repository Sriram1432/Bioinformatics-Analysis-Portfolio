# Validate transcript identifiers, missingness and expression structure.  Creates cancer-specific
# discovery objects, documents structured missingness and checks residual expression patterns
# before survival modelling.  Run with: Rscript 02_expression_integrity.R

options(stringsAsFactors = FALSE)

resolve_project_root <- function() {
    configured <- Sys.getenv("PROJECT3_ROOT", unset = "")
    if (nzchar(configured))
        return(normalizePath(configured, winslash = "/", mustWork = TRUE))
    command <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(command))
        return(normalizePath(dirname(sub("^--file=", "", command[[1]])), winslash = "/", mustWork = TRUE))
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- resolve_project_root()


# ---------------------------------------------------------------------------- Feature identity and
# cancer-specific eligibility
# ----------------------------------------------------------------------------

options(width = 185, scipen = 5)
root <- PROJECT_ROOT
out <- file.path(root, "Results_Rebuild", "Stage03")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "matrixStats", "ggplot2", "ragg")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ",
    pkg)
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
inputs <- list.files(file.path(root, "Data"), recursive = TRUE, full.names = TRUE)
getfile <- function(name) {
    z <- inputs[basename(inputs) == name]
    if (length(z) != 1L)
        stop("Exactly one file required: ", name, "; got ", length(z))
    z
}
expfile <- getfile("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
clifile <- getfile("Survival_SupplementalTable_S1_20171025_xena_sp")
phefile <- getfile("TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
cli <- fread(clifile, showProgress = FALSE)
phe <- fread(phefile, showProgress = FALSE)
req <- c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time")
if (length(setdiff(req, names(cli)))) stop("Missing clinical fields: ", paste(setdiff(req, names(cli)),
    collapse = ", "))
if (!all(c("sample", "sample_type_id") %in% names(phe))) stop("Missing phenotype fields")
setnames(cli, req, c("sample_id", "patient_id", "cancer", "OS", "OS_time"))
setnames(phe, "sample", "sample_id")
cli[, `:=`(sample_id = as.character(sample_id), patient_id = as.character(patient_id), cancer = as.character(cancer),
    OS = suppressWarnings(as.numeric(OS)), OS_time = suppressWarnings(as.numeric(OS_time)))]
phe[, sample_id := as.character(sample_id)]
# Phenotype is decisive for primary-specimen assignment; mismatched barcodes are excluded.
primary_ids <- unique(phe[!is.na(sample_type_id) & sample_type_id == 1L, sample_id])
header <- names(fread(expfile, nrows = 0L, check.names = FALSE, showProgress = FALSE))
if (header[1L] != "sample") stop("Unexpected expression feature field")
expnames <- header[-1L]
dup_expr <- unique(expnames[duplicated(expnames) | duplicated(expnames, fromLast = TRUE)])
selected <- cli[cancer %chin% cancers & sample_id %chin% primary_ids]
selected[, valid := !is.na(sample_id) & !is.na(patient_id) & patient_id == substr(sample_id, 1, 12) &
    OS %in% c(0, 1) & is.finite(OS_time) & OS_time > 0]
repeated <- selected[, .N, by = patient_id][N > 1L, patient_id]
matched <- selected[valid == TRUE & !patient_id %chin% repeated & !sample_id %chin% dup_expr & sample_id %chin%
    expnames]
if (anyDuplicated(matched$patient_id) || anyDuplicated(matched$sample_id)) stop("Nonunique analytic identities")
if (nrow(matched) != 4718L || sum(matched$OS == 1L) != 1519L) stop("Cohort differs from checked Stage02: ",
    nrow(matched), " patients / ", sum(matched$OS == 1L), " events")
message("Stage03: importing 4,718 expression profiles ...")
colpos <- c(1L, match(matched$sample_id, header))
if (anyNA(colpos)) stop("Expression matching failed")
dt <- fread(expfile, select = sort(colpos), check.names = FALSE, showProgress = FALSE)
feature <- as.character(dt[[1L]])
colids <- setdiff(names(dt), names(dt)[1L])
if (!setequal(colids, matched$sample_id)) stop("Expression columns mismatched")
if (!all(vapply(dt[, ..colids], is.numeric, logical(1)))) stop("Non-numeric expression column")
expr <- as.matrix(dt[, ..colids])
storage.mode(expr) <- "double"
rownames(expr) <- feature
rm(dt)
invisible(gc(FALSE))
matched <- matched[match(colnames(expr), sample_id)]
if (!identical(as.character(matched$sample_id), colnames(expr))) stop("Clinical-expression order mismatch")
# Never collapse duplicated historical names without a verified one-to-one map.
dupfeat <- unique(feature[duplicated(feature) | duplicated(feature, fromLast = TRUE)])
badid <- is.na(feature) | !nzchar(trimws(feature)) | feature %chin% dupfeat
idtab <- data.table(Feature_ID = feature, ID_class = ifelse(grepl("^[0-9]+$", feature), "Numeric_Entrez_like",
    "Symbol_like"), Ambiguous_or_empty = badid)
fwrite(idtab, file.path(out, "Stage03_source_feature_id_inventory.csv"))
fwrite(data.table(Duplicate_feature_ID = dupfeat), file.path(out, "Stage03_ambiguous_feature_ids.csv"))
validexpr <- expr[!badid, , drop = FALSE]
rm(expr)
invisible(gc(FALSE))
# Optional annotation: only exact unique ENTREZID -> SYMBOL via installed versioned Bioconductor
# database.
numids <- unique(rownames(validexpr)[grepl("^[0-9]+$", rownames(validexpr))])
mapping <- data.table(Source_ID = numids, HGNC_symbol = NA_character_, Mapping_status = "UNMAPPED_NO_LOCAL_DATABASE")
annotation_version <- "Not installed"
if (length(numids) && requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi",
    quietly = TRUE)) {
    annotation_version <- as.character(utils::packageVersion("org.Hs.eg.db"))
    mm <- tryCatch(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = numids, columns = "SYMBOL",
        keytype = "ENTREZID"), error = function(e) NULL)
    if (!is.null(mm)) {
        mm <- as.data.table(mm)
        setnames(mm, c("ENTREZID", "SYMBOL"), c("Source_ID", "HGNC_symbol"))
        mm <- mm[!is.na(HGNC_symbol) & nzchar(HGNC_symbol)]
        counts <- mm[, .(Unique_symbols = uniqueN(HGNC_symbol), HGNC_symbol = if (uniqueN(HGNC_symbol) ==
            1L)
            unique(HGNC_symbol)[1L] else NA_character_), by = Source_ID]
        mapping <- merge(mapping[, .(Source_ID)], counts, by = "Source_ID", all.x = TRUE, sort = FALSE)
        mapping[, Mapping_status := fifelse(Unique_symbols == 1L, "UNIQUE_LOCAL_MAP", "UNMAPPED_OR_AMBIGUOUS")]
    }
}
fwrite(mapping, file.path(out, "Stage03_numeric_identifier_mapping_diagnostic.csv"))
# Mapping table is diagnostic only; do NOT rename rows or merge with existing symbols.  Cohort
# filtering is feature-specific and cancer-specific: no imputation.
stats <- list()
patterns <- list()
sets <- list()
tech <- list()
for (cc in cancers) {
    message("Stage03: ", cc)
    idx <- which(matched$cancer == cc)
    cl <- copy(matched[idx])
    cl[, valid := NULL]
    mat <- validexpr[, idx, drop = FALSE]
    if (!identical(colnames(mat), as.character(cl$sample_id)))
        stop(cc, ": alignment failure")
    missing <- !is.finite(mat)
    missing_per_sample <- colSums(missing)
    patterns[[cc]] <- data.table(Cancer = cc, Patient = cl$patient_id, Event = cl$OS, Missing_features = missing_per_sample)
    complete <- rowSums(missing) == 0L
    m <- mat[complete, , drop = FALSE]
    vary <- matrixStats::rowMaxs(m) > matrixStats::rowMins(m)
    eligible <- m[vary, , drop = FALSE]
    if (any(!is.finite(eligible)) || anyDuplicated(rownames(eligible)))
        stop(cc, ": invalid final matrix")
    obj <- list(cancer = cc, expression = eligible, clinical = cl, feature_ids = rownames(eligible),
        scale = "UCSC Xena PanCancer processed expression (transformed/batch-adjusted); not raw counts",
        selection = "Phenotype-primary, unique patient/sample, OS_time>0, complete finite variable features",
        annotation = "Original source feature IDs retained; numeric identifiers not renamed without verified unique map")
    saveRDS(obj, file.path(out, paste0("TCGA_", cc, "_discovery.rds")), compress = FALSE)
    stats[[cc]] <- data.table(Cancer = cc, Patients = nrow(cl), Deaths = sum(cl$OS == 1), Censored = sum(cl$OS ==
        0), Features_source = nrow(mat), Features_complete = sum(complete), Features_eligible = nrow(eligible),
        Features_incomplete = sum(!complete), Features_constant = sum(!vary), Samples_with_missing = sum(missing_per_sample >
            0), Median_missing_per_sample = median(missing_per_sample))
    sets[[cc]] <- rownames(eligible)
    # Identify structured missingness pattern without claiming why it arose.
    if (any(missing_per_sample > 0)) {
        missgroup <- missing_per_sample > 0
        complete_outcomes <- table(factor(cl$OS, levels = c(0, 1)), factor(missgroup, levels = c(FALSE,
            TRUE)))
        tech[[cc]] <- data.table(Cancer = cc, Missingness_group = c("No_missing_features", "Any_missing_features"),
            N = as.integer(colSums(complete_outcomes)), Events = as.integer(complete_outcomes["1", ]))
    }
    rm(mat, m, eligible, missing, cl)
    invisible(gc(FALSE))
}
stat <- rbindlist(stats)
samp <- rbindlist(patterns)
tgroups <- if (length(tech)) rbindlist(tech) else data.table()
shared <- Reduce(intersect, sets)
if (length(shared) != 17465L) warning("Common feature count differs from prior 17465: ", length(shared))
fwrite(stat, file.path(out, "Stage03_cohort_and_feature_eligibility.csv"))
fwrite(samp, file.path(out, "Stage03_sample_missingness.csv"))
fwrite(tgroups, file.path(out, "Stage03_structured_missingness_outcomes.csv"))
fwrite(data.table(Feature_ID = shared), file.path(out, "Stage03_common_feature_universe.csv"))
# Modern, standalone figures: each addresses a specific data-integrity question.
pal <- c(Complete = "#236B8E", Incomplete = "#E4984B")
plotdf <- melt(stat[, .(Cancer, Features_complete, Features_incomplete)], id.vars = "Cancer", variable.name = "Class",
    value.name = "N")
plotdf[, Class := factor(ifelse(Class == "Features_complete", "Complete", "Incomplete"), levels = c("Complete",
    "Incomplete"))]
plotdf[, Cancer := factor(Cancer, levels = rev(cancers))]
p1 <- ggplot2::ggplot(plotdf, ggplot2::aes(Cancer, N, fill = Class)) + ggplot2::geom_col(width = 0.67) +
    ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = pal) + ggplot2::scale_y_continuous(labels = function(x) format(x,
    big.mark = ",", trim = TRUE), expand = ggplot2::expansion(mult = c(0, 0.04))) + ggplot2::labs(title = "Feature completeness differs by TCGA cancer",
    subtitle = "Nonfinite source measurements before per-cancer complete-feature selection", x = NULL,
    y = "Number of features", fill = NULL, caption = "Source identifiers with missing or duplicated names were excluded before calculation.") +
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 20), axis.text = ggplot2::element_text(colour = "black", face = "bold"),
    axis.title = ggplot2::element_text(colour = "black", face = "bold"), legend.position = "top", panel.grid.major.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(colour = "black"))
ragg::agg_webp(file.path(figdir, "Stage03_feature_completeness_by_cancer.webp"), width = 3200, height = 2000,
    res = 300, quality = 95)
print(p1)
grDevices::dev.off()
# Sample-level structured missingness: show exact groups and number of patients.
groups <- samp[, .(Patients = .N, Deaths = sum(Event == 1)), by = .(Cancer, Missing_features)]
groups[, Cancer := factor(Cancer, levels = rev(cancers))]
p2 <- ggplot2::ggplot(groups, ggplot2::aes(x = Missing_features, y = Cancer, size = Patients, colour = Deaths/Patients)) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::scale_size_area(max_size = 17, breaks = c(2, 50, 150,
    300, 500, 1000)) + ggplot2::scale_colour_gradient(low = "#5CC9C0", high = "#173E75", limits = c(0,
    1), name = "Observed death fraction") + ggplot2::scale_x_continuous(breaks = sort(unique(groups$Missing_features)),
    labels = function(x) format(x, big.mark = ",")) + ggplot2::labs(title = "Structured expression missingness across specimens",
    subtitle = "Circle area reflects patients sharing the same number of unmeasured features", x = "Nonfinite features per specimen",
    y = NULL, size = "Patients", caption = "Observed death fractions are descriptive; neither the cause nor outcome mechanism of missingness is established.") +
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 20), axis.text = ggplot2::element_text(colour = "black", face = "bold"),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right", plot.caption = ggplot2::element_text(colour = "black"))
ragg::agg_webp(file.path(figdir, "Stage03_specimen_missingness_bubble_plot.webp"), width = 3500, height = 1900,
    res = 300, quality = 95)
print(p2)
grDevices::dev.off()
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 03 \342\200\224 REBUILD SCIENTIFIC CHECKPOINT\n====================================================================\n")
    cat("1. FINAL ANALYTIC COHORTS / FEATURE ELIGIBILITY\n")
    print(stat)
    cat("\n2. STRUCTURED MISSINGNESS AND EVENT DISTRIBUTION\n")
    print(tgroups)
    cat("\n3. MISSINGNESS GROUPS BY CANCER\n")
    print(groups[order(Cancer, Missing_features)])
    cat("\n4. FEATURE IDENTIFIER AUDIT\n")
    print(idtab[, .(N = .N), by = ID_class])
    cat("Duplicated source IDs: ", paste(dupfeat, collapse = ", "), "\n")
    cat("Source numeric IDs: ", length(numids), " | local annotation version: ", annotation_version,
        "\n", sep = "")
    print(mapping)
    cat("\n5. PAN-CANCER CORRESPONDENCE\n")
    cat("Shared feature IDs: ", length(shared), " | patients: ", sum(stat$Patients), " | deaths: ", sum(stat$Deaths),
        "\n", sep = "")
    cat("\n6. SCIENTIFIC LIMITATIONS / NEXT CHECKPOINT\n")
    cat("No imputation, no raw-count assumptions, no silent ID remapping. Feature missingness is structured in some cancers.\n")
    cat("Gene identity mapping is DIAGNOSTIC; source names retained in output; numeric IDs must not be assumed to be symbols.\n")
    cat("Whole-cancer complete-feature filtering may remove informative genes when a subset of samples has technical missingness.\n")
    cat("Any differences in event proportions by missingness group motivate an explicit technical confounding sensitivity analysis.\n")
    cat("Expression scale and batch effects require the Stage04 assessment; do not infer batch origin from these patterns alone.\n")
    cat("\n7. FIGURES\n", file.path(figdir, "Stage03_feature_completeness_by_cancer.webp"), "\n", file.path(figdir,
        "Stage03_specimen_missingness_bubble_plot.webp"), "\n")
    cat("\n8. OUTPUT\n", out, "\n====================================================================\n")
})
writeLines(report, file.path(out, "Stage03_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Structured
# missingness assessment
# ----------------------------------------------------------------------------

options(width = 180, scipen = 5)
root <- PROJECT_ROOT
out <- file.path(root, "Results_Rebuild", "Stage03B")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "survival", "ggplot2", "ragg")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install required package: ",
    pkg)
suppressPackageStartupMessages(library(data.table))
set.seed(20261102)
cancers <- c("COAD", "OV", "UCEC")
require_file <- function(p) {
    if (!file.exists(p))
        stop("Missing required file: ", p)
    p
}
files <- list.files(file.path(root, "Data"), recursive = TRUE, full.names = TRUE)
expfile <- files[basename(files) == "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"]
if (length(expfile) != 1L) stop("Exactly one original expression file expected")
header <- names(fread(expfile, nrows = 0L, check.names = FALSE, showProgress = FALSE))
dup <- unique(header[-1L][duplicated(header[-1L]) | duplicated(header[-1L], fromLast = TRUE)])
# Stage03 has already guaranteed unique primary tumour identifiers and OS > 0.
cl <- list()
for (cc in cancers) {
    p <- require_file(file.path(root, "Results_Rebuild", "Stage03", paste0("TCGA_", cc, "_discovery.rds")))
    d <- readRDS(p)
    cl[[cc]] <- as.data.table(d$clinical)
    if (!identical(as.character(cl[[cc]]$sample_id), colnames(d$expression)))
        stop("Stage03 clinical\342\200\223expression mismatch: ", cc)
    if (anyDuplicated(cl[[cc]]$sample_id) || any(cl[[cc]]$sample_id %in% dup))
        stop("Ambiguous sample in ", cc)
    rm(d)
}
all_ids <- unlist(lapply(cl, function(x) as.character(x$sample_id)), use.names = FALSE)
if (anyDuplicated(all_ids)) stop("Duplicate samples across selected cohorts")
pos <- match(all_ids, header)
if (anyNA(pos)) stop("Some Stage03 samples are missing from original expression")
message("Stage03B: reading original expression for ", length(all_ids), " samples once ...")
dt <- fread(expfile, select = sort(c(1L, pos)), check.names = FALSE, showProgress = FALSE)
ids <- as.character(dt[[1L]])
if (anyNA(ids)) stop("Missing source gene identifier")
dup_feature <- unique(ids[duplicated(ids) | duplicated(ids, fromLast = TRUE)])
keep <- !(ids %in% dup_feature) & nzchar(trimws(ids))
colids <- names(dt)[-1L]
m <- as.matrix(dt[, ..colids])
storage.mode(m) <- "double"
rownames(m) <- ids
rm(dt)
invisible(gc(FALSE))
m <- m[keep, all_ids, drop = FALSE]
if (!identical(colnames(m), all_ids)) stop("Columns not aligned after raw import")
# Unambiguous original identifiers only, and actual source measurement missingness.
stat <- list()
patterns <- list()
genesets <- list()
byfeature <- list()
groups <- list()
for (cc in cancers) {
    x <- m[, as.character(cl[[cc]]$sample_id), drop = FALSE]
    bad <- !is.finite(x)
    nbad <- colSums(bad)
    tab <- table(nbad)
    patterns[[cc]] <- data.table(Cancer = cc, Missing_features = as.integer(names(tab)), Patients = as.integer(tab))
    groups[[cc]] <- nbad > 0
    # Which genes are missing in ANY of the samples, and in what fraction?
    rate <- rowMeans(bad)
    genesets[[cc]] <- rownames(x)[rate > 0]
    byfeature[[cc]] <- data.table(Feature_ID = rownames(x), Cancer = cc, Missing_rate = rate)
    stat[[cc]] <- data.table(Cancer = cc, Patients = ncol(x), Genes = nrow(x), Genes_any_missing = sum(rate >
        0), Genes_all_missing = sum(rate == 1), Genes_complete = sum(rate == 0), Genes_ge90pct_observed = sum(rate <=
        0.1), Genes_ge80pct_observed = sum(rate <= 0.2), Affected_patients = sum(nbad > 0), Distinct_sample_missingness_counts = length(unique(nbad)),
        Distinct_nonzero_missingness_patterns = length(unique(nbad[nbad > 0])))
    rm(x, bad)
    invisible(gc(FALSE))
}
stat <- rbindlist(stat)
pat <- rbindlist(patterns)
rates <- rbindlist(byfeature)
setorder(stat, Cancer)
# Compare EXACT feature sets, not just the number of missing genes.
set_comps <- rbindlist(lapply(combn(cancers, 2, simplify = FALSE), function(pair) {
    a <- genesets[[pair[1]]]
    b <- genesets[[pair[2]]]
    data.table(Cancer_A = pair[1], Cancer_B = pair[2], Missing_A = length(a), Missing_B = length(b),
        Shared_missing = length(intersect(a, b)), Only_A = length(setdiff(a, b)), Only_B = length(setdiff(b,
            a)), Jaccard = if (length(union(a, b)))
            length(intersect(a, b))/length(union(a, b)) else NA_real_)
}))
# Exact sample-level missing signatures (hash-free: group by TRUE/FALSE pattern equality).  For
# this dataset the number of patterns should be small. Identify and compare.
signatures <- list()
for (cc in cancers) {
    x <- m[, as.character(cl[[cc]]$sample_id), drop = FALSE]
    miss_counts <- colSums(!is.finite(x))
    ucounts <- sort(unique(miss_counts))
    sigs <- lapply(ucounts, function(ct) {
        ix <- which(miss_counts == ct)
        ref <- !is.finite(x[, ix[1]])
        same <- if (length(ix) > 1)
            all(vapply(ix[-1], function(j) identical(!is.finite(x[, j]), ref), logical(1))) else TRUE
        data.table(Cancer = cc, Missing_features = ct, Patients = length(ix), Identical_missing_gene_set_within_group = same)
    })
    signatures[[cc]] <- rbindlist(sigs)
}
signatures <- rbindlist(signatures)
# Compare old conservative eligibility with per-gene 90%/80% observation.
ret <- melt(stat[, .(Cancer, Genes, Genes_complete, Genes_ge90pct_observed, Genes_ge80pct_observed)],
    id.vars = c("Cancer", "Genes"), variable.name = "Strategy", value.name = "Genes_retained")
ret[, Strategy := factor(Strategy, levels = c("Genes_complete", "Genes_ge90pct_observed", "Genes_ge80pct_observed"),
    labels = c("100% measured", ">=90% measured", ">=80% measured"))]
ret[, Pct_of_source := 100 * Genes_retained/Genes]
# Clinical association: adjusted missingness indicator, NEVER impute source expression.
clean_chr <- function(x) {
    z <- toupper(trimws(as.character(x)))
    z[z %in% c("", "NA", "N/A", "NULL", "[NOT AVAILABLE]", "[NOT EVALUATED]", "[DISCREPANCY]", "UNKNOWN",
        "NOT REPORTED")] <- NA_character_
    z
}
major_stage <- function(x) {
    z <- gsub("^STAGE[[:space:]]*", "", clean_chr(x))
    y <- rep(NA_character_, length(z))
    y[grepl("^IV", z)] <- "IV"
    y[grepl("^III", z) & is.na(y)] <- "III"
    y[grepl("^II", z) & is.na(y)] <- "II"
    y[grepl("^I($|A|B|C|[[:space:]])", z) & is.na(y)] <- "I"
    y
}
clinical_rows <- list()
cov_rows <- list()
km_rows <- list()
for (cc in cancers) {
    cdat <- cl[[cc]]
    n <- nrow(cdat)
    bad <- groups[[cc]]
    dat <- data.frame(time = as.numeric(cdat$OS_time), event = as.integer(cdat$OS), Missing_group = as.integer(bad))
    if (any(!is.finite(dat$time) | dat$time <= 0) | !all(dat$event %in% c(0, 1)))
        stop("Invalid survival in ", cc)
    age <- if ("age_at_initial_pathologic_diagnosis" %in% names(cdat))
        suppressWarnings(as.numeric(as.character(cdat$age_at_initial_pathologic_diagnosis))) else rep(NA_real_, n)
    age[age < 18 | age > 105] <- NA_real_
    sex <- if ("gender" %in% names(cdat))
        clean_chr(cdat$gender) else rep(NA_character_, n)
    sex[!sex %in% c("MALE", "FEMALE")] <- NA_character_
    aa <- if ("ajcc_pathologic_tumor_stage" %in% names(cdat))
        major_stage(cdat$ajcc_pathologic_tumor_stage) else rep(NA_character_, n)
    bb <- if ("clinical_stage" %in% names(cdat))
        major_stage(cdat$clinical_stage) else rep(NA_character_, n)
    stage <- if (sum(!is.na(aa)) >= sum(!is.na(bb)))
        aa else bb
    dat$age10 <- age/10
    dat$sex <- sex
    dat$stage <- stage
    # UCEC histology is a potentially major confounder, match Stage07 grouping.
    if (cc == "UCEC" && "histological_type" %in% names(cdat)) {
        h <- clean_chr(cdat$histological_type)
        hg <- rep(NA_character_, n)
        hg[!is.na(h) & grepl("ENDOMETRIOID", h)] <- "Endometrioid"
        hg[!is.na(h) & grepl("SEROUS|MIXED", h)] <- "Serous_or_mixed"
        dat$histology <- hg
    }
    groups_table <- rbindlist(lapply(c(0L, 1L), function(g) {
        idx <- which(dat$Missing_group == g)
        data.table(Cancer = cc, Missing_group = g, N = length(idx), Events = sum(dat$event[idx]), Event_fraction = if (length(idx))
            mean(dat$event[idx]) else NA_real_, Median_observed_time = if (length(idx))
            median(dat$time[idx]) else NA_real_, Median_age = if (length(idx))
            median(age[idx], na.rm = TRUE) else NA_real_)
    }))
    clinical_rows[[cc]] <- groups_table
    # Compare covariate distributions; do not claim bias is absent if adjusted model is estimable.
    cand <- c("age10", "sex", "stage", if (cc == "UCEC") "histology")
    for (v in cand) {
        vals <- dat[[v]]
        obs <- !is.na(vals)
        cov_rows[[length(cov_rows) + 1L]] <- data.table(Cancer = cc, Covariate = v, Coverage_pct = 100 *
            mean(obs), Available_nonmissing = sum(obs), Missinggroup_prevalence_pct = if (sum(obs))
            100 * mean(dat$Missing_group[obs]) else NA_real_)
    }
    if (length(unique(dat$Missing_group)) < 2 || min(table(dat$Missing_group)) < 10) {
        km_rows[[cc]] <- data.table(Cancer = cc, Analysis = "UNADJUSTED_KM_AND_LOGRANK", N = n, Events = sum(dat$event),
            HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_, Status = "NOT_ESTIMABLE_GROUP_LT10")
        next
    }
    safe_fit <- function(formula, d) {
        w <- character()
        fit <- withCallingHandlers(tryCatch(survival::coxph(formula, data = d, ties = "efron", control = survival::coxph.control(iter.max = 50)),
            error = function(e) e), warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        })
        list(fit = fit, warning = paste(unique(w), collapse = " | "))
    }
    for (whichmodel in c("Unadjusted", "Adjusted")) {
        sub <- dat[, c("time", "event", "Missing_group"), drop = FALSE]
        cov <- character()
        if (whichmodel == "Adjusted") {
            # Covariates must be measurable in >=80% of the cohort and have >=2 observed levels.
            for (v in cand) {
                values <- dat[[v]]
                if (mean(!is.na(values)) >= 0.8 && length(unique(na.omit(values))) >= 2) {
                  sub[[v]] <- values
                  cov <- c(cov, v)
                }
            }
        }
        sub <- sub[complete.cases(sub), , drop = FALSE]
        for (v in cov) if (is.character(sub[[v]]))
            sub[[v]] <- droplevels(factor(sub[[v]]))
        if (length(unique(sub$Missing_group)) != 2 || min(table(sub$Missing_group)) < 10) {
            km_rows[[length(km_rows) + 1L]] <- data.table(Cancer = cc, Analysis = whichmodel, N = nrow(sub),
                Events = sum(sub$event), HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                Status = "SPARSE_OR_SINGLE_GROUP_AFTER_COMPLETE_CASE")
            next
        }
        fm <- as.formula(paste("survival::Surv(time,event) ~ Missing_group", if (length(cov))
            paste0(" + ", paste(cov, collapse = " + ")) else ""))
        ans <- safe_fit(fm, sub)
        if (inherits(ans$fit, "error") || !("Missing_group" %in% names(stats::coef(ans$fit)))) {
            status <- if (inherits(ans$fit, "error"))
                conditionMessage(ans$fit) else "MISSING_COEFFICIENT"
            km_rows[[length(km_rows) + 1L]] <- data.table(Cancer = cc, Analysis = whichmodel, N = nrow(sub),
                Events = sum(sub$event), HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
                Status = status)
            next
        }
        fit <- ans$fit
        b <- unname(coef(fit)["Missing_group"])
        se <- sqrt(unname(vcov(fit)["Missing_group", "Missing_group"]))
        ok <- is.finite(b) && is.finite(se) && se > 0
        km_rows[[length(km_rows) + 1L]] <- data.table(Cancer = cc, Analysis = whichmodel, N = nrow(sub),
            Events = sum(sub$event), HR = if (ok)
                exp(b) else NA_real_, CI_low = if (ok)
                exp(b - 1.96 * se) else NA_real_, CI_high = if (ok)
                exp(b + 1.96 * se) else NA_real_, P = if (ok)
                2 * pnorm(-abs(b/se)) else NA_real_, Status = if (!ok)
                "INVALID_FIT" else if (nzchar(ans$warning))
                paste("WARNING:", ans$warning) else "OK")
    }
}
clin <- rbindlist(clinical_rows)
cov <- rbindlist(cov_rows)
surv <- rbindlist(km_rows, fill = TRUE)
# Only report crude OR as descriptive, never label it a hazard ratio.
source_num <- unique(rates[grepl("^[0-9]+$", Feature_ID), Feature_ID])
map_input <- file.path(root, "Results_Rebuild", "Stage03", "Stage03_numeric_identifier_mapping_diagnostic.csv")
mapcheck <- if (file.exists(map_input)) fread(map_input) else data.table(Source_ID = source_num, Mapping_status = "NOT_AVAILABLE")
# Prevent accidental use of this diagnostic as a newly corrected cohort.
fwrite(stat, file.path(out, "Stage03B_feature_missingness_summary.csv"))
fwrite(pat, file.path(out, "Stage03B_specimen_missingness_patterns.csv"))
fwrite(set_comps, file.path(out, "Stage03B_missing_feature_set_overlap.csv"))
fwrite(signatures, file.path(out, "Stage03B_missingness_pattern_identity.csv"))
fwrite(ret, file.path(out, "Stage03B_gene_coverage_thresholds.csv"))
fwrite(clin, file.path(out, "Stage03B_clinical_group_summary.csv"))
fwrite(cov, file.path(out, "Stage03B_covariate_coverage.csv"))
fwrite(surv, file.path(out, "Stage03B_missingness_group_survival_models.csv"))
fwrite(mapcheck, file.path(out, "Stage03B_numeric_id_mapping_status.csv"))
# Figure 1 \342\200\224 exact cross-cancer gene-set membership for missing features.
member <- data.table(Feature_ID = unique(unlist(genesets, use.names = FALSE)))
for (cc in cancers) member[, (cc) := Feature_ID %in% genesets[[cc]]]
comb <- member[, .(Gene_count = .N), by = .(COAD, OV, UCEC)]
comb[, Pattern := paste0("COAD ", as.integer(COAD), "  |  OV ", as.integer(OV), "  |  UCEC ", as.integer(UCEC))]
comb <- comb[order(-Gene_count)]
p1 <- ggplot2::ggplot(comb, ggplot2::aes(x = reorder(Pattern, Gene_count), y = Gene_count, fill = Gene_count)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() + ggplot2::scale_fill_gradient(low = "#63B6C9",
    high = "#153F64", guide = "none") + ggplot2::labs(title = "Overlap of structurally unmeasured gene sets",
    subtitle = "Patterns indicate whether any specimen in a cancer has a nonfinite measurement", x = NULL,
    y = "Genes in missingness-set intersection", caption = "Pattern 1 = a feature is nonfinite in at least one specimen; 0 = no nonfinite values.") +
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 19), axis.text = ggplot2::element_text(face = "bold", colour = "black"),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), plot.caption = ggplot2::element_text(colour = "black"),
    panel.grid.major.y = ggplot2::element_blank())
ragg::agg_webp(file.path(figdir, "Stage03B_missing_gene_set_overlap.webp"), width = 3500, height = 1900,
    res = 300, quality = 95)
print(p1)
grDevices::dev.off()
# Figure 2 \342\200\224 retainable features by observation coverage; DOES NOT impute.
p2 <- ggplot2::ggplot(ret, ggplot2::aes(x = Strategy, y = Genes_retained, colour = Cancer, group = Cancer)) +
    ggplot2::geom_line(linewidth = 1.2) + ggplot2::geom_point(size = 3.8) + ggplot2::scale_y_continuous(labels = function(x) format(x,
    big.mark = ",")) + ggplot2::labs(title = "Feature availability under alternative completeness rules",
    subtitle = "Relaxing completeness preserves more genes but requires per-feature missing-data handling",
    x = NULL, y = "Features meeting minimum observed-sample fraction", colour = "Cancer", caption = "Exploratory accounting only. Original complete-feature cohorts remain unchanged.") +
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 18), axis.text = ggplot2::element_text(face = "bold", colour = "black"),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), plot.caption = ggplot2::element_text(colour = "black"))
ragg::agg_webp(file.path(figdir, "Stage03B_gene_retention_completeness_thresholds.webp"), width = 3500,
    height = 2100, res = 300, quality = 95)
print(p2)
grDevices::dev.off()
# Figure 3 \342\200\224 conditional descriptive survival association, with 95% CI and null line.
p3dat <- surv[is.finite(HR) & is.finite(CI_low) & is.finite(CI_high)]
if (nrow(p3dat)) {
    p3dat[, Label := paste(Cancer, Analysis, sep = " \302\267 ")]
    p3dat[, Label := factor(Label, levels = rev(Label))]
    p3 <- ggplot2::ggplot(p3dat, ggplot2::aes(x = HR, y = Label, colour = Analysis)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = "#777777") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low,
        xmax = CI_high), orientation = "y", height = 0.18, linewidth = 1.1) + ggplot2::geom_point(size = 3.8) +
        ggplot2::scale_x_log10() + ggplot2::labs(title = "Missingness-group association with overall survival",
        subtitle = "Cox hazard ratios for specimens with any nonfinite feature versus none", x = "Hazard ratio (log scale), with 95% confidence interval",
        y = NULL, colour = NULL, caption = "Exploratory observational association; missingness-group status is not a biological exposure. OV excluded: comparison group n=2.") +
        ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
        colour = "black", size = 18), axis.text = ggplot2::element_text(face = "bold", colour = "black"),
        axis.title = ggplot2::element_text(face = "bold", colour = "black"), plot.caption = ggplot2::element_text(colour = "black"))
    ragg::agg_webp(file.path(figdir, "Stage03B_missingness_survival_association_95CI.webp"), width = 3500,
        height = 1900, res = 300, quality = 95)
    print(p3)
    grDevices::dev.off()
}
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 03B \342\200\224 SCIENTIFIC INVESTIGATION\n====================================================================\n")
    cat("1. STRUCTURED MISSINGNESS AND RETAINED GENE COUNTS\n")
    print(stat)
    cat("\n2. EXACT OVERLAP OF UNMEASURED FEATURE SETS\n")
    print(set_comps)
    cat("\n3. IDENTICAL MISSINGNESS PATTERN WITHIN EACH GROUP?\n")
    print(signatures)
    cat("\n4. GENE RETENTION UNDER ALTERNATIVE OBSERVATION THRESHOLDS\n")
    print(ret)
    cat("\n5. MISSINGNESS GROUPS, OUTCOMES, AGE\n")
    print(clin)
    cat("\n6. COVARIATE AVAILABILITY\n")
    print(cov)
    cat("\n7. SURVIVAL ASSOCIATION OF MISSINGNESS GROUP (EXPLORATORY; 95% CI)\n")
    print(surv)
    cat("\n8. 29 NUMERIC SOURCE IDENTIFIERS: ANNOTATION STATUS\n")
    print(mapcheck)
    cat("\n9. SCIENTIFIC DECISION RULES\n")
    cat("Do not interpret missingness as proven laboratory batch without source processing metadata.\n")
    cat("Survival group associations are exploratory observational comparisons, not evidence of causality.\n")
    cat("OV has only 2 samples with complete features: group survival comparison is not estimable.\n")
    cat("90%/80% thresholds measure prospective gene availability; they do NOT imply valid imputation or retrained models.\n")
    cat("Source-scale TCGA expression is transformed and batch adjusted; not raw counts.\n")
    cat("Do not silently use numeric ID mappings, or overwrite original Stage03 saved cohorts.\n")
    cat("\n10. FIGURES\n")
    print(list.files(figdir, pattern = "\\.webp$", full.names = TRUE))
    cat("\nOUTPUT: ", out, "\n====================================================================\n",
        sep = "")
})
writeLines(report, file.path(out, "Stage03B_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Expression-scale and
# PCA diagnostics ----------------------------------------------------------------------------

options(width = 190, scipen = 5)
root <- PROJECT_ROOT
source_dir <- file.path(root, "Results_Rebuild", "Stage03")
out <- file.path(root, "Results_Rebuild", "Stage04")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (p in c("data.table", "matrixStats", "ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Required package missing: ",
    p)
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
files <- setNames(file.path(source_dir, paste0("TCGA_", cancers, "_discovery.rds")), cancers)
if (any(!file.exists(files))) stop("Missing Stage03 cohort RDS: ", paste(names(files)[!file.exists(files)],
    collapse = ", "))
missfile <- file.path(source_dir, "Stage03_sample_missingness.csv")
if (!file.exists(missfile)) stop("Missing Stage03_sample_missingness.csv")
miss <- fread(missfile)
required_miss <- c("Cancer", "Patient", "Missing_features")
if (!all(required_miss %in% names(miss))) stop("Stage03 missingness schema differs")
if (anyDuplicated(miss[, .(Cancer, Patient)])) stop("Duplicate Stage03 patient/missingness records")
check_finite <- function(x, lab) {
    if (any(!is.finite(x)))
        stop(lab, ": unexpectedly nonfinite expression")
}
all_stats <- list()
all_samp <- list()
all_pc <- list()
pca_summary <- list()
warns <- character()
set.seed(20261102)
for (cc in cancers) {
    message("Stage04: auditing ", cc, " ...")
    obj <- readRDS(files[[cc]])
    if (!is.list(obj) || !all(c("expression", "clinical") %in% names(obj)))
        stop(cc, ": invalid Stage03 RDS format")
    x <- obj$expression
    cl <- as.data.table(obj$clinical)
    if (!is.matrix(x) || !is.numeric(x) || !all(c("patient_id", "sample_id", "OS", "OS_time") %in% names(cl)))
        stop(cc, ": unexpected matrix/clinical format")
    if (ncol(x) != nrow(cl) || !identical(colnames(x), as.character(cl$sample_id)))
        stop(cc, ": sample alignment failure")
    if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)))
        stop(cc, ": duplicate identifiers")
    check_finite(x, cc)
    if (nrow(x) < 100L || ncol(x) < 30L)
        stop(cc, ": unexpectedly small input")
    n <- ncol(x)
    g <- nrow(x)
    med <- matrixStats::colMedians(x)
    q25 <- matrixStats::colQuantiles(x, probs = 0.25, drop = TRUE)
    q75 <- matrixStats::colQuantiles(x, probs = 0.75, drop = TRUE)
    iqr <- as.numeric(q75 - q25)
    neg <- colMeans(x < 0) * 100
    zero <- colMeans(x == 0) * 100
    # Robust *diagnostic* flags only. No specimen removed automatically.
    robust_z <- function(v) {
        d <- stats::mad(v, constant = 1.4826, na.rm = TRUE)
        if (!is.finite(d) || d < 1e-08)
            return(rep(NA_real_, length(v)))
        (v - stats::median(v))/d
    }
    zmed <- robust_z(med)
    ziqr <- robust_z(iqr)
    samp <- data.table(Cancer = cc, Patient = as.character(cl$patient_id), Sample = as.character(cl$sample_id),
        Expression_median = as.numeric(med), Expression_IQR = iqr, Negative_pct = neg, Exactly_zero_pct = zero,
        Robust_Z_median = zmed, Robust_Z_IQR = ziqr, Flag_median = abs(zmed) > 4, Flag_IQR = abs(ziqr) >
            4, Event = as.numeric(cl$OS), Survival_days = as.numeric(cl$OS_time))
    setkey(miss, Cancer, Patient)
    setkey(samp, Cancer, Patient)
    samp <- miss[samp, on = .(Cancer, Patient)]
    if (anyNA(samp$Missing_features) || nrow(samp) != n)
        stop(cc, ": missing Stage03 sample QC mapping")
    all_samp[[cc]] <- samp
    all_stats[[cc]] <- data.table(Cancer = cc, Patients = n, Genes = g, Deaths = sum(cl$OS == 1), Matrix_min = min(x),
        Matrix_median = stats::median(med), Matrix_max = max(x), Patient_median_Q05 = as.numeric(quantile(med,
            0.05)), Patient_median_Q95 = as.numeric(quantile(med, 0.95)), Patient_IQR_median = stats::median(iqr),
        Negative_value_pct = mean(x < 0) * 100, Zero_value_pct = mean(x == 0) * 100, Flagged_sample_median = sum(samp$Flag_median,
            na.rm = TRUE), Flagged_sample_IQR = sum(samp$Flag_IQR, na.rm = TRUE))
    # PCA is *diagnostic*: select 1000 highest-variance nonconstant genes within each cancer.
    # Center and scale each gene to unit SD; PC loadings and sample PCs are NOT survival
    # predictors.
    variances <- matrixStats::rowVars(x)
    ord <- order(variances, decreasing = TRUE)
    take <- ord[seq_len(min(1000L, length(ord)))]
    xx <- t(x[take, , drop = FALSE])
    pc <- stats::prcomp(xx, center = TRUE, scale. = TRUE, rank. = 3L)
    pvar <- 100 * pc$sdev^2/sum(pc$sdev^2)
    # IMPORTANT: rank.=3 can truncate total variance calculation; compute total standardized
    # variance exactly = number of selected scaled features, times sample variance n/(n-1) in the
    # covariance convention.
    pvar <- 100 * pc$sdev^2/ncol(xx)
    coords <- data.table(Cancer = cc, Patient = as.character(cl$patient_id), Sample = as.character(cl$sample_id),
        PC1 = as.numeric(pc$x[, 1L]), PC2 = as.numeric(pc$x[, 2L]), PC1_variance_pct = pvar[1L], PC2_variance_pct = pvar[2L])
    coords <- merge(coords, miss[, .(Cancer, Patient, Missing_features)], by = c("Cancer", "Patient"),
        all.x = TRUE, sort = FALSE)
    if (nrow(coords) != n || anyNA(coords$Missing_features))
        stop(cc, ": failed PCA metadata alignment")
    coords[, Missingness_group := ifelse(Missing_features > 0, "Structured missingness", "No missing features")]
    all_pc[[cc]] <- coords
    pc1_by_group <- NA_real_
    if (length(unique(coords$Missingness_group)) == 2L) {
        a <- coords[Missingness_group == "No missing features", PC1]
        b <- coords[Missingness_group == "Structured missingness", PC1]
        if (length(a) >= 10L && length(b) >= 10L) {
            pooled <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b))/(length(a) +
                length(b) - 2))
            if (is.finite(pooled) && pooled > 0)
                pc1_by_group <- (mean(b) - mean(a))/pooled
        }
    }
    pca_summary[[cc]] <- data.table(Cancer = cc, Genes_used = length(take), PC1_variance_pct = pvar[1L],
        PC2_variance_pct = pvar[2L], PC3_variance_pct = pvar[3L], PCA_PC1_missingness_group_standardized_mean_diff = pc1_by_group)
    if (any(abs(c(zmed, ziqr)) > 4, na.rm = TRUE))
        warns <- c(warns, paste0(cc, ": sample distribution outlier flags require later inspection; NOT excluded"))
    rm(obj, x, cl, xx, pc, coords, samp)
    invisible(gc(FALSE))
}
smry <- rbindlist(all_stats)
samples <- rbindlist(all_samp)
pcs <- rbindlist(all_pc)
pstat <- rbindlist(pca_summary)
if (sum(smry$Patients) != 4718L || sum(smry$Deaths) != 1519L) stop("Unexpected cohort drift")
if (any(!is.finite(smry$Matrix_min)) || any(smry$Genes < 17465)) stop("Invalid expression summary")
fwrite(smry, file.path(out, "Stage04_expression_scale_summary.csv"))
fwrite(samples, file.path(out, "Stage04_patient_level_expression_diagnostics.csv"))
fwrite(pstat, file.path(out, "Stage04_PCA_variance_and_missingness_association.csv"))
fwrite(pcs, file.path(out, "Stage04_PCA_scores_diagnostic.csv"))
fwrite(data.table(Statement = c("Input consists of processed/batch-adjusted pan-cancer expression, not raw counts.",
    "PCA and robust distribution flags are diagnostics; neither establishes a processing batch nor identifies mislabeled specimens.",
    "No sample exclusion, imputation, log transformation, batch correction or outcome-based adjustment performed.",
    "Cancer-wise gene z-scoring for Cox models must be TRAINING-SCOPE aware if subsequent prediction uses gene selection or fitting.",
    "Residual batch effects cannot be ruled out without independently verified batch/aliquot metadata.",
    "The original Stage03 RDS objects remain the model input; Stage04 has not overwritten or regenerated them.")),
    file.path(out, "Stage04_methodological_boundaries.csv"))
# Standalone figures \342\200\224 each with distinct diagnostic content, no multipanels.
theme_pub <- function() {
    ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
        size = 19, colour = "black"), plot.subtitle = ggplot2::element_text(colour = "black", size = 11),
        axis.title = ggplot2::element_text(face = "bold", colour = "black"), axis.text = ggplot2::element_text(face = "bold",
            colour = "black"), legend.title = ggplot2::element_text(face = "bold", colour = "black"),
        plot.caption = ggplot2::element_text(colour = "black", size = 9), panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank())
}
save_webp <- function(p, name, width = 3600, height = 2100) {
    path <- file.path(figdir, name)
    ragg::agg_webp(path, width = width, height = height, res = 300, quality = 96)
    tryCatch(print(p), finally = grDevices::dev.off())
    path
}
samples[, Cancer := factor(Cancer, levels = cancers)]
p1 <- ggplot2::ggplot(samples, ggplot2::aes(Cancer, Expression_median)) + ggplot2::geom_violin(ggplot2::aes(fill = Cancer),
    scale = "width", trim = TRUE, alpha = 0.7, colour = NA, show.legend = FALSE) + ggplot2::geom_boxplot(width = 0.15,
    fill = "white", outlier.shape = NA, linewidth = 0.45) + ggplot2::labs(title = "Sample-wise expression location varies across cancers",
    subtitle = "Within-sample median of eligible processed genes | each distribution represents patients",
    x = NULL, y = "Median processed expression (source units)", caption = "Processed, batch-adjusted expression; distributions are not raw RNA-seq library-size measures.") +
    theme_pub()
fig1 <- save_webp(p1, "Stage04_sample_expression_median_distributions.webp")
p2 <- ggplot2::ggplot(samples, ggplot2::aes(Cancer, Expression_IQR)) + ggplot2::geom_violin(ggplot2::aes(fill = Cancer),
    scale = "width", trim = TRUE, alpha = 0.7, colour = NA, show.legend = FALSE) + ggplot2::geom_boxplot(width = 0.15,
    fill = "white", outlier.shape = NA, linewidth = 0.45) + ggplot2::labs(title = "Within-sample expression dispersion across cancers",
    subtitle = "Interquartile range of eligible processed genes per patient", x = NULL, y = "Within-sample expression IQR (source units)",
    caption = "Differences may reflect biology, preprocessing or residual technical variation; do not infer batch causality.") +
    theme_pub()
fig2 <- save_webp(p2, "Stage04_sample_expression_dispersion.webp")
figs <- c(fig1, fig2)
for (cc in c("COAD", "OV", "UCEC")) {
    dd <- copy(pcs[Cancer == cc])
    if (!nrow(dd))
        stop("No PCA rows for ", cc)
    pc1 <- round(dd$PC1_variance_pct[1], 1)
    pc2 <- round(dd$PC2_variance_pct[1], 1)
    dd[, Missingness_group := factor(Missingness_group, levels = c("No missing features", "Structured missingness"))]
    pp <- ggplot2::ggplot(dd, ggplot2::aes(PC1, PC2, colour = Missingness_group)) + ggplot2::geom_point(alpha = 0.8,
        size = 2.3) + ggplot2::scale_colour_manual(values = c(`No missing features` = "#267D88", `Structured missingness` = "#D9674E"),
        drop = FALSE) + ggplot2::labs(title = paste0(cc, " | Expression PCA and source missingness"),
        subtitle = "Top 1,000 variable eligible genes; genes centred and scaled within cancer", x = sprintf("PC1 (%.1f%% variance)",
            pc1), y = sprintf("PC2 (%.1f%% variance)", pc2), colour = "Source data completeness", caption = if (cc ==
            "OV")
            "Only 2 specimens lack structured missingness; group inference is not valid." else "Missingness-group separation is descriptive; metadata are insufficient to establish technical batch origin.") +
        theme_pub() + ggplot2::theme(legend.position = "bottom", panel.grid.major.x = ggplot2::element_line(colour = "#eeeeee"))
    figs <- c(figs, save_webp(pp, paste0("Stage04_", cc, "_PCA_structured_missingness.webp"), width = 3200,
        height = 2300))
}
report <- capture.output({
    cat("=====================================================================\nPROJECT 3 | STAGE 04 \342\200\224 EXPRESSION SCIENTIFIC CHECKPOINT\n=====================================================================\n")
    cat("1. EXPRESSION DISTRIBUTIONS (TRANSFORMED/ADJUSTED; NOT RAW COUNTS)\n")
    print(smry)
    cat("\n2. PCA VARIANCE AND STRUCTURED-MISSINGNESS EFFECT-SIZE DIAGNOSTICS\n")
    print(pstat)
    cat("\n3. EXPLORATORY OUTLIER FLAGS (NO AUTOMATIC EXCLUSIONS)\n")
    print(head(samples[Flag_median == TRUE | Flag_IQR == TRUE, .(Cancer, Patient, Sample, Expression_median,
        Expression_IQR, Robust_Z_median, Robust_Z_IQR, Missing_features)][order(Cancer, -abs(Robust_Z_median))],
        40L))
    cat("\n4. ANALYTIC CORRESPONDENCE\n")
    cat("Patients: ", sum(smry$Patients), " | deaths: ", sum(smry$Deaths), " | cancers: ", nrow(smry),
        "\n", sep = "")
    cat("5. CAUTIONS\n")
    cat("PCA is not a batch-effect test without independent batch labels. Unusual PCA grouping may reflect genuine biological heterogeneity.\n")
    cat("Cancer-specific source-scale distributions are not required to be identical; no additional transformation or correction has been imposed.\n")
    cat("No imputation, gene renaming, survival modelling, or arbitrary sample exclusions.\n")
    cat("Source clinical covariates and aliquot/batch information remain incompletely characterised.\n")
    cat("Gene-wise standardisation should be within cancer; full-cohort z-scoring must not precede predictive train/test splits.\n")
    if (length(warns)) {
        cat("\n6. OUTLIER NOTES\n")
        cat(paste(unique(warns), collapse = "\n"), "\n")
    }
    cat("\n7. FIGURES WRITTEN\n")
    cat(paste(figs, collapse = "\n"), "\n")
    cat("\n8. OUTPUT\n", out, "\n=====================================================================\n")
})
writeLines(report, file.path(out, "Stage04_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
