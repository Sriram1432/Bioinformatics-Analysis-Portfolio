# Define the ten TCGA cohorts and validate patient-level survival records.  Inventories the source
# files, resolves duplicate specimens and constructs the established 4,718-patient analysis
# population without altering the original data.  Run with: Rscript 01_prepare_cohorts.R

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


# ---------------------------------------------------------------------------- Source inventory and
# study coverage ----------------------------------------------------------------------------

options(width = 180, scipen = 6)
root <- PROJECT_ROOT
data_dir <- file.path(root, "Data")
out <- file.path(root, "Results_Rebuild", "Stage01")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "ggplot2", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Install the required package: ", pkg, " (install.packages(\"", pkg, "\"))")
}
suppressPackageStartupMessages(library(data.table))
set.seed(20261102)
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")

files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
files <- files[!dir.exists(files)]
find_unique <- function(name) {
    hits <- files[basename(files) == name]
    if (length(hits) != 1L)
        stop("Expected exactly one source named ", name, "; found ", length(hits))
    hits
}
paths <- c(expression = find_unique("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"), survival = find_unique("Survival_SupplementalTable_S1_20171025_xena_sp"),
    phenotype = find_unique("TCGA_phenotype_denseDataOnlyDownload.tsv.gz"))
info <- file.info(paths)
inventory <- data.table(Source = names(paths), Filename = basename(paths), Size_MiB = round(info$size/1024^2,
    3), Modified_local = as.character(info$mtime), MD5 = unname(tools::md5sum(paths)))
# MD5 fingerprints detect file changes between stages, not source authenticity.
fwrite(inventory, file.path(out, "Stage01_source_inventory.csv"))

# Expression: HEADER ONLY + eight rows. This is not full-gene validation.
header <- names(fread(paths[["expression"]], nrows = 0L, check.names = FALSE, showProgress = FALSE))
if (!length(header) || header[[1]] != "sample") stop("Expression identifier column differs from expected \"sample\"")
exids <- header[-1L]
if (length(exids) < 100L) stop("Suspiciously small expression sample column count")
ex_preview <- fread(paths[["expression"]], nrows = 8L, select = seq_len(min(7L, length(header))), check.names = FALSE,
    showProgress = FALSE)
clinical <- fread(paths[["survival"]], showProgress = FALSE)
phenotype <- fread(paths[["phenotype"]], showProgress = FALSE)
reqclin <- c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time", "age_at_initial_pathologic_diagnosis",
    "gender", "ajcc_pathologic_tumor_stage", "clinical_stage", "histological_type")
reqphe <- c("sample", "sample_type_id", "sample_type")
if (length(setdiff(reqclin, names(clinical)))) stop("Missing clinical columns: ", paste(setdiff(reqclin,
    names(clinical)), collapse = ", "))
if (length(setdiff(reqphe, names(phenotype)))) stop("Missing phenotype columns: ", paste(setdiff(reqphe,
    names(phenotype)), collapse = ", "))

# Barcodes and identifiers: exploratory metadata integrity; no patient selection yet.
cl <- copy(clinical)
setnames(cl, c("sample", "_PATIENT", "cancer type abbreviation"), c("sample_id", "patient_id", "cancer"))
cl[, `:=`(sample_id = as.character(sample_id), patient_id = as.character(patient_id), cancer = as.character(cancer))]
ph <- copy(phenotype)
setnames(ph, "sample", "sample_id")
ph[, sample_id := as.character(sample_id)]
ph[, sample_type_id := suppressWarnings(as.integer(sample_type_id))]
cl10 <- cl[cancer %chin% cancers]
ph1 <- ph[sample_type_id == 1L & !is.na(sample_id), .(sample_id)]
primary_ids <- unique(ph1$sample_id)
expr_dup <- unique(exids[duplicated(exids) | duplicated(exids, fromLast = TRUE)])
expr_unique <- unique(exids[!exids %chin% expr_dup])
cl10[, `:=`(primary_annotated = sample_id %chin% primary_ids, expression_present = sample_id %chin% expr_unique,
    barcode_consistent = !is.na(sample_id) & !is.na(patient_id) & substr(sample_id, 1L, 12L) == patient_id)]
# Do not call these final analytical counts: Stage02 decides survival and duplicate rules.
coverage <- cl10[, .(Clinical_rows = .N, Primary_annotated = sum(primary_annotated, na.rm = TRUE), Primary_with_unique_expression = sum(primary_annotated &
    expression_present, na.rm = TRUE), Barcode_mismatches = sum(!barcode_consistent, na.rm = TRUE)),
    by = cancer]
coverage[, cancer := factor(cancer, levels = cancers)]
setorder(coverage, cancer)
coverage[, cancer := as.character(cancer)]
fwrite(coverage, file.path(out, "Stage01_sample_coverage_diagnostic.csv"))

issues <- data.table(Check = c("Expression header duplicates", "Clinical rows: ten cancers", "Phenotype duplicate sample rows",
    "Clinical duplicate sample rows (ten cancers)", "Clinical barcode mismatches (ten cancers)", "Expression IDs matching any clinical sample (ten cancers)",
    "Expression header sample columns"), Value = c(length(expr_dup), nrow(cl10), sum(duplicated(ph$sample_id)),
    sum(duplicated(cl10$sample_id)), sum(!cl10$barcode_consistent), sum(unique(cl10$sample_id) %chin%
        expr_unique), length(exids)))
fwrite(issues, file.path(out, "Stage01_identifier_diagnostics.csv"))

# One scientifically informative, standalone inventory figure (not a biological result).
plot_data <- melt(coverage[, .(cancer, Clinical_rows, Primary_annotated, Primary_with_unique_expression)],
    id.vars = "cancer", variable.name = "Availability", value.name = "Records")
plot_data[, Availability := factor(Availability, levels = c("Clinical_rows", "Primary_annotated",
    "Primary_with_unique_expression"), labels = c("Clinical records", "Primary tissue annotation", "Primary + unique RNA profile"))]
plot_data[, cancer := factor(cancer, levels = cancers)]
p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cancer, y = Records, fill = Availability)) + ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.82),
    width = 0.72) + ggplot2::scale_fill_manual(values = c("#778899", "#368AA0", "#1E4B71")) + ggplot2::labs(title = "TCGA input coverage across ten cancer types",
    subtitle = "Diagnostic sample availability; not a final survival-analysis cohort", x = "TCGA cancer type",
    y = "Number of clinical records", fill = NULL) + ggplot2::theme_classic(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 17), plot.subtitle = ggplot2::element_text(colour = "black", size = 11),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), axis.text = ggplot2::element_text(colour = "black"),
    legend.text = ggplot2::element_text(colour = "black"), legend.position = "bottom")
figpath <- file.path(figdir, "Stage01_TCGA_input_sample_coverage.webp")
ggplot2::ggsave(figpath, p, device = ragg::agg_webp, width = 10, height = 6, units = "in", res = 300)
if (!file.exists(figpath) || file.info(figpath)$size < 5000) stop("WebP figure was not saved correctly")

# The original Xena dataset's preprocessing cannot be established from these filenames or first
# eight rows alone. Provenance remains a required checkpoint.
flags <- data.table(Topic = c("Clinical input schema", "Expression header", "Phenotype schema", "Expression preprocessing/provenance",
    "Historical/mixed gene identifier mapping", "Residual batch effects", "Patient eligibility / survival status"),
    Status = c("CHECKED", "CHECKED", "CHECKED", "NOT YET VERIFIED", "NOT YET VERIFIED", "NOT YET ASSESSED",
        "DEFER TO STAGE 02"), Meaning = c("Required clinical variables present", "Sample columns identified",
        "Sample classification fields available", "Source transformations and original build require source documentation",
        "Previewed feature IDs do not establish modern gene-symbol mapping", "Requires documented batches / full expression assessment",
        "Stage 01 deliberately does not define the analytic cohort"))
fwrite(flags, file.path(out, "Stage01_scientific_checkpoint.csv"))

report <- capture.output({
    cat("====================================================================\n")
    cat("PROJECT 3 | STAGE 01 \342\200\224 REBUILD SCIENTIFIC CHECKPOINT\n")
    cat("====================================================================\n")
    cat("1. SOURCE FINGERPRINTS\n")
    print(inventory)
    cat("\n2. SOURCE SCHEMAS\n")
    print(data.table(Source = c("expression", "survival", "phenotype"), Rows_loaded = c(8L, nrow(clinical),
        nrow(phenotype)), Columns = c(length(header), ncol(clinical), ncol(phenotype))))
    cat("\n3. EXPRESSION IDENTIFIER PREVIEW (NOT A COMPLETE ANNOTATION)\n")
    print(ex_preview[, 1, with = FALSE])
    cat("\n4. CANCER INPUT COVERAGE \342\200\224 NOT FINAL COHORT\n")
    print(coverage)
    cat("\n5. IDENTIFIER FLAGS\n")
    print(issues)
    if (length(expr_dup))
        cat("Duplicate expression sample names: ", paste(head(expr_dup, 30L), collapse = ", "), "\n")
    cat("\n6. SCIENTIFIC CHECKPOINTS\n")
    print(flags)
    cat("\n7. FIGURE CREATED\n", figpath, "\n")
    cat("\n8. INTERPRETATION\n")
    cat("Counts are diagnostic availability, not matched survival-eligible patients.\n")
    cat("No genome-wide expression import, no survival endpoint filtering, no models.\n")
    cat("MD5 verifies file identity between runs; it does not authenticate dataset origin.\n")
    cat("Stage 01 remains REVIEW until expression provenance and identifier semantics are checked.\n")
    cat("Original data and existing Results/ outputs were NOT modified.\n")
    cat("\nOutput: ", out, "\n", sep = "")
    cat("====================================================================\n")
})
writeLines(report, file.path(out, "Stage01_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Identifier and
# provenance checks ----------------------------------------------------------------------------

options(width = 180, scipen = 5)
root <- PROJECT_ROOT
out <- file.path(root, "Results_Rebuild", "Stage01B")
fig <- file.path(out, "Figures")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
for (p in c("data.table", "ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ",
    p)
suppressPackageStartupMessages(library(data.table))
data_dir <- file.path(root, "Data")
find_one <- function(nm) {
    p <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
    p <- p[basename(p) == nm]
    if (length(p) != 1L)
        stop("Expected one file ", nm)
    p
}
expfile <- find_one("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
clfile <- find_one("Survival_SupplementalTable_S1_20171025_xena_sp")
phfile <- find_one("TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
headings <- names(fread(expfile, nrows = 0L, check.names = FALSE, showProgress = FALSE))
ids <- headings[-1L]
dupids <- unique(ids[duplicated(ids)])
# Read just source identifier column, all gene rows. Treat as string: do not let fread coerce
# Entrez-like numeric identifiers to floating-point.
gene_column <- fread(expfile, select = 1L, colClasses = "character", showProgress = FALSE)[[1L]]
feature <- trimws(as.character(gene_column))
feature_info <- data.table(Feature_ID = feature, Category = fifelse(is.na(feature) | !nzchar(feature),
    "Blank", fifelse(grepl("^[0-9]+$", feature), "Numeric/Entrez-like", fifelse(grepl("^ENS[A-Z]*G[0-9]+(\\.[0-9]+)?$",
        feature), "Ensembl-like", fifelse(grepl("^[[:alpha:]][[:alnum:]._-]*$", feature), "Symbol-like",
        "Other")))))
feature_counts <- feature_info[, .(N = .N), by = Category][order(-N)]
dupgene <- feature_info[, .N, by = Feature_ID][N > 1L][order(-N)]
fwrite(feature_counts, file.path(out, "Stage01B_identifier_types.csv"))
fwrite(dupgene, file.path(out, "Stage01B_duplicate_feature_identifiers.csv"))
# Inspect duplicate sample columns as *independent physical columns*. If source columns differ,
# choosing one arbitrarily would be scientifically unsafe.
pairs <- list()
discrepancies <- list()
for (s in dupids) {
    pos <- which(headings == s)  # positions in full file, including gene ID field
    dt <- fread(expfile, select = c(1L, pos), check.names = FALSE, showProgress = FALSE)
    yy <- as.matrix(dt[, -1L, with = FALSE])
    storage.mode(yy) <- "double"
    if (ncol(yy) != length(pos))
        stop("Duplicate selection altered columns for ", s)
    ref <- yy[, 1L]
    for (j in 2:ncol(yy)) {
        other <- yy[, j]
        equal_na <- (is.na(ref) & is.na(other)) | (!is.na(ref) & !is.na(other) & ref == other)
        dis <- sum(!equal_na)
        pd <- data.table(Sample_ID = s, Source_column_first = pos[1L], Source_column_other = pos[j],
            N_features = length(ref), N_discordant = dis, Discordance_pct = round(100 * dis/length(ref),
                5), Max_absolute_difference = if (any(is.finite(ref) & is.finite(other)))
                max(abs(ref[is.finite(ref) & is.finite(other)] - other[is.finite(ref) & is.finite(other)])) else NA_real_)
        pairs[[length(pairs) + 1L]] <- pd
        if (dis)
            discrepancies[[length(discrepancies) + 1L]] <- data.table(Sample_ID = s, First_10_discordant_feature_ids = paste(head(feature[!equal_na],
                10), collapse = ";"))
    }
}
dup_results <- if (length(pairs)) rbindlist(pairs) else data.table(Sample_ID = character(), N_discordant = integer())
fwrite(dup_results, file.path(out, "Stage01B_duplicate_expression_column_comparison.csv"))
if (length(discrepancies)) fwrite(rbindlist(discrepancies), file.path(out, "Stage01B_duplicate_column_discrepancy_examples.csv"))
# Cohort-level metadata overlap; do not select analytical patients here.
cl <- fread(clfile, showProgress = FALSE)
ph <- fread(phfile, showProgress = FALSE)
stopifnot(all(c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time") %in% names(cl)), all(c("sample",
    "sample_type_id") %in% names(ph)))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
sub <- cl[get("cancer type abbreviation") %chin% cancers]
clids <- unique(as.character(sub$sample))
phids <- unique(as.character(ph[sample_type_id == 1L]$sample))
exp_unique <- setdiff(unique(ids), dupids)
# Explicitly flag duplicates with primary sample type and whether they occur in ten selected cancer
# clinical records.
dup_meta <- data.table(Sample_ID = dupids, In_10_cancer_clinical = dupids %in% clids, Phenotype_primary = dupids %in%
    phids, N_expression_columns = vapply(dupids, function(s) sum(ids == s), integer(1)))
fwrite(dup_meta, file.path(out, "Stage01B_duplicate_sample_metadata.csv"))
# Check phenotype sample types for mismatching TCGA vial positions.
ph_check <- ph[!is.na(sample), .(sample = as.character(sample), sample_type_id)]
ph_check[, Barcode_type := suppressWarnings(as.integer(substr(sample, 14L, 15L)))]
ph_check[, Expected_primary := Barcode_type == 1L]
ph_mismatch <- ph_check[!is.na(Barcode_type) & !is.na(sample_type_id) & ((sample_type_id == 1L) != Expected_primary)]
fwrite(ph_mismatch, file.path(out, "Stage01B_primary_barcode_phenotype_discordance.csv"))
# Figure: informative identifier composition, not decorative biology.
p <- ggplot2::ggplot(feature_counts, ggplot2::aes(x = reorder(Category, N), y = N)) + ggplot2::geom_col(fill = "#22618A",
    width = 0.62) + ggplot2::geom_text(ggplot2::aes(label = format(N, big.mark = ",")), hjust = -0.08,
    fontface = "bold", size = 4) + ggplot2::coord_flip(clip = "off") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
    0.2))) + ggplot2::labs(title = "Gene identifier composition in the TCGA expression source", subtitle = "Classification by identifier syntax; database mappings require authoritative annotation",
    x = NULL, y = "Number of source feature rows") + ggplot2::theme_classic(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 16), axis.title = ggplot2::element_text(face = "bold", colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"))
plot_path <- file.path(fig, "Stage01B_feature_identifier_classes.webp")
ggplot2::ggsave(plot_path, p, device = ragg::agg_webp, width = 11, height = 5.5, units = "in", res = 300)
# Provenance: source URL is cited but source transformation metadata is not independently
# authenticated by an MD5 or file name. Do not infer from values.
source_note <- c("SOURCE: UCSC Xena TCGA PanCanAtlas, EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena",
    "UCSC Xena FAQ: TCGA PanCan Atlas expression was generated by the PanCan Atlas project and normalised for batch effects.",
    "https://ucsc-xena.gitbook.io/project/faq/advanced-data-and-datasets", "Dataset landing page: https://xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net",
    "CAUTION: Specific normalisation units, algorithm/version, probe map, and processing batch labels require original dataset metadata; not verified by this local analysis.",
    "Historical feature IDs must be mapped against an explicitly versioned gene annotation before enrichment/protein integration.",
    "Do not impute structured missing measurements or treat transformed expression as raw RNA-seq counts.")
writeLines(source_note, file.path(out, "Stage01B_source_provenance_notes.txt"))
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 01B \342\200\224 SCIENTIFIC RESOLUTION REPORT\n====================================================================\n")
    cat("1. ALL SOURCE FEATURE ID CLASSES\n")
    print(feature_counts)
    cat("\n2. DUPLICATED FEATURE IDENTIFIERS\n")
    print(dupgene)
    cat("\n3. DUPLICATED EXPRESSION SAMPLE COMPARISONS (ACTUAL VALUES)\n")
    print(dup_results)
    cat("\n4. DUPLICATE SAMPLES: COHORT RELEVANCE\n")
    print(dup_meta)
    cat("\n5. PHENOTYPE VS BARCODE SAMPLE-TYPE VALIDATION\n")
    cat("Evaluated rows: ", nrow(ph_check), " | Discordant rows: ", nrow(ph_mismatch), "\n", sep = "")
    if (nrow(ph_mismatch))
        print(head(ph_mismatch, 20))
    cat("\n6. SCIENTIFIC DECISION CRITERIA\n")
    cat("Duplicate expression columns with differing values -> investigate before automatic exclusion or retention.\n")
    cat("Duplicated feature symbols -> exclude ambiguous rows until verified mapping.\n")
    cat("Numeric/Entrez-like IDs are not inherently invalid; use versioned mapping in Stage03.\n")
    cat("Provenance source is documented, but detailed transformation and probe mapping remain external metadata checks.\n")
    cat("A barcode type mismatch requires specimen-level review.\n")
    cat("\n7. FIGURE\n")
    cat(plot_path, "\n")
    cat("\n8. REPORT DIRECTORY\n")
    cat(out, "\n====================================================================\n")
})
writeLines(report, file.path(out, "Stage01B_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Duplicate-specimen
# resolution ----------------------------------------------------------------------------

options(width = 175, scipen = 5)
root <- PROJECT_ROOT
out <- file.path(root, "Results_Rebuild", "Stage01C")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (p in c("data.table", "ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ",
    p)
suppressPackageStartupMessages(library(data.table))
find_one <- function(filename) {
    f <- list.files(file.path(root, "Data"), recursive = TRUE, full.names = TRUE)
    h <- f[basename(f) == filename]
    if (length(h) != 1L)
        stop("Expected exactly one ", filename)
    h
}
xfile <- find_one("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
cfile <- find_one("Survival_SupplementalTable_S1_20171025_xena_sp")
pfile <- find_one("TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
headers <- names(fread(xfile, nrows = 0L, check.names = FALSE, showProgress = FALSE))
duplicated_ids <- unique(headers[-1L][duplicated(headers[-1L])])
clinical <- fread(cfile, showProgress = FALSE)
phenotype <- fread(pfile, showProgress = FALSE)
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
ct <- clinical[get("cancer type abbreviation") %chin% cancers]
rows <- list()
details <- list()
for (s in duplicated_ids) {
    pos <- which(headers == s)
    # fread select by positional indices, preserving independently measured columns
    dat <- fread(xfile, select = c(1L, pos), check.names = FALSE, showProgress = FALSE)
    if (ncol(dat) != length(pos) + 1L)
        stop("Column selection failed: ", s)
    a <- as.numeric(dat[[2L]])
    b <- as.numeric(dat[[3L]])
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 100L)
        stop("Too few comparable features for ", s)
    pearson <- suppressWarnings(cor(a[ok], b[ok], method = "pearson"))
    spearman <- suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
    dd <- a[ok] - b[ok]
    ph <- phenotype[sample == s]
    cl <- ct[sample == s]
    if (nrow(ph) > 1L || nrow(cl) > 1L)
        warning("Multiple metadata rows: ", s)
    ph_type <- if (nrow(ph) == 1L)
        as.character(ph$sample_type_id) else NA_character_
    cl_cancer <- if (nrow(cl) == 1L)
        as.character(cl[["cancer type abbreviation"]]) else NA_character_
    cl_os <- if (nrow(cl) == 1L)
        as.character(cl$OS) else NA_character_
    cl_time <- if (nrow(cl) == 1L)
        as.character(cl[["OS.time"]]) else NA_character_
    barcode <- suppressWarnings(as.integer(substr(s, 14L, 15L)))
    source_type <- if (barcode == 1L)
        "01 primary tumour" else if (barcode == 2L)
        "02 recurrent tumour" else paste("Type", barcode)
    rows[[length(rows) + 1L]] <- data.table(Sample_ID = s, Cancer = cl_cancer, In_selected_cancers = nrow(cl) ==
        1L, Source_position_1 = pos[1L], Source_position_2 = pos[2L], Barcode_type = barcode, Barcode_interpretation = source_type,
        Phenotype_type = ph_type, Phenotype_disagrees_barcode = !is.na(ph_type) && is.finite(barcode) &&
            as.integer(ph_type) != barcode, Event = cl_os, Survival_days = cl_time, Finite_pairs = sum(ok),
        Pearson_r = pearson, Spearman_r = spearman, Median_abs_difference = median(abs(dd)), P95_abs_difference = as.numeric(quantile(abs(dd),
            0.95)), Pct_abs_difference_gt_1 = 100 * mean(abs(dd) > 1), Pct_exact_equal = 100 * mean(dd ==
            0), Scientific_action = if (nrow(cl) == 1L && nrow(ph) == 1L && as.integer(ph_type) == 1L)
            "PRIMARY: retain exclusion until aliquot identity and provenance are resolved" else if (nrow(cl) == 1L && nrow(ph) == 1L && as.integer(ph_type) != barcode)
            "BARCODE/PHENOTYPE MISMATCH: quarantine, confirm specimen classification" else "Not eligible as an unambiguous primary tumour in selected cancer cohort")
    details[[length(details) + 1L]] <- data.table(Sample_ID = s, Clinical_rows = nrow(cl), Phenotype_rows = nrow(ph),
        Clinical_patient = if (nrow(cl) == 1L)
            as.character(cl[["_PATIENT"]]) else NA_character_, Phenotype_sample_type = if (nrow(ph) == 1L)
            as.character(ph$sample_type) else NA_character_)
}
res <- rbindlist(rows)
extra <- rbindlist(details)
fwrite(res, file.path(out, "Stage01C_duplicate_expression_similarity_and_metadata.csv"))
fwrite(extra, file.path(out, "Stage01C_duplicate_metadata_detail.csv"))
# Crucial: quantify impact of duplicates on eligible Stage03 cohort, not the entire TCGA source.
primary_impacted <- res[In_selected_cancers == TRUE & Phenotype_type == "1"]
needs_resolution <- res[In_selected_cancers == TRUE & (Phenotype_disagrees_barcode | Phenotype_type ==
    "1")]
# A real data-quality figure: systematic similarity of independent source columns.
p <- ggplot2::ggplot(res, ggplot2::aes(x = reorder(Sample_ID, Pearson_r), y = Pearson_r)) + ggplot2::geom_segment(ggplot2::aes(xend = Sample_ID,
    y = 0, yend = Pearson_r), linewidth = 0.85, colour = "#A2B9C9") + ggplot2::geom_point(ggplot2::aes(fill = In_selected_cancers),
    shape = 21, size = 4.7, stroke = 0.75, colour = "black") + ggplot2::scale_fill_manual(values = c(`FALSE` = "#B9C3CC",
    `TRUE` = "#087E8B"), labels = c(`FALSE` = "Outside selected cancer cohorts", `TRUE` = "In selected cancer cohorts")) +
    ggplot2::coord_flip(clip = "off") + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0,
    1, 0.2)) + ggplot2::labs(title = "Technical similarity of duplicated TCGA sample columns", subtitle = "Pearson correlation across finite gene-level measurements; similarity does not establish specimen equivalence",
    x = NULL, y = "Gene-expression Pearson correlation", fill = NULL) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "black", size = 15), axis.title = ggplot2::element_text(face = "bold",
        colour = "black"), axis.text = ggplot2::element_text(colour = "black"), legend.position = "bottom",
        panel.grid.major.y = ggplot2::element_blank())
figure <- file.path(figdir, "Stage01C_duplicate_sample_expression_similarity.webp")
ggplot2::ggsave(figure, p, device = ragg::agg_webp, width = 12, height = 6.5, units = "in", res = 300,
    bg = "white")
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 01C \342\200\224 SPECIMEN RESOLUTION CHECKPOINT\n====================================================================\n")
    cat("1. DUPLICATE SOURCE COLUMNS: ACTUAL EXPRESSION SIMILARITY\n")
    print(res[, .(Sample_ID, Cancer, Barcode_type, Phenotype_type, Pearson_r, Spearman_r, Median_abs_difference,
        P95_abs_difference, Pct_abs_difference_gt_1, Pct_exact_equal)])
    cat("\n2. METADATA AND COHORT RELEVANCE\n")
    print(res[, .(Sample_ID, In_selected_cancers, Barcode_interpretation, Phenotype_type, Phenotype_disagrees_barcode,
        Event, Survival_days, Scientific_action)], nrows = nrow(res))
    cat("\n3. ADDITIONAL METADATA\n")
    print(extra)
    cat("\n4. PRIMARY SAMPLES EXCLUDED BY DUPLICATE-HEADER POLICY\n")
    print(primary_impacted[, .(Sample_ID, Cancer, Event, Survival_days)])
    cat("\n5. OUTSTANDING IDENTITY/TYPE ANOMALIES\n")
    print(needs_resolution[, .(Sample_ID, Cancer, Barcode_type, Phenotype_type, Pearson_r)])
    cat("\n6. SCIENTIFIC DECISION\n")
    cat("Different expression columns under the SAME truncated TCGA sample identifier cannot be safely merged or automatically selected.\n")
    cat("Correlation identifies similarity but cannot verify aliquot, sample identity, laboratory batch or technical replication.\n")
    cat("Without externally authenticated aliquot-level metadata, continue conservative exclusion of ambiguous primary profiles.\n")
    cat("A phenotype/barcode disagreement is quarantined until source annotation is reconciled.\n")
    cat("These results do not alter original Stage01\342\200\22309 datasets or establish downstream biological conclusions.\n")
    cat("\nFIGURE: ", figure, "\nOUTPUT: ", out, "\n====================================================================\n")
})
writeLines(report, file.path(out, "Stage01C_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Patient cohort and
# survival quality control
# ----------------------------------------------------------------------------

options(width = 180, scipen = 6)
root <- PROJECT_ROOT
indir <- file.path(root, "Data")
out <- file.path(root, "Results_Rebuild", "Stage02")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (p in c("data.table", "ggplot2", "survival", "ragg")) {
    if (!requireNamespace(p, quietly = TRUE))
        stop("Please install missing package: ", p)
}
suppressPackageStartupMessages(library(data.table))
ccs <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
files <- list.files(indir, full.names = TRUE, recursive = TRUE)
getfile <- function(nm) {
    hits <- files[basename(files) == nm]
    if (length(hits) != 1L)
        stop("Expected one input named ", nm, "; found ", length(hits))
    hits
}
clinical_file <- getfile("Survival_SupplementalTable_S1_20171025_xena_sp")
phe_file <- getfile("TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
exp_file <- getfile("EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz")
cli <- data.table::fread(clinical_file, showProgress = FALSE)
phe <- data.table::fread(phe_file, showProgress = FALSE)
check <- function(d, cols, nm) if (length(setdiff(cols, names(d)))) stop(nm, " missing: ", paste(setdiff(cols,
    names(d)), collapse = ", "))
check(cli, c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time"), "Clinical")
check(phe, c("sample", "sample_type_id"), "Phenotype")
data.table::setnames(cli, c("sample", "_PATIENT", "cancer type abbreviation", "OS", "OS.time"), c("sample_id",
    "patient_id", "cancer", "OS", "time"))
data.table::setnames(phe, "sample", "sample_id")
cli[, `:=`(sample_id = as.character(sample_id), patient_id = as.character(patient_id), cancer = as.character(cancer),
    OS = suppressWarnings(as.numeric(as.character(OS))), time = suppressWarnings(as.numeric(as.character(time))))]
phe[, sample_id := as.character(sample_id)]
header <- names(data.table::fread(exp_file, nrows = 0, check.names = FALSE, showProgress = FALSE))
if (!length(header) || header[1] != "sample") stop("Unexpected expression first column")
exp_cols <- header[-1]
ambig <- unique(exp_cols[duplicated(exp_cols) | duplicated(exp_cols, fromLast = TRUE)])
# Restrict interpretation to phenotype type 01, not the barcode suffix alone.
primary_ids <- unique(phe[!is.na(sample_type_id) & sample_type_id == 1L, sample_id])
selected <- cli[cancer %chin% ccs]
primary <- selected[sample_id %chin% primary_ids]
primary[, `:=`(barcode_patient = substr(sample_id, 1, 12), barcode_type = suppressWarnings(as.integer(substr(sample_id,
    14, 15))), valid_endpoint = !is.na(OS) & OS %in% c(0, 1) & is.finite(time) & time >= 0, positive_time = is.finite(time) &
    time > 0)]
# Critical: choose no ambiguous profiles; report exclusions rather than guessing.
primary[, barcode_valid := !is.na(sample_id) & !is.na(patient_id) & patient_id == barcode_patient]
primary[, expression_present := sample_id %chin% exp_cols]
primary[, ambiguous_expression := sample_id %chin% ambig]
# Repeated patient records are excluded conservatively, before filtering time.
dup_pat <- primary[, .N, by = .(cancer, patient_id)][N > 1L]
primary[, duplicate_patient := patient_id %chin% dup_pat$patient_id]
primary[, eligible_pre_time := valid_endpoint & barcode_valid & expression_present & !ambiguous_expression &
    !duplicate_patient]
eligible <- primary[eligible_pre_time == TRUE & positive_time == TRUE]
if (anyDuplicated(eligible$patient_id) || anyDuplicated(eligible$sample_id)) stop("Eligible identifiers are not unique")
if (!all(eligible$barcode_type == 1L)) stop("Unexpected non-primary barcode in phenotype primary cohort")
# Per-cancer characteristics before and after endpoint and expression exclusions.
cohort <- selected[, .(Clinical = .N), by = cancer]
cohort <- merge(cohort, primary[, .(Phenotype_primary = .N, Expression_present = sum(expression_present),
    Ambiguous_profiles = sum(ambiguous_expression), Barcode_mismatch = sum(!barcode_valid), Invalid_endpoint = sum(!valid_endpoint),
    Zero_day = sum(valid_endpoint & time == 0)), by = cancer], by = "cancer", all.x = TRUE)
cohort <- merge(cohort, eligible[, .(Eligible = .N, Deaths = sum(OS == 1), Censored = sum(OS == 0), Followup_median_days = median(time),
    Deaths_under_30_days = sum(OS == 1 & time < 30)), by = cancer], by = "cancer", all.x = TRUE)
for (nm in setdiff(names(cohort), "cancer")) data.table::set(cohort, which(is.na(cohort[[nm]])), nm,
    0)
cohort[, cancer := factor(cancer, levels = ccs)]
data.table::setorder(cohort, cancer)
# Detailed flow with explicit distinct units (row counts, not uniquely matched patients).
flow <- data.table::data.table(Order = 1:7, Stage = c("Clinical records: all cancers", "Clinical records: selected 10 cancers",
    "Phenotype-annotated primary records", "Valid endpoint + matched identifiers + unique expression",
    "Positive survival duration (final eligible)", "Deaths in eligible cohort", "Censored in eligible cohort"),
    N = c(nrow(cli), nrow(selected), nrow(primary), sum(primary$eligible_pre_time), nrow(eligible), sum(eligible$OS ==
        1), sum(eligible$OS == 0)))
# Diagnostic exclusions, potentially overlapping; never sum to infer a mutually exclusive flow.
exclusions <- data.table::data.table(Check = c("Invalid endpoint (NA/event coding/negative time)", "Zero-day endpoint among valid",
    "Missing expression column", "Duplicate expression sample name", "Clinical/patient barcode mismatch",
    "Repeated clinical patient"), N = c(sum(!primary$valid_endpoint), sum(primary$valid_endpoint & primary$time ==
    0), sum(!primary$expression_present), sum(primary$ambiguous_expression), sum(!primary$barcode_valid),
    sum(primary$duplicate_patient)))
zero <- primary[valid_endpoint == TRUE & time == 0, .(Cancer = cancer, Sample_ID = sample_id, Patient_ID = patient_id,
    Event = OS)]
# Short follow-up and observed survival describe the cohort, not model performance.
fit_km <- survival::survfit(survival::Surv(time, OS) ~ cancer, data = as.data.frame(eligible), conf.type = "log-log")
ss <- summary(fit_km)
km <- data.table::data.table(Stratum = as.character(ss$strata), Time_years = ss$time/365.25, S = ss$surv,
    Low = ss$lower, High = ss$upper, At_risk = ss$n.risk, Deaths = ss$n.event)
km[, Cancer := sub("^cancer=", "", Stratum)]
# Add S(0)=1 independently per stratum; no interpolation and no artificial CI.
km <- data.table::rbindlist(list(data.table::data.table(Cancer = ccs, Time_years = 0, S = 1, Low = 1,
    High = 1, At_risk = as.integer(cohort$Eligible), Deaths = 0L, Stratum = paste0("cancer=", ccs)),
    km), fill = TRUE)
data.table::setorder(km, Cancer, Time_years)
# Optional survival landmarks; explicitly NA if beyond maximum observed event/censor time.
landmarks <- data.table::rbindlist(lapply(ccs, function(cc) {
    d <- eligible[cancer == cc]
    tmax <- max(d$time)
    f <- survival::survfit(survival::Surv(time, OS) ~ 1, data = as.data.frame(d), conf.type = "log-log")
    vals <- lapply(c(365.25, 3 * 365.25, 5 * 365.25), function(t) {
        if (t > tmax)
            return(c(S = NA_real_, Low = NA_real_, High = NA_real_, At_risk = 0))
        x <- summary(f, times = t, extend = FALSE)
        c(S = x$surv, Low = x$lower, High = x$upper, At_risk = x$n.risk)
    })
    data.table::data.table(Cancer = cc, Horizon_years = c(1, 3, 5), S = vapply(vals, `[[`, numeric(1),
        "S"), CI_low = vapply(vals, `[[`, numeric(1), "Low"), CI_high = vapply(vals, `[[`, numeric(1),
        "High"), N_at_risk = as.integer(vapply(vals, `[[`, numeric(1), "At_risk")))
}))
# Save all tables first.
data.table::fwrite(cohort, file.path(out, "Stage02_cohort_by_cancer.csv"))
data.table::fwrite(flow, file.path(out, "Stage02_cohort_flow.csv"))
data.table::fwrite(exclusions, file.path(out, "Stage02_nonexclusive_exclusion_flags.csv"))
data.table::fwrite(zero, file.path(out, "Stage02_zero_day_records.csv"))
data.table::fwrite(km, file.path(out, "Stage02_KM_curves_descriptive.csv"))
data.table::fwrite(landmarks, file.path(out, "Stage02_KM_landmarks.csv"))
# Journal-style graphics: one standalone graphic per scientific question.
ink <- "#111827"
muted <- "#526175"
teal <- "#0E7490"
rust <- "#CF563C"
blue <- "#2563A6"
theme_pub <- function(base = 13) ggplot2::theme_minimal(base_size = base, base_family = "sans") + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = ink, size = base + 4, margin = ggplot2::margin(b = 8)), plot.subtitle = ggplot2::element_text(colour = muted,
    size = base - 1, margin = ggplot2::margin(b = 15)), axis.title = ggplot2::element_text(face = "bold",
    colour = ink), axis.text = ggplot2::element_text(colour = ink), strip.text = ggplot2::element_text(face = "bold",
    colour = ink), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(face = "bold", colour = ink), plot.margin = ggplot2::margin(22,
        24, 22, 22))
saveplot <- function(plot, name, w = 12, h = 7) {
    fp <- file.path(figdir, paste0(name, ".webp"))
    ggplot2::ggsave(fp, plot, width = w, height = h, units = "in", dpi = 320, device = ragg::agg_webp,
        bg = "white")
    fp
}
# Fig 1: sorted dot-and-range visual; primary and final count endpoints avoid a basic bar chart.
d1 <- data.table::melt(as.data.table(cohort)[, .(cancer, Phenotype_primary, Eligible)], id.vars = "cancer",
    measure.vars = c("Phenotype_primary", "Eligible"), variable.name = "Cohort", value.name = "Patients")
d1[, Cohort := factor(Cohort, levels = c("Phenotype_primary", "Eligible"), labels = c("Primary annotation",
    "Final eligible"))]
p1 <- ggplot2::ggplot(d1, ggplot2::aes(x = Patients, y = cancer, group = cancer)) + ggplot2::geom_line(colour = "#B6C5D3",
    linewidth = 1.7) + ggplot2::geom_point(ggplot2::aes(colour = Cohort), size = 4.2) + ggplot2::scale_colour_manual(values = c(`Primary annotation` = "#9BAFC1",
    `Final eligible` = teal)) + ggplot2::labs(title = "From primary-tumour annotation to analytic eligibility",
    subtitle = "Diagnostic attrition across ten TCGA cancers \302\267 one primary-tumour profile per patient",
    x = "Number of specimens / patients", y = NULL, colour = NULL) + theme_pub(15) + ggplot2::theme(legend.position = "top")
f1 <- saveplot(p1, "Stage02_cohort_attrition_dotplot", 13, 7)
# Fig 2: event fractions with per-cancer N labelled.
d2 <- copy(cohort)
d2[, Event_pct := 100 * Deaths/Eligible]
p2 <- ggplot2::ggplot(d2, ggplot2::aes(x = Event_pct, y = reorder(as.character(cancer), Event_pct))) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Event_pct, yend = reorder(as.character(cancer),
        Event_pct)), colour = "#D5E2E7", linewidth = 2.5) + ggplot2::geom_point(colour = rust, size = 5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%   (n=%d)", Event_pct, Eligible)), hjust = -0.12,
        size = 4, fontface = "bold", colour = ink) + ggplot2::scale_x_continuous(limits = c(0, 100),
    breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = c(0.01, 0.03))) + ggplot2::labs(title = "Observed mortality differs substantially across cancer cohorts",
    subtitle = "Event fraction is not a mortality rate: follow-up durations and censoring differ", x = "Patients with observed death (%)",
    y = NULL) + theme_pub(14)
f2 <- saveplot(p2, "Stage02_observed_event_fraction", 12, 7)
# Fig 3: ten individually faceted descriptive KM curves, with valid 95% CI ribbons.  Cap plotting
# horizon at 5y for readability; not an estimate of 5y survival in sparse tails.
kplot <- km[Time_years <= 5 & is.finite(S) & is.finite(Low) & is.finite(High)]
kplot[, Cancer := factor(Cancer, levels = ccs)]
p3 <- ggplot2::ggplot(kplot, ggplot2::aes(Time_years, S, group = Cancer)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = Low,
    ymax = High), fill = teal, alpha = 0.18, colour = NA) + ggplot2::geom_step(colour = teal, linewidth = 0.9) +
    ggplot2::facet_wrap(~Cancer, ncol = 5) + ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 *
    x), "%")) + ggplot2::coord_cartesian(xlim = c(0, 5)) + ggplot2::labs(title = "Overall survival across the ten analytic TCGA cohorts",
    subtitle = "Descriptive Kaplan\342\200\223Meier estimates with 95% confidence intervals \302\267 no cross-cancer causal comparison",
    x = "Time since baseline (years)", y = "Estimated overall survival") + theme_pub(12) + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#EEF2F5"))
f3 <- saveplot(p3, "Stage02_overall_survival_10_cancers_95CI", 17, 8)
report <- capture.output({
    cat("=====================================================================\nPROJECT 3 | STAGE 02 \342\200\224 REBUILD SCIENTIFIC CHECKPOINT\n=====================================================================\n")
    cat("1. COHORT FLOW (COUNTS NOT ALL THE SAME UNIT)\n")
    print(flow)
    cat("\n2. SELECTED CANCER COHORTS\n")
    print(cohort)
    cat("\n3. NONEXCLUSIVE EXCLUSION FLAGS (DO NOT ADD COUNTS)\n")
    print(exclusions)
    cat("\n4. ZERO-DAY RECORDS (COUNTS BY CANCER/EVENT)\n")
    print(zero[, .(N = .N), by = .(Cancer, Event)][order(Cancer, Event)])
    cat("\n5. OVERALL SURVIVAL LANDMARKS (95% CI; NA WHEN UNESTIMABLE)\n")
    print(landmarks)
    cat("\n6. PATIENT / BARCODE AND SURVIVAL CHECKPOINTS\n")
    cat("Eligible patients: ", nrow(eligible), " | deaths: ", sum(eligible$OS == 1), " | censored: ",
        sum(eligible$OS == 0), "\n", sep = "")
    cat("Duplicate eligible patient IDs: ", anyDuplicated(eligible$patient_id), " | duplicate eligible sample IDs: ",
        anyDuplicated(eligible$sample_id), "\n", sep = "")
    cat("Duplicate expression HEADER names: ", length(ambig), " (ambiguous records excluded)\n", sep = "")
    cat("Primary specimen type 01 agreement: ", all(eligible$barcode_type == 1L), "\n", sep = "")
    cat("N death before 30 days in eligible cohort: ", sum(eligible$OS == 1 & eligible$time < 30), "\n",
        sep = "")
    cat("\n7. FIGURES GENERATED\n", paste(c(f1, f2, f3), collapse = "\n"), "\n", sep = "")
    cat("\n8. SCIENTIFIC INTERPRETATION\n")
    cat("Survival time zero is excluded as previously specified; zero-day events are reported for sensitivity planning.\n")
    cat("Cohort event fraction is descriptive, not risk-adjusted; KM curves do not adjust for cancer stage, grade, or age.\n")
    cat("KM 95% CIs reflect uncertainty under independent censoring assumptions; late follow-up may be sparse.\n")
    cat("Truncated TCGA identifiers with multiple distinct expression columns remain excluded pending aliquot metadata.\n")
    cat("Stage 02 is a scientific checkpoint: no gene-expression Cox models or significance testing performed.\n")
    cat("Original Data and Results were not modified.\n")
    cat("\nOUTPUT: ", out, "\n=====================================================================\n",
        sep = "")
})
writeLines(report, file.path(out, "Stage02_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
