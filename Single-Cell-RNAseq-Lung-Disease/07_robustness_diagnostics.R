# Project 1 | Human lung single-cell RNA sequencing Run sample-level sensitivity analyses and
# targeted biological diagnostics.  Source dataset: GEO GSE136831. Clinical samples, rather than
# individual cells, define biological replication for inferential analyses.

options(stringsAsFactors = FALSE)

resolve_project_root <- function() {
    configured <- Sys.getenv("PROJECT1_ROOT", unset = "")
    if (nzchar(configured))
        return(normalizePath(configured, winslash = "/", mustWork = TRUE))
    command <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(command))
        return(normalizePath(dirname(sub("^--file=", "", command[[1]])), winslash = "/", mustWork = TRUE))
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

PROJECT_ROOT <- resolve_project_root()

# Shared validation and modelling utilities

# Shared utilities for the revised GSE136831 portfolio workflow.  These functions deliberately keep
# file handling, statistical fitting, and export validation separate from the biological questions
# in later scripts.  Every run writes to a new directory, so an old figure cannot be mistaken for
# output produced by the current code.

suppressPackageStartupMessages(library(data.table))

assert_columns <- function(x, required, label = deparse(substitute(x))) {
    absent <- setdiff(required, names(x))
    if (length(absent)) {
        stop(label, " is missing columns: ", paste(absent, collapse = ", "))
    }
    invisible(TRUE)
}

write_table <- function(x, name) {
    destination <- file.path(TABLE_DIR, name)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(x, destination)
    if (!file.exists(destination) || file.info(destination)$size == 0) {
        stop("Failed to write table: ", destination)
    }
    destination
}

read_table <- function(name, columns = NULL) {
    path <- file.path(TABLE_DIR, name)
    if (!file.exists(path))
        stop("Missing revised result: ", path)
    x <- data.table::fread(path)
    if (!is.null(columns))
        assert_columns(x, columns, name)
    x
}

# Create a clinical-sample-by-cell-class table that includes true zero observations.  Zero
# completion matters because omission of absent classes biases disease comparisons toward clinical
# samples in which a class happened to be recovered.
complete_donor_composition <- function(metadata, class_column, subset_rows = NULL) {
    x <- if (is.null(subset_rows))
        copy(metadata) else copy(metadata[subset_rows])
    classes <- sort(unique(as.character(x[[class_column]])))
    donor_map <- unique(metadata[, .(Subject_Identity, Disease_Identity)])
    observed <- x[, .(Cells = .N), by = c("Subject_Identity", "Disease_Identity", class_column)]
    grid <- data.table::CJ(Subject_Identity = donor_map$Subject_Identity, Class = classes, unique = TRUE)
    grid <- merge(grid, donor_map, by = "Subject_Identity", all.x = TRUE)
    setnames(observed, class_column, "Class")
    grid <- merge(grid, observed[, .(Subject_Identity, Class, Cells)], by = c("Subject_Identity", "Class"),
        all.x = TRUE)
    grid[is.na(Cells), Cells := 0L]
    totals <- x[, .(TotalCells = .N), by = Subject_Identity]
    grid <- merge(grid, totals, by = "Subject_Identity", all.x = TRUE)
    grid[is.na(TotalCells), TotalCells := 0L]
    grid[, Percentage := fifelse(TotalCells > 0, 100 * Cells/TotalCells, 0)]
    grid[]
}

# Aggregate raw counts by clinical sample after enforcing a minimum number of cells.  The function
# returns unfiltered gene counts because harmonized sensitivity analyses need to define one common
# gene universe across several thresholds.
make_pseudobulk <- function(count_matrix, metadata, cell_selector, disease_groups, min_cells, label) {
    selected <- which(cell_selector & metadata$Disease_Identity %in% disease_groups)
    if (!length(selected))
        stop("No cells selected for ", label)

    md <- copy(metadata[selected])
    donor_cells <- md[, .(Cells = .N), by = Subject_Identity]
    eligible <- donor_cells[Cells >= min_cells, Subject_Identity]
    selected <- selected[metadata$Subject_Identity[selected] %in% eligible]
    md <- copy(metadata[selected])
    if (!nrow(md))
        stop("No clinical samples passed the cell threshold for ", label)

    raw <- BPCells::pseudobulk_matrix(mat = count_matrix[, selected, drop = FALSE], cell_groups = as.character(md$Subject_Identity),
        method = "sum")
    raw <- as.matrix(raw)
    samples <- unique(md[, .(Subject_Identity, Disease_Identity)])
    samples <- samples[match(colnames(raw), Subject_Identity)]
    cell_totals <- md[, .(Cells = .N), by = Subject_Identity]
    samples[, Cells := cell_totals$Cells[match(Subject_Identity, cell_totals$Subject_Identity)]]
    samples[, Disease_Identity := factor(as.character(Disease_Identity), levels = disease_groups)]

    if (anyNA(samples$Subject_Identity) || anyNA(samples$Disease_Identity)) {
        stop("Clinical-sample or disease alignment failed for ", label)
    }
    if (!identical(as.character(samples$Subject_Identity), colnames(raw))) {
        stop("Pseudobulk column order failed for ", label)
    }
    group_sizes <- samples[, .N, by = Disease_Identity]
    if (!all(disease_groups %in% as.character(group_sizes$Disease_Identity)) || any(group_sizes$N < 3L)) {
        stop("Fewer than three eligible clinical samples in a disease group for ", label)
    }

    list(raw = raw, samples = samples, donor_cells = cell_totals, selected_cells = selected)
}

# Fit an edgeR quasi-likelihood model on clinical-sample pseudobulks. Disease-only models are
# primary. The optional recovered-cell covariate is a sensitivity analysis, because cell recovery
# may be both technical and biologically disease-related.
fit_pseudobulk <- function(pb, contrasts, gene_symbols, gene_universe = NULL, adjust_cell_count = FALSE,
    label = "pseudobulk") {
    dge <- edgeR::DGEList(counts = pb$raw)
    disease_design <- model.matrix(~0 + Disease_Identity, data = pb$samples)
    colnames(disease_design) <- levels(pb$samples$Disease_Identity)

    if (is.null(gene_universe)) {
        keep <- edgeR::filterByExpr(dge, design = disease_design)
    } else {
        keep <- rownames(dge) %in% gene_universe
    }
    if (sum(keep) < 100L)
        stop("Too few expressed genes for ", label)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)

    design <- disease_design
    if (adjust_cell_count) {
        design <- cbind(design, log2_cells = log2(pb$samples$Cells))
    }
    if (qr(design)$rank != ncol(design)) {
        stop("Non-identifiable design for ", label)
    }

    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    results <- lapply(names(contrasts), function(comparison) {
        contrast <- contrasts[[comparison]]
        if (length(contrast) != ncol(design)) {
            stop("Contrast length mismatch for ", label, ": ", comparison)
        }
        test <- edgeR::glmQLFTest(fit, contrast = contrast)
        tab <- as.data.table(edgeR::topTags(test, n = Inf)$table, keep.rownames = "gene_id")
        tab[, Gene := gene_symbols[match(gene_id, names(gene_symbols))]]
        tab[, `:=`(Comparison = comparison, EffectClass = fifelse(FDR < 0.05 & abs(logFC) >= 1, "FDR < 0.05 and |logFC| >= 1",
            fifelse(FDR < 0.05, "FDR < 0.05, small effect", "Not FDR significant")))]
        list(table = tab, test = test)
    })
    names(results) <- names(contrasts)
    list(dge = dge, fit = fit, samples = pb$samples, results = results, logCPM = edgeR::cpm(dge, log = TRUE,
        prior.count = 2), design = design)
}

run_hallmark_gsea <- function(test, pathways, comparison) {
    tab <- test$table
    statistic <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
    names(statistic) <- rownames(tab)
    statistic <- sort(statistic[is.finite(statistic)], decreasing = TRUE)
    set.seed(42)
    # A single process is intentional: on Windows, the default BiocParallel backend repeatedly
    # launches workers and can vary across installations.
    ans <- as.data.table(fgsea::fgsea(pathways = pathways, stats = statistic, minSize = 15, maxSize = 500,
        eps = 0, nproc = 1))
    ans[, Comparison := comparison]
    ans
}

# Pairwise quasibinomial models retain the clinical sample as the observation while using its
# evaluated-cell count as the denominator. They are preferable to an unweighted rank test when
# clinical-sample cell recovery varies.
pairwise_quasibinomial <- function(x, success, total, group = "Disease", family_label = "All") {
    pairs <- combn(c("Control", "IPF", "COPD"), 2, simplify = FALSE)
    out <- rbindlist(lapply(pairs, function(pair) {
        z <- copy(x[get(group) %in% pair])
        z[, GroupModel := factor(get(group), levels = pair)]
        if (uniqueN(z$GroupModel) != 2L || any(table(z$GroupModel) < 3L)) {
            return(data.table(Family = family_label, Group1 = pair[1], Group2 = pair[2], N1 = sum(z$GroupModel ==
                pair[1]), N2 = sum(z$GroupModel == pair[2]), LogOdds = NA_real_, OddsRatio = NA_real_,
                CI_low = NA_real_, CI_high = NA_real_, PValue = NA_real_))
        }
        fit <- glm(cbind(z[[success]], z[[total]] - z[[success]]) ~ GroupModel, family = quasibinomial(),
            data = z)
        co <- summary(fit)$coefficients[2, ]
        estimate <- unname(co[1])
        standard_error <- unname(co[2])
        data.table(Family = family_label, Group1 = pair[1], Group2 = pair[2], N1 = sum(z$GroupModel ==
            pair[1]), N2 = sum(z$GroupModel == pair[2]), LogOdds = estimate, OddsRatio = exp(estimate),
            CI_low = exp(estimate - 1.96 * standard_error), CI_high = exp(estimate + 1.96 * standard_error),
            PValue = unname(co[4]))
    }), fill = TRUE)
    out[, FDR := p.adjust(PValue, method = "BH")]
    out[]
}

save_webp <- function(plot, filename, width = 11, height = 7, quality = 94, dpi = 190) {
    destination <- file.path(FIGURE_DIR, filename)
    temporary_png <- tempfile(fileext = ".png")
    on.exit(unlink(temporary_png), add = TRUE)
    ggplot2::ggsave(temporary_png, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white",
        limitsize = FALSE)
    image <- png::readPNG(temporary_png)
    webp::write_webp(image, destination, quality = quality)
    if (!file.exists(destination) || file.info(destination)$size < 1000) {
        stop("Figure export failed: ", destination)
    }
    destination
}

portfolio_theme <- function(base_size = 11) {
    ggplot2::theme_minimal(base_size = base_size) + ggplot2::theme(text = ggplot2::element_text(family = "sans",
        colour = "#111111"), plot.title = ggplot2::element_text(size = 16, face = "bold", colour = "#000000"),
        plot.subtitle = ggplot2::element_blank(), plot.caption = ggplot2::element_blank(), axis.title = ggplot2::element_text(face = "bold",
            colour = "#000000"), axis.text = ggplot2::element_text(colour = "#111111"), strip.background = ggplot2::element_rect(fill = "#F1F5F9",
            colour = NA), strip.text = ggplot2::element_text(face = "bold"), panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_line(colour = "#E7EDF2",
            linewidth = 0.3), plot.background = ggplot2::element_rect(fill = "white", colour = NA), panel.background = ggplot2::element_rect(fill = "white",
            colour = NA))
}

# Reproducible output locations and checkpoint handling

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(BPCells)
    library(edgeR)
    library(fgsea)
})


RUN_DIR <- file.path(PROJECT_ROOT, "outputs")
TABLE_DIR <- file.path(RUN_DIR, "tables")
FIGURE_DIR <- file.path(RUN_DIR, "figures")
CACHE_DIR <- file.path(RUN_DIR, "checkpoints")
for (p in c(RUN_DIR, TABLE_DIR, FIGURE_DIR, CACHE_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

load_checkpoint <- function(name) {
    path <- file.path(CACHE_DIR, paste0(name, ".rds"))
    if (!file.exists(path))
        stop("Missing upstream checkpoint: ", path)
    readRDS(path)
}
store_checkpoint <- function(name, values) {
    path <- file.path(CACHE_DIR, paste0(name, ".rds"))
    temp <- tempfile(tmpdir = CACHE_DIR, fileext = ".rds")
    on.exit(unlink(temp), add = TRUE)
    saveRDS(values, temp, compress = FALSE)
    if (!file.rename(temp, path))
        stop("Cannot save checkpoint: ", path)
    invisible(path)
}
restore <- function(x) list2env(x, envir = .GlobalEnv)

# ----------------------------------------------------------------------------- Analysis
# -----------------------------------------------------------------------------

state <- load_checkpoint("01_data")
metadata <- fread(file.path(TABLE_DIR, "cell_metadata.csv"), showProgress = FALSE)
metadata[, Disease_Identity := factor(Disease_Identity, levels = c("Control", "IPF", "COPD"))]
genes <- fread(file.path(PROJECT_ROOT, "Data", "GSE136831_AllCells.GeneIDs.txt.gz"), showProgress = FALSE)
setnames(genes, c("Ensembl_GeneID", "HGNC_EnsemblAlt_GeneID"), c("gene_id", "gene_symbol"))
gene_symbols <- setNames(as.character(genes$gene_symbol), genes$gene_id)
counts <- BPCells::open_matrix_dir(state$cache)[, state$keep, drop = FALSE]
if (!identical(colnames(counts), as.character(metadata$CellBarcode_Identity))) stop("Count/meta barcode mismatch")
restore(load_checkpoint("03_primary_models"))
restore(load_checkpoint("04_gene_sets"))
restore(load_checkpoint("05_macrophage_models"))
restore(load_checkpoint("06_atii_models"))
metadata[, MacrophageClass := fifelse(Manuscript_Identity == "Macrophage_Alveolar", fifelse(grepl("CellCycle",
    Subclass_Cell_Identity), "Cell-cycle", fifelse(grepl("Outlier", Subclass_Cell_Identity), "Outlier",
    "Core")), NA_character_)]
macrophage_rows <- which(metadata$Manuscript_Identity == "Macrophage_Alveolar")
macrophage_representation <- complete_donor_composition(metadata, "MacrophageClass", macrophage_rows)
# ------------------------------------------------------------ Fibroblast clinical-sample influence
# ------------------------------------------------------------ Research question: Do the leading
# fibroblast effects retain their direction and magnitude after each clinical sample is omitted,
# rather than merely crossing an FDR threshold?
top_fibroblast <- fib_fit$results$IPF_vs_Control$table[order(FDR)][1:20, gene_id]
fib_loo <- rbindlist(lapply(seq_len(ncol(fib_fit$dge)), function(excluded) {
    retain <- seq_len(ncol(fib_fit$dge)) != excluded
    dge <- fib_fit$dge[, retain, keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge)
    samples <- fib_fit$samples[retain]
    design <- model.matrix(~0 + Disease_Identity, data = samples)
    colnames(design) <- c("Control", "IPF")
    dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
    test <- edgeR::glmQLFTest(fit, contrast = c(-1, 1))
    tab <- as.data.table(edgeR::topTags(test, n = Inf)$table, keep.rownames = "gene_id")
    tab[gene_id %in% top_fibroblast, .(ExcludedDonor = colnames(fib_fit$dge)[excluded], gene_id, Gene = gene_symbols[gene_id],
        logFC, FDR)]
}), fill = TRUE)
write_table(fib_loo, "sensitivity/fibroblast_top20_leave_one_out.csv")

# ------------------------------------------------------------ Clinical-sample-level macrophage
# cell-cycle evidence ------------------------------------------------------------ Research
# question: Does MKI67/TOP2A co-detection differ among diseases after accounting for the number of
# macrophages evaluated per clinical sample?  A quasibinomial model treats each clinical sample as
# an observation and uses its evaluated macrophage count as the denominator. Co-detection is still
# an RNA phenotype and is not proof of completed cell division.
gene_vector <- function(symbol, cell_rows) {
    ids <- genes[gene_symbol == symbol, gene_id]
    if (length(ids) != 1L)
        stop("Expected one Ensembl ID for ", symbol)
    as.numeric(as.matrix(counts[ids, cell_rows, drop = FALSE]))
}
mki67 <- gene_vector("MKI67", macrophage_rows)
top2a <- gene_vector("TOP2A", macrophage_rows)
cycle_detection <- metadata[macrophage_rows, .(Subject_Identity, Disease = as.character(Disease_Identity))]
cycle_detection[, Positive := as.integer(mki67 > 0 & top2a > 0)]
cycle_detection <- cycle_detection[, .(Positive = sum(Positive), Total = .N), by = .(Subject_Identity,
    Disease)]
cycle_detection[, Percentage := 100 * Positive/Total]
write_table(cycle_detection, "macrophage/cell_cycle_codetection_by_donor.csv")
cycle_models <- pairwise_quasibinomial(cycle_detection, "Positive", "Total")
write_table(cycle_models, "macrophage/cell_cycle_quasibinomial.csv")

# The source cycling-subclass representation is evaluated separately. Agreement with transcript
# co-detection is supportive but not independent validation, because both definitions derive from
# transcriptomic information.
cycle_subclass <- macrophage_representation[Class == "Cell-cycle", .(Subject_Identity, Disease = as.character(Disease_Identity),
    Positive = Cells, Total = TotalCells, Percentage)]
cycle_subclass_models <- pairwise_quasibinomial(cycle_subclass, "Positive", "Total", family_label = "Cycling subclass")
write_table(cycle_subclass_models, "macrophage/cycling_subclass_quasibinomial.csv")

# ------------------------------------------------------------ Macrophage identity and possible
# mixed-profile signals ------------------------------------------------------------ Research
# question: Are PPBP, PF4, MMP7, and KRT19 detections widespread within core macrophages, or
# concentrated in cycling/outlier subclasses and particular clinical samples?  Detection alone
# cannot distinguish ambient RNA, residual multiplets, phagocytosed material, or authentic
# transcription. Results are therefore contamination diagnostics rather than automatic evidence for
# removal or biological synthesis.
diagnostic_genes <- c("PPBP", "PF4", "MMP7", "KRT19")
detection_matrix <- sapply(diagnostic_genes, function(symbol) {
    gene_vector(symbol, macrophage_rows) > 0
})
diagnostic_cell_metadata <- metadata[macrophage_rows, .(Subject_Identity, Disease = as.character(Disease_Identity),
    Class = MacrophageClass)]
marker_prevalence <- rbindlist(lapply(seq_along(diagnostic_genes), function(i) {
    z <- copy(diagnostic_cell_metadata)
    z[, Detected := detection_matrix[, i]]
    z[, .(Detected = sum(Detected), Total = .N, Percent = 100 * mean(Detected)), by = .(Subject_Identity,
        Disease, Class)][, Gene := diagnostic_genes[i]][]
}), fill = TRUE)
write_table(marker_prevalence, "macrophage/diagnostic_marker_prevalence.csv")

core_marker_models <- rbindlist(lapply(diagnostic_genes, function(symbol) {
    x <- marker_prevalence[Gene == symbol & Class == "Core"]
    pairwise_quasibinomial(x, "Detected", "Total", family_label = symbol)
}), fill = TRUE)
core_marker_models[, GlobalFDR := p.adjust(PValue, method = "BH")]
write_table(core_marker_models, "macrophage/core_marker_quasibinomial.csv")

# ------------------------------------------------------------ Plot-ready marker summaries
# ------------------------------------------------------------ Balanced deterministic samples bound
# memory. These summaries are descriptive visual checks and are not used for
# differential-expression inference.
marker_panel <- c("SFTPC", "SFTPA1", "AGER", "PDPN", "KRT5", "KRT17", "KRT8", "COL1A1", "DCN", "ACTA2",
    "PECAM1", "VWF", "C1QA", "MARCO", "LYZ", "CD3D", "MS4A1", "MKI67", "TOP2A", "PPBP", "PF4", "MMP7",
    "KRT19", "IGHM", "CCNA2")
marker_summary <- function(indices, groups, group_name) {
    gene_rows <- which(genes$gene_symbol %in% marker_panel)
    out <- rbindlist(lapply(unique(as.character(groups)), function(group_value) {
        selected <- indices[as.character(groups) == group_value]
        mat <- as.matrix(counts[gene_rows, selected, drop = FALSE])
        library_size <- metadata$nUMI[selected]
        data.table(Group = group_value, Gene = genes$gene_symbol[gene_rows], pct_exp = 100 * rowMeans(mat >
            0), avg_exp = rowMeans(log1p(t(t(mat)/pmax(library_size, 1) * 10000))))
    }), fill = TRUE)
    out[, avg_exp_scaled := {
        standard_deviation <- sd(avg_exp)
        if (is.finite(standard_deviation) && standard_deviation > 0) {
            pmax(-2.5, pmin(2.5, (avg_exp - mean(avg_exp))/standard_deviation))
        } else rep(0, .N)
    }, by = Gene]
    setnames(out, "Group", group_name)
    out[]
}
set.seed(42)
marker_indices <- unlist(lapply(split(seq_len(nrow(metadata)), metadata$Manuscript_Identity), function(x) sample(x,
    min(length(x), 1500L))), use.names = FALSE)
marker_table <- marker_summary(marker_indices, metadata$Manuscript_Identity[marker_indices], "CellType")
write_table(marker_table, "plot_data/celltype_markers.csv")

set.seed(42)
macro_marker_indices <- unlist(lapply(split(macrophage_rows, interaction(metadata$Disease_Identity[macrophage_rows],
    metadata$MacrophageClass[macrophage_rows], drop = TRUE)), function(x) sample(x, min(length(x), 1500L))),
    use.names = FALSE)
macro_marker_groups <- paste(metadata$Disease_Identity[macro_marker_indices], metadata$MacrophageClass[macro_marker_indices],
    sep = " | ")
macro_marker_table <- marker_summary(macro_marker_indices, macro_marker_groups, "DiseaseClass")
write_table(macro_marker_table, "plot_data/macrophage_markers.csv")


store_checkpoint("07_analysis_models", list(fib_fit = fib_fit, atii_fit = atii_fit, macro_fit = macro_fit,
    core_fit = core_fit, core_adjusted = core_adjusted, cycling_fit = cycling_fit, atii_subclass_fits = atii_subclass_fits,
    pathways = pathways))
