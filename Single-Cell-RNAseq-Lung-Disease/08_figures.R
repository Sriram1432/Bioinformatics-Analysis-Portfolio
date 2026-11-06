# Project 1 | Human lung single-cell RNA sequencing Regenerate the approved scientific figure atlas
# from saved analytical outputs.  Source dataset: GEO GSE136831. Clinical samples, rather than
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

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(scales)
})
if (!file.exists(file.path(CACHE_DIR, "07_analysis_models.rds"))) {
    stop("Run 07_robustness_diagnostics.R before generating figures")
}
umap_path <- file.path(TABLE_DIR, "umap_plot_data.csv")
if (!file.exists(umap_path)) stop("Required UMAP embedding missing: ", umap_path)
# Generate the 20 approved scientific figures.  Every image is generated from analytical result
# tables and source-derived objects.  The manifest stores a SHA-256 digest, so stale images cannot
# pass validation.

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
})

disease_palette <- c(Control = "#10B5D1", IPF = "#FF527B", COPD = "#FFB52E")
figure_records <- list()

register_figure <- function(number, title, plot, width = 11, height = 7) {
    if (number %in% c(7L, 12L))
        return(invisible(NULL))
    filename <- sprintf("Figure%02d_%s.webp", number, gsub("_+", "_", gsub("[^A-Za-z0-9]+", "_", title)))
    path <- save_webp(plot, filename, width, height)
    figure_records[[length(figure_records) + 1L]] <<- data.table(Figure = number, Title = title, File = basename(path),
        SHA256 = digest::digest(file = path, algo = "sha256"), Bytes = file.info(path)$size, CreatedUTC = format(Sys.time(),
            tz = "UTC", usetz = TRUE))

    invisible(path)
}

volcano_plot <- function(tab, title) {
    x <- copy(tab)
    x[, Category := fifelse(FDR < 0.05 & abs(logFC) >= 1, fifelse(logFC > 0, "Higher in IPF", "Lower in IPF"),
        fifelse(FDR < 0.05, "FDR only (<2-fold)", "Not significant"))]
    x[, LogFDR := -log10(pmax(FDR, .Machine$double.xmin))]
    label <- x[Category %in% c("Higher in IPF", "Lower in IPF") & !is.na(Gene) & nzchar(Gene) & !grepl("^ENSG",
        Gene)][order(FDR)][1:min(.N, 12)]
    ggplot(x, aes(logFC, LogFDR)) + geom_point(aes(colour = Category), size = 0.85, alpha = 0.6) + geom_vline(xintercept = c(-1,
        1), linetype = "dotted", colour = "#94A3B8") + geom_hline(yintercept = -log10(0.05), linetype = "dotted",
        colour = "#94A3B8") + geom_text_repel(data = label, aes(label = Gene), colour = "#17273B", seed = 42,
        size = 3, max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0) + scale_colour_manual(values = c(`Higher in IPF` = disease_palette[["IPF"]],
        `Lower in IPF` = disease_palette[["Control"]], `FDR only (<2-fold)` = "#8B60C9", `Not significant` = "#D4DCE5")) +
        labs(title = title, x = expression(log[2](" fold change")), y = expression(-log[10]("FDR")),
            colour = NULL) + portfolio_theme() + theme(legend.position = "top")
}

gsea_lollipop <- function(tab, title, top_n = 14L) {
    x <- copy(tab)[is.finite(NES) & !is.na(padj)][order(padj)][1:min(.N, top_n)]
    x[, Pathway := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
    x[, Pathway := factor(Pathway, levels = Pathway[order(NES)])]
    ggplot(x, aes(NES, Pathway)) + geom_segment(aes(x = 0, xend = NES, yend = Pathway), colour = "#D9E2EA",
        linewidth = 0.8) + geom_point(aes(fill = -log10(pmax(padj, 1e-300))), shape = 21, size = 5, colour = "white") +
        geom_text(aes(label = sprintf("%+.2f", NES), hjust = ifelse(NES > 0, -0.25, 1.25)), size = 3) +
        scale_fill_gradient(low = "#B7DDEA", high = "#A21D46", name = expression(-log[10]("adjusted P"))) +
        scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) + labs(title = title, x = "Normalized enrichment score",
        y = NULL, fill = expression(-log[10]("adjusted P"))) + portfolio_theme()
}

# 01 | Cohort structure
cohort <- read_table("cohort.csv")
p1 <- ggplot(cohort, aes(Disease_Identity, Cells, fill = Disease_Identity)) + geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%s cells
%s clinical samples", scales::comma(Cells), Donors)), vjust = -0.35,
        fontface = "bold", lineheight = 1.2) + scale_fill_manual(values = disease_palette) + scale_y_continuous(labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.2))) + labs(title = "Clinical Sample and Cell Recovery Across Disease Groups",
    x = NULL, y = "Retained cells") + portfolio_theme(12) + theme(legend.position = "none")
register_figure(1, "Clinical Sample and Cell Recovery Across Disease Groups", p1, 9.5, 6)

# 02 | Clinical-sample-level quality landscape
qc <- read_table("qc_by_donor.csv")
qc_long <- melt(qc, id.vars = c("Subject_Identity", "Disease_Identity"), measure.vars = c("MedianGenes",
    "MedianUMI", "MedianMT"), variable.name = "Metric", value.name = "Value")
qc_long[, Metric := factor(Metric, levels = c("MedianGenes", "MedianUMI", "MedianMT"), labels = c("Median genes",
    "Median UMIs", "Median mitochondrial %"))]
p2 <- ggplot(qc_long, aes(Disease_Identity, Value, fill = Disease_Identity)) + geom_violin(alpha = 0.18,
    trim = FALSE) + geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12,
    seed = 42), shape = 21, size = 1.7, alpha = 0.8) + facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = disease_palette) + labs(title = "Sequencing Quality Across Clinical Samples",
    x = NULL, y = NULL) + portfolio_theme() + theme(legend.position = "none")
register_figure(2, "Sequencing Quality Across Clinical Samples", p2, 13, 5.5)

# 03 | Source cell-type landscape
umap <- read_table("umap_plot_data.csv")
assert_columns(umap, c("UMAP_1", "UMAP_2", "Manuscript_Identity"), "UMAP export")
set.seed(42)
umap_draw <- umap[sample.int(nrow(umap), min(65000L, nrow(umap)))]
umap_centres <- umap_draw[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = Manuscript_Identity]
p3 <- ggplot(umap_draw, aes(UMAP_1, UMAP_2, colour = Manuscript_Identity)) + geom_point(size = 0.12,
    alpha = 0.62) + coord_equal() + geom_text_repel(data = umap_centres, aes(label = Manuscript_Identity),
    colour = "#111111", size = 2.5, fontface = "bold", seed = 42, box.padding = 0.25, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE) + labs(title = "Human Lung Cellular Atlas", x = "UMAP 1",
    y = "UMAP 2") + portfolio_theme(9) + theme(panel.grid = element_blank(), legend.position = "none")
register_figure(3, "Human Lung Cellular Atlas", p3, 14, 9)

# 04 | Disease distribution on shared embedding
umap[, Disease_Identity := factor(Disease_Identity, levels = names(disease_palette))]
p4 <- ggplot(umap, aes(UMAP_1, UMAP_2, colour = Disease_Identity)) + geom_point(size = 0.12, alpha = 0.6) +
    facet_wrap(~Disease_Identity, nrow = 1) + coord_equal() + scale_colour_manual(values = disease_palette) +
    labs(title = "Disease Groups Across the Lung Transcriptional Landscape", x = "UMAP 1", y = "UMAP 2") +
    portfolio_theme() + theme(panel.grid = element_blank(), legend.position = "none")
register_figure(4, "Disease Groups Across the Lung Transcriptional Landscape", p4, 14, 5.8)

# 05 | Canonical marker coherence
markers <- read_table("plot_data/celltype_markers.csv")
p5 <- ggplot(markers, aes(Gene, CellType)) + geom_point(aes(size = pct_exp, fill = avg_exp_scaled), shape = 21,
    colour = "#17273B", stroke = 0.18) + scale_size_continuous(range = c(0.4, 8), name = "Cells detected (%)") +
    scale_fill_gradient2(low = "#3E7CB1", mid = "white", high = "#B23A48", midpoint = 0, name = "Scaled mean") +
    labs(title = "Canonical Markers Support Lung Cell-Type Annotations", x = NULL, y = NULL) + portfolio_theme(10) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1))
register_figure(5, "Canonical Markers Support Lung Cell-Type Annotations", p5, 14, 12)

# 06 | Broad recovered composition
composition <- read_table("broad_cell_composition.csv")
focus <- c("ATII", "Fibroblast", "Myofibroblast", "Aberrant_Basaloid")
composition <- composition[CellType %in% focus]
composition[, CellType := factor(CellType, levels = focus)]
p6 <- ggplot(composition, aes(Disease_Identity, Percentage, fill = Disease_Identity)) + geom_boxplot(width = 0.52,
    outlier.shape = NA, alpha = 0.45) + geom_point(position = position_jitter(width = 0.12, seed = 42),
    shape = 21, size = 1.6, alpha = 0.8) + facet_wrap(~CellType, scales = "free_y", ncol = 2) + scale_fill_manual(values = disease_palette) +
    labs(title = "Recovered Lung Cell Fractions Across Disease Groups", x = NULL, y = "Recovered cell fraction (%)") +
    portfolio_theme() + theme(legend.position = "none")
register_figure(6, "Recovered Lung Cell Fractions Across Disease Groups", p6, 11, 8)

# 07 | Pairwise recovered-composition effects
prop <- read_table("audit/broad_composition_pairwise_propeller.csv", c("CellType", "Comparison", "DifferencePercentagePoints",
    "FDRWithinContrast"))
selected_types <- prop[, .(MinimumFDR = min(FDRWithinContrast)), by = CellType][order(MinimumFDR)][1:min(.N,
    18L), CellType]
prop <- prop[CellType %in% selected_types]
prop[, Comparison := factor(Comparison, levels = c("IPF_vs_Control", "COPD_vs_Control", "COPD_vs_IPF"),
    labels = c("IPF vs Control", "COPD vs Control", "COPD vs IPF"))]
prop[, CellType := factor(CellType, levels = rev(selected_types))]
p7 <- ggplot(prop, aes(Comparison, CellType, fill = DifferencePercentagePoints)) + geom_tile(colour = "white",
    linewidth = 0.5) + geom_point(data = prop[FDRWithinContrast < 0.05], shape = 21, fill = "white",
    colour = "#111111", size = 2, stroke = 0.5) + scale_fill_gradient2(low = disease_palette[["Control"]],
    mid = "white", high = disease_palette[["IPF"]], midpoint = 0, name = "Group 2 - Group 1
(percentage points)") +
    labs(title = "Differential Cell Representation Across Lung Diseases", x = NULL, y = NULL) + portfolio_theme(10) +
    theme(panel.grid = element_blank(), legend.position = "right")
register_figure(7, "Differential Cell Representation Across Lung Diseases", p7, 11, 8)

# 08 | Broad fibroblast differential expression
fib <- read_table("de/fibroblast_broad_IPF_vs_Control.csv")
p8 <- volcano_plot(fib, "IPF-Associated Fibroblast Transcriptional Remodeling")
register_figure(8, "IPF-Associated Fibroblast Transcriptional Remodeling", p8, 11, 7.2)

# 09 | Fibroblast leave-one-clinical-sample-out coefficients
loo <- read_table("sensitivity/fibroblast_top20_leave_one_out.csv")
loo_summary <- loo[, .(MedianLogFC = median(logFC), Low = min(logFC), High = max(logFC), WorstFDR = max(FDR),
    DirectionStable = uniqueN(sign(logFC)) == 1L), by = Gene][order(MedianLogFC)]
loo_summary[, Gene := factor(Gene, levels = Gene)]
p9 <- ggplot(loo_summary, aes(MedianLogFC, Gene)) + geom_errorbar(aes(xmin = Low, xmax = High), orientation = "y",
    width = 0.18, linewidth = 0.8, colour = "#2A9D8F") + geom_point(size = 3, colour = "#17273B") + geom_vline(xintercept = 0,
    linetype = "dotted", colour = "#94A3B8") + labs(title = "Fibroblast Expression Stability Across Clinical Sample Exclusions",
    x = expression(log[2](" fold change (IPF vs Control)")), y = NULL) + portfolio_theme()
register_figure(9, "Fibroblast Expression Stability Across Clinical Sample Exclusions", p9, 11, 8)

# 10 | Fibroblast Hallmark programs
fib_gsea <- read_table("gsea/fibroblast_broad.csv")
p10 <- gsea_lollipop(fib_gsea, "IPF-Associated Transcriptional Programs in Fibroblasts")
register_figure(10, "IPF-Associated Transcriptional Programs in Fibroblasts", p10, 11.5, 8)

# 11 | Broad ATII differential expression
atii <- read_table("de/atii_broad_IPF_vs_Control.csv")
p11 <- volcano_plot(atii, "IPF-Associated Transcriptional Changes in ATII Cells")
register_figure(11, "IPF-Associated Transcriptional Changes in ATII Cells", p11, 11, 7.2)

# 12 | Harmonized ATII threshold stability
stability <- read_table("comparisons/atii_harmonized_threshold_stability.csv")
stable_genes <- stability[StableFDRAndTwofold == TRUE][order(FDR_10), Gene]
selected_genes <- stable_genes
stability_long <- melt(stability[Gene %in% selected_genes], id.vars = c("gene_id", "Gene", "StableDirection",
    "StableFDRAndTwofold"), measure.vars = c("logFC_10", "logFC_15", "logFC_20"), variable.name = "Threshold",
    value.name = "logFC")
stability_long[, Threshold := factor(Threshold, levels = c("logFC_10", "logFC_15", "logFC_20"), labels = c("10",
    "15", "20"))]
end_labels <- stability_long[Threshold == "20"]
p12 <- ggplot(stability_long, aes(Threshold, logFC, group = Gene, colour = Gene)) + geom_hline(yintercept = 0,
    colour = "#CBD5E1") + geom_hline(yintercept = c(-1, 1), colour = "#94A3B8", linetype = "dotted") +
    geom_line(alpha = 0.7, linewidth = 0.75) + geom_point(size = 2) + geom_text_repel(data = end_labels,
    aes(label = Gene), direction = "y", hjust = -0.1, size = 2.7, seed = 42, force = 5, box.padding = 0.35,
    max.time = 5, min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE) + scale_x_discrete(expand = expansion(add = c(0.2,
    1.25))) + labs(title = "Robustness of ATII Gene Expression Findings", x = "Minimum ATII cells required per clinical sample",
    y = expression(log[2](" fold change (IPF vs Control)"))) + portfolio_theme(10) + theme(legend.position = "none")
register_figure(12, "Robustness of ATII Gene Expression Findings", p12, 13, 8)

# 13 | ATII subclass representation
atii_rep <- read_table("atii_subclass_representation.csv")
atii_rep <- atii_rep[TotalCells > 0 & Class %in% c("High surfactant", "Low surfactant")]
p13 <- ggplot(atii_rep, aes(Disease_Identity, Percentage, fill = Disease_Identity)) + geom_boxplot(width = 0.5,
    outlier.shape = NA, alpha = 0.45) + geom_point(position = position_jitter(width = 0.12, height = 0,
    seed = 42), shape = 21, size = 1.7) + facet_wrap(~Class, nrow = 1) + scale_fill_manual(values = disease_palette) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) + labs(title = "High- and Low-Surfactant ATII Cell Representation",
    x = NULL, y = "Fraction of recovered ATII cells (%)") + portfolio_theme() + theme(legend.position = "none")
register_figure(13, "High- and Low-Surfactant ATII Cell Representation", p13, 11, 6)

# 14 | Broad and subclass ATII coefficients
atii_compare <- read_table("comparisons/atii_broad_vs_subclasses.csv")
top_atii <- atii_compare[order(BroadFDR)][1:min(.N, 250)]
atii_compare_long <- rbindlist(list(top_atii[, .(Gene, BroadLogFC, SubclassLogFC = HighLogFC, SubclassFDR = HighFDR,
    Subclass = "High surfactant")], top_atii[, .(Gene, BroadLogFC, SubclassLogFC = LowLogFC, SubclassFDR = LowFDR,
    Subclass = "Low surfactant")]))
p14 <- ggplot(atii_compare_long, aes(BroadLogFC, SubclassLogFC)) + geom_hline(yintercept = 0, colour = "#D8E2EB") +
    geom_vline(xintercept = 0, colour = "#D8E2EB") + geom_abline(slope = 1, linetype = "dotted", colour = "#94A3B8") +
    geom_point(aes(colour = fifelse(SubclassFDR < 0.05, "Subclass BH-FDR < 0.05", "Not significant")),
        alpha = 0.65, size = 1.6) + facet_wrap(~Subclass, nrow = 1) + scale_colour_manual(values = c(`Subclass BH-FDR < 0.05` = disease_palette[["IPF"]],
    `Not significant` = "#B9C4CF")) + labs(title = "Gene Expression Differences Across ATII Subclasses",
    x = expression("Broad ATII " * log[2] * " fold change"), y = expression("Subclass " * log[2] * " fold change"),
    colour = NULL) + portfolio_theme()
register_figure(14, "Gene Expression Differences Across ATII Subclasses", p14, 12, 6.5)

# 15 | Broad macrophage IPF effect
macro_ipf <- read_table("de/macrophage_broad_IPF_vs_Control.csv")
p15 <- volcano_plot(macro_ipf, "IPF-Associated Transcriptional Changes in Alveolar Macrophages")
register_figure(15, "IPF-Associated Transcriptional Changes in Alveolar Macrophages", p15, 11, 7.2)

# 16 | Correct direct IPF-COPD evidence
macro_copd <- read_table("de/macrophage_broad_COPD_vs_Control.csv")
macro_direct <- read_table("de/macrophage_broad_IPF_vs_COPD.csv")
direct <- Reduce(function(a, b) merge(a, b, by = c("gene_id", "Gene")), list(macro_ipf[, .(gene_id, Gene,
    IPF_logFC = logFC, IPF_FDR = FDR)], macro_copd[, .(gene_id, Gene, COPD_logFC = logFC, COPD_FDR = FDR)],
    macro_direct[, .(gene_id, Gene, Direct_logFC = logFC, Direct_FDR = FDR)]))
direct[, Evidence := fifelse(Direct_FDR < 0.05 & abs(Direct_logFC) >= 1, fifelse(Direct_logFC > 0,
    "Greater in IPF: direct support", "Greater in COPD: direct support"), fifelse(IPF_FDR < 0.05 & COPD_FDR <
    0.05 & sign(IPF_logFC) == sign(COPD_logFC), "Concordant versus Control", "No direct twofold evidence"))]
p16 <- ggplot(direct, aes(COPD_logFC, IPF_logFC, colour = Evidence)) + geom_hline(yintercept = 0, colour = "#D8E2EB") +
    geom_vline(xintercept = 0, colour = "#D8E2EB") + geom_abline(slope = 1, linetype = "dotted", colour = "#94A3B8") +
    geom_point(size = 0.8, alpha = 0.5) + scale_colour_manual(values = c(`Greater in IPF: direct support` = disease_palette[["IPF"]],
    `Greater in COPD: direct support` = disease_palette[["COPD"]], `Concordant versus Control` = "#7B61A8",
    `No direct twofold evidence` = "#CBD5E1")) + coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 7)) +
    labs(title = "Direct Macrophage Expression Differences Between IPF and COPD", x = expression(log[2](" fold change (COPD vs Control)")),
        y = expression(log[2](" fold change (IPF vs Control)")), colour = NULL) + portfolio_theme(10) +
    theme(legend.position = "bottom")
register_figure(16, "Direct Macrophage Expression Differences Between IPF and COPD", p16, 11, 8)

# 17 | Macrophage subclass representation
macro_rep <- read_table("macrophage_subclass_representation.csv")
macro_rep[, Class := factor(Class, levels = c("Core", "Cell-cycle", "Outlier"), labels = c("Core macrophage transcriptomic group",
    "Cell-cycle-associated macrophage group", "Atypical macrophage transcriptomic group"))]
p17 <- ggplot(macro_rep, aes(Disease_Identity, Percentage, fill = Disease_Identity)) + geom_boxplot(width = 0.5,
    outlier.shape = NA, alpha = 0.4) + geom_point(position = position_jitter(width = 0.12, height = 0,
    seed = 42), shape = 21, size = 1.6, alpha = 0.8) + facet_wrap(~Class, scales = "fixed", nrow = 1) +
    scale_fill_manual(values = disease_palette) + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0,
    0.08))) + labs(title = "Macrophage Transcriptional States Across Disease Groups", x = NULL, y = "Fraction of recovered alveolar macrophages (%)") +
    portfolio_theme() + theme(legend.position = "none")
register_figure(17, "Macrophage Transcriptional States Across Disease Groups", p17, 14, 6)

# 18 | Core macrophage DE
core_ipf <- read_table("de/macrophage_core_IPF_vs_Control.csv")
p18 <- volcano_plot(core_ipf, "IPF-Associated Transcriptional Changes in Core Macrophages")
register_figure(18, "IPF-Associated Transcriptional Changes in Core Macrophages", p18, 11, 7.2)

# 19 | Broad versus core effects
broad_core <- read_table("comparisons/macrophage_broad_vs_core.csv")
bc <- broad_core[Comparison == "IPF_vs_Control"]
bc[, DisplayPattern := factor(Pattern, levels = c("Supported in broad and core", "Broad only at threshold",
    "Core only at threshold", "Neither"), labels = c("Meets criteria in both", "Meets criteria in broad only",
    "Meets criteria in core only", "Below joint criteria"))]
p19 <- ggplot(bc, aes(BroadLogFC, CoreLogFC, colour = DisplayPattern)) + geom_hline(yintercept = 0, colour = "#D8E2EB") +
    geom_vline(xintercept = 0, colour = "#D8E2EB") + geom_abline(slope = 1, linetype = "dotted", colour = "#94A3B8") +
    geom_point(size = 0.85, alpha = 0.55) + scale_colour_manual(values = c(`Meets criteria in both` = "#2A9D8F",
    `Meets criteria in broad only` = disease_palette[["IPF"]], `Meets criteria in core only` = "#3E7CB1",
    `Below joint criteria` = "#D3DBE4")) + labs(title = "Concordance of Broad and Core Macrophage Effects",
    x = expression("Broad macrophage " * log[2] * " fold change"), y = expression("Core macrophage " *
        log[2] * " fold change"), colour = NULL) + portfolio_theme(10) + theme(legend.position = "bottom")
register_figure(19, "Concordance of Broad and Core Macrophage Effects", p19, 11, 8)

# 20 | Primary versus cell-count-adjusted core macrophage pathways
primary_gsea <- read_table("gsea/macrophage_core_primary.csv")
adjusted_gsea <- read_table("gsea/macrophage_core_cellcount_adjusted.csv")
primary_gsea[, Model := "Primary: disease only"]
adjusted_gsea[, Model := "Sensitivity: + recovered cells"]
gsea_both <- rbindlist(list(primary_gsea, adjusted_gsea), fill = TRUE)
top_paths <- primary_gsea[Comparison == "IPF_vs_Control"][order(padj)][1:14, pathway]
gsea_both <- gsea_both[Comparison == "IPF_vs_Control" & pathway %in% top_paths]
gsea_both[, Pathway := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
path_levels <- rev(gsub("_", " ", sub("^HALLMARK_", "", top_paths)))
gsea_both[, Pathway := factor(Pathway, levels = path_levels)]
p20 <- ggplot(gsea_both, aes(NES, Pathway, colour = Model)) + geom_vline(xintercept = 0, colour = "#D8E2EB") +
    geom_line(aes(group = Pathway), colour = "#C7D1DB", linewidth = 0.6) + geom_point(aes(size = -log10(pmax(padj,
    1e-300))), alpha = 0.85) + scale_colour_manual(values = c(`Primary: disease only` = "#17273B", `Sensitivity: + recovered cells` = disease_palette[["IPF"]])) +
    scale_size_continuous(range = c(2.5, 7), name = expression(-log[10]("FDR"))) + labs(title = "Robustness of Macrophage Pathway Enrichment",
    x = "Normalized enrichment score", y = NULL, colour = "Model") + portfolio_theme(10)
register_figure(20, "Robustness of Macrophage Pathway Enrichment", p20, 12, 8)

# 21 | Clinical-sample-weighted cell-cycle evidence
cycle <- read_table("macrophage/cell_cycle_codetection_by_donor.csv")
cycle[, Disease := factor(Disease, levels = names(disease_palette))]
cycle_tests <- read_table("macrophage/cell_cycle_quasibinomial.csv")
ipf_cycle <- cycle_tests[Group1 == "Control" & Group2 == "IPF"]
cycle_annotation <- sprintf("IPF vs Control: OR %.2f; 95%% CI %.2f-%.2f; FDR %.2e", ipf_cycle$OddsRatio,
    ipf_cycle$CI_low, ipf_cycle$CI_high, ipf_cycle$FDR)
p21 <- ggplot(cycle, aes(Disease, Percentage, fill = Disease)) + geom_boxplot(width = 0.42, outlier.shape = NA,
    alpha = 0.35) + geom_point(aes(size = Total), position = position_jitter(width = 0.12, height = 0,
    seed = 42), shape = 21, alpha = 0.75) + scale_fill_manual(values = disease_palette) + scale_size_continuous(range = c(1.5,
    5)) + scale_y_sqrt(limits = c(0, NA), breaks = c(0, 1, 4, 9, 16, 25), expand = expansion(mult = c(0,
    0.18))) + annotate("label", x = 2, y = Inf, label = cycle_annotation, vjust = 1.3, size = 3.6, fontface = "bold",
    fill = "white") + labs(title = "MKI67 and TOP2A Co-detection in Lung Macrophages", x = NULL, y = "Cells with both transcripts detected (%)",
    size = "Macrophages evaluated") + portfolio_theme() + theme(legend.position = "right")
register_figure(21, "MKI67 and TOP2A Co-detection in Lung Macrophages", p21, 11, 7)

# 22 | Possible contamination and mixed-profile diagnostics
diagnostic <- read_table("macrophage/diagnostic_marker_prevalence.csv")
diagnostic[, Disease := factor(Disease, levels = names(disease_palette))]
diagnostic[, Class := factor(Class, levels = c("Core", "Cell-cycle", "Outlier"), labels = c("Core group",
    "Cell-cycle-associated group", "Atypical group"))]
p22 <- ggplot(diagnostic, aes(Disease, Percent, colour = Class)) + geom_boxplot(aes(group = interaction(Disease,
    Class)), outlier.shape = NA, width = 0.6, alpha = 0.12, position = position_dodge(width = 0.68)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.68, seed = 42), size = 1.15,
        alpha = 0.6) + facet_wrap(~Gene, scales = "free_y", ncol = 2) + scale_colour_manual(values = c(`Core group` = "#17273B",
    `Cell-cycle-associated group` = "#8B60C9", `Atypical group` = "#2A9D8F")) + scale_y_continuous(limits = c(0,
    NA), expand = expansion(mult = c(0, 0.08))) + labs(title = "Epithelial- and Platelet-Associated Transcripts in Macrophages",
    x = NULL, y = "Cells with detected transcript (%)", colour = NULL) + portfolio_theme(9) + theme(legend.position = "top")
register_figure(22, "Epithelial- and Platelet-Associated Transcripts in Macrophages", p22, 13, 9)

manifest <- rbindlist(figure_records)
setorder(manifest, Figure)
if (!identical(as.integer(manifest$Figure), setdiff(1:22, c(7L, 12L)))) {
    stop("Figure numbering is incomplete")
}
fwrite(manifest, file.path(RUN_DIR, "figure_manifest.csv"))
manifest <- fread(file.path(RUN_DIR, "figure_manifest.csv"), showProgress = FALSE)
if (nrow(manifest) != 20L || any(!file.exists(file.path(FIGURE_DIR, manifest$File)))) stop("Expected 20 scientifically approved figures")
