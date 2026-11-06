# Project 1 | Human lung single-cell RNA sequencing Build the lung cellular atlas and test
# donor-level recovered-cell composition.  Source dataset: GEO GSE136831. Clinical samples, rather
# than individual cells, define biological replication for inferential analyses.

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

# ------------------------------------------------------------ Broad cell composition
# ------------------------------------------------------------ Research question: Which
# source-annotated populations differ in recovered representation among disease groups when
# clinical samples, rather than individual cells, define replication?  The result concerns
# recovered dissociated cells, not absolute abundance in intact tissue. Zero-count
# clinical-sample/population combinations remain represented.
composition <- complete_donor_composition(metadata, "Manuscript_Identity")
setnames(composition, "Class", "CellType")
write_table(composition, "broad_cell_composition.csv")
propeller <- speckle::propeller(clusters = metadata$Manuscript_Identity, sample = metadata$Subject_Identity,
    group = metadata$Disease_Identity, transform = "asin", robust = TRUE)
propeller <- as.data.table(propeller, keep.rownames = "CellType")
write_table(propeller, "broad_composition_propeller.csv")

# ------------------------------------------------------------ Source subclasses and
# clinical-sample representation ------------------------------------------------------------
# Research question: Are broad macrophage and ATII signals associated with changing mixtures of
# transcriptionally defined source subclasses?  Subclasses are not independent ground truth: they
# were themselves defined from RNA profiles. Stratification reduces one form of mixture confounding
# but cannot prove cell identity or eliminate technical and clinical confounding.
metadata[, MacrophageClass := NA_character_]
metadata[Manuscript_Identity == "Macrophage_Alveolar", MacrophageClass := fifelse(grepl("CellCycle",
    Subclass_Cell_Identity), "Cell-cycle", fifelse(grepl("Outlier", Subclass_Cell_Identity), "Outlier",
    "Core"))]
metadata[, ATIIClass := NA_character_]
metadata[Manuscript_Identity == "ATII", ATIIClass := fifelse(grepl("High-Surfactants", Subclass_Cell_Identity),
    "High surfactant", fifelse(grepl("Low-Surfactants", Subclass_Cell_Identity), "Low surfactant", "Other"))]

macrophage_rows <- which(metadata$Manuscript_Identity == "Macrophage_Alveolar")
atii_rows <- which(metadata$Manuscript_Identity == "ATII")
macrophage_representation <- complete_donor_composition(metadata, "MacrophageClass", macrophage_rows)
atii_representation <- complete_donor_composition(metadata, "ATIIClass", atii_rows)
write_table(macrophage_representation, "macrophage_subclass_representation.csv")
write_table(atii_representation, "atii_subclass_representation.csv")

summarise_subclass_support <- function(rows, class_column, population) {
    # Materialise the aliased columns before grouping. This avoids data.table's non-standard
    # evaluation treating an aliased by-expression differently when the parent table is
    # simultaneously subset with an external integer vector.
    donor_cells <- metadata[rows, .(Subclass = as.character(get(class_column)), Disease = as.character(Disease_Identity),
        Subject_Identity = as.character(Subject_Identity))][, .(Cells = .N), by = .(Subclass, Disease,
        Subject_Identity)]
    result <- donor_cells[, .(DonorsObserved = .N, Donors10 = sum(Cells >= 10), Donors20 = sum(Cells >=
        20), MedianCells = as.numeric(median(Cells)), MaxCells = max(Cells)), by = .(Subclass, Disease)]
    result[, Population := population]
    setcolorder(result, c("Population", "Subclass", "Disease", setdiff(names(result), c("Population",
        "Subclass", "Disease"))))
    result[]
}
representation_support <- rbindlist(list(summarise_subclass_support(macrophage_rows, "MacrophageClass",
    "Macrophage"), summarise_subclass_support(atii_rows, "ATIIClass", "ATII")), fill = TRUE)
write_table(representation_support, "subclass_donor_support.csv")


write_table(metadata[, .(CellBarcode_Identity, Subject_Identity, Disease_Identity, Manuscript_Identity,
    Subclass_Cell_Identity, MacrophageClass, ATIIClass)], "annotated_cell_metadata.csv")

pairs <- list(c("Control", "IPF"), c("Control", "COPD"), c("IPF", "COPD"))
pairwise <- rbindlist(lapply(pairs, function(pair) {
    z <- metadata[Disease_Identity %in% pair]
    ans <- as.data.table(speckle::propeller(clusters = z$Manuscript_Identity, sample = z$Subject_Identity,
        group = factor(z$Disease_Identity, levels = pair), transform = "asin", robust = TRUE), keep.rownames = "CellType")
    data.table(CellType = ans$CellType, Comparison = paste0(pair[2L], "_vs_", pair[1L]), Group1 = pair[1L],
        Group2 = pair[2L], MeanPercentage1 = 100 * ans[[paste0("PropMean.", pair[1L])]], MeanPercentage2 = 100 *
            ans[[paste0("PropMean.", pair[2L])]], DifferencePercentagePoints = 100 * (ans[[paste0("PropMean.",
            pair[2L])]] - ans[[paste0("PropMean.", pair[1L])]]), TStatisticGroup1MinusGroup2 = ans$Tstatistic,
        PValue = ans$P.Value, FDRWithinContrast = ans$FDR)
}), fill = TRUE)
pairwise[, FDRGlobal := p.adjust(PValue, method = "BH")]
if (nrow(pairwise) != 114L) stop("Unexpected number of pairwise cellular composition comparisons")
write_table(pairwise, "audit/broad_composition_pairwise_propeller.csv")

set.seed(42)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required for atlas embedding")
seurat_metadata <- as.data.frame(metadata)
rownames(seurat_metadata) <- as.character(seurat_metadata$CellBarcode_Identity)
if (!identical(rownames(seurat_metadata), colnames(counts))) stop("Seurat metadata misalignment")
obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = seurat_metadata, min.cells = 0, min.features = 0)
obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- Seurat::SketchData(obj, ncells = 30000, method = "LeverageScore", sketched.assay = "sketch", verbose = FALSE)
Seurat::DefaultAssay(obj) <- "sketch"
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = 1:30, return.model = TRUE, verbose = FALSE, seed.use = 42)
emb <- Seurat::Embeddings(obj, "umap")
if (ncol(emb) != 2L || anyDuplicated(rownames(emb))) stop("UMAP coordinate integrity failure")
umap <- as.data.table(emb, keep.rownames = "CellBarcode_Identity")
setnames(umap, names(umap)[2:3], c("UMAP_1", "UMAP_2"))
umap <- merge(umap, metadata[, .(CellBarcode_Identity, Disease_Identity, Manuscript_Identity)], by = "CellBarcode_Identity",
    all.x = TRUE, sort = FALSE)
if (anyNA(umap$Disease_Identity) || anyNA(umap$Manuscript_Identity) || any(!is.finite(umap$UMAP_1)) ||
    any(!is.finite(umap$UMAP_2))) stop("Invalid embedding metadata")
write_table(umap, "umap_plot_data.csv")
