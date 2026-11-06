# Analyse directional pathway overrepresentation and ranked enrichment.  Uses the cancer-specific
# tested universe, signed Cox rankings and fixed MSigDB package version used by the completed
# analysis.  Run with: Rscript 07_pathway_analysis.R

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

if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr 26.1.1 before running this script.")
if (as.character(utils::packageVersion("msigdbr")) != "26.1.1") stop("This analysis requires msigdbr 26.1.1 to preserve the validated gene-set definitions.")


# ---------------------------------------------------------------------------- Directional pathway
# overrepresentation ----------------------------------------------------------------------------

options(width = 180, scipen = 5)
root <- PROJECT_ROOT
in_dir <- file.path(root, "Results_Rebuild", "Stage05")
out_dir <- file.path(root, "Results_Rebuild", "Stage11A")
fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "ggplot2", "msigdbr", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing package ", pkg, ". Install with install.packages(\"", pkg, "\"). ", "For msigdbr, internet may be needed on first installation.")
    }
}
suppressPackageStartupMessages(library(data.table))
ccs <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
paths <- file.path(in_dir, paste0("Stage05_", ccs, "_all_gene_Cox.csv"))
if (any(!file.exists(paths))) stop("Missing Stage05 files: ", paste(basename(paths[!file.exists(paths)]),
    collapse = ", "))

message("Stage11A 1/5: Load MSigDB Hallmark and Reactome gene-set definitions")
ms <- tryCatch(data.table::as.data.table(msigdbr::msigdbr(species = "Homo sapiens")), error = function(e) stop("MSigDB retrieval failed: ",
    conditionMessage(e), ". Install/update msigdbr and ensure gene-set data are available."))
collection <- if ("gs_collection" %in% names(ms)) "gs_collection" else "gs_cat"
subcoll <- if ("gs_subcollection" %in% names(ms)) "gs_subcollection" else "gs_subcat"
if (!all(c("gs_name", "gene_symbol", collection, subcoll) %in% names(ms))) stop("Unknown msigdbr schema: ",
    paste(names(ms), collapse = ", "))
ms[, Library := fifelse(get(collection) == "H", "Hallmark", fifelse(get(collection) == "C2" & grepl("REACTOME",
    get(subcoll), ignore.case = TRUE), "Reactome", NA_character_))]
ms <- unique(ms[!is.na(Library) & !is.na(gene_symbol) & nzchar(gene_symbol), .(Library, Pathway = gs_name,
    Gene = gene_symbol)])
if (!nrow(ms)) stop("No Hallmark or Reactome pathways found. Check msigdbr release/column names.")
gene_sets <- split(ms$Gene, paste(ms$Library, ms$Pathway, sep = "|||"))
set_metadata <- unique(ms[, .(Library, Pathway)])
setkey(set_metadata, Library, Pathway)

message("Stage11A 2/5: Directional over-representation against per-cancer tested genes")
# The per-cancer tested feature universe, not all annotated human genes, is used.  Exclude
# Entrez-like and deprecated non-symbol identifiers conservatively.
is_symbol <- function(x) !is.na(x) & grepl("^[A-Za-z][A-Za-z0-9._-]*$", x) & !grepl("^[0-9]+$", x)
# Main selection fixed a priori: top <=300 FDR<0.05 genes per direction.  This prevents huge
# KIRC/PAAD lists from dominating solely via large selection sets.  Also run top100 / top500 as
# explicitly labeled threshold sensitivity; never use sensitivity to select a more favourable
# threshold.
max_list_sizes <- c(100L, 300L, 500L)
min_set_size <- 10L
max_set_size <- 500L
analyses <- list()
cohort <- list()
selections <- list()
for (cc in ccs) {
    message("  ", cc)
    x <- fread(file.path(in_dir, paste0("Stage05_", cc, "_all_gene_Cox.csv")))
    req <- c("Feature_ID", "Beta", "P", "FDR", "Fit_status")
    if (!all(req %in% names(x)))
        stop("Stage05 schema mismatch for ", cc)
    x <- x[Fit_status == "ok" & is.finite(Beta) & is.finite(P) & is.finite(FDR) & is_symbol(Feature_ID)]
    setorder(x, P, Feature_ID)
    if (anyDuplicated(x$Feature_ID)) {
        warning(cc, ": duplicated feature symbols; keep lowest-P per symbol, report duplicates.")
        x <- unique(x, by = "Feature_ID")
    }
    universe <- unique(x$Feature_ID)
    cohort[[cc]] <- data.table(Cancer = cc, Tested_symbols = length(universe), Discovery_FDR_005 = sum(x$FDR <
        0.05), Adverse_FDR_005 = sum(x$FDR < 0.05 & x$Beta > 0), Favourable_FDR_005 = sum(x$FDR < 0.05 &
        x$Beta < 0))
    for (topn in max_list_sizes) for (direction in c("Adverse", "Favourable")) {
        z <- x[FDR < 0.05 & (if (direction == "Adverse")
            Beta > 0 else Beta < 0)]
        setorder(z, P, Feature_ID)
        query <- head(z$Feature_ID, topn)
        selections[[length(selections) + 1L]] <- data.table(Cancer = cc, Direction = direction, Maximum_genes = topn,
            Gene = query)
        K <- length(query)
        N <- length(universe)
        if (K < 5)
            next
        for (key in names(gene_sets)) {
            present <- intersect(gene_sets[[key]], universe)
            M <- length(present)
            if (M < min_set_size || M > max_set_size)
                next
            hit <- intersect(query, present)
            k <- length(hit)
            # p(X>=k), universe is the valid tested human gene symbols.
            p <- stats::phyper(k - 1, M, N - M, K, lower.tail = FALSE)
            pair <- strsplit(key, "|||", fixed = TRUE)[[1]]
            analyses[[length(analyses) + 1L]] <- data.table(Cancer = cc, Direction = direction, Maximum_genes = topn,
                Library = pair[1], Pathway = pair[2], Universe_N = N, Query_K = K, Pathway_in_universe_M = M,
                Overlap_k = k, Expected_overlap = K * M/N, Fold_enrichment = if (M > 0)
                  k/(K * M/N) else NA_real_, P = p, Overlapping_genes = paste(sort(hit), collapse = ";"))
        }
    }
}
cohort_dt <- rbindlist(cohort, fill = TRUE)
res <- if (length(analyses)) rbindlist(analyses, fill = TRUE) else data.table()
if (!nrow(res)) stop("No enrichment tests: likely no FDR-significant candidates. Review cohort summary first.")
# Separate BH correction within cancer, direction, threshold and library.
res[, FDR := p.adjust(P, method = "BH"), by = .(Cancer, Direction, Maximum_genes, Library)]
setorder(res, Maximum_genes, FDR, P)
fwrite(res, file.path(out_dir, "Stage11A_all_ORA_tests.csv"))
fwrite(cohort_dt, file.path(out_dir, "Stage11A_discovery_universe_summary.csv"))
fwrite(rbindlist(selections, fill = TRUE), file.path(out_dir, "Stage11A_selected_gene_lists.csv"))
main <- res[Maximum_genes == 300L]
main_sig <- main[FDR < 0.05 & Overlap_k >= 3L]
fwrite(main_sig, file.path(out_dir, "Stage11A_main_significant_pathways.csv"))

message("Stage11A 3/5: Threshold consistency (100/300/500; exploratory)")
sensitivity <- res[FDR < 0.05 & Overlap_k >= 3, .(Significant_pathways = uniqueN(Pathway)), by = .(Cancer,
    Direction, Library, Maximum_genes)]
all_grid <- CJ(Cancer = ccs, Direction = c("Adverse", "Favourable"), Library = c("Hallmark", "Reactome"),
    Maximum_genes = max_list_sizes)
sensitivity <- merge(all_grid, sensitivity, by = c("Cancer", "Direction", "Library", "Maximum_genes"),
    all.x = TRUE)
sensitivity[is.na(Significant_pathways), Significant_pathways := 0L]
fwrite(sensitivity, file.path(out_dir, "Stage11A_threshold_sensitivity.csv"))

message("Stage11A 4/5: Standalone scientific figures")
figures <- character()
save_plot <- function(g, name, width = 12, height = 8) {
    f <- file.path(fig_dir, name)
    ragg::agg_webp(filename = f, width = width, height = height, units = "in", res = 300, quality = 95)
    tryCatch(print(g), finally = grDevices::dev.off())
    figures <<- c(figures, f)
    invisible(f)
}
plot_dt <- main_sig[Library == "Hallmark" & is.finite(Fold_enrichment)]
if (nrow(plot_dt)) {
    # At most five pathways per cancer and direction (no duplication inflation).
    setorder(plot_dt, Cancer, Direction, FDR, -Fold_enrichment)
    plot_dt <- plot_dt[, head(.SD, 5L), by = .(Cancer, Direction)]
    plot_dt[, Label := gsub("^HALLMARK_", "", Pathway)]
    plot_dt[, Label := gsub("_", " ", Label, fixed = TRUE)]
    plot_dt[, Panel := paste(Cancer, Direction, sep = " | ")]
    g <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = Panel, y = Label)) + ggplot2::geom_point(ggplot2::aes(size = Overlap_k,
        colour = pmin(-log10(pmax(FDR, 1e-300)), 35))) + ggplot2::scale_colour_viridis_c(name = "\342\210\222log10(FDR)",
        option = "C") + ggplot2::scale_size_continuous(name = "Overlapping genes", range = c(2, 9)) +
        ggplot2::labs(title = "Hallmark pathways enriched among survival-associated genes", subtitle = "Directional unadjusted discovery; top 300 FDR<0.05 genes per cancer and direction",
            x = NULL, y = NULL, caption = "BH correction per cancer \303\227 direction \303\227 library | exploratory ORA") +
        ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(text = ggplot2::element_text(colour = "black"),
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, colour = "black"), axis.text.y = ggplot2::element_text(colour = "black"),
        plot.title = ggplot2::element_text(face = "bold"), panel.grid.major.y = ggplot2::element_blank())
    save_plot(g, "Stage11A_Hallmark_directional_ORA_bubble.webp", width = 19, height = 12)
}
count_dt <- main_sig[, .(Significant_pathways = .N), by = .(Cancer, Direction, Library)]
count_grid <- CJ(Cancer = ccs, Direction = c("Adverse", "Favourable"), Library = c("Hallmark", "Reactome"))
count_dt <- merge(count_grid, count_dt, by = c("Cancer", "Direction", "Library"), all.x = TRUE)
count_dt[is.na(Significant_pathways), Significant_pathways := 0L]
g2 <- ggplot2::ggplot(count_dt, ggplot2::aes(x = Cancer, y = Significant_pathways, fill = Direction)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) + ggplot2::facet_wrap(~Library, scales = "free_y") +
    ggplot2::labs(title = "Discovery pathway-enrichment yield by cancer", subtitle = "Count of BH-significant pathways; top 300 direction-specific genes (not independent biology)",
        x = NULL, y = "Significant pathways (FDR < 0.05; overlap \342\211\245 3)") + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(text = ggplot2::element_text(colour = "black"), plot.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(colour = "black"), axis.text.y = ggplot2::element_text(colour = "black"),
        legend.position = "top")
save_plot(g2, "Stage11A_pathway_yield_by_cancer.webp", width = 13, height = 7)

message("Stage11A 5/5: Consolidated checkpoint")
key_genes <- c("SERPINE1", "MET", "STX1A", "PGK1", "DKK1", "PSMD2", "PLK1", "FLNC", "LY6D", "UTY")
key_pathways <- main_sig[grepl(paste(key_genes, collapse = "|"), Overlapping_genes)]
setorder(key_pathways, FDR)
fwrite(key_pathways, file.path(out_dir, "Stage11A_key_gene_containing_pathways.csv"))
manifest <- data.table(Figure_file = figures)
fwrite(manifest, file.path(out_dir, "Stage11A_figures.csv"))
fmt <- function(dt, n = 15) if (nrow(dt)) capture.output(print(dt[seq_len(min(n, .N))])) else "No rows."
report <- c("===============================================================", "PROJECT 3 | STAGE 11A \342\200\224 SURVIVAL-ASSOCIATED PATHWAY ENRICHMENT",
    "===============================================================", "1. COHORT DISCOVERY SUMMARIES",
    fmt(cohort_dt, 10), "2. MAIN SIGNIFICANT PATHWAYS BY CANCER / DIRECTION / LIBRARY", fmt(count_dt[order(Cancer,
        Library, Direction)], 40), "3. TOP ENRICHMENTS (FDR<0.05; overlap>=3)", fmt(main_sig[order(FDR),
        .(Cancer, Direction, Library, Pathway, Query_K, Overlap_k, Fold_enrichment, FDR)], 25), "4. FOCAL GENE-CONTAINING ENRICHMENTS",
    fmt(key_pathways[, .(Cancer, Direction, Library, Pathway, FDR, Overlapping_genes)], 20), "5. THRESHOLD SENSITIVITY",
    fmt(sensitivity, 35), "6. FIGURES", if (length(figures)) figures else "No Hallmark figure (no significant Hallmark pathways)",
    "7. SCIENTIFIC INTERPRETATION RULES", "* Cox associations are UNADJUSTED and discovered in TCGA; ORA is hypothesis-generating.",
    "* ORA uses the successfully tested, identifiable per-cancer gene universe.", "* Selection is capped at top 300 FDR<0.05 genes per direction; 100/500 are sensitivities.",
    "* Gene-set overlaps, correlation, model confounding and gene-set redundancy prevent causal interpretation.",
    "* No enrichment in a cancer can indicate limited statistical power, not biological absence.", "* No independent pathway validation, no protein-level confirmation, no prediction claim.",
    "RESULTS: ", out_dir, "===============================================================")
writeLines(report, file.path(out_dir, "Stage11A_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Rank-based gene-set
# enrichment ----------------------------------------------------------------------------

options(width = 180, scipen = 4)
root <- PROJECT_ROOT
input <- file.path(root, "Results_Rebuild", "Stage05")
output <- file.path(root, "Results_Rebuild", "Stage11B")
figure_dir <- file.path(output, "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
needed <- c("data.table", "ggplot2", "fgsea", "msigdbr", "ragg")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), ". Install CRAN packages with install.packages(c(\"data.table\",\"ggplot2\",\"msigdbr\",\"ragg\")); ",
    "install fgsea with BiocManager::install(\"fgsea\", ask=FALSE, update=FALSE).")
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
files <- setNames(file.path(input, paste0("Stage05_", cancers, "_all_gene_Cox.csv")), cancers)
if (any(!file.exists(files))) stop("Missing Stage05 files: ", paste(basename(files[!file.exists(files)]),
    collapse = ", "))

message("11B 1/5 \342\200\224 Load fixed MSigDB Hallmark and Reactome libraries")
ms <- tryCatch(as.data.table(msigdbr::msigdbr(species = "Homo sapiens")), error = function(e) stop("msigdbr failure: ",
    conditionMessage(e)))
coll <- if ("gs_collection" %in% names(ms)) "gs_collection" else "gs_cat"
subcoll <- if ("gs_subcollection" %in% names(ms)) "gs_subcollection" else "gs_subcat"
if (!all(c(coll, subcoll, "gs_name", "gene_symbol") %in% names(ms))) stop("Unsupported msigdbr schema")
ms[, Library := fifelse(get(coll) == "H", "Hallmark", fifelse(get(coll) == "C2" & grepl("REACTOME",
    get(subcoll), ignore.case = TRUE), "Reactome", NA_character_))]
ms <- unique(ms[!is.na(Library) & !is.na(gene_symbol) & nzchar(gene_symbol), .(Library, Pathway = gs_name,
    Gene = gene_symbol)])
if (!nrow(ms)) stop("No matching MSigDB pathways: verify library version")
sets <- lapply(split(ms, ms$Library), function(dt) split(dt$Gene, dt$Pathway))
lib_info <- data.table(Library = names(sets), N_pathways = vapply(sets, length, integer(1)), msigdbr_package_version = as.character(utils::packageVersion("msigdbr")))
fwrite(lib_info, file.path(output, "Stage11B_gene_set_inventory.csv"))

# Cox P is two-sided. Rank by signed Z; maintain the same adverse/favourable interpretation used in
# Stage 11A. Never rank by FDR (FDR induces ties).  Restrict to successful models and
# gene-symbol-like IDs. This does NOT convert unadjusted Cox results to clinical-adjusted
# estimates.
is_symbol <- function(s) !is.na(s) & grepl("^[A-Za-z][A-Za-z0-9._-]*$", s) & !grepl("^[0-9]+$", s)
all_results <- list()
summary_rows <- list()
problems <- list()
rank_catalog <- list()
message("11B 2/5 \342\200\224 Genome-wide signed Z ranking and fgseaMultilevel")
for (cc in cancers) {
    message("  ", cc)
    dt <- fread(files[[cc]])
    required <- c("Feature_ID", "Beta", "P", "Fit_status")
    if (!all(required %in% names(dt)))
        stop(cc, ": missing ", paste(setdiff(required, names(dt)), collapse = ", "))
    total_rows <- nrow(dt)
    dt <- dt[Fit_status == "ok" & is.finite(Beta) & is.finite(P) & P >= 0 & P <= 1 & is_symbol(Feature_ID) &
        Beta != 0]
    # P=0 is possible due to numerical underflow; cap at machine minimum.
    dt[, Signed_Z := sign(Beta) * stats::qnorm(pmax(P, .Machine$double.xmin)/2, lower.tail = FALSE)]
    dt <- dt[is.finite(Signed_Z)]
    setorder(dt, Feature_ID)
    dt <- dt[order(-abs(Signed_Z), Feature_ID)]
    duplicates <- sum(duplicated(dt$Feature_ID))
    dt <- unique(dt, by = "Feature_ID")
    rank <- dt$Signed_Z
    names(rank) <- dt$Feature_ID
    rank <- sort(rank, decreasing = TRUE)
    fwrite(dt[, .(Feature_ID, Beta, P, Signed_Z)], file.path(output, paste0("Stage11B_", cc, "_signed_gene_ranking.csv")))
    rank_catalog[[cc]] <- rank
    sr <- data.table(Cancer = cc, Stage05_rows = total_rows, Ranked_genes = length(rank), Duplicate_symbols_collapsed = duplicates,
        Positive_ranked = sum(rank > 0), Negative_ranked = sum(rank < 0), Zero_capped_P = sum(dt$P ==
            0), Median_abs_Z = median(abs(rank)), Status = "OK")
    if (length(rank) < 1000) {
        problems[[length(problems) + 1L]] <- data.table(Cancer = cc, Problem = "Too few ranking features")
        sr[, Status := "FAILED"]
        summary_rows[[cc]] <- sr
        next
    }
    summary_rows[[cc]] <- sr
    for (lib in names(sets)) {
        valid_sets <- lapply(sets[[lib]], function(g) intersect(g, names(rank)))
        valid_sets <- valid_sets[lengths(valid_sets) >= 15L & lengths(valid_sets) <= 500L]
        if (!length(valid_sets))
            next
        ans <- tryCatch({
            # Stable independent seed per cancer+library, no post-hoc selection.
            set.seed(20261102L + match(cc, cancers) * 100L + match(lib, names(sets)))
            suppressWarnings(fgsea::fgseaMultilevel(pathways = valid_sets, stats = rank, minSize = 15,
                maxSize = 500, eps = 1e-30, scoreType = "std", nproc = 1))
        }, error = function(e) {
            problems[[length(problems) + 1L]] <<- data.table(Cancer = cc, Problem = paste(lib, conditionMessage(e)))
            NULL
        })
        if (is.null(ans) || !nrow(ans))
            next
        ans <- as.data.table(ans)
        ans[, `:=`(Cancer = cc, Library = lib, Direction = fifelse(NES > 0, "Adverse", "Favourable"),
            Leading_edge = paste(vapply(leadingEdge, function(v) paste(v, collapse = ";"), character(1))))]
        ans[, leadingEdge := NULL]
        # fgsea padj is BH-adjusted within each cancer/library, across both NES directions.
        setnames(ans, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size"), c("Pathway", "P",
            "FDR", "Log2Error", "ES", "NES", "Set_size"), skip_absent = TRUE)
        all_results[[length(all_results) + 1L]] <- ans[, .(Cancer, Library, Pathway, Direction, Set_size,
            ES, NES, P, FDR, Log2Error, Leading_edge)]
    }
}
summary_dt <- rbindlist(summary_rows, fill = TRUE)
fwrite(summary_dt, file.path(output, "Stage11B_ranked_universe_diagnostics.csv"))
problems_dt <- if (length(problems)) rbindlist(problems, fill = TRUE) else data.table(Cancer = character(),
    Problem = character())
fwrite(problems_dt, file.path(output, "Stage11B_problems.csv"))
if (!length(all_results)) stop("All GSEA runs failed; see Stage11B_problems.csv")
res <- rbindlist(all_results, fill = TRUE)
res <- res[order(Cancer, Library, FDR, -abs(NES))]
fwrite(res, file.path(output, "Stage11B_all_ranked_pathway_tests.csv"))
sig <- res[is.finite(FDR) & FDR < 0.05]
fwrite(sig, file.path(output, "Stage11B_significant_ranked_pathways.csv"))

message("11B 3/5 \342\200\224 Compare GSEA patterns with Stage 11A without treating as replication")
ora_path <- file.path(root, "Results_Rebuild", "Stage11A", "Stage11A_main_significant_pathways.csv")
if (file.exists(ora_path)) {
    ora <- fread(ora_path)
    ora <- unique(ora[, .(Cancer, Library, Pathway, ORA_Direction = Direction, ORA_FDR = FDR)])
    concord <- merge(sig[, .(Cancer, Library, Pathway, GSEA_Direction = Direction, GSEA_NES = NES, GSEA_FDR = FDR)],
        ora, by = c("Cancer", "Library", "Pathway"), all = FALSE)
    concord[, Direction_agrees := GSEA_Direction == ORA_Direction]
} else {
    concord <- data.table(Cancer = character(), Library = character(), Pathway = character(), GSEA_Direction = character(),
        GSEA_NES = numeric(), GSEA_FDR = numeric(), ORA_Direction = character(), ORA_FDR = numeric(),
        Direction_agrees = logical())
    problems_dt <- rbind(problems_dt, data.table(Cancer = "ALL", Problem = "Stage11A main ORA absent; overlap not calculated"))
}
fwrite(concord, file.path(output, "Stage11B_overlap_with_Stage11A_ORA.csv"))
counts <- res[, .(Tested_pathways = .N, Significant_pathways = sum(FDR < 0.05, na.rm = TRUE)), by = .(Cancer,
    Library, Direction)]
all_grid <- CJ(Cancer = cancers, Library = names(sets), Direction = c("Adverse", "Favourable"))
counts <- merge(all_grid, counts, by = c("Cancer", "Library", "Direction"), all.x = TRUE)
counts[is.na(Significant_pathways), Significant_pathways := 0L]
counts[is.na(Tested_pathways), Tested_pathways := 0L]
fwrite(counts, file.path(output, "Stage11B_pathway_yield_by_cancer.csv"))

message("11B 4/5 \342\200\224 Publication-oriented standalone figures")
figures <- character()
save_webp <- function(plot, filename, w = 13, h = 8) {
    path <- file.path(figure_dir, filename)
    ragg::agg_webp(filename = path, width = w, height = h, units = "in", res = 300, quality = 95)
    tryCatch(print(plot), finally = grDevices::dev.off())
    figures <<- c(figures, path)
}
base_theme <- ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(text = ggplot2::element_text(colour = "black"),
    plot.title = ggplot2::element_text(face = "bold", colour = "black"), axis.text.x = ggplot2::element_text(colour = "black"),
    axis.text.y = ggplot2::element_text(colour = "black"), panel.grid.major.y = ggplot2::element_blank())
# Fixed biologically interpretable Hallmarks, not opportunistically selected genes.
hallmarks <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_HYPOXIA",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS", "HALLMARK_APOPTOSIS", "HALLMARK_ANGIOGENESIS",
    "HALLMARK_MTORC1_SIGNALING")
heat <- res[Library == "Hallmark" & Pathway %in% hallmarks]
if (nrow(heat)) {
    heat[, Pathway_label := factor(gsub("_", " ", sub("HALLMARK_", "", Pathway, fixed = TRUE), fixed = TRUE),
        levels = rev(gsub("_", " ", sub("HALLMARK_", "", hallmarks, fixed = TRUE), fixed = TRUE)))]
    heat[, Cancer := factor(Cancer, levels = cancers)]
    g <- ggplot2::ggplot(heat, ggplot2::aes(x = Cancer, y = Pathway_label)) + ggplot2::geom_tile(ggplot2::aes(fill = NES),
        colour = "white", linewidth = 0.6) + ggplot2::geom_point(data = heat[FDR < 0.05], shape = 21,
        fill = "black", colour = "white", size = 2.9, stroke = 0.5) + ggplot2::scale_fill_gradient2(low = "#2b6cb0",
        mid = "#f8fafc", high = "#b83236", midpoint = 0, name = "Normalised enrichment\nscore (NES)",
        na.value = "grey93") + ggplot2::labs(title = "Ranked survival associations reveal cancer-dependent pathways",
        subtitle = "Signed Cox Z ranking; positive NES = adverse-associated genes; dot = FDR < 0.05",
        x = NULL, y = NULL, caption = "Exploratory unadjusted TCGA associations | BH within cancer and library") +
        base_theme + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
    save_webp(g, "Stage11B_Hallmark_directional_NES_heatmap.webp", w = 15, h = 9)
}
g2 <- ggplot2::ggplot(counts, ggplot2::aes(x = Cancer, y = Significant_pathways, fill = Direction)) +
    ggplot2::geom_col(position = "dodge", width = 0.72) + ggplot2::facet_wrap(~Library, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(Adverse = "#b83236", Favourable = "#2b6cb0")) + ggplot2::labs(title = "Directional pathway enrichment by cancer",
    subtitle = "Genome-wide signed-Z GSEA; overlapping gene sets are not independent findings", x = NULL,
    y = "Pathways with FDR < 0.05") + base_theme + ggplot2::theme(legend.position = "top")
save_webp(g2, "Stage11B_ranked_GSEA_pathway_yield.webp", w = 13, h = 7)
# A separate bubble figure for the most significant Reactome programs, when present.
react <- sig[Library == "Reactome" & is.finite(NES)]
if (nrow(react)) {
    react <- react[order(Cancer, FDR, -abs(NES))]
    react <- react[, head(.SD, 3L), by = Cancer]
    react[, Pathway_label := gsub("_", " ", sub("REACTOME_", "", Pathway, fixed = TRUE), fixed = TRUE)]
    react[, Pathway_label := ifelse(nchar(Pathway_label) > 62, paste0(substr(Pathway_label, 1, 59),
        "..."), Pathway_label)]
    react[, Panel := paste(Cancer, Direction, sep = " | ")]
    g3 <- ggplot2::ggplot(react, ggplot2::aes(x = Panel, y = Pathway_label)) + ggplot2::geom_point(ggplot2::aes(size = abs(NES),
        colour = pmin(-log10(pmax(FDR, 1e-300)), 40))) + ggplot2::scale_colour_viridis_c(name = "\342\210\222log10(FDR)",
        option = "C") + ggplot2::scale_size_continuous(name = "Absolute NES", range = c(2.2, 8)) + ggplot2::labs(title = "Reactome pathways identified by ranked GSEA",
        subtitle = "Top 3 FDR-ranked pathways per cancer; descriptive selected view", x = NULL, y = NULL) +
        base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, colour = "black"))
    save_webp(g3, "Stage11B_Reactome_ranked_enrichment_bubbles.webp", w = 19, h = 12)
}
fwrite(data.table(Figure = figures), file.path(output, "Stage11B_figure_manifest.csv"))

message("11B 5/5 \342\200\224 Consolidated scientific checkpoint (all report output at end)")
fmt <- function(z, n = 25L) if (nrow(z)) capture.output(print(head(z, n))) else "None"
top <- sig[order(FDR, -abs(NES)), .(Cancer, Library, Direction, Pathway, NES, FDR, Set_size)][1:min(25,
    .N)]
recurrent <- sig[, .(N_cancers = uniqueN(Cancer), Cancers = paste(sort(unique(Cancer)), collapse = ";")),
    by = .(Library, Pathway, Direction)][N_cancers >= 2][order(-N_cancers, Library, Pathway)]
fwrite(recurrent, file.path(output, "Stage11B_recurrent_significant_pathways.csv"))
report <- c("================================================================", "PROJECT 3 | STAGE 11B \342\200\224 GENOME-WIDE RANKED PATHWAY ENRICHMENT",
    "================================================================", "1. RANKING UNIVERSE / DIAGNOSTICS",
    fmt(summary_dt, 10), "2. PATHWAY ENRICHMENT COUNTS BY CANCER, LIBRARY, DIRECTION", fmt(counts, 40),
    "3. TOP SIGNIFICANT ENRICHMENTS", fmt(top, 25), "4. RECURRENT SIGNIFICANT PATHWAYS", fmt(recurrent,
        25), "5. OVERLAP WITH STAGE 11A ORA (NOT INDEPENDENT REPLICATION)", paste0("Matching significant cancer-library-pathway results: ",
        nrow(concord)), paste0("Same sign: ", sum(concord$Direction_agrees, na.rm = TRUE), " of ", nrow(concord)),
    fmt(concord[order(GSEA_FDR)], 12), "6. PROBLEMS", fmt(problems_dt, 25), "7. FIGURES", if (length(figures)) figures else "None",
    "8. SCIENTIFIC CAUTIONS", "* ALL rank statistics derive from UNADJUSTED cancer-specific Cox discovery; no biological causation asserted.",
    "* Positive NES: positive Cox log-HR rank end (adverse). Negative NES: favourable rank end.", "* Rank-based GSEA avoids Stage11A top-300 selection but shares the SAME TCGA cohort, genes and biases.",
    "* Gene\342\200\223gene correlation, tumour cell composition, sex, clinical covariates and pathway redundancy may confound interpretations.",
    "* BH correction performed within each cancer and gene-set library across both enrichment directions.",
    "* No external protein/survival validation; negative enrichment findings are not evidence of biological absence.",
    "RESULTS:", output, "================================================================")
writeLines(report, file.path(output, "Stage11B_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
