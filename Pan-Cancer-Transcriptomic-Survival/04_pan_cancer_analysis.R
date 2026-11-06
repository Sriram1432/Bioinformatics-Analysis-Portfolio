# Integrate cancer-specific effects and quantify cross-cancer heterogeneity.  Combines stratified
# Cox estimates, recurrence summaries and heterogeneity statistics while retaining cancer-specific
# effect directions.  Run with: Rscript 04_pan_cancer_analysis.R

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


# ---------------------------------------------------------------------------- Pan-cancer
# integration and heterogeneity
# ----------------------------------------------------------------------------

options(scipen = 999, width = 190)
root <- PROJECT_ROOT
in_dir <- file.path(root, "Results_Rebuild", "Stage03")
stage05 <- file.path(root, "Results_Rebuild", "Stage05")
out_dir <- file.path(root, "Results_Rebuild", "Stage06")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (p in c("data.table", "survival", "matrixStats", "ggplot2", "ragg")) {
    if (!requireNamespace(p, quietly = TRUE))
        stop("Install required package: ", p)
}
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")

# 1. Validate upstream inputs and establish the COMMON feature universe.
input_files <- setNames(file.path(in_dir, paste0("TCGA_", cancers, "_discovery.rds")), cancers)
prior_files <- setNames(file.path(stage05, paste0("Stage05_", cancers, "_all_gene_Cox.csv")), cancers)
if (any(!file.exists(input_files))) stop("Missing Stage03 RDS: ", paste(names(input_files)[!file.exists(input_files)],
    collapse = ", "))
if (any(!file.exists(prior_files))) stop("Missing Stage05 gene results: ", paste(names(prior_files)[!file.exists(prior_files)],
    collapse = ", "))
meta <- vector("list", length(cancers))
names(meta) <- cancers
shared <- NULL
for (cc in cancers) {
    d <- readRDS(input_files[[cc]])
    if (!is.list(d) || !all(c("expression", "clinical") %in% names(d)))
        stop("Stage03 schema mismatch: ", cc)
    x <- d$expression
    cl <- as.data.table(d$clinical)
    if (!is.matrix(x) || !is.numeric(x) || is.null(rownames(x)) || !all(c("sample_id", "patient_id",
        "OS", "OS_time") %in% names(cl)) || !identical(colnames(x), as.character(cl$sample_id)) || anyDuplicated(rownames(x)) ||
        anyDuplicated(cl$patient_id) || anyNA(cl$OS) || anyNA(cl$OS_time) || !all(cl$OS %in% c(0, 1)) ||
        !all(is.finite(cl$OS_time) & cl$OS_time > 0))
        stop("Cohort integrity failed: ", cc)
    shared <- if (is.null(shared))
        rownames(x) else intersect(shared, rownames(x))
    meta[[cc]] <- data.table(Cancer = cc, Patients = nrow(cl), Events = sum(cl$OS == 1L), Features = nrow(x))
    rm(d, x, cl)
}
cohorts <- rbindlist(meta)
if (sum(cohorts$Patients) != 4718L || sum(cohorts$Events) != 1519L || length(shared) != 17465L) stop("Upstream cohort/feature checkpoint changed; review before proceeding")

# 2. Load common genes only; standardise EACH gene ACROSS patients WITHIN each cancer.  This makes
# HR units: one within-cancer SD increase; not an inter-tumour absolute expression comparison and
# not additional batch correction.
zmat <- vector("list", length(cancers))
names(zmat) <- cancers
clinical_parts <- vector("list", length(cancers))
for (j in seq_along(cancers)) {
    cc <- cancers[j]
    message("Stage06: preparing ", cc, " (", j, "/10)")
    d <- readRDS(input_files[[cc]])
    x <- d$expression[shared, , drop = FALSE]
    cl <- as.data.table(d$clinical)
    if (anyNA(x) || any(!is.finite(x)))
        stop("Nonfinite shared expression: ", cc)
    mu <- matrixStats::rowMeans2(x)
    ss <- matrixStats::rowSds(x)
    if (any(!is.finite(ss) | ss < 1e-10))
        stop("Invariant/invalid common feature: ", cc)
    x <- sweep(x, 1L, mu, "-")
    x <- sweep(x, 1L, ss, "/")
    storage.mode(x) <- "double"
    if (!identical(colnames(x), as.character(cl$sample_id)))
        stop("Alignment mismatch: ", cc)
    zmat[[cc]] <- x
    clinical_parts[[j]] <- data.table(Cancer = cc, sample_id = cl$sample_id, OS_time = as.numeric(cl$OS_time),
        OS = as.integer(cl$OS))
    rm(d, x, cl, mu, ss)
    invisible(gc(verbose = FALSE))
}
pan_clin <- rbindlist(clinical_parts)
if (!identical(unlist(lapply(zmat, colnames), use.names = FALSE), pan_clin$sample_id)) stop("Pan-cancer specimen order mismatch")
pan_clin[, Cancer := factor(Cancer, levels = cancers)]
fit_data <- data.frame(time = pan_clin$OS_time, event = pan_clin$OS, Cancer = pan_clin$Cancer, z = numeric(nrow(pan_clin)))

# 3. Cancer-stratified Cox model, common association coefficient; separate baseline hazard per
# malignancy. Efron tie correction.
fit_one <- function(z) {
    if (any(!is.finite(z)))
        return(list(beta = NA_real_, se = NA_real_, p = NA_real_, status = "nonfinite_input"))
    fit_data$z <- z
    notes <- character()
    model <- withCallingHandlers(tryCatch(survival::coxph(survival::Surv(time, event) ~ z + strata(Cancer),
        data = fit_data, ties = "efron", control = survival::coxph.control(iter.max = 30), model = FALSE,
        x = FALSE, y = FALSE), error = function(e) structure(list(message = conditionMessage(e)), class = "cox_error")),
        warning = function(w) {
            notes <<- c(notes, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    if (inherits(model, "cox_error"))
        return(list(beta = NA_real_, se = NA_real_, p = NA_real_, status = "fit_error"))
    if (length(notes))
        return(list(beta = NA_real_, se = NA_real_, p = NA_real_, status = "fit_warning"))
    b <- unname(coef(model)[1L])
    V <- model$var[1L, 1L]
    if (!is.finite(b) || !is.finite(V) || V <= 0 || abs(b) > 15)
        return(list(beta = NA_real_, se = NA_real_, p = NA_real_, status = "unstable_estimate"))
    se <- sqrt(V)
    list(beta = b, se = se, p = 2 * pnorm(-abs(b/se)), status = "ok")
}

message("Stage06: fitting ", length(shared), " cancer-stratified Cox models...")
results <- vector("list", length(shared))
for (i in seq_along(shared)) {
    # Unlist preserves the exact cancer block order of pan_clin.
    zz <- unlist(lapply(zmat, function(x) x[i, ]), use.names = FALSE)
    results[[i]] <- fit_one(zz)
    if (i%%2000L == 0L)
        message("  completed ", i, "/", length(shared))
}
pan <- data.table(Feature_ID = shared, Beta = vapply(results, `[[`, numeric(1), "beta"), SE = vapply(results,
    `[[`, numeric(1), "se"), P = vapply(results, `[[`, numeric(1), "p"), Fit_status = vapply(results,
    `[[`, character(1), "status"))
pan[, `:=`(HR = exp(Beta), CI_low = exp(Beta - qnorm(0.975) * SE), CI_high = exp(Beta + qnorm(0.975) *
    SE))]
pan[, FDR := p.adjust(P, method = "BH", n = .N)]
pan[, Direction := fifelse(is.finite(HR) & HR > 1, "Adverse", fifelse(is.finite(HR) & HR < 1, "Favourable",
    NA_character_))]
setorder(pan, FDR, P, Feature_ID, na.last = TRUE)
fwrite(pan, file.path(out_dir, "Stage06_pan_cancer_stratified_Cox_all_features.csv"))
rm(zmat, results)
invisible(gc(verbose = FALSE))

# 4. Cross-cancer evidence from Stage05 (do NOT refit cancer-specific models).
prior <- rbindlist(lapply(cancers, function(cc) {
    x <- fread(prior_files[[cc]])
    required <- c("Cancer", "Feature_ID", "Beta", "SE", "HR", "P", "FDR", "Fit_status")
    if (!all(required %in% names(x)))
        stop("Stage05 schema changed: ", cc)
    if (!all(x$Cancer == cc) || anyDuplicated(x$Feature_ID))
        stop("Stage05 identifiers invalid: ", cc)
    x[Feature_ID %chin% shared, ..required]
}))
# Every shared identifier must have exactly one row in each cancer.
if (nrow(prior) != length(shared) * length(cancers)) stop("Incomplete Stage05 cross-cancer feature correspondence")

# Report recurrence among cancer-specific FDR discoveries; this alone is NOT a formal test that
# effects differ between cancers.
recurrence <- prior[FDR < 0.05 & is.finite(FDR), .(N_significant_cancers = .N, N_adverse = sum(HR > 1),
    N_favourable = sum(HR < 1), Significant_in = paste(Cancer, collapse = ", ")), by = Feature_ID]
setorder(recurrence, -N_significant_cancers, Feature_ID)
fwrite(recurrence, file.path(out_dir, "Stage06_cross_cancer_FDR_recurrence.csv"))

# 5. Exploratory inverse-variance effect meta-analysis and heterogeneity.  Cochran's Q / I2
# quantify dispersion; because samples are distinct across cancers, cancer-specific estimates are
# treated as independent. This is NOT an adjusted/random-effects estimate and does not validate the
# shared-effect assumption of the stratified Cox model.
good <- prior[Fit_status == "ok" & is.finite(Beta) & is.finite(SE) & SE > 0, .(Feature_ID, Cancer, Beta,
    SE, Weight = 1/SE^2)]
meta <- good[, {
    k <- .N
    if (k < 2L)
        .(K = k, Meta_HR = NA_real_, Meta_P = NA_real_, Q = NA_real_, I2_pct = NA_real_, Q_P = NA_real_) else {
        b <- sum(Weight * Beta)/sum(Weight)
        q <- sum(Weight * (Beta - b)^2)
        .(K = k, Meta_HR = exp(b), Meta_P = 2 * pnorm(-abs(b * sqrt(sum(Weight)))), Q = q, I2_pct = if (q >
            0)
            max(0, (q - (k - 1))/q) * 100 else 0, Q_P = pchisq(q, df = k - 1, lower.tail = FALSE))
    }
}, by = Feature_ID]
# Include absent/failed genes in each correction family.
meta <- merge(data.table(Feature_ID = shared), meta, by = "Feature_ID", all.x = TRUE, sort = FALSE)
meta[, `:=`(Q_FDR = p.adjust(Q_P, method = "BH", n = .N), Meta_FDR = p.adjust(Meta_P, method = "BH",
    n = .N))]
fwrite(meta, file.path(out_dir, "Stage06_cross_cancer_meta_heterogeneity.csv"))

joined <- merge(pan, meta, by = "Feature_ID", all.x = TRUE, sort = FALSE)
joined <- merge(joined, recurrence, by = "Feature_ID", all.x = TRUE, sort = FALSE)
setorder(joined, FDR, P, Feature_ID, na.last = TRUE)
fwrite(joined, file.path(out_dir, "Stage06_integrated_pan_cancer_evidence.csv"))


# 6. PRESENTATION-QUALITY FIGURES \342\200\224 independent, no lettered panels.  Statistical
# design: point AREA shows |log2(HR)|, point FILL shows -log10 adjusted FDR, and point OUTLINE
# shows HR direction. CI is shown in an accompanying forest plot, not incorrectly encoded as a 2D
# bubble attribute.  These are exploratory UNADJUSTED estimates and cannot establish causality.
for (p in c("ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing ", p)
fig_dir <- file.path(out_dir, "Figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
figure_paths <- character()
publication_theme <- function(base = 13) {
    ggplot2::theme_minimal(base_size = base) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
        colour = "black", size = 18, margin = ggplot2::margin(b = 8)), plot.subtitle = ggplot2::element_text(colour = "black",
        size = 11, margin = ggplot2::margin(b = 12)), plot.caption = ggplot2::element_text(colour = "black",
        size = 9, hjust = 0), axis.title = ggplot2::element_text(colour = "black", face = "bold"), axis.text = ggplot2::element_text(colour = "black",
        face = "bold"), strip.text = ggplot2::element_text(colour = "black", face = "bold"), legend.title = ggplot2::element_text(colour = "black",
        face = "bold"), panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = "#E9EDF1",
        linewidth = 0.25), plot.margin = ggplot2::margin(14, 22, 16, 14))
}
save_fig <- function(g, filename, w = 4200, h = 3000) {
    f <- file.path(fig_dir, filename)
    ragg::agg_webp(filename = f, width = w, height = h, res = 300, quality = 96)
    tryCatch(print(g), finally = grDevices::dev.off())
    figure_paths <<- c(figure_paths, f)
    invisible(f)
}
# Candidate selection is descriptive, not confirmatory. Force broadly informative candidates, not
# just those with the smallest pan-cancer FDR.  KIRC/PAAD findings alone cannot dominate the
# display.
joined[, Recurrence := data.table::fifelse(is.na(N_significant_cancers), 0L, as.integer(N_significant_cancers))]
joined[, Hetero_rank := data.table::fifelse(is.finite(I2_pct), I2_pct, 100)]
valid <- joined[Fit_status == "ok" & is.finite(FDR) & is.finite(HR) & HR > 0 & is.finite(CI_low) & CI_low >
    0]
if (nrow(valid) < 40) stop("Insufficient stable pan-cancer results to design heatmap")
# Select a deterministic balanced set: recurrence, homogeneous associations, heterogeneous
# associations, strongest pan results, and clear biologic exemplars.
choose <- function(x, n) head(unique(x), n)
setorder(valid, -Recurrence, FDR, Feature_ID)
selected <- choose(valid$Feature_ID, 16)
homo <- valid[I2_pct <= 25 & FDR < 0.05 & is.finite(I2_pct)][order(FDR)]$Feature_ID
hetero <- valid[I2_pct >= 75 & Q_FDR < 0.05 & is.finite(Q_FDR)][order(FDR)]$Feature_ID
selected <- unique(c(selected, choose(homo, 9), choose(hetero, 9)))
selected <- unique(c(selected, c("SERPINE1", "STX1A", "DKK1", "PGK1", "MET", "DONSON", "SRD5A1")))
selected <- selected[selected %in% valid$Feature_ID]
selected <- unique(c(selected, choose(valid[order(FDR)]$Feature_ID, 45)))
selected <- head(selected, 45L)
# Recheck Stage05 associations from long-format file for every matrix position.
sel_specific <- prior[Feature_ID %chin% selected, .(Feature_ID, Cancer, Beta, SE, HR, P, FDR, Fit_status)]
sel_pan <- pan[Feature_ID %chin% selected, .(Feature_ID, Cancer = "PAN-CANCER", Beta, SE, HR, P, FDR,
    Fit_status)]
bubble <- rbindlist(list(sel_specific, sel_pan), use.names = TRUE, fill = TRUE)
if (nrow(bubble) != length(selected) * 11L) stop("Master plot incomplete: selected gene/cancer positions missing")
bubble[, `:=`(effect_size = abs(log2(HR)), fdr_mlog10 = -log10(pmax(FDR, .Machine$double.xmin)), Direction = fifelse(is.finite(HR) &
    HR > 1, "HR > 1 (adverse)", fifelse(is.finite(HR) & HR < 1, "HR < 1 (favourable)", "Unavailable")))]
bubble[!is.finite(effect_size), effect_size := NA_real_]
bubble[!is.finite(fdr_mlog10), fdr_mlog10 := NA_real_]
bubble[, capped_colour := pmin(fdr_mlog10, 12)]
bubble[!is.finite(FDR) | FDR >= 0.05, capped_colour := 0]
# Cluster rows using signed log-HRs, with zero for failed fits ONLY for clustering (never as a
# measured effect); gene results preserve NA.
cluster_wide <- data.table::dcast(bubble[Cancer != "PAN-CANCER"], Feature_ID ~ Cancer, value.var = "Beta")
ord <- match(selected, cluster_wide$Feature_ID)
mat <- as.matrix(cluster_wide[ord, ..cancers])
mat[!is.finite(mat)] <- 0
# Per-cancer centered/scaled column patterns prevent one cancer dominating row order; fill-in zero
# applies ONLY to unsupervised display ordering.
if (nrow(mat) > 1L && ncol(mat) > 1L) {
    z <- scale(mat)
    z[!is.finite(z)] <- 0
    hc <- hclust(dist(z), method = "ward.D2")
    row_order <- selected[hc$order]
} else row_order <- selected
bubble[, Feature_ID := factor(Feature_ID, levels = rev(row_order))]
bubble[, Cancer := factor(Cancer, levels = c(cancers, "PAN-CANCER"))]
bubble[, Direction := factor(Direction, levels = c("HR > 1 (adverse)", "HR < 1 (favourable)", "Unavailable"))]
master <- ggplot2::ggplot(bubble, ggplot2::aes(x = Cancer, y = Feature_ID)) + ggplot2::geom_tile(fill = "#F7F9FB",
    colour = "white", linewidth = 0.35) + ggplot2::geom_point(ggplot2::aes(size = effect_size, fill = capped_colour,
    colour = Direction), shape = 21, stroke = 0.75, alpha = 0.96, na.rm = TRUE) + ggplot2::scale_size_continuous(name = "|log2(HR)|",
    range = c(0.7, 10), breaks = c(0.25, 0.5, 1, 1.5, 2), limits = c(0, max(2, max(bubble$effect_size,
        na.rm = TRUE)))) + ggplot2::scale_fill_gradientn(name = expression(-log[10]("FDR") ~ "(capped at 12)"),
    colours = c("#E4EAF0", "#84C9CE", "#2D7FA8", "#263C7B", "#8D2C65"), limits = c(0, 12), oob = scales::squish,
    breaks = c(0, 2, 4, 8, 12)) + ggplot2::scale_colour_manual(name = "Direction", values = c(`HR > 1 (adverse)` = "#A92A37",
    `HR < 1 (favourable)` = "#176F76", Unavailable = "#BAC1C9"), drop = FALSE) + ggplot2::labs(title = "Shared and cancer-specific prognostic associations",
    subtitle = "Selected 45 source genes | 10 cancer cohorts + cancer-stratified pan-cancer model | hierarchical gene ordering",
    x = NULL, y = NULL, caption = "Bubble area: magnitude of log2(HR), per 1 within-cancer SD expression increase. Outline: HR direction.\nFill: BH-FDR significance; non-significant positions set to zero colour. Exploratory, unadjusted associations; CIs in companion forest plot.") +
    publication_theme(13) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
    vjust = 1, size = 11), axis.text.y = ggplot2::element_text(size = 10), panel.grid = ggplot2::element_blank(),
    legend.position = "right", plot.caption = ggplot2::element_text(size = 9))
save_fig(master, "Stage06_MASTER_clustered_gene_by_cancer_pan_cancer_bubble.webp", w = 5900, h = 5100)
fwrite(bubble, file.path(out_dir, "Stage06_master_heatmap_long_data.csv"))
fwrite(data.table(Feature_ID = row_order, Gene_cluster_order = seq_along(row_order)), file.path(out_dir,
    "Stage06_master_heatmap_gene_order.csv"))

# Pan-cancer association with the exact 95% Wald CIs as an interpretable forest.
forest <- valid[order(FDR)][1:min(22L, .N), .(Feature_ID, HR, CI_low, CI_high, FDR, I2_pct)]
forest[, Feature_ID := factor(Feature_ID, levels = rev(Feature_ID))]
fp <- ggplot2::ggplot(forest, ggplot2::aes(x = HR, y = Feature_ID)) + ggplot2::geom_vline(xintercept = 1,
    linetype = "dashed", colour = "#697986") + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high),
    orientation = "y", width = 0.22, linewidth = 0.8, colour = "#315373") + ggplot2::geom_point(ggplot2::aes(colour = I2_pct),
    size = 3.8) + ggplot2::scale_x_log10() + ggplot2::scale_colour_gradient(low = "#297E8A", high = "#B63D54",
    na.value = "#93999F", name = "Heterogeneity I\302\262 (%)") + ggplot2::labs(title = "Pan-cancer prognostic effect sizes and uncertainty",
    subtitle = "Cancer-stratified Cox model | 22 strongest exploratory associations | 95% Wald confidence intervals",
    x = "Hazard ratio (logarithmic axis; 95% CI)", y = NULL, caption = "HR per 1 within-cancer SD expression increase. Shared coefficient, cancer-specific baseline hazards.\nColour indicates cross-cancer heterogeneity; clinical adjustment and PH diagnostics are pending.") +
    publication_theme()
save_fig(fp, "Stage06_pan_cancer_top22_forest_95CI.webp", w = 4200, h = 3600)

# Heterogeneity landscape: important competing explanation for pooled signal.
het_plot <- meta[is.finite(I2_pct) & is.finite(Meta_FDR) & is.finite(Q_FDR)]
if (nrow(het_plot)) {
    het_plot[, Discovery := factor(ifelse(Meta_FDR < 0.05, "Meta FDR < 0.05", "Meta FDR >= 0.05"),
        levels = c("Meta FDR < 0.05", "Meta FDR >= 0.05"))]
    hp <- ggplot2::ggplot(het_plot, ggplot2::aes(x = I2_pct, y = -log10(pmax(Q_FDR, .Machine$double.xmin)))) +
        ggplot2::geom_hex(bins = 48) + ggplot2::scale_fill_gradient(low = "#E3EDF2", high = "#284B70",
        name = "Gene count") + ggplot2::geom_vline(xintercept = c(25, 50, 75), linetype = "dashed", colour = "#8B98A5",
        linewidth = 0.3) + ggplot2::labs(title = "Effect heterogeneity across cancer types", subtitle = "Cochran Q BH-FDR and I\302\262 from univariable cancer-specific Cox estimates",
        x = "I\302\262 heterogeneity (%)", y = "\342\210\222log10(Q FDR)", caption = "Exploratory heterogeneity, not a causal or mechanism test. Horizontal interpretation depends on sample/event sizes.") +
        publication_theme()
    # geom_hex needs hexbin; scatter fallback when absent
    if (!requireNamespace("hexbin", quietly = TRUE)) {
        hp <- ggplot2::ggplot(het_plot, ggplot2::aes(x = I2_pct, y = -log10(pmax(Q_FDR, .Machine$double.xmin)))) +
            ggplot2::geom_bin2d(bins = 45) + ggplot2::scale_fill_gradient(low = "#E3EDF2", high = "#284B70",
            name = "Gene count") + ggplot2::geom_vline(xintercept = c(25, 50, 75), linetype = "dashed",
            colour = "#8B98A5", linewidth = 0.3) + ggplot2::labs(title = "Effect heterogeneity across cancer types",
            subtitle = "Cochran Q BH-FDR and I\302\262 from cancer-specific Cox estimates", x = "I\302\262 heterogeneity (%)",
            y = "\342\210\222log10(Q FDR)", caption = "Exploratory heterogeneity only; heterogeneity tests do not establish biological mechanism.") +
            publication_theme()
    }
    save_fig(hp, "Stage06_cross_cancer_effect_heterogeneity_density.webp", w = 4200, h = 3000)
}

# Strong vs weak cross-cancer recurrence, not a false impression of validation.
rr <- recurrence[, .N, by = N_significant_cancers][order(N_significant_cancers)]
if (nrow(rr)) {
    rr[, N_significant_cancers := factor(N_significant_cancers, levels = sort(unique(N_significant_cancers)))]
    rp <- ggplot2::ggplot(rr, ggplot2::aes(x = N_significant_cancers, y = N)) + ggplot2::geom_col(fill = "#257D91",
        width = 0.72) + ggplot2::geom_text(ggplot2::aes(label = format(N, big.mark = ",")), vjust = -0.5,
        fontface = "bold", size = 4) + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
        0.1))) + ggplot2::labs(title = "How often are prognostic associations reproduced across cancers?",
        subtitle = "Counts of genes with within-cancer BH-FDR < 0.05; not independent validation", x = "Number of cancers showing FDR < 0.05",
        y = "Number of source genes", caption = "Cancer type changes patient mix and tissue context; recurrence does not guarantee equal effect direction.") +
        publication_theme()
    save_fig(rp, "Stage06_cross_cancer_FDR_recurrence.webp", w = 3800, h = 2650)
}

# 7. SCIENTIFIC CHECKPOINT
summary_pan <- data.table(Patients = nrow(pan_clin), Events = sum(pan_clin$OS == 1), Cancers = length(cancers),
    Features_tested = nrow(pan), Successful_fits = sum(pan$Fit_status == "ok"), Failed_or_flagged = sum(pan$Fit_status !=
        "ok"), Nominal_P_005 = sum(pan$P < 0.05, na.rm = TRUE), FDR_005 = sum(pan$FDR < 0.05, na.rm = TRUE),
    FDR_adverse = sum(pan$FDR < 0.05 & pan$HR > 1, na.rm = TRUE), FDR_favourable = sum(pan$FDR < 0.05 &
        pan$HR < 1, na.rm = TRUE))
fit_status <- pan[, .N, by = Fit_status][order(-N)]
rec_summary <- recurrence[, .N, by = N_significant_cancers][order(-N_significant_cancers)]
hetero_summary <- meta[, .(Tested = sum(!is.na(Q_P)), Q_FDR_005 = sum(Q_FDR < 0.05, na.rm = TRUE), Median_I2_pct = median(I2_pct,
    na.rm = TRUE), All_10_valid = sum(K == 10L, na.rm = TRUE))]
top_pan <- joined[order(FDR, P)][1:min(20L, .N), .(Feature_ID, HR, CI_low, CI_high, P, FDR, K, I2_pct,
    Q_P, Q_FDR, N_significant_cancers)]
top_rec <- recurrence[1:min(20L, .N)]
for (nm in c("summary_pan", "fit_status", "rec_summary", "hetero_summary", "top_pan", "top_rec")) fwrite(get(nm),
    file.path(out_dir, paste0("Stage06_", nm, ".csv")))
report <- capture.output({
    cat("=====================================================================\n")
    cat("PROJECT 3 | STAGE 06 \342\200\224 PAN-CANCER SCIENTIFIC CHECKPOINT\n")
    cat("=====================================================================\n")
    cat("1. STRATIFIED COX DISCOVERY\n")
    print(summary_pan)
    cat("\n2. FIT STATUSES\n")
    print(fit_status)
    cat("\n3. TOP 20 PAN-CANCER GENES (95% CI; Q/I\302\262)\n")
    print(top_pan)
    cat("\n4. CROSS-CANCER FDR RECURRENCE\n")
    print(rec_summary)
    cat("\n5. RECURRENT FEATURE EXAMPLES\n")
    print(top_rec)
    cat("\n6. HETEROGENEITY SUMMARY\n")
    print(hetero_summary)
    cat("\n7. FIGURES GENERATED\n")
    print(figure_paths)
    cat("\n8. SCIENTIFIC LIMITATIONS\n")
    cat("Unadjusted Cox associations, not validated clinical biomarkers.\n")
    cat("One within-cancer SD expression increase; cancer-stratified baseline hazard.\n")
    cat("Shared pan-cancer coefficient does not imply homogeneous cancer effects.\n")
    cat("Q/I\302\262 uses unadjusted per-cancer beta/SE; effect heterogeneity may arise from case-mix.\n")
    cat("Gene selection for figures is outcome-informed and purely descriptive.\n")
    cat("Selection, PH assumption, influence and clinical confounding require later stages.\n")
    cat("BH includes failed/NA tests in total feature family through explicit n.\n")
    cat("95% CI is represented in forest plot, not implied by bubble size or colour.\n")
    cat("Master heatmap cells with unavailable/failed estimates are not interpreted.\n")
    cat("Upstream source gene annotations remain historical; no silent gene remapping.\n")
    cat("\nOUTPUT: ", out_dir, "\n", sep = "")
    cat("=====================================================================\n")
})
writeLines(report, file.path(out_dir, "Stage06_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
