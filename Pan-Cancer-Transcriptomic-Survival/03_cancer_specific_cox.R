# Fit cancer-specific genome-wide Cox proportional-hazards models.  Expression is standardised
# within cancer and false-discovery rates are controlled within each cancer-specific testing
# family.  Run with: Rscript 03_cancer_specific_cox.R

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


# ---------------------------------------------------------------------------- Genome-wide
# cancer-specific Cox modelling
# ----------------------------------------------------------------------------

options(scipen = 999, width = 180)
root <- PROJECT_ROOT
input_dir <- file.path(root, "Results_Rebuild", "Stage03")
output_dir <- file.path(root, "Results_Rebuild", "Stage05")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (p in c("survival", "data.table", "matrixStats", "ggplot2", "ragg")) {
    if (!requireNamespace(p, quietly = TRUE))
        stop("Missing package: ", p, ". Install once using install.packages('", p, "').")
}
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")

# A failed/unstable fit remains in the multiple-testing family, with P = NA.  Do not silently
# interpret non-convergent or separated models as discoveries.
empty_fit <- function(reason) {
    list(beta = NA_real_, SE = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, P = NA_real_,
        status = reason)
}

fit_one <- function(z, fit_data) {
    if (any(!is.finite(z)) || !is.finite(stats::sd(z)) || stats::sd(z) < 1e-10)
        return(empty_fit("zero_or_near_zero_variance"))
    fit_data$z <- z
    notes <- character()
    fit <- withCallingHandlers(tryCatch(survival::coxph(survival::Surv(OS_time, OS) ~ z, data = fit_data,
        ties = "efron", control = survival::coxph.control(iter.max = 30), model = FALSE, x = FALSE, y = FALSE),
        error = function(e) structure(list(message = conditionMessage(e)), class = "cox_error")), warning = function(w) {
        notes <<- c(notes, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (inherits(fit, "cox_error"))
        return(empty_fit("fit_error"))
    if (length(notes))
        return(empty_fit("fit_warning"))
    b <- unname(stats::coef(fit)[1L])
    V <- fit$var[1L, 1L]
    if (!is.finite(b) || !is.finite(V) || V <= 0)
        return(empty_fit("invalid_estimate"))
    se <- sqrt(V)
    # Large coefficients and unstable confidence intervals are not ranked as discoveries.
    if (abs(b) > 15 || !is.finite(se) || abs(b/se) > 1e+05)
        return(empty_fit("extreme_or_unstable_estimate"))
    p <- 2 * stats::pnorm(-abs(b/se))
    list(beta = b, SE = se, HR = exp(b), CI_low = exp(b - 1.95996398454 * se), CI_high = exp(b + 1.95996398454 *
        se), P = p, status = "ok")
}

all_summaries <- vector("list", length(cancers))
names(all_summaries) <- cancers
all_tops <- vector("list", length(cancers))
names(all_tops) <- cancers
all_status <- vector("list", length(cancers))
names(all_status) <- cancers
all_feature_files <- character()
figure_dir <- file.path(output_dir, "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

for (cancer in cancers) {
    input_file <- file.path(input_dir, paste0("TCGA_", cancer, "_discovery.rds"))
    if (!file.exists(input_file))
        stop("Missing Stage03 cohort: ", input_file)
    message("Stage05: fitting Cox models for ", cancer, " ...")
    cohort <- readRDS(input_file)
    if (!is.list(cohort) || !all(c("expression", "clinical") %in% names(cohort)))
        stop("Unexpected Stage03 object structure for ", cancer)
    expr <- cohort$expression
    clin <- data.table::as.data.table(cohort$clinical)
    if (!is.matrix(expr) || !is.numeric(expr) || is.null(rownames(expr)) || !all(c("sample_id", "patient_id",
        "OS", "OS_time") %in% names(clin)))
        stop("Invalid expression/clinical schema for ", cancer)
    if (!identical(colnames(expr), as.character(clin$sample_id)) || anyDuplicated(rownames(expr)) ||
        anyDuplicated(clin$patient_id) || anyNA(expr) || any(!is.finite(expr)) || anyNA(clin$OS) || anyNA(clin$OS_time) ||
        !all(clin$OS %in% c(0, 1)) || !all(is.finite(clin$OS_time) & clin$OS_time > 0))
        stop("Cohort validation failed for ", cancer)
    if (sum(clin$OS == 1L) < 10L)
        stop("Insufficient survival events: ", cancer)

    # Genes are rows; patients are columns. Preserve native Stage03 matrices.
    mu <- matrixStats::rowMeans2(expr)
    sdv <- matrixStats::rowSds(expr)
    eligible <- is.finite(sdv) & sdv > 1e-10 & is.finite(mu)
    gene_ids <- rownames(expr)[eligible]
    expr <- expr[eligible, , drop = FALSE]
    mu <- mu[eligible]
    sdv <- sdv[eligible]
    # row-wise scale; each coefficient is per one cancer-specific expression SD.
    expr <- sweep(expr, 1L, mu, "-")
    expr <- sweep(expr, 1L, sdv, "/")
    if (any(!is.finite(expr)))
        stop("Standardisation produced nonfinite data in ", cancer)
    fit_data <- data.frame(OS_time = as.numeric(clin$OS_time), OS = as.integer(clin$OS))
    n <- nrow(expr)
    out <- vector("list", n)
    for (j in seq_len(n)) {
        out[[j]] <- fit_one(as.numeric(expr[j, ]), fit_data)
        if (j%%2500L == 0L)
            message("  ", cancer, ": ", j, "/", n, " features")
    }
    res <- data.table(Cancer = cancer, Feature_ID = gene_ids, Beta = vapply(out, `[[`, numeric(1), "beta"),
        SE = vapply(out, `[[`, numeric(1), "SE"), HR = vapply(out, `[[`, numeric(1), "HR"), CI_low = vapply(out,
            `[[`, numeric(1), "CI_low"), CI_high = vapply(out, `[[`, numeric(1), "CI_high"), P = vapply(out,
            `[[`, numeric(1), "P"), Fit_status = vapply(out, `[[`, character(1), "status"))
    # n = total eligible features, including fits that failed numerically.  Setting n explicitly
    # avoids under-adjustment when P-values are missing.
    res[, FDR := stats::p.adjust(P, method = "BH", n = .N)]
    res[, Direction := fifelse(is.finite(HR) & HR > 1, "Adverse", fifelse(is.finite(HR) & HR < 1,
        "Favourable", NA_character_))]
    res[, `:=`(Patients = nrow(clin), Events = sum(clin$OS == 1L))]
    data.table::setorder(res, FDR, P, Feature_ID, na.last = TRUE)
    outfile <- file.path(output_dir, paste0("Stage05_", cancer, "_all_gene_Cox.csv"))
    data.table::fwrite(res, outfile)
    all_feature_files <- c(all_feature_files, outfile)
    all_summaries[[cancer]] <- data.table(Cancer = cancer, Patients = nrow(clin), Events = sum(clin$OS ==
        1L), Features_tested = nrow(res), Successful_models = sum(res$Fit_status == "ok"), Failed_or_flagged = sum(res$Fit_status !=
        "ok"), Nominal_P_lt_0_05 = sum(res$P < 0.05, na.rm = TRUE), FDR_lt_0_05 = sum(res$FDR < 0.05,
        na.rm = TRUE), FDR_adverse = sum(res$FDR < 0.05 & res$HR > 1, na.rm = TRUE), FDR_favourable = sum(res$FDR <
        0.05 & res$HR < 1, na.rm = TRUE))
    all_tops[[cancer]] <- res[seq_len(min(10L, .N)), .(Cancer, Feature_ID, HR, CI_low, CI_high, P, FDR,
        Fit_status)]
    all_status[[cancer]] <- res[, .N, by = .(Cancer, Fit_status)]
    rm(cohort, expr, clin, res, out)
    invisible(gc(verbose = FALSE))
}
summary_dt <- data.table::rbindlist(all_summaries)
top_dt <- data.table::rbindlist(all_tops)
status_dt <- data.table::rbindlist(all_status)
if (sum(summary_dt$Patients) != 4718L || sum(summary_dt$Events) != 1519L) stop("Cross-cancer patient/event totals differ from validated Stage03 totals")
data.table::fwrite(summary_dt, file.path(output_dir, "Stage05_cancer_specific_summary.csv"))
data.table::fwrite(top_dt, file.path(output_dir, "Stage05_top10_per_cancer.csv"))
data.table::fwrite(status_dt, file.path(output_dir, "Stage05_model_status.csv"))


# JOURNAL-ORIENTED STANDALONE FIGURES; do not use unadjusted associations to claim independence.
# All full genome-wide CSV results are retained, including unsuccessful model fits.
make_theme <- function() ggplot2::theme_minimal(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 18), plot.subtitle = ggplot2::element_text(colour = "black", size = 11),
    plot.caption = ggplot2::element_text(colour = "black", size = 9), axis.title = ggplot2::element_text(face = "bold",
        colour = "black"), axis.text = ggplot2::element_text(face = "bold", colour = "black"), legend.title = ggplot2::element_text(face = "bold",
        colour = "black"), panel.grid.minor = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
save_publication_webp <- function(plot, filename, w = 3600L, h = 2600L) {
    path <- file.path(figure_dir, filename)
    ragg::agg_webp(filename = path, width = w, height = h, res = 300, quality = 96)
    tryCatch(print(plot), finally = grDevices::dev.off())
    path
}
figures <- character()
# Top-line landscape: discovery fractions are not comparable across cancers without considering
# event counts, gene variation, biological composition and study design.
sumplot <- data.table::copy(summary_dt)
sumplot[, Pct_significant := 100 * FDR_lt_0_05/Features_tested]
sumplot[, Cancer := factor(Cancer, levels = Cancer[order(Pct_significant)])]
p_overview <- ggplot2::ggplot(sumplot, ggplot2::aes(x = Pct_significant, y = Cancer)) + ggplot2::geom_segment(ggplot2::aes(x = 0,
    xend = Pct_significant, yend = Cancer), colour = "#D0D7DF", linewidth = 1.2) + ggplot2::geom_point(ggplot2::aes(size = Events,
    colour = Pct_significant), alpha = 0.94) + ggplot2::scale_colour_gradient(low = "#3D7D9C", high = "#B34A58") +
    ggplot2::scale_size_continuous(range = c(4, 11)) + ggplot2::scale_x_continuous(labels = function(x) paste0(round(x,
    1), "%"), expand = ggplot2::expansion(mult = c(0, 0.15))) + ggplot2::labs(title = "Genome-wide survival associations vary considerably across cancers",
    subtitle = "Univariable Cox models | BH-adjusted FDR < 0.05 within each cancer", x = "Proportion of eligible features associated with overall survival",
    y = NULL, size = "Observed deaths", colour = "FDR-positive (%)", caption = "Descriptive discovery burden; not evidence that one cancer is biologically more prognostic than another.") +
    make_theme()
figures <- c(figures, save_publication_webp(p_overview, "Stage05_genomewide_FDR_discovery_landscape.webp"))

# A separate large volcano figure for each cancer; all points correspond to an actual successful
# fit.  X-axis clipped for legibility, with exact untruncated HR/beta preserved in CSV.
for (cc in cancers) {
    d <- data.table::fread(file.path(output_dir, paste0("Stage05_", cc, "_all_gene_Cox.csv")))
    d <- d[Fit_status == "ok" & is.finite(Beta) & is.finite(P) & is.finite(FDR)]
    if (!nrow(d))
        next
    d[, Status := data.table::fifelse(FDR < 0.05 & HR > 1, "Adverse association", data.table::fifelse(FDR <
        0.05 & HR < 1, "Favourable association", "Not FDR-significant"))]
    d[, NegLogP := -log10(pmax(P, .Machine$double.xmin))]
    xx <- stats::quantile(abs(d$Beta), 0.995, na.rm = TRUE, names = FALSE)
    xx <- max(0.75, as.numeric(xx))
    d[, Plot_beta := pmax(-xx, pmin(xx, Beta))]
    d[, Status := factor(Status, levels = c("Not FDR-significant", "Favourable association", "Adverse association"))]
    # Label only 6 high-confidence symbol-like features per cancer.
    labs <- d[Status != "Not FDR-significant" & grepl("^[A-Za-z][A-Za-z0-9._-]+$", Feature_ID)]
    data.table::setorder(labs, FDR, P)
    labs <- head(labs, 6L)
    pp <- ggplot2::ggplot(d, ggplot2::aes(Plot_beta, NegLogP, colour = Status)) + ggplot2::geom_point(size = 0.65,
        alpha = 0.48) + ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#555555",
        linewidth = 0.55) + ggplot2::scale_colour_manual(values = c(`Not FDR-significant` = "#BDC6CF",
        `Favourable association` = "#277E91", `Adverse association` = "#B34A58"), drop = FALSE) + ggplot2::geom_text(data = labs,
        ggplot2::aes(label = Feature_ID), colour = "black", fontface = "bold", size = 3.2, check_overlap = TRUE,
        vjust = -0.8, show.legend = FALSE) + ggplot2::labs(title = paste0(cc, " | Genome-wide unadjusted prognostic associations"),
        subtitle = sprintf("%s patients | %s deaths | %s/%s FDR-positive features", summary_dt[Cancer ==
            cc, Patients], summary_dt[Cancer == cc, Events], summary_dt[Cancer == cc, FDR_lt_0_05], summary_dt[Cancer ==
            cc, Features_tested]), x = "Cox log(HR) per 1 within-cancer SD (extreme values clipped)",
        y = expression(-log[10](italic(p))), colour = "Within-cancer BH FDR", caption = "Y-axis uses nominal p-values; point categories use BH FDR. Exploratory and not clinically adjusted.") +
        make_theme() + ggplot2::theme(legend.position = "bottom")
    figures <- c(figures, save_publication_webp(pp, paste0("Stage05_", cc, "_genomewide_volcano.webp"),
        w = 3900, h = 2700))
    rm(d, pp, labs)
    invisible(gc(FALSE))
}
# Prioritised illustration: top 2 genes per selected cancer, with accurate 95% CI bars.  These are
# discovery-selected illustrations, NOT independent validation.
for (cc in c("KIRC", "LUAD", "PAAD")) {
    d <- data.table::fread(file.path(output_dir, paste0("Stage05_", cc, "_all_gene_Cox.csv")))
    d <- head(d[Fit_status == "ok" & is.finite(HR) & is.finite(CI_low) & is.finite(CI_high) & CI_low >
        0 & CI_high > 0 & grepl("^[A-Za-z][A-Za-z0-9._-]+$", Feature_ID)][order(FDR, P)], 10L)
    if (nrow(d) < 2L)
        next
    d[, Feature_ID := factor(Feature_ID, levels = rev(Feature_ID))]
    pp <- ggplot2::ggplot(d, ggplot2::aes(x = HR, y = Feature_ID)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = "#666666") + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low,
        xmax = CI_high), orientation = "y", width = 0.24, linewidth = 0.85, colour = "#4E6572") + ggplot2::geom_point(ggplot2::aes(colour = HR >
        1), size = 3.7) + ggplot2::scale_colour_manual(values = c(`TRUE` = "#B34A58", `FALSE` = "#277E91"),
        guide = "none") + ggplot2::scale_x_log10() + ggplot2::labs(title = paste0(cc, " | Top discovery associations"),
        subtitle = "Hazard ratios per within-cancer expression SD; horizontal bars = 95% Wald confidence intervals",
        x = "Hazard ratio (log scale)", y = NULL, caption = "Top 10 features selected using these same outcomes; effects are unadjusted and selection-biased.") +
        make_theme()
    figures <- c(figures, save_publication_webp(pp, paste0("Stage05_", cc, "_top10_Cox_forest_95CI.webp"),
        w = 3400, h = 2700))
}
data.table::fwrite(data.table(Figure_file = figures), file.path(output_dir, "Stage05_figure_manifest.csv"))

# A single reproducible final console section, also saved verbatim as text.
report <- capture.output({
    cat("====================================================================\n")
    cat("PROJECT 3 | STAGE 05 \342\200\224 CONSOLIDATED RESULTS\n")
    cat("====================================================================\n")
    cat("1. CANCER-SPECIFIC GENOME-WIDE PROGNOSTIC ASSOCIATIONS\n")
    print(summary_dt)
    cat("\n2. FIT STATUS BY CANCER\n")
    print(status_dt)
    cat("\n3. TOP 10 FEATURES PER CANCER (EXPLORATORY; NOT CLINICALLY ADJUSTED)\n")
    for (cancer in cancers) {
        cat("\n--- ", cancer, " ---\n", sep = "")
        print(top_dt[Cancer == cancer])
    }
    cat("\n4. ANALYTICAL INTERPRETATION\n")
    cat("HR represents one within-cancer standard deviation increase in expression.\n")
    cat("FDR is Benjamini-Hochberg correction within each cancer over all eligible features.\n")
    cat("Gene IDs are original source identifiers, not yet a fully harmonised annotation.\n")
    cat("These are UNADJUSTED associations, not validated independent biomarkers.\n")
    cat("PH assumption, confounding, extreme expression and stability require later evaluation.\n")
    cat("Stage 06: pan-cancer cancer-stratified Cox models and effect heterogeneity.\n")
    cat("\n5. FIGURES GENERATED (WEBP)\n")
    print(figures)
    cat("\n6. OUTPUT DIRECTORY\n", output_dir, "\n", sep = "")
    cat("====================================================================\n")
})
writeLines(report, file.path(output_dir, "Stage05_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
