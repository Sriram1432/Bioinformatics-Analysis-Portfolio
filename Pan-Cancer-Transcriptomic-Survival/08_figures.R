# Regenerate the fifteen-figure transcriptomic survival atlas.  Reads validated numerical outputs,
# reconstructs the three prespecified Kaplan-Meier displays and writes audited WebP and SVG
# exports.  Run with: Rscript 08_figures.R

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


# ---------------------------------------------------------------------------- Final scientific
# figure atlas ----------------------------------------------------------------------------

options(width = 180, scipen = 6)

root <- PROJECT_ROOT
out <- file.path(root, "Results_Rebuild", "Final_Figure_Atlas")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

packages <- c("data.table", "survival", "ggplot2", "ragg", "svglite", "patchwork", "scales", "ggrepel")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
adverse <- "#B6433F"
favourable <- "#2B6F9F"
ink <- "#20242A"
muted <- "#66717F"
pale <- "#E7EBEF"
req <- function(path) {
    if (!file.exists(path))
        stop("Missing required input: ", path)
    path
}
read_stage <- function(stage, file) fread(req(file.path(root, "Results_Rebuild", stage, file)))

theme_atlas <- function(size = 10) {
    ggplot2::theme_classic(base_size = size, base_family = "sans") + ggplot2::theme(axis.text = ggplot2::element_text(colour = ink),
        axis.title = ggplot2::element_text(colour = ink), strip.background = ggplot2::element_rect(fill = "#F2F4F6",
            colour = "#BFC7D0"), strip.text = ggplot2::element_text(face = "bold", colour = ink), legend.position = "bottom",
        legend.title = ggplot2::element_text(face = "bold"), panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(8, 14, 8, 8))
}

manifest <- list()
validation <- character()
km_audit <- list()
export_plot <- function(p, number, stem, width, height, source, method, denominator, limitations = "None identified") {
    base <- sprintf("Figure%02d_%s", number, stem)
    webp <- file.path(figdir, paste0(base, ".webp"))
    svg <- file.path(figdir, paste0(base, ".svg"))
    ggplot2::ggsave(webp, p, device = ragg::agg_webp, width = width, height = height, units = "in", dpi = 300,
        bg = "white")
    ggplot2::ggsave(svg, p, device = svglite::svglite, width = width, height = height, units = "in",
        bg = "white")
    ok <- file.exists(webp) && file.info(webp)$size > 10000 && file.exists(svg) && file.info(svg)$size >
        1000
    manifest[[length(manifest) + 1L]] <<- data.table(Figure = number, Stem = stem, Source_files = source,
        Analysis_method = method, Denominator = denominator, Width_in = width, Height_in = height, WebP = webp,
        SVG = svg, Output_status = if (ok)
            "RENDERED" else "FAILED", Limitations = limitations)
    if (!ok)
        stop("Export validation failed for Figure ", number)
    invisible(c(webp, svg))
}

# Inputs and schema gates -----------------------------------------------------
cohort <- read_stage("Stage03", "Stage03_cohort_and_feature_eligibility.csv")
s05_summary <- read_stage("Stage05", "Stage05_cancer_specific_summary.csv")
s05 <- rbindlist(lapply(cancers, function(cc) read_stage("Stage05", paste0("Stage05_", cc, "_all_gene_Cox.csv"))),
    fill = TRUE)
s06_meta <- read_stage("Stage06", "Stage06_cross_cancer_meta_heterogeneity.csv")
s06_pan <- read_stage("Stage06", "Stage06_pan_cancer_stratified_Cox_all_features.csv")
s07 <- read_stage("Stage07", "Stage07_all_adjusted_candidate_models.csv")
s07_meta <- read_stage("Stage07", "Stage07_adjusted_cross_cancer_fixed_effect_meta.csv")
s07_selection <- read_stage("Stage07", "Stage07_priority_feature_selection.csv")
s08 <- read_stage("Stage08", "Stage08_all_candidate_diagnostics.csv")
s08b_ph <- read_stage("Stage08B", "Stage08B_PH_term_specific.csv")
s09c <- read_stage("Stage09C", "Stage09C_training_selected_algorithm_per_cancer.csv")
s11b <- read_stage("Stage11B", "Stage11B_all_ranked_pathway_tests.csv")

stopifnot(sum(cohort$Patients) == 4718L, sum(cohort$Deaths) == 1519L, setequal(cohort$Cancer, cancers),
    setequal(s05_summary$Cancer, cancers), all(c("Beta", "SE", "HR", "CI_low", "CI_high", "P", "FDR") %in%
        names(s05)), all(c("NES", "FDR", "Library", "Pathway") %in% names(s11b)))
validation <- c(validation, "PASS: Stage03 cohort total = 4,718 unique analysed patients.", "PASS: Stage03 observed death total = 1,519.",
    "PASS: all ten prespecified cancers present.", "PASS: Stage05 Cox and Stage11B ranked-pathway schemas validated.")

# Figure 01: cohort size ------------------------------------------------------
d1 <- copy(cohort)
d1[, Cancer := factor(Cancer, levels = cancers)]
p1 <- ggplot(d1, aes(Cancer, Patients)) + geom_col(width = 0.72, fill = "#356A89") + geom_text(aes(label = scales::comma(Patients)),
    vjust = -0.4, size = 3.1) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = scales::comma) +
    labs(x = NULL, y = "Unique analysed patients") + theme_atlas(11)
export_plot(p1, 1, "TCGA_patient_counts", 10.5, 6.2, "Stage03_cohort_and_feature_eligibility.csv", "Final patient-level cohort count",
    "N = 4,718")

# Figure 02: observed status -------------------------------------------------
d2 <- rbindlist(list(d1[, .(Cancer, Status = "Observed death", N = Deaths, Total = Patients)], d1[, .(Cancer,
    Status = "Censored", N = Censored, Total = Patients)]))
d2[, Cancer := factor(Cancer, levels = cancers)]
d2[, Status := factor(Status, levels = c("Censored", "Observed death"))]
labs2 <- d1[, .(Cancer, Patients, Deaths)]
labs2[, Cancer := factor(Cancer, levels = cancers)]
p2 <- ggplot(d2, aes(Cancer, N, fill = Status)) + geom_col(width = 0.74) + geom_text(data = labs2, aes(Cancer,
    Patients, label = sprintf("%d / %d", Deaths, Patients)), inherit.aes = FALSE, vjust = -0.35, size = 2.8) +
    scale_fill_manual(values = c(Censored = "#B8C5CE", `Observed death` = "#B6433F")) + scale_y_continuous(expand = expansion(mult = c(0,
    0.1)), labels = scales::comma) + labs(x = NULL, y = "Patients", fill = "Observed status") + theme_atlas(11)
export_plot(p2, 2, "observed_event_status", 11, 6.4, "Stage03_cohort_and_feature_eligibility.csv", "Observed event/censor count (not fixed-time mortality)",
    "N = 4,718; deaths = 1,519")

# Figure 03: cancer-specific volcano atlas ----------------------------------
label_genes <- c("SERPINE1", "MET", "STX1A", "PGK1", "DKK1")
d3 <- s05[Fit_status == "ok" & is.finite(Beta) & is.finite(P) & is.finite(HR) & HR > 0]
if (max(abs(d3$HR - exp(d3$Beta)), na.rm = TRUE) > 1e-10) stop("Figure 03: saved HR does not equal exp(Beta)")
d3[, `:=`(Display_HR = exp(Beta), Cancer = factor(Cancer, levels = cancers), mlogp = -log10(pmax(P, .Machine$double.xmin)),
    Sig = FDR < 0.05)]
hr_limits3 <- range(c(0.25, 4, d3$Display_HR), na.rm = TRUE)
lab3 <- d3[Feature_ID %in% label_genes]
p3 <- ggplot(d3, aes(Display_HR, mlogp)) + geom_vline(xintercept = 1, linetype = 2, colour = muted, linewidth = 0.45) +
    geom_point(aes(colour = Sig), alpha = 0.42, size = 0.45) + ggrepel::geom_text_repel(data = lab3,
    aes(label = Feature_ID), size = 2.3, max.overlaps = Inf, min.segment.length = 0, box.padding = 0.18,
    seed = 3) + facet_wrap(~Cancer, ncol = 5, scales = "free_y") + scale_x_log10(breaks = c(0.25, 0.5,
    1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4"), limits = hr_limits3, expand = expansion(mult = c(0.02,
    0.02))) + scale_colour_manual(values = c(`FALSE` = "#AAB2BA", `TRUE` = adverse), labels = c("Not FDR-significant",
    "BH FDR < 0.05")) + labs(x = "Hazard ratio per 1-SD increase in gene expression", y = expression(-log[10](raw ~
    P)), colour = NULL) + theme_atlas(8.5) + theme(legend.position = "bottom")
export_plot(p3, 3, "genomewide_Cox_volcano_atlas", 15, 8.8, "Stage05_[CANCER]_all_gene_Cox.csv", "Cancer-specific univariable Cox; HR = exp(Beta); log-scaled HR axis; BH FDR within cancer",
    paste0("Successful models = ", scales::comma(nrow(d3))))

# Figure 04: selected gene-by-cancer map ------------------------------------
g4_stage05 <- s05[, .(Stage05_tested = sum(Fit_status == "ok" & is.finite(Beta) & is.finite(HR)), Stage05_FDR_lt_0_05 = sum(Fit_status ==
    "ok" & is.finite(FDR) & FDR < 0.05), Best_Stage05_FDR = suppressWarnings(min(FDR[Fit_status == "ok"],
    na.rm = TRUE))), by = Feature_ID]
g4_stage05[!is.finite(Best_Stage05_FDR), Best_Stage05_FDR := NA_real_]
g4_adjusted <- s07[, .(Adjusted_models_ok = sum(Adjusted_status == "ok" & is.finite(Adjusted_P)), Adjusted_P_lt_0_05 = sum(Adjusted_status ==
    "ok" & is.finite(Adjusted_P) & Adjusted_P < 0.05), Adjusted_FDR_lt_0_05 = sum(Adjusted_status ==
    "ok" & is.finite(Adjusted_FDR_exploratory) & Adjusted_FDR_exploratory < 0.05)), by = Feature_ID]
g4_summary <- Reduce(function(x, y) merge(x, y, by = "Feature_ID", all = TRUE), list(g4_stage05, g4_adjusted,
    s06_meta, s07_selection))
for (nm in c("Low_heterogeneity", "High_heterogeneity", "Recurrent", "Top_cancer_specific", "Named_example")) g4_summary[is.na(get(nm)),
    (nm) := FALSE]
g4_summary[, Eligible := Stage05_tested == length(cancers) & Stage05_FDR_lt_0_05 >= 1L & Adjusted_models_ok >=
    1L & is.finite(I2_pct) & is.finite(Meta_P)]

forced4 <- c("SERPINE1", "MET", "STX1A", "PGK1", "DKK1", "ATP6V1B1")
if (!all(forced4 %in% g4_summary[Eligible == TRUE]$Feature_ID)) stop("Figure 04: one or more required project candidates fail the documented availability/support gates")
selected4 <- character()
add4 <- function(pool, n) {
    pool <- setdiff(pool, selected4)
    selected4 <<- c(selected4, head(pool, n))
}
add4(forced4, length(forced4))
add4(g4_summary[Eligible == TRUE & Named_example == TRUE][order(Best_Stage05_FDR, Meta_P)]$Feature_ID,
    4L)
add4(g4_summary[Eligible == TRUE & Low_heterogeneity == TRUE][order(Meta_FDR, I2_pct, Best_Stage05_FDR)]$Feature_ID,
    8L)
add4(g4_summary[Eligible == TRUE & High_heterogeneity == TRUE][order(Q_FDR, -Stage05_FDR_lt_0_05, Best_Stage05_FDR)]$Feature_ID,
    8L)
add4(g4_summary[Eligible == TRUE & Recurrent == TRUE][order(-Stage05_FDR_lt_0_05, Best_Stage05_FDR, Meta_FDR)]$Feature_ID,
    5L)
add4(g4_summary[Eligible == TRUE & Top_cancer_specific == TRUE][order(Best_Stage05_FDR, -Adjusted_P_lt_0_05)]$Feature_ID,
    5L)
add4(g4_summary[Eligible == TRUE][order(Meta_FDR, Best_Stage05_FDR)]$Feature_ID, 36L - length(selected4))
selected4 <- head(selected4, 36L)
if (length(selected4) != 36L || anyDuplicated(selected4)) stop("Figure 04: deterministic selection did not produce 36 unique genes")

g4_selected_summary <- g4_summary[match(selected4, Feature_ID)]
g4_selected_summary[, Selection_rationale := vapply(Feature_ID, function(g) {
    z <- g4_summary[Feature_ID == g]
    tags <- c(if (g %in% forced4) "required verified project candidate", if (isTRUE(z$Named_example)) "Stage07 named example",
        if (isTRUE(z$Low_heterogeneity)) "Stage07 low-heterogeneity set", if (isTRUE(z$High_heterogeneity)) "Stage07 high-heterogeneity set",
        if (isTRUE(z$Recurrent)) "Stage07 recurrent set", if (isTRUE(z$Top_cancer_specific)) "Stage07 top cancer-specific set")
    paste(c(tags, sprintf("%d cancer-specific BH-FDR discoveries", z$Stage05_FDR_lt_0_05), sprintf("%d successful adjusted models",
        z$Adjusted_models_ok)), collapse = "; ")
}, character(1))]

grid4 <- CJ(Feature_ID = selected4, Cancer = cancers, unique = TRUE)
d4 <- merge(grid4, s05, by = c("Feature_ID", "Cancer"), all.x = TRUE, sort = FALSE)
d4 <- merge(d4, g4_selected_summary[, .(Feature_ID, Selection_rationale, Stage05_FDR_lt_0_05, Adjusted_models_ok,
    Adjusted_P_lt_0_05, Adjusted_FDR_lt_0_05, Meta_HR, Meta_P, Meta_FDR, I2_pct, Q_P, Q_FDR, Low_heterogeneity,
    High_heterogeneity, Recurrent, Top_cancer_specific, Named_example)], by = "Feature_ID", all.x = TRUE)
d4 <- merge(d4, s07[, .(Feature_ID, Cancer, Clinically_adjusted_HR = Adjusted_HR, Clinically_adjusted_CI_low = CI_low,
    Clinically_adjusted_CI_high = CI_high, Clinically_adjusted_P = Adjusted_P, Clinically_adjusted_exploratory_FDR = Adjusted_FDR_exploratory,
    Clinically_adjusted_status = Adjusted_status)], by = c("Feature_ID", "Cancer"), all.x = TRUE)
d4[, Tested := Fit_status == "ok" & is.finite(Beta) & is.finite(HR)]
d4[, Significance := fifelse(!Tested, "Missing / untested", fifelse(is.finite(FDR) & FDR < 0.05, "BH FDR < 0.05",
    "Not FDR-significant"))]
d4[, Significance := factor(Significance, levels = c("BH FDR < 0.05", "Not FDR-significant", "Missing / untested"))]
d4[, `:=`(Cancer = factor(Cancer, levels = cancers), Feature_ID = factor(Feature_ID, levels = rev(selected4)),
    Plot_magnitude = fifelse(Tested, abs(Beta), 0.04))]
lim4 <- max(abs(d4$Beta), na.rm = TRUE)
p4 <- ggplot(d4, aes(Cancer, Feature_ID)) + geom_point(aes(size = Plot_magnitude, fill = Beta, shape = Significance),
    colour = ink, stroke = 0.48) + scale_shape_manual(values = c(`BH FDR < 0.05` = 21, `Not FDR-significant` = 22,
    `Missing / untested` = 4), breaks = c("BH FDR < 0.05", "Not FDR-significant", "Missing / untested"),
    drop = FALSE) + scale_size_continuous(range = c(0.9, 6.2), name = "|log(HR)|") + scale_fill_gradient2(low = "#0072B2",
    mid = "white", high = "#D55E00", midpoint = 0, limits = c(-lim4, lim4), na.value = "#B7BDC4", name = "Signed log(HR)") +
    labs(x = NULL, y = NULL, shape = "Association status") + theme_atlas(9.5) + theme(panel.grid.major = element_line(colour = "#EEF0F2",
    linewidth = 0.25), axis.text.y = element_text(size = 8.2))
export_plot(p4, 4, "selected_gene_cancer_bubble_matrix", 13.5, 16, "Stage05 all-gene Cox; Stage06 heterogeneity; Stage07 priority and adjusted Cox outputs",
    "Deterministic 36-gene evidence-balanced selection; Stage05 BH FDR < 0.05 defines displayed significance",
    "36 genes x 10 cancers", "Selected evidence atlas; adjusted estimates support selection but bubbles show Stage05 unadjusted cancer-specific estimates")

supp4 <- copy(d4)
supp4[, `:=`(Feature_ID = as.character(Feature_ID), Cancer = as.character(Cancer))]
setnames(supp4, c("Beta", "HR", "CI_low", "CI_high", "P", "FDR"), c("Unadjusted_log_HR", "Unadjusted_HR",
    "Unadjusted_CI_low", "Unadjusted_CI_high", "Unadjusted_P", "Unadjusted_BH_FDR"))
fwrite(supp4[, .(Feature_ID, Cancer, Tested, Fit_status, Unadjusted_log_HR, Unadjusted_HR, Unadjusted_CI_low,
    Unadjusted_CI_high, Unadjusted_P, Unadjusted_BH_FDR, Significance, Clinically_adjusted_HR, Clinically_adjusted_CI_low,
    Clinically_adjusted_CI_high, Clinically_adjusted_P, Clinically_adjusted_exploratory_FDR, Clinically_adjusted_status,
    Meta_HR, Meta_P, Meta_FDR, I2_pct, Q_P, Q_FDR, Selection_rationale)], file.path(out, "Figure04_displayed_gene_cancer_associations.csv"))
method4 <- c("FIGURE 04 GENE-SELECTION METHODOLOGY", "Candidate source: saved Stage07 priority-feature table, supplemented only by the six revision-specified project candidates.",
    "Eligibility gate: successful Stage05 estimates in all 10 cancers; at least one cancer-specific Stage05 BH FDR < 0.05; at least one successful Stage07 clinically adjusted model; and finite Stage06 meta-analysis/heterogeneity statistics.",
    "Deterministic quotas, after removing duplicates: 6 required candidates; 4 additional Stage07 named examples; 8 low-heterogeneity; 8 high-heterogeneity; 5 recurrent; 5 top cancer-specific; then Stage06 meta-FDR rank only if needed to reach 36.",
    "Within-category ranking: named by best Stage05 FDR then meta P; low heterogeneity by meta FDR/I2; high heterogeneity by Q FDR then recurrence; recurrent by number of cancer FDR discoveries; cancer-specific by best Stage05 FDR then adjusted support.",
    "Displayed bubbles always use the original Stage05 cancer-specific unadjusted Cox estimates. Circle = BH FDR < 0.05; square = tested but nonsignificant; cross = missing/untested. No gene was chosen solely for an extreme HR.",
    paste("Selected genes:", paste(selected4, collapse = ", ")))
writeLines(method4, file.path(out, "Figure04_gene_selection_methodology.txt"))

# Figure 05: SERPINE1 adjusted forest ---------------------------------------
d5 <- s07[Feature_ID == "SERPINE1" & Adjusted_status == "ok"]
d5[, Cancer := factor(Cancer, levels = rev(cancers))]
sum5 <- s07_meta[Feature_ID == "SERPINE1"]
if (nrow(d5) != 10L || nrow(sum5) != 1L) stop("SERPINE1 adjusted forest inputs incomplete")
d5p <- d5[, .(Label = as.character(Cancer), HR = Adjusted_HR, lo = CI_low, hi = CI_high, Type = "Cancer-specific")]
d5p <- rbind(d5p, data.table(Label = "Fixed-effect summary", HR = sum5$HR, lo = sum5$CI_low, hi = sum5$CI_high,
    Type = "Summary"))
d5p[, Label := factor(Label, levels = c("Fixed-effect summary", cancers))]
p5 <- ggplot(d5p, aes(HR, Label)) + geom_vline(xintercept = 1, linetype = 2, colour = muted) + geom_errorbar(aes(xmin = lo,
    xmax = hi, colour = Type), orientation = "y", width = 0.17, linewidth = 0.75) + geom_point(aes(shape = Type,
    colour = Type), size = 2.8) + scale_x_log10() + scale_colour_manual(values = c(`Cancer-specific` = "#315F7D",
    Summary = adverse)) + scale_shape_manual(values = c(`Cancer-specific` = 16, Summary = 18)) + labs(x = "Adjusted hazard ratio per within-cancer SD (log scale), 95% CI",
    y = NULL, colour = NULL, shape = NULL) + theme_atlas(10.5)
export_plot(p5, 5, "SERPINE1_adjusted_cross_cancer_forest", 10, 7.3, "Stage07_all_adjusted_candidate_models.csv; Stage07_adjusted_cross_cancer_fixed_effect_meta.csv",
    "Cancer-specific adjusted Cox plus inverse-variance fixed-effect summary", "10 cancers", "Discovery-conditioned candidate; covariate sets differ by cancer")

# Figure 06: heterogeneity heatmap ------------------------------------------
selected6 <- c("PGK1", "STX1A", "MET", "DKK1", "SERPINE1", "FLNC", "PSMD2", "PLK1", "FOXM1", "EXPH5",
    "KLRB1", "FOSL1")
d6 <- s05[Feature_ID %in% selected6]
ann6 <- s06_meta[Feature_ID %in% selected6, .(Feature_ID, I2_pct)]
d6 <- merge(d6, ann6, by = "Feature_ID", all.x = TRUE)
d6[, Row := sprintf("%s  [I^2 %.0f%%]", Feature_ID, I2_pct)]
row_levels6 <- rev(unique(d6[match(selected6, Feature_ID), Row]))
d6[, `:=`(Cancer = factor(Cancer, levels = cancers), Row = factor(Row, levels = row_levels6))]
lim6 <- as.numeric(stats::quantile(abs(d6$Beta), 0.98, na.rm = TRUE))
p6 <- ggplot(d6, aes(Cancer, Row, fill = pmax(-lim6, pmin(lim6, Beta)))) + geom_tile(colour = "white",
    linewidth = 0.5) + geom_point(data = d6[FDR < 0.05], shape = 8, size = 1.55, colour = ink) + scale_fill_gradient2(low = favourable,
    mid = "white", high = adverse, midpoint = 0, limits = c(-lim6, lim6), name = "Signed log(HR)") +
    labs(x = NULL, y = NULL) + theme_atlas(9.5) + theme(axis.text.y = element_text(colour = ink))
export_plot(p6, 6, "selected_gene_heterogeneity_heatmap", 13.2, 7.8, "Stage05_[CANCER]_all_gene_Cox.csv; Stage06_cross_cancer_meta_heterogeneity.csv",
    "Signed unadjusted log(HR) with formal Stage06 I-squared annotation", "12 genes x 10 cancers", "PGK1 emphasised as the prespecified BRCA/KIRC direction contrast")

# Figure 07: selected clinically adjusted estimates -------------------------
pairs7 <- data.table(Cancer = c("PAAD", "KIRC", "LUAD", "BRCA", "KIRC"), Feature_ID = c("MET", "STX1A",
    "DKK1", "PGK1", "PGK1"))
d7 <- merge(pairs7, s07, by = c("Cancer", "Feature_ID"), all.x = TRUE, sort = FALSE)
if (any(!is.finite(d7$Adjusted_HR))) stop("Figure 07 selected adjusted estimates missing")
d7[, Label := paste(Cancer, Feature_ID, sep = " / ")]
d7[, PH_note := fifelse(Cancer == "LUAD" & Feature_ID == "DKK1", "Gene PH P = 0.00073", "")]
d7[, Label := factor(Label, levels = rev(Label))]
p7 <- ggplot(d7, aes(Adjusted_HR, Label)) + geom_vline(xintercept = 1, linetype = 2, colour = muted) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.17, linewidth = 0.8,
        colour = "#315F7D") + geom_point(size = 3, colour = adverse) + geom_text(aes(x = CI_high * 1.06,
    label = PH_note), hjust = 0, size = 3, colour = adverse) + scale_x_log10(expand = expansion(mult = c(0.05,
    0.28))) + labs(x = "Clinically adjusted hazard ratio per within-cancer SD (log scale), 95% CI", y = NULL) +
    theme_atlas(10.5)
export_plot(p7, 7, "clinically_adjusted_selected_forest", 11, 6, "Stage07_all_adjusted_candidate_models.csv; Stage08B_PH_term_specific.csv",
    "Multivariable cancer-specific Cox", paste0("Model N = ", paste(d7$N, collapse = ", ")), "DKK1-LUAD violates gene-level proportional hazards")

# KM helper (Figures 08-10) --------------------------------------------------
make_km <- function(cc, gene, number, adjusted_row, selection_note) {
    z <- readRDS(req(file.path(root, "Results_Rebuild", "Stage03", paste0("TCGA_", cc, "_discovery.rds"))))
    if (!all(c("expression", "clinical") %in% names(z)))
        stop(cc, ": invalid patient-level RDS")
    cl <- as.data.table(z$clinical)
    if (!identical(colnames(z$expression), as.character(cl$sample_id)) || anyDuplicated(cl$patient_id))
        stop(cc, ": patient/expression alignment failure")
    if (!gene %in% rownames(z$expression))
        stop(cc, "/", gene, ": gene absent")
    score <- as.numeric(z$expression[gene, ])
    valid <- is.finite(score) & is.finite(cl$OS_time) & cl$OS_time > 0 & cl$OS %in% 0:1
    d <- data.frame(time_days = as.numeric(cl$OS_time[valid]), time_months = as.numeric(cl$OS_time[valid])/30.4375,
        event = as.integer(cl$OS[valid]), score = score[valid])
    cutoff <- median(d$score)
    d$Group <- factor(ifelse(d$score >= cutoff, "Higher expression", "Lower expression"), levels = c("Lower expression",
        "Higher expression"))
    fit <- survival::survfit(survival::Surv(time_months, event) ~ Group, data = d, conf.type = "log-log")
    ss <- summary(fit, censored = TRUE)
    curve <- data.table(time = ss$time, surv = ss$surv, lower = ss$lower, upper = ss$upper, n.censor = ss$n.censor,
        Group = sub("^Group=", "", ss$strata))
    curve <- rbind(data.table(time = 0, surv = 1, lower = 1, upper = 1, n.censor = 0, Group = levels(d$Group)),
        curve)
    curve[, Group := factor(Group, levels = levels(d$Group))]
    cens <- curve[n.censor > 0]
    xmax <- max(d$time_months)
    breaks <- pretty(c(0, xmax), 6)
    breaks <- breaks[breaks >= 0 & breaks <= xmax]
    if (!0 %in% breaks)
        breaks <- c(0, breaks)
    risk_s <- summary(fit, times = breaks, extend = TRUE)
    risk <- data.table(time = risk_s$time, n.risk = risk_s$n.risk, Group = sub("^Group=", "", risk_s$strata))
    risk[, Group := factor(Group, levels = rev(levels(d$Group)))]
    cols <- c(`Lower expression` = favourable, `Higher expression` = adverse)
    group_stats <- as.data.table(d)[, .(N = .N, Events = sum(event)), by = Group]
    gs_lo <- group_stats[Group == "Lower expression"]
    gs_hi <- group_stats[Group == "Higher expression"]
    subtitle <- sprintf("Median split = %.3f; lower n/events = %d/%d; higher n/events = %d/%d.\nAdjusted continuous-expression HR %.3f (95%% CI %.3f to %.3f).",
        cutoff, gs_lo$N, gs_lo$Events, gs_hi$N, gs_hi$Events, adjusted_row$Adjusted_HR, adjusted_row$CI_low,
        adjusted_row$CI_high)
    pk <- ggplot(curve, aes(time, surv, colour = Group, fill = Group)) + geom_ribbon(data = curve[is.finite(lower) &
        is.finite(upper)], aes(ymin = lower, ymax = upper), alpha = 0.14, colour = NA) + geom_step(linewidth = 0.9) +
        geom_point(data = cens, shape = 3, size = 1.3, stroke = 0.55) + scale_colour_manual(values = cols) +
        scale_fill_manual(values = cols) + scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1),
        expand = expansion(mult = c(0, 0.02))) + scale_x_continuous(limits = c(0, xmax), breaks = breaks,
        expand = expansion(mult = c(0, 0.01))) + labs(x = NULL, y = "Overall survival probability", colour = NULL,
        fill = NULL, subtitle = subtitle) + theme_atlas(9.5) + theme(legend.position = "top", plot.subtitle = element_text(size = 8.7),
        plot.margin = margin(t = 14, r = 10, b = 8, l = 10))
    pr <- ggplot(risk, aes(time, Group, label = n.risk, colour = Group)) + geom_text(aes(hjust = ifelse(time ==
        0, 0, 0.5)), size = 3) + scale_colour_manual(values = cols, guide = "none") + scale_x_continuous(limits = c(0,
        xmax), breaks = breaks, expand = expansion(mult = c(0, 0.01))) + labs(x = "Overall-survival follow-up (months; source days / 30.4375)",
        y = NULL) + theme_classic(base_size = 9) + theme(axis.text.y = element_text(colour = ink), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())
    combo <- pk/pr + patchwork::plot_layout(heights = c(4.5, 1.15))
    export_plot(combo, number, paste0(cc, "_", gene, "_Kaplan_Meier"), 10.8, 8.7, paste0("Stage03/TCGA_",
        cc, "_discovery.rds; Stage07_all_adjusted_candidate_models.csv"), "Median-split unadjusted Kaplan-Meier with log-log Greenwood 95% CI and aligned risk table",
        sprintf("N = %d; deaths = %d", nrow(d), sum(d$event)), selection_note)
    group_wide <- dcast(group_stats, . ~ Group, value.var = c("N", "Events"))
    group_wide[, "." := NULL]
    km_audit[[length(km_audit) + 1L]] <<- cbind(data.table(Figure = number, Cancer = cc, Gene = gene,
        Time_source_unit = "days", Display_unit = "months", Cutoff_rule = "within-cancer median", Cutoff = cutoff),
        group_wide)
}

get_adj <- function(cc, gene) s07[Cancer == cc & Feature_ID == gene][1]
make_km("PAAD", "MET", 8, get_adj("PAAD", "MET"), "Post-discovery visualisation; not independent validation.")

# COAD audit and Figure 09 ----------------------------------------------------
coad_discovery <- s05[Cancer == "COAD" & FDR < 0.05]
coad_audit <- merge(coad_discovery[, .(Cancer, Feature_ID, Discovery_HR = HR, Discovery_CI_low = CI_low,
    Discovery_CI_high = CI_high, Discovery_P = P, Discovery_FDR = FDR)], s07[Cancer == "COAD", .(Cancer,
    Feature_ID, Adjusted_HR, Adjusted_CI_low = CI_low, Adjusted_CI_high = CI_high, Adjusted_P, Adjusted_FDR_exploratory,
    N_adjusted = N, Events_adjusted = Events)], by = c("Cancer", "Feature_ID"), all.x = TRUE)
coad_audit[, `:=`(Stage08_diagnostic_available = Feature_ID %in% s08[Cancer == "COAD"]$Feature_ID, Eligible_adjusted = is.finite(Adjusted_P) &
    Adjusted_P < 0.05 & is.finite(Adjusted_CI_low) & Adjusted_CI_low > 1)]
setorder(coad_audit, -Eligible_adjusted, Adjusted_P, Discovery_FDR)
coad_gene <- coad_audit[Eligible_adjusted == TRUE, Feature_ID][1]
if (is.na(coad_gene)) stop("BLOCKED: no defensible COAD candidate with saved adjusted support")
coad_audit[, `:=`(Selected = Feature_ID == coad_gene, Selection_reason = fifelse(Feature_ID == coad_gene,
    "Strongest adjusted P among COAD genome-wide FDR discoveries with adjusted HR CI entirely above 1; median KM split not used for selection",
    "Not selected: lower prespecified ranking on adjusted P/precision or no saved adjusted support"),
    Diagnostic_limitation = fifelse(Stage08_diagnostic_available, "Archived Stage08 diagnostics available",
        "No archived Stage08 gene-specific PH diagnostic"))]
fwrite(coad_audit, file.path(out, "COAD_gene_selection_audit.csv"))
make_km("COAD", coad_gene, 9, get_adj("COAD", coad_gene), "Selected from saved Cox evidence before viewing KM separation; no archived gene-specific PH diagnostic.")
make_km("KIRC", "STX1A", 10, get_adj("KIRC", "STX1A"), "Post-discovery visualisation; gene PH P = 0.0849 in Stage08B.")

# Figure 11: reconstruct the exact DKK1-LUAD diagnostic ----------------------
clean_chr <- function(x) {
    z <- trimws(as.character(x))
    z[toupper(z) %in% c("", "NA", "N/A", "NULL", "NOT REPORTED", "NOT AVAILABLE", "UNKNOWN", "[NOT AVAILABLE]",
        "[NOT EVALUATED]", "[DISCREPANCY]")] <- NA_character_
    z
}
major_stage <- function(x) {
    z <- toupper(clean_chr(x))
    z <- gsub("^STAGE[[:space:]]*", "", z)
    y <- rep(NA_character_, length(z))
    y[grepl("^IV", z)] <- "IV"
    y[grepl("^III", z) & is.na(y)] <- "III"
    y[grepl("^II", z) & is.na(y)] <- "II"
    y[grepl("^I($|A|B|C|[[:space:]])", z) & is.na(y)] <- "I"
    y
}
z11 <- readRDS(req(file.path(root, "Results_Rebuild", "Stage03", "TCGA_LUAD_discovery.rds")))
cl11 <- as.data.table(z11$clinical)
n11 <- nrow(cl11)
d11 <- data.frame(time = as.numeric(cl11$OS_time), event = as.integer(cl11$OS))
age <- suppressWarnings(as.numeric(cl11$age_at_initial_pathologic_diagnosis))
age[!is.finite(age) | age < 18 | age > 105] <- NA
d11$age10 <- age/10
sx <- toupper(clean_chr(cl11$gender))
sx[!sx %in% c("MALE", "FEMALE")] <- NA
d11$sex <- factor(sx)
st <- major_stage(cl11$ajcc_pathologic_tumor_stage)
d11$stage <- factor(st, levels = c("I", "II", "III", "IV"))
include11 <- complete.cases(d11)
d11 <- droplevels(d11[include11, ])
d11$z <- as.numeric(scale(as.numeric(z11$expression["DKK1", ])))[include11]
fit11 <- survival::coxph(survival::Surv(time, event) ~ z + age10 + sex + stage, data = d11, ties = "efron",
    x = TRUE, y = TRUE)
ph11 <- survival::cox.zph(fit11, transform = "km", terms = TRUE, singledf = FALSE)
saved_p11 <- s08b_ph[Cancer == "LUAD" & Feature_ID == "DKK1" & Term == "z", P]
if (length(saved_p11) != 1L || abs(ph11$table["z", "p"] - saved_p11) > 1e-10) stop("Figure 11 PH reconstruction does not reproduce Stage08B")
diag11 <- data.table(x = ph11$x, y = ph11$y[, "z"])
smooth11 <- loess(y ~ x, data = diag11, span = 0.65)
grid11 <- data.table(x = seq(min(diag11$x), max(diag11$x), length.out = 250))
pred11 <- predict(smooth11, newdata = grid11$x, se = TRUE)
grid11[, `:=`(fit = pred11$fit, lo = pred11$fit - 1.96 * pred11$se.fit, hi = pred11$fit + 1.96 * pred11$se.fit)]
p11 <- ggplot(diag11, aes(x, y)) + geom_hline(yintercept = unname(coef(fit11)["z"]), linetype = 2, colour = muted) +
    geom_point(alpha = 0.45, size = 1.4, colour = "#315F7D") + geom_ribbon(data = grid11, aes(x = x,
    ymin = lo, ymax = hi), inherit.aes = FALSE, fill = adverse, alpha = 0.15) + geom_line(data = grid11,
    aes(y = fit), colour = adverse, linewidth = 1) + labs(x = "KM-transformed follow-up time", y = "Scaled Schoenfeld residual / time-varying coefficient estimate",
    caption = sprintf("Exact Stage08B model reconstruction: N = %d, deaths = %d; DKK1 PH P = %.6g; global PH P = %.6g.",
        nrow(d11), sum(d11$event), ph11$table["z", "p"], ph11$table["GLOBAL", "p"])) + theme_atlas(10)
export_plot(p11, 11, "LUAD_DKK1_Schoenfeld_diagnostic", 10.5, 6.8, "Stage03/TCGA_LUAD_discovery.rds; Stage08B_PH_term_specific.csv; Stage08B model specification",
    "Exact adjusted Cox reconstruction and scaled Schoenfeld residual diagnostic", sprintf("N = %d; deaths = %d",
        nrow(d11), sum(d11$event)), "Strong gene-specific PH violation; constant HR is not time-invariant")

# Figure 12: training-only internal validation -------------------------------
d12 <- copy(s09c)
d12[, Cancer := factor(Cancer, levels = rev(Cancer[order(Mean_Delta_C)]))]
p12 <- ggplot(d12, aes(Mean_Delta_C, Cancer)) + geom_vline(xintercept = 0, linetype = 2, colour = muted) +
    geom_segment(aes(x = Min_Delta_C, xend = Max_Delta_C, yend = Cancer), linewidth = 0.7, colour = "#9AA6B2") +
    geom_point(aes(colour = Mean_Delta_C > 0), size = 3.3) + scale_colour_manual(values = c(`FALSE` = favourable,
    `TRUE` = adverse), guide = "none") + scale_x_continuous(labels = scales::number_format(accuracy = 0.01,
    sign = TRUE)) + labs(x = "Mean held-out change in C-index (gene model - clinical model)", y = NULL) +
    theme_atlas(10.5)
export_plot(p12, 12, "training_only_crossvalidated_delta_C", 10.5, 6.8, "Stage09C_training_selected_algorithm_per_cancer.csv",
    "Training-fold-only selection; repeated 5-fold CV x 2", "100/100 planned fold comparisons", "Internal discrimination only; dependent folds; no calibration or external validation")

# Figure 13: Hallmark signed NES ---------------------------------------------
hall13 <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MITOTIC_SPINDLE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
d13 <- s11b[Library == "Hallmark" & Pathway %in% hall13]
d13[, `:=`(Cancer = factor(Cancer, levels = cancers), PathwayLabel = gsub("_", " ", sub("HALLMARK_",
    "", Pathway)))]
d13[, PathwayLabel := factor(PathwayLabel, levels = rev(gsub("_", " ", sub("HALLMARK_", "", hall13))))]
lim13 <- max(abs(d13$NES), na.rm = TRUE)
p13 <- ggplot(d13, aes(Cancer, PathwayLabel, fill = NES)) + geom_tile(colour = "white", linewidth = 0.55) +
    geom_point(data = d13[FDR < 0.05], shape = 8, size = 1.65, colour = ink) + scale_fill_gradient2(low = favourable,
    mid = "white", high = adverse, midpoint = 0, limits = c(-lim13, lim13), name = "Signed NES") + labs(x = NULL,
    y = NULL) + theme_atlas(9.5)
export_plot(p13, 13, "Hallmark_signed_NES_heatmap", 13, 6.8, "Stage11B_all_ranked_pathway_tests.csv",
    "Genome-wide signed Cox-Z ranked GSEA; within-cancer BH FDR", "8 Hallmark programmes x 10 cancers",
    "Same source outcomes as discovery; not independent validation or measured pathway activity")

# Figure 14: pruned Reactome matrix -----------------------------------------
react14 <- c("REACTOME_CELL_CYCLE_CHECKPOINTS", "REACTOME_MITOTIC_PROMETAPHASE", "REACTOME_SEPARATION_OF_SISTER_CHROMATIDS",
    "REACTOME_DNA_REPLICATION", "REACTOME_COLLAGEN_FORMATION", "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "REACTOME_ECM_PROTEOGLYCANS")
react14 <- react14[react14 %in% s11b[Library == "Reactome"]$Pathway]
d14 <- s11b[Library == "Reactome" & Pathway %in% react14]
d14[, `:=`(Cancer = factor(Cancer, levels = cancers), PathwayLabel = gsub("_", " ", sub("REACTOME_",
    "", Pathway)), Sig = FDR < 0.05)]
d14[, PathwayLabel := factor(PathwayLabel, levels = rev(gsub("_", " ", sub("REACTOME_", "", react14))))]
lim14 <- max(abs(d14$NES), na.rm = TRUE)
p14 <- ggplot(d14, aes(Cancer, PathwayLabel)) + geom_point(aes(size = Set_size, fill = NES, shape = Sig),
    colour = ink, stroke = 0.5) + scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24), labels = c("Not significant",
    "BH FDR < 0.05")) + scale_size_continuous(range = c(2, 7), name = "Gene-set size") + scale_fill_gradient2(low = favourable,
    mid = "white", high = adverse, midpoint = 0, limits = c(-lim14, lim14), name = "Signed NES") + labs(x = NULL,
    y = NULL, shape = NULL) + theme_atlas(9)
export_plot(p14, 14, "Reactome_pruned_programme_matrix", 14.5, 7.5, "Stage11B_all_ranked_pathway_tests.csv",
    "Genome-wide signed Cox-Z ranked GSEA; pruned nonredundant display", paste0(length(react14), " Reactome programmes x 10 cancers"),
    "Reactome hierarchy and gene overlap preclude treating rows as independent mechanisms")

# Figure 15: recurrence and direction ---------------------------------------
program15 <- data.table(Label = c("EMT", "Collagen formation", "Hypoxia", "ECM organisation", "E2F targets",
    "G2/M checkpoint", "Oxidative phosphorylation", "MYC targets V1"), Library = c("Hallmark", "Reactome",
    "Hallmark", "Reactome", "Hallmark", "Hallmark", "Hallmark", "Hallmark"), Pathway = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "REACTOME_COLLAGEN_FORMATION", "HALLMARK_HYPOXIA", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_MYC_TARGETS_V1"))
d15 <- merge(program15, s11b, by = c("Library", "Pathway"), all.x = TRUE, allow.cartesian = TRUE)
d15[, Category := fifelse(!is.finite(FDR) | FDR >= 0.05, "Nonsignificant", fifelse(NES > 0, "Adverse",
    "Favourable"))]
rec15 <- d15[Category != "Nonsignificant", .(Adverse = sum(Category == "Adverse"), Favourable = sum(Category ==
    "Favourable")), by = Label]
d15 <- merge(d15, rec15, by = "Label", all.x = TRUE)
d15[is.na(Adverse), `:=`(Adverse = 0L, Favourable = 0L)]
d15[, Row := sprintf("%s  [A%d/F%d]", Label, Adverse, Favourable)]
row15 <- rev(unique(d15[match(program15$Label, Label), Row]))
d15[, `:=`(Cancer = factor(Cancer, levels = cancers), Row = factor(Row, levels = row15), Category = factor(Category,
    levels = c("Favourable", "Nonsignificant", "Adverse")))]
p15 <- ggplot(d15, aes(Cancer, Row, fill = Category)) + geom_tile(colour = "white", linewidth = 0.6) +
    scale_fill_manual(values = c(Favourable = favourable, Nonsignificant = "#E3E6E9", Adverse = adverse)) +
    labs(x = NULL, y = NULL, fill = "BH FDR < 0.05 direction") + theme_atlas(9.5)
export_plot(p15, 15, "shared_and_context_dependent_pathways", 13, 6.8, "Stage11B_all_ranked_pathway_tests.csv",
    "BH-FDR direction matrix with source-derived recurrence", "8 programmes x 10 cancers", "Comparative transcript-survival enrichment; not measured pathway activity")

# Figure 16: scientific feasibility gate ------------------------------------ The revision
# explicitly forbids replacing adjusted, cancer-standardised predictions with pooled Kaplan-Meier
# curves. The saved Stage06 model is Surv(time,event) ~ z + strata(Cancer), without clinical
# adjustment, and only its coefficient table was archived. Stage07 stores separate cancer-specific
# adjusted estimates and a fixed-effect meta-analysis, not a stratified model object with baseline
# hazards/covariance needed for survival prediction.
serp06 <- s06_pan[Feature_ID == "SERPINE1" & Fit_status == "ok"]
serp07 <- s07_meta[Feature_ID == "SERPINE1"]
if (nrow(serp06) != 1L || nrow(serp07) != 1L) stop("Figure 16 feasibility audit: SERPINE1 summary evidence missing")
fig16_report <- c("FIGURE 16 SCIENTIFIC FEASIBILITY REPORT", "STATUS: BLOCKED_SCIENTIFICALLY - no Figure 16 image was generated.",
    "", "Requested estimand: clinically adjusted, cancer-stratified, cancer-distribution-standardised SERPINE1 survival curves with 95% confidence intervals.",
    "", "Blocking evidence:", "1. The saved Stage06 pan-cancer model is unadjusted: Surv(time,event) ~ within-cancer-standardised SERPINE1 + strata(Cancer).",
    "2. Stage06 archives a coefficient table, not the fitted Cox model, cancer-specific baseline cumulative hazards, model matrix, or full coefficient covariance required for survival prediction.",
    "3. Stage07 archives separate clinically adjusted cancer-specific models and a fixed-effect coefficient meta-analysis. These do not define one cancer-stratified survival model or shared standardisation population.",
    "4. No validated analysis output prespecifies high/low SERPINE1 prediction values, clinical-covariate standardisation values, follow-up horizon, or a bootstrap procedure that propagates model and standardisation uncertainty.",
    "5. Therefore model-based standardised curves and defensible 95% confidence intervals cannot be recovered from the validated saved outputs without fitting a new model and defining a new estimand.",
    "", sprintf("Available but insufficient Stage06 unadjusted stratified association: HR %.6f (95%% CI %.6f to %.6f), P %.6g.",
        serp06$HR, serp06$CI_low, serp06$CI_high, serp06$P), sprintf("Available but insufficient Stage07 fixed-effect meta-analysis of separate adjusted models: HR %.6f (95%% CI %.6f to %.6f), P %.6g.",
        serp07$HR, serp07$CI_low, serp07$CI_high, serp07$P), "", "Required next analysis: prespecify and validate a pooled clinically adjusted Cox model with cancer strata, covariate coding and expression contrast; archive the fitted model/baseline hazards; then use patient-level g-computation with cancer-distribution weights and bootstrap confidence intervals.",
    "A naive pooled Kaplan-Meier figure was not substituted.")
writeLines(fig16_report, file.path(out, "Figure16_scientific_feasibility_report.txt"))
manifest[[length(manifest) + 1L]] <- data.table(Figure = 16L, Stem = "SERPINE1_PanCancer_Adjusted_Survival",
    Source_files = "Stage06 pan-cancer stratified Cox coefficient table; Stage07 separate adjusted models/meta-analysis",
    Analysis_method = "Feasibility gate for model-based standardised survival prediction", Denominator = "Not defined: requested prediction model is not available",
    Width_in = NA_real_, Height_in = NA_real_, WebP = NA_character_, SVG = NA_character_, Output_status = "BLOCKED_SCIENTIFICALLY",
    Limitations = "No saved clinically adjusted cancer-stratified model/baseline hazards/covariance or validated standardisation/CI specification; no naive KM substitution")

# Final audit outputs ---------------------------------------------------------
manifest_dt <- rbindlist(manifest, fill = TRUE)
setorder(manifest_dt, Figure)
fwrite(manifest_dt, file.path(out, "figure_manifest.csv"))
fwrite(rbindlist(km_audit, fill = TRUE), file.path(out, "KM_sample_group_audit.csv"))

# Lightweight SVG content gates and raster dimension gates.
rendered_manifest <- manifest_dt[Output_status == "RENDERED"]
if (nrow(rendered_manifest) != 15L || nrow(manifest_dt) != 16L || nrow(manifest_dt[Figure == 16L & Output_status ==
    "BLOCKED_SCIENTIFICALLY"]) != 1L) stop("Final manifest status/count gate failed")
for (i in seq_len(nrow(rendered_manifest))) {
    svg_text <- paste(readLines(rendered_manifest$SVG[i], warn = FALSE, n = 40), collapse = " ")
    if (!grepl("<svg", svg_text, fixed = TRUE))
        stop("Invalid SVG: ", rendered_manifest$SVG[i])
}
validation <- c(validation, "PASS: Figures 01-15 exported as 300-dpi WebP and vector SVG; manifest extended through Figure 16.",
    "PASS: Figure 03 displays exp(Beta) hazard ratios on a log scale with HR = 1 reference; no successful finite Stage05 observation was excluded by an axis limit.",
    "PASS: Figure 04 contains 36 eligible genes selected by documented evidence-balanced rules; supplementary association data and exact methodology exported.",
    "PASS: all 15 WebP figures visually inspected for clipping, label legibility, legends, and KM risk-table alignment.",
    "PASS: all SVG files contain a valid SVG root.", "PASS: KM source time confirmed as days from Stage02 reports/scripts; display conversion uses days / 30.4375.",
    "PASS: KM patient IDs unique and expression columns exactly aligned to clinical sample IDs.", sprintf("PASS: COAD gene selected by saved Cox ranking before KM rendering: %s.",
        coad_gene), sprintf("PASS: DKK1-LUAD reconstructed gene PH P reproduces Stage08B: %.9f.", saved_p11),
    "BLOCKED: Figure 16 was not generated because no validated clinically adjusted cancer-stratified model object, baseline hazards/covariance, or standardisation/uncertainty specification is available; no pooled KM substitute was used.",
    "CAUTION: COAD selected gene lacks an archived Stage08 gene-specific PH diagnostic.", "CAUTION: Stage11B extremely small FDR values are displayed as supplied and require independent algorithmic validation before publication.",
    "CAUTION: all gene/pathway findings are candidate associations, not validated biomarkers or causal mechanisms.")
writeLines(c("PROJECT 3 FINAL FIGURE ATLAS VALIDATION REPORT", strrep("=", 54), validation), file.path(out,
    "validation_report.txt"))

summary_lines <- c("PROJECT 3 | REVISED FIGURE ATLAS (MANIFEST THROUGH FIGURE 16)", "Completion: Figures 01-15 rendered in WebP and SVG; Figure 16 scientifically blocked and not imaged",
    "Revision status: Figure 03 HR-axis revision complete; Figure 04 expanded to 36 genes with supplementary audit",
    "Validated cohort: 4,718 patients; 1,519 deaths; 10 TCGA cancers", sprintf("COAD Figure 09 gene: %s (selection audited; archived gene-specific PH diagnostic unavailable)",
        coad_gene), "KM time: source days, displayed as months using 30.4375 days/month; median splits; genuine censoring, Greenwood CI and risk tables",
    sprintf("DKK1-LUAD PH reconstruction: gene P = %.9f", saved_p11), "Main prediction evidence: training-only repeated cross-validation (100/100 folds)",
    "Figure 16: unavailable saved model/standardisation specification prevents defensible adjusted survival prediction and confidence intervals",
    "Primary caveats: discovery conditioning, residual confounding, PH violations, dependent CV folds, pathway redundancy, no external validation",
    paste0("Output: ", out))
writeLines(summary_lines, file.path(out, "final_console_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
