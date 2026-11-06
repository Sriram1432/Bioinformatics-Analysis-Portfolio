# Fit prioritised clinically adjusted models and evaluate model assumptions.  Reports
# matched-sample estimates, proportional-hazards diagnostics and influence sensitivity without
# presenting selected candidates as validated biomarkers.  Run with: Rscript
# 05_adjusted_models_diagnostics.R

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


# ---------------------------------------------------------------------------- Clinically adjusted
# candidate models ----------------------------------------------------------------------------

options(width = 190, scipen = 5)
root <- PROJECT_ROOT
in03 <- file.path(root, "Results_Rebuild", "Stage03")
in05 <- file.path(root, "Results_Rebuild", "Stage05")
in06 <- file.path(root, "Results_Rebuild", "Stage06")
out <- file.path(root, "Results_Rebuild", "Stage07")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (pkg in c("data.table", "survival", "matrixStats", "ggplot2", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Install required package: ", pkg)
}
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")

# ---------- 1. Read prior discovery outputs; do not repeat genome-wide Cox. ---
reqfile <- function(path) {
    if (!file.exists(path))
        stop("Required input absent: ", path)
    path
}
int <- fread(reqfile(file.path(in06, "Stage06_integrated_pan_cancer_evidence.csv")))
needed_int <- c("Feature_ID", "FDR", "P", "I2_pct", "Q_FDR", "N_significant_cancers")
if (!all(needed_int %in% names(int))) stop("Stage06 integrated schema differs: ", paste(setdiff(needed_int,
    names(int)), collapse = ", "))
prior <- rbindlist(lapply(cancers, function(cc) {
    x <- fread(reqfile(file.path(in05, paste0("Stage05_", cc, "_all_gene_Cox.csv"))))
    needed <- c("Cancer", "Feature_ID", "P", "FDR", "Beta", "SE", "Fit_status")
    if (!all(needed %in% names(x)))
        stop("Stage05 schema differs in ", cc)
    x[, ..needed]
}))

# Prioritisation is deliberately modest, prespecified from AVAILABLE discoveries.  High and low
# heterogeneity both matter; recurrence includes opposite effects.
setorder(int, FDR, P, na.last = TRUE)
lowhet <- int[is.finite(FDR) & FDR < 0.05 & is.finite(I2_pct) & I2_pct < 25, head(Feature_ID, 12)]
highhet <- int[is.finite(FDR) & FDR < 0.05 & is.finite(I2_pct) & I2_pct >= 75, head(Feature_ID, 12)]
recur <- int[is.finite(N_significant_cancers) & N_significant_cancers >= 4][order(-N_significant_cancers,
    FDR), head(Feature_ID, 10)]
# Strongest 5 FDR hits per malignancy, including those without any FDR hits:
per_cancer <- unique(unlist(lapply(cancers, function(cc) {
    z <- prior[Cancer == cc & Fit_status == "ok" & is.finite(FDR)][order(FDR, P)]
    head(z$Feature_ID, 5)
}), use.names = FALSE))
# Explicit biologically interesting representatives from reported outputs.
a_priori_examples <- c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1")
priority <- unique(c(lowhet, highhet, recur, per_cancer, a_priori_examples))
priority <- priority[!is.na(priority) & nzchar(priority)]
selection <- data.table(Feature_ID = priority, Low_heterogeneity = priority %chin% lowhet, High_heterogeneity = priority %chin%
    highhet, Recurrent = priority %chin% recur, Top_cancer_specific = priority %chin% per_cancer, Named_example = priority %chin%
    a_priori_examples)
fwrite(selection, file.path(out, "Stage07_priority_feature_selection.csv"))

# ---------- 2. Harmonise clinical covariates (no outcome-driven selection). ---
clean_chr <- function(x) {
    z <- trimws(as.character(x))
    z[toupper(z) %in% c("", "NA", "N/A", "NULL", "NOT REPORTED", "NOT AVAILABLE", "UNKNOWN", "[NOT AVAILABLE]",
        "[NOT EVALUATED]", "[DISCREPANCY]")] <- NA_character_
    z
}
major_stage <- function(x) {
    z <- toupper(clean_chr(x))
    z <- gsub("^STAGE[[:space:]]*", "", z)
    out <- rep(NA_character_, length(z))
    out[grepl("^IV", z)] <- "IV"
    out[grepl("^III", z)] <- "III"
    out[grepl("^II", z) & is.na(out)] <- "II"
    out[grepl("^I($|A|B|C|[[:space:]])", z) & is.na(out)] <- "I"
    out
}
prepare_clinical <- function(cl, cc) {
    n <- nrow(cl)
    d <- data.frame(time = as.numeric(cl$OS_time), event = as.integer(cl$OS))
    # Age in years; do not silently impute absent covariates.
    if ("age_at_initial_pathologic_diagnosis" %in% names(cl)) {
        age <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
        age[!is.finite(age) | age < 18 | age > 105] <- NA_real_
        d$age10 <- age/10
    }
    if ("gender" %in% names(cl)) {
        sx <- toupper(clean_chr(cl$gender))
        sx[!sx %in% c("MALE", "FEMALE")] <- NA_character_
        if (sum(!is.na(sx)) >= 0.8 * n && min(table(sx), na.rm = TRUE) >= 10)
            d$sex <- factor(sx)
    }
    # Use AJCC pathologic stage where available, otherwise clinical_stage.
    a <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        major_stage(cl$ajcc_pathologic_tumor_stage) else rep(NA_character_, n)
    b <- if ("clinical_stage" %in% names(cl))
        major_stage(cl$clinical_stage) else rep(NA_character_, n)
    stage <- if (sum(!is.na(a)) >= sum(!is.na(b)))
        a else b
    if (sum(!is.na(stage)) >= 0.8 * n && length(unique(na.omit(stage))) >= 2) {
        d$stage <- factor(stage, levels = c("I", "II", "III", "IV"))
    }
    # UCEC histology has meaningful within-cancer subtype variation. Group serous and mixed-serous
    # vs endometrioid; omit only if sparsity/coverage precludes it.
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- toupper(clean_chr(cl$histological_type))
        group <- rep(NA_character_, n)
        group[grepl("ENDOMETRIOID", h)] <- "Endometrioid"
        group[grepl("SEROUS|MIXED", h)] <- "Serous_or_mixed"
        if (sum(!is.na(group)) >= 0.8 * n && all(table(group) >= 10))
            d$histology <- factor(group)
    }
    d
}

# ---------- 3. Fit matched-sample unadjusted and adjusted models. ----------
run_fit <- function(dat, adjust = TRUE) {
    terms <- if (adjust)
        setdiff(names(dat), c("time", "event", "z")) else character()
    f <- as.formula(paste("survival::Surv(time, event) ~ z", if (length(terms))
        paste0(" + ", paste(terms, collapse = " + ")) else ""))
    warnings <- character()
    mod <- withCallingHandlers(tryCatch(survival::coxph(f, data = dat, ties = "efron", control = survival::coxph.control(iter.max = 40),
        x = FALSE, y = FALSE), error = function(e) e), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (inherits(mod, "error"))
        return(list(beta = NA_real_, se = NA_real_, p = NA_real_, status = paste0("error: ", conditionMessage(mod)),
            ph = NA_real_))
    bet <- unname(coef(mod)["z"])
    se <- tryCatch(sqrt(unname(vcov(mod)["z", "z"])), error = function(e) NA_real_)
    ok <- is.finite(bet) && is.finite(se) && se > 0 && abs(bet) < 15
    status <- if (!ok)
        "invalid_estimate" else if (length(warnings))
        "fit_warning" else "ok"
    p <- if (ok)
        2 * pnorm(-abs(bet/se)) else NA_real_
    list(beta = if (ok) bet else NA_real_, se = if (ok) se else NA_real_, p = p, status = status, ph = NA_real_)
}

messages <- list()
result <- list()
summary <- list()
ph_results <- list()
for (cc in cancers) {
    message("Stage07: clinical adjustment ", cc)
    d <- readRDS(reqfile(file.path(in03, paste0("TCGA_", cc, "_discovery.rds"))))
    x <- d$expression
    cl <- as.data.table(d$clinical)
    if (!is.matrix(x) || !is.numeric(x) || !identical(colnames(x), as.character(cl$sample_id)))
        stop("Expression/clinical mismatch: ", cc)
    if (anyNA(cl$OS) || anyNA(cl$OS_time) || any(cl$OS_time <= 0))
        stop("Invalid survival: ", cc)
    clinical <- prepare_clinical(cl, cc)
    if (!any(c("age10", "stage", "histology", "sex") %in% names(clinical)))
        warning("No adjustment variables in ", cc)
    terms <- setdiff(names(clinical), c("time", "event"))
    # Use complete cases for the clinical baseline, equally for both gene fits.
    include <- complete.cases(clinical)
    model_cl <- clinical[include, , drop = FALSE]
    # Drop unused factor levels to avoid singular dummy terms.
    model_cl[] <- lapply(model_cl, function(z) if (is.factor(z))
        droplevels(z) else z)
    terms <- Filter(function(k) length(unique(model_cl[[k]])) >= 2L, terms)
    model_cl <- model_cl[, c("time", "event", terms), drop = FALSE]
    if (nrow(model_cl) < 30L || sum(model_cl$event) < 20L)
        stop("Insufficient clinical complete cases/events for ", cc)
    baseline_warnings <- character()
    baseline <- withCallingHandlers(tryCatch(survival::coxph(as.formula(paste("survival::Surv(time,event)~",
        if (length(terms))
            paste(terms, collapse = "+") else "1")), data = model_cl, ties = "efron"), error = function(e) e), warning = function(w) {
        baseline_warnings <<- c(baseline_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    baseline_status <- if (inherits(baseline, "error"))
        conditionMessage(baseline) else "ok"
    # No silent use of a failed clinical baseline.
    if (inherits(baseline, "error"))
        stop("Clinical baseline cannot be fitted for ", cc, ": ", baseline_status)
    # Note coefficients & event information rather than assuming all factors are stable.
    if (length(baseline_warnings))
        warning("Clinical baseline warning in ", cc, ": ", paste(unique(baseline_warnings), collapse = "; "))
    if (any(!is.finite(coef(baseline))))
        stop("Nonfinite clinical baseline coefficient for ", cc)
    original <- as.integer(sum(cl$OS == 1))
    summary[[cc]] <- data.table(Cancer = cc, N_total = nrow(cl), Events_total = original, N_adjusted = nrow(model_cl),
        Events_adjusted = sum(model_cl$event), Covariates = if (length(terms))
            paste(terms, collapse = " + ") else "none", Baseline_status = baseline_status, Candidate_features = sum(priority %chin% rownames(x)))

    for (feature in priority[priority %chin% rownames(x)]) {
        raw <- as.numeric(x[feature, ])
        ss <- sd(raw)
        if (!is.finite(ss) || ss < 1e-10 || any(!is.finite(raw)))
            next
        z <- (raw - mean(raw))/ss  # global cancer-wise SD, preserved on complete-case subset
        dat <- model_cl
        dat$z <- z[include]
        unadjusted <- run_fit(dat, adjust = FALSE)
        adjusted <- run_fit(dat, adjust = TRUE)
        result[[length(result) + 1L]] <- data.table(Cancer = cc, Feature_ID = feature, N = nrow(dat),
            Events = sum(dat$event), Covariates = if (length(terms))
                paste(terms, collapse = "+") else "none", Unadjusted_HR = exp(unadjusted$beta), Unadjusted_P = unadjusted$p, Adjusted_Beta = adjusted$beta,
            Adjusted_SE = adjusted$se, Adjusted_HR = exp(adjusted$beta), CI_low = exp(adjusted$beta -
                1.96 * adjusted$se), CI_high = exp(adjusted$beta + 1.96 * adjusted$se), Adjusted_P = adjusted$p,
            Unadjusted_status = unadjusted$status, Adjusted_status = adjusted$status)
    }
    rm(d, x, cl)
    invisible(gc(verbose = FALSE))
}
res <- rbindlist(result, use.names = TRUE, fill = TRUE)
# FDR across ALL prespecified priority-feature \303\227 cancer adjusted tests.
res[, Adjusted_FDR_exploratory := p.adjust(Adjusted_P, method = "BH", n = .N)]
# Also combine adjusted estimates by cancer (effect heterogeneity may remain).
valid <- res[Adjusted_status == "ok" & is.finite(Adjusted_Beta) & is.finite(Adjusted_SE) & Adjusted_SE >
    0]
meta <- valid[, {
    w <- 1/Adjusted_SE^2
    k <- .N
    b <- sum(w * Adjusted_Beta)/sum(w)
    se <- sqrt(1/sum(w))
    q <- sum(w * (Adjusted_Beta - b)^2)
    .(K = k, HR = exp(b), CI_low = exp(b - 1.96 * se), CI_high = exp(b + 1.96 * se), P = 2 * pnorm(-abs(b/se)),
        Q = q, Q_P = if (k > 1)
            pchisq(q, df = k - 1, lower.tail = FALSE) else NA_real_, I2_pct = if (k > 1 && q > 0)
            max(0, (q - (k - 1))/q) * 100 else NA_real_, N_adverse = sum(Adjusted_Beta > 0), N_favourable = sum(Adjusted_Beta < 0))
}, by = Feature_ID]
# DerSimonian-Laird random-effects sensitivity: descriptive when K is small; not a substitute for a
# multilevel hierarchical model or individual-level adjusted pooling.
random_meta <- valid[, {
    k <- .N
    w <- 1/Adjusted_SE^2
    b0 <- sum(w * Adjusted_Beta)/sum(w)
    q <- sum(w * (Adjusted_Beta - b0)^2)
    cval <- sum(w) - sum(w^2)/sum(w)
    tau2 <- if (k >= 2 && is.finite(cval) && cval > 0)
        max(0, (q - (k - 1))/cval) else NA_real_
    wr <- 1/(Adjusted_SE^2 + if (is.finite(tau2))
        tau2 else 0)
    br <- sum(wr * Adjusted_Beta)/sum(wr)
    ser <- sqrt(1/sum(wr))
    .(Tau2_DL = tau2, Random_HR = exp(br), Random_CI_low = exp(br - 1.96 * ser), Random_CI_high = exp(br +
        1.96 * ser), Random_P = 2 * pnorm(-abs(br/ser)))
}, by = Feature_ID]
meta <- merge(meta, random_meta, by = "Feature_ID", all.x = TRUE, sort = FALSE)
meta[, Random_FDR_exploratory := p.adjust(Random_P, method = "BH", n = length(priority))]
meta[, Exploratory_FDR := p.adjust(P, method = "BH", n = length(priority))]
setorder(meta, Exploratory_FDR, P, na.last = TRUE)
setorder(res, Cancer, Adjusted_P, na.last = TRUE)
clinical_summary <- rbindlist(summary)
fwrite(res, file.path(out, "Stage07_all_adjusted_candidate_models.csv"))
fwrite(meta, file.path(out, "Stage07_adjusted_cross_cancer_fixed_effect_meta.csv"))
fwrite(clinical_summary, file.path(out, "Stage07_clinical_model_cohorts.csv"))

per_cancer_summary <- res[, .(Models = .N, Successful = sum(Adjusted_status == "ok"), Adjusted_P_005 = sum(Adjusted_P <
    0.05, na.rm = TRUE), Exploratory_FDR_005 = sum(Adjusted_FDR_exploratory < 0.05, na.rm = TRUE)), by = Cancer]
status_counts <- res[, .N, by = .(Cancer, Adjusted_status)]
top_adjusted <- res[is.finite(Adjusted_P)][order(Adjusted_P)][seq_len(min(25L, .N)), .(Cancer, Feature_ID,
    N, Events, Adjusted_HR, CI_low, CI_high, Adjusted_P, Adjusted_FDR_exploratory)]
lead <- res[Feature_ID %chin% c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1"),
    .(Cancer, Feature_ID, Adjusted_HR, CI_low, CI_high, Adjusted_P, Adjusted_status)][order(Feature_ID,
    Cancer)]
top_meta <- meta[seq_len(min(20L, .N)), .(Feature_ID, K, HR, CI_low, CI_high, P, Exploratory_FDR, I2_pct,
    Q_P, N_adverse, N_favourable)]
# ---------- 4. Editorial-quality, standalone diagnostic association figures. --- All graphics are
# descriptive: selected using the same survival outcomes.
black <- "#171717"
muted <- "#5A6573"
base_theme <- function(size = 13) ggplot2::theme_minimal(base_size = size, base_family = "sans") + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = black, size = 16), axis.title = ggplot2::element_text(face = "bold", colour = black), axis.text = ggplot2::element_text(colour = black),
    legend.title = ggplot2::element_text(face = "bold", colour = black), panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(15, 22, 15, 15))
export_plot <- function(p, file, width = 12, height = 8) {
    ggplot2::ggsave(file, p, device = ragg::agg_webp, width = width, height = height, units = "in", dpi = 300,
        bg = "white")
}
figs <- character()
# Clinical inclusion proportions. Disclosure: covariate selection based on availability.
cs <- copy(clinical_summary)
cs[, Fraction := N_adjusted/N_total]
p0 <- ggplot2::ggplot(cs, ggplot2::aes(x = Fraction, y = reorder(Cancer, Fraction))) + ggplot2::geom_segment(ggplot2::aes(x = 0,
    xend = Fraction, yend = Cancer), linewidth = 1.1, colour = "#D7DDE7") + ggplot2::geom_point(size = 5,
    colour = "#164E77") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%d / %d", N_adjusted, N_total)),
    hjust = -0.15, size = 4, colour = black) + ggplot2::scale_x_continuous(labels = function(z) sprintf("%d%%",
    round(z * 100)), limits = c(0, 1.23), breaks = c(0, 0.25, 0.5, 0.75, 1)) + ggplot2::labs(title = "Clinical adjustment: available patient cohorts",
    subtitle = "Complete-case inclusion; event counts and covariates reported in accompanying table",
    x = "Patients retained for clinical models", y = NULL) + base_theme()
f0 <- file.path(figdir, "Stage07_clinical_adjustment_cohort_retention.webp")
export_plot(p0, f0, 11, 7)
figs <- c(figs, f0)
# Per-cancer ranked forest plots, select top 10 by adjusted P, report FDR exploratory.
for (cc in cancers) {
    z <- copy(res[Cancer == cc & Adjusted_status == "ok" & is.finite(Adjusted_HR) & is.finite(CI_low) &
        is.finite(CI_high) & CI_low > 0 & CI_high > 0])
    if (!nrow(z))
        next
    setorder(z, Adjusted_P)
    z <- head(z, 10L)
    z[, `:=`(Gene = factor(Feature_ID, levels = rev(Feature_ID)), Direction = ifelse(Adjusted_HR > 1,
        "Higher hazard", "Lower hazard"))]
    p <- ggplot2::ggplot(z, ggplot2::aes(x = Adjusted_HR, y = Gene, colour = Direction)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = muted) + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high),
        orientation = "y", width = 0.2, linewidth = 0.8) + ggplot2::geom_point(size = 3) + ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(values = c(`Higher hazard` = "#B34039", `Lower hazard` = "#237A74")) +
        ggplot2::labs(title = paste0(cc, " | clinically adjusted candidate associations"), subtitle = "Exploratory, outcome-selected genes; HR per within-cancer SD; 95% Wald CI",
            x = "Adjusted hazard ratio (log scale)", y = NULL, colour = NULL) + base_theme(12) + ggplot2::theme(legend.position = "bottom")
    f <- file.path(figdir, paste0("Stage07_", cc, "_adjusted_top10_forest_95CI.webp"))
    export_plot(p, f, 11, 7.2)
    figs <- c(figs, f)
}
# One multi-cancer gene forest per biologically motivated gene with estimates.
for (g in c("SERPINE1", "STX1A", "MET", "DKK1", "PGK1", "PSMD2")) {
    z <- copy(res[Feature_ID == g & Adjusted_status == "ok" & is.finite(Adjusted_HR) & is.finite(CI_low) &
        is.finite(CI_high) & CI_low > 0 & CI_high > 0])
    if (nrow(z) < 3)
        next
    z[, Cancer := factor(Cancer, levels = rev(cancers))]
    z[, Direction := ifelse(Adjusted_HR > 1, "Higher hazard", "Lower hazard")]
    p <- ggplot2::ggplot(z, ggplot2::aes(x = Adjusted_HR, y = Cancer, colour = Direction)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = muted) + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high),
        orientation = "y", width = 0.19, linewidth = 0.8) + ggplot2::geom_point(size = 3.4) + ggplot2::scale_x_log10() +
        ggplot2::scale_colour_manual(values = c(`Higher hazard` = "#B34039", `Lower hazard` = "#237A74")) +
        ggplot2::labs(title = paste0(g, " | adjusted association across cancers"), subtitle = "Different cancer-specific covariates; 95% CI; exploratory; no assumption of shared effects",
            x = "Cancer-specific adjusted HR per 1 SD (log scale)", y = NULL, colour = NULL) + base_theme(12) +
        ggplot2::theme(legend.position = "bottom")
    f <- file.path(figdir, paste0("Stage07_", g, "_cross_cancer_adjusted_forest_95CI.webp"))
    export_plot(p, f, 11, 7.5)
    figs <- c(figs, f)
}
# Fixed vs random model comparison: label both as exploratory sensitivity estimates.
mm <- copy(meta[K >= 5 & is.finite(Random_HR) & is.finite(HR) & is.finite(I2_pct)])
if (nrow(mm)) {
    setorder(mm, Exploratory_FDR, P)
    mm <- head(mm, 18L)
    u <- rbindlist(list(mm[, .(Feature_ID, Method = "Fixed effect", HR, low = CI_low, high = CI_high)],
        mm[, .(Feature_ID, Method = "Random effects (DL)", HR = Random_HR, low = Random_CI_low, high = Random_CI_high)]))
    u[, Feature_ID := factor(Feature_ID, levels = rev(mm$Feature_ID))]
    p <- ggplot2::ggplot(u, ggplot2::aes(x = HR, y = Feature_ID, colour = Method)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = muted) + ggplot2::geom_errorbar(ggplot2::aes(xmin = low, xmax = high),
        orientation = "y", width = 0.22, position = ggplot2::position_dodge(width = 0.65)) + ggplot2::geom_point(size = 2.7,
        position = ggplot2::position_dodge(width = 0.65)) + ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = c(`Fixed effect` = "#164E77",
        `Random effects (DL)` = "#CC7440")) + ggplot2::labs(title = "Adjusted cross-cancer effects: model sensitivity",
        subtitle = "Discovery-conditioned candidates; fixed vs DerSimonian\342\200\223Laird random effects; 95% CI",
        x = "Pooled HR per within-cancer SD (log scale)", y = NULL, colour = NULL) + base_theme(11) +
        ggplot2::theme(legend.position = "bottom")
    f <- file.path(figdir, "Stage07_adjusted_fixed_vs_random_meta_forest_95CI.webp")
    export_plot(p, f, 12, 9)
    figs <- c(figs, f)
}
fwrite(data.table(Figure_path = figs), file.path(out, "Stage07_figure_index.csv"))

report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 07 \342\200\224 CONSOLIDATED RESULTS\n====================================================================\n")
    cat("1. CANDIDATE SELECTION (DISCOVERY-CONDITIONED; EXPLORATORY)\n")
    print(data.table(Candidates = length(priority), Low_heterogeneity = length(lowhet), High_heterogeneity = length(highhet),
        Recurrent = length(recur), Top_cancer_specific = length(per_cancer)))
    cat("\n2. CLINICAL ADJUSTMENT DESIGN / COMPLETE CASES\n")
    print(clinical_summary)
    cat("\n3. ADJUSTED ASSOCIATIONS BY CANCER (PRIORITISED ONLY)\n")
    print(per_cancer_summary)
    cat("\n4. MODEL STATUS\n")
    print(status_counts)
    cat("\n5. TOP 25 CLINICALLY ADJUSTED CANDIDATE ASSOCIATIONS\n")
    print(top_adjusted)
    cat("\n6. SELECTED BIOLOGICAL EXAMPLES ACROSS CANCERS\n")
    print(lead)
    cat("\n7. ADJUSTED CROSS-CANCER META-ANALYSIS (TOP 20)\n")
    print(top_meta)
    cat("\n8. INTERPRETATION AND NEXT DECISION\n")
    cat("All tests are exploratory: candidate selection used the same TCGA outcome data.\n")
    cat("Adjusted-FDR is only across the selected family, NOT genome-wide confirmation.\n")
    cat("Within-cancer gene-wise z-scores from full original cohort; all comparisons use identical clinical complete cases within cancer.\n")
    cat("Covariates selected by availability/biological relevance, not outcome. Stage/sex/subtype harmonisation is approximate.\n")
    cat("Fixed-effect meta assumes shared log-HR; random-effects DL is a descriptive sensitivity analysis, not confirmation.\n")
    cat("Stage07 FDR is selected-panel exploratory only; no independent validation.\n")
    cat("Do not call associations independent biomarkers until model diagnostics, robustness and validation.\n")
    cat("\n9. SCIENTIFIC FIGURES (WEBP)\n")
    print(figs)
    cat("\n10. OUTPUT DIRECTORY\n", out, "\n====================================================================\n")
})
writeLines(report, file.path(out, "Stage07_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Cox diagnostics and
# influence analysis ----------------------------------------------------------------------------

options(width = 190, scipen = 5)
root <- PROJECT_ROOT
in03 <- file.path(root, "Results_Rebuild", "Stage03")
in07 <- file.path(root, "Results_Rebuild", "Stage07")
out <- file.path(root, "Results_Rebuild", "Stage08")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "survival", "ggplot2", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Missing package: ", pkg, ". Install with install.packages('", pkg, "')")
}
suppressPackageStartupMessages(library(data.table))
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
reqfile <- function(p) {
    if (!file.exists(p))
        stop("Required input missing: ", p)
    p
}
prior <- fread(reqfile(file.path(in07, "Stage07_all_adjusted_candidate_models.csv")))
meta <- fread(reqfile(file.path(in07, "Stage07_adjusted_cross_cancer_fixed_effect_meta.csv")))
stopifnot(all(c("Cancer", "Feature_ID", "Adjusted_status", "Adjusted_P") %in% names(prior)))
stopifnot(all(c("Feature_ID", "P", "K") %in% names(meta)))
# Fixed biological exemplars + strongest additional adjusted meta candidates.
anchors <- c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1", "DKK1", "DONSON")
meta_candidates <- meta[is.finite(P) & K >= 8][order(P)]$Feature_ID
extras <- head(setdiff(meta_candidates, anchors), 5)
features <- unique(c(anchors, extras))
selection <- data.table(Feature_ID = features, Selection = ifelse(features %in% anchors, "Prespecified_example",
    "Top_adjusted_meta_exploratory"))
fwrite(selection, file.path(out, "Stage08_selected_features.csv"))

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
prepare_clinical <- function(cl, cc) {
    n <- nrow(cl)
    d <- data.frame(time = as.numeric(cl$OS_time), event = as.integer(cl$OS))
    if ("age_at_initial_pathologic_diagnosis" %in% names(cl)) {
        age <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
        age[!is.finite(age) | age < 18 | age > 105] <- NA_real_
        d$age10 <- age/10
    }
    if ("gender" %in% names(cl)) {
        sx <- toupper(clean_chr(cl$gender))
        sx[!sx %in% c("MALE", "FEMALE")] <- NA_character_
        if (sum(!is.na(sx)) >= 0.8 * n && length(table(sx)) > 1 && min(table(sx)) >= 10)
            d$sex <- factor(sx)
    }
    a <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        major_stage(cl$ajcc_pathologic_tumor_stage) else rep(NA_character_, n)
    b <- if ("clinical_stage" %in% names(cl))
        major_stage(cl$clinical_stage) else rep(NA_character_, n)
    st <- if (sum(!is.na(a)) >= sum(!is.na(b)))
        a else b
    if (sum(!is.na(st)) >= 0.8 * n && length(unique(na.omit(st))) >= 2)
        d$stage <- factor(st, levels = c("I", "II", "III", "IV"))
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- toupper(clean_chr(cl$histological_type))
        group <- rep(NA_character_, n)
        group[grepl("ENDOMETRIOID", h)] <- "Endometrioid"
        group[grepl("SEROUS|MIXED", h)] <- "Serous_or_mixed"
        if (sum(!is.na(group)) >= 0.8 * n && length(table(group)) > 1 && all(table(group) >= 10))
            d$histology <- factor(group)
    }
    d
}

# Reconstruct Stage07 clinical specification exactly; do not redo genome-wide discovery.
rows <- list()
errors <- list()
ix <- 0L
for (cc in cancers) {
    message("Stage08: diagnostics for ", cc)
    dat <- readRDS(reqfile(file.path(in03, paste0("TCGA_", cc, "_discovery.rds"))))
    x <- dat$expression
    cl <- as.data.table(dat$clinical)
    if (!is.matrix(x) || !is.numeric(x) || !identical(colnames(x), as.character(cl$sample_id)))
        stop("Invalid alignment in ", cc)
    model_cl <- prepare_clinical(cl, cc)
    include <- complete.cases(model_cl)
    model_cl <- model_cl[include, , drop = FALSE]
    model_cl[] <- lapply(model_cl, function(z) if (is.factor(z))
        droplevels(z) else z)
    terms <- setdiff(names(model_cl), c("time", "event"))
    terms <- Filter(function(k) length(unique(model_cl[[k]])) >= 2L, terms)
    model_cl <- model_cl[, c("time", "event", terms), drop = FALSE]
    base_formula <- as.formula(paste("survival::Surv(time,event) ~", if (length(terms))
        paste(terms, collapse = " + ") else "1"))
    gene_formula <- as.formula(paste("survival::Surv(time,event) ~ z", if (length(terms))
        paste0(" + ", paste(terms, collapse = " + ")) else ""))
    # Baseline is fixed for all genes in this cancer.
    baseline <- tryCatch(survival::coxph(base_formula, data = model_cl, ties = "efron", x = TRUE, y = TRUE),
        error = function(e) e)
    if (inherits(baseline, "error") || any(!is.finite(coef(baseline))))
        stop("Clinical baseline unsuccessful for ", cc)
    for (feature in features) {
        if (!(feature %in% rownames(x)))
            next
        if (!any(prior$Cancer == cc & prior$Feature_ID == feature & prior$Adjusted_status == "ok"))
            next
        raw <- as.numeric(x[feature, ])
        sd_raw <- sd(raw)
        if (!is.finite(sd_raw) || sd_raw < 1e-10 || anyNA(raw))
            next
        df <- model_cl
        df$z <- ((raw - mean(raw))/sd_raw)[include]
        warns <- character()
        fit <- withCallingHandlers(tryCatch(survival::coxph(gene_formula, data = df, ties = "efron",
            x = TRUE, y = TRUE, control = survival::coxph.control(iter.max = 40)), error = function(e) e),
            warning = function(w) {
                warns <<- c(warns, conditionMessage(w))
                invokeRestart("muffleWarning")
            })
        if (inherits(fit, "error") || length(warns) > 0 || any(!is.finite(coef(fit)))) {
            errors[[length(errors) + 1L]] <- data.table(Cancer = cc, Feature_ID = feature, Problem = if (inherits(fit,
                "error"))
                conditionMessage(fit) else paste(warns, collapse = " | "))
            next
        }
        beta <- unname(coef(fit)["z"])
        se <- sqrt(unname(vcov(fit)["z", "z"]))
        if (!is.finite(beta) || !is.finite(se) || se <= 0)
            next
        ph <- tryCatch(survival::cox.zph(fit, transform = "km"), error = function(e) e)
        ph_gene <- if (!inherits(ph, "error") && "z" %in% rownames(ph$table))
            unname(ph$table["z", "p"]) else NA_real_
        ph_global <- if (!inherits(ph, "error") && "GLOBAL" %in% rownames(ph$table))
            unname(ph$table["GLOBAL", "p"]) else NA_real_
        # DFBETA captures each individual's approximate influence on the z coefficient.
        db <- tryCatch(residuals(fit, type = "dfbeta"), error = function(e) e)
        # residuals.coxph(type='dfbeta') may return an unnamed matrix.  Its columns follow
        # coef(fit) order; select gene by coefficient INDEX, not by residual matrix column name.
        z_index <- match("z", names(coef(fit)))
        if (inherits(db, "error") || is.na(z_index)) {
            max_dfbeta <- NA_real_
            influential_n <- NA_integer_
            drop <- integer()
        } else {
            if (is.null(dim(db)))
                db <- matrix(db, ncol = length(coef(fit)))
            if (ncol(db) != length(coef(fit)) || nrow(db) != nrow(df)) {
                max_dfbeta <- NA_real_
                influential_n <- NA_integer_
                drop <- integer()
                errors[[length(errors) + 1L]] <- data.table(Cancer = cc, Feature_ID = feature, Problem = "DFBETA dimensions do not match model")
            } else {
                zdb <- abs(as.numeric(db[, z_index]))
                max_dfbeta <- if (any(is.finite(zdb)))
                  max(zdb[is.finite(zdb)]) else NA_real_
                influential_n <- sum(zdb > 2/sqrt(sum(df$event)), na.rm = TRUE)
                # Sensitivity: remove the top 1% most influential observations.
                n_drop <- min(sum(is.finite(zdb)), max(1L, ceiling(0.01 * nrow(df))))
                drop <- if (n_drop)
                  head(order(zdb, decreasing = TRUE, na.last = NA), n_drop) else integer()
            }
        }
        sensitivity_beta <- NA_real_
        sensitivity_p <- NA_real_
        if (length(drop) > 0 && sum(df$event[-drop]) >= 20) {
            sf <- tryCatch(suppressWarnings(survival::coxph(gene_formula, data = df[-drop, , drop = FALSE],
                ties = "efron", control = survival::coxph.control(iter.max = 40))), error = function(e) e)
            if (!inherits(sf, "error") && is.finite(unname(coef(sf)["z"]))) {
                sensitivity_beta <- unname(coef(sf)["z"])
                sx <- tryCatch(summary(sf)$coefficients["z", "Pr(>|z|)"], error = function(e) NA_real_)
                sensitivity_p <- as.numeric(sx)
            }
        }
        # Nested 1-df partial likelihood ratio test and APPARENT concordance.
        lrt <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(baseline)))
        lrt_p <- if (is.finite(lrt))
            pchisq(max(lrt, 0), df = 1, lower.tail = FALSE) else NA_real_
        c_base <- unname(summary(baseline)$concordance[1])
        c_gene <- unname(summary(fit)$concordance[1])
        ix <- ix + 1L
        rows[[ix]] <- data.table(Cancer = cc, Feature_ID = feature, N = nrow(df), Events = sum(df$event),
            Adjusted_HR = exp(beta), CI_low = exp(beta - 1.96 * se), CI_high = exp(beta + 1.96 * se),
            Gene_PH_P = ph_gene, Global_PH_P = ph_global, Max_abs_DFBETA = max_dfbeta, Influential_patients = influential_n,
            Excluded_in_sensitivity = length(drop), Sensitivity_HR = exp(sensitivity_beta), Sensitivity_P = sensitivity_p,
            LogHR_change = if (is.finite(sensitivity_beta))
                sensitivity_beta - beta else NA_real_, LRT_P = lrt_p, Apparent_C_clinical = c_base, Apparent_C_plus_gene = c_gene,
            Apparent_Delta_C = c_gene - c_base)
    }
    rm(dat, x, cl, model_cl, baseline)
    invisible(gc(verbose = FALSE))
}
res <- if (length(rows)) rbindlist(rows, fill = TRUE) else stop("No diagnostics completed. Check Stage07 inputs.")
res[, PH_FDR_exploratory := p.adjust(Gene_PH_P, method = "BH")]
res[, LRT_FDR_exploratory := p.adjust(LRT_P, method = "BH")]
res[, HR_change_pct := 100 * (exp(LogHR_change) - 1)]
setorder(res, Cancer, Feature_ID)
fwrite(res, file.path(out, "Stage08_all_candidate_diagnostics.csv"))
err <- if (length(errors)) rbindlist(errors, fill = TRUE) else data.table(Cancer = character(), Feature_ID = character(),
    Problem = character())
fwrite(err, file.path(out, "Stage08_model_problems.csv"))
summary_cancer <- res[, .(Models = .N, Gene_PH_P_005 = sum(Gene_PH_P < 0.05, na.rm = TRUE), Global_PH_P_005 = sum(Global_PH_P <
    0.05, na.rm = TRUE), LRT_P_005 = sum(LRT_P < 0.05, na.rm = TRUE), Median_apparent_Delta_C = median(Apparent_Delta_C,
    na.rm = TRUE), Sensitivity_change_gt_20pct = sum(abs(HR_change_pct) > 20, na.rm = TRUE)), by = Cancer]

# ----- Figure atlas: single-panel, publication-oriented diagnostics. -----
ink <- "#17191E"
base_theme <- function() ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = ink, size = 17), plot.subtitle = ggplot2::element_text(colour = "#4C5663", size = 11), axis.title = ggplot2::element_text(face = "bold",
    colour = ink), axis.text = ggplot2::element_text(colour = ink), legend.title = ggplot2::element_text(face = "bold",
    colour = ink), panel.grid.minor = ggplot2::element_blank(), plot.margin = ggplot2::margin(16, 24,
    16, 16))
export <- function(p, name, w = 12, h = 8) {
    target <- file.path(figdir, name)
    ggplot2::ggsave(target, p, device = ragg::agg_webp, width = w, height = h, units = "in", dpi = 300,
        bg = "white")
    target
}
figs <- character()
# Heatmap: gene-specific proportional hazards diagnostic significance; missing remains NA.
phd <- data.table::copy(res)
phd[, PH_category := factor(ifelse(!is.finite(Gene_PH_P), "Not estimable", ifelse(Gene_PH_P < 0.05,
    "Potential violation (p<0.05)", "No evidence (p>=0.05)")), levels = c("Potential violation (p<0.05)",
    "No evidence (p>=0.05)", "Not estimable"))]
phd[, Feature_ID := factor(Feature_ID, levels = rev(features))]
phd[, Cancer := factor(Cancer, levels = cancers)]
p <- ggplot2::ggplot(phd, ggplot2::aes(x = Cancer, y = Feature_ID)) + ggplot2::geom_tile(ggplot2::aes(fill = PH_category),
    colour = "white", linewidth = 0.85, width = 0.94, height = 0.91) + ggplot2::scale_fill_manual(values = c(`Potential violation (p<0.05)` = "#A83736",
    `No evidence (p>=0.05)` = "#1A7879", `Not estimable` = "#BAC1C9"), drop = FALSE) + ggplot2::labs(title = "Gene-level proportional hazards diagnostics",
    subtitle = "Schoenfeld residual tests; p < 0.05 indicates potential time-varying association, not proof of failure",
    x = "TCGA cancer", y = NULL, fill = "Test outcome") + base_theme() + ggplot2::theme(legend.position = "bottom")
figs <- c(figs, export(p, "Stage08_gene_level_PH_diagnostic_heatmap.webp", 13, 9))
# Stability forest: before vs after top 1% DFBETA observation exclusion.
focus <- data.table::copy(res[(Cancer == "PAAD" & Feature_ID == "MET") | (Cancer == "KIRC" & Feature_ID ==
    "STX1A") | (Cancer == "LUAD" & Feature_ID == "DKK1") | (Cancer == "BRCA" & Feature_ID == "PGK1") |
    (Feature_ID == "SERPINE1")])
focus <- focus[is.finite(Adjusted_HR) & is.finite(Sensitivity_HR) & Sensitivity_HR > 0]
if (nrow(focus) > 0) {
    focus[, Label := paste(Cancer, Feature_ID, sep = "  |  ")]
    focus[, Label := factor(Label, levels = rev(unique(Label)))]
    lo <- data.table::rbindlist(list(focus[, .(Label, Model = "Primary adjusted", HR = Adjusted_HR)],
        focus[, .(Label, Model = "Excluding top 1% influence", HR = Sensitivity_HR)]))
    p <- ggplot2::ggplot(lo, ggplot2::aes(x = HR, y = Label, colour = Model, shape = Model)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = "#66717B") + ggplot2::geom_point(size = 3.2, position = ggplot2::position_dodge(width = 0.55)) +
        ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = c(`Primary adjusted` = "#234E70",
        `Excluding top 1% influence` = "#C05B38")) + ggplot2::labs(title = "Influential-case sensitivity of selected hazard ratios",
        subtitle = "Refitting after removal of top 1% DFBETA observations; diagnostic, not a replacement estimator",
        x = "Adjusted hazard ratio per within-cancer SD (log scale)", y = NULL, colour = NULL, shape = NULL) +
        base_theme() + ggplot2::theme(legend.position = "bottom")
    figs <- c(figs, export(p, "Stage08_selected_gene_influence_sensitivity.webp", 13, 8.5))
}
# Apparent concordance gain: labelled explicitly as in-sample, not predictive validation.
cc <- data.table::copy(res[is.finite(Apparent_Delta_C)])
if (nrow(cc)) {
    cc <- cc[, .(Median_Delta_C = median(Apparent_Delta_C), Q1 = as.numeric(quantile(Apparent_Delta_C,
        0.25)), Q3 = as.numeric(quantile(Apparent_Delta_C, 0.75)), N = .N), by = Cancer]
    cc[, Cancer := factor(Cancer, levels = rev(cancers))]
    p <- ggplot2::ggplot(cc, ggplot2::aes(x = Median_Delta_C, y = Cancer)) + ggplot2::geom_vline(xintercept = 0,
        linetype = "dashed", colour = "#7A8490") + ggplot2::geom_errorbar(ggplot2::aes(xmin = Q1, xmax = Q3),
        orientation = "y", width = 0.21, linewidth = 1, colour = "#8BA4BC") + ggplot2::geom_point(size = 3.6,
        colour = "#1C5678") + ggplot2::labs(title = "Apparent discrimination gain from adding a gene",
        subtitle = "Median and interquartile range across selected genes within cancer; IN-SAMPLE and selection-biased",
        x = "Apparent \316\224 concordance (gene + clinical minus clinical)", y = NULL) + base_theme()
    figs <- c(figs, export(p, "Stage08_apparent_concordance_gain_by_cancer.webp", 12, 7.4))
}
# Combined per-gene PH 95% CI figure with explicit p-values for PH diagnostic.
q <- data.table::copy(res[Feature_ID %in% c("SERPINE1", "MET", "STX1A", "PGK1", "DKK1") & is.finite(CI_low) &
    is.finite(CI_high) & CI_low > 0])
if (nrow(q)) {
    for (gene in unique(q$Feature_ID)) {
        t <- data.table::copy(q[Feature_ID == gene])
        t[, Cancer := factor(Cancer, levels = rev(cancers))]
        t[, Flag := ifelse(is.finite(Gene_PH_P) & Gene_PH_P < 0.05, "Potential gene PH violation",
            "No detected gene PH violation")]
        p <- ggplot2::ggplot(t, ggplot2::aes(x = Adjusted_HR, y = Cancer, colour = Flag)) + ggplot2::geom_vline(xintercept = 1,
            linetype = "dashed", colour = "#7A8490") + ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low,
            xmax = CI_high), orientation = "y", width = 0.2, linewidth = 0.85) + ggplot2::geom_point(size = 3) +
            ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = c(`Potential gene PH violation` = "#B4433D",
            `No detected gene PH violation` = "#237A74")) + ggplot2::labs(title = paste0(gene, " | adjusted HR and proportional hazards diagnostic"),
            subtitle = "95% Wald CI; gene-level Schoenfeld test; observational, outcome-selected candidates",
            x = "Adjusted HR per within-cancer SD (log scale)", y = NULL, colour = NULL) + base_theme() +
            ggplot2::theme(legend.position = "bottom")
        figs <- c(figs, export(p, paste0("Stage08_", gene, "_PH_flagged_forest_95CI.webp"), 12, 8))
    }
}
data.table::fwrite(data.table::data.table(Figure_path = figs), file.path(out, "Stage08_figure_index.csv"))

report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 08 \342\200\224 CONSOLIDATED RESULTS\n====================================================================\n")
    cat("1. ANALYSES COMPLETED\n")
    print(data.table(Cancers = uniqueN(res$Cancer), Genes_selected = length(features), Models_completed = nrow(res),
        Models_problem = nrow(err)))
    cat("\n2. SELECTED FEATURES\n")
    print(selection)
    cat("\n3. MODEL DIAGNOSTICS BY CANCER\n")
    print(summary_cancer)
    cat("\n4. PH VIOLATIONS (GENE TEST P<0.05): TOP 25\n")
    print(head(res[is.finite(Gene_PH_P) & Gene_PH_P < 0.05][order(Gene_PH_P), .(Cancer, Feature_ID, N,
        Events, Adjusted_HR, Gene_PH_P, Global_PH_P)], 25))
    cat("\n5. HIGH-INFLUENCE / SENSITIVITY EFFECTS: TOP 25 BY ABS CHANGE IN HR\n")
    print(head(res[is.finite(HR_change_pct)][order(-abs(HR_change_pct)), .(Cancer, Feature_ID, Adjusted_HR,
        Sensitivity_HR, HR_change_pct, Max_abs_DFBETA, Influential_patients)], 25))
    cat("\n6. INCREMENTAL ASSOCIATION AND APPARENT DISCRIMINATION: TOP 25\n")
    print(head(res[is.finite(LRT_P)][order(LRT_P), .(Cancer, Feature_ID, LRT_P, Apparent_C_clinical,
        Apparent_C_plus_gene, Apparent_Delta_C, Gene_PH_P)], 25))
    cat("\n7. KEY CROSS-CANCER EXAMPLES\n")
    print(res[Feature_ID %in% anchors, .(Cancer, Feature_ID, Adjusted_HR, Gene_PH_P, Global_PH_P, Sensitivity_HR,
        HR_change_pct, LRT_P, Apparent_Delta_C)][order(Feature_ID, Cancer)])
    cat("\n8. MODEL PROBLEMS\n")
    print(err)
    cat("\n9. INTERPRETATION NOTES\n")
    cat("All analyses remain exploratory because selection used TCGA outcomes. P-values and exploratory FDR are NOT external validation.\n")
    cat("Schoenfeld tests examine the PH assumption; low P indicates potential time-varying effect, not proof of model invalidity.\n")
    cat("DFBETA threshold 2/sqrt(events) is a screening heuristic, not a deletion criterion.\n")
    cat("Sensitivity excludes the top 1% most influential observations for the gene coefficient; report shifts, do not replace primary estimates.\n")
    cat("Delta-C values are APPARENT (in-sample), optimistic, and NOT validated predictive gains.\n")
    cat("LRT compares nested clinical vs clinical+gene models on the SAME patients; candidate selection makes inference exploratory.\n")
    cat("\n10. JOURNAL-ORIENTED FIGURES (WebP; actual visual quality requires inspection)\n")
    print(figs)
    cat("\n11. OUTPUT DIRECTORY\n", out, "\n====================================================================\n")
})
writeLines(report, file.path(out, "Stage08_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Proportional-hazards
# and influence resolution
# ----------------------------------------------------------------------------

options(width = 180, scipen = 6)
root <- PROJECT_ROOT
in03 <- file.path(root, "Results_Rebuild", "Stage03")
in08 <- file.path(root, "Results_Rebuild", "Stage08")
out <- file.path(root, "Results_Rebuild", "Stage08B")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("data.table", "survival", "ggplot2", "ragg")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ",
    pkg)
suppressPackageStartupMessages(library(data.table))
checkfile <- function(x) {
    if (!file.exists(x))
        stop("Missing required file: ", x)
    x
}
prior <- fread(checkfile(file.path(in08, "Stage08_all_candidate_diagnostics.csv")))
stopifnot(all(c("Cancer", "Feature_ID", "Adjusted_HR", "Sensitivity_HR", "Excluded_in_sensitivity") %in%
    names(prior)))
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
prepare_clinical <- function(cl, cc) {
    n <- nrow(cl)
    d <- data.frame(time = as.numeric(cl$OS_time), event = as.integer(cl$OS))
    if ("age_at_initial_pathologic_diagnosis" %in% names(cl)) {
        age <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
        age[!is.finite(age) | age < 18 | age > 105] <- NA_real_
        d$age10 <- age/10
    }
    if ("gender" %in% names(cl)) {
        sx <- toupper(clean_chr(cl$gender))
        sx[!sx %in% c("MALE", "FEMALE")] <- NA_character_
        if (sum(!is.na(sx)) >= 0.8 * n && length(table(sx)) > 1 && min(table(sx)) >= 10)
            d$sex <- factor(sx)
    }
    a <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        major_stage(cl$ajcc_pathologic_tumor_stage) else rep(NA_character_, n)
    b <- if ("clinical_stage" %in% names(cl))
        major_stage(cl$clinical_stage) else rep(NA_character_, n)
    st <- if (sum(!is.na(a)) >= sum(!is.na(b)))
        a else b
    if (sum(!is.na(st)) >= 0.8 * n && length(unique(na.omit(st))) >= 2)
        d$stage <- factor(st, levels = c("I", "II", "III", "IV"))
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- toupper(clean_chr(cl$histological_type))
        group <- rep(NA_character_, n)
        group[grepl("ENDOMETRIOID", h)] <- "Endometrioid"
        group[grepl("SEROUS|MIXED", h)] <- "Serous_or_mixed"
        if (sum(!is.na(group)) >= 0.8 * n && length(table(group)) > 1 && all(table(group) >= 10))
            d$histology <- factor(group)
    }
    d
}


# Focus on global PH failures, gene-level PH failure, and influential refits.  Two genes per
# global-PH cancer give a check that clinical drivers recur.
checks <- data.table(Cancer = c("BRCA", "BRCA", "KIRC", "KIRC", "OV", "OV", "UCEC", "UCEC", "LUAD", "PAAD"),
    Feature_ID = c("PGK1", "SERPINE1", "STX1A", "DONSON", "SERPINE1", "MET", "DKK1", "SERPINE1", "DKK1",
        "MET"))
# Values below are descriptive flags, NOT selection by significance for validation.
term_rows <- list()
influence_rows <- list()
loo_rows <- list()
model_rows <- list()
figpaths <- character()
problems <- character()
for (cc in unique(checks$Cancer)) {
    message("Stage08B: ", cc)
    dat <- readRDS(checkfile(file.path(in03, paste0("TCGA_", cc, "_discovery.rds"))))
    x <- dat$expression
    cl <- as.data.table(dat$clinical)
    if (!identical(colnames(x), as.character(cl$sample_id)))
        stop("Expression/clinical sample mismatch: ", cc)
    original <- prepare_clinical(cl, cc)
    include <- complete.cases(original)
    d0 <- original[include, , drop = FALSE]
    d0[] <- lapply(d0, function(z) if (is.factor(z))
        droplevels(z) else z)
    terms <- setdiff(names(d0), c("time", "event"))
    terms <- Filter(function(k) length(unique(d0[[k]])) >= 2, terms)
    d0 <- d0[, c("time", "event", terms), drop = FALSE]
    for (gene in checks[Cancer == cc]$Feature_ID) {
        if (!gene %in% rownames(x)) {
            problems <- c(problems, paste(cc, gene, "gene absent"))
            next
        }
        raw <- as.numeric(x[gene, ])
        if (anyNA(raw) || sd(raw) <= 0) {
            problems <- c(problems, paste(cc, gene, "invalid measurements"))
            next
        }
        df <- d0
        df$z <- as.numeric(scale(raw))[include]  # SAME full cancer z-scale; do not restandardise after deletion
        f <- as.formula(paste("survival::Surv(time,event) ~ z", if (length(terms))
            paste0(" + ", paste(terms, collapse = " + ")) else ""))
        fit <- tryCatch(survival::coxph(f, data = df, ties = "efron", x = TRUE, y = TRUE), error = function(e) e)
        if (inherits(fit, "error")) {
            problems <- c(problems, paste(cc, gene, conditionMessage(fit)))
            next
        }
        ph <- tryCatch(survival::cox.zph(fit, transform = "km", terms = TRUE, singledf = FALSE), error = function(e) e)
        if (inherits(ph, "error")) {
            problems <- c(problems, paste(cc, gene, "PH:", conditionMessage(ph)))
            next
        }
        tb <- as.data.table(as.data.frame(ph$table), keep.rownames = "Term")
        if (ncol(tb) < 4)
            stop("Unexpected cox.zph table schema")
        setnames(tb, names(tb)[2:4], c("Chisq", "Df", "P"))
        tb[, `:=`(Cancer = cc, Feature_ID = gene, N = nrow(df), Events = sum(df$event))]
        term_rows[[length(term_rows) + 1L]] <- tb[, .(Cancer, Feature_ID, Term, Chisq, Df, P, N, Events)]
        # Explicit comparison: unscaled DFBETA versus scaled DFBETAS.
        db <- residuals(fit, type = "dfbeta")
        ds <- residuals(fit, type = "dfbetas")
        if (is.null(dim(db)))
            db <- matrix(db, ncol = length(coef(fit)))
        if (is.null(dim(ds)))
            ds <- matrix(ds, ncol = length(coef(fit)))
        idx <- match("z", names(coef(fit)))
        stopifnot(!is.na(idx), nrow(db) == nrow(df), ncol(db) == length(coef(fit)), ncol(ds) == ncol(db))
        zd <- abs(as.numeric(db[, idx]))
        zs <- abs(as.numeric(ds[, idx]))
        # 2/sqrt(N_events) is intended as heuristic on scaled DFBETAS; original Stage08 applied it
        # to unscaled DFBETA, and hence reported no above-threshold observations.
        threshold <- 2/sqrt(sum(df$event))
        k <- max(1, ceiling(0.01 * nrow(df)))
        drop <- head(order(zd, decreasing = TRUE), k)
        refit <- survival::coxph(f, data = df[-drop, , drop = FALSE], ties = "efron")
        beta <- unname(coef(fit)["z"])
        b2 <- unname(coef(refit)["z"])
        # LOO for top 5 most influential cases; interpretable stress test.
        for (i in head(order(zd, decreasing = TRUE), 5)) {
            lf <- tryCatch(survival::coxph(f, data = df[-i, , drop = FALSE], ties = "efron"), error = function(e) e)
            if (!inherits(lf, "error"))
                loo_rows[[length(loo_rows) + 1]] <- data.table(Cancer = cc, Feature_ID = gene, Row_in_adjusted_cohort = i,
                  Sample = as.character(cl$sample_id[which(include)[i]]), Event = df$event[i], Abs_DFBETA = zd[i],
                  Abs_DFBETAS = zs[i], Leave_one_out_HR = exp(unname(coef(lf)["z"])), LOO_HR_change_pct = 100 *
                    (exp(unname(coef(lf)["z"]) - beta) - 1))
        }
        old <- prior[Cancer == cc & Feature_ID == gene]
        influence_rows[[length(influence_rows) + 1]] <- data.table(Cancer = cc, Feature_ID = gene, N = nrow(df),
            Events = sum(df$event), Adjusted_HR = exp(beta), Sensitivity_HR = exp(b2), HR_change_pct = 100 *
                (exp(b2 - beta) - 1), DFBETA_max = max(zd), DFBETAS_max = max(zs), Threshold_2_per_sqrt_events = threshold,
            DFBETA_above_threshold = sum(zd > threshold), DFBETAS_above_threshold = sum(zs > threshold),
            Dropped_top1pct = length(drop), Stage08_reported_primary_HR = if (nrow(old))
                old$Adjusted_HR[1] else NA_real_, Stage08_reported_sensitivity_HR = if (nrow(old))
                old$Sensitivity_HR[1] else NA_real_, Primary_HR_reproduction_error = if (nrow(old))
                abs(exp(beta) - old$Adjusted_HR[1]) else NA_real_, Sensitivity_HR_reproduction_error = if (nrow(old))
                abs(exp(b2) - old$Sensitivity_HR[1]) else NA_real_)
        model_rows[[length(model_rows) + 1]] <- data.table(Cancer = cc, Feature_ID = gene, Term_count = length(terms),
            Covariates = paste(terms, collapse = " + "), Gene_PH_P = tb[Term == "z", P][1], Global_PH_P = tb[Term ==
                "GLOBAL", P][1])
        # A genuine scaled-Schoenfeld residual plot with smoothed coefficient estimate and
        # uncertainty.  Plot is a single figure for priority model; labels bold/black via base
        # graphics.
        if ((cc == "LUAD" && gene == "DKK1") || (cc == "PAAD" && gene == "MET") || (cc == "KIRC" && gene ==
            "STX1A")) {
            path <- file.path(figdir, paste0("Stage08B_", cc, "_", gene, "_Schoenfeld_timevarying.webp"))
            ragg::agg_webp(path, width = 2400, height = 1650, res = 300, background = "white")
            op <- par(mar = c(5.4, 5.4, 3.6, 1.3), cex.axis = 1.1, font.axis = 2, font.lab = 2, col.axis = "black",
                col.lab = "black")
            tryCatch({
                plot(ph, var = "z", main = paste(cc, gene, "| Schoenfeld residual diagnostic"), xlab = "Follow-up time (days; KM-transformed axis)",
                  ylab = "Estimated log hazard ratio over time", col = "#174E73", lwd = 2)
                abline(h = beta, lty = 2, col = "#BA513D", lwd = 2)
            }, error = function(e) message("Plot failure: ", conditionMessage(e)))
            par(op)
            grDevices::dev.off()
            figpaths <- c(figpaths, path)
        }
    }
    rm(dat, x, cl)
    gc(verbose = FALSE)
}
terms_dt <- rbindlist(term_rows, fill = TRUE)
infl <- rbindlist(influence_rows, fill = TRUE)
loo <- rbindlist(loo_rows, fill = TRUE)
mods <- rbindlist(model_rows, fill = TRUE)
for (nm in c("terms_dt", "infl", "loo", "mods")) if (!nrow(get(nm))) stop("No output rows: ", nm)
fwrite(terms_dt, file.path(out, "Stage08B_PH_term_specific.csv"))
fwrite(infl, file.path(out, "Stage08B_DFBETA_vs_DFBETAS_and_refits.csv"))
fwrite(loo, file.path(out, "Stage08B_leave_one_out_top5.csv"))
fwrite(mods, file.path(out, "Stage08B_model_specifications.csv"))
# Modern term-level heatmap. P values exploratory and repeated across genes within cancer.
pt <- copy(terms_dt[Term != "GLOBAL"])
pt[, Flag := fifelse(is.na(P), "Not estimable", fifelse(P < 0.05, "PH concern (P < 0.05)", "No detected departure"))]
pt[, Model := paste(Cancer, Feature_ID, sep = " / ")]
pt[, Model := factor(Model, levels = rev(unique(Model)))]
p <- ggplot2::ggplot(pt, ggplot2::aes(x = Term, y = Model, fill = Flag)) + ggplot2::geom_tile(colour = "white",
    linewidth = 0.85, width = 0.94, height = 0.88) + ggplot2::scale_fill_manual(values = c(`PH concern (P < 0.05)` = "#C84B42",
    `No detected departure` = "#2A807A", `Not estimable` = "#BEC7D0"), drop = FALSE) + ggplot2::labs(title = "Which model terms violate proportional hazards?",
    subtitle = "Term-level Schoenfeld tests; exploratory diagnostics (not multiple-testing confirmation)",
    x = "Model term", y = NULL, fill = NULL) + ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 17), axis.text = ggplot2::element_text(face = "bold", colour = "black"),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), panel.grid = ggplot2::element_blank(),
    legend.position = "bottom")
fig <- file.path(figdir, "Stage08B_clinical_vs_gene_PH_heatmap.webp")
ggplot2::ggsave(fig, p, device = ragg::agg_webp, width = 13, height = 9, units = "in", dpi = 300, bg = "white")
figpaths <- c(figpaths, fig)
# Influence contrast dot plot, labels make clear the different residual scaling.
ll <- melt(infl[, .(Cancer, Feature_ID, DFBETA_above_threshold, DFBETAS_above_threshold)], id.vars = c("Cancer",
    "Feature_ID"), variable.name = "Residual_type", value.name = "Patients_flagged")
ll[, Label := paste(Cancer, Feature_ID)]
p2 <- ggplot2::ggplot(ll, ggplot2::aes(x = Patients_flagged, y = reorder(Label, Patients_flagged), shape = Residual_type,
    colour = Residual_type)) + ggplot2::geom_point(size = 3.3, position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::labs(title = "Influence diagnostics depend on residual scaling", subtitle = "Same heuristic threshold shown for unscaled DFBETA vs scaled DFBETAS; neither automatically warrants exclusion",
        x = "Patients above screening threshold", y = NULL, colour = NULL, shape = NULL) + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", colour = "black", size = 17), axis.text = ggplot2::element_text(face = "bold",
        colour = "black"), axis.title = ggplot2::element_text(face = "bold", colour = "black"), legend.position = "bottom")
fig2 <- file.path(figdir, "Stage08B_DFBETAS_influence_diagnostic.webp")
ggplot2::ggsave(fig2, p2, device = ragg::agg_webp, width = 13, height = 8, units = "in", dpi = 300, bg = "white")
figpaths <- c(figpaths, fig2)
fwrite(data.table(Figure = figpaths), file.path(out, "Stage08B_figure_index.csv"))
report <- capture.output({
    cat("==============================================================\nPROJECT 3 | STAGE 08B \342\200\224 PH AND INFLUENCE RESOLUTION\n==============================================================\n")
    cat("1. TERM-SPECIFIC PH DIAGNOSTICS (GENE + CLINICAL)\n")
    print(terms_dt)
    cat("\n2. GLOBAL VS GENE PH CHECK\n")
    print(mods)
    cat("\n3. DFBETA / DFBETAS COMPARISON AND TOP 1% REFIT\n")
    print(infl)
    cat("\n4. LARGEST SINGLE-PATIENT LEAVE-ONE-OUT CHANGES\n")
    print(loo[order(-abs(LOO_HR_change_pct))][1:min(.N, 20)])
    cat("\n5. AUTOMATED REPRODUCTION CHECK\n")
    print(infl[, .(Cancer, Feature_ID, Primary_HR_reproduction_error, Sensitivity_HR_reproduction_error)])
    cat("\n6. SCIENTIFIC NOTES\n")
    cat("Term-specific PH tests diagnose whether gene or clinical covariates contribute to nonproportional hazards.\n")
    cat("DFBETA is unscaled; DFBETAS is scaled by coefficient uncertainty. Original zero flags used a heuristic on unscaled DFBETA.\n")
    cat("Top-1-percent influence exclusion is an exploratory stress test and is NOT justified as the primary estimator.\n")
    cat("An individual gene may satisfy PH even when stage/histology/age violates the global PH test.\n")
    cat("Time-varying coefficient plots are diagnostic; changing model strategy requires a targeted biological/statistical decision.\n")
    cat("Do not interpret selected-gene P values as external validation. No source file overwritten.\n")
    cat("\n7. FIGURES\n")
    print(figpaths)
    cat("\n8. PROBLEMS\n")
    print(problems)
    cat("\nOUTPUT: ", out, "\n==============================================================\n")
})
writeLines(report, file.path(out, "Stage08B_CONSOLIDATED_REPORT.txt"))
cat(paste(report, collapse = "\n"), "\n")
