# Evaluate incremental survival discrimination under repeated cross-validation.  Separates
# outcome-conditioned candidate assessment from training-only feature selection and uses the
# corrected BRCA fold handling.  Run with: Rscript 06_internal_validation.R

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


# ----------------------------------------------------------------------------
# Candidate-conditioned internal validation
# ----------------------------------------------------------------------------

options(width = 185, scipen = 5)
root <- PROJECT_ROOT
in03 <- file.path(root, "Results_Rebuild", "Stage03")
out <- file.path(root, "Results_Rebuild", "Stage09")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("survival", "data.table", "ggplot2", "ragg")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop("Install required R package: ", pkg)
}
suppressPackageStartupMessages(library(data.table))
set.seed(20261102)
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
genes <- c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1", "DKK1")
K <- 5L
R <- 2L
required <- function(f) {
    if (!file.exists(f))
        stop("Required validated file missing: ", f)
    f
}
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
# Reproduce the Stage07 clinical design; do not select covariates by survival P value.
clinical_design <- function(cl, cc) {
    n <- nrow(cl)
    d <- data.frame(time = as.numeric(cl$OS_time), event = as.integer(cl$OS))
    if ("age_at_initial_pathologic_diagnosis" %in% names(cl)) {
        a <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
        a[!is.finite(a) | a < 18 | a > 105] <- NA_real_
        d$age10 <- a/10
    }
    if ("gender" %in% names(cl)) {
        s <- toupper(clean_chr(cl$gender))
        s[!s %in% c("MALE", "FEMALE")] <- NA_character_
        if (sum(!is.na(s)) >= 0.8 * n && length(table(s)) > 1 && min(table(s)) >= 10)
            d$sex <- factor(s)
    }
    a <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        major_stage(cl$ajcc_pathologic_tumor_stage) else rep(NA_character_, n)
    b <- if ("clinical_stage" %in% names(cl))
        major_stage(cl$clinical_stage) else rep(NA_character_, n)
    s <- if (sum(!is.na(a)) >= sum(!is.na(b)))
        a else b
    if (sum(!is.na(s)) >= 0.8 * n && length(unique(na.omit(s))) >= 2)
        d$stage <- factor(s, levels = c("I", "II", "III", "IV"))
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- toupper(clean_chr(cl$histological_type))
        g <- rep(NA_character_, n)
        g[grepl("ENDOMETRIOID", h)] <- "Endometrioid"
        g[grepl("SEROUS|MIXED", h)] <- "Serous_or_mixed"
        if (sum(!is.na(g)) >= 0.8 * n && length(table(g)) > 1 && all(table(g) >= 10))
            d$histology <- factor(g)
    }
    d
}
# In a test fold with a previously unseen factor level, score only if estimable; record failure
# transparently, never silently turn failures into successes.
safe_cox <- function(form, df) {
    wrn <- character()
    fit <- withCallingHandlers(tryCatch(survival::coxph(form, data = df, ties = "efron", control = survival::coxph.control(iter.max = 50),
        model = TRUE), error = function(e) e), warning = function(w) {
        wrn <<- c(wrn, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    if (inherits(fit, "error"))
        return(list(fit = NULL, problem = paste("fit error:", conditionMessage(fit))))
    if (length(wrn) || any(!is.finite(stats::coef(fit))))
        return(list(fit = NULL, problem = paste("fit warning/unstable:", paste(unique(wrn), collapse = "; "))))
    list(fit = fit, problem = "")
}
# Holdout Harrell C; reverse=TRUE because larger Cox linear predictor means greater hazard.
score_c <- function(fit, test) {
    pred <- tryCatch(as.numeric(stats::predict(fit, newdata = test, type = "lp", reference = "zero")),
        error = function(e) e)
    if (inherits(pred, "error") || length(pred) != nrow(test) || any(!is.finite(pred)))
        return(NA_real_)
    if (sum(test$event) < 2L || length(unique(pred)) < 2L)
        return(NA_real_)
    tryCatch(as.numeric(survival::concordance(survival::Surv(test$time, test$event) ~ pred, reverse = TRUE)$concordance),
        error = function(e) NA_real_)
}
# Stratify folds by event indicator within cancer. One shared fold map for all genes.
cv_rows <- list()
issue_rows <- list()
cohort_rows <- list()
phi_rows <- list()
idx <- 0L
ii <- 0L
hi <- 0L
for (cc in cancers) {
    message("Stage09: validating ", cc)
    o <- readRDS(required(file.path(in03, paste0("TCGA_", cc, "_discovery.rds"))))
    x <- o$expression
    cl <- as.data.table(o$clinical)
    stopifnot(identical(colnames(x), as.character(cl$sample_id)))
    d <- clinical_design(cl, cc)
    complete <- complete.cases(d) & is.finite(d$time) & d$time > 0 & d$event %in% c(0, 1)
    df <- d[complete, , drop = FALSE]
    df[] <- lapply(df, function(v) if (is.factor(v))
        droplevels(v) else v)
    covs <- setdiff(names(df), c("time", "event"))
    covs <- Filter(function(v) length(unique(df[[v]])) > 1, covs)
    df <- df[, c("time", "event", covs), drop = FALSE]
    cohort_rows[[cc]] <- data.table(Cancer = cc, N = nrow(df), Events = sum(df$event), Clinical_covariates = paste(covs,
        collapse = " + "))
    f0 <- as.formula(paste("survival::Surv(time,event) ~", if (length(covs))
        paste(covs, collapse = " + ") else "1"))
    f1 <- as.formula(paste("survival::Surv(time,event) ~ z", if (length(covs))
        paste0(" + ", paste(covs, collapse = " + ")) else ""))
    # Assign event-stratified fold ids; repeated partitions are not statistically independent.
    folds <- matrix(NA_integer_, nrow(df), R)
    for (r in seq_len(R)) for (e in c(0, 1)) {
        ind <- which(df$event == e)
        folds[ind, r] <- sample(rep(seq_len(K), length.out = length(ind)))
    }
    # Training and test cases must have the same clinical eligibility, for paired comparison.
    for (g in genes) {
        if (!g %in% rownames(x)) {
            ii <- ii + 1L
            issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = NA_integer_, Fold = NA_integer_,
                Reason = "Gene absent")
            next
        }
        values <- as.numeric(x[g, ])[complete]
        if (any(!is.finite(values))) {
            ii <- ii + 1L
            issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = NA_integer_, Fold = NA_integer_,
                Reason = "Nonfinite gene expression")
            next
        }
        for (r in seq_len(R)) for (k in seq_len(K)) {
            test <- which(folds[, r] == k)
            train <- which(folds[, r] != k)
            tr <- df[train, , drop = FALSE]
            te <- df[test, , drop = FALSE]
            # Preprocessing moments calculated ONLY on training patients, then applied to held-out
            # patients.
            mu <- mean(values[train])
            s <- stats::sd(values[train])
            if (!is.finite(s) || s < 1e-09) {
                ii <- ii + 1L
                issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Reason = "Near-constant gene in training")
                next
            }
            tr$z <- (values[train] - mu)/s
            te$z <- (values[test] - mu)/s
            # Prevent a test-only categorical level causing an unknown dummy coefficient.
            issue <- NULL
            for (v in covs) if (is.factor(tr[[v]])) {
                tr[[v]] <- droplevels(tr[[v]])
                te[[v]] <- factor(as.character(te[[v]]), levels = levels(tr[[v]]))
                if (anyNA(te[[v]]))
                  issue <- paste0("Unknown test factor level in ", v)
                if (nlevels(tr[[v]]) < 2L)
                  issue <- paste0("Collapsed training factor ", v)
            }
            if (!is.null(issue)) {
                ii <- ii + 1L
                issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Reason = issue)
                next
            }
            a <- safe_cox(f0, tr)
            b <- safe_cox(f1, tr)
            if (is.null(a$fit) || is.null(b$fit)) {
                ii <- ii + 1L
                issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Reason = paste(a$problem,
                  b$problem))
                next
            }
            c0 <- score_c(a$fit, te)
            c1 <- score_c(b$fit, te)
            if (!is.finite(c0) || !is.finite(c1)) {
                ii <- ii + 1L
                issue_rows[[ii]] <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Reason = "Holdout concordance not estimable")
                next
            }
            idx <- idx + 1L
            cv_rows[[idx]] <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, N_train = length(train),
                N_test = length(test), Deaths_test = sum(te$event), C_clinical = c0, C_plus_gene = c1,
                Delta_C = c1 - c0)
        }
        # Stage-stratified PH sensitivity is a model sensitivity, not held-out validation.
        if ("stage" %in% covs && g %in% c("SERPINE1", "STX1A", "MET", "DKK1", "PGK1")) {
            full <- df
            full$z <- (values - mean(values))/stats::sd(values)
            of <- safe_cox(f1, full)
            scovs <- setdiff(covs, "stage")
            ff <- as.formula(paste("survival::Surv(time,event) ~ z", if (length(scovs))
                paste0(" + ", paste(scovs, collapse = " + ")) else "", "+ strata(stage)"))
            sf <- safe_cox(ff, full)
            hi <- hi + 1L
            phi_rows[[hi]] <- data.table(Cancer = cc, Gene = g, Ordinary_HR = if (is.null(of$fit))
                NA_real_ else exp(unname(stats::coef(of$fit)["z"])), Stratified_HR = if (is.null(sf$fit))
                NA_real_ else exp(unname(stats::coef(sf$fit)["z"])))
        }
    }
    rm(o, x, cl)
    invisible(gc(verbose = FALSE))
}
cv <- if (length(cv_rows)) rbindlist(cv_rows) else stop("All cross-validation folds failed; inspect cohort/covariates.")
issues <- if (length(issue_rows)) rbindlist(issue_rows, fill = TRUE) else data.table(Cancer = character(),
    Gene = character(), Repeat = integer(), Fold = integer(), Reason = character())
cohorts <- rbindlist(cohort_rows)
ph <- if (length(phi_rows)) rbindlist(phi_rows) else data.table()
# Completeness counts are explicit; a missing fold does not disappear from the denominator.
summary <- cv[, .(Successful_folds = .N, Mean_C_clinical = mean(C_clinical), Mean_C_plus_gene = mean(C_plus_gene),
    Mean_Delta_C = mean(Delta_C), Median_Delta_C = median(Delta_C), Positive_folds = sum(Delta_C > 0)),
    by = .(Cancer, Gene)]
summary[, `:=`(Expected_folds = K * R, Failed_folds = K * R - Successful_folds)]
setorder(summary, -Mean_Delta_C)
if (nrow(ph)) ph[, HR_change_pct := 100 * (Stratified_HR/Ordinary_HR - 1)]
fwrite(cv, file.path(out, "Stage09_fold_level_validation.csv"))
fwrite(summary, file.path(out, "Stage09_conditional_validation_summary.csv"))
fwrite(issues, file.path(out, "Stage09_failures_and_warnings.csv"))
fwrite(cohorts, file.path(out, "Stage09_eligible_clinical_cohorts.csv"))
if (nrow(ph)) fwrite(ph, file.path(out, "Stage09_stage_stratification_sensitivity.csv"))
# Journal-style chart theme and full-raster WebP. Every chart is separate, not a multipanel.
theme_science <- function() ggplot2::theme_minimal(base_size = 16) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 19), plot.subtitle = ggplot2::element_text(colour = "#414A52", size = 12),
    axis.title = ggplot2::element_text(face = "bold", colour = "black"), axis.text = ggplot2::element_text(colour = "black"),
    panel.grid.minor = ggplot2::element_blank(), legend.title = ggplot2::element_text(face = "bold"))
figure_paths <- character()
save_plot <- function(p, name, width = 13, height = 8) {
    path <- file.path(figdir, name)
    ggplot2::ggsave(path, plot = p, device = ragg::agg_webp, width = width, height = height, units = "in",
        res = 300, background = "white")
    figure_paths <<- c(figure_paths, path)
}
# Report mean delta-C for every gene-cancer pair; missing results remain unplotted.
p <- ggplot2::ggplot(summary, ggplot2::aes(x = Mean_Delta_C, y = reorder(paste(Cancer, Gene, sep = "  |  "),
    Mean_Delta_C), colour = Mean_Delta_C > 0)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
    colour = "#60676E") + ggplot2::geom_segment(ggplot2::aes(x = 0, xend = Mean_Delta_C, yend = reorder(paste(Cancer,
    Gene, sep = "  |  "), Mean_Delta_C)), alpha = 0.45, linewidth = 0.55) + ggplot2::geom_point(size = 2.9) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#087F72", `FALSE` = "#B44747"), guide = "none") +
    ggplot2::labs(title = "Held-out discrimination: gene versus clinical model", subtitle = "Mean paired difference in Harrell C across 5 folds \303\227 2 repetitions; discovery-conditioned",
        x = "Mean held-out \316\224C (gene + clinical minus clinical)", y = NULL) + theme_science()
save_plot(p, "Stage09_conditional_heldout_deltaC_all_genes.webp", width = 16, height = 19)
# Focused gene-by-cancer heatmap: numerical \316\224C, with each cell annotated and failures
# identified separately.
p <- ggplot2::ggplot(summary, ggplot2::aes(x = Cancer, y = Gene, fill = Mean_Delta_C)) + ggplot2::geom_tile(colour = "white",
    linewidth = 1.2) + ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", Mean_Delta_C)), size = 4.2,
    colour = "black", fontface = "bold") + ggplot2::scale_fill_gradient2(low = "#B85A61", mid = "#FAFAFA",
    high = "#46A99D", midpoint = 0, name = "Mean \316\224C") + ggplot2::labs(title = "Out-of-fold prognostic discrimination across ten cancers",
    subtitle = "Positive values favour gene-augmented clinical models; cells are not external validation",
    x = "TCGA cancer cohort", y = "Discovery-selected gene") + theme_science() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
    hjust = 1))
save_plot(p, "Stage09_conditional_deltaC_gene_by_cancer_heatmap.webp", width = 15, height = 9)
# Visualize fold instability for the three primary examples, rather than hiding spread behind
# averages.
focal <- cv[(Cancer == "PAAD" & Gene == "MET") | (Cancer == "KIRC" & Gene == "STX1A") | (Cancer == "LUAD" &
    Gene == "DKK1")]
if (nrow(focal)) {
    focal[, Association := paste0(Cancer, " \302\267 ", Gene)]
    p <- ggplot2::ggplot(focal, ggplot2::aes(x = Association, y = Delta_C, colour = Association)) + ggplot2::geom_hline(yintercept = 0,
        linetype = "dashed", colour = "#63666A") + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12,
        height = 0), size = 3, alpha = 0.8) + ggplot2::stat_summary(fun = mean, geom = "point", shape = 18,
        size = 6, colour = "black") + ggplot2::labs(title = "Internal validation variability for focal associations",
        subtitle = "Dots: held-out folds; black diamonds: mean. Repeated folds are correlated, not independent replicates.",
        x = NULL, y = "Held-out \316\224C") + theme_science() + ggplot2::theme(legend.position = "none")
    save_plot(p, "Stage09_focal_genes_fold_variability.webp", width = 12, height = 8)
}
if (nrow(ph)) {
    q <- melt(ph, id.vars = c("Cancer", "Gene"), measure.vars = c("Ordinary_HR", "Stratified_HR"), variable.name = "Model",
        value.name = "HR")
    q <- q[is.finite(HR) & HR > 0]
    q[, Association := paste(Cancer, Gene, sep = " | ")]
    p <- ggplot2::ggplot(q, ggplot2::aes(x = HR, y = reorder(Association, HR), colour = Model)) + ggplot2::geom_vline(xintercept = 1,
        linetype = "dashed", colour = "#63666A") + ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.45)) +
        ggplot2::scale_x_log10() + ggplot2::scale_colour_manual(values = c(Ordinary_HR = "#9D4C55", Stratified_HR = "#117A72"),
        labels = c("Clinical stage as covariate", "Stratified baseline hazard by stage")) + ggplot2::labs(title = "Sensitivity to clinical-stage proportional hazards",
        subtitle = "Exploratory full-cohort gene effects; estimates, not independently validated", x = "Adjusted hazard ratio (log scale)",
        y = NULL, colour = "Model") + theme_science()
    save_plot(p, "Stage09_clinical_stage_stratification_sensitivity.webp", width = 15, height = 15)
}
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 09 \342\200\224 INTERNAL VALIDATION CHECKPOINT\n====================================================================\n")
    cat("1. DESIGN: cancer-specific 5-fold \303\227 2 repeats; 9 discovered genes; 10 cancers.\n")
    cat("   Expression centre/SD estimated within TRAINING folds only; same split for clinical and augmented model.\n")
    cat("   Gene identities were selected using the FULL TCGA outcomes in Stage05\342\200\22307 (selection leakage remains).\n")
    print(data.table(Cancers = length(cancers), Genes = length(genes), Expected_folds_per_pair = K *
        R, Expected_comparisons = length(cancers) * length(genes) * K * R, Successful = nrow(cv), Issues = nrow(issues)))
    cat("\n2. CLINICAL COMPLETE-CASE COHORTS\n")
    print(cohorts)
    cat("\n3. HELD-OUT PERFORMANCE: FULL SUMMARY (sorted by mean \316\224C)\n")
    print(summary)
    cat("\n4. FOCAL ASSOCIATIONS\n")
    print(summary[(Cancer == "PAAD" & Gene == "MET") | (Cancer == "KIRC" & Gene == "STX1A") | (Cancer ==
        "LUAD" & Gene == "DKK1")])
    cat("\n5. FAILURES / WARNINGS: by cancer and gene\n")
    print(if (nrow(issues))
        issues[, .(N = .N), by = .(Cancer, Gene, Reason)][order(-N)] else "None")
    cat("\n6. CLINICAL STAGE STRATIFICATION SENSITIVITY\n")
    print(if (nrow(ph))
        ph[order(-abs(HR_change_pct))] else "No estimable comparisons")
    cat("\n7. FIGURES\n")
    print(figure_paths)
    cat("\n8. SCIENTIFIC LIMITATIONS\n")
    cat("OUTCOME-INFORMED SELECTION LEAKAGE: genes selected on entire TCGA before cross-validation. This is NOT unbiased prospective assessment of gene discovery.\n")
    cat("Fold-level preprocessing is training-only; clinical covariate design follows earlier analysis rather than being outcome-optimised.\n")
    cat("Harrell C can be affected by censoring; no calibration or fixed-time IPCW C estimated. Repeats are correlated.\n")
    cat("No external validation. No proof of clinical utility. PH violations documented in Stage08B; time-varying gene effects require separate modelling.\n")
    cat("Models with estimation errors/warnings are recorded and not treated as successful; missing folds are not silently ignored.\n")
    cat("Do not interpret descriptive spread of repeated folds as a confidence interval.\n")
    cat("\nOUTPUT: ", out, "\n====================================================================\n",
        sep = "")
})
writeLines(report, file.path(out, "Stage09_scientific_checkpoint.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Training-only nested
# validation ----------------------------------------------------------------------------

options(width = 175, scipen = 5)
root <- PROJECT_ROOT
input <- file.path(root, "Results_Rebuild", "Stage03")
out <- file.path(root, "Results_Rebuild", "Stage09B")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (p in c("survival", "data.table", "ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ",
    p)
suppressPackageStartupMessages(library(data.table))
set.seed(20261102)
cancers <- c("BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "PAAD", "UCEC")
fixed_genes <- c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1", "DKK1")
K <- 5L
repeats <- 2L
N_SCREEN <- 150L
# Deliberately parsimonious PREDECLARED clinical encoding to prevent separation: age per 10 years;
# sex binary; advanced stage III/IV vs I/II; UCEC histology binary.  This differs from Stage07's
# four-level stage model and is a predictive sensitivity model.
clean <- function(x) {
    z <- toupper(trimws(as.character(x)))
    z[z %in% c("", "NA", "N/A", "NULL", "UNKNOWN", "NOT AVAILABLE", "[NOT AVAILABLE]", "[NOT EVALUATED]")] <- NA_character_
    z
}
stage_numeric <- function(x) {
    z <- gsub("^STAGE[[:space:]]*", "", clean(x))
    v <- rep(NA_real_, length(z))
    v[grepl("^IV", z)] <- 4
    v[grepl("^III", z)] <- 3
    v[is.na(v) & grepl("^II", z)] <- 2
    v[is.na(v) & grepl("^I($|A|B|C|[[:space:]])", z)] <- 1
    v
}
make_clinical <- function(cl, cc) {
    d <- data.frame(time = as.numeric(as.character(cl$OS_time)), event = as.integer(as.character(cl$OS)))
    a <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
    a[!is.finite(a) | a < 18 | a > 105] <- NA_real_
    d$age10 <- a/10
    s <- clean(cl$gender)
    s[!s %in% c("FEMALE", "MALE")] <- NA_character_
    if (cc != "OV" && cc != "UCEC" && sum(!is.na(s)) >= 0.8 * nrow(d) && length(unique(na.omit(s))) ==
        2)
        d$male <- as.numeric(s == "MALE")
    p <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        stage_numeric(cl$ajcc_pathologic_tumor_stage) else rep(NA_real_, nrow(d))
    q <- if ("clinical_stage" %in% names(cl))
        stage_numeric(cl$clinical_stage) else rep(NA_real_, nrow(d))
    st <- if (sum(is.finite(p)) >= sum(is.finite(q)))
        p else q
    if (sum(is.finite(st)) >= 0.8 * nrow(d) && length(unique(na.omit(st))) >= 2)
        d$advanced <- as.numeric(st >= 3)
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- clean(cl$histological_type)
        v <- rep(NA_real_, length(h))
        v[grepl("ENDOMETRIOID", h)] <- 0
        v[grepl("SEROUS|MIXED", h)] <- 1
        if (sum(is.finite(v)) >= 0.8 * nrow(d) && length(unique(na.omit(v))) == 2)
            d$serous_mixed <- v
    }
    d
}
get_fit <- function(form, data) {
    w <- character()
    fit <- withCallingHandlers(tryCatch(survival::coxph(form, data = data, ties = "efron", control = survival::coxph.control(iter.max = 60),
        model = TRUE), error = function(e) e), warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
    })
    if (inherits(fit, "error"))
        return(list(fit = NULL, warning = paste(unique(w), collapse = " | "), reason = conditionMessage(fit)))
    b <- coef(fit)
    vv <- tryCatch(diag(vcov(fit)), error = function(e) NA_real_)
    ok <- length(b) > 0 && all(is.finite(b)) && all(is.finite(vv)) && all(vv > 0) && all(abs(b) < 12)
    # Cox separation warnings are NOT accepted as valid fits, even if coefficients finite.
    sep <- any(grepl("infinite|did not converge|singular|ran out of iterations", w, ignore.case = TRUE))
    list(fit = if (ok && !sep) fit else NULL, warning = paste(unique(w), collapse = " | "), reason = if (!ok) "Invalid/singular coefficient or variance" else if (sep) "Suspected separation/nonconvergence" else "")
}
concord <- function(fit, test) {
    if (is.null(fit))
        return(NA_real_)
    lp <- tryCatch(as.numeric(predict(fit, newdata = test, type = "lp", reference = "zero")), error = function(e) rep(NA_real_,
        nrow(test)))
    if (length(lp) != nrow(test) || any(!is.finite(lp)) || sum(test$event) < 2 || length(unique(lp)) <
        2)
        return(NA_real_)
    tryCatch(as.numeric(survival::concordance(survival::Surv(test$time, test$event) ~ lp, reverse = TRUE)$concordance),
        error = function(e) NA_real_)
}
rows <- list()
issues <- list()
nested <- list()
cohorts <- list()
ir <- 0L
ii <- 0L
ni <- 0L
for (cc in cancers) {
    message("Stage09B: ", cc, " | correcting CV and nested selection")
    path <- file.path(input, paste0("TCGA_", cc, "_discovery.rds"))
    if (!file.exists(path))
        stop("Missing: ", path)
    o <- readRDS(path)
    x <- o$expression
    cl <- as.data.table(o$clinical)
    stopifnot(is.matrix(x), identical(colnames(x), as.character(cl$sample_id)))
    d0 <- make_clinical(cl, cc)
    keep <- complete.cases(d0) & is.finite(d0$time) & d0$time > 0 & d0$event %in% c(0, 1)
    d <- d0[keep, , drop = FALSE]
    covars <- setdiff(names(d), c("time", "event"))
    # Drop covariates if zero variation within THIS cohort (not based on outcomes).
    covars <- covars[vapply(d[covars], function(z) length(unique(z)) > 1, logical(1))]
    d <- d[c("time", "event", covars)]
    cohorts[[cc]] <- data.table(Cancer = cc, N = nrow(d), Events = sum(d$event), Covariates = paste(covars,
        collapse = " + "))
    x <- x[, keep, drop = FALSE]
    f0 <- as.formula(paste("survival::Surv(time,event)~", paste(covars, collapse = "+")))
    f1 <- as.formula(paste("survival::Surv(time,event)~", paste(c("z", covars), collapse = "+")))
    foldid <- matrix(NA_integer_, nrow(d), repeats)
    for (r in seq_len(repeats)) for (ev in 0:1) {
        ix <- which(d$event == ev)
        foldid[ix, r] <- sample(rep(seq_len(K), length.out = length(ix)))
    }
    for (r in seq_len(repeats)) for (k in seq_len(K)) {
        train <- which(foldid[, r] != k)
        test <- which(foldid[, r] == k)
        tr <- d[train, , drop = FALSE]
        te <- d[test, , drop = FALSE]
        # Sparse training binary covariates can separate: predeclared ridge is not silently
        # substituted.
        bfit <- get_fit(f0, tr)
        bc <- concord(bfit$fit, te)
        # Prespecified-gene arm: explicitly still affected by prior outcome-based selection.
        for (g in fixed_genes) {
            result <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Analysis = "Discovered_gene_conditional",
                C_clinical = bc, C_gene = NA_real_, Delta_C = NA_real_, Status = "")
            if (!g %in% rownames(x))
                result$Status <- "gene_absent" else {
                v <- as.numeric(x[g, ])
                m <- mean(v[train])
                s <- sd(v[train])
                if (any(!is.finite(v)) || !is.finite(s) || s < 1e-09)
                  result$Status <- "invalid_expression" else {
                  train_df <- tr
                  test_df <- te
                  train_df$z <- (v[train] - m)/s
                  test_df$z <- (v[test] - m)/s
                  gf <- get_fit(f1, train_df)
                  gc <- concord(gf$fit, test_df)
                  result$C_gene <- gc
                  result$Delta_C <- gc - bc
                  result$Status <- if (is.finite(result$Delta_C))
                    "ok" else paste0("fit_or_score_failed: ", gf$reason, "; ", bfit$reason)
                  if (nzchar(gf$warning) || nzchar(bfit$warning)) {
                    ii <- ii + 1L
                    issues[[ii]] <- data.table(Cancer = cc, Repeat = r, Fold = k, Analysis = "fixed",
                      Gene = g, Warnings = paste(gf$warning, bfit$warning, sep = " | "))
                  }
                }
            }
            ir <- ir + 1L
            rows[[ir]] <- result
        }
        # Nested arm: within the TRAIN fold only, choose 150 most variable genes, then select the
        # ONE strongest training-adjusted gene (lowest valid Wald P).  Both screening steps must be
        # repeated independently for EVERY test fold.
        nn <- data.table(Cancer = cc, Gene = NA_character_, Repeat = r, Fold = k, Analysis = "Nested_train_selected",
            C_clinical = bc, C_gene = NA_real_, Delta_C = NA_real_, Status = "", Candidates_screened = 0L)
        if (!is.finite(bc))
            nn$Status <- paste0("clinical_failure: ", bfit$reason) else {
            m <- rowMeans(x[, train, drop = FALSE])
            s2 <- rowMeans(x[, train, drop = FALSE]^2) - m^2
            s2[pmax(s2, 0) < 1e-12 | !is.finite(s2) | !is.finite(m)] <- NA_real_
            eligible <- which(is.finite(s2))
            eligible <- eligible[!grepl("^[0-9]+$", rownames(x)[eligible])]
            eligible <- eligible[!duplicated(rownames(x)[eligible])]
            top <- head(eligible[order(s2[eligible], decreasing = TRUE)], N_SCREEN)
            bestp <- Inf
            bestg <- NA_character_
            bestfit <- NULL
            ntested <- 0L
            bestmu <- NA_real_
            bestsd <- NA_real_
            for (j in top) {
                v <- as.numeric(x[j, ])
                ss <- sqrt(s2[j])
                tt <- tr
                tt$z <- (v[train] - m[j])/ss
                fit <- get_fit(f1, tt)
                if (is.null(fit$fit))
                  next
                sm <- tryCatch(summary(fit$fit)$coefficients, error = function(e) NULL)
                if (is.null(sm) || !("z" %in% rownames(sm)))
                  next
                p <- unname(sm["z", "Pr(>|z|)"])
                ntested <- ntested + 1L
                if (is.finite(p) && p < bestp) {
                  bestp <- p
                  bestg <- rownames(x)[j]
                  bestfit <- fit$fit
                  bestmu <- m[j]
                  bestsd <- ss
                }
            }
            nn$Candidates_screened <- ntested
            if (is.null(bestfit))
                nn$Status <- "no_valid_training_gene" else {
                te$z <- (as.numeric(x[bestg, test]) - bestmu)/bestsd
                cc1 <- concord(bestfit, te)
                nn$Gene <- bestg
                nn$C_gene <- cc1
                nn$Delta_C <- cc1 - bc
                nn$Status <- if (is.finite(cc1))
                  "ok" else "nested_test_score_failed"
            }
        }
        ni <- ni + 1L
        nested[[ni]] <- nn
    }
    rm(o, x, cl)
    invisible(gc(verbose = FALSE))
}
fixed <- rbindlist(rows, fill = TRUE)
nested <- rbindlist(nested, fill = TRUE)
all <- rbindlist(list(fixed, nested), fill = TRUE)
cohorts <- rbindlist(cohorts)
summary <- all[, .(Expected = .N, Successful = sum(Status == "ok"), Mean_C_clinical = if (any(Status ==
    "ok")) mean(C_clinical[Status == "ok"]) else NA_real_, Mean_C_plus_gene = if (any(Status == "ok")) mean(C_gene[Status ==
    "ok"]) else NA_real_, Mean_Delta_C = if (any(Status == "ok")) mean(Delta_C[Status == "ok"]) else NA_real_,
    Positive_folds = sum(Delta_C[Status == "ok"] > 0)), by = .(Cancer, Gene, Analysis)]
setorder(summary, Analysis, -Mean_Delta_C)
issues <- if (length(issues)) rbindlist(issues) else data.table()
fwrite(all, file.path(out, "Stage09B_all_validation_attempts.csv"))
fwrite(summary, file.path(out, "Stage09B_validation_summary.csv"))
fwrite(cohorts, file.path(out, "Stage09B_clinical_cohorts.csv"))
if (nrow(issues)) fwrite(issues, file.path(out, "Stage09B_fit_warnings.csv"))
# Results with failed folds are reported WITH their original expected denominator.
fail <- all[Status != "ok", .(Failures = .N), by = .(Cancer, Analysis, Status)]
fwrite(fail, file.path(out, "Stage09B_failure_details.csv"))
theme_pub <- function() ggplot2::theme_minimal(base_size = 16) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    size = 19, colour = "black"), axis.title = ggplot2::element_text(face = "bold", colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"), panel.grid.minor = ggplot2::element_blank())
figures <- character()
savefig <- function(p, n, w = 13, h = 8) {
    f <- file.path(figdir, n)
    ggplot2::ggsave(f, plot = p, device = ragg::agg_webp, width = w, height = h, units = "in", res = 300,
        background = "white")
    figures <<- c(figures, f)
}
ss <- summary[Analysis == "Discovered_gene_conditional" & is.finite(Mean_Delta_C)]
if (nrow(ss)) {
    p <- ggplot2::ggplot(ss, ggplot2::aes(x = Cancer, y = Gene, fill = Mean_Delta_C)) + ggplot2::geom_tile(colour = "white",
        linewidth = 1) + ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f", Mean_Delta_C)), colour = "black",
        fontface = "bold", size = 4) + ggplot2::scale_fill_gradient2(low = "#BF535B", mid = "white",
        high = "#168C82", midpoint = 0, name = "Held-out \316\224C") + ggplot2::labs(title = "Discovery-conditioned prognostic discrimination",
        subtitle = "Preselected genes; held-out folds do NOT remove original gene-selection bias", x = "Cancer",
        y = "Candidate gene") + theme_pub()
    savefig(p, "Stage09B_discovered_gene_conditional_deltaC_heatmap.webp", 15, 9)
}
ss <- summary[Analysis == "Nested_train_selected" & is.finite(Mean_Delta_C)]
if (nrow(ss)) {
    p <- ggplot2::ggplot(ss, ggplot2::aes(x = Mean_Delta_C, y = reorder(Cancer, Mean_Delta_C), colour = Mean_Delta_C >
        0)) + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") + ggplot2::geom_segment(ggplot2::aes(x = 0,
        xend = Mean_Delta_C, yend = Cancer), linewidth = 1) + ggplot2::geom_point(size = 5) + ggplot2::scale_colour_manual(values = c(`TRUE` = "#087F72",
        `FALSE` = "#C35259"), guide = "none") + ggplot2::labs(title = "Training-fold-selected gene: internal predictive value",
        subtitle = "150 training-high-variance genes; one selected by training-adjusted Cox P per fold",
        x = "Mean held-out \316\224C", y = NULL) + theme_pub()
    savefig(p, "Stage09B_nested_gene_discovery_deltaC.webp", 13, 8)
}
focal <- all[Analysis == "Discovered_gene_conditional" & ((Cancer == "PAAD" & Gene == "MET") | (Cancer ==
    "KIRC" & Gene == "STX1A") | (Cancer == "LUAD" & Gene == "DKK1")) & Status == "ok"]
if (nrow(focal)) {
    focal[, Label := paste(Cancer, Gene, sep = " \302\267 ")]
    p <- ggplot2::ggplot(focal, ggplot2::aes(x = Label, y = Delta_C, colour = Label)) + ggplot2::geom_hline(yintercept = 0,
        linetype = "dashed") + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.12,
        height = 0), size = 3) + ggplot2::stat_summary(fun = mean, geom = "point", shape = 18, size = 6,
        colour = "black") + ggplot2::labs(title = "Focal genes: cross-validation fold variability", subtitle = "Points are correlated folds, not independent 95% confidence intervals",
        x = NULL, y = "Held-out \316\224C") + theme_pub() + ggplot2::theme(legend.position = "none")
    savefig(p, "Stage09B_focal_fold_variability.webp", 12, 8)
}
report <- capture.output({
    cat("====================================================================\nPROJECT 3 | STAGE 09B \342\200\224 CORRECTED AND NESTED INTERNAL VALIDATION\n====================================================================\n")
    cat("1. DESIGN\n")
    cat("5-fold x 2 repeats, ten cancers, fixed preselected 9-gene arm AND separate training-only nested-selection arm.\n")
    cat("Nested candidate pool: top ", N_SCREEN, " training-variable complete genes per fold; selection uses training adjusted Cox only.\n",
        sep = "")
    cat("Clinical simplification: age per decade, sex binary, advanced stage III/IV vs I/II, UCEC histology binary (when available).\n")
    cat("\n2. COMPLETENESS BY ANALYSIS\n")
    print(all[, .(Expected = .N, Success = sum(Status == "ok"), Failed = sum(Status != "ok")), by = Analysis])
    cat("\n3. CLINICAL COHORTS\n")
    print(cohorts)
    cat("\n4. DISCOVERY-CONDITIONED RESULTS (TOP 20)\n")
    print(head(summary[Analysis == "Discovered_gene_conditional"], 20))
    cat("\n5. NESTED TRAINING-SELECTED GENE PERFORMANCE (ALL CANCERS)\n")
    print(summary[Analysis == "Nested_train_selected"])
    cat("\n6. NESTED SELECTED GENES BY FOLD\n")
    print(nested[, .(Cancer, Repeat, Fold, Gene, Candidates_screened, Delta_C, Status)])
    cat("\n7. FAILED FOLDS (NEVER HIDDEN)\n")
    print(fail)
    cat("\n8. FIGURES\n")
    print(figures)
    cat("\n9. INTERPRETATION\n")
    cat("Nested arm assesses a PROSPECTIVE training-only selection ALGORITHM, not the prior nominated 9 genes.\n")
    cat("Discovery-conditioned candidate arm remains outcome-selection biased. Both arms are INTERNAL CV, not external validation.\n")
    cat("Only valid folds are averaged; missingness may bias mean. Interpret completion counts. Repeated folds are correlated.\n")
    cat("Stage used a deliberately simpler clinical model to reduce separation; not identical to Stage07 adjustment.\n")
    cat("C-index is discrimination, not calibration or clinical utility. No patient-level predicted survival calibration performed.\n")
    cat("\nOUTPUT: ", out, "\n====================================================================\n",
        sep = "")
})
writeLines(report, file.path(out, "Stage09B_scientific_checkpoint.txt"))
cat(paste(report, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------- Corrected BRCA folds
# and algorithm summary
# ----------------------------------------------------------------------------

options(width = 160, scipen = 5)
root <- PROJECT_ROOT
input <- file.path(root, "Results_Rebuild", "Stage03")
previous <- file.path(root, "Results_Rebuild", "Stage09B")
out <- file.path(root, "Results_Rebuild", "Stage09C")
figdir <- file.path(out, "Figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
for (p in c("survival", "data.table", "ggplot2", "ragg")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ",
    p)
suppressPackageStartupMessages(library(data.table))
set.seed(20261102)
fixed_genes <- c("SERPINE1", "FLNC", "STX1A", "PSMD2", "EXPH5", "PLK1", "MET", "PGK1", "DKK1")
K <- 5L
repeats <- 2L
N_SCREEN <- 150L
for (f in c("Stage09B_all_validation_attempts.csv", "Stage09B_validation_summary.csv")) if (!file.exists(file.path(previous,
    f))) stop("Missing Stage09B output: ", f)
clean <- function(x) {
    z <- toupper(trimws(as.character(x)))
    z[z %in% c("", "NA", "N/A", "NULL", "UNKNOWN", "NOT AVAILABLE", "[NOT AVAILABLE]", "[NOT EVALUATED]")] <- NA_character_
    z
}
stage_numeric <- function(x) {
    z <- gsub("^STAGE[[:space:]]*", "", clean(x))
    v <- rep(NA_real_, length(z))
    v[grepl("^IV", z)] <- 4
    v[grepl("^III", z)] <- 3
    v[is.na(v) & grepl("^II", z)] <- 2
    v[is.na(v) & grepl("^I($|A|B|C|[[:space:]])", z)] <- 1
    v
}
make_clinical <- function(cl, cc) {
    d <- data.frame(time = as.numeric(as.character(cl$OS_time)), event = as.integer(as.character(cl$OS)))
    a <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
    a[!is.finite(a) | a < 18 | a > 105] <- NA_real_
    d$age10 <- a/10
    s <- clean(cl$gender)
    s[!s %in% c("FEMALE", "MALE")] <- NA_character_
    if (cc != "OV" && cc != "UCEC" && sum(!is.na(s)) >= 0.8 * nrow(d) && length(unique(na.omit(s))) ==
        2)
        d$male <- as.numeric(s == "MALE")
    p <- if ("ajcc_pathologic_tumor_stage" %in% names(cl))
        stage_numeric(cl$ajcc_pathologic_tumor_stage) else rep(NA_real_, nrow(d))
    q <- if ("clinical_stage" %in% names(cl))
        stage_numeric(cl$clinical_stage) else rep(NA_real_, nrow(d))
    st <- if (sum(is.finite(p)) >= sum(is.finite(q)))
        p else q
    if (sum(is.finite(st)) >= 0.8 * nrow(d) && length(unique(na.omit(st))) >= 2)
        d$advanced <- as.numeric(st >= 3)
    if (cc == "UCEC" && "histological_type" %in% names(cl)) {
        h <- clean(cl$histological_type)
        v <- rep(NA_real_, length(h))
        v[grepl("ENDOMETRIOID", h)] <- 0
        v[grepl("SEROUS|MIXED", h)] <- 1
        if (sum(is.finite(v)) >= 0.8 * nrow(d) && length(unique(na.omit(v))) == 2)
            d$serous_mixed <- v
    }
    d
}
get_fit <- function(form, data) {
    w <- character()
    fit <- withCallingHandlers(tryCatch(survival::coxph(form, data = data, ties = "efron", control = survival::coxph.control(iter.max = 60),
        model = TRUE), error = function(e) e), warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
    })
    if (inherits(fit, "error"))
        return(list(fit = NULL, warning = paste(unique(w), collapse = " | "), reason = conditionMessage(fit)))
    b <- coef(fit)
    vv <- tryCatch(diag(vcov(fit)), error = function(e) NA_real_)
    ok <- length(b) > 0 && all(is.finite(b)) && all(is.finite(vv)) && all(vv > 0) && all(abs(b) < 12)
    # Cox separation warnings are NOT accepted as valid fits, even if coefficients finite.
    sep <- any(grepl("infinite|did not converge|singular|ran out of iterations", w, ignore.case = TRUE))
    list(fit = if (ok && !sep) fit else NULL, warning = paste(unique(w), collapse = " | "), reason = if (!ok) "Invalid/singular coefficient or variance" else if (sep) "Suspected separation/nonconvergence" else "")
}
concord <- function(fit, test) {
    if (is.null(fit))
        return(NA_real_)
    lp <- tryCatch(as.numeric(predict(fit, newdata = test, type = "lp", reference = "zero")), error = function(e) rep(NA_real_,
        nrow(test)))
    if (length(lp) != nrow(test) || any(!is.finite(lp)) || sum(test$event) < 2 || length(unique(lp)) <
        2)
        return(NA_real_)
    tryCatch(as.numeric(survival::concordance(survival::Surv(test$time, test$event) ~ lp, reverse = TRUE)$concordance),
        error = function(e) NA_real_)
}

# Targeted recovery rule: when unpenalized clinical fit is invalid, remove only binary predictors
# with absent event or nonevent strata in TRAINING, one at a time.  Model choice is based
# exclusively on TRAINING outcomes, not held-out outcomes.  Reuse identical reduced covariate set
# for baseline, gene fits and training selection.
repair_covars <- function(tr, covars) {
    changes <- character()
    attempted <- list()
    cur <- covars
    repeat {
        f <- as.formula(paste("survival::Surv(time,event)~", if (length(cur))
            paste(cur, collapse = "+") else "1"))
        model <- get_fit(f, tr)
        if (!is.null(model$fit))
            return(list(covars = cur, model = model, removed = changes, status = "ok"))
        if (!length(cur))
            return(list(covars = cur, model = model, removed = changes, status = "failed"))
        tab <- lapply(cur, function(v) {
            z <- tr[[v]]
            count <- table(factor(z, levels = sort(unique(z))), factor(tr$event, levels = 0:1))
            score <- if (length(unique(z)) == 2)
                min(count) else Inf
            list(score = score, counts = paste(apply(count, 1, paste, collapse = "/"), collapse = ";"))
        })
        names(tab) <- cur
        scores <- vapply(tab, `[[`, numeric(1), "score")
        # First remove binary predictors with the clearest training separation.
        candidates <- names(scores)[is.finite(scores) & scores == 0]
        if (!length(candidates)) {
            # If clinical fit still fails, choose sparsest binary predictor, but flag adaptation.
            candidates <- names(scores)[is.finite(scores)]
        }
        if (!length(candidates))
            return(list(covars = cur, model = model, removed = changes, status = "failed_no_binary_fix"))
        chosen <- candidates[which.min(scores[candidates])]
        changes <- c(changes, paste0(chosen, " (training min count=", scores[chosen], ")"))
        cur <- setdiff(cur, chosen)
    }
}

cc <- "BRCA"
path <- file.path(input, paste0("TCGA_", cc, "_discovery.rds"))
if (!file.exists(path)) stop("Missing: ", path)
o <- readRDS(path)
x <- o$expression
cl <- as.data.table(o$clinical)
stopifnot(is.matrix(x), identical(colnames(x), as.character(cl$sample_id)))
d0 <- make_clinical(cl, cc)
keep <- complete.cases(d0) & is.finite(d0$time) & d0$time > 0 & d0$event %in% c(0, 1)
d <- d0[keep, , drop = FALSE]
covars <- setdiff(names(d), c("time", "event"))
covars <- covars[vapply(d[covars], function(z) length(unique(z)) > 1, logical(1))]
d <- d[c("time", "event", covars)]
x <- x[, keep, drop = FALSE]
foldid <- matrix(NA_integer_, nrow(d), repeats)
for (r in seq_len(repeats)) for (ev in 0:1) {
    ix <- which(d$event == ev)
    foldid[ix, r] <- sample(rep(seq_len(K), length.out = length(ix)))
}
original <- fread(file.path(previous, "Stage09B_all_validation_attempts.csv"))
require_cols <- c("Cancer", "Gene", "Repeat", "Fold", "Analysis", "C_clinical", "C_gene", "Delta_C",
    "Status")
if (!all(require_cols %in% names(original))) stop("Stage09B attempts missing columns: ", paste(setdiff(require_cols,
    names(original)), collapse = ", "))
# Audit intact BRCA folds to ensure exactly the same random allocation and estimators.
audit <- list()
repaired <- list()
repaired_nested <- list()
design <- list()
for (r in seq_len(repeats)) for (k in seq_len(K)) {
    tridx <- which(foldid[, r] != k)
    teidx <- which(foldid[, r] == k)
    tr <- d[tridx, , drop = FALSE]
    te <- d[teidx, , drop = FALSE]
    repinfo <- repair_covars(tr, covars)
    cur <- repinfo$covars
    f1 <- as.formula(paste("survival::Surv(time,event)~", paste(c("z", cur), collapse = "+")))
    bc <- concord(repinfo$model$fit, te)
    prior <- original[Cancer == "BRCA" & Analysis == "Discovered_gene_conditional" & Repeat == r & Fold ==
        k]
    if (nrow(prior) != length(fixed_genes))
        stop("Expected 9 earlier BRCA gene results in fold ", r, "/", k)
    priorvalid <- all(prior$Status == "ok")
    # On previously successful folds, verify that the same clinical C results reproduce.
    if (priorvalid) {
        err <- max(abs(prior$C_clinical - bc))
        if (!is.finite(err) || err > 1e-08)
            stop("BRCA fold reproduction mismatch: ", r, "/", k, "; max abs C error: ", err)
    }
    design[[length(design) + 1L]] <- data.table(Repeat = r, Fold = k, Train_n = length(tridx), Train_events = sum(tr$event),
        Test_n = length(teidx), Test_events = sum(te$event), Original_covars = paste(covars, collapse = "+"),
        Final_covars = paste(cur, collapse = "+"), Removed = paste(repinfo$removed, collapse = "; "),
        Clinical_C = bc, Original_success = priorvalid, Clinical_status = repinfo$status)
    if (priorvalid)
        next  # retain the saved estimates
    for (g in fixed_genes) {
        row <- data.table(Cancer = cc, Gene = g, Repeat = r, Fold = k, Analysis = "Discovered_gene_conditional",
            C_clinical = bc, C_gene = NA_real_, Delta_C = NA_real_, Status = "")
        if (!g %in% rownames(x))
            row$Status <- "gene_absent" else {
            v <- as.numeric(x[g, ])
            mu <- mean(v[tridx])
            sdv <- sd(v[tridx])
            if (any(!is.finite(v)) || !is.finite(sdv) || sdv < 1e-09)
                row$Status <- "invalid_expression" else {
                tr1 <- tr
                te1 <- te
                tr1$z <- (v[tridx] - mu)/sdv
                te1$z <- (v[teidx] - mu)/sdv
                fit <- get_fit(f1, tr1)
                gc <- concord(fit$fit, te1)
                row$C_gene <- gc
                row$Delta_C <- gc - bc
                row$Status <- if (is.finite(row$Delta_C))
                  "ok" else paste0("failed: ", fit$reason, "; clinical=", repinfo$status)
            }
        }
        repaired[[length(repaired) + 1L]] <- row
    }
    nn <- data.table(Cancer = cc, Gene = NA_character_, Repeat = r, Fold = k, Analysis = "Nested_train_selected",
        C_clinical = bc, C_gene = NA_real_, Delta_C = NA_real_, Status = "", Candidates_screened = 0L)
    if (is.finite(bc)) {
        mu <- rowMeans(x[, tridx, drop = FALSE])
        vv <- rowMeans(x[, tridx, drop = FALSE]^2) - mu^2
        vv[vv < 1e-12 | !is.finite(vv) | !is.finite(mu)] <- NA_real_
        eligible <- which(is.finite(vv))
        eligible <- eligible[!grepl("^[0-9]+$", rownames(x)[eligible])]
        eligible <- eligible[!duplicated(rownames(x)[eligible])]
        top <- head(eligible[order(vv[eligible], decreasing = TRUE)], N_SCREEN)
        bestp <- Inf
        bestg <- NA_character_
        bestfit <- NULL
        bestmu <- NA_real_
        bestsd <- NA_real_
        ntested <- 0L
        for (j in top) {
            v <- as.numeric(x[j, ])
            ss <- sqrt(vv[j])
            tt <- tr
            tt$z <- (v[tridx] - mu[j])/ss
            fit <- get_fit(f1, tt)
            if (is.null(fit$fit))
                next
            sm <- tryCatch(summary(fit$fit)$coefficients, error = function(e) NULL)
            if (is.null(sm) || !("z" %in% rownames(sm)))
                next
            p <- unname(sm["z", "Pr(>|z|)"])
            ntested <- ntested + 1L
            if (is.finite(p) && p < bestp) {
                bestp <- p
                bestg <- rownames(x)[j]
                bestfit <- fit$fit
                bestmu <- mu[j]
                bestsd <- ss
            }
        }
        nn$Candidates_screened <- ntested
        if (!is.null(bestfit)) {
            te$z <- (as.numeric(x[bestg, teidx]) - bestmu)/bestsd
            gc <- concord(bestfit, te)
            nn$Gene <- bestg
            nn$C_gene <- gc
            nn$Delta_C <- gc - bc
            nn$Status <- if (is.finite(gc))
                "ok" else "score_failed"
        } else nn$Status <- "no_valid_training_gene"
    } else nn$Status <- "clinical_failure"
    repaired_nested[[length(repaired_nested) + 1L]] <- nn
}
design <- rbindlist(design)
# Keep every successful original Stage09B evaluation untouched; replace only failed BRCA rows.
orig_bad <- original[Cancer == "BRCA" & Analysis == "Discovered_gene_conditional" & Status != "ok"]
orig_bad_nested <- original[Cancer == "BRCA" & Analysis == "Nested_train_selected" & Status != "ok"]
if (nrow(orig_bad) != 18L || nrow(orig_bad_nested) != 2L) stop("Unexpected Stage09B BRCA failures: fixed=",
    nrow(orig_bad), "; training-selected=", nrow(orig_bad_nested))
newfixed <- if (length(repaired)) rbindlist(repaired) else original[0]
if (nrow(newfixed) != nrow(orig_bad)) stop("Targeted repair count mismatch: expected ", nrow(orig_bad),
    " got ", nrow(newfixed))
fixed_all <- rbindlist(list(original[Analysis == "Discovered_gene_conditional" & !(Cancer == "BRCA" &
    Status != "ok")], newfixed), fill = TRUE)
previous_nested <- original[Analysis == "Nested_train_selected" & !(Cancer == "BRCA" & Status != "ok")]
newnested <- if (length(repaired_nested)) rbindlist(repaired_nested, fill = TRUE) else original[0]
nested_all <- rbindlist(list(previous_nested, newnested), fill = TRUE)
all <- rbindlist(list(fixed_all, nested_all), fill = TRUE)
setorder(all, Analysis, Cancer, Repeat, Fold, Gene)
# Exactly 900 fixed and 100 training-selected attempts, with explicit failures.
if (nrow(fixed_all) != 900L || nrow(nested_all) != 100L) stop("Unexpected output row counts")
fwrite(all, file.path(out, "Stage09C_combined_attempts.csv"))
fwrite(design, file.path(out, "Stage09C_BRCA_fold_design_audit.csv"))
# Aggregate nested selection by FOLD first, never by selected gene: all fold results, including
# negative results and gene switches, contribute once per cancer/fold.
algorithm <- nested_all[, .(Expected = .N, Successful = sum(Status == "ok"), Failed = sum(Status != "ok"),
    Clinical_C = if (any(Status == "ok")) mean(C_clinical[Status == "ok"]) else NA_real_, Gene_C = if (any(Status ==
        "ok")) mean(C_gene[Status == "ok"]) else NA_real_, Mean_Delta_C = if (any(Status == "ok")) mean(Delta_C[Status ==
        "ok"]) else NA_real_, Median_Delta_C = if (any(Status == "ok")) median(Delta_C[Status == "ok"]) else NA_real_,
    Positive_folds = sum(Delta_C[Status == "ok"] > 0), Min_Delta_C = if (any(Status == "ok")) min(Delta_C[Status ==
        "ok"]) else NA_real_, Max_Delta_C = if (any(Status == "ok")) max(Delta_C[Status == "ok"]) else NA_real_),
    by = Cancer]
setorder(algorithm, -Mean_Delta_C)
fwrite(algorithm, file.path(out, "Stage09C_training_selected_algorithm_per_cancer.csv"))
fixed_summary <- fixed_all[, .(Expected = .N, Successful = sum(Status == "ok"), Failed = sum(Status !=
    "ok"), Clinical_C = if (any(Status == "ok")) mean(C_clinical[Status == "ok"]) else NA_real_, Gene_C = if (any(Status ==
    "ok")) mean(C_gene[Status == "ok"]) else NA_real_, Mean_Delta_C = if (any(Status == "ok")) mean(Delta_C[Status ==
    "ok"]) else NA_real_, Positive_folds = sum(Delta_C[Status == "ok"] > 0)), by = .(Cancer, Gene)]
setorder(fixed_summary, -Mean_Delta_C)
fwrite(fixed_summary, file.path(out, "Stage09C_discovered_genes_per_cancer.csv"))
selected <- nested_all[Status == "ok", .(Fold_selections = .N), by = .(Cancer, Gene)]
selected[, Selection_frequency := Fold_selections/10]
setorder(selected, Cancer, -Fold_selections)
fwrite(selected, file.path(out, "Stage09C_training_selected_gene_frequency.csv"))
fail <- all[Status != "ok", .(Failures = .N), by = .(Analysis, Cancer, Status)]
fwrite(fail, file.path(out, "Stage09C_remaining_failures.csv"))
# Deliberate original-vs-repaired audit for comparable folds, not artificial performance
# improvement.
original_summary <- original[Analysis == "Nested_train_selected" & Status == "ok", .(Mean_Delta_C_09B = mean(Delta_C)),
    by = Cancer]
algorithm <- merge(algorithm, original_summary, by = "Cancer", all.x = TRUE)
setorder(algorithm, -Mean_Delta_C)
fwrite(algorithm, file.path(out, "Stage09C_training_selected_algorithm_per_cancer.csv"))

pub <- function() ggplot2::theme_minimal(base_size = 16) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold",
    colour = "black", size = 19), axis.title = ggplot2::element_text(face = "bold", colour = "black"),
    axis.text = ggplot2::element_text(face = "bold", colour = "black"), panel.grid.minor = ggplot2::element_blank())
p <- ggplot2::ggplot(algorithm, ggplot2::aes(x = Mean_Delta_C, y = reorder(Cancer, Mean_Delta_C))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey45") + ggplot2::geom_segment(ggplot2::aes(x = 0,
    xend = Mean_Delta_C, yend = Cancer), linewidth = 1.3, colour = "#527D8B") + ggplot2::geom_point(size = 4,
    colour = "#007C79") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.3f  (%d/%d)", Mean_Delta_C,
    Successful, Expected)), hjust = ifelse(algorithm$Mean_Delta_C >= 0, -0.09, 1.09), size = 4, fontface = "bold") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.3, 0.35))) + ggplot2::labs(title = "Training-selected prognostic algorithm by cancer",
    subtitle = "Mean held-out \316\224C across all folds; genes reselected within training folds", x = "Mean held-out concordance gain (\316\224C)",
    y = "Cancer") + pub()
fig <- file.path(figdir, "Stage09C_training_selected_algorithm_deltaC.webp")
ggplot2::ggsave(fig, p, device = ragg::agg_webp, width = 13, height = 8, units = "in", res = 300, background = "white")
report <- capture.output({
    cat("==================================================================\nPROJECT 3 | STAGE 09C \342\200\224 BRCA REPAIR & ALGORITHM VALIDATION\n==================================================================\n")
    cat("1. BRCA TRAINING-ONLY DESIGN AUDIT\n")
    print(design)
    cat("\n2. COMPLETENESS\n")
    print(all[, .(Expected = .N, Successful = sum(Status == "ok"), Failed = sum(Status != "ok")), by = Analysis])
    cat("\n3. TRAINING-SELECTED ALGORITHM \342\200\224 PER CANCER (all selected genes combined)\n")
    print(algorithm)
    cat("\n4. TOP DISCOVERY-CONDITIONED ASSOCIATIONS\n")
    print(head(fixed_summary, 15))
    cat("\n5. SELECTION STABILITY (top 20)\n")
    print(head(selected, 20))
    cat("\n6. UNRESOLVED FOLDS\n")
    print(fail)
    cat("\n7. FIGURE\n")
    print(fig)
    cat("\n8. INTERPRETATION\n")
    cat("09C repaired BRCA only, reusing all successful Stage09B folds without alteration.\n")
    cat("Adaptive clinical covariate removal, when needed, uses TRAINING fold event counts; models differ by fold and differ from Stage07.\n")
    cat("Training-only gene discovery: one selected among 150 highest TRAINING-variance eligible genes; not genome-wide screening.\n")
    cat("Report the ALGORITHM mean across all folds, NOT the mean for each repeatedly selected gene as overall predictive performance.\n")
    cat("Preselected 9-gene associations remain discovery-conditioned; outcome-informed selection bias persists.\n")
    cat("These are correlated, repeated held-out CV folds, not external validation; no clinical calibration or prospective validation.\n")
    cat("09B output and source datasets were not modified.\n")
    cat("\nOUTPUT: ", out, "\n==================================================================\n",
        sep = "")
})
writeLines(report, file.path(out, "Stage09C_scientific_checkpoint.txt"))
cat(paste(report, collapse = "\n"), "\n")
