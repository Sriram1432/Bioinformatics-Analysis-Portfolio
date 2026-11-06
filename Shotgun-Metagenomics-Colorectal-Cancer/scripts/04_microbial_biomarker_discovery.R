# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "metafor", "sandwich"))
suppressPackageStartupMessages(library(data.table))

# Discovery of CRC-associated microbial species

library(data.table)


project_dir <- project_root()

data_dir <- file.path(project_dir, "data")
cohort_dir <- file.path(project_dir, "results", "cohort")
profiling_dir <- file.path(project_dir, "results", "profiling")
output_dir <- file.path(project_dir, "results", "microbial_associations")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


metadata <- fread(file.path(
  cohort_dir, "metadata_eligible_primary.csv"
))

species <- fread(file.path(
  data_dir, "species_signal_2340_CRC_cohort_20240322.tab"
))

taxonomy <- fread(file.path(
  data_dir,
  "MSP_set_ref_oss_gss_status_612_ref_paper_20250220_gtdb_220.tab"
))


metadata <- metadata[
  `class` %in% c("healthy", "CRC")
]

sample_ids <- metadata$sample

stopifnot(
  !anyDuplicated(sample_ids),
  all(sample_ids %in% names(species)),
  !anyDuplicated(species$msp_id)
)

X <- as.matrix(species[, ..sample_ids])
storage.mode(X) <- "double"
rownames(X) <- species$msp_id

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0)
)


prevalence <- rowMeans(X > 0)

# Define filtering independently of disease labels.
keep <- prevalence >= 0.10

X <- X[keep, , drop = FALSE]


# Start with all original MSP features.
# The retained subset is used for association testing only.
# CLR references the complete measured community.

X_all <- as.matrix(species[, ..sample_ids])
storage.mode(X_all) <- "double"

P_all <- sweep(X_all, 2, colSums(X_all), "/")

# Multiplicative zero replacement:
# Assign a small positive fraction to zeros within each sample,
# scaling nonzero components to preserve closure.

replace_zeros <- function(p) {
  
  zeros <- p == 0
  z <- sum(zeros)
  
  if (z == 0) return(p)
  
  positive <- p[p > 0]
  delta <- min(positive) * 0.5
  
  # Ensure replacement mass is less than total composition.
  delta <- min(delta, 0.5 / z)
  
  replacement_mass <- z * delta
  
  p[zeros] <- delta
  p[!zeros] <- p[!zeros] *
    (1 - replacement_mass)
  
  p
}

P_replaced <- apply(
  P_all,
  2,
  replace_zeros
)

stopifnot(
  all(P_replaced > 0),
  max(abs(colSums(P_replaced) - 1)) < 1e-8
)

logP <- log(P_replaced)

CLR_all <- sweep(
  logP,
  2,
  colMeans(logP),
  "-"
)

rownames(CLR_all) <- species$msp_id
colnames(CLR_all) <- sample_ids
CLR <- CLR_all[keep, , drop = FALSE]


metadata[, Disease := as.integer(`class` == "CRC")]
metadata[, Study := factor(study_accession)]
metadata[, Age := suppressWarnings(as.numeric(age))]
metadata[, BMI := suppressWarnings(as.numeric(bmi))]

metadata[, Sex := tolower(trimws(as.character(gender)))]
metadata[!Sex %in% c("male", "female"), Sex := NA_character_]
metadata[, Sex := factor(Sex)]

model_vars <- c("Disease", "Study", "Age", "Sex", "BMI")

complete <- complete.cases(metadata[, ..model_vars])

analysis_metadata <- metadata[complete]
analysis_metadata[, Study := droplevels(Study)]
analysis_metadata[, Sex := droplevels(Sex)]

CLR <- CLR[, complete, drop = FALSE]


design <- model.matrix(
  ~ Disease + Study + Age + Sex + BMI,
  data = analysis_metadata
)

stopifnot(
  qr(design)$rank == ncol(design),
  "Disease" %in% colnames(design)
)

# Features x samples -> samples x features
Y <- t(CLR)

fit <- lm.fit(
  x = design,
  y = Y
)

df_residual <- nrow(design) - fit$rank

stopifnot(df_residual > 0)

coefficient_index <- match(
  "Disease",
  colnames(design)
)

beta <- fit$coefficients[coefficient_index, ]

# Homoskedastic OLS covariance:
# Used only for preliminary association screening.
R <- fit$residuals

sigma2 <- colSums(R^2) / df_residual

XtX_inverse <- chol2inv(chol(crossprod(design)))

se <- sqrt(
  sigma2 * XtX_inverse[coefficient_index, coefficient_index]
)

t_stat <- beta / se

p_values <- 2 * pt(
  -abs(t_stat),
  df = df_residual
)


results <- data.table(
  Species_ID = rownames(CLR),
  Prevalence_percent = prevalence[keep] * 100,
  CLR_difference = as.numeric(beta),
  Standard_error = as.numeric(se),
  CI_lower = as.numeric(
    beta - qt(0.975, df_residual) * se
  ),
  CI_upper = as.numeric(
    beta + qt(0.975, df_residual) * se
  ),
  P_value = as.numeric(p_values)
)

results[, FDR := p.adjust(P_value, method = "BH")]

results[, Direction := fifelse(
  CLR_difference > 0,
  "CRC-enriched",
  "CRC-depleted"
)]


annotations <- taxonomy[
  ,
  .(
    Species_ID = msp_id,
    Scientific_name = species,
    Genus = genus,
    Phylum = phylum,
    Domain = superkingdom
  )
]

results <- merge(
  results,
  annotations,
  by = "Species_ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(results, FDR, P_value)


significant <- results[
  !is.na(FDR) & FDR < 0.05
]

discovery_summary <- data.table(
  Metric = c(
    "CRC samples",
    "Healthy samples",
    "Complete-case samples",
    "Features before filtering",
    "Features tested",
    "FDR significant features",
    "CRC-enriched significant features",
    "CRC-depleted significant features"
  ),
  Value = c(
    sum(metadata$Disease == 1),
    sum(metadata$Disease == 0),
    nrow(analysis_metadata),
    nrow(X_all),
    nrow(results),
    nrow(significant),
    sum(significant$Direction == "CRC-enriched"),
    sum(significant$Direction == "CRC-depleted")
  )
)


fwrite(
  results,
  file.path(output_dir, "CRC_species_discovery_all.csv")
)

fwrite(
  significant,
  file.path(output_dir, "CRC_species_discovery_FDR05.csv")
)

fwrite(
  discovery_summary,
  file.path(output_dir, "CRC_species_discovery_summary.csv")
)

# Cross-project synthesis of CRC-associated MSPs
#
# These source projects contributed to discovery and are not treated as
# independent external replication.

packages <- c("data.table", "sandwich", "metafor")

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

library(data.table)
library(sandwich)
library(metafor)


project_dir <- project_root()

output_dir <- file.path(
  project_dir, "results", "microbial_associations"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Select 20 strongest exploratory candidates
candidates <- head(
  results[order(FDR, P_value)],
  20
)

# Preserve the original analysis sample order
stopifnot(
  identical(
    colnames(CLR),
    analysis_metadata$sample
  )
)


fit_candidate <- function(feature_id, study_id) {
  
  d <- copy(
    analysis_metadata[study_accession == study_id]
  )
  
  counts <- table(
    factor(d$Disease, levels = c(0, 1))
  )
  
  if (any(counts < 10)) return(NULL)
  
  sample_index <- match(d$sample, colnames(CLR))
  
  d[, Abundance_CLR := as.numeric(
    CLR[feature_id, sample_index]
  )]
  
  d[, Sex := droplevels(Sex)]
  
  # All included observations already passed
  # the common complete-case criteria.
  
  fit <- lm(
    Abundance_CLR ~ Disease + Age + Sex + BMI,
    data = d
  )
  
  if (anyNA(coef(fit))) return(NULL)
  
  V <- sandwich::vcovHC(fit, type = "HC3")
  
  estimate <- unname(coef(fit)["Disease"])
  
  se <- sqrt(V["Disease", "Disease"])
  
  if (!is.finite(se) || se <= 0) return(NULL)
  
  data.table(
    Species_ID = feature_id,
    Study = study_id,
    N_healthy = unname(counts["0"]),
    N_CRC = unname(counts["1"]),
    Estimate = estimate,
    SE = se,
    CI_lower = estimate - 1.96 * se,
    CI_upper = estimate + 1.96 * se
  )
}

studies <- unique(analysis_metadata$study_accession)

cohort_results <- rbindlist(
  lapply(candidates$Species_ID, function(id) {
    
    rbindlist(
      lapply(studies, function(study) {
        fit_candidate(id, study)
      }),
      fill = TRUE
    )
    
  }),
  fill = TRUE
)


meta_results <- rbindlist(
  lapply(candidates$Species_ID, function(id) {
    
    d <- cohort_results[Species_ID == id]
    
    if (nrow(d) < 3) return(NULL)
    
    fit <- metafor::rma(
      yi = d$Estimate,
      sei = d$SE,
      method = "REML",
      test = "knha"
    )
    
    data.table(
      Species_ID = id,
      Studies = nrow(d),
      Positive_studies = sum(d$Estimate > 0),
      Negative_studies = sum(d$Estimate < 0),
      Meta_estimate = as.numeric(fit$b[1, 1]),
      CI_lower = fit$ci.lb,
      CI_upper = fit$ci.ub,
      P_value = fit$pval,
      I2 = fit$I2,
      Tau2 = fit$tau2
    )
    
  }),
  fill = TRUE
)

meta_results[
  ,
  FDR := p.adjust(P_value, method = "BH")
]

# Attach organism names
meta_results <- merge(
  meta_results,
  unique(candidates[, .(
    Species_ID,
    Scientific_name
  )]),
  by = "Species_ID",
  all.x = TRUE
)

setorder(meta_results, FDR)


fwrite(
  cohort_results,
  file.path(
    output_dir,
    "CRC_candidate_cohort_effects.csv"
  )
)

fwrite(
  meta_results,
  file.path(
    output_dir,
    "CRC_candidate_meta_analysis.csv"
  )
)

# Rapid robustness testing of CRC-associated microorganisms
#
#
# Reuses existing results. No abundance recalculation.

suppressPackageStartupMessages(library(metafor))


output_dir <- file.path(project_dir, "results", "microbial_associations")

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

candidate_ids <- as.character(candidates$Species_ID)

stopifnot(
  length(candidate_ids) == 20,
  !anyDuplicated(candidate_ids),
  all(candidate_ids %in% species$msp_id),
  all(candidate_ids %in% cohort_results$Species_ID)
)


loo_details <- rbindlist(
  lapply(candidate_ids, function(id) {
    
    d <- cohort_results[Species_ID == id]
    
    stopifnot(
      nrow(d) >= 4,
      !anyDuplicated(d$Study)
    )
    
    rbindlist(
      lapply(d$Study, function(excluded_study) {
        
        sub <- d[Study != excluded_study]
        
        fit <- metafor::rma(
          yi = sub$Estimate,
          sei = sub$SE,
          method = "REML",
          test = "knha"
        )
        
        data.table(
          Species_ID = id,
          Excluded_study = excluded_study,
          Studies_remaining = nrow(sub),
          Estimate = as.numeric(fit$b[1, 1]),
          CI_lower = fit$ci.lb,
          CI_upper = fit$ci.ub,
          P_value = fit$pval,
          I2 = fit$I2
        )
        
      }),
      fill = TRUE
    )
    
  }),
  fill = TRUE
)


loo_summary <- loo_details[
  ,
  .(
    Studies_removed = .N,
    
    Minimum_effect = min(Estimate),
    Maximum_effect = max(Estimate),
    
    Maximum_P = max(P_value),
    
    Significant_removals = sum(
      CI_lower > 0 | CI_upper < 0
    ),
    
    Consistent_direction = all(Estimate > 0) |
      all(Estimate < 0),
    
    Maximum_I2 = max(I2)
  ),
  by = Species_ID
]

loo_summary[
  ,
  Robust_LOO := (
    Consistent_direction &
      Significant_removals == Studies_removed
  )
]


# Use the original CRC/healthy metadata.
# This includes samples omitted from complete-case regression.

meta_prev <- copy(metadata)

meta_prev <- meta_prev[
  `class` %in% c("healthy", "CRC")
]

stopifnot(
  !anyDuplicated(meta_prev$sample),
  all(meta_prev$sample %in% colnames(X_all))
)

# Explicit sample alignment
sample_indices <- match(
  meta_prev$sample,
  colnames(X_all)
)

stopifnot(
  !anyNA(sample_indices),
  identical(
    colnames(X_all)[sample_indices],
    meta_prev$sample
  )
)

# MSP identifiers correspond to original species row order
feature_indices <- match(
  candidate_ids,
  species$msp_id
)

stopifnot(
  !anyNA(feature_indices),
  nrow(X_all) == nrow(species)
)

# Detection/non-detection matrix
detected <- X_all[
  feature_indices,
  sample_indices,
  drop = FALSE
] > 0

rownames(detected) <- candidate_ids
colnames(detected) <- meta_prev$sample

# Retain studies containing both groups
study_counts <- dcast(
  meta_prev[
    ,
    .N,
    by = .(study_accession, `class`)
  ],
  study_accession ~ class,
  value.var = "N",
  fill = 0
)

informative_studies <- study_counts[
  healthy >= 5 & CRC >= 5,
  study_accession
]

meta_prev <- meta_prev[
  study_accession %in% informative_studies
]

detected <- detected[
  ,
  match(meta_prev$sample, colnames(detected)),
  drop = FALSE
]

stopifnot(
  identical(
    colnames(detected),
    meta_prev$sample
  )
)


prevalence_test <- function(id) {
  
  present <- as.logical(detected[id, ])
  
  # Overall descriptive detection frequencies
  healthy <- meta_prev$`class` == "healthy"
  crc <- meta_prev$`class` == "CRC"
  
  healthy_prevalence <- mean(present[healthy]) * 100
  crc_prevalence <- mean(present[crc]) * 100
  
  # Study-specific 2 x 2 tables:
  # Rows = CRC / healthy
  # Columns = detected / not detected
  
  tables <- lapply(
    informative_studies,
    function(study) {
      
      idx <- meta_prev$study_accession == study
      
      group <- factor(
        meta_prev$`class`[idx],
        levels = c("CRC", "healthy")
      )
      
      status <- factor(
        present[idx],
        levels = c(TRUE, FALSE)
      )
      
      table(group, status)
    }
  )
  
  array_3d <- array(
    NA_integer_,
    dim = c(2, 2, length(tables))
  )
  
  for (i in seq_along(tables)) {
    array_3d[, , i] <- tables[[i]]
  }
  
  dimnames(array_3d) <- list(
    Disease = c("CRC", "healthy"),
    Detection = c("Detected", "Not_detected"),
    Study = informative_studies
  )
  
  mh <- tryCatch(
    mantelhaen.test(
      array_3d,
      correct = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(mh)) {
    
    return(data.table(
      Species_ID = id,
      Healthy_prevalence = healthy_prevalence,
      CRC_prevalence = crc_prevalence,
      Prevalence_difference = crc_prevalence -
        healthy_prevalence,
      Odds_ratio = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      P_value = NA_real_
    ))
    
  }
  
  data.table(
    Species_ID = id,
    
    Healthy_prevalence = healthy_prevalence,
    CRC_prevalence = crc_prevalence,
    
    Prevalence_difference =
      crc_prevalence - healthy_prevalence,
    
    Odds_ratio = unname(mh$estimate),
    
    CI_lower = mh$conf.int[1],
    CI_upper = mh$conf.int[2],
    
    P_value = mh$p.value
  )
}

prevalence_results <- rbindlist(
  lapply(candidate_ids, prevalence_test),
  fill = TRUE
)

prevalence_results[
  ,
  FDR := p.adjust(P_value, method = "BH")
]


summary_table <- merge(
  candidates[
    ,
    .(
      Species_ID,
      Scientific_name,
      Discovery_CLR = CLR_difference,
      Discovery_FDR = FDR
    )
  ],
  loo_summary,
  by = "Species_ID"
)

summary_table <- merge(
  summary_table,
  prevalence_results,
  by = "Species_ID"
)

summary_table[
  ,
  Prevalence_direction_agrees :=
    sign(Prevalence_difference) ==
    sign(Discovery_CLR)
]

summary_table[
  ,
  Prevalence_significant :=
    !is.na(FDR) & FDR < 0.05
]

# Restore discovery ranking
summary_table[
  ,
  Candidate_order := match(
    Species_ID,
    candidate_ids
  )
]

setorder(
  summary_table,
  Candidate_order
)


fwrite(
  loo_details,
  file.path(
    output_dir,
    "CRC_leave_one_study_out_all.csv"
  )
)

fwrite(
  loo_summary,
  file.path(
    output_dir,
    "CRC_leave_one_study_out_summary.csv"
  )
)

fwrite(
  prevalence_results,
  file.path(
    output_dir,
    "CRC_prevalence_associations.csv"
  )
)

fwrite(
  summary_table,
  file.path(
    output_dir,
    "CRC_candidate_robustness_summary.csv"
  )
)
