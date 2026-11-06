# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "metafor", "sandwich"))
suppressPackageStartupMessages(library(data.table))

# CRC-ASSOCIATED MICROBIAL FUNCTIONAL MODULE DISCOVERY
#
# Functional potential, not measured metabolic activity.
#
# HC3 heteroskedasticity-robust regression
# Study + age + sex + BMI + mapped sequencing reads
# BH correction across all tested modules

library(data.table)


project_dir <- project_root()

data_dir <- file.path(project_dir, "data")

output_dir <- file.path(
  project_dir,
  "results",
  "functional_associations"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


metadata <- fread(file.path(
  project_dir,
  "results",
  "cohort",
  "metadata_eligible_primary.csv"
))

functional <- fread(file.path(
  data_dir,
  "functional_signal_2340_CRC_cohort_20240329.tab"
))

definitions <- fread(file.path(
  data_dir,
  "modules_definition_GMM_GBM_KEGG_92_20240329.tab"
))

metadata <- metadata[
  `class` %in% c("healthy", "CRC")
]

stopifnot(
  !anyDuplicated(metadata$sample),
  !anyDuplicated(functional$modules_id),
  all(metadata$sample %in% names(functional))
)


sample_ids <- metadata$sample

X <- as.matrix(functional[, ..sample_ids])

storage.mode(X) <- "double"

rownames(X) <- functional$modules_id

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0)
)


prevalence <- rowMeans(X > 0)

# Retain modules detected in at least 10% of samples.
# Filtering is independent of disease labels.

keep <- prevalence >= 0.10


P <- sweep(
  X,
  2,
  colSums(X),
  "/"
)

replace_zeros <- function(p) {
  
  zeros <- p == 0
  z <- sum(zeros)
  
  if (z == 0) return(p)
  
  positive <- p[p > 0]
  
  delta <- min(positive) * 0.5
  
  delta <- min(delta, 0.5 / z)
  
  replacement_mass <- z * delta
  
  p[zeros] <- delta
  
  p[!zeros] <- p[!zeros] *
    (1 - replacement_mass)
  
  p
}

P_replaced <- apply(
  P,
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

rownames(CLR_all) <- functional$modules_id
colnames(CLR_all) <- sample_ids
CLR <- CLR_all[keep, , drop = FALSE]


metadata[
  ,
  Disease := as.integer(`class` == "CRC")
]

metadata[
  ,
  Study := factor(study_accession)
]

metadata[
  ,
  Age := suppressWarnings(as.numeric(age))
]

metadata[
  ,
  BMI := suppressWarnings(as.numeric(bmi))
]

metadata[
  ,
  Sex := tolower(trimws(as.character(gender)))
]

metadata[
  !Sex %in% c("male", "female"),
  Sex := NA_character_
]

metadata[
  ,
  Sex := factor(Sex)
]

metadata[
  ,
  Mapped_reads :=
    suppressWarnings(
      as.numeric(gut_mapped_read_count)
    )
]

metadata[
  ,
  Log_mapped_reads := fifelse(
    is.finite(Mapped_reads) & Mapped_reads > 0,
    log10(Mapped_reads),
    NA_real_
  )
]

model_vars <- c(
  "Disease",
  "Study",
  "Age",
  "Sex",
  "BMI",
  "Log_mapped_reads"
)

complete <- complete.cases(
  metadata[, ..model_vars]
)

analysis_metadata <- copy(metadata[complete])

analysis_metadata[
  ,
  Study := droplevels(Study)
]

analysis_metadata[
  ,
  Sex := droplevels(Sex)
]

CLR <- CLR[, complete, drop = FALSE]

stopifnot(
  identical(
    colnames(CLR),
    analysis_metadata$sample
  )
)


design <- model.matrix(
  ~ Disease + Study + Age + Sex +
    BMI + Log_mapped_reads,
  data = analysis_metadata
)

stopifnot(
  qr(design)$rank == ncol(design),
  "Disease" %in% colnames(design)
)

Y <- t(CLR)

stopifnot(
  nrow(Y) == nrow(design)
)


fit <- lm.fit(
  x = design,
  y = Y
)

df_residual <- nrow(design) - fit$rank

disease_index <- match(
  "Disease",
  colnames(design)
)

beta <- as.numeric(
  fit$coefficients[disease_index, ]
)


XtX_inverse <- solve(crossprod(design))

# Leverage
hat_values <- rowSums(
  (design %*% XtX_inverse) * design
)

stopifnot(
  all(is.finite(hat_values)),
  all(hat_values < 1 - 1e-10)
)

# Influence weights for the disease coefficient
unit <- rep(0, ncol(design))
unit[disease_index] <- 1

weights <- as.vector(
  design %*% XtX_inverse %*% unit
)

# HC3 standard errors for all module regressions
adjusted_weights <- weights / (1 - hat_values)

residuals <- as.matrix(fit$residuals)

robust_variance <- colSums(
  sweep(
    residuals,
    1,
    adjusted_weights,
    "*"
  )^2
)

robust_se <- sqrt(robust_variance)

stopifnot(
  length(beta) == ncol(Y),
  all(is.finite(robust_se)),
  all(robust_se > 0)
)

# Large-sample HC3 Wald inference
z_value <- beta / robust_se

p_values <- 2 * pnorm(-abs(z_value))


results <- data.table(
  Module_ID = rownames(CLR),
  Prevalence_percent =
    as.numeric(prevalence[keep]) * 100,
  CLR_difference = beta,
  Robust_SE = robust_se,
  CI_lower = beta - 1.96 * robust_se,
  CI_upper = beta + 1.96 * robust_se,
  P_value = p_values
)

results[
  ,
  FDR := p.adjust(P_value, method = "BH")
]

results[
  ,
  Direction := fifelse(
    CLR_difference > 0,
    "CRC-enriched",
    "CRC-depleted"
  )
]


annotations <- definitions[
  ,
  .(
    Module_ID = modules_id,
    Category = category,
    Biological_implication = implication
  )
]

stopifnot(
  !anyDuplicated(annotations$Module_ID)
)

results <- merge(
  results,
  annotations,
  by = "Module_ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  results,
  FDR,
  P_value
)

significant <- results[
  !is.na(FDR) & FDR < 0.05
]


category_summary <- results[
  ,
  .(
    Modules_tested = .N,
    Significant = sum(FDR < 0.05, na.rm = TRUE),
    CRC_enriched = sum(
      FDR < 0.05 &
        Direction == "CRC-enriched",
      na.rm = TRUE
    ),
    CRC_depleted = sum(
      FDR < 0.05 &
        Direction == "CRC-depleted",
      na.rm = TRUE
    )
  ),
  by = Category
]

setorder(
  category_summary,
  -Significant
)


fwrite(
  results,
  file.path(
    output_dir,
    "CRC_functional_modules_all.csv"
  )
)

fwrite(
  significant,
  file.path(
    output_dir,
    "CRC_functional_modules_FDR05.csv"
  )
)

fwrite(
  category_summary,
  file.path(
    output_dir,
    "CRC_functional_category_summary.csv"
  )
)

# CROSS-COHORT FUNCTIONAL MODULE ASSOCIATIONS
#
# No repeated abundance processing.

library(data.table)

library(metafor)


output_dir <- file.path(project_dir, "results", "functional_associations")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


top_enriched <- head(
  results[
    Direction == "CRC-enriched"
  ][order(FDR)],
  8
)

top_depleted <- head(
  results[
    Direction == "CRC-depleted"
  ][order(FDR)],
  8
)

# Include biologically relevant metabolic modules
# independently of their discovery significance.

metabolic_modules <- results[
  grepl(
    "butyrate|butyryl|propionate|indole|tryptophan",
    Biological_implication,
    ignore.case = TRUE
  )
]

candidates <- unique(
  rbindlist(
    list(
      top_enriched,
      top_depleted,
      metabolic_modules
    ),
    fill = TRUE
  ),
  by = "Module_ID"
)

candidate_ids <- as.character(candidates$Module_ID)

# Verify feature identifiers and sample alignment
expected_ids <- functional$modules_id[keep]

stopifnot(
  nrow(CLR) == length(expected_ids),
  ncol(CLR) == nrow(analysis_metadata),
  !anyDuplicated(expected_ids)
)

if (is.null(rownames(CLR))) {
  rownames(CLR) <- expected_ids
}

stopifnot(
  identical(rownames(CLR), expected_ids),
  all(candidate_ids %in% rownames(CLR)),
  identical(
    colnames(CLR),
    analysis_metadata$sample
  )
)


fit_module <- function(module_id, study_id) {
  
  d <- copy(
    analysis_metadata[
      study_accession == study_id
    ]
  )
  
  counts <- table(
    factor(d$Disease, levels = c(0, 1))
  )
  
  # Require both comparison groups
  if (any(counts < 10)) return(NULL)
  
  idx <- match(
    d$sample,
    colnames(CLR)
  )
  
  stopifnot(!anyNA(idx))
  
  d[
    ,
    Module_CLR := as.numeric(
      CLR[module_id, idx]
    )
  ]
  
  d[, Sex := droplevels(Sex)]
  
  # Retain only estimable clinical covariates.
  potential <- c(
    "Age",
    "Sex",
    "BMI",
    "Log_mapped_reads"
  )
  
  usable <- potential[
    vapply(
      potential,
      function(v) {
        
        x <- d[[v]]
        
        if (is.factor(x)) {
          nlevels(x) > 1
        } else {
          uniqueN(x) > 1
        }
        
      },
      logical(1)
    )
  ]
  
  formula <- reformulate(
    c("Disease", usable),
    response = "Module_CLR"
  )
  
  fit <- tryCatch(
    lm(formula, data = d),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  if (!"Disease" %in% names(coef(fit))) {
    return(NULL)
  }
  
  if (is.na(coef(fit)["Disease"])) {
    return(NULL)
  }
  
  # HC3 heteroskedasticity-robust standard error
  V <- tryCatch(
    sandwich::vcovHC(fit, type = "HC3"),
    error = function(e) NULL
  )
  
  if (is.null(V)) return(NULL)
  
  estimate <- unname(coef(fit)["Disease"])
  
  se <- sqrt(V["Disease", "Disease"])
  
  if (!is.finite(se) || se <= 0) {
    return(NULL)
  }
  
  data.table(
    Module_ID = module_id,
    Study = study_id,
    N_healthy = unname(counts["0"]),
    N_CRC = unname(counts["1"]),
    Estimate = estimate,
    SE = se,
    CI_lower = estimate - 1.96 * se,
    CI_upper = estimate + 1.96 * se
  )
}

studies <- unique(
  analysis_metadata$study_accession
)

cohort_results <- rbindlist(
  lapply(candidate_ids, function(id) {
    
    rbindlist(
      lapply(studies, function(study) {
        fit_module(id, study)
      }),
      fill = TRUE
    )
    
  }),
  fill = TRUE
)


meta_results <- rbindlist(
  lapply(candidate_ids, function(id) {
    
    d <- cohort_results[
      Module_ID == id
    ]
    
    if (nrow(d) < 3) return(NULL)
    
    fit <- metafor::rma(
      yi = d$Estimate,
      sei = d$SE,
      method = "REML",
      test = "knha"
    )
    
    data.table(
      Module_ID = id,
      Studies = nrow(d),
      
      Positive_studies = sum(
        d$Estimate > 0
      ),
      
      Negative_studies = sum(
        d$Estimate < 0
      ),
      
      Meta_estimate = as.numeric(
        fit$b[1, 1]
      ),
      
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


meta_results <- merge(
  meta_results,
  candidates[
    ,
    .(
      Module_ID,
      Category,
      Biological_implication,
      Discovery_CLR = CLR_difference,
      Discovery_FDR = FDR
    )
  ],
  by = "Module_ID",
  all.x = TRUE,
  sort = FALSE
)

setorder(
  meta_results,
  FDR
)


fwrite(
  cohort_results,
  file.path(
    output_dir,
    "CRC_functional_cohort_effects.csv"
  )
)

fwrite(
  meta_results,
  file.path(
    output_dir,
    "CRC_functional_meta_analysis.csv"
  )
)
