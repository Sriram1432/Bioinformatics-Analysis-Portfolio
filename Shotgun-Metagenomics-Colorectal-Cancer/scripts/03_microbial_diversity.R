# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "metafor", "vegan", "permute", "sandwich", "lmtest"))
suppressPackageStartupMessages(library(data.table))

# Alpha diversity and preliminary cohort-adjusted associations

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})


project_dir <- project_root()

data_dir <- file.path(project_dir, "data")

cohort_dir <- file.path(project_dir, "results", "cohort")

output_dir <- file.path(project_dir, "results", "diversity")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


metadata <- fread(
  file.path(cohort_dir, "metadata_eligible_primary.csv")
)

species <- fread(
  file.path(
    data_dir,
    "species_signal_2340_CRC_cohort_20240322.tab"
  )
)

stopifnot(
  nrow(metadata) == 2138L,
  !anyDuplicated(metadata$sample),
  all(metadata$sample %in% names(species))
)


sample_ids <- metadata$sample

X <- as.matrix(species[, ..sample_ids])

storage.mode(X) <- "double"

rownames(X) <- species$msp_id

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0)
)


# Observed MSP richness
richness <- colSums(X > 0)

# Within-sample proportional MSP signal
P <- sweep(X, 2, colSums(X), "/")

# Shannon diversity
# Zero entries contribute zero by convention.
logP <- matrix(
  0,
  nrow = nrow(P),
  ncol = ncol(P)
)

positive <- P > 0

logP[positive] <- P[positive] * log(P[positive])

shannon <- -colSums(logP)

# Gini-Simpson diversity
simpson <- 1 - colSums(P^2)


dt <- copy(metadata)

dt[, Richness := as.numeric(richness)]
dt[, Shannon := shannon]
dt[, Simpson := simpson]

dt[, Disease := factor(
  `class`,
  levels = c("healthy", "adenoma", "CRC")
)]

dt[, Study := factor(study_accession)]

dt[, Age := suppressWarnings(as.numeric(age))]
dt[, BMI := suppressWarnings(as.numeric(bmi))]

dt[, Sex := tolower(trimws(as.character(gender)))]

dt[
  !Sex %in% c("male", "female"),
  Sex := NA_character_
]

dt[, Sex := factor(Sex)]

dt[, Mapped_reads := as.numeric(gut_mapped_read_count)]

dt[
  ,
  Log_mapped_reads := ifelse(
    is.finite(Mapped_reads) & Mapped_reads > 0,
    log10(Mapped_reads),
    NA_real_
  )
]


diversity_summary <- dt[
  ,
  .(
    N = .N,
    
    Richness_median = median(Richness),
    Richness_Q1 = quantile(Richness, 0.25),
    Richness_Q3 = quantile(Richness, 0.75),
    
    Shannon_mean = mean(Shannon),
    Shannon_SD = sd(Shannon),
    
    Simpson_mean = mean(Simpson),
    Simpson_SD = sd(Simpson)
  ),
  by = Disease
]


depth_associations <- data.table(
  Metric = c("Richness", "Shannon", "Simpson"),
  
  Spearman_rho = c(
    cor(dt$Richness, dt$Mapped_reads,
        method = "spearman", use = "complete.obs"),
    
    cor(dt$Shannon, dt$Mapped_reads,
        method = "spearman", use = "complete.obs"),
    
    cor(dt$Simpson, dt$Mapped_reads,
        method = "spearman", use = "complete.obs")
  )
)


model_variables <- c(
  "Richness",
  "Shannon",
  "Disease",
  "Study",
  "Log_mapped_reads",
  "Age",
  "Sex",
  "BMI"
)

model_data <- dt[
  complete.cases(dt[, ..model_variables])
]

model_data[, Disease := droplevels(Disease)]
model_data[, Study := droplevels(Study)]
model_data[, Sex := droplevels(Sex)]

stopifnot(
  nlevels(model_data$Disease) == 3L,
  nlevels(model_data$Sex) == 2L,
  all(table(model_data$Disease) > 0)
)

model_group_counts <- model_data[
  ,
  .N,
  by = Disease
]


model_formula <- function(outcome) {
  as.formula(
    paste0(
      outcome,
      " ~ Disease + Study + Log_mapped_reads + Age + Sex + BMI"
    )
  )
}

fit_richness <- lm(
  model_formula("Richness"),
  data = model_data
)

fit_shannon <- lm(
  model_formula("Shannon"),
  data = model_data
)


extract_effects <- function(model, outcome) {
  
  V <- vcovHC(model, type = "HC3")
  
  coefficients <- coef(model)
  
  terms <- c(
    "Diseaseadenoma",
    "DiseaseCRC"
  )
  
  stopifnot(
    all(terms %in% names(coefficients)),
    all(is.finite(coefficients[terms])),
    all(is.finite(diag(V)[terms]))
  )
  
  estimates <- coefficients[terms]
  
  standard_errors <- sqrt(diag(V))[terms]
  
  z <- estimates / standard_errors
  
  p_values <- 2 * pnorm(-abs(z))
  
  data.table(
    Outcome = outcome,
    Comparison = c(
      "Adenoma vs healthy",
      "CRC vs healthy"
    ),
    Estimate = as.numeric(estimates),
    Robust_SE = as.numeric(standard_errors),
    CI_lower = as.numeric(
      estimates - qnorm(0.975) * standard_errors
    ),
    CI_upper = as.numeric(
      estimates + qnorm(0.975) * standard_errors
    ),
    P_value = as.numeric(p_values)
  )
}

regression_results <- rbindlist(list(
  extract_effects(fit_richness, "Richness"),
  extract_effects(fit_shannon, "Shannon")
))

# Exploratory BH correction across four comparisons
regression_results[
  ,
  FDR := p.adjust(P_value, method = "BH")
]


cohort_groups <- dcast(
  dt[, .N, by = .(study_accession, Disease)],
  study_accession ~ Disease,
  value.var = "N",
  fill = 0
)

complete_case_studies <- dcast(
  model_data[, .N, by = .(study_accession, Disease)],
  study_accession ~ Disease,
  value.var = "N",
  fill = 0
)


fwrite(
  dt[, .(
    sample,
    study_accession,
    country,
    Disease,
    Richness,
    Shannon,
    Simpson,
    Mapped_reads
  )],
  file.path(output_dir, "alpha_diversity_samples.csv")
)

fwrite(
  diversity_summary,
  file.path(output_dir, "alpha_diversity_group_summary.csv")
)

fwrite(
  regression_results,
  file.path(output_dir, "alpha_diversity_adjusted_models.csv")
)

fwrite(
  cohort_groups,
  file.path(output_dir, "alpha_diversity_cohort_groups.csv")
)

saveRDS(
  list(
    richness = fit_richness,
    shannon = fit_shannon
  ),
  file.path(output_dir, "alpha_diversity_models.rds")
)

# Cross-cohort alpha diversity and random-effects meta-analysis

packages <- c("data.table", "metafor")

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

library(data.table)
library(metafor)

set.seed(20261027)


project_dir <- project_root()

output_dir <- file.path(project_dir, "results", "diversity")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(file.path(
  project_dir,
  "results",
  "cohort",
  "metadata_eligible_primary.csv"
))

diversity <- fread(file.path(
  output_dir,
  "alpha_diversity_samples.csv"
))

stopifnot(
  !anyDuplicated(metadata$sample),
  !anyDuplicated(diversity$sample),
  setequal(metadata$sample, diversity$sample)
)


dt <- merge(
  metadata,
  diversity[, .(sample, Richness, Shannon)],
  by = "sample",
  sort = FALSE
)

dt[, Disease := factor(
  `class`,
  levels = c("healthy", "CRC")
)]

dt[, Age := suppressWarnings(as.numeric(age))]
dt[, BMI := suppressWarnings(as.numeric(bmi))]

dt[, Sex := tolower(trimws(as.character(gender)))]
dt[!Sex %in% c("male", "female"), Sex := NA_character_]
dt[, Sex := factor(Sex)]

dt[, Mapped_reads := as.numeric(gut_mapped_read_count)]
dt[, Log_mapped_reads := fifelse(
  is.finite(Mapped_reads) & Mapped_reads > 0,
  log10(Mapped_reads),
  NA_real_
)]

dt <- dt[`class` %in% c("healthy", "CRC")]


fit_study <- function(data, outcome) {
  
  study <- unique(data$study_accession)
  
  counts <- table(factor(
    data$`class`,
    levels = c("healthy", "CRC")
  ))
  
  result <- data.table(
    Study = study,
    Outcome = outcome,
    N_healthy = unname(counts["healthy"]),
    N_CRC = unname(counts["CRC"]),
    N_model = NA_integer_,
    Estimate = NA_real_,
    SE = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    P_value = NA_real_,
    Covariates = NA_character_,
    Status = "Not estimable"
  )
  
  # Minimum exploratory group size
  if (any(counts < 10)) {
    result[, Status := "Insufficient group size"]
    return(result)
  }
  
  # Begin with the same potential covariates.
  # Remove covariates that are constant or unavailable
  # within a study rather than creating singular models.
  
  candidate <- c(
    "Age", "Sex", "BMI", "Log_mapped_reads"
  )
  
  use <- c(outcome, "Disease", candidate)
  
  d <- data[
    complete.cases(data[, ..use])
  ]
  
  d[, Disease := droplevels(Disease)]
  d[, Sex := droplevels(Sex)]
  
  if (nrow(d) < 20 ||
      nlevels(d$Disease) != 2 ||
      any(table(d$Disease) < 5)) {
    
    result[, Status := "Insufficient complete cases"]
    result[, N_model := nrow(d)]
    return(result)
  }
  
  valid_covariates <- candidate[
    vapply(candidate, function(v) {
      x <- d[[v]]
      if (is.factor(x)) {
        nlevels(x) > 1
      } else {
        uniqueN(x) > 1
      }
    }, logical(1))
  ]
  
  formula <- reformulate(
    c("Disease", valid_covariates),
    response = outcome
  )
  
  model <- tryCatch(
    lm(formula, data = d),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    result[, Status := "Model failed"]
    return(result)
  }
  
  term <- "DiseaseCRC"
  
  if (!term %in% names(coef(model)) ||
      is.na(coef(model)[term])) {
    result[, Status := "Disease effect not estimable"]
    return(result)
  }
  
  sm <- summary(model)$coefficients
  
  estimate <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  p <- unname(sm[term, "Pr(>|t|)"])
  
  if (!is.finite(se) || se <= 0) {
    result[, Status := "Invalid standard error"]
    return(result)
  }
  
  critical <- qt(0.975, df.residual(model))
  
  result[, `:=`(
    N_model = nobs(model),
    Estimate = estimate,
    SE = se,
    CI_lower = estimate - critical * se,
    CI_upper = estimate + critical * se,
    P_value = p,
    Covariates = paste(valid_covariates, collapse = ", "),
    Status = "Estimated"
  )]
  
  result
}

study_results <- rbindlist(
  lapply(c("Richness", "Shannon"), function(outcome) {
    
    rbindlist(
      lapply(
        split(dt, by = "study_accession"),
        fit_study,
        outcome = outcome
      ),
      fill = TRUE
    )
    
  }),
  fill = TRUE
)


meta_models <- list()

meta_summary <- rbindlist(
  lapply(c("Richness", "Shannon"), function(outcome) {
    
    d <- study_results[
      Outcome == outcome & Status == "Estimated"
    ]
    
    if (nrow(d) < 3) {
      return(data.table(
        Outcome = outcome,
        Studies = nrow(d),
        Estimate = NA_real_,
        CI_lower = NA_real_,
        CI_upper = NA_real_,
        P_value = NA_real_,
        Tau2 = NA_real_,
        I2 = NA_real_,
        Q_pvalue = NA_real_
      ))
    }
    
    model <- rma(
      yi = d$Estimate,
      sei = d$SE,
      method = "REML",
      test = "knha",
      slab = d$Study
    )
    
    meta_models[[outcome]] <<- model
    
    data.table(
      Outcome = outcome,
      Studies = model$k,
      Estimate = as.numeric(model$b[1, 1]),
      CI_lower = model$ci.lb,
      CI_upper = model$ci.ub,
      P_value = model$pval,
      Tau2 = model$tau2,
      I2 = model$I2,
      Q_pvalue = model$QEp
    )
  }),
  fill = TRUE
)

meta_summary[
  ,
  FDR := p.adjust(P_value, method = "BH")
]


fwrite(
  study_results,
  file.path(output_dir, "alpha_diversity_cohort_effects.csv")
)

fwrite(
  meta_summary,
  file.path(output_dir, "alpha_diversity_meta_analysis.csv")
)

saveRDS(
  meta_models,
  file.path(output_dir, "alpha_diversity_meta_models.rds")
)

# CORRECTED DISCOVERY BETA DIVERSITY
#
# Correct study grouping: study_accession
# Original discovery MSP profiles
#

library(data.table)
library(vegan)
library(permute)

set.seed(20261027)

project_dir <- project_root()

output_dir <- file.path(
  project_dir,
  "results",
  "community_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- fread(
  file.path(
    project_dir,
    "results",
    "cohort",
    "metadata_eligible_primary.csv"
  )
)

metadata <- metadata[
  class %in% c("CRC", "healthy")
]


species_data <- fread(
  file.path(
    project_dir,
    "data",
    "species_signal_2340_CRC_cohort_20240322.tab"
  ),
  showProgress = FALSE
)

stopifnot(
  !anyDuplicated(metadata$sample),
  all(metadata$sample %in% names(species_data)),
  !anyDuplicated(species_data$msp_id)
)

sample_ids <- metadata$sample

community <- t(
  as.matrix(species_data[, ..sample_ids])
)

storage.mode(community) <- "double"

rownames(community) <- sample_ids
colnames(community) <- species_data$msp_id

stopifnot(
  all(is.finite(community)),
  all(community >= 0),
  all(rowSums(community) > 0)
)

community <- decostand(
  community,
  method = "total"
)


md <- data.table(
  Sample = metadata$sample,
  Study = metadata$study_accession,
  Group = metadata$class,
  Age = suppressWarnings(as.numeric(metadata$age)),
  BMI = suppressWarnings(as.numeric(metadata$bmi)),
  Sex = tolower(trimws(as.character(metadata$gender)))
)

md[
  !Sex %in% c("male", "female"),
  Sex := NA_character_
]

md[, Study := factor(Study)]
md[, Group := factor(Group, levels = c("healthy", "CRC"))]

stopifnot(
  identical(rownames(community), md$Sample),
  !anyNA(md$Study)
)


study_counts <- md[
  ,
  .(
    CRC = sum(Group == "CRC"),
    Healthy = sum(Group == "healthy")
  ),
  by = Study
]

eligible_studies <- study_counts[
  CRC > 0 & Healthy > 0,
  Study
]

keep <- md$Study %in% eligible_studies

md <- droplevels(md[keep])

community <- community[
  keep,
  ,
  drop = FALSE
]


distance <- vegdist(
  community,
  method = "bray"
)

permutations <- permute::how(
  nperm = 999,
  blocks = md$Study
)

primary_fit <- adonis2(
  distance ~ Study + Group,
  data = as.data.frame(md),
  permutations = permutations,
  by = "terms"
)


complete_index <- complete.cases(
  md[, .(Age, Sex, BMI)]
)

md_adjusted <- droplevels(
  as.data.frame(md[complete_index])
)

community_adjusted <- community[
  complete_index,
  ,
  drop = FALSE
]

adjusted_counts <- table(
  md_adjusted$Study,
  md_adjusted$Group
)

# Only retain studies with both groups after adjustment.

adjusted_eligible <- rownames(adjusted_counts)[
  rowSums(adjusted_counts > 0) == 2
]

adjusted_keep <- md_adjusted$Study %in% adjusted_eligible

md_adjusted <- droplevels(
  md_adjusted[adjusted_keep, , drop = FALSE]
)

community_adjusted <- community_adjusted[
  adjusted_keep,
  ,
  drop = FALSE
]

distance_adjusted <- vegdist(
  community_adjusted,
  method = "bray"
)

adjusted_permutations <- permute::how(
  nperm = 999,
  blocks = md_adjusted$Study
)

adjusted_fit <- adonis2(
  distance_adjusted ~ Study + Age + Sex + BMI + Group,
  data = md_adjusted,
  permutations = adjusted_permutations,
  by = "terms"
)


study_results <- list()

for (study in levels(md$Study)) {
  
  idx <- md$Study == study
  
  study_md <- droplevels(
    as.data.frame(md[idx])
  )
  
  study_matrix <- community[
    idx,
    ,
    drop = FALSE
  ]
  
  study_distance <- vegdist(
    study_matrix,
    method = "bray"
  )
  
  fit <- adonis2(
    study_distance ~ Group,
    data = study_md,
    permutations = 999
  )
  
  dispersion <- betadisper(
    study_distance,
    study_md$Group,
    type = "median",
    bias.adjust = TRUE
  )
  
  dispersion_test <- permutest(
    dispersion,
    permutations = 999
  )
  
  study_results[[study]] <- data.table(
    Study = study,
    N = nrow(study_md),
    CRC_N = sum(study_md$Group == "CRC"),
    Healthy_N = sum(study_md$Group == "healthy"),
    R2 = fit$R2[1],
    PERMANOVA_p = fit$`Pr(>F)`[1],
    Dispersion_p =
      dispersion_test$tab[1, "Pr(>F)"]
  )
}

study_results <- rbindlist(study_results)

study_results[
  ,
  PERMANOVA_FDR := p.adjust(
    PERMANOVA_p,
    method = "BH"
  )
]

study_results[
  ,
  Dispersion_FDR := p.adjust(
    Dispersion_p,
    method = "BH"
  )
]


fwrite(
  study_results,
  file.path(
    output_dir,
    "community_validation_4_corrected_study_results.csv"
  )
)

saveRDS(
  list(
    primary = primary_fit,
    adjusted = adjusted_fit,
    study_results = study_results,
    metadata = md,
    adjusted_metadata = md_adjusted
  ),
  file.path(
    output_dir,
    "community_validation_4_corrected_beta_diversity.rds"
  )
)
