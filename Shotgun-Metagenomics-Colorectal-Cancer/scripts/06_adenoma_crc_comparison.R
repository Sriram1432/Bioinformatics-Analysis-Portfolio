# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table"))
suppressPackageStartupMessages(library(data.table))

# HEALTHY vs ADENOMA vs CRC
#
# Targeted analysis of previously identified CRC-associated MSPs
#
# Primary adjustment:
#   Study + log10(mapped sequencing reads)
#
# Three pairwise contrasts:
#   Adenoma vs healthy
#   CRC vs healthy
#   CRC vs adenoma
#
# HC3 robust SE and BH correction

library(data.table)


project_dir <- project_root()

output_dir <- file.path(
  project_dir,
  "results",
  "clinical_comparisons"
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

species <- fread(file.path(
  project_dir,
  "data",
  "species_signal_2340_CRC_cohort_20240322.tab"
))

previous_results <- fread(file.path(
  project_dir,
  "results",
  "microbial_associations",
  "CRC_species_discovery_all.csv"
))

stopifnot(
  !anyDuplicated(metadata$sample),
  !anyDuplicated(species$msp_id),
  !anyDuplicated(previous_results$Species_ID)
)


groups <- c("healthy", "adenoma", "CRC")

study_counts <- dcast(
  metadata[
    `class` %in% groups,
    .N,
    by = .(study_accession, `class`)
  ],
  study_accession ~ class,
  value.var = "N",
  fill = 0
)

three_group_studies <- study_counts[
  healthy > 0 & adenoma > 0 & CRC > 0,
  study_accession
]

stopifnot(
  length(three_group_studies) == 5
)

dt <- copy(
  metadata[
    study_accession %in% three_group_studies &
      `class` %in% groups
  ]
)

dt[, Disease := factor(
  `class`,
  levels = c("healthy", "adenoma", "CRC")
)]

dt[, Study := factor(study_accession)]

dt[, Mapped_reads := suppressWarnings(
  as.numeric(gut_mapped_read_count)
)]

dt[, Log_mapped_reads := fifelse(
  is.finite(Mapped_reads) & Mapped_reads > 0,
  log10(Mapped_reads),
  NA_real_
)]


required <- c(
  "Disease",
  "Study",
  "Log_mapped_reads"
)

complete <- complete.cases(dt[, ..required])

dt <- dt[complete]

dt[, Study := droplevels(Study)]
dt[, Disease := droplevels(Disease)]

stopifnot(
  nlevels(dt$Disease) == 3,
  nlevels(dt$Study) == 5,
  !anyDuplicated(dt$sample)
)


candidate_ids <- head(
  previous_results[
    order(FDR, P_value),
    Species_ID
  ],
  20
)

candidate_annotations <- previous_results[
  match(candidate_ids, Species_ID),
  .(
    Species_ID,
    Scientific_name,
    Discovery_CLR = CLR_difference,
    Discovery_FDR = FDR
  )
]

stopifnot(
  length(candidate_ids) == 20,
  !anyDuplicated(candidate_ids),
  all(candidate_ids %in% species$msp_id),
  all(dt$sample %in% names(species))
)


sample_ids <- dt$sample

X <- as.matrix(species[, ..sample_ids])

storage.mode(X) <- "double"

rownames(X) <- species$msp_id

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0),
  identical(colnames(X), dt$sample)
)

#
# Important:
# Use all 2,543 MSP features to define the CLR reference.
# Then extract the 20 prespecified candidate features.

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
  
  delta <- min(p[p > 0]) * 0.5
  
  delta <- min(delta, 0.5 / z)
  
  mass <- z * delta
  
  p[zeros] <- delta
  
  p[!zeros] <- p[!zeros] * (1 - mass)
  
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

rownames(CLR_all) <- species$msp_id
colnames(CLR_all) <- dt$sample

CLR <- CLR_all[
  candidate_ids,
  ,
  drop = FALSE
]

stopifnot(
  identical(rownames(CLR), candidate_ids),
  identical(colnames(CLR), dt$sample)
)


design <- model.matrix(
  ~ Disease + Study + Log_mapped_reads,
  data = dt
)

stopifnot(
  qr(design)$rank == ncol(design)
)

Y <- t(CLR)

fit <- lm.fit(
  x = design,
  y = Y
)

stopifnot(
  fit$rank == ncol(design)
)


XtX_inverse <- solve(
  crossprod(design)
)

hat_values <- rowSums(
  (design %*% XtX_inverse) * design
)

stopifnot(
  all(hat_values < 1 - 1e-10)
)

residuals <- as.matrix(fit$residuals)

# Contrast vectors
columns <- colnames(design)

contrast_adenoma <- setNames(
  rep(0, length(columns)),
  columns
)

contrast_crc <- contrast_adenoma
contrast_crc_adenoma <- contrast_adenoma

contrast_adenoma["Diseaseadenoma"] <- 1

contrast_crc["DiseaseCRC"] <- 1

contrast_crc_adenoma["DiseaseCRC"] <- 1
contrast_crc_adenoma["Diseaseadenoma"] <- -1

contrasts <- list(
  "Adenoma vs healthy" = contrast_adenoma,
  "CRC vs healthy" = contrast_crc,
  "CRC vs adenoma" = contrast_crc_adenoma
)


get_contrast <- function(name, contrast) {
  
  beta <- as.numeric(
    crossprod(
      contrast,
      fit$coefficients
    )
  )
  
  influence_weights <- as.vector(
    design %*% XtX_inverse %*% contrast
  )
  
  adjusted_weights <- influence_weights /
    (1 - hat_values)
  
  robust_variance <- colSums(
    sweep(
      residuals,
      1,
      adjusted_weights,
      "*"
    )^2
  )
  
  se <- sqrt(robust_variance)
  
  z <- beta / se
  
  p <- 2 * pnorm(-abs(z))
  
  data.table(
    Species_ID = candidate_ids,
    Comparison = name,
    CLR_difference = beta,
    Robust_SE = se,
    CI_lower = beta - 1.96 * se,
    CI_upper = beta + 1.96 * se,
    P_value = p
  )
}

association_results <- rbindlist(
  Map(
    get_contrast,
    names(contrasts),
    contrasts
  )
)

#
# Correct across all 60 prespecified contrasts.

association_results[
  ,
  FDR := p.adjust(
    P_value,
    method = "BH"
  )
]

association_results <- merge(
  association_results,
  candidate_annotations,
  by = "Species_ID",
  all.x = TRUE,
  sort = FALSE
)

association_results[
  ,
  Direction := fifelse(
    CLR_difference > 0,
    "Higher in first group",
    "Lower in first group"
  )
]


wide <- dcast(
  association_results,
  Species_ID + Scientific_name ~ Comparison,
  value.var = "CLR_difference"
)

wide[
  ,
  Association_pattern := fcase(
    `Adenoma vs healthy` > 0 &
      `CRC vs healthy` > 0,
    "Higher in adenoma and CRC",
    
    `Adenoma vs healthy` < 0 &
      `CRC vs healthy` < 0,
    "Lower in adenoma and CRC",
    
    default = "Opposing directions"
  )
]


detected <- X[
  candidate_ids,
  ,
  drop = FALSE
] > 0

detection_summary <- rbindlist(
  lapply(candidate_ids, function(id) {
    
    data.table(
      Species_ID = id,
      
      Healthy_percent = 100 * mean(
        detected[id, dt$Disease == "healthy"]
      ),
      
      Adenoma_percent = 100 * mean(
        detected[id, dt$Disease == "adenoma"]
      ),
      
      CRC_percent = 100 * mean(
        detected[id, dt$Disease == "CRC"]
      )
    )
    
  })
)

detection_summary <- merge(
  detection_summary,
  candidate_annotations[
    ,
    .(Species_ID, Scientific_name)
  ],
  by = "Species_ID",
  all.x = TRUE
)


fwrite(
  association_results,
  file.path(
    output_dir,
    "three_group_species_associations.csv"
  )
)

fwrite(
  wide,
  file.path(
    output_dir,
    "three_group_species_patterns.csv"
  )
)

fwrite(
  detection_summary,
  file.path(
    output_dir,
    "three_group_species_detection.csv"
  )
)

# HEALTHY vs ADENOMA vs CRC
# TARGETED FUNCTIONAL MODULE ANALYSIS

library(data.table)

project_dir <- project_root()

output_dir <- file.path(
  project_dir, "results", "clinical_comparisons"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


stopifnot(
  exists("dt"),
  all(c(
    "sample", "Disease", "Study",
    "Log_mapped_reads"
  ) %in% names(dt)),
  !anyDuplicated(dt$sample),
  nlevels(dt$Disease) == 3,
  uniqueN(dt$study_accession) == 5
)


module_ids <- c(
  "MGB017",  # IPA synthesis
  "MF0009",  # Tryptophan degradation
  "MGB005",  # Tryptophan synthesis
  "MF0121",  # Propionate, acrylate pathway
  "MGB055",  # Propionate synthesis III
  "MF0116",  # Butyrate, transferase pathway
  "MGB053",  # Butyrate synthesis II
  "MF0001",  # Ethanol production
  "MF0058",  # Fructose degradation
  "M00017",  # Methionine biosynthesis
  "M00632",  # Galactose degradation
  "M00028"   # Ornithine biosynthesis
)


functional <- fread(file.path(
  project_dir,
  "data",
  "functional_signal_2340_CRC_cohort_20240329.tab"
))

definitions <- fread(file.path(
  project_dir,
  "data",
  "modules_definition_GMM_GBM_KEGG_92_20240329.tab"
))

stopifnot(
  all(module_ids %in% functional$modules_id),
  all(dt$sample %in% names(functional)),
  !anyDuplicated(functional$modules_id)
)

sample_ids <- as.character(dt$sample)
sample_ids <- as.character(dt$sample)
sample_ids <- as.character(dt$sample)
X <- as.matrix(
  functional[, ..sample_ids]
)

storage.mode(X) <- "double"

rownames(X) <- functional$modules_id

stopifnot(
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0)
)


P <- sweep(X, 2, colSums(X), "/")

replace_zeros <- function(p) {
  
  zero <- p == 0
  n_zero <- sum(zero)
  
  if (!n_zero) return(p)
  
  delta <- min(
    min(p[!zero]) * 0.5,
    0.5 / n_zero
  )
  
  p[zero] <- delta
  
  p[!zero] <- p[!zero] *
    (1 - n_zero * delta)
  
  p
}

P <- apply(P, 2, replace_zeros)

stopifnot(
  all(P > 0),
  max(abs(colSums(P) - 1)) < 1e-8
)

L <- log(P)

CLR_all <- sweep(
  L,
  2,
  colMeans(L),
  "-"
)

rownames(CLR_all) <- functional$modules_id
colnames(CLR_all) <- dt$sample

CLR_modules <- CLR_all[
  module_ids,
  ,
  drop = FALSE
]


design <- model.matrix(
  ~ Disease + Study + Log_mapped_reads,
  data = dt
)

stopifnot(
  qr(design)$rank == ncol(design),
  identical(colnames(CLR_modules), dt$sample)
)

Y <- t(CLR_modules)

fit <- lm.fit(design, Y)

XtX_inverse <- chol2inv(
  chol(crossprod(design))
)

h <- rowSums(
  (design %*% XtX_inverse) * design
)

stopifnot(
  all(h < 1 - 1e-10)
)

residual_matrix <- as.matrix(fit$residuals)


make_contrast <- function(
    adenoma = 0,
    crc = 0
) {
  
  v <- setNames(
    numeric(ncol(design)),
    colnames(design)
  )
  
  v["Diseaseadenoma"] <- adenoma
  v["DiseaseCRC"] <- crc
  
  v
}

contrasts <- list(
  "Adenoma vs healthy" = make_contrast(
    adenoma = 1
  ),
  
  "CRC vs healthy" = make_contrast(
    crc = 1
  ),
  
  "CRC vs adenoma" = make_contrast(
    adenoma = -1,
    crc = 1
  )
)

estimate_contrast <- function(label, contrast) {
  
  estimate <- as.numeric(
    crossprod(contrast, fit$coefficients)
  )
  
  weights <- as.vector(
    design %*% XtX_inverse %*% contrast
  ) / (1 - h)
  
  se <- sqrt(
    colSums(
      sweep(
        residual_matrix,
        1,
        weights,
        "*"
      )^2
    )
  )
  
  stopifnot(
    all(is.finite(se)),
    all(se > 0)
  )
  
  p <- 2 * pnorm(
    -abs(estimate / se)
  )
  
  data.table(
    Module_ID = module_ids,
    Comparison = label,
    CLR_difference = estimate,
    Robust_SE = se,
    CI_lower = estimate - 1.96 * se,
    CI_upper = estimate + 1.96 * se,
    P_value = p
  )
}

module_associations <- rbindlist(
  Map(
    estimate_contrast,
    names(contrasts),
    contrasts
  )
)

# Correct across all 36 comparisons
module_associations[
  ,
  FDR := p.adjust(
    P_value,
    method = "BH"
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

module_associations <- merge(
  module_associations,
  annotations,
  by = "Module_ID",
  all.x = TRUE,
  sort = FALSE
)


descriptive <- rbindlist(
  lapply(module_ids, function(id) {
    
    data.table(
      Module_ID = id,
      Group = as.character(dt$Disease),
      CLR = as.numeric(CLR_modules[id, ])
    )[
      ,
      .(
        Mean_CLR = mean(CLR),
        Median_CLR = median(CLR),
        N = .N
      ),
      by = .(Module_ID, Group)
    ]
  })
)


fwrite(
  module_associations,
  file.path(
    output_dir,
    "three_group_functional_associations.csv"
  )
)

fwrite(
  descriptive,
  file.path(
    output_dir,
    "three_group_functional_descriptive.csv"
  )
)
