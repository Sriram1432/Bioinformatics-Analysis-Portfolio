# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table"))
suppressPackageStartupMessages(library(data.table))

# Microbial community landscape

library(data.table)


project_dir <- project_root()

species_file <- file.path(
  project_dir,
  "data",
  "species_signal_2340_CRC_cohort_20240322.tab"
)

metadata_file <- file.path(
  project_dir,
  "results",
  "cohort",
  "metadata_eligible_primary.csv"
)

output_dir <- file.path(
  project_dir,
  "results",
  "profiling"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


metadata <- fread(metadata_file)

species <- fread(species_file)


sample_ids <- metadata$sample

stopifnot(
  length(sample_ids) == 2138,
  !anyDuplicated(sample_ids),
  all(sample_ids %in% names(species))
)

# Create features x samples matrix
X <- as.matrix(
  species[
    ,
    ..sample_ids
  ]
)

rownames(X) <- species$msp_id

storage.mode(X) <- "double"

stopifnot(
  !anyNA(X),
  all(is.finite(X)),
  all(X >= 0),
  all(colSums(X) > 0)
)


presence <- X > 0

species_prevalence <- rowMeans(presence)


sample_richness <- colSums(presence)


sample_totals <- colSums(X)

P <- sweep(
  X,
  2,
  sample_totals,
  FUN = "/"
)


mean_proportion <- rowMeans(P)


species_sd <- apply(P, 1, sd)

species_cv <- species_sd / pmax(
  mean_proportion,
  .Machine$double.eps
)


species_results <- data.table(
  Species_ID = rownames(X),
  Prevalence_percent = 100 * species_prevalence,
  Mean_proportion_percent = 100 * mean_proportion,
  SD_proportion_percent = 100 * species_sd,
  CV = species_cv
)

setorder(
  species_results,
  -Prevalence_percent
)


sample_results <- data.table(
  sample = colnames(X),
  Richness = as.integer(sample_richness),
  Original_signal_total = sample_totals
)

sample_results <- merge(
  sample_results,
  metadata[
    ,
    .(
      sample,
      study_accession,
      country,
      class
    )
  ],
  by = "sample",
  sort = FALSE
)


richness_summary <- data.table(
  Metric = c(
    "Samples",
    "Species features",
    "Median species richness",
    "Mean species richness",
    "Minimum richness",
    "25th percentile richness",
    "75th percentile richness",
    "Maximum richness"
  ),
  
  Value = c(
    ncol(X),
    nrow(X),
    median(sample_richness),
    mean(sample_richness),
    min(sample_richness),
    as.numeric(quantile(sample_richness, 0.25)),
    as.numeric(quantile(sample_richness, 0.75)),
    max(sample_richness)
  )
)

prevalence_summary <- data.table(
  Prevalence_threshold = c(
    "At least 1%",
    "At least 5%",
    "At least 10%",
    "At least 25%",
    "At least 50%",
    "At least 75%",
    "At least 90%"
  ),
  
  Species_count = c(
    sum(species_prevalence >= 0.01),
    sum(species_prevalence >= 0.05),
    sum(species_prevalence >= 0.10),
    sum(species_prevalence >= 0.25),
    sum(species_prevalence >= 0.50),
    sum(species_prevalence >= 0.75),
    sum(species_prevalence >= 0.90)
  )
)


top_abundant <- copy(species_results)

setorder(
  top_abundant,
  -Mean_proportion_percent
)

top_abundant <- head(top_abundant, 20)

top_prevalent <- head(species_results, 20)


cohort_richness <- sample_results[
  ,
  .(
    Samples = .N,
    Median_richness = median(Richness),
    Mean_richness = mean(Richness),
    SD_richness = sd(Richness)
  ),
  by = .(
    study_accession,
    country
  )
][order(-Samples)]


fwrite(
  species_results,
  file.path(output_dir, "species_landscape.csv")
)

fwrite(
  sample_results,
  file.path(output_dir, "sample_richness.csv")
)

fwrite(
  cohort_richness,
  file.path(output_dir, "cohort_richness.csv")
)

fwrite(
  richness_summary,
  file.path(output_dir, "richness_summary.csv")
)

fwrite(
  prevalence_summary,
  file.path(output_dir, "prevalence_summary.csv")
)

# Biological identification of prominent microorganisms

library(data.table)

project_dir <- project_root()

data_dir <- file.path(project_dir, "data")

output_dir <- file.path(
  project_dir,
  "results",
  "profiling"
)

taxonomy <- fread(
  file.path(
    data_dir,
    "MSP_set_ref_oss_gss_status_612_ref_paper_20250220_gtdb_220.tab"
  )
)

landscape <- fread(
  file.path(output_dir, "species_landscape.csv")
)

taxonomy_selected <- taxonomy[
  ,
  .(
    Species_ID = msp_id,
    Scientific_name = species,
    Genus = genus,
    Family = family,
    Phylum = phylum,
    Domain = superkingdom,
    Catalogue = catalogue,
    MSP_class = MSP_class
  )
]

annotated <- merge(
  landscape,
  taxonomy_selected,
  by = "Species_ID",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  nrow(annotated) == nrow(landscape),
  !anyDuplicated(annotated$Species_ID)
)

top_prevalent <- copy(annotated)

setorder(
  top_prevalent,
  -Prevalence_percent
)

top_prevalent <- head(top_prevalent, 20)

top_abundant <- copy(annotated)

setorder(
  top_abundant,
  -Mean_proportion_percent
)

top_abundant <- head(top_abundant, 20)

annotation_summary <- data.table(
  Metric = c(
    "Total MSP features",
    "Matched taxonomy records",
    "Features with species names",
    "Features with genus names",
    "Features without species names",
    "Unique species names",
    "Unique genera"
  ),
  
  Value = c(
    nrow(annotated),
    sum(!is.na(annotated$Domain)),
    sum(!is.na(annotated$Scientific_name)),
    sum(!is.na(annotated$Genus)),
    sum(is.na(annotated$Scientific_name)),
    uniqueN(
      annotated$Scientific_name,
      na.rm = TRUE
    ),
    uniqueN(
      annotated$Genus,
      na.rm = TRUE
    )
  )
)

domain_summary <- annotated[
  ,
  .N,
  by = Domain
][order(-N)]

phylum_summary <- annotated[
  ,
  .N,
  by = .(Domain, Phylum)
][order(-N)]

fwrite(
  annotated,
  file.path(
    output_dir,
    "species_landscape_annotated.csv"
  )
)
