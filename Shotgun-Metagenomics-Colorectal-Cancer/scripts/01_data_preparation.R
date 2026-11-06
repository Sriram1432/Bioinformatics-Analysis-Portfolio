# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "readxl"))
suppressPackageStartupMessages(library(data.table))

# Clinical sample eligibility and study population
suppressPackageStartupMessages(library(readxl))

project_dir <- project_root()

metadata_file <- file.path(
  project_dir,
  "data",
  "metadata_extracted",
  "metadata_2340_CRC_cohort_20240704.xlsx"
)

output_dir <- file.path(
  project_dir,
  "results",
  "cohort"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
metadata_archive <- file.path(project_root(), "data", "metadata_2340_CRC_cohort_20240704.tar.gz")
if (!file.exists(metadata_file)) {
  if (!file.exists(metadata_archive)) stop("Metadata workbook or archive is missing")
  dir.create(dirname(metadata_file), recursive = TRUE, showWarnings = FALSE)
  utils::untar(metadata_archive, exdir = dirname(metadata_file))
}
dt <- as.data.table(
  read_excel(
    metadata_file,
    sheet = "metadata_2340_CRC_cohort"
  )
)


normalise_flag <- function(x) {
  tolower(trimws(as.character(x)))
}

dt[, contamination_flag :=
     normalise_flag(excluded) == "yes"]

dt[, metadata_exclusion_flag :=
     normalise_flag(to_exclude) == "yes"]

dt[is.na(contamination_flag), contamination_flag := FALSE]
dt[is.na(metadata_exclusion_flag), metadata_exclusion_flag := FALSE]


flag_overlap <- dt[
  ,
  .N,
  by = .(
    contamination_flag,
    metadata_exclusion_flag
  )
][order(contamination_flag, metadata_exclusion_flag)]

exclusion_by_timepoint <- dt[
  ,
  .N,
  by = .(
    timepoint,
    contamination_flag,
    metadata_exclusion_flag
  )
][order(timepoint)]


dt[, exclude_primary :=
     contamination_flag |
     metadata_exclusion_flag |
     timepoint != "t1"]

dt[is.na(exclude_primary), exclude_primary := TRUE]

eligible <- copy(
  dt[
    exclude_primary == FALSE &
      class %in% c("healthy", "adenoma", "CRC")
  ]
)

setorder(
  eligible,
  study_accession,
  individual_id,
  sample
)

eligible <- unique(
  eligible,
  by = c("study_accession", "individual_id")
)


eligibility_summary <- data.table(
  Metric = c(
    "Original samples",
    "Explicit contamination flags",
    "Metadata exclusion flags",
    "Follow-up samples (t2)",
    "Total excluded (unique samples)",
    "Eligible independent samples",
    "Healthy controls",
    "Adenoma",
    "CRC",
    "Studies represented",
    "Countries represented"
  ),
  
  Value = c(
    nrow(dt),
    sum(dt$contamination_flag),
    sum(dt$metadata_exclusion_flag),
    sum(dt$timepoint == "t2"),
    sum(
      dt$exclude_primary |
        !(dt$class %in% c("healthy", "adenoma", "CRC"))
    ),
    nrow(eligible),
    sum(eligible$class == "healthy"),
    sum(eligible$class == "adenoma"),
    sum(eligible$class == "CRC"),
    uniqueN(eligible$study_accession),
    uniqueN(eligible$country)
  )
)


study_summary <- dcast(
  eligible[
    ,
    .N,
    by = .(study_accession, class)
  ],
  study_accession ~ class,
  value.var = "N",
  fill = 0
)


independence_summary <- data.table(
  Metric = c(
    "Unique sample identifiers",
    "Unique study-participant combinations",
    "Duplicated study-participant combinations",
    "Non-baseline observations"
  ),
  
  Value = c(
    uniqueN(eligible$sample),
    uniqueN(
      eligible,
      by = c("study_accession", "individual_id")
    ),
    sum(
      duplicated(
        eligible,
        by = c("study_accession", "individual_id")
      )
    ),
    sum(eligible$timepoint != "t1")
  )
)


stopifnot(
  nrow(eligible) > 0,
  uniqueN(eligible$sample) == nrow(eligible),
  all(eligible$timepoint == "t1")
)

fwrite(
  eligible,
  file.path(
    output_dir,
    "metadata_eligible_primary.csv"
  )
)

fwrite(
  study_summary,
  file.path(
    output_dir,
    "eligible_clinical_groups_by_study.csv"
  )
)
