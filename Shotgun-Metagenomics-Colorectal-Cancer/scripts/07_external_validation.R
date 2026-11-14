# Colorectal cancer shotgun metagenomics
# Source profiles are harmonised signals; zero entries denote nondetection,
# not confirmed biological absence. All tests are observational.
source(file.path(Sys.getenv("PROJECT4_ROOT", unset = getwd()), "scripts", "utils.R"), local = TRUE)
require_packages(c("data.table", "vegan", "permute", "brglm2", "metafor", "sandwich"))
suppressPackageStartupMessages(library(data.table))

# INDEPENDENT VALIDATION DATASET INSPECTION
#
# Source: Piccinno et al., 2025
# DOI: 10.5281/zenodo.15069069
#
library(data.table)

project_dir <- project_root()

data_dir <- file.path(
  project_dir,
  "data",
  "validation"
)

output_dir <- file.path(
  project_dir,
  "results",
  "taxonomic_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


archives <- c(
  "metadata.tar.gz",
  "metaphlan4_profiles.tar.gz",
  "humann3.6_pathways.tar.gz"
)

archive_paths <- file.path(data_dir, archives)

if (!all(file.exists(archive_paths))) {
  stop(
    "Missing files: ",
    paste(
      archives[!file.exists(archive_paths)],
      collapse = ", "
    )
  )
}

file_info <- data.table(
  Archive = archives,
  Size_MB = round(
    file.size(archive_paths) / 1024^2,
    2
  )
)


archive_contents <- lapply(
  archive_paths,
  function(path) {
    utils::untar(
      path,
      list = TRUE
    )
  }
)

names(archive_contents) <- archives

archive_summary <- rbindlist(
  lapply(seq_along(archives), function(i) {
    data.table(
      Archive = archives[i],
      Entries = length(archive_contents[[i]]),
      First_entries = paste(
        head(archive_contents[[i]], 3),
        collapse = " | "
      )
    )
  }),
  fill = TRUE
)


metadata_dir <- file.path(
  data_dir,
  "metadata_extracted"
)

dir.create(
  metadata_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

utils::untar(
  archive_paths[1],
  exdir = metadata_dir
)

metadata_files <- list.files(
  metadata_dir,
  recursive = TRUE,
  full.names = TRUE
)

metadata_files <- metadata_files[
  !dir.exists(metadata_files)
]

metadata_file_info <- data.table(
  File = basename(metadata_files),
  Path = metadata_files,
  Size_KB = round(
    file.size(metadata_files) / 1024,
    2
  )
)


tabular_files <- metadata_files[
  grepl(
    "\\.(tsv|tab|txt|csv)(\\.gz)?$",
    metadata_files,
    ignore.case = TRUE
  )
]

metadata_previews <- list()

for (f in tabular_files) {
  
  result <- tryCatch({
    
    dt <- fread(
      f,
      nrows = 5,
      showProgress = FALSE
    )
    
    list(
      file = basename(f),
      columns = names(dt),
      preview = dt
    )
    
  }, error = function(e) {
    
    list(
      file = basename(f),
      error = conditionMessage(e)
    )
    
  })
  
  metadata_previews[[length(metadata_previews) + 1]] <- result
}


fwrite(
  archive_summary,
  file.path(
    output_dir,
    "external_validation_archive_summary.csv"
  )
)

fwrite(
  metadata_file_info,
  file.path(
    output_dir,
    "external_validation_metadata_files.csv"
  )
)

# INDEPENDENT VALIDATION COHORT ELIGIBILITY

library(data.table)

project_dir <- project_root()

metadata_dir <- file.path(
  project_dir,
  "data",
  "validation",
  "metadata_extracted",
  "metadata"
)

output_dir <- file.path(project_dir, "results", "taxonomic_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


files <- list.files(
  metadata_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

all_metadata <- rbindlist(
  lapply(files, function(f) {
    dt <- fread(f)
    dt[, Cohort := tools::file_path_sans_ext(basename(f))]
    dt
  }),
  fill = TRUE
)


all_metadata[, Group := fcase(
  `Study condition` == "CRC", "CRC",
  `Study condition` == "Control", "healthy",
  `Study condition` == "Adenoma", "adenoma",
  default = "other"
)]

all_metadata[, Cohort_type := fifelse(
  grepl("^This_study_", Cohort),
  "New cohort",
  "Previously published cohort"
)]


counts <- dcast(
  all_metadata[
    ,
    .N,
    by = .(Cohort_type, Cohort, Group)
  ],
  Cohort_type + Cohort ~ Group,
  value.var = "N",
  fill = 0
)

for (column in c("CRC", "healthy", "adenoma", "other")) {
  if (!column %in% names(counts)) {
    counts[, (column) := 0L]
  }
}

counts[, Total := CRC + healthy + adenoma + other]

counts[, CRC_validation_eligible :=
         CRC > 0 & healthy > 0
]

setorder(counts, Cohort_type, Cohort)


new_cohorts <- counts[
  Cohort_type == "New cohort"
]

eligible_new <- new_cohorts[
  CRC_validation_eligible == TRUE
]


duplicate_ids <- all_metadata[
  ,
  .(
    Duplicate_IDs = sum(duplicated(Name)),
    Missing_IDs = sum(is.na(Name) | !nzchar(Name))
  ),
  by = Cohort
]


fwrite(
  counts,
  file.path(output_dir, "validation_all_cohort_counts.csv")
)

fwrite(
  all_metadata[
    Cohort_type == "New cohort"
  ],
  file.path(output_dir, "validation_new_cohort_metadata.csv")
)

# EXTERNAL VALIDATION — TAXONOMIC DATA PREPARATION
#
library(data.table)

project_dir <- project_root()

data_dir <- file.path(project_dir, "data", "validation")

output_dir <- file.path(project_dir, "results", "taxonomic_validation")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


eligible_cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)

metadata_validation <- copy(
  all_metadata[Cohort %in% eligible_cohorts]
)


archive <- file.path(
  data_dir,
  "metaphlan4_profiles.tar.gz"
)

stopifnot(file.exists(archive))

archive_entries <- untar(archive, list = TRUE)

selected_entries <- archive_entries[
  basename(archive_entries) %in%
    paste0(eligible_cohorts, ".tsv")
]

if (length(selected_entries) != length(eligible_cohorts)) {
  stop(
    "Expected four taxonomic profile files; found ",
    length(selected_entries),
    ". Check archive filenames."
  )
}

extract_dir <- file.path(
  data_dir,
  "taxonomic_validation"
)

dir.create(
  extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

untar(
  archive,
  files = selected_entries,
  exdir = extract_dir
)


profile_inspection <- list()

for (cohort in eligible_cohorts) {
  
  profile_file <- file.path(
    extract_dir,
    "metaphlan4_profiles",
    paste0(cohort, ".tsv")
  )
  
  stopifnot(file.exists(profile_file))
  
  profile <- fread(
    profile_file,
    showProgress = FALSE
  )
  
  cohort_metadata <- metadata_validation[
    Cohort == cohort
  ]
  
  # Candidate sample identifiers are metadata Name values.
  # Check whether they occur in the profile column names.
  
  sample_matches <- intersect(
    cohort_metadata$Name,
    names(profile)
  )
  
  profile_inspection[[cohort]] <- list(
    Cohort = cohort,
    Profile_rows = nrow(profile),
    Profile_columns = ncol(profile),
    Metadata_samples = nrow(cohort_metadata),
    Matching_sample_columns = length(sample_matches),
    Column_names = head(names(profile), 8),
    Example_rows = head(profile, 3),
    Sample_matches = head(sample_matches, 5)
  )
}


profile_summary <- rbindlist(
  lapply(profile_inspection, function(x) {
    
    data.table(
      Cohort = x$Cohort,
      Profile_rows = x$Profile_rows,
      Profile_columns = x$Profile_columns,
      Metadata_samples = x$Metadata_samples,
      Matching_sample_columns =
        x$Matching_sample_columns
    )
    
  })
)


fwrite(
  profile_summary,
  file.path(
    output_dir,
    "validation_taxonomic_profile_inspection.csv"
  )
)

# EXTERNAL VALIDATION — SPECIES COMPATIBILITY
#
# Checks:
#
# No statistical testing yet

library(data.table)

project_dir <- project_root()

profile_dir <- file.path(
  project_dir,
  "data",
  "validation",
  "taxonomic_validation",
  "metaphlan4_profiles"
)

output_dir <- file.path(
  project_dir,
  "results",
  "taxonomic_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


target_species <- c(
  "Parvimonas micra",
  "Peptostreptococcus stomatis",
  "Allisonella pneumosinta",
  "Gemella morbillorum",
  "Porphyromonas asaccharolytica",
  "Bulleidia moorei",
  "Fusobacterium animalis",
  "Roseburia intestinalis",
  "Agathobacter faecis",
  "Hominilimicola fabiformis",
  "Streptococcus salivarius",
  "Butyribacter hominis"
)

expected_direction <- c(
  rep("Higher in CRC", 7),
  rep("Lower in CRC", 5)
)

targets <- data.table(
  Species = target_species,
  Expected_direction = expected_direction
)

cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)


extract_species <- function(x) {
  
  # Match the terminal species-level identifier.
  # Excludes genus, family, strain and other ranks.
  
  is_species <- grepl(
    "(^|\\|)s__[^|]+$",
    x
  )
  
  result <- sub(
    "^.*s__",
    "",
    x[is_species]
  )
  
  list(
    indices = which(is_species),
    species = gsub("_", " ", result, fixed = TRUE)
  )
}


coverage_results <- list()
alternative_results <- list()
species_counts <- list()

for (cohort in cohorts) {
  
  file <- file.path(
    profile_dir,
    paste0(cohort, ".tsv")
  )
  
  stopifnot(file.exists(file))
  
  dt <- fread(file, showProgress = FALSE)
  
  stopifnot(
    "clade_name" %in% names(dt)
  )
  
  parsed <- extract_species(dt$clade_name)
  
  species_dt <- dt[parsed$indices]
  
  species_dt[, Species := parsed$species]
  
  # All other columns are sample abundances.
  sample_columns <- setdiff(
    names(species_dt),
    c("clade_name", "Species")
  )
  
  abundance <- as.matrix(
    species_dt[, ..sample_columns]
  )
  
  storage.mode(abundance) <- "double"
  
  detected <- rowSums(
    abundance > 0,
    na.rm = TRUE
  )
  
  prevalence <- detected / ncol(abundance)
  
  species_info <- data.table(
    Species = species_dt$Species,
    Clade = species_dt$clade_name,
    Detected_samples = detected,
    Prevalence = prevalence
  )
  
  # Exact candidate matches
  matched <- merge(
    targets,
    species_info,
    by = "Species",
    all.x = TRUE,
    sort = FALSE
  )
  
  matched[, Cohort := cohort]
  
  matched[, Exact_match := !is.na(Clade)]
  
  coverage_results[[cohort]] <- matched
  
  # Potential alternative genus/species classifications
  search_terms <- c(
    "Parvimonas",
    "Peptostreptococcus",
    "Allisonella",
    "Gemella",
    "Porphyromonas",
    "Bulleidia",
    "Fusobacterium",
    "Roseburia",
    "Agathobacter",
    "Hominilimicola",
    "Streptococcus",
    "Butyribacter"
  )
  
  alternatives <- species_info[
    grepl(
      paste(search_terms, collapse = "|"),
      Species,
      ignore.case = TRUE
    )
  ]
  
  alternatives[, Cohort := cohort]
  
  alternative_results[[cohort]] <- alternatives
  
  species_counts[[cohort]] <- data.table(
    Cohort = cohort,
    Total_taxonomic_rows = nrow(dt),
    Species_rows = nrow(species_dt),
    Samples = length(sample_columns)
  )
}


coverage <- rbindlist(
  coverage_results,
  fill = TRUE
)

alternatives <- rbindlist(
  alternative_results,
  fill = TRUE
)

species_summary <- rbindlist(species_counts)

# Compact cross-cohort summary
coverage_summary <- coverage[
  ,
  .(
    Cohorts_with_exact_match = sum(Exact_match),
    Mean_prevalence = if (
      all(is.na(Prevalence))
    ) NA_real_ else mean(Prevalence, na.rm = TRUE),
    Expected_direction = Expected_direction[1]
  ),
  by = Species
]


fwrite(
  coverage,
  file.path(
    output_dir,
    "validation_species_coverage_by_cohort.csv"
  )
)

fwrite(
  coverage_summary,
  file.path(
    output_dir,
    "validation_species_coverage_summary.csv"
  )
)

fwrite(
  alternatives,
  file.path(
    output_dir,
    "validation_possible_taxonomic_alternatives.csv"
  )
)

# INDEPENDENT VALIDATION OF CRC MICROBIAL SIGNATURES
#
# Primary question:
# Are discovery-associated microorganisms differentially
# detected in independent CRC versus healthy samples?
#
# Analysis:

library(data.table)

project_dir <- project_root()

profile_dir <- file.path(
  project_dir,
  "data",
  "validation",
  "taxonomic_validation",
  "metaphlan4_profiles"
)

output_dir <- file.path(
  project_dir,
  "results",
  "taxonomic_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


targets <- data.table(
  Species = c(
    "Parvimonas micra",
    "Peptostreptococcus stomatis",
    "Gemella morbillorum",
    "Porphyromonas asaccharolytica",
    "Roseburia intestinalis",
    "Streptococcus salivarius"
  ),
  Expected = c(
    "Higher",
    "Higher",
    "Higher",
    "Higher",
    "Lower",
    "Lower"
  )
)

cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)


validation_records <- list()

for (cohort in cohorts) {
  
  path <- file.path(
    profile_dir,
    paste0(cohort, ".tsv")
  )
  
  dt <- fread(path, showProgress = FALSE)
  
  meta <- copy(
    all_metadata[
      Cohort == cohort &
        `Study condition` %in% c("CRC", "Control")
    ]
  )
  
  stopifnot(
    !anyDuplicated(meta$Name),
    all(meta$Name %in% names(dt))
  )
  
  # Use terminal species-level taxonomic entries only.
  
  species_rows <- grepl(
    "(^|\\|)s__[^|]+$",
    dt$clade_name
  )
  
  species_dt <- dt[species_rows]
  
  species_dt[, Species := gsub(
    "_",
    " ",
    sub("^.*s__", "", clade_name),
    fixed = TRUE
  )]
  
  # Require unique species identification.
  
  for (species in targets$Species) {
    
    hit <- species_dt[
      Species == species
    ]
    
    if (nrow(hit) != 1L) {
      stop(
        "Expected one row for ",
        species,
        " in ",
        cohort,
        "; found ",
        nrow(hit)
      )
    }
    
    sample_ids <- meta$Name
    
    abundance <- as.numeric(
      unlist(
        hit[, ..sample_ids],
        use.names = FALSE
      )
    )
    
    if (
      length(abundance) != nrow(meta) ||
      anyNA(abundance) ||
      any(abundance < 0)
    ) {
      stop("Invalid abundance values: ", cohort, " / ", species)
    }
    
    validation_records[[length(validation_records) + 1L]] <-
      data.table(
        Cohort = cohort,
        Sample = sample_ids,
        Group = ifelse(
          meta[["Study condition"]] == "CRC",
          "CRC",
          "healthy"
        ),
        Species = species,
        Abundance_percent = abundance,
        Detected = abundance > 0
      )
  }
}

validation_long <- rbindlist(validation_records)


study_results <- list()
pooled_results <- list()

for (species in targets$Species) {
  
  species_data <- validation_long[
    Species == species
  ]
  
  tables <- list()
  
  for (cohort in cohorts) {
    
    sub <- species_data[
      Cohort == cohort
    ]
    
    # Rows: CRC / healthy
    # Columns: detected / not detected
    
    tab <- table(
      factor(
        sub$Group,
        levels = c("CRC", "healthy")
      ),
      factor(
        sub$Detected,
        levels = c(TRUE, FALSE)
      )
    )
    
    tables[[cohort]] <- tab
    
    crc_yes <- tab[1, 1]
    crc_no <- tab[1, 2]
    
    healthy_yes <- tab[2, 1]
    healthy_no <- tab[2, 2]
    
    # Haldane-Anscombe correction for finite
    # descriptive odds ratios in sparse tables.
    
    odds_ratio <- (
      (crc_yes + 0.5) *
        (healthy_no + 0.5)
    ) / (
      (crc_no + 0.5) *
        (healthy_yes + 0.5)
    )
    
    fisher_p <- fisher.test(tab)$p.value
    
    study_results[[length(study_results) + 1L]] <-
      data.table(
        Species = species,
        Cohort = cohort,
        CRC_N = sum(tab[1, ]),
        Healthy_N = sum(tab[2, ]),
        CRC_detected = crc_yes,
        Healthy_detected = healthy_yes,
        CRC_prevalence = crc_yes / sum(tab[1, ]),
        Healthy_prevalence = healthy_yes / sum(tab[2, ]),
        OR_corrected = odds_ratio,
        Fisher_p = fisher_p
      )
  }
  
  
  tab_array <- array(
    0,
    dim = c(2, 2, length(cohorts)),
    dimnames = list(
      Group = c("CRC", "healthy"),
      Detection = c("Detected", "Not_detected"),
      Cohort = cohorts
    )
  )
  
  for (i in seq_along(cohorts)) {
    tab_array[, , i] <- tables[[cohorts[i]]]
  }
  
  mh <- mantelhaen.test(
    tab_array,
    correct = FALSE
  )
  
  expected <- targets[
    Species == species,
    Expected
  ][1]
  
  or_estimate <- unname(
    mh$estimate
  )
  
  observed_direction <- if (
    is.finite(or_estimate) &&
    or_estimate > 1
  ) {
    "Higher"
  } else if (
    is.finite(or_estimate) &&
    or_estimate < 1
  ) {
    "Lower"
  } else {
    "Undetermined"
  }
  
  pooled_results[[length(pooled_results) + 1L]] <-
    data.table(
      Species = species,
      Expected = expected,
      MH_OR = or_estimate,
      CI_lower = mh$conf.int[1],
      CI_upper = mh$conf.int[2],
      P_value = mh$p.value,
      Observed = observed_direction,
      Direction_agreement =
        observed_direction == expected
    )
}

study_results <- rbindlist(study_results)

pooled_results <- rbindlist(pooled_results)

pooled_results[
  ,
  FDR := p.adjust(
    P_value,
    method = "BH"
  )
]

setorder(
  pooled_results,
  FDR,
  Species
)


fwrite(
  validation_long,
  file.path(
    output_dir,
    "taxonomic_validation_8_validation_detection_data.csv"
  )
)

fwrite(
  study_results,
  file.path(
    output_dir,
    "taxonomic_validation_8_cohort_detection_results.csv"
  )
)

fwrite(
  pooled_results,
  file.path(
    output_dir,
    "taxonomic_validation_8_MH_validation_results.csv"
  )
)

# INDEPENDENT EXTERNAL VALIDATION
#
# A. Covariate-adjusted detection associations
# B. Covariate-adjusted relative-abundance associations
# C. Random-effects cross-cohort meta-analysis
#
# Six prespecified species.

library(data.table)

packages <- c(
  "brglm2",
  "metafor",
  "sandwich"
)

project_dir <- project_root()

output_dir <- file.path(project_dir, "results", "taxonomic_validation")

profile_dir <- file.path(
  project_dir,
  "data",
  "validation",
  "taxonomic_validation",
  "metaphlan4_profiles"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

validation_file <- file.path(
  output_dir,
  "taxonomic_validation_8_validation_detection_data.csv"
)

stopifnot(
  file.exists(validation_file),
  exists("all_metadata")
)

vdat <- fread(validation_file)


target_species <- data.table(
  Species = c(
    "Parvimonas micra",
    "Peptostreptococcus stomatis",
    "Gemella morbillorum",
    "Porphyromonas asaccharolytica",
    "Roseburia intestinalis",
    "Streptococcus salivarius"
  ),
  
  Expected = c(
    1, 1, 1, 1, -1, -1
  )
)

cohorts <- unique(vdat$Cohort)


meta <- copy(
  all_metadata[
    Cohort %in% cohorts &
      `Study condition` %in% c("CRC", "Control")
  ]
)

meta <- meta[
  ,
  .(
    Cohort,
    Sample = Name,
    Age = suppressWarnings(as.numeric(Age)),
    BMI = suppressWarnings(as.numeric(BMI)),
    Sex = as.character(Sex)
  )
]

meta[
  Sex %in% c("", "nd", "NA", "Unknown"),
  Sex := NA_character_
]

stopifnot(
  !anyDuplicated(meta[, .(Cohort, Sample)])
)

vdat <- merge(
  vdat,
  meta,
  by = c("Cohort", "Sample"),
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  !anyNA(vdat$Group),
  !anyNA(vdat$Abundance_percent)
)

vdat[
  ,
  CRC := as.integer(Group == "CRC")
]

vdat[
  ,
  Detected := as.integer(Abundance_percent > 0)
]

#
# Important:
# MetaPhlAn tables contain overlapping taxonomic ranks.
# Only terminal species rows are used.
#
# Zero values are replaced by a fixed small value
# after within-sample species-level closure.
#
# CLR is calculated across all reported species,
# not only the six validation candidates.

clr_records <- list()

for (cohort in cohorts) {
  
  path <- file.path(
    profile_dir,
    paste0(cohort, ".tsv")
  )
  
  dt <- fread(path)
  
  keep <- grepl(
    "(^|\\|)s__[^|]+$",
    dt$clade_name
  )
  
  sp <- dt[keep]
  
  sp[
    ,
    Species := gsub(
      "_",
      " ",
      sub("^.*s__", "", clade_name),
      fixed = TRUE
    )
  ]
  
  ids <- unique(
    vdat[Cohort == cohort, Sample]
  )
  
  stopifnot(
    all(ids %in% names(sp)),
    !anyDuplicated(sp$Species)
  )
  
  abundance_matrix <- as.matrix(
    sp[, ..ids]
  )
  
  storage.mode(abundance_matrix) <- "double"
  
  stopifnot(
    all(is.finite(abundance_matrix)),
    all(abundance_matrix >= 0)
  )
  
  # Profiles are supplied as percentages.
  # Close species-level abundances to proportions.
  
  sample_totals <- colSums(abundance_matrix)
  
  stopifnot(all(sample_totals > 0))
  
  rel <- sweep(
    abundance_matrix,
    2,
    sample_totals,
    "/"
  )
  
  # Explicit fixed zero replacement for comparison.
  # A sensitivity analysis can assess this choice later.
  
  log_rel <- log(rel + 1e-6)
  
  clr <- sweep(
    log_rel,
    2,
    colMeans(log_rel),
    "-"
  )
  
  rownames(clr) <- sp$Species
  
  for (species in target_species$Species) {
    
    stopifnot(species %in% rownames(clr))
    
    clr_records[[length(clr_records) + 1L]] <-
      data.table(
        Cohort = cohort,
        Sample = ids,
        Species = species,
        CLR = as.numeric(clr[species, ])
      )
  }
}

clr_data <- rbindlist(clr_records)

vdat <- merge(
  vdat,
  clr_data,
  by = c("Cohort", "Sample", "Species"),
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  !anyNA(vdat$CLR),
  nrow(vdat) == nrow(clr_data)
)

#
# Keep age, sex and BMI when available.
#
# If missing covariates would remove too many subjects,
# remove the covariate with the greatest missingness.
#
# Preserve at least 8 CRC and 8 healthy individuals
# for the initial adjusted cohort analysis.
#
# A covariate with no within-cohort variation is excluded.

prepare_model <- function(d) {
  
  covars <- c("Age", "Sex", "BMI")
  
  covars <- covars[
    vapply(
      covars,
      function(z) {
        x <- d[[z]]
        x <- x[!is.na(x)]
        length(unique(x)) > 1
      },
      logical(1)
    )
  ]
  
  repeat {
    
    variables <- c(
      "CRC", "Detected", "CLR", covars
    )
    
    cc <- complete.cases(
      d[, ..variables]
    )
    
    n_crc <- sum(d$CRC[cc] == 1)
    n_healthy <- sum(d$CRC[cc] == 0)
    
    if (
      (n_crc >= 8 && n_healthy >= 8) ||
      length(covars) == 0L
    ) {
      break
    }
    
    missing_fraction <- vapply(
      covars,
      function(z) mean(is.na(d[[z]])),
      numeric(1)
    )
    
    drop_var <- covars[
      which.max(missing_fraction)
    ]
    
    covars <- setdiff(covars, drop_var)
  }
  
  variables <- c(
    "CRC", "Detected", "CLR", covars
  )
  
  d <- d[
    complete.cases(d[, ..variables])
  ]
  
  # Recheck variation after complete-case selection.
  
  covars <- covars[
    vapply(
      covars,
      function(z) {
        length(unique(d[[z]])) > 1
      },
      logical(1)
    )
  ]
  
  if ("Sex" %in% covars) {
    d[, Sex := factor(Sex)]
  }
  
  list(
    data = d,
    covars = covars
  )
}


model_results <- list()
model_diagnostics <- list()

for (species in target_species$Species) {
  
  for (cohort in cohorts) {
    
    d <- copy(
      vdat[
        Species == species &
          Cohort == cohort
      ]
    )
    
    prepared <- prepare_model(d)
    
    d <- prepared$data
    covars <- prepared$covars
    
    stopifnot(
      sum(d$CRC == 1) >= 8,
      sum(d$CRC == 0) >= 8
    )
    
    formula_terms <- c("CRC", covars)
    
    model_formula <- reformulate(
      formula_terms,
      response = "Detected"
    )
    
    abundance_formula <- reformulate(
      formula_terms,
      response = "CLR"
    )
    
    # A. Bias-reduced logistic regression
    
    detection_fit <- glm(
      model_formula,
      data = d,
      family = binomial(),
      method = brglm2::brglmFit,
      type = "AS_mean"
    )
    
    detection_summary <- coef(
      summary(detection_fit)
    )
    
    detection_beta <- detection_summary[
      "CRC", "Estimate"
    ]
    
    detection_se <- detection_summary[
      "CRC", "Std. Error"
    ]
    
    # B. Linear association with robust HC3 SE
    
    abundance_fit <- lm(
      abundance_formula,
      data = d
    )
    
    abundance_beta <- unname(
      coef(abundance_fit)["CRC"]
    )
    
    abundance_cov <- sandwich::vcovHC(
      abundance_fit,
      type = "HC3"
    )
    
    abundance_se <- sqrt(
      abundance_cov["CRC", "CRC"]
    )
    
    # Store effect estimates
    
    model_results[[length(model_results) + 1L]] <-
      data.table(
        Species = species,
        Cohort = cohort,
        N = nrow(d),
        CRC_N = sum(d$CRC == 1),
        Healthy_N = sum(d$CRC == 0),
        Covariates = if (
          length(covars) == 0
        ) "none" else paste(covars, collapse = "+"),
        Detection_logOR = detection_beta,
        Detection_SE = detection_se,
        Detection_OR = exp(detection_beta),
        CLR_beta = abundance_beta,
        CLR_SE = abundance_se
      )
    
    model_diagnostics[[length(model_diagnostics) + 1L]] <-
      data.table(
        Species = species,
        Cohort = cohort,
        Detection_CRC = sum(
          d$Detected[d$CRC == 1]
        ),
        Detection_healthy = sum(
          d$Detected[d$CRC == 0]
        ),
        Covariates = paste(
          covars,
          collapse = "+"
        )
      )
  }
}

cohort_models <- rbindlist(model_results)

diagnostics <- rbindlist(model_diagnostics)

stopifnot(
  all(is.finite(cohort_models$Detection_logOR)),
  all(is.finite(cohort_models$Detection_SE)),
  all(cohort_models$Detection_SE > 0),
  all(is.finite(cohort_models$CLR_beta)),
  all(is.finite(cohort_models$CLR_SE)),
  all(cohort_models$CLR_SE > 0)
)

#
# REML estimation with Hartung-Knapp inference.
#
# Small number of cohorts: interpret I2 and tau2
# cautiously, and do not interpret a non-significant
# heterogeneity test as evidence of homogeneity.

meta_results <- list()

for (species in target_species$Species) {
  
  d <- cohort_models[
    Species == species
  ]
  
  for (outcome in c("Detection", "CLR")) {
    
    if (outcome == "Detection") {
      
      yi <- d$Detection_logOR
      sei <- d$Detection_SE
      
    } else {
      
      yi <- d$CLR_beta
      sei <- d$CLR_SE
      
    }
    
    fit <- metafor::rma(
      yi = yi,
      sei = sei,
      method = "REML",
      test = "knha"
    )
    
    expected <- target_species[
      Species == species,
      Expected
    ][1]
    
    beta <- as.numeric(fit$b[1, 1])
    
    meta_results[[length(meta_results) + 1L]] <-
      data.table(
        Species = species,
        Outcome = outcome,
        K = fit$k,
        Beta = beta,
        CI_lower = fit$ci.lb,
        CI_upper = fit$ci.ub,
        P_value = fit$pval,
        I2 = fit$I2,
        Tau2 = fit$tau2,
        Expected_direction = expected,
        Direction_agreement =
          sign(beta) == expected,
        Cohorts_expected_direction =
          sum(sign(yi) == expected),
        Cohorts_evaluated = length(yi)
      )
  }
}

meta_summary <- rbindlist(meta_results)

# Separate BH correction by outcome.
# Each outcome evaluates the same six
# prespecified microbial species.

meta_summary[
  ,
  FDR := p.adjust(P_value, method = "BH"),
  by = Outcome
]

meta_summary[
  ,
  Association_scale := fifelse(
    Outcome == "Detection",
    "Log odds ratio",
    "CLR difference"
  )
]

meta_summary[
  ,
  OR := fifelse(
    Outcome == "Detection",
    exp(Beta),
    NA_real_
  )
]

meta_summary[
  ,
  OR_lower := fifelse(
    Outcome == "Detection",
    exp(CI_lower),
    NA_real_
  )
]

meta_summary[
  ,
  OR_upper := fifelse(
    Outcome == "Detection",
    exp(CI_upper),
    NA_real_
  )
]

setorder(
  meta_summary,
  Outcome,
  FDR
)


fwrite(
  cohort_models,
  file.path(
    output_dir,
    "taxonomic_validation_9_adjusted_cohort_models.csv"
  )
)

fwrite(
  diagnostics,
  file.path(
    output_dir,
    "taxonomic_validation_9_model_diagnostics.csv"
  )
)

fwrite(
  meta_summary,
  file.path(
    output_dir,
    "taxonomic_validation_9_random_effects_meta.csv"
  )
)

saveRDS(
  list(
    cohort_models = cohort_models,
    meta_summary = meta_summary
  ),
  file.path(
    output_dir,
    "taxonomic_validation_9_validation_results.rds"
  )
)

# SENSITIVITY ANALYSIS EXCLUDING NHSII
#
# No repeat of species processing or model fitting.

library(data.table)
library(metafor)

project_dir <- project_root()

output_dir <- file.path(
  project_dir,
  "results",
  "taxonomic_validation"
)

cohort_models <- fread(
  file.path(
    output_dir,
    "taxonomic_validation_9_adjusted_cohort_models.csv"
  )
)

original_meta <- fread(
  file.path(
    output_dir,
    "taxonomic_validation_9_random_effects_meta.csv"
  )
)


sensitivity_data <- cohort_models[
  !grepl("NSHII", Cohort)
]

stopifnot(
  uniqueN(sensitivity_data$Cohort) == 3
)

expected_direction <- data.table(
  Species = c(
    "Parvimonas micra",
    "Peptostreptococcus stomatis",
    "Gemella morbillorum",
    "Porphyromonas asaccharolytica",
    "Roseburia intestinalis",
    "Streptococcus salivarius"
  ),
  Expected = c(1, 1, 1, 1, -1, -1)
)


results <- list()

for (species in expected_direction$Species) {
  
  d <- sensitivity_data[
    Species == species
  ]
  
  expected <- expected_direction[
    Species == species,
    Expected
  ][1]
  
  for (outcome in c("Detection", "CLR")) {
    
    if (outcome == "Detection") {
      
      effect <- d$Detection_logOR
      se <- d$Detection_SE
      
    } else {
      
      effect <- d$CLR_beta
      se <- d$CLR_SE
    }
    
    fit <- metafor::rma(
      yi = effect,
      sei = se,
      method = "REML",
      test = "knha"
    )
    
    beta <- as.numeric(fit$b[1, 1])
    
    results[[length(results) + 1L]] <-
      data.table(
        Species = species,
        Outcome = outcome,
        Cohorts = fit$k,
        Beta = beta,
        CI_lower = fit$ci.lb,
        CI_upper = fit$ci.ub,
        P_value = fit$pval,
        I2 = fit$I2,
        Direction_agreement = sign(beta) == expected,
        Cohorts_expected_direction =
          sum(sign(effect) == expected)
      )
  }
}

sensitivity_results <- rbindlist(results)

sensitivity_results[
  ,
  FDR := p.adjust(P_value, method = "BH"),
  by = Outcome
]


sensitivity_results[
  ,
  OR := fifelse(
    Outcome == "Detection",
    exp(Beta),
    NA_real_
  )
]

sensitivity_results[
  ,
  OR_lower := fifelse(
    Outcome == "Detection",
    exp(CI_lower),
    NA_real_
  )
]

sensitivity_results[
  ,
  OR_upper := fifelse(
    Outcome == "Detection",
    exp(CI_upper),
    NA_real_
  )
]


comparison <- merge(
  original_meta[
    ,
    .(
      Species,
      Outcome,
      Original_Beta = Beta,
      Original_P = P_value,
      Original_FDR = FDR
    )
  ],
  sensitivity_results[
    ,
    .(
      Species,
      Outcome,
      Sensitivity_Beta = Beta,
      Sensitivity_P = P_value,
      Sensitivity_FDR = FDR,
      Direction_agreement
    )
  ],
  by = c("Species", "Outcome")
)


fwrite(
  sensitivity_results,
  file.path(
    output_dir,
    "taxonomic_validation_10_excluding_NHSII.csv"
  )
)

fwrite(
  comparison,
  file.path(
    output_dir,
    "taxonomic_validation_10_original_vs_sensitivity.csv"
  )
)

# INDEPENDENT FUNCTIONAL VALIDATION
# HUMAnN 3.6 PATHWAY COMPATIBILITY

library(data.table)

project_dir <- project_root()

data_dir <- file.path(
  project_dir, "data", "validation"
)

output_dir <- file.path(
  project_dir, "results", "pathway_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)

archive <- file.path(
  data_dir,
  "humann3.6_pathways.tar.gz"
)

stopifnot(file.exists(archive))


entries <- untar(archive, list = TRUE)

selected <- entries[
  basename(entries) %in% paste0(cohorts, ".tsv")
]

stopifnot(length(selected) == 4L)

extract_dir <- file.path(
  data_dir, "functional_validation"
)

dir.create(
  extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

untar(
  archive,
  files = selected,
  exdir = extract_dir
)


keywords <- c(
  "tryptophan",
  "indole",
  "propionate",
  "propionic",
  "acrylate",
  "butyrate",
  "butyric",
  "fructose",
  "galactose",
  "methionine",
  "ornithine"
)

keyword_pattern <- paste(
  keywords,
  collapse = "|"
)

profile_summary <- list()
pathway_matches <- list()
pathway_examples <- list()

for (cohort in cohorts) {
  
  path <- file.path(
    extract_dir,
    "humann3.6_pathways",
    paste0(cohort, ".tsv")
  )
  
  stopifnot(file.exists(path))
  
  dt <- fread(
    path,
    showProgress = FALSE
  )
  
  feature_col <- names(dt)[1]
  
  features <- as.character(
    dt[[feature_col]]
  )
  
  # HUMAnN profiles may contain taxon-stratified
  # rows separated by "|".
  # Retain only community-level pathway rows.
  
  community_rows <- !grepl(
    "|",
    features,
    fixed = TRUE
  )
  
  community_features <- features[
    community_rows
  ]
  
  matched <- community_features[
    grepl(
      keyword_pattern,
      community_features,
      ignore.case = TRUE
    )
  ]
  
  pathway_matches[[cohort]] <- data.table(
    Cohort = cohort,
    Pathway = matched
  )
  
  profile_summary[[cohort]] <- data.table(
    Cohort = cohort,
    Total_rows = nrow(dt),
    Community_rows = sum(community_rows),
    Stratified_rows = sum(!community_rows),
    Profile_samples = ncol(dt) - 1L,
    Matching_pathways = length(matched),
    Feature_column = feature_col
  )
  
  pathway_examples[[cohort]] <- data.table(
    Cohort = cohort,
    Example = head(community_features, 8)
  )
}


functional_summary <- rbindlist(profile_summary)

functional_matches <- unique(
  rbindlist(pathway_matches)
)

pathway_coverage <- functional_matches[
  ,
  .(
    Cohorts_present = uniqueN(Cohort)
  ),
  by = Pathway
][order(-Cohorts_present, Pathway)]


fwrite(
  functional_summary,
  file.path(
    output_dir,
    "functional_profile_summary.csv"
  )
)

fwrite(
  pathway_coverage,
  file.path(
    output_dir,
    "candidate_functional_pathways.csv"
  )
)

# INDEPENDENT FUNCTIONAL VALIDATION
#
# Two biologically corresponding pathways:
# PWY-5494 — Propionate production (acrylate pathway)
# PWY-6317 — Galactose degradation (Leloir pathway)
#
# Outcome: community-relative HUMAnN pathway abundance
# Transformation: CLR across all community pathways
# Models: cohort-specific HC3 + random-effects meta-analysis
#
# Print all outputs at the end.

library(data.table)

project_dir <- project_root()

profile_dir <- file.path(
  project_dir,
  "data",
  "validation",
  "functional_validation",
  "humann3.6_pathways"
)

output_dir <- file.path(
  project_dir,
  "results",
  "pathway_validation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


targets <- data.table(
  Pathway_ID = c("PWY-5494", "PWY-6317"),
  Original_module = c("MF0121", "M00632"),
  Expected = c(1L, -1L)
)

cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)


metadata_files <- list.files(
  file.path(
    project_dir,
    "data",
    "validation",
    "metadata_extracted",
    "metadata"
  ),
  pattern = "\\.tsv$",
  full.names = TRUE
)

metadata_files <- metadata_files[
  tools::file_path_sans_ext(
    basename(metadata_files)
  ) %in% cohorts
]

stopifnot(length(metadata_files) == 4L)

metadata <- rbindlist(
  lapply(metadata_files, function(f) {
    x <- fread(f)
    x[, Cohort := tools::file_path_sans_ext(basename(f))]
    x
  }),
  fill = TRUE
)

metadata <- metadata[
  `Study condition` %in% c("CRC", "Control")
]

metadata[, CRC := as.integer(`Study condition` == "CRC")]

metadata[, Age := suppressWarnings(as.numeric(Age))]
metadata[, BMI := suppressWarnings(as.numeric(BMI))]

metadata[, Sex := trimws(as.character(Sex))]
metadata[!Sex %in% c("Male", "Female"), Sex := NA_character_]

stopifnot(
  !anyDuplicated(metadata[, .(Cohort, Name)])
)


records <- list()
profile_checks <- list()

for (cohort in cohorts) {
  
  path <- file.path(profile_dir, paste0(cohort, ".tsv"))
  stopifnot(file.exists(path))
  
  dt <- fread(path, showProgress = FALSE)
  
  pathway_col <- names(dt)[1]
  
  pathway_names <- as.character(dt[[pathway_col]])
  
  # Community-level pathways only.
  keep <- !grepl("|", pathway_names, fixed = TRUE)
  
  dt <- dt[keep]
  
  pathway_names <- pathway_names[keep]
  
  pathway_ids <- sub(":.*$", "", pathway_names)
  
  # Exclude bookkeeping features if present.
  valid <- !pathway_ids %in% c(
    "UNMAPPED",
    "UNINTEGRATED"
  )
  
  dt <- dt[valid]
  pathway_ids <- pathway_ids[valid]
  
  stopifnot(!anyDuplicated(pathway_ids))
  
  meta <- copy(metadata[Cohort == cohort])
  
  sample_ids <- meta$Name
  
  stopifnot(all(sample_ids %in% names(dt)))
  
  abundance <- as.matrix(dt[, ..sample_ids])
  
  storage.mode(abundance) <- "double"
  
  stopifnot(
    all(is.finite(abundance)),
    all(abundance >= 0)
  )
  
  # Community-relative pathway abundance.
  # HUMAnN pathway abundances are not necessarily
  # fractions that sum to one before closure.
  
  totals <- colSums(abundance)
  
  stopifnot(all(totals > 0))
  
  relative <- sweep(abundance, 2, totals, "/")
  
  # Fixed pseudocount for initial validation.
  # The choice will be documented as a limitation.
  
  log_relative <- log(relative + 1e-6)
  
  clr <- sweep(
    log_relative,
    2,
    colMeans(log_relative),
    "-"
  )
  
  rownames(clr) <- pathway_ids
  
  for (i in seq_len(nrow(targets))) {
    
    id <- targets$Pathway_ID[i]
    
    if (!id %in% rownames(clr)) {
      stop("Missing pathway: ", id, " / ", cohort)
    }
    
    records[[length(records) + 1L]] <- data.table(
      Cohort = cohort,
      Sample = sample_ids,
      Pathway_ID = id,
      CLR = as.numeric(clr[id, ])
    )
  }
  
  profile_checks[[cohort]] <- data.table(
    Cohort = cohort,
    Samples = length(sample_ids),
    Community_pathways = nrow(dt),
    Target_pathways = sum(
      targets$Pathway_ID %in% pathway_ids
    )
  )
}

functional_data <- rbindlist(records)

functional_data <- merge(
  functional_data,
  metadata[
    ,
    .(
      Cohort,
      Sample = Name,
      CRC,
      Age,
      Sex,
      BMI
    )
  ],
  by = c("Cohort", "Sample"),
  all.x = TRUE,
  sort = FALSE
)

stopifnot(
  !anyNA(functional_data$CRC),
  !anyNA(functional_data$CLR)
)


model_results <- list()

for (id in targets$Pathway_ID) {
  
  for (cohort in cohorts) {
    
    d <- copy(
      functional_data[
        Pathway_ID == id &
          Cohort == cohort
      ]
    )
    
    # Use cohort-specific available covariates.
    # Avoid losing the majority of samples because
    # one demographic variable is poorly recorded.
    
    candidate_covariates <- c("Age", "Sex", "BMI")
    
    covars <- candidate_covariates[
      vapply(
        candidate_covariates,
        function(v) {
          observed <- d[[v]][!is.na(d[[v]])]
          length(unique(observed)) > 1
        },
        logical(1)
      )
    ]
    
    repeat {
      
      needed <- c("CRC", "CLR", covars)
      
      cc <- complete.cases(d[, ..needed])
      
      if (
        (sum(d$CRC[cc] == 1) >= 8 &&
         sum(d$CRC[cc] == 0) >= 8) ||
        length(covars) == 0L
      ) {
        break
      }
      
      missing_fraction <- vapply(
        covars,
        function(v) mean(is.na(d[[v]])),
        numeric(1)
      )
      
      covars <- setdiff(
        covars,
        covars[which.max(missing_fraction)]
      )
    }
    
    needed <- c("CRC", "CLR", covars)
    
    d <- d[complete.cases(d[, ..needed])]
    
    covars <- covars[
      vapply(
        covars,
        function(v) length(unique(d[[v]])) > 1,
        logical(1)
      )
    ]
    
    if ("Sex" %in% covars) {
      d[, Sex := droplevels(factor(Sex))]
    }
    
    stopifnot(
      sum(d$CRC == 1) >= 8,
      sum(d$CRC == 0) >= 8
    )
    
    model_formula <- reformulate(
      c("CRC", covars),
      response = "CLR"
    )
    
    fit <- lm(model_formula, data = d)
    
    robust_cov <- sandwich::vcovHC(fit, type = "HC3")
    
    beta <- unname(coef(fit)["CRC"])
    
    se <- sqrt(robust_cov["CRC", "CRC"])
    
    stopifnot(
      is.finite(beta),
      is.finite(se),
      se > 0
    )
    
    model_results[[length(model_results) + 1L]] <-
      data.table(
        Pathway_ID = id,
        Cohort = cohort,
        N = nrow(d),
        CRC_N = sum(d$CRC == 1),
        Healthy_N = sum(d$CRC == 0),
        Covariates = if (
          length(covars)
        ) paste(covars, collapse = "+") else "none",
        Beta = beta,
        SE = se,
        Expected_direction = targets[
          Pathway_ID == id,
          Expected
        ][1],
        Direction_agreement =
          sign(beta) == targets[
            Pathway_ID == id,
            Expected
          ][1]
      )
  }
}

cohort_results <- rbindlist(model_results)


meta_results <- list()

for (id in targets$Pathway_ID) {
  
  d <- cohort_results[Pathway_ID == id]
  
  fit <- metafor::rma(
    yi = d$Beta,
    sei = d$SE,
    method = "REML",
    test = "knha"
  )
  
  beta <- as.numeric(fit$b[1, 1])
  
  expected <- targets[
    Pathway_ID == id,
    Expected
  ][1]
  
  meta_results[[id]] <- data.table(
    Pathway_ID = id,
    K = fit$k,
    Beta = beta,
    CI_lower = fit$ci.lb,
    CI_upper = fit$ci.ub,
    P_value = fit$pval,
    I2 = fit$I2,
    Tau2 = fit$tau2,
    Expected_direction = expected,
    Direction_agreement = sign(beta) == expected,
    Cohorts_expected_direction =
      sum(d$Direction_agreement)
  )
}

meta_results <- rbindlist(meta_results)

meta_results[
  ,
  FDR := p.adjust(P_value, method = "BH")
]


fwrite(
  functional_data,
  file.path(
    output_dir,
    "pathway_validation_2_functional_validation_data.csv"
  )
)

fwrite(
  cohort_results,
  file.path(
    output_dir,
    "pathway_validation_2_cohort_models.csv"
  )
)

fwrite(
  meta_results,
  file.path(
    output_dir,
    "pathway_validation_2_random_effects_meta.csv"
  )
)

# FINAL ANALYSIS — MICROBIAL COMMUNITY COMPOSITION
#
# Independent cohorts: Czechia, Italy, NHSII, Turkey
#
# Primary analysis:
# Bray-Curtis PERMANOVA, controlling for cohort
# and restricting permutations within cohort.
#
# Secondary:
# Age + sex + BMI adjustment
# Cohort-specific PERMANOVA
# Cohort-specific dispersion testing
#
# Biological replication: individual participants

library(data.table)

library(vegan)

set.seed(20261027)

project_dir <- project_root()

profile_dir <- file.path(
  project_dir, "data", "validation",
  "taxonomic_validation", "metaphlan4_profiles"
)

metadata_dir <- file.path(
  project_dir, "data", "validation",
  "metadata_extracted", "metadata"
)

output_dir <- file.path(project_dir, "results", "community_validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c(
  "This_study_cohort3__ONCOBIOME_IIGM_CZ",
  "This_study_cohort4__ONCOBIOME_IIGM_IT",
  "This_study_cohort5__NSHII",
  "This_study_cohort6__IIGM_TU"
)


abundance_list <- list()
metadata_list <- list()

for (cohort in cohorts) {
  
  profile <- fread(
    file.path(profile_dir, paste0(cohort, ".tsv")),
    showProgress = FALSE
  )
  
  meta <- fread(
    file.path(metadata_dir, paste0(cohort, ".tsv"))
  )
  
  meta <- meta[
    `Study condition` %in% c("CRC", "Control")
  ]
  
  stopifnot(
    !anyDuplicated(meta$Name),
    all(meta$Name %in% names(profile))
  )
  
  species_rows <- grepl(
    "(^|\\|)s__[^|]+$",
    profile$clade_name
  )
  
  sp <- profile[species_rows]
  
  species <- gsub(
    "_", " ",
    sub("^.*s__", "", sp$clade_name),
    fixed = TRUE
  )
  
  stopifnot(!anyDuplicated(species))
  
  ids <- meta$Name
  
  matrix_sp <- t(as.matrix(sp[, ..ids]))
  storage.mode(matrix_sp) <- "double"
  
  rownames(matrix_sp) <- paste(cohort, ids, sep = "::")
  colnames(matrix_sp) <- species
  
  stopifnot(
    all(is.finite(matrix_sp)),
    all(matrix_sp >= 0),
    all(rowSums(matrix_sp) > 0)
  )
  
  # Close species-level abundances to proportions.
  matrix_sp <- decostand(matrix_sp, method = "total")
  
  abundance_list[[cohort]] <- matrix_sp
  
  metadata_list[[cohort]] <- data.table(
    Sample_key = rownames(matrix_sp),
    Cohort = cohort,
    Group = ifelse(
      meta[["Study condition"]] == "CRC",
      "CRC", "healthy"
    ),
    Age = suppressWarnings(as.numeric(meta$Age)),
    BMI = suppressWarnings(as.numeric(meta$BMI)),
    Sex = trimws(as.character(meta$Sex))
  )
}


all_species <- sort(unique(unlist(
  lapply(abundance_list, colnames)
)))

community_parts <- lapply(abundance_list, function(m) {
  
  aligned <- matrix(
    0,
    nrow = nrow(m),
    ncol = length(all_species),
    dimnames = list(rownames(m), all_species)
  )
  
  aligned[, colnames(m)] <- m
  aligned
})

community <- do.call(rbind, community_parts)

metadata <- rbindlist(metadata_list)

stopifnot(
  identical(rownames(community), metadata$Sample_key),
  !anyDuplicated(metadata$Sample_key)
)

metadata[
  !Sex %in% c("Male", "Female"),
  Sex := NA_character_
]

metadata[, Cohort := factor(Cohort)]
metadata[, Group := factor(Group, levels = c("healthy", "CRC"))]


distance <- vegdist(
  community,
  method = "bray"
)

#
# The model includes cohort before CRC.
# The CRC term is evaluated conditionally.
# Permutations are restricted to within-cohort samples.

permutation_design <- permute::how(
  nperm = 999,
  blocks = metadata$Cohort
)

primary_fit <- adonis2(
  distance ~ Cohort + Group,
  data = metadata,
  permutations = permutation_design,
  by = "terms"
)


adjusted_index <- complete.cases(
  metadata[, .(Age, Sex, BMI)]
)

meta_adjusted <- droplevels(
  as.data.frame(metadata[adjusted_index])
)

community_adjusted <- community[adjusted_index, , drop = FALSE]

distance_adjusted <- vegdist(
  community_adjusted,
  method = "bray"
)

stopifnot(
  all(table(
    meta_adjusted$Cohort,
    meta_adjusted$Group
  ) > 0)
)

adjusted_permutations <- permute::how(
  nperm = 999,
  blocks = meta_adjusted$Cohort
)

adjusted_fit <- adonis2(
  distance_adjusted ~ Cohort + Age + Sex + BMI + Group,
  data = meta_adjusted,
  permutations = adjusted_permutations,
  by = "terms"
)


cohort_results <- list()

for (cohort in cohorts) {
  
  idx <- as.character(metadata$Cohort) == cohort
  
  m <- community[idx, , drop = FALSE]
  
  md <- droplevels(
    as.data.frame(metadata[idx])
  )
  
  d <- vegdist(m, method = "bray")
  
  fit <- adonis2(
    d ~ Group,
    data = md,
    permutations = 999
  )
  
  # Dispersion checks whether within-group
  # variability differs between CRC and controls.
  
  dispersion <- betadisper(
    d,
    group = md$Group,
    type = "median",
    bias.adjust = TRUE
  )
  
  dispersion_test <- permutest(
    dispersion,
    permutations = 999
  )
  
  cohort_results[[cohort]] <- data.table(
    Cohort = cohort,
    N = nrow(md),
    CRC_N = sum(md$Group == "CRC"),
    Healthy_N = sum(md$Group == "healthy"),
    PERMANOVA_R2 = fit$R2[1],
    PERMANOVA_p = fit$`Pr(>F)`[1],
    Dispersion_p =
      dispersion_test$tab[1, "Pr(>F)"]
  )
}

cohort_results <- rbindlist(cohort_results)

cohort_results[
  ,
  PERMANOVA_FDR := p.adjust(
    PERMANOVA_p,
    method = "BH"
  )
]

cohort_results[
  ,
  Dispersion_FDR := p.adjust(
    Dispersion_p,
    method = "BH"
  )
]


fwrite(
  cohort_results,
  file.path(
    output_dir,
    "cohort_beta_diversity_results.csv"
  )
)

saveRDS(
  list(
    primary = primary_fit,
    adjusted = adjusted_fit,
    cohorts = cohort_results,
    sample_metadata = metadata,
    adjusted_sample_metadata = meta_adjusted,
    permutation_count = 999L
  ),
  file.path(
    output_dir,
    "beta_diversity_results.rds"
  )
)
