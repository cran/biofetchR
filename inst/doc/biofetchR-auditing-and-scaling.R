## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = TRUE,
  warning = TRUE
)

# Helper for calculating retention rates safely.
safe_rate <- function(numerator, denominator) {
  ifelse(
    !is.na(denominator) & denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

# Helper for compact display of selected columns.
select_existing <- function(x, columns) {
  x[
    ,
    intersect(columns, names(x)),
    drop = FALSE
  ]
}

# Vignette-safe example paths. These keep all example outputs and caches under
# tempdir(), which is safe for CRAN examples/vignettes/tests.
example_output_dir <- file.path(tempdir(), "biofetchR-vignette-outputs")
example_cache_dir <- file.path(tempdir(), "biofetchR-vignette-cache")

## ----load-package-------------------------------------------------------------
library(biofetchR)

## ----summary-column-reference, echo = FALSE-----------------------------------
summary_columns <- data.frame(
  column_group = c(
    "Input identifiers",
    "Spatial identifiers",
    "Record counts",
    "Status fields",
    "Output fields",
    "Evidence fields"
  ),
  example_columns = c(
    "species, iso2c, region_id, overlay",
    "region_source, gadm_unit, overlay, marine_region_name",
    "n_total, n_cleaned, n_thinned, n_with_overlay, n_without_overlay",
    "status, fail_stage, fail_reason",
    "output_file, output_dir",
    "native_filter_mode, griis_filter_mode, evidence_action"
  ),
  why_check = c(
    "Confirm which row was attempted.",
    "Confirm which spatial framework was used.",
    "Check how many records survived each processing step.",
    "Identify rows that succeeded, failed, or returned no usable records.",
    "Confirm where exported occurrence files were written.",
    "Confirm how native-range or invasive-status evidence was handled."
  ),
  stringsAsFactors = FALSE
)

summary_columns

## ----example-summary----------------------------------------------------------
summary_example <- data.frame(
  species = c(
    "Xenopus laevis",
    "Trachemys scripta",
    "Rattus rattus",
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida"
  ),
  workflow = c(
    "terrestrial_gadm",
    "terrestrial_gadm",
    "terrestrial_gadm",
    "marine_eez",
    "marine_eez",
    "marine_lme"
  ),
  region_id = c(
    "FR",
    "GB",
    "GB",
    NA,
    NA,
    NA
  ),
  overlay = c(
    NA,
    NA,
    NA,
    "eez",
    "eez",
    "lme"
  ),
  n_total = c(
    1340,
    220,
    0,
    980,
    640,
    410
  ),
  n_cleaned = c(
    1309,
    198,
    0,
    944,
    612,
    389
  ),
  n_thinned = c(
    861,
    174,
    0,
    702,
    455,
    301
  ),
  n_with_overlay = c(
    861,
    174,
    0,
    690,
    430,
    292
  ),
  n_without_overlay = c(
    0,
    0,
    0,
    12,
    25,
    9
  ),
  status = c(
    "success",
    "success",
    "no_records",
    "success",
    "success",
    "success"
  ),
  fail_stage = c(
    NA,
    NA,
    "import",
    NA,
    NA,
    NA
  ),
  fail_reason = c(
    NA,
    NA,
    "No records returned",
    NA,
    NA,
    NA
  ),
  output_file = c(
    file.path(
      example_output_dir,
      "terrestrial_gadm_batch",
      "occurrence_files",
      "Xenopus_laevis_FR.csv"
    ),
    file.path(
      example_output_dir,
      "terrestrial_gadm_batch",
      "occurrence_files",
      "Trachemys_scripta_GB.csv"
    ),
    NA,
    file.path(
      example_output_dir,
      "marine_eez_batch",
      "occurrence_files",
      "Carcinus_maenas_eez.csv"
    ),
    file.path(
      example_output_dir,
      "marine_eez_batch",
      "occurrence_files",
      "Mnemiopsis_leidyi_eez.csv"
    ),
    file.path(
      example_output_dir,
      "marine_lme_batch",
      "occurrence_files",
      "Undaria_pinnatifida_lme.csv"
    )
  ),
  stringsAsFactors = FALSE
)

summary_example

## ----read-summary-template, eval = FALSE--------------------------------------
# summary_tbl <- read.csv(
#   file.path(example_output_dir, "terrestrial_gadm_batch", "gbif_summary.csv"),
#   stringsAsFactors = FALSE
# )

## ----status-check-------------------------------------------------------------
table(summary_example$status, useNA = "ifany")

## ----failed-rows--------------------------------------------------------------
failed_rows <- summary_example[
  !summary_example$status %in% "success",
  ,
  drop = FALSE
]

failed_rows

## ----failure-diagnostic-fields------------------------------------------------
select_existing(
  failed_rows,
  c(
    "species",
    "workflow",
    "region_id",
    "overlay",
    "status",
    "fail_stage",
    "fail_reason"
  )
)

## ----retention-checks---------------------------------------------------------
summary_example$retained_after_cleaning <- safe_rate(
  summary_example$n_cleaned,
  summary_example$n_total
)

summary_example$retained_after_thinning <- safe_rate(
  summary_example$n_thinned,
  summary_example$n_cleaned
)

summary_example$overlay_assignment_rate <- safe_rate(
  summary_example$n_with_overlay,
  summary_example$n_thinned
)

select_existing(
  summary_example,
  c(
    "species",
    "workflow",
    "n_total",
    "n_cleaned",
    "n_thinned",
    "n_with_overlay",
    "n_without_overlay",
    "retained_after_cleaning",
    "retained_after_thinning",
    "overlay_assignment_rate"
  )
)

## ----review-flags-------------------------------------------------------------
summary_example$needs_review <- with(
  summary_example,
  status != "success" |
    retained_after_cleaning < 0.5 |
    retained_after_thinning < 0.5 |
    overlay_assignment_rate < 0.9
)

select_existing(
  summary_example[summary_example$needs_review, , drop = FALSE],
  c(
    "species",
    "workflow",
    "status",
    "n_total",
    "n_cleaned",
    "n_thinned",
    "n_with_overlay",
    "n_without_overlay",
    "retained_after_cleaning",
    "retained_after_thinning",
    "overlay_assignment_rate",
    "fail_stage",
    "fail_reason"
  )
)

## ----output-path-checks-------------------------------------------------------
summary_example$output_file_listed <- !is.na(summary_example$output_file) &
  nzchar(summary_example$output_file)

summary_example$output_file_exists <- ifelse(
  summary_example$output_file_listed,
  file.exists(summary_example$output_file),
  FALSE
)

select_existing(
  summary_example,
  c(
    "species",
    "status",
    "output_file_listed",
    "output_file_exists",
    "output_file"
  )
)

## ----rows-to-rerun------------------------------------------------------------
rows_to_rerun <- summary_example[
  !summary_example$status %in% "success",
  ,
  drop = FALSE
]

terrestrial_rows_to_rerun <- rows_to_rerun[
  grepl("^terrestrial|^freshwater", rows_to_rerun$workflow),
  ,
  drop = FALSE
]

marine_rows_to_rerun <- rows_to_rerun[
  grepl("^marine", rows_to_rerun$workflow),
  ,
  drop = FALSE
]

select_existing(
  rows_to_rerun,
  c(
    "species",
    "workflow",
    "region_id",
    "overlay",
    "fail_stage",
    "fail_reason"
  )
)

## ----rerun-failed-terrestrial-template, eval = FALSE--------------------------
# failed_terrestrial_input <- data.frame(
#   species = terrestrial_rows_to_rerun$species,
#   iso2c = terrestrial_rows_to_rerun$region_id,
#   stringsAsFactors = FALSE
# )
# 
# process_gbif_terrestrial_freshwater_pipeline(
#   df = failed_terrestrial_input,
# 
#   output_dir = file.path(example_output_dir, "terrestrial_gadm_rerun"),
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   region_source = "gadm",
#   gadm_unit = 1,
# 
#   batch_size = 1,
#   apply_cleaning = TRUE,
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE
# )

## ----rerun-failed-marine-template, eval = FALSE-------------------------------
# failed_marine_input <- data.frame(
#   species = marine_rows_to_rerun$species,
#   stringsAsFactors = FALSE
# )
# 
# process_gbif_marine_pipeline(
#   df = failed_marine_input,
# 
#   output_dir = file.path(example_output_dir, "marine_eez_rerun"),
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   overlay = "eez",
#   overlay_cache_dir = file.path(example_cache_dir, "marine_overlays"),
# 
#   batch_size = 1,
#   apply_cleaning = TRUE,
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE
# )

## ----origin-audit-example-----------------------------------------------------
origin_audit_example <- data.frame(
  species = c(
    "Xenopus laevis",
    "Trachemys scripta",
    "Carcinus maenas",
    "Undaria pinnatifida"
  ),
  recipient_region = c(
    "FR",
    "GB",
    "GB",
    "NZ"
  ),
  evidence_source = c(
    "native_web_and_griis",
    "native_web_and_griis",
    "griis_country_and_native_web",
    "griis_country_and_native_web"
  ),
  native_range_decision = c(
    "outside_native_range",
    "outside_native_range",
    "outside_native_range",
    "outside_native_range"
  ),
  griis_decision = c(
    "listed_invasive",
    "listed_invasive",
    "listed_alien_or_invasive",
    "listed_alien_or_invasive"
  ),
  evidence_action = c(
    "retain",
    "retain",
    "retain",
    "review"
  ),
  stringsAsFactors = FALSE
)

origin_audit_example

## ----spatial-attribution-audit-example----------------------------------------
spatial_audit_example <- data.frame(
  species = c(
    "Xenopus laevis",
    "Carcinus maenas",
    "Carcinus maenas",
    "Undaria pinnatifida"
  ),
  workflow = c(
    "terrestrial_gadm",
    "marine_eez",
    "marine_lme",
    "marine_ecs"
  ),
  spatial_framework = c(
    "GADM level 1",
    "Exclusive Economic Zones",
    "Large Marine Ecosystems",
    "Extended continental shelves"
  ),
  records_with_assignment = c(
    861,
    690,
    675,
    288
  ),
  records_without_assignment = c(
    0,
    12,
    27,
    13
  ),
  stringsAsFactors = FALSE
)

spatial_audit_example

## ----second-stage-terrestrial-input-------------------------------------------
terrestrial_test_batch <- data.frame(
  species = c(
    "Xenopus laevis",
    "Trachemys scripta",
    "Sturnus vulgaris",
    "Rattus rattus",
    "Myocastor coypus"
  ),
  iso2c = c(
    "FR",
    "GB",
    "DE",
    "GB",
    "IT"
  ),
  stringsAsFactors = FALSE
)

terrestrial_test_batch

## ----second-stage-marine-input------------------------------------------------
marine_test_batch <- data.frame(
  species = c(
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida",
    "Didemnum vexillum"
  ),
  stringsAsFactors = FALSE
)

marine_test_batch

## ----project-scale-template, eval = FALSE-------------------------------------
# project_output <- file.path(example_output_dir, "terrestrial_gadm_batch")
# project_cache <- example_cache_dir
# 
# batch_result <- process_gbif_terrestrial_freshwater_pipeline(
#   df = terrestrial_test_batch,
# 
#   output_dir = project_output,
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   region_source = "gadm",
#   gadm_unit = 1,
#   cache_dir_gadm = file.path(project_cache, "gadm"),
# 
#   batch_size = 1,
# 
#   apply_cleaning = TRUE,
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE,
#   quiet = FALSE
# )

## ----provider-metadata-template-----------------------------------------------
provider_metadata_template <- data.frame(
  provider = c(
    "GBIF",
    "GADM",
    "Marine Regions",
    "GRIIS",
    "SInAS",
    "WoRMS"
  ),
  used_for = c(
    "Occurrence records",
    "Administrative attribution",
    "Marine overlay attribution",
    "Species-country alien or invasive status evidence",
    "Native and alien distribution evidence",
    "Marine taxonomy or distribution evidence"
  ),
  record_with_outputs = c(
    "Download DOI or download key",
    "Version or access date",
    "Product name, version or access date, and overlay name",
    "Version or access date",
    "Dataset version and access date",
    "Dataset DOI, version, or access date"
  ),
  stringsAsFactors = FALSE
)

provider_metadata_template

## ----session-info-------------------------------------------------------------
sessionInfo()

## ----save-session-info-template, eval = FALSE---------------------------------
# dir.create(example_output_dir, recursive = TRUE, showWarnings = FALSE)
# 
# writeLines(
#   capture.output(sessionInfo()),
#   file.path(example_output_dir, "sessionInfo.txt")
# )

