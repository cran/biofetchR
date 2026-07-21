params <-
list(run_live = FALSE)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = TRUE,
  warning = TRUE
)

# Live examples are skipped by default.
# Set params$run_live to TRUE when rendering the vignette locally if you want to
# contact GBIF, GRIIS, SInAS, WoRMS, or other external provider sources.
run_live <- isTRUE(params$run_live)

# Read GBIF credentials from environment variables.
# This avoids writing credentials directly into the vignette or project files.
gbif_user <- Sys.getenv("GBIF_USER")
gbif_pwd <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")

# GBIF examples are attempted only when all credentials are available.
gbif_ready <- all(nzchar(c(gbif_user, gbif_pwd, gbif_email)))

# Example outputs are written to a temporary folder.
# In a real project, replace tempdir() with a stable project directory.
output_root <- file.path(tempdir(), "biofetchR_origin_evidence_vignette")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Run live examples without stopping the whole vignette if a provider request,
# download, or external service fails.
run_live_safely <- function(expr) {
  tryCatch(
    force(expr),
    error = function(e) {
      message("Live example did not complete: ", conditionMessage(e))
      NULL
    }
  )
}

# Print a clear message when a live example is skipped.
skip_live_message <- function(reason) {
  message("Live example skipped: ", reason)
  invisible(NULL)
}

## ----load-package-------------------------------------------------------------
library(biofetchR)

## ----render-live-version, eval = FALSE----------------------------------------
# rmarkdown::render(
#   "vignettes/biofetchR-origin-evidence.Rmd",
#   params = list(run_live = TRUE)
# )

## ----evidence-workflow-overview, echo = FALSE---------------------------------
workflow_steps <- data.frame(
  step = 1:7,
  stage = c(
    "Start with species-region records",
    "Prepare taxonomy",
    "Attach native-range evidence",
    "Attach GRIIS evidence",
    "Assign origin-status fields",
    "Inspect audit outputs",
    "Apply optional filtering"
  ),
  purpose = c(
    "Define the taxon and recipient country or region being evaluated.",
    "Clean names, apply manual corrections, and reject obvious non-taxonomic strings.",
    "Identify countries or regions treated as part of the taxon's native range.",
    "Identify species-country combinations listed as introduced or invasive.",
    "Classify rows using the available evidence and keep unresolved cases visible.",
    "Check taxonomy, evidence sources, status fields, and retained or excluded rows.",
    "Use stricter filters only when they match the biological question."
  ),
  stringsAsFactors = FALSE
)

workflow_steps

## ----taxonomy-correction-input------------------------------------------------
# Example input with deliberately mixed taxonomic quality.
# Some names are clean, some contain authority strings, some need manual
# correction, and some are not usable taxon names.
taxonomy_input <- data.frame(
  species = c(
    "Rattus rattus (Linnaeus, 1758)",
    "Sturnus vulgaris",
    "Xenopus laevis Daudin, 1802",
    "Trachemys scripta elegans",
    "Rattus sp.",
    "unknown",
    "water temperature"
  ),
  iso2c = c("GB", "DE", "FR", "GB", "GB", "GB", "FR"),
  stringsAsFactors = FALSE
)

taxonomy_input

## ----taxonomy-correction-example----------------------------------------------
# First, prepare the input names automatically.
# This handles routine cleaning and creates an audit trail.
taxonomy_prepared_auto <- bf_prepare_taxa_for_gbif(
  df = taxonomy_input,
  name_col = "species",
  iso2_col = "iso2c",
  require_species_level = TRUE,
  deduplicate = TRUE
)

# Inspect the automatic taxonomy audit.
taxonomy_prepared_auto$audit

# Optional: define manual overrides only for names the user has reviewed.
# These are not required for every workflow.
manual_taxonomy_fixes <- c(
  "Xenopus laevis Daudin, 1802" = "Xenopus laevis",
  "Trachemys scripta elegans" = "Trachemys scripta"
)

# Rerun taxonomy preparation with the reviewed manual overrides.
taxonomy_prepared_reviewed <- bf_prepare_taxa_for_gbif(
  df = taxonomy_input,
  name_col = "species",
  iso2_col = "iso2c",
  manual_fixes = manual_taxonomy_fixes,
  require_species_level = TRUE,
  deduplicate = TRUE
)

# Compare the reviewed audit table with the automatic audit table.
taxonomy_prepared_reviewed$audit

## ----local-recipient-input----------------------------------------------------
recipient_records <- data.frame(
  species = c(
    "Rattus rattus",
    "Sturnus vulgaris",
    "Xenopus laevis",
    "Trachemys scripta",
    "Carcinus maenas"
  ),
  iso2c = c("GB", "DE", "FR", "GB", "GB"),
  stringsAsFactors = FALSE
)

recipient_records

## ----status-category-table, echo = FALSE--------------------------------------
status_categories <- data.frame(
  status = c("native", "non_native", "invasive", "unresolved"),
  meaning = c(
    "The recipient region is supported as part of the taxon's native range.",
    "The taxon is interpreted as introduced or non-native, but invasive-status evidence is absent in the selected sources.",
    "The taxon is introduced or non-native and has invasive-status evidence, such as impact, rapid spread, or high abundance.",
    "The available evidence is not sufficient to assign a confident status."
  ),
  stringsAsFactors = FALSE
)

status_categories

## ----local-native-evidence----------------------------------------------------
# Example project native-range evidence table.
# Each row gives the countries treated as part of the taxon's native range.
# Multiple ISO3 country codes are separated with semicolons.
native_ranges <- data.frame(
  species = c(
    "Rattus rattus",
    "Sturnus vulgaris",
    "Xenopus laevis",
    "Trachemys scripta",
    "Carcinus maenas"
  ),
  native_origin_iso3 = c(
    "IND;PAK",
    "DEU;FRA;GBR",
    "ZAF;LSO;SWZ",
    "USA;MEX",
    "GBR;IRL;FRA;ESP;PRT"
  ),
  native_origin_flag = "HAS_ORIGIN",
  native_sources_used = "example_project_table",
  stringsAsFactors = FALSE
)

native_ranges

## ----attach-native-local------------------------------------------------------
native_status <- bf_attach_native_status(
  df = recipient_records,
  species_col = "species",
  iso2c_col = "iso2c",
  native_ranges = native_ranges,
  native_filter_mode = "audit_only",
  quiet = TRUE
)

native_status

## ----native-summary-local-----------------------------------------------------
bf_native_status_summary(native_status)

## ----inspect-native-audit-columns---------------------------------------------
native_audit_cols <- unique(c(
  intersect(c("species", "iso2c"), names(native_status)),
  grep("native|origin|source|status|decision|filter", names(native_status), value = TRUE)
))

native_status[, native_audit_cols, drop = FALSE]

## ----filtering-mode-table, echo = FALSE---------------------------------------
filtering_modes <- data.frame(
  mode = c(
    "audit_only",
    "non_native_only",
    "non_native_or_unknown",
    "native_only"
  ),
  behaviour = c(
    "Keep all rows and add evidence fields.",
    "Keep only rows interpreted as non-native in the recipient country.",
    "Keep rows interpreted as non-native plus unresolved rows.",
    "Keep only rows interpreted as native in the recipient country."
  ),
  typical_use = c(
    "First-pass review.",
    "Strict non-native-only analysis.",
    "Precautionary workflows where unresolved rows need later review.",
    "Native-range checks or validation workflows."
  ),
  stringsAsFactors = FALSE
)

filtering_modes

## ----native-filter-modes------------------------------------------------------
# Strict non-native-only table.
non_native_only <- bf_attach_native_status(
  df = recipient_records,
  species_col = "species",
  iso2c_col = "iso2c",
  native_ranges = native_ranges,
  native_filter_mode = "non_native_only",
  quiet = TRUE
)

# Less strict table that keeps unresolved rows for later review.
non_native_or_unknown <- bf_attach_native_status(
  df = recipient_records,
  species_col = "species",
  iso2c_col = "iso2c",
  native_ranges = native_ranges,
  native_filter_mode = "non_native_or_unknown",
  quiet = TRUE
)

non_native_only
non_native_or_unknown

## ----native-web-source-registry-----------------------------------------------
# Show native-range web/API sources supported by the installed package version.
bf_available_native_web_sources()

## ----native-web-live, eval = run_live-----------------------------------------
# # Retrieve native-range evidence from the selected web/API sources.
# # This is wrapped in run_live_safely() so that a temporary provider failure
# # does not stop the whole vignette from rendering.
# native_web <- run_live_safely(
# 
#   # Compile native-origin evidence for the species in recipient_records.
#   bf_fetch_native_ranges_web(
# 
#     # Species-country table containing the taxa to query.
#     species = recipient_records,
# 
#     # Column containing the taxon names.
#     species_col = "species",
# 
#     # Native-range evidence sources to query.
#     # These currently include SInAS, GBIF Species API distribution records,
#     # and WoRMS Aphia distribution records.
#     sources = c("sinas", "gbif", "worms"),
# 
#     # Cache provider results so repeated runs do not need to re-query
#     # the same sources unnecessarily.
#     cache_dir = file.path(output_root, "cache", "native_web"),
# 
#     # Use cached results when available.
#     # Set this to TRUE only when the provider evidence should be refreshed.
#     force_refresh = FALSE,
# 
#     # Pause briefly between provider requests to avoid sending rapid repeated
#     # queries to external services.
#     sleep_sec = 0.25,
# 
#     # Show provider messages during the live example.
#     quiet = FALSE,
# 
#     # Return a structured list containing summary, species-level,
#     # long-format, and unmapped evidence outputs.
#     return = "list"
#   )
# )
# 
# # If the live query completed successfully, inspect the returned evidence.
# if (!is.null(native_web)) {
# 
#   # Summary of the native-web retrieval.
#   native_web$summary
# 
#   # Species-level native-range evidence table.
#   # This is usually the table used for attachment to species-country records.
#   native_web$species
# 
#   # Long-format provider evidence.
#   # This is useful for checking which source returned which evidence.
#   native_web$long
# 
#   # Provider strings or regions that could not be mapped cleanly.
#   # These should be reviewed if many species remain unresolved.
#   native_web$unmapped
# }

## ----write-native-web-live, eval = run_live-----------------------------------
# if (exists("native_web") && !is.null(native_web)) {
#   bf_write_native_web_outputs(
#     native_web,
#     output_dir = file.path(output_root, "native_web_audit"),
#     prefix = "native_web"
#   )
# } else {
#   skip_live_message("native web evidence was not available from the previous live chunk")
# }

## ----attach-native-web-live, eval = run_live----------------------------------
# # Attach the species-level native-range evidence returned by the web/API
# # retrieval step to the species-country input table.
# # This chunk depends on the previous native_web object being available.
# if (exists("native_web") && !is.null(native_web)) {
# 
#   # Add native-status fields to the recipient species-country table.
#   native_web_status <- bf_attach_native_status(
# 
#     # Species-country table to classify.
#     df = recipient_records,
# 
#     # Column containing the taxon name in the species-country table.
#     species_col = "species",
# 
#     # Column containing the recipient country as an ISO2 code.
#     iso2c_col = "iso2c",
# 
#     # Species-level native-range evidence returned by bf_fetch_native_ranges_web().
#     # This table combines the selected web/API sources into a format that can be
#     # used by bf_attach_native_status().
#     native_ranges = native_web$species,
# 
#     # Column containing the taxon name in the native-range evidence table.
#     native_species_col = "species",
# 
#     # Keep all rows and add native-range evidence fields for inspection.
#     # No records are removed at this stage.
#     native_filter_mode = "audit_only",
# 
#     # Suppress routine messages from the attachment function.
#     quiet = TRUE
#   )
# 
#   # Display the species-country table with attached native-range evidence.
#   native_web_status
# 
# } else {
# 
#   # If the previous live retrieval step did not run or failed, print a clear
#   # skip message rather than causing the vignette to fail.
#   skip_live_message("native web evidence was not available from the previous live chunk")
# }

## ----attach-native-web-direct-live, eval = run_live---------------------------
# # Attach native-range evidence directly from supported web/API sources.
# # This route is useful when the workflow should retrieve provider evidence
# # during the native-status attachment step rather than supplying a prebuilt
# # native_ranges table.
# native_web_direct <- run_live_safely(
# 
#   # Attach native-status fields to the species-country input table.
#   bf_attach_native_status(
# 
#     # Species-country table to classify.
#     df = recipient_records,
# 
#     # Column containing the taxon name.
#     species_col = "species",
# 
#     # Column containing the recipient country as an ISO2 code.
#     iso2c_col = "iso2c",
# 
#     # No local native-range table is supplied here.
#     # Instead, evidence is retrieved from the selected web/API sources below.
#     native_ranges = NULL,
# 
#     # Turn on native-range evidence retrieval from supported provider sources.
#     use_native_web = TRUE,
# 
#     # Query all selected native-range web/API sources.
#     # SInAS, GBIF distribution metadata, and WoRMS Aphia distributions are used
#     # as complementary evidence sources.
#     native_web_sources = c("sinas", "gbif", "worms"),
# 
#     # Store downloaded or retrieved provider evidence in a cache directory.
#     # This avoids repeating the same provider requests unnecessarily.
#     native_web_cache_dir = file.path(output_root, "cache", "native_web_direct"),
# 
#     # Write native-web audit outputs so the provider evidence can be inspected.
#     export_native_web_audit = TRUE,
# 
#     # Directory where native-web audit files will be written.
#     native_web_audit_dir = file.path(output_root, "native_web_direct_audit"),
# 
#     # Keep all rows and attach evidence fields for review.
#     # No records are removed at this stage.
#     native_filter_mode = "audit_only",
# 
#     # Suppress routine messages from the attachment function.
#     quiet = TRUE
#   )
# )
# 
# # If the live provider query completed successfully, display the attached
# # native-status table.
# if (!is.null(native_web_direct)) {
#   native_web_direct
# }

## ----griis-live, eval = run_live----------------------------------------------
# # Attach GRIIS introduced/invasive-status evidence to the species-country table.
# # This is a live example because the function may download or retrieve GRIIS
# # evidence when a local GRIIS table is not supplied.
# griis_status <- run_live_safely(
# 
#   # Add GRIIS status fields to the input table.
#   bf_attach_griis_status(
# 
#     # Species-country table to classify.
#     df = recipient_records,
# 
#     # Column containing the taxon name.
#     species_col = "species",
# 
#     # Column containing the recipient country as an ISO2 code.
#     iso2c_col = "iso2c",
# 
#     # No preloaded GRIIS table is supplied here.
#     # When griis = NULL, the function uses the configured GRIIS retrieval/cache
#     # workflow to obtain the evidence.
#     griis = NULL,
# 
#     # Directory used to cache GRIIS files or processed evidence.
#     # This avoids repeating the same retrieval step unnecessarily.
#     cache_dir = file.path(output_root, "cache", "griis"),
# 
#     # Use cached GRIIS evidence when available.
#     # Set this to TRUE only when you deliberately want to refresh the cache.
#     force_refresh = FALSE,
# 
#     # Require the species-country combination to match GRIIS evidence.
#     # This is stricter than checking whether the species appears anywhere in GRIIS.
#     require_country = TRUE,
# 
#     # Show progress messages for the live example.
#     quiet = FALSE
#   )
# )
# 
# # If the GRIIS attachment step completed successfully, display the output table
# # with the added GRIIS evidence fields.
# if (!is.null(griis_status)) {
#   griis_status
# }

## ----griis-filter-live, eval = run_live---------------------------------------
# # Create a stricter invasive-only table from the already attached GRIIS fields.
# # The previous chunk created griis_status with columns such as griis_listed,
# # griis_invasive, and griis_status. Here we simply filter that existing audit
# # table rather than re-running the GRIIS attachment step.
# if (exists("griis_status") && !is.null(griis_status)) {
# 
#   # Keep only rows that were flagged as invasive by the attached GRIIS evidence.
#   invasive_only <- griis_status[
#     griis_status$griis_invasive %in% TRUE,
#     ,
#     drop = FALSE
#   ]
# 
#   # Display the stricter invasive-only table.
#   invasive_only
# 
# } else {
# 
#   # If the GRIIS attachment chunk did not run or failed, skip cleanly.
#   skip_live_message("GRIIS status was not available from the previous live chunk")
# }

## ----origin-status-audit-example----------------------------------------------
# Example origin-status audit table.
# This is not live evidence; it illustrates how evidence fields can support
# one final status assignment for each species-country row.
origin_status_audit_example <- data.frame(
  species = c(
    "Rattus rattus",
    "Sturnus vulgaris",
    "Xenopus laevis",
    "Trachemys scripta",
    "Carcinus maenas"
  ),
  iso2c = c("GB", "DE", "FR", "GB", "GB"),
  native_range_result = c(
    "outside_native_range",
    "inside_native_range",
    "outside_native_range",
    "outside_native_range",
    "inside_native_range"
  ),
  griis_result = c(
    "listed_invasive",
    "not_listed_or_unresolved",
    "listed_invasive",
    "listed_introduced",
    "not_listed_or_unresolved"
  ),
  assigned_status = c(
    "invasive",
    "native",
    "invasive",
    "non_native",
    "native"
  ),
  review_action = c(
    "retain_for_invasive_workflow",
    "exclude_from_non_native_workflow",
    "retain_for_invasive_workflow",
    "retain_for_non_native_workflow",
    "exclude_from_non_native_workflow"
  ),
  stringsAsFactors = FALSE
)

origin_status_audit_example

## ----recommended-defaults-table, echo = FALSE---------------------------------
recommended_defaults <- data.frame(
  setting = c(
    "prepare_taxonomy",
    "export_taxonomy_audit",
    "use_native_range",
    "use_native_web",
    "native_web_sources",
    "native_filter_mode",
    "use_griis",
    "filter_griis_invasive",
    "reconcile_origin_evidence",
    "export_origin_audit",
    "store_in_memory"
  ),
  recommended_value = c(
    "TRUE",
    "TRUE",
    "TRUE when origin status is needed",
    "TRUE to supplement project native-range evidence or fill gaps",
    "c(\"sinas\", \"gbif\", \"worms\")",
    "audit_only",
    "TRUE when GRIIS evidence is needed",
    "FALSE",
    "TRUE",
    "TRUE",
    "FALSE for larger batches"
  ),
  reason = c(
    "Clean and correct names before evidence matching or GBIF submission.",
    "Preserve raw, fixed, cleaned, accepted, and rejected names.",
    "Attach native-origin fields for status assignment.",
    "Compile native-origin evidence from supported provider helpers.",
    "Use all supported native-range web/API sources unless the project requires a subset.",
    "Keep all rows while evidence is being checked.",
    "Attach introduced and invasive-status fields.",
    "Avoid invasive-only filtering before evidence review.",
    "Create combined origin-evidence status fields where possible.",
    "Write files that explain status assignment and filtering decisions.",
    "Write occurrence outputs to disk rather than keeping every table in memory."
  ),
  stringsAsFactors = FALSE
)

recommended_defaults

## ----origin-in-pipeline-live, eval = run_live---------------------------------
# # Run the live pipeline example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   # Skip cleanly rather than failing the vignette.
#   skip_live_message("set GBIF credentials to run the origin-evidence pipeline example")
# 
# } else {
# 
#   # Small species-country input table for the live example.
#   origin_pipeline_input <- data.frame(
#     species = c(
#       "Xenopus laevis",
#       "Trachemys scripta"
#     ),
#     iso2c = c("FR", "GB"),
#     stringsAsFactors = FALSE
#   )
# 
#   # Run the terrestrial/freshwater pipeline with taxonomy preparation,
#   # native-range evidence, GRIIS evidence, and audit outputs enabled.
#   origin_pipeline_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
# 
#       # Species-country input table.
#       df = origin_pipeline_input,
# 
#       # Output folder for this vignette example.
#       output_dir = file.path(output_root, "origin_pipeline"),
# 
#       # GBIF credentials read from environment variables.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Use GADM level 1 for spatial attribution.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Prepare taxon names and write taxonomy audit files.
#       prepare_taxonomy = TRUE,
#       export_taxonomy_audit = TRUE,
# 
#       # Attach GRIIS introduced/invasive-status evidence.
#       # Do not filter to GRIIS-invasive rows in this audit-mode example.
#       use_griis = TRUE,
#       filter_griis_invasive = FALSE,
#       griis_cache_dir = file.path(output_root, "cache", "griis"),
# 
#       # Attach native-range evidence from supported web/API sources.
#       # No local native-range table is supplied in this example.
#       use_native_range = TRUE,
#       native_ranges = NULL,
#       use_native_web = TRUE,
#       native_web_sources = c("sinas", "gbif", "worms"),
#       native_web_cache_dir = file.path(output_root, "cache", "native_web"),
#       export_native_web_audit = TRUE,
# 
#       # Keep all rows while attaching native-range evidence.
#       native_filter_mode = "audit_only",
# 
#       # Write origin-status audit outputs.
#       reconcile_origin_evidence = TRUE,
#       export_origin_audit = TRUE,
#       export_summary = TRUE,
# 
#       # Keep the vignette example small.
#       batch_size = 1,
# 
#       # Clean and thin occurrence records before export.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Write detailed outputs to disk rather than storing every table in memory.
#       return_all_results = FALSE,
#       store_in_memory = FALSE,
# 
#       # Show progress messages during the live example.
#       quiet = FALSE
#     )
#   )
# 
#   # Display the returned pipeline result object.
#   origin_pipeline_result
# }

## ----inspect-origin-outputs-live, eval = run_live-----------------------------
# # Path to the pipeline output directory.
# pipeline_dir <- file.path(output_root, "origin_pipeline")
# 
# # Path to the main summary file.
# summary_path <- file.path(pipeline_dir, "gbif_summary.csv")
# 
# # Read and display the most useful summary fields if the file exists.
# if (file.exists(summary_path)) {
# 
#   origin_summary <- read.csv(summary_path, stringsAsFactors = FALSE)
# 
#   origin_summary[
#     ,
#     intersect(
#       c(
#         "species",
#         "region_id",
#         "region_type",
#         "n_total",
#         "n_cleaned",
#         "n_thinned",
#         "status",
#         "fail_stage",
#         "fail_reason",
#         "output_file"
#       ),
#       names(origin_summary)
#     ),
#     drop = FALSE
#   ]
# 
# } else {
#   skip_live_message("summary file was not available from the live pipeline run")
# }

## ----list-origin-audit-files-live, eval = run_live----------------------------
# # List audit and summary files written by the pipeline.
# # These files should be kept with the exported occurrence data.
# if (dir.exists(pipeline_dir)) {
#   list.files(
#     pipeline_dir,
#     pattern = "taxonomy|origin|native_web|griis|summary",
#     full.names = TRUE
#   )
# } else {
#   skip_live_message("pipeline output directory was not available")
# }

