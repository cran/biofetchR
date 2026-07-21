params <-
list(run_live = FALSE)

## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = TRUE,
  warning = TRUE
)

# Read the vignette parameter that controls whether live examples are run.
# When run_live is FALSE, live GBIF and provider-download examples are skipped.
# When run_live is TRUE, those examples are executed during rendering.
run_live <- isTRUE(params$run_live)

# Read GBIF credentials from environment variables.
# This avoids writing usernames, passwords, or email addresses directly into
# scripts, vignettes, or shared project files.
gbif_user <- Sys.getenv("GBIF_USER")
gbif_pwd <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")

# Check whether all three GBIF credential fields are available.
# Live GBIF examples are only attempted when this is TRUE.
gbif_ready <- all(nzchar(c(gbif_user, gbif_pwd, gbif_email)))

# Run a live example without stopping the whole vignette if an external service,
# provider download, or GBIF request fails. The error message is printed so the
# user can see what happened, and NULL is returned for that example.
run_live_safely <- function(expr) {
  tryCatch(
    force(expr),
    error = function(e) {
      message("Live example did not complete: ", conditionMessage(e))
      NULL
    }
  )
}

# Print a clear message when a live example is skipped because credentials or
# another required setup step are missing.
skip_live_message <- function(reason) {
  message("Live example skipped: ", reason)
  invisible(NULL)
}

## ----load-package-------------------------------------------------------------
library(biofetchR)

## ----render-live-version, eval = FALSE----------------------------------------
# rmarkdown::render(
#   "vignettes/biofetchR-batch-pipelines.Rmd",
#   params = list(run_live = TRUE)
# )

## ----terrestrial-batch-input--------------------------------------------------
terrestrial_batch <- data.frame(
  species = c(
    "Xenopus laevis",
    "Trachemys scripta",
    "Rattus rattus",
    "Sturnus vulgaris"
  ),
  iso2c = c(
    "FR",
    "GB",
    "GB",
    "DE"
  ),
  stringsAsFactors = FALSE
)

terrestrial_batch

## ----terrestrial-batch-input-with-metadata------------------------------------
terrestrial_batch_metadata <- data.frame(
  species = c(
    "Xenopus laevis",
    "Trachemys scripta"
  ),
  iso2c = c(
    "FR",
    "GB"
  ),
  taxon_group = c(
    "Amphibian",
    "Reptile"
  ),
  notes = c(
    "low-memory smoke test",
    "small batch example"
  ),
  stringsAsFactors = FALSE
)

terrestrial_batch_metadata

## ----marine-batch-input-------------------------------------------------------
marine_batch <- data.frame(
  species = c(
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida"
  ),
  stringsAsFactors = FALSE
)

marine_batch

## ----credentials-and-output---------------------------------------------------
# Create a temporary output folder for this vignette run.
# In a real project, replace tempdir() with a stable project directory,
# such as "outputs/gbif_batches".
output_root <- file.path(tempdir(), "biofetchR_batch_vignette")

# Create the output folder if it does not already exist.
# recursive = TRUE also creates any missing parent folders.
# showWarnings = FALSE avoids unnecessary messages if the folder already exists.
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Create a small credential-check table.
# This does not print the credential values themselves; it only reports whether
# each required environment variable is available.
credential_status <- data.frame(
  variable = c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"),
  available = nzchar(c(gbif_user, gbif_pwd, gbif_email)),
  stringsAsFactors = FALSE
)

# Display the credential check and the output folder being used for this run.
credential_status
output_root

## ----cleaning-thinning-template, eval = FALSE---------------------------------
# process_gbif_terrestrial_freshwater_pipeline(
#   df = terrestrial_batch,
# 
#   output_dir = file.path(output_root, "cleaned_thinned_batch"),
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   region_source = "gadm",
#   gadm_unit = 1,
# 
#   apply_cleaning = TRUE,
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE
# )

## ----low-memory-terrestrial-live, eval = run_live-----------------------------
# if (!gbif_ready) {
#   skip_live_message(
#     "set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run the GBIF smoke test"
#   )
# } else {
# 
#   # Use one species-country pair for the first test run.
#   # This is enough to check credentials, GBIF access, spatial attribution,
#   # cleaning, thinning, evidence attachment, and output writing.
#   low_memory_input <- data.frame(
#     species = "Xenopus laevis",
#     iso2c = "FR",
#     stringsAsFactors = FALSE
#   )
# 
#   low_memory_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = low_memory_input,
# 
#       # Write all example outputs into the vignette output folder.
#       output_dir = file.path(output_root, "terrestrial_low_memory"),
# 
#       # GBIF credentials are read from environment variables above.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Assign records to GADM level-1 administrative units.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Keep the first run deliberately small and conservative.
#       batch_size = 1,
# 
#       # Apply coordinate quality control and then thin clustered records.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Prepare taxon names before submitting the GBIF download.
#       prepare_taxonomy = TRUE,
# 
#       # Attach GRIIS evidence in audit mode.
#       # Audit mode keeps the records and adds evidence fields for inspection.
#       use_griis_status = TRUE,
#       griis_filter_mode = "audit_only",
#       griis_cache_dir = file.path(output_root, "cache", "griis"),
# 
#       # Attach native-range evidence from selected web sources.
#       # For marine workflows, WoRMS can be useful alongside SInAS and GBIF.
#       use_native_status = TRUE,
#       use_native_web = TRUE,
#       native_web_sources = c("sinas", "gbif", "worms"),
#       native_web_cache_dir = file.path(output_root, "cache", "native_web"),
#       native_filter_mode = "audit_only",
# 
#       # Write origin-evidence and processing summaries.
#       reconcile_origin_evidence = TRUE,
#       export_origin_audit = TRUE,
#       export_summary = TRUE,
# 
#       # Write occurrence outputs to disk, but also return a result object
#       # so the smoke test can be inspected immediately in the R session.
#       return_all_results = TRUE,
#       store_in_memory = FALSE,
# 
#       # Show progress messages during the live run.
#       quiet = FALSE
#     )
#   )
# 
#   low_memory_result
# }

## ----low-memory-summary-live, eval = run_live---------------------------------
# # Build the path to the summary file written by the smoke-test run.
# # The summary is stored inside the output directory used for this example.
# summary_path <- file.path(
#   output_root,
#   "terrestrial_low_memory",
#   "gbif_summary.csv"
# )
# 
# # Read and inspect the summary only if the file was created.
# # This avoids an error if the previous live example was skipped or did not finish.
# if (file.exists(summary_path)) {
# 
#   # Import the summary table as a plain data frame.
#   low_memory_summary <- read.csv(summary_path, stringsAsFactors = FALSE)
# 
#   # Display the most useful audit columns.
#   # intersect() is used so the code still works if some optional columns are not
#   # present in a particular pipeline run or package version.
#   low_memory_summary[
#     ,
#     intersect(
#       c(
#         "species",     # taxon processed
#         "region_id",   # spatial unit or recipient-region identifier
#         "region_type", # spatial framework used, such as GADM
#         "n_total",     # records before cleaning
#         "n_cleaned",   # records retained after coordinate cleaning
#         "n_thinned",   # records retained after spatial thinning
#         "status",      # success, failure, or no-record status
#         "fail_stage",  # pipeline stage where a failure occurred, if any
#         "fail_reason", # short explanation of the failure, if any
#         "output_file"  # exported occurrence file written for the row
#       ),
#       names(low_memory_summary)
#     )
#   ]
# } else {
# 
#   # This message is shown when the live example was skipped or no summary was
#   # written. It keeps the vignette readable without stopping execution.
#   skip_live_message("summary file was not created by the previous live chunk")
# }

## ----synthetic-terrestrial-summary--------------------------------------------
# In a real workflow, this information is read from gbif_summary.csv.
summary_example <- data.frame(

  # Taxon processed by the pipeline.
  species = c("Xenopus laevis", "Trachemys scripta", "Rattus rattus"),

  # Recipient region or spatial unit associated with each row.
  # For this terrestrial example, these are ISO2 country codes.
  region_id = c("FR", "GB", "GB"),

  # Spatial framework used for the exported records.
  # Here, records are assumed to have been assigned to GADM level-1 units.
  region_type = c("GADM_1", "GADM_1", "GADM_1"),

  # Number of records returned before coordinate cleaning.
  n_total = c(1340, 220, 0),

  # Number of records retained after coordinate quality control.
  n_cleaned = c(1309, 198, 0),

  # Number of records retained after spatial thinning.
  n_thinned = c(861, 174, 0),

  # Processing outcome for each species-region row.
  status = c("success", "success", "no_records"),

  # Pipeline stage where a row failed or stopped, if applicable.
  fail_stage = c(NA, NA, "import"),

  # Short explanation of the failure or stopping reason, if applicable.
  fail_reason = c(NA, NA, "No records returned"),

  # Path to the exported occurrence file.
  # Rows with no usable records may not produce an output file.
  output_file = c(
    "Xenopus_laevis_FR_gadm1.csv",
    "Trachemys_scripta_GB_gadm1.csv",
    NA
  ),

  stringsAsFactors = FALSE
)

# Display the example summary table.
summary_example

## ----inspect-synthetic-terrestrial-summary------------------------------------
table(summary_example$status, useNA = "ifany")

summary_example[summary_example$status != "success", ]

## ----terrestrial-batch-live, eval = run_live----------------------------------
# if (!gbif_ready) {
#   skip_live_message("set GBIF credentials to run the terrestrial batch example")
# } else {
# 
#   # Run the same terrestrial/freshwater workflow on a slightly larger
#   # species-country table. This is the next step after the one-row smoke test.
#   terrestrial_batch_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = terrestrial_batch,
# 
#       # Write this batch into its own output folder so it is easy to inspect
#       # separately from the smoke-test output.
#       output_dir = file.path(output_root, "terrestrial_batch"),
# 
#       # GBIF credentials are read from environment variables defined earlier.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Assign records to GADM level-1 administrative units.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Keep batches small in this example.
#       # For larger projects, increase batch_size gradually after checking
#       # that the first few rows complete successfully.
#       batch_size = 1,
# 
#       # Apply coordinate cleaning and 5 km spatial thinning.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Standardise and resolve taxon names before GBIF download submission.
#       prepare_taxonomy = TRUE,
# 
#       # Attach GRIIS evidence without filtering rows out of this example.
#       # This keeps the evidence available for inspection in the audit outputs.
#       use_griis = TRUE,
#       filter_griis_invasive = FALSE,
#       griis_cache_dir = file.path(output_root, "cache", "griis"),
# 
#       # Attach native-range evidence from selected web sources.
#       # audit_only keeps all records and writes evidence for later checking.
#       use_native_range = TRUE,
#       use_native_web = TRUE,
#       native_web_sources = c("sinas", "gbif"),
#       native_web_cache_dir = file.path(output_root, "cache", "native_web"),
#       native_filter_mode = "audit_only",
# 
#       # Reconcile GRIIS and native-origin evidence and write audit outputs.
#       reconcile_origin_evidence = TRUE,
#       export_origin_audit = TRUE,
# 
#       # Write gbif_summary.csv for the batch.
#       export_summary = TRUE,
# 
#       # For batch workflows, write occurrence tables to disk rather than
#       # keeping all processed records in memory.
#       store_in_memory = FALSE,
# 
#       # Return only the main run object because outputs are written to disk.
#       return_all_results = FALSE,
# 
#       # Print progress messages during the live run.
#       quiet = FALSE
#     )
#   )
# 
#   # Display the returned run object.
#   terrestrial_batch_result
# }

## ----terrestrial-overlay-live, eval = run_live--------------------------------
# if (!gbif_ready) {
#   skip_live_message("set GBIF credentials to run the overlay example")
# } else {
# 
#   # Run a small overlay example using one species-country pair.
#   # This keeps the example manageable while showing how to combine
#   # administrative and ecological spatial attribution.
#   feow_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
# 
#       # Input table for the example.
#       # The species is processed in France as the recipient country.
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "FR",
#         stringsAsFactors = FALSE
#       ),
# 
#       # Write the overlay example into its own output folder.
#       output_dir = file.path(output_root, "terrestrial_feow"),
# 
#       # GBIF credentials are read from environment variables defined earlier.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Use both GADM and FEOW as spatial frameworks.
#       # GADM provides administrative attribution.
#       # FEOW provides freshwater ecoregion attribution.
#       region_source = c("gadm", "feow"),
# 
#       # Keep the GADM and FEOW joins as separate outputs.
#       # This is useful when the two spatial frameworks should be inspected
#       # independently rather than merged into one combined unit.
#       overlay_mode = "dual_separate",
# 
#       # Use GADM level 1 for subnational administrative attribution.
#       gadm_unit = 1,
# 
#       # Cache downloaded or loaded spatial layers so they can be reused.
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
#       feow_cache_dir = file.path(output_root, "cache", "feow"),
# 
#       # Keep the example to one GBIF request.
#       batch_size = 1,
# 
#       # Apply coordinate cleaning and 5 km spatial thinning before export.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Write gbif_summary.csv for the overlay run.
#       export_summary = TRUE,
# 
#       # Write processed occurrence tables to disk rather than storing
#       # all records in the active R session.
#       store_in_memory = FALSE,
# 
#       # Return only the main run object because detailed outputs are written
#       # to the output directory.
#       return_all_results = FALSE,
# 
#       # Show progress messages during the live run.
#       quiet = FALSE
#     )
#   )
# 
#   # Display the returned run object.
#   feow_result
# }

## ----marine-batch-live, eval = run_live---------------------------------------
# if (!gbif_ready) {
#   skip_live_message("set GBIF credentials to run the marine batch example")
# } else {
# 
#   # Run a small marine batch example.
#   # Unlike the terrestrial workflow, the marine pipeline starts from a species
#   # list and assigns records to marine regions after the GBIF records are imported.
#   marine_batch_result <- run_live_safely(
#     process_gbif_marine_pipeline(
#       df = marine_batch,
# 
#       # Write the marine batch outputs into their own folder.
#       output_dir = file.path(output_root, "marine_batch"),
# 
#       # GBIF credentials are read from environment variables defined earlier.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Assign imported occurrence records to Exclusive Economic Zones.
#       # Other marine overlays can be used in the dedicated marine workflow.
#       overlay = "eez",
# 
#       # Cache the marine overlay so it can be reused in later runs.
#       overlay_cache_dir = file.path(output_root, "cache", "marine_eez"),
# 
#       # Allow the example to continue more flexibly if a remote overlay source
#       # is temporarily unavailable or needs to be resolved from cache.
#       strict_overlay_loading = FALSE,
# 
#       # Keep the example small by processing one species at a time.
#       batch_size = 1,
# 
#       # Apply coordinate cleaning and 5 km spatial thinning before export.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Standardise and resolve taxon names before submitting GBIF downloads.
#       prepare_taxonomy = TRUE,
# 
#       # Attach GRIIS evidence in audit mode.
#       # Audit mode keeps the records and adds evidence fields for inspection.
#       use_griis_status = TRUE,
#       griis_filter_mode = "audit_only",
# 
#       # Attach native-range evidence from selected web sources.
#       # For marine workflows, WoRMS can be useful alongside SInAS and GBIF.
#       use_native_status = TRUE,
#       use_native_web = TRUE,
#       native_web_sources = c("sinas", "gbif", "worms"),
#       native_filter_mode = "audit_only",
# 
#       # Compare and reconcile GRIIS and native-origin evidence where available.
#       reconcile_griis_native = TRUE,
# 
#       # Write gbif_summary.csv for the marine batch.
#       export_summary = TRUE,
# 
#       # Write processed occurrence tables to disk rather than storing all
#       # records in the active R session.
#       store_in_memory = FALSE,
# 
#       # Return only the main run object because detailed outputs are written
#       # to the output directory.
#       return_all_results = FALSE,
# 
#       # Show progress messages during the live run.
#       quiet = FALSE
#     )
#   )
# 
#   # Display the returned run object.
#   marine_batch_result
# }

## ----output-audit-example-----------------------------------------------------
# Create a compact audit view from the summary table.
# This keeps the columns most useful for checking whether each row completed
# successfully and where the exported occurrence file was written.
audit_view <- summary_example[
  ,
  intersect(
    c(
      "species",     # taxon processed
      "region_id",   # recipient region or spatial unit
      "region_type", # spatial framework used
      "status",      # processing outcome for the row
      "fail_stage",  # stage where the row failed or stopped, if applicable
      "fail_reason", # explanation of the failure or stopping reason
      "output_file"  # exported occurrence file path, if one was written
    ),
    names(summary_example)
  )
]

# Display the audit view.
audit_view

