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
# Set params$run_live to TRUE when rendering locally if you want to contact
# GBIF or download external spatial-provider resources.
run_live <- isTRUE(params$run_live)

# Read GBIF credentials from environment variables.
# This avoids writing usernames, passwords, or email addresses directly into
# vignettes, scripts, or shared project files.
gbif_user <- Sys.getenv("GBIF_USER")
gbif_pwd <- Sys.getenv("GBIF_PWD")
gbif_email <- Sys.getenv("GBIF_EMAIL")

# GBIF live examples are attempted only when all credentials are available.
gbif_ready <- all(nzchar(c(gbif_user, gbif_pwd, gbif_email)))

# Write vignette example outputs to a temporary folder.
# In a real project, replace tempdir() with a stable project directory.
output_root <- file.path(tempdir(), "biofetchR_spatial_overlays_vignette")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

# Run live examples without stopping the whole vignette if a provider request,
# spatial download, or GBIF request fails.
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
#   "vignettes/biofetchR-spatial-overlays.Rmd",
#   params = list(run_live = TRUE)
# )

## ----spatial-decision-table, echo = FALSE-------------------------------------
spatial_decision <- data.frame(
  analysis_goal = c(
    "Summarise records by administrative units",
    "Describe terrestrial ecological setting",
    "Describe freshwater ecological setting",
    "Add freshwater or aquatic-infrastructure context",
    "Add protected-area or wetland context",
    "Assign offshore records to marine regions",
    "Attach local environmental or human-pressure values"
  ),
  suitable_context = c(
    "GADM",
    "TEOW or RESOLVE 2017",
    "FEOW, basins, rivers, or lakes",
    "Rivers, lakes, basins, GDW barriers, HydroWASTE, or GDW reservoirs",
    "WDPA or Ramsar",
    "EEZ, LME, MEOW, or IHO",
    "Raster context"
  ),
  typical_output = c(
    "Administrative identifiers and names",
    "Terrestrial ecoregion, biome, or realm labels",
    "Freshwater ecoregion or hydrological identifiers",
    "Hydrological or infrastructure-context fields",
    "Protected-area or wetland-context fields",
    "Marine region identifiers and names",
    "Raster-derived covariate columns"
  ),
  stringsAsFactors = FALSE
)

spatial_decision

## ----example-occurrence-table-------------------------------------------------
example_occurrences <- data.frame(
  species = c(
    "Xenopus laevis",
    "Xenopus laevis",
    "Carcinus maenas"
  ),
  decimalLongitude = c(2.35, 1.44, -5.10),
  decimalLatitude = c(48.86, 43.60, 50.15),
  countryCode = c("FR", "FR", "GB"),
  year = c(2019, 2021, 2020),
  stringsAsFactors = FALSE
)

example_occurrences

## ----cleaning-options-table, echo = FALSE-------------------------------------
cleaning_options <- data.frame(
  option = c(
    "Coordinate screening",
    "Coordinate-uncertainty filtering",
    "High-uncertainty flagging",
    "Cleaning-only mode",
    "Distance-based thinning",
    "Density-based thinning",
    "Priority-based retention",
    "Diagnostic output"
  ),
  argument = c(
    "automatic",
    "filter_uncertain, max_uncertainty_m",
    "warn_uncertainty_m",
    "dist_km = 0",
    "dist_km > 0",
    "max_density",
    "priority_col",
    "return_all = TRUE"
  ),
  purpose = c(
    "Remove missing, non-finite, impossible, or exact zero-zero coordinates.",
    "Remove records with coordinate uncertainty above a chosen threshold.",
    "Flag retained records with high coordinate uncertainty.",
    "Apply coordinate cleaning without thinning.",
    "Reduce spatial clustering by retaining points separated by a minimum distance.",
    "Limit the approximate number of retained points per km².",
    "Preferentially retain higher-priority records during thinning.",
    "Return all post-cleaning rows with a `.thinned` flag showing which rows were retained."
  ),
  stringsAsFactors = FALSE
)

cleaning_options

## ----cleaning-example-input---------------------------------------------------
# Example occurrence table for demonstrating cleaning and thinning.
# This is deliberately small and includes several common coordinate problems.
cleaning_example <- data.frame(
  species = c(
    "Xenopus laevis",
    "Xenopus laevis",
    "Xenopus laevis",
    "Xenopus laevis",
    "Xenopus laevis",
    "Xenopus laevis",
    "Xenopus laevis"
  ),

  # Longitude values include valid coordinates, an impossible longitude,
  # an exact zero-zero point, and a missing value.
  decimalLongitude = c(
    2.3500,
    2.3505,
    2.3520,
    181.0000,
    0.0000,
    NA,
    1.4400
  ),

  # Latitude values include valid coordinates and an exact zero-zero point.
  decimalLatitude = c(
    48.8600,
    48.8605,
    48.8620,
    48.8600,
    0.0000,
    43.6000,
    43.6000
  ),

  # Coordinate uncertainty is optional in GBIF-style data.
  # Here, one otherwise valid record has high uncertainty.
  coordinateUncertaintyInMeters = c(
    30,
    50,
    12000,
    20,
    10,
    100,
    300
  ),

  # individualCount is used here as an example priority column.
  # Higher values are retained first during thinning when priority_col is used.
  individualCount = c(
    10,
    2,
    1,
    1,
    1,
    1,
    5
  ),

  year = c(2019, 2019, 2020, 2020, 2021, 2021, 2022),

  stringsAsFactors = FALSE
)

cleaning_example

## ----cleaning-only-example----------------------------------------------------
# Apply coordinate cleaning only.
# Records with missing, impossible, or exact zero-zero coordinates are removed.
# Records above max_uncertainty_m are also removed because filter_uncertain = TRUE.
cleaned_only <- thin_spatial_points(
  cleaning_example,
  dist_km = 0,
  filter_uncertain = TRUE,
  max_uncertainty_m = 10000,
  quiet = TRUE
)

# Drop geometry for compact display in the vignette.
sf::st_drop_geometry(cleaned_only)

## ----uncertainty-flag-example-------------------------------------------------
# Retain high-uncertainty records but flag them for review.
cleaned_with_uncertainty_flags <- thin_spatial_points(
  cleaning_example,
  dist_km = 0,
  filter_uncertain = FALSE,
  warn_uncertainty_m = 5000,
  quiet = TRUE
)

# Inspect the retained records and the uncertainty flag.
sf::st_drop_geometry(cleaned_with_uncertainty_flags)[
  ,
  intersect(
    c(
      "species",
      "decimalLongitude",
      "decimalLatitude",
      "coordinateUncertaintyInMeters",
      ".high_uncertainty"
    ),
    names(cleaned_with_uncertainty_flags)
  ),
  drop = FALSE
]

## ----distance-thinning-example------------------------------------------------
# Apply coordinate cleaning and then distance-based thinning.
# return_all = TRUE keeps all cleaned records and adds a `.thinned` diagnostic
# column showing which rows were retained by thinning.
distance_thinned <- thin_spatial_points(
  cleaning_example,
  dist_km = 5,
  priority_col = "individualCount",
  return_all = TRUE,
  filter_uncertain = FALSE,
  warn_uncertainty_m = 5000,
  quiet = TRUE
)

# Inspect thinning diagnostics.
sf::st_drop_geometry(distance_thinned)[
  ,
  intersect(
    c(
      "species",
      "decimalLongitude",
      "decimalLatitude",
      "individualCount",
      "coordinateUncertaintyInMeters",
      ".high_uncertainty",
      ".thinned"
    ),
    names(distance_thinned)
  ),
  drop = FALSE
]

## ----distance-thinning-example-overlay----------------------------------------
# Apply coordinate cleaning and then distance-based thinning.
# return_all = TRUE keeps all cleaned records and adds a `.thinned` diagnostic
# column showing which rows were retained by thinning.
distance_thinned <- thin_spatial_points(
  cleaning_example,
  dist_km = 5,
  priority_col = "individualCount",
  return_all = TRUE,
  filter_uncertain = FALSE,
  warn_uncertainty_m = 5000,
  quiet = TRUE
)

# Inspect thinning diagnostics.
sf::st_drop_geometry(distance_thinned)[
  ,
  intersect(
    c(
      "species",
      "decimalLongitude",
      "decimalLatitude",
      "individualCount",
      "coordinateUncertaintyInMeters",
      ".high_uncertainty",
      ".thinned"
    ),
    names(distance_thinned)
  ),
  drop = FALSE
]

## ----pipeline-cleaning-settings-example, eval = FALSE-------------------------
# process_gbif_terrestrial_freshwater_pipeline(
#   df = data.frame(
#     species = "Xenopus laevis",
#     iso2c = "FR",
#     stringsAsFactors = FALSE
#   ),
#   output_dir = file.path(output_root, "example_outputs"),
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   region_source = "gadm",
#   gadm_unit = 1,
#   cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#   # Apply coordinate quality control.
#   apply_cleaning = TRUE,
# 
#   # Apply spatial thinning after cleaning.
#   apply_thinning = TRUE,
# 
#   # Retain records separated by at least 5 km where possible.
#   dist_km = 5,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE
# )

## ----primary-secondary-table, echo = FALSE------------------------------------
primary_secondary <- data.frame(
  workflow = c(
    "Subnational terrestrial reporting",
    "Freshwater invasion analysis",
    "Protected-area screening",
    "Marine jurisdictional analysis",
    "Environmental-covariate analysis"
  ),
  primary_context = c(
    "GADM level 1",
    "GADM level 1 or basin unit",
    "GADM or ecological region",
    "EEZ",
    "GADM, ecoregion, basin, or marine region"
  ),
  secondary_context = c(
    "TEOW or RESOLVE 2017",
    "FEOW, rivers, lakes, dams, reservoirs",
    "WDPA or Ramsar",
    "LME, MEOW, IHO, or high seas",
    "Raster context"
  ),
  reason = c(
    "Administrative outputs can be interpreted alongside ecological regions.",
    "Reporting units and freshwater ecological units answer different questions.",
    "Protected-area context is usually an additional attribute, not the main grouping unit.",
    "Marine records need marine spatial units rather than land boundaries.",
    "Raster values describe point-level covariates rather than polygon membership."
  ),
  stringsAsFactors = FALSE
)

primary_secondary

## ----overlay-reference-table--------------------------------------------------
overlay_reference <- data.frame(
  overlay = c(
    "gadm",
    "teow",
    "feow",
    "lakes",
    "rivers",
    "basins",
    "ne_admin1",
    "resolve2017",
    "wdpa",
    "ramsar",
    "gdw_barriers",
    "hydrowaste",
    "gdw_reservoirs",
    "global_mining",
    "eez",
    "lme",
    "meow",
    "iho"
  ),
  typical_use = c(
    "Administrative attribution for terrestrial and freshwater records",
    "Terrestrial ecoregion context",
    "Freshwater ecoregion context",
    "Lake-associated occurrence context",
    "River-associated occurrence context",
    "Hydrological basin context",
    "Natural Earth administrative-unit context",
    "RESOLVE ecoregion and biome context",
    "Protected-area context",
    "Ramsar wetland context",
    "Barrier and dam context",
    "Wastewater-treatment context",
    "Reservoir context",
    "Mining or extraction-pressure context",
    "Marine jurisdictional attribution",
    "Large Marine Ecosystem attribution",
    "Marine ecoregion attribution",
    "Sea-area attribution"
  ),
  broad_domain = c(
    "terrestrial/freshwater",
    "terrestrial",
    "freshwater",
    "freshwater",
    "freshwater",
    "freshwater",
    "terrestrial/freshwater",
    "terrestrial",
    "terrestrial/freshwater/marine depending on source",
    "wetland/freshwater/coastal",
    "freshwater",
    "freshwater/infrastructure",
    "freshwater",
    "terrestrial/infrastructure",
    "marine",
    "marine",
    "marine",
    "marine"
  ),
  stringsAsFactors = FALSE
)

overlay_reference

## ----load-gadm-live, eval = run_live------------------------------------------
# # Load GADM level 1 polygons for Portugal.
# # This is useful for checking that GADM can be retrieved and cached before
# # running a larger GBIF workflow.
# gadm_prt_l1 <- run_live_safely(
#   load_gadm(
#     iso2c = "PT",
#     level = 1,
#     cache_dir = file.path(output_root, "cache", "gadm"),
#     quiet = FALSE
#   )
# )
# 
# # Inspect the returned GADM object if loading succeeded.
# if (!is.null(gadm_prt_l1)) {
#   gadm_prt_l1
# }

## ----gadm-pipeline-live, eval = run_live--------------------------------------
# # Run the live GADM pipeline example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   # Skip cleanly rather than failing the vignette.
#   skip_live_message("set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run this example")
# 
# } else {
# 
#   # Run a minimal terrestrial/freshwater pipeline with GADM attribution.
#   gadm_pipeline_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
# 
#       # Species-country input table.
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "FR",
#         stringsAsFactors = FALSE
#       ),
# 
#       # Output folder for this example.
#       output_dir = file.path(output_root, "gadm_pipeline"),
# 
#       # GBIF credentials read from environment variables.
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Use GADM level 1 as the primary spatial framework.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Clean coordinates and thin retained records.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Keep the example small and write outputs to disk.
#       batch_size = 1,
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
# 
#       # Show progress messages.
#       quiet = FALSE
#     )
#   )
# 
#   gadm_pipeline_result
# }

## ----inspect-gadm-summary-live, eval = run_live-------------------------------
# # Path to the summary file written by the GADM example.
# summary_path <- file.path(
#   output_root,
#   "gadm_pipeline",
#   "gbif_summary.csv"
# )
# 
# # Read the summary if it exists.
# if (file.exists(summary_path)) {
#   gadm_summary <- read.csv(summary_path, stringsAsFactors = FALSE)
#   gadm_summary
# } else {
#   skip_live_message("summary file was not available from the GADM live example")
# }

## ----feow-overlay-pipeline-live, eval = run_live------------------------------
# # Run the live GADM + FEOW example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   skip_live_message("set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run this example")
# 
# } else {
# 
#   # Run a pipeline with both administrative and freshwater ecological context.
#   feow_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "FR",
#         stringsAsFactors = FALSE
#       ),
# 
#       output_dir = file.path(output_root, "feow_pipeline"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Use GADM and FEOW together.
#       region_source = c("gadm", "feow"),
# 
#       # Keep the two spatial contexts as separate output fields.
#       overlay_mode = "dual_separate",
# 
#       # GADM and FEOW cache directories.
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
#       feow_cache_dir = file.path(output_root, "cache", "feow"),
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       batch_size = 1,
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   )
# 
#   feow_result
# }

## ----provider-loaders-live, eval = run_live-----------------------------------
# # Load Natural Earth administrative level 1 context for Great Britain and Ireland.
# ne_admin1 <- run_live_safely(
#   bf_load_ne_admin1(
#     iso2c = c("GB", "IE"),
#     cache_dir = file.path(output_root, "cache", "ne_admin1"),
#     quiet = FALSE
#   )
# )
# 
# # Load a RESOLVE 2017 terrestrial ecoregion subset.
# resolve2017 <- run_live_safely(
#   bf_load_resolve2017(
#     iso2c = c("GB", "IE"),
#     cache_dir = file.path(output_root, "cache", "resolve2017"),
#     quiet = FALSE
#   )
# )
# 
# if (!is.null(ne_admin1)) {
#   ne_admin1
# }
# 
# if (!is.null(resolve2017)) {
#   resolve2017
# }

## ----custom-overlay-template, eval = FALSE------------------------------------
# # Example only: start from an existing occurrence table with coordinates.
# existing_occurrences <- example_occurrences
# 
# # Convert the occurrence table to an sf point object.
# occ_sf <- sf::st_as_sf(
#   existing_occurrences,
#   coords = c("decimalLongitude", "decimalLatitude"),
#   crs = 4326,
#   remove = FALSE
# )
# 
# # Load a provider layer.
# resolve2017 <- bf_load_resolve2017(
#   iso2c = c("GB", "FR"),
#   cache_dir = file.path(output_root, "cache", "resolve2017"),
#   quiet = FALSE
# )
# 
# # Join points to polygons.
# # In a real workflow, inspect unmatched points and check CRS/geometry validity.
# occ_with_resolve <- sf::st_join(
#   occ_sf,
#   resolve2017,
#   left = TRUE
# )
# 
# occ_with_resolve

## ----raster-context-live, eval = run_live-------------------------------------
# # Run the raster-context example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   skip_live_message("set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run this example")
# 
# } else {
# 
#   # Run a small pipeline with GADM attribution and raster context.
#   raster_context_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "FR",
#         stringsAsFactors = FALSE
#       ),
# 
#       output_dir = file.path(output_root, "raster_context"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Keep GADM as the primary vector framework.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Add gridded environmental or human-pressure context.
#       raster_context = c("worldcover", "soilgrids", "human_footprint"),
#       raster_cache_dir = file.path(output_root, "cache", "raster"),
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       batch_size = 1,
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   )
# 
#   raster_context_result
# }

## ----raster-context-live-example, eval = run_live-----------------------------
# # Run the raster-context example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   skip_live_message("set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run this example")
# 
# } else {
# 
#   # Run a small pipeline with GADM attribution and raster context.
#   raster_context_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "FR",
#         stringsAsFactors = FALSE
#       ),
# 
#       output_dir = file.path(output_root, "raster_context"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Keep GADM as the primary vector framework.
#       region_source = "gadm",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # Add gridded environmental or human-pressure context.
#       raster_context = c("worldcover", "soilgrids", "human_footprint"),
#       raster_cache_dir = file.path(output_root, "cache", "raster"),
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       batch_size = 1,
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   )
# 
#   raster_context_result
# }

## ----protected-wetland-loaders-live, eval = run_live--------------------------
# # WDPA and Ramsar loading is opt-in because these datasets have
# # external licensing and attribution requirements.
# #
# # To run this chunk locally, set the option below in the console before
# # rendering the vignette:
# #
# # options(biofetchR.wdpa_opt_in = TRUE)
# 
# protected_overlay_opt_in <- isTRUE(
#   getOption("biofetchR.wdpa_opt_in", FALSE)
# )
# 
# if (!protected_overlay_opt_in) {
# 
#   skip_live_message(
#     "set options(biofetchR.wdpa_opt_in = TRUE) after checking provider terms"
#   )
# 
# } else {
# 
#   # Load protected-area polygons for Great Britain.
#   wdpa_gb <- run_live_safely(
#     bf_load_wdpa(
#       iso2c = "GB",
#       cache_dir = file.path(output_root, "cache", "wdpa"),
#       opt_in = TRUE,
#       quiet = FALSE
#     )
#   )
# 
#   # Load Ramsar wetland polygons for Great Britain.
#   ramsar_gb <- run_live_safely(
#     bf_load_ramsar(
#       iso2c = "GB",
#       cache_dir = file.path(output_root, "cache", "ramsar"),
#       opt_in = TRUE,
#       quiet = FALSE
#     )
#   )
# 
#   list(
#     wdpa = wdpa_gb,
#     ramsar = ramsar_gb
#   )
# }

## ----protected-overlay-pipeline-live, eval = run_live-------------------------
# # Run the live protected-area overlay example only when GBIF credentials
# # are available and WDPA opt-in has been confirmed.
# if (!gbif_ready) {
# 
#   skip_live_message("set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run this example")
# 
# } else if (!isTRUE(getOption("biofetchR.wdpa_opt_in", FALSE))) {
# 
#   skip_live_message(
#     "set options(biofetchR.wdpa_opt_in = TRUE) after checking provider terms"
#   )
# 
# } else {
# 
#   protected_overlay_result <- run_live_safely(
#     process_gbif_terrestrial_freshwater_pipeline(
#       df = data.frame(
#         species = "Xenopus laevis",
#         iso2c = "GB",
#         stringsAsFactors = FALSE
#       ),
# 
#       output_dir = file.path(output_root, "protected_overlay_pipeline"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Keep GADM as the main administrative framework and add WDPA.
#       region_source = c("gadm", "wdpa"),
#       overlay_mode = "dual_separate",
#       gadm_unit = 1,
#       cache_dir_gadm = file.path(output_root, "cache", "gadm"),
# 
#       # WDPA is opt-in because it is an externally licensed dataset.
#       wdpa_cache_dir = file.path(output_root, "cache", "wdpa"),
#       wdpa_opt_in = TRUE,
#       wdpa_exclude_marine = TRUE,
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       batch_size = 1,
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       quiet = FALSE
#     )
#   )
# 
#   protected_overlay_result
# }

## ----spatial-summary-example--------------------------------------------------
spatial_summary_example <- data.frame(
  species = c(
    "Xenopus laevis",
    "Xenopus laevis",
    "Carcinus maenas"
  ),
  workflow = c(
    "gadm",
    "gadm_plus_feow",
    "marine_eez"
  ),
  primary_spatial_unit = c(
    "GADM level 1",
    "GADM level 1",
    "EEZ"
  ),
  secondary_spatial_context = c(
    NA,
    "FEOW",
    NA
  ),
  n_total = c(1340, 1340, 980),
  n_cleaned = c(1309, 1309, 944),
  n_thinned = c(861, 861, 702),
  status = c("success", "success", "success"),
  output_file = c(
    "Xenopus_laevis_FR_gadm1.csv",
    "Xenopus_laevis_FR_gadm1_feow.csv",
    "Carcinus_maenas_eez.csv"
  ),
  stringsAsFactors = FALSE
)

spatial_summary_example

## ----inspect-all-summary-template, eval = FALSE-------------------------------
# summary_path <- file.path("path", "to", "gbif_summary.csv")
# 
# summary_tbl <- read.csv(summary_path, stringsAsFactors = FALSE)
# 
# summary_tbl[
#   ,
#   intersect(
#     c(
#       "species",
#       "region_id",
#       "region_type",
#       "n_total",
#       "n_cleaned",
#       "n_thinned",
#       "status",
#       "fail_stage",
#       "fail_reason",
#       "output_file"
#     ),
#     names(summary_tbl)
#   ),
#   drop = FALSE
# ]

## ----troubleshooting-table, echo = FALSE--------------------------------------
troubleshooting <- data.frame(
  symptom = c(
    "No records after cleaning",
    "Many records have no spatial unit",
    "Provider layer fails to load",
    "Raster context fails",
    "Marine records are missing",
    "Output files are not where expected",
    "Very large memory use"
  ),
  possible_cause = c(
    "Coordinate cleaning removed records with invalid, missing, or suspicious coordinates.",
    "Coordinates may fall outside the overlay extent, offshore, near boundaries, or in the wrong CRS.",
    "Provider service, URL, dependency, cache, or licence/opt-in requirement may be unavailable.",
    "Raster file may be large, unavailable, in an unexpected format, or incompatible with local memory.",
    "Marine species may have been submitted through a terrestrial/country-filtered workflow.",
    "The pipeline may have written split outputs by spatial unit or failed before export.",
    "Large GBIF downloads, high-resolution overlays, raster extraction, or in-memory storage may be too heavy."
  ),
  first_check = c(
    "Inspect n_total, n_cleaned, n_thinned, fail_stage, and fail_reason in gbif_summary.csv.",
    "Check coordinates, overlay choice, CRS, coastal/offshore points, and unmatched rows.",
    "Try the provider loader directly with a small country or overlay subset.",
    "Test one species-region pair and inspect raster_cache_dir.",
    "Use process_gbif_marine_pipeline() and a marine overlay such as EEZ, LME, MEOW, or IHO.",
    "Inspect output_dir, output_file, status, and fail_reason in gbif_summary.csv.",
    "Set store_in_memory = FALSE and add overlays or raster context in stages."
  ),
  stringsAsFactors = FALSE
)

troubleshooting

