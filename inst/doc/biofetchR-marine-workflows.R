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
# GBIF or download external marine-region resources.
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
output_root <- file.path(tempdir(), "biofetchR_marine_workflows_vignette")
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

# Resolve output paths listed in summary files.
# Some workflows may write absolute paths, while others may write paths relative
# to the output directory.
resolve_output_file <- function(path, output_dir) {
  if (is.na(path) || !nzchar(path)) {
    return(NA_character_)
  }

  if (file.exists(path)) {
    return(path)
  }

  candidate <- file.path(output_dir, path)

  if (file.exists(candidate)) {
    return(candidate)
  }

  path
}

## ----load-package-------------------------------------------------------------
library(biofetchR)

## ----render-live-version, eval = FALSE----------------------------------------
# rmarkdown::render(
#   "vignettes/biofetchR-marine-workflows.Rmd",
#   params = list(run_live = TRUE)
# )

## ----marine-input-table-------------------------------------------------------
marine_species <- data.frame(
  species = c(
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida",
    "Didemnum vexillum"
  ),
  stringsAsFactors = FALSE
)

marine_species

## ----marine-input-with-metadata-----------------------------------------------
marine_species_with_metadata <- data.frame(
  species = c(
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida"
  ),
  broad_group = c(
    "crustacean",
    "ctenophore",
    "macroalga"
  ),
  expected_habitat = c(
    "coastal_benthic",
    "pelagic_coastal",
    "coastal_benthic"
  ),
  project_id = c(
    "marine_001",
    "marine_002",
    "marine_003"
  ),
  stringsAsFactors = FALSE
)

marine_species_with_metadata

## ----marine-overlay-decision-table, echo = FALSE------------------------------
marine_overlay_decision <- data.frame(
  overlay = c(
    "EEZ",
    "LME",
    "MEOW",
    "IHO"
  ),
  represents = c(
    "Exclusive Economic Zones",
    "Large Marine Ecosystems",
    "Marine Ecoregions of the World",
    "International Hydrographic Organization sea areas"
  ),
  use_when = c(
    "The analysis needs country-associated marine jurisdiction or offshore management areas.",
    "The analysis needs broad marine ecosystem regions for large-scale ecological summaries.",
    "The analysis focuses on coastal marine biogeography or nearshore ecoregional patterns.",
    "The analysis needs named seas, basins, or hydrographic sea-area labels."
  ),
  typical_question = c(
    "Which country-associated marine area contains these records?",
    "Which broad marine ecosystem contains these records?",
    "Which coastal marine ecoregion contains these records?",
    "Which named sea area contains these records?"
  ),
  stringsAsFactors = FALSE
)

marine_overlay_decision

## ----available-marine-overlays------------------------------------------------
# Core marine overlays.
bf_available_marine_overlays()

## ----available-marine-overlays-aliases----------------------------------------
# Core overlays plus aliases.
bf_available_marine_overlays(include_aliases = TRUE)

## ----available-extended-marine-overlays---------------------------------------
bf_available_marine_overlays(include_extended = TRUE)

## ----extended-marine-overlay-table, echo = FALSE------------------------------
extended_marine_overlays <- data.frame(
  overlay = c(
    "high_seas",
    "territorial_seas",
    "contiguous_zone",
    "internal_waters",
    "archipelagic_waters",
    "ecs",
    "goas",
    "eez_iho",
    "longhurst"
  ),
  spatial_meaning = c(
    "Marine areas beyond national jurisdiction",
    "Coastal waters close to national coastlines",
    "Zones beyond territorial seas but still linked to coastal-state jurisdiction",
    "Waters landward of coastal baselines",
    "Waters enclosed by archipelagic baselines",
    "Extended continental shelves",
    "Global Oceans and Seas regions",
    "Intersection of EEZ and IHO sea-area layers",
    "Biogeochemical ocean provinces"
  ),
  useful_when = c(
    "The analysis needs to identify records outside EEZ-style national marine areas.",
    "The analysis focuses on nearshore or coastal-state waters.",
    "The analysis needs a broader coastal jurisdictional zone than territorial seas.",
    "The analysis includes records in enclosed coastal or port-associated waters.",
    "The analysis includes island or archipelagic marine regions.",
    "The analysis needs continental shelf context for offshore records.",
    "The analysis needs broad named ocean or sea regions.",
    "The analysis needs both jurisdictional and named sea-area context together.",
    "The analysis needs broad oceanographic or productivity-region context."
  ),
  stringsAsFactors = FALSE
)

extended_marine_overlays

## ----complementary-marine-overlay-template, eval = FALSE----------------------
# # Example overlays for a project interested in complementary offshore coverage.
# # high_seas is useful for open-ocean areas beyond national jurisdiction.
# # ecs is useful for extended continental shelf context.
# marine_overlays_to_run <- c(
#   "high_seas",
#   "ecs"
# )
# 
# complementary_overlay_results <- lapply(
#   marine_overlays_to_run,
#   function(this_overlay) {
# 
#     process_gbif_marine_pipeline(
#       df = marine_species,
# 
#       output_dir = file.path(
#         output_root,
#         paste0("marine_", this_overlay)
#       ),
# 
#       user = Sys.getenv("GBIF_USER"),
#       pwd = Sys.getenv("GBIF_PWD"),
#       email = Sys.getenv("GBIF_EMAIL"),
# 
#       overlay = this_overlay,
#       overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
# 
#       batch_size = 1,
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       prepare_taxonomy = TRUE,
# 
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   }
# )
# 
# names(complementary_overlay_results) <- marine_overlays_to_run

## ----load-eez-overlay-live, eval = run_live-----------------------------------
# # Load an EEZ overlay from the package-managed marine-region workflow.
# # This tests provider access, cache writing, and local overlay loading before
# # any GBIF records are downloaded.
# eez_overlay <- run_live_safely(
#   bf_load_marine_regions_overlay(
#     region_source = "eez",
#     cache_dir = file.path(output_root, "cache", "marine_overlays"),
#     quiet = FALSE
#   )
# )
# 
# # Inspect basic overlay diagnostics if loading succeeded.
# if (!is.null(eez_overlay)) {
# 
#   data.frame(
#     overlay = "eez",
#     n_features = nrow(eez_overlay),
#     n_columns = ncol(eez_overlay),
#     geometry_type = paste(
#       unique(as.character(sf::st_geometry_type(eez_overlay))),
#       collapse = ", "
#     ),
#     crs = as.character(sf::st_crs(eez_overlay)$input),
#     stringsAsFactors = FALSE
#   )
# 
# }

## ----inspect-eez-overlay-fields-live, eval = run_live-------------------------
# if (!is.null(eez_overlay)) {
#   names(eez_overlay)
# }

## ----low-memory-marine-live, eval = run_live----------------------------------
# # Run the marine smoke test only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   skip_live_message(
#     "set GBIF_USER, GBIF_PWD, and GBIF_EMAIL to run the marine smoke test"
#   )
# 
# } else {
# 
#   # Keep the smoke-test input deliberately small.
#   low_memory_marine <- data.frame(
#     species = "Carcinus maenas",
#     stringsAsFactors = FALSE
#   )
# 
#   low_memory_marine_result <- run_live_safely(
#     process_gbif_marine_pipeline(
#       df = low_memory_marine,
# 
#       output_dir = file.path(output_root, "marine_low_memory"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Assign retained records to EEZ polygons.
#       overlay = "eez",
#       overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
# 
#       # Allow the vignette example to fail softly if the provider is unavailable.
#       strict_overlay_loading = FALSE,
# 
#       # Keep memory use low.
#       batch_size = 1,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
# 
#       # Apply basic occurrence-quality controls.
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       # Prepare taxonomy before export.
#       prepare_taxonomy = TRUE,
# 
#       export_summary = TRUE,
#       quiet = FALSE
#     )
#   )
# 
#   low_memory_marine_result
# }

## ----inspect-low-memory-marine-summary-live, eval = run_live------------------
# marine_summary_path <- file.path(
#   output_root,
#   "marine_low_memory",
#   "gbif_summary.csv"
# )
# 
# if (file.exists(marine_summary_path)) {
#   marine_summary <- read.csv(marine_summary_path, stringsAsFactors = FALSE)
#   marine_summary
# } else {
#   skip_live_message("marine summary file was not created by the live smoke test")
# }

## ----marine-summary-example---------------------------------------------------
marine_summary_example <- data.frame(
  species = c(
    "Carcinus maenas",
    "Mnemiopsis leidyi",
    "Undaria pinnatifida"
  ),
  overlay = c(
    "eez",
    "eez",
    "eez"
  ),
  n_total = c(
    980,
    640,
    410
  ),
  n_cleaned = c(
    944,
    612,
    389
  ),
  n_thinned = c(
    702,
    455,
    301
  ),
  n_with_overlay = c(
    690,
    430,
    292
  ),
  n_without_overlay = c(
    12,
    25,
    9
  ),
  status = c(
    "success",
    "success",
    "success"
  ),
  output_file = c(
    "Carcinus_maenas_eez.csv",
    "Mnemiopsis_leidyi_eez.csv",
    "Undaria_pinnatifida_eez.csv"
  ),
  stringsAsFactors = FALSE
)

marine_summary_example

## ----inspect-marine-summary-example-------------------------------------------
# Count successful and failed rows.
table(marine_summary_example$status, useNA = "ifany")

# Inspect species with at least one retained record lacking overlay assignment.
marine_summary_example[
  marine_summary_example$n_without_overlay > 0,
  c(
    "species",
    "overlay",
    "n_thinned",
    "n_with_overlay",
    "n_without_overlay"
  )
]

## ----marine-multiple-overlay-template, eval = FALSE---------------------------
# # Choose the overlays to compare.
# marine_overlays_to_run <- c(
#   "eez",
#   "lme",
#   "meow",
#   "iho"
# )
# 
# # Run each overlay separately.
# overlay_results <- lapply(
#   marine_overlays_to_run,
#   function(this_overlay) {
# 
#     process_gbif_marine_pipeline(
#       df = marine_species,
# 
#       output_dir = file.path(
#         output_root,
#         paste0("marine_", this_overlay, "_batch")
#       ),
# 
#       user = Sys.getenv("GBIF_USER"),
#       pwd = Sys.getenv("GBIF_PWD"),
#       email = Sys.getenv("GBIF_EMAIL"),
# 
#       overlay = this_overlay,
#       overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
# 
#       batch_size = 1,
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       prepare_taxonomy = TRUE,
# 
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   }
# )
# 
# names(overlay_results) <- marine_overlays_to_run

## ----compare-marine-overlay-summaries-template, eval = FALSE------------------
# # Build paths to each overlay summary file.
# summary_paths <- file.path(
#   output_root,
#   paste0("marine_", marine_overlays_to_run, "_batch"),
#   "gbif_summary.csv"
# )
# 
# # Read available summary files and add the requested overlay name.
# summary_list <- Map(
#   function(path, overlay_name) {
# 
#     if (!file.exists(path)) {
#       return(NULL)
#     }
# 
#     x <- read.csv(path, stringsAsFactors = FALSE)
#     x$overlay_requested <- overlay_name
#     x
#   },
#   summary_paths,
#   marine_overlays_to_run
# )
# 
# # Remove missing summaries.
# summary_list <- summary_list[
#   !vapply(summary_list, is.null, logical(1))
# ]
# 
# # Combine summaries if at least one was found.
# overlay_summary <- do.call(rbind, summary_list)
# 
# # Inspect the main audit fields.
# overlay_summary[
#   ,
#   intersect(
#     c(
#       "species",
#       "overlay_requested",
#       "n_total",
#       "n_cleaned",
#       "n_thinned",
#       "n_with_overlay",
#       "n_without_overlay",
#       "status",
#       "fail_stage",
#       "fail_reason",
#       "output_file"
#     ),
#     names(overlay_summary)
#   ),
#   drop = FALSE
# ]

## ----marine-cleaning-thinning-template, eval = FALSE--------------------------
# process_gbif_marine_pipeline(
#   df = marine_species,
# 
#   output_dir = file.path(output_root, "marine_cleaned_thinned"),
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   overlay = "eez",
#   overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
# 
#   # Remove invalid, missing, suspicious, or highly uncertain coordinates
#   # according to the pipeline cleaning settings.
#   apply_cleaning = TRUE,
# 
#   # Reduce dense spatial clustering after cleaning.
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   batch_size = 1,
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE
# )

## ----marine-origin-evidence-live, eval = run_live-----------------------------
# # Run the marine origin-evidence example only when GBIF credentials are available.
# if (!gbif_ready) {
# 
#   skip_live_message("set GBIF credentials to run the marine origin-evidence example")
# 
# } else {
# 
#   marine_origin_result <- run_live_safely(
#     process_gbif_marine_pipeline(
#       df = data.frame(
#         species = c(
#           "Carcinus maenas",
#           "Undaria pinnatifida"
#         ),
#         stringsAsFactors = FALSE
#       ),
# 
#       output_dir = file.path(output_root, "marine_origin"),
# 
#       user = gbif_user,
#       pwd = gbif_pwd,
#       email = gbif_email,
# 
#       # Assign retained records to a country-associated marine overlay.
#       # EEZ attribution is useful here because GRIIS-style evidence is
#       # interpreted at species-by-country or species-by-territory level.
#       overlay = "eez",
#       overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
#       strict_overlay_loading = FALSE,
# 
#       # Standardise taxonomy before attaching evidence.
#       prepare_taxonomy = TRUE,
# 
#       # Attach species-country invasive or alien-status evidence without
#       # removing records.
#       use_griis_status = TRUE,
#       griis_filter_mode = "audit_only",
#       griis_cache_dir = file.path(output_root, "cache", "griis"),
# 
#       # Attach native-range evidence without removing records.
#       use_native_status = TRUE,
#       use_native_web = TRUE,
#       native_web_sources = c(
#         "sinas",
#         "gbif",
#         "worms"
#       ),
#       native_web_cache_dir = file.path(output_root, "cache", "native_web"),
#       native_filter_mode = "audit_only",
# 
#       # Add reconciled audit fields that help compare GRIIS-style
#       # species-country evidence with native-range evidence.
#       reconcile_griis_native = TRUE,
# 
#       batch_size = 1,
# 
#       apply_cleaning = TRUE,
#       apply_thinning = TRUE,
#       dist_km = 5,
# 
#       export_summary = TRUE,
#       store_in_memory = FALSE,
#       return_all_results = FALSE,
#       quiet = FALSE
#     )
#   )
# 
#   marine_origin_result
# }

## ----inspect-exported-marine-file-live, eval = run_live-----------------------
# summary_dir <- file.path(output_root, "marine_low_memory")
# 
# summary_path <- file.path(
#   summary_dir,
#   "gbif_summary.csv"
# )
# 
# if (file.exists(summary_path)) {
# 
#   marine_summary <- read.csv(summary_path, stringsAsFactors = FALSE)
# 
#   if ("output_file" %in% names(marine_summary)) {
# 
#     candidate_files <- marine_summary$output_file[
#       !is.na(marine_summary$output_file) & nzchar(marine_summary$output_file)
#     ]
# 
#     if (length(candidate_files) > 0) {
# 
#       candidate_file <- resolve_output_file(
#         path = candidate_files[1],
#         output_dir = summary_dir
#       )
# 
#       if (file.exists(candidate_file)) {
#         marine_occ <- read.csv(candidate_file, stringsAsFactors = FALSE)
#         utils::head(marine_occ)
#       } else {
#         skip_live_message("no readable exported occurrence file was listed")
#       }
# 
#     } else {
#       skip_live_message("summary file did not list any exported occurrence files")
#     }
# 
#   } else {
#     skip_live_message("summary file did not contain an output_file column")
#   }
# 
# } else {
#   skip_live_message("summary file was not available from the live marine run")
# }

## ----marine-scaling-template, eval = FALSE------------------------------------
# marine_batch_result <- process_gbif_marine_pipeline(
#   df = marine_species,
# 
#   output_dir = file.path(output_root, "marine_eez_batch"),
# 
#   user = Sys.getenv("GBIF_USER"),
#   pwd = Sys.getenv("GBIF_PWD"),
#   email = Sys.getenv("GBIF_EMAIL"),
# 
#   overlay = "eez",
#   overlay_cache_dir = file.path(output_root, "cache", "marine_overlays"),
# 
#   batch_size = 1,
# 
#   apply_cleaning = TRUE,
#   apply_thinning = TRUE,
#   dist_km = 5,
# 
#   prepare_taxonomy = TRUE,
# 
#   use_griis_status = TRUE,
#   griis_filter_mode = "audit_only",
#   griis_cache_dir = file.path(output_root, "cache", "griis"),
# 
#   use_native_status = TRUE,
#   use_native_web = TRUE,
#   native_web_sources = c("sinas", "gbif", "worms"),
#   native_web_cache_dir = file.path(output_root, "cache", "native_web"),
#   native_filter_mode = "audit_only",
#   reconcile_griis_native = TRUE,
# 
#   export_summary = TRUE,
#   store_in_memory = FALSE,
#   return_all_results = FALSE,
#   quiet = FALSE
# )

## ----marine-troubleshooting-table, echo = FALSE-------------------------------
marine_troubleshooting <- data.frame(
  symptom = c(
    "GBIF query does not run",
    "No records after cleaning",
    "Many records lack overlay assignment",
    "Overlay fails to load",
    "Output files are not where expected",
    "Taxon names need review",
    "Origin or invasive-status evidence is unresolved",
    "Run uses too much memory"
  ),
  possible_cause = c(
    "GBIF credentials may be missing, incorrect, or unavailable in the R session.",
    "Coordinates may be missing, impossible, exact zero-zero points, or too uncertain.",
    "Records may fall outside the selected overlay, near boundaries, in provider gaps, or in a different marine spatial domain.",
    "Provider service, URL, cache, optional dependency, or local file may be unavailable.",
    "The run may have failed before export, or output_file may be relative to the output directory.",
    "GBIF names may require synonym resolution, spelling review, or manual checking.",
    "Native-range or GRIIS evidence may not resolve cleanly for the selected taxon and spatial unit.",
    "Large GBIF downloads, dense coastal records, multiple overlays, or in-memory storage may be too heavy."
  ),
  first_check = c(
    "Check GBIF_USER, GBIF_PWD, GBIF_EMAIL, and the live-run messages.",
    "Inspect n_total, n_cleaned, n_thinned, fail_stage, and fail_reason in gbif_summary.csv.",
    "Inspect unmatched coordinates and consider whether EEZ, LME, MEOW, or IHO is the right overlay.",
    "Load the overlay directly with bf_load_marine_regions_overlay().",
    "Inspect output_dir, status, fail_stage, fail_reason, and output_file in gbif_summary.csv.",
    "Run taxonomy preparation and inspect the taxonomy audit outputs.",
    "Use audit_only mode first and inspect evidence fields before filtering.",
    "Use batch_size = 1, store_in_memory = FALSE, and add overlays or evidence layers in stages."
  ),
  stringsAsFactors = FALSE
)

marine_troubleshooting

## ----provider-citation-table, echo = FALSE------------------------------------
provider_citations <- data.frame(
  source = c(
    "GBIF occurrence data",
    "Marine Regions",
    "GRIIS",
    "WoRMS",
    "SInAS",
    "mregions2"
  ),
  cite_when_used_for = c(
    "Occurrence downloads from GBIF.",
    "EEZ, LME, MEOW, IHO, high-seas, continental-shelf, or other Marine Regions products.",
    "Species-by-country or species-by-territory alien/invasive status evidence.",
    "Marine taxonomic names or WoRMS-derived native/distribution evidence.",
    "Native and alien distribution evidence from SInAS.",
    "R access to Marine Regions gazetteer or data products."
  ),
  citation_note = c(
    "Cite the GBIF download DOI generated for the exact occurrence download.",
    "Cite MarineRegions.org and, where relevant, the specific Marine Regions product used.",
    "Cite Pagad et al. (2018), Scientific Data, 5, 170202, and any checklist-specific source where used.",
    "Cite the WoRMS dataset DOI and/or the recommended WoRMS citation from the provider.",
    "Cite the SInAS dataset version used and the associated Scientific Data paper.",
    "Cite mregions2 if the R package was used directly in the workflow."
  ),
  stringsAsFactors = FALSE
)

provider_citations

