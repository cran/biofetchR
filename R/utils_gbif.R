################################################################################
# utils_gbif.R
# -----------------------------------------------------------------------------
# biofetchR: lightweight GBIF lookup and import helpers
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide small GBIF-facing utilities used by the retained biofetchR download
#   and processing workflows. These helpers cover three narrow tasks:
#
#   1. checking whether GBIF has coordinate-bearing records for a species in a
#      supplied recipient region;
#   2. resolving a scientific name against the GBIF backbone taxonomy to retrieve
#      a taxonKey; and
#   3. polling submitted GBIF occurrence-download keys, importing completed
#      downloads, and coercing returned records to WGS84 sf point objects.
#
# DESIGN NOTES
#   - This file is intentionally small and does not submit GBIF downloads itself;
#     submission is handled by the dedicated GBIF download-backend script.
#   - `check_gbif_presence()` uses the GBIF occurrence-search API for a quick
#     presence/absence pre-check. API search results do not themselves create a
#     single reproducible download DOI.
#   - `wait_and_import_gbif()` imports completed occurrence downloads and keeps
#     the user-supplied download key as the reproducibility anchor in downstream
#     summaries.
#   - `harmonize_column_types()` is retained for backwards compatibility with
#     older scripts that bind heterogeneous GBIF sf objects.
#
# GBIF DATA USE AND CITATION
#   GBIF-mediated data remain subject to GBIF data-use terms and the licences of
#   the original publishing datasets. When occurrence-download data are used in
#   analyses, retain the GBIF download key/DOI and cite the relevant GBIF
#   occurrence download or derived dataset. For direct occurrence-search API
#   results, users should track dataset keys and create an appropriate citation
#   route where the records are used analytically.
#
################################################################################

#' Check whether GBIF has coordinate-based records for a species
#'
#' Query the GBIF occurrence-search API through [rgbif::occ_search()] and return
#' a simple `TRUE`/`FALSE` flag indicating whether at least one coordinate-bearing
#' record is available for a species. The helper is intended as a lightweight
#' pre-check before running heavier GBIF download workflows.
#'
#' @details
#' By default, `region_id` is passed to GBIF as an ISO2 country filter. When
#' `use_eez_for_marine = TRUE` and the package-level `is_marine_species()` helper
#' identifies the species as marine, the country predicate is skipped so the
#' query behaves as a global marine pre-check.
#'
#' This function only checks record availability. It does not submit or import a
#' GBIF occurrence download and therefore does not generate a download DOI.
#'
#' @param species Character scalar. Scientific name to query.
#' @param region_id Character scalar or `NULL`. ISO2 country code used for the
#'   GBIF country predicate. Use `NULL` for a global check.
#' @param use_eez_for_marine Logical. If `TRUE`, skip the country predicate for
#'   species detected as marine by `is_marine_species()`.
#' @param verbose Logical. If `TRUE`, print a short CLI status message.
#'
#' @return A logical value of length one. The value is `TRUE` when the GBIF
#'   occurrence-search query returns at least one coordinate-bearing record for
#'   the supplied species and region settings. The value is `FALSE` when no
#'   matching coordinate-bearing record is returned or when the query fails. This
#'   result is an availability pre-check only and does not represent a GBIF
#'   occurrence download or download DOI.
#'
#' @section GBIF data use:
#' This helper uses GBIF's occurrence-search API for a quick availability check.
#' Records used in analyses should still be cited through an appropriate GBIF
#' download DOI, derived dataset, or dataset-level citation workflow.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   has_records <- check_gbif_presence(
#'     species = "Carcinus maenas"
#'   )
#'
#'   has_records
#' }
#' }
#'
#' @family GBIF helpers
#' @md
#' @export
check_gbif_presence <- function(species, region_id = NULL, use_eez_for_marine = FALSE, verbose = FALSE) {
  out <- tryCatch({
    args <- list(scientificName = species, hasCoordinate = TRUE, limit = 1)
    if (!(use_eez_for_marine && is_marine_species(species))) {
      args$country <- region_id
    }

    res <- do.call(rgbif::occ_search, args)
    has_records <- is.data.frame(res$data) && nrow(res$data) > 0

    if (verbose) {
      msg <- if (has_records) "[OK] Records found" else "[x] No records"
      cli::cli_alert_info("{.val {species}} in {.val {region_id %||% 'global'}}: {msg}")
    }

    has_records
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{.val {species}}: Error - {e$message}")
    FALSE
  })

  return(out)
}

#' Retrieve a GBIF taxonKey for a scientific name
#'
#' Resolve a supplied name against the GBIF backbone taxonomy using
#' [rgbif::name_backbone()] and return the matched `usageKey` when available.
#'
#' @details
#' This helper is a thin, fail-soft wrapper around `rgbif::name_backbone()`. It
#' returns `NA_integer_` rather than stopping when GBIF name matching fails, which
#' makes it suitable for batch pipelines and audit tables.
#'
#' @param species Character scalar. Scientific name to resolve.
#' @param verbose Logical. If `TRUE`, print a CLI status message.
#'
#' @return An integer value of length one containing the matched GBIF backbone
#'   `usageKey` for `species`. Returns `NA_integer_` when no usage key is found
#'   or when the GBIF backbone lookup fails. The value can be used as a GBIF
#'   taxon identifier in downstream occurrence-download predicates and audit
#'   tables.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   key <- get_taxon_key("Carcinus maenas")
#'
#'   key
#' }
#' }
#'
#' @family GBIF helpers
#' @md
#' @export
get_taxon_key <- function(species, verbose = FALSE) {
  key <- tryCatch({
    backbone <- rgbif::name_backbone(name = species)
    if (!is.null(backbone$usageKey)) {
      if (verbose) cli::cli_alert_success("{.val {species}}: taxonKey = {.val {backbone$usageKey}}")
      backbone$usageKey
    } else {
      if (verbose) cli::cli_alert_warning("{.val {species}}: [x] No taxonKey found")
      NA_integer_
    }
  }, error = function(e) {
    if (verbose) cli::cli_alert_danger("{.val {species}}: Error - {e$message}")
    NA_integer_
  })

  return(key)
}

#' Wait for GBIF occurrence downloads and import completed records
#'
#' Poll GBIF download metadata until each supplied download key succeeds, then
#' download, import and convert the returned occurrence records to WGS84 `sf`
#' points.
#'
#' @details
#' `keys` should normally be the named list returned by a GBIF download-submission
#' helper, where names represent the species, country, region or batch label used
#' by the calling pipeline. Each successful import is stored under the same name.
#'
#' Completed downloads are imported with [rgbif::occ_download_import()], filtered
#' to records with non-missing `decimalLatitude` and `decimalLongitude`, and then
#' converted to an `sf` point object with EPSG:4326 coordinates. Failed,
#' cancelled, killed or timed-out downloads are skipped with a CLI message.
#'
#' @param keys Named list or named character vector of GBIF occurrence download
#'   keys. Names are used as the returned list names.
#' @param wait_time Numeric. Seconds to wait between polling attempts.
#' @param max_tries Integer. Maximum number of polling attempts per key.
#'
#' @return A named list of `sf` point objects in WGS84 longitude/latitude
#'   coordinates. Each list element contains imported GBIF occurrence records for
#'   one completed download key, filtered to records with non-missing
#'   `decimalLongitude` and `decimalLatitude` values and converted to EPSG:4326
#'   point geometry. List names are taken from `names(keys)`, so they usually
#'   represent the species, country, region or batch label supplied by the
#'   calling pipeline. Downloads that fail, are cancelled, time out or cannot be
#'   imported are skipped, so the returned list may be shorter than `keys` or
#'   empty.
#'
#' @section GBIF data use and citation:
#' This function imports GBIF occurrence-download data. Retain the GBIF download
#' key/DOI in downstream summaries and cite the corresponding GBIF occurrence
#' download, or an appropriate derived dataset, when using the records in
#' analyses or publications.
#'
#' @examples
#' \donttest{
#' if (interactive() && exists("keys")) {
#'   imported <- wait_and_import_gbif(keys)
#'
#'   names(imported)
#' }
#' }
#'
#' @family GBIF helpers
#' @md
#' @export
wait_and_import_gbif <- function(keys, wait_time = 30, max_tries = 60) {
  results <- list()

  for (region_id in names(keys)) {
    key <- keys[[region_id]]
    success <- FALSE

    for (i in seq_len(max_tries)) {
      status <- tryCatch({
        rgbif::occ_download_meta(key)$status
      }, error = function(e) NA_character_)

      if (!is.na(status)) {
        if (status == "SUCCEEDED") {
          success <- TRUE
          break
        } else if (status %in% c("KILLED", "CANCELLED")) {
          cli::cli_alert_danger("[x] Download {key} failed with status: {status}")
          break
        }
      }

      Sys.sleep(wait_time)
    }

    if (!success) {
      cli::cli_alert_warning("[timeout] Timeout waiting for download {.val {key}} - skipping.")
      next
    }

    occ_data <- tryCatch({
      zipfile <- rgbif::occ_download_get(key, overwrite = TRUE, path = tempdir())
      df <- rgbif::occ_download_import(zipfile)
      df <- dplyr::filter(df, !is.na(decimalLatitude), !is.na(decimalLongitude))
      sf::st_as_sf(df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    }, error = function(e) {
      cli::cli_alert_danger("[x] Failed to import download {.val {key}}: {e$message}")
      NULL
    })

    if (!is.null(occ_data)) {
      results[[region_id]] <- occ_data
    }
  }

  return(results)
}

