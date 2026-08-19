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
#   - This file does not submit GBIF downloads itself; submission and queue
#     management are handled by gbif_download_backends.R.
#   - `check_gbif_presence()` uses the GBIF occurrence-search API for a lightweight
#     availability pre-check. Search API calls do not themselves create a
#     reproducible GBIF occurrence-download DOI.
#   - `wait_and_import_gbif()` treats request labels and GBIF download keys as
#     separate concepts. Several request labels can intentionally share one GBIF
#     download key.
#   - Shared download keys are deduplicated before polling/import, so one GBIF
#     archive is downloaded and imported only once.
#   - When a download object carries `biofetchR_split_field`, the imported archive
#     is split back to its original request labels after import. The
#     terrestrial/freshwater backend uses `countryCode` for this purpose.
#   - Polling timeouts remain explicit unresolved states. They retain the GBIF key
#     and are never interpreted as evidence that GBIF contained zero records.
#   - `wait_and_import_gbif()` attaches a request-level retrieval audit to the
#     returned list so successful, empty, failed, and unresolved requests remain
#     distinguishable.
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
#' Poll submitted GBIF occurrence-download keys, import successfully completed
#' SIMPLE_CSV archives, convert coordinate-bearing records to WGS84 `sf` points,
#' and preserve one explicit retrieval outcome for every original request label.
#'
#' @details
#' `keys` can be a named list or named character vector. Names represent logical
#' request labels such as species names or ISO2 countries, while values are GBIF
#' asynchronous download keys.
#'
#' Multiple labels are allowed to share the same GBIF key. Before polling,
#' biofetchR constructs the request-label-to-key mapping and reduces the supplied
#' keys to their unique non-missing values. Each unique GBIF download is therefore
#' polled, downloaded and imported only once.
#'
#' For each unique key, GBIF metadata are polled up to `max_tries` times with
#' `wait_time` seconds between attempts. `SUCCEEDED` proceeds to import.
#' `FAILED`, `KILLED`, `CANCELLED`, and `FILE_ERASED` are treated as terminal
#' failure states. If no terminal state is reached before the polling limit, the
#' request is recorded as `pending_timeout`; the GBIF key is retained for later
#' review or recovery.
#'
#' Completed archives are imported with [rgbif::occ_download_import()]. Records
#' lacking `decimalLatitude` or `decimalLongitude` are removed before conversion
#' to EPSG:4326 point geometry.
#'
#' @section Shared download keys and country splitting:
#' The terrestrial/freshwater multi-country backend deliberately returns several
#' country labels pointing to one shared GBIF download key and attaches
#' `biofetchR_split_field = "countryCode"`.
#'
#' When a split field is present, the completed archive is imported once and then
#' partitioned locally for each original label. A country with no matching rows
#' after a successful import is retained as a zero-row WGS84 `sf` object and is
#' audited as `success_zero`, rather than disappearing from the returned result.
#'
#' If no split field is supplied, each request label associated with a successful
#' key receives the complete imported `sf` object.
#'
#' @section Retrieval-status audit:
#' The returned list carries a `gbif_status` attribute containing one row per
#' original request label. Its columns are:
#'
#' \itemize{
#'   \item `label`: original request label;
#'   \item `gbif_key`: associated GBIF download key;
#'   \item `gbif_status`: terminal or last observed GBIF download status;
#'   \item `retrieval_status`: biofetchR retrieval interpretation;
#'   \item `n_records`: number of coordinate-bearing records assigned to the
#'     request label when known;
#'   \item `message`: diagnostic or recovery message where relevant.
#' }
#'
#' Retrieval statuses include `success`, `success_zero`, `invalid_key`,
#' `pending_timeout`, `download_failed`, `killed`, `cancelled`, `file_erased`,
#' `import_failed`, and `split_failed`.
#'
#' Importantly, `pending_timeout` means that the download remained unresolved
#' within the configured polling window. It is not interpreted as a biological
#' zero and the associated GBIF key remains available in the audit.
#'
#' @param keys Named list or named character vector of GBIF occurrence-download
#'   keys. Names are logical request labels. Several names may intentionally map
#'   to the same key.
#' @param wait_time Numeric. Seconds between GBIF metadata polling attempts.
#' @param max_tries Integer. Maximum number of metadata polling attempts for each
#'   UNIQUE GBIF download key.
#'
#' @return A named list of WGS84 `sf` point objects for successfully resolved
#'   request labels. For a normal species-level download, one imported object is
#'   returned under the corresponding species label. For a shared multi-country
#'   download, the archive is imported once and separate country-labelled `sf`
#'   objects are returned after splitting by the configured field.
#'
#'   Successfully resolved labels with zero matching records are retained as
#'   zero-row `sf` objects. Failed or unresolved requests need not appear as list
#'   elements, but every original request is represented in the attached
#'   `gbif_status` audit attribute.
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
#'   attr(imported, "gbif_status")
#' }
#' }
#'
#' @family GBIF helpers
#' @md
#' @export
wait_and_import_gbif <- function(keys, wait_time = 30, max_tries = 60) {
  split_field <- attr(keys, "biofetchR_split_field", exact = TRUE)
  requested_labels_attr <- attr(keys, "biofetchR_requested_labels", exact = TRUE)

  if (is.list(keys)) {
    key_vec <- vapply(
      keys,
      function(x) {
        x <- as.character(x)
        x <- x[!is.na(x) & nzchar(x)]
        if (length(x)) x[[1L]] else NA_character_
      },
      character(1)
    )
  } else {
    key_vec <- as.character(keys)
  }

  labels <- names(keys)
  if (is.null(labels) || length(labels) != length(key_vec) || any(is.na(labels) | !nzchar(labels))) {
    if (!is.null(requested_labels_attr) && length(requested_labels_attr) == length(key_vec)) {
      labels <- as.character(requested_labels_attr)
    } else {
      labels <- paste0("request_", seq_along(key_vec))
    }
  }

  request_map <- data.frame(
    label = as.character(labels),
    gbif_key = as.character(key_vec),
    stringsAsFactors = FALSE
  )

  audit <- data.frame(
    label = request_map$label,
    gbif_key = request_map$gbif_key,
    gbif_status = NA_character_,
    retrieval_status = ifelse(
      is.na(request_map$gbif_key) | !nzchar(request_map$gbif_key),
      "invalid_key",
      "pending"
    ),
    n_records = NA_integer_,
    message = NA_character_,
    stringsAsFactors = FALSE
  )

  results <- list()

  .set_audit <- function(labels_for_key,
                         gbif_status = NA_character_,
                         retrieval_status,
                         n_records = NA_integer_,
                         message = NA_character_) {
    idx <- which(audit$label %in% labels_for_key)
    if (!length(idx)) return(invisible(NULL))

    audit$gbif_status[idx] <<- gbif_status
    audit$retrieval_status[idx] <<- retrieval_status
    audit$n_records[idx] <<- as.integer(n_records)
    audit$message[idx] <<- message
    invisible(NULL)
  }

  .to_sf <- function(df) {
    if (!is.data.frame(df)) return(NULL)
    if (!all(c("decimalLatitude", "decimalLongitude") %in% names(df))) return(NULL)

    df <- dplyr::filter(
      df,
      !is.na(.data$decimalLatitude),
      !is.na(.data$decimalLongitude)
    )

    if (!nrow(df)) {
      return(
        sf::st_sf(
          df,
          geometry = sf::st_sfc(crs = 4326)
        )
      )
    }

    sf::st_as_sf(
      df,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326
    )
  }

  valid_keys <- unique(
    request_map$gbif_key[
      !is.na(request_map$gbif_key) & nzchar(request_map$gbif_key)
    ]
  )

  for (key in valid_keys) {
    labels_for_key <- request_map$label[request_map$gbif_key == key]
    final_status <- NA_character_
    success <- FALSE

    for (i in seq_len(max_tries)) {
      status <- tryCatch(
        rgbif::occ_download_meta(key)$status,
        error = function(e) NA_character_
      )

      if (!is.na(status) && nzchar(status)) {
        final_status <- toupper(as.character(status[[1L]]))

        if (identical(final_status, "SUCCEEDED")) {
          success <- TRUE
          break
        }

        if (final_status %in% c("FAILED", "KILLED", "CANCELLED", "FILE_ERASED")) {
          break
        }
      }

      Sys.sleep(wait_time)
    }

    if (!isTRUE(success)) {
      retrieval_status <- if (identical(final_status, "FAILED")) {
        "download_failed"
      } else if (identical(final_status, "KILLED")) {
        "killed"
      } else if (identical(final_status, "CANCELLED")) {
        "cancelled"
      } else if (identical(final_status, "FILE_ERASED")) {
        "file_erased"
      } else {
        "pending_timeout"
      }

      msg <- if (identical(retrieval_status, "pending_timeout")) {
        paste0(
          "Polling limit reached before GBIF download completed; key retained for recovery: ",
          key
        )
      } else {
        paste0("GBIF download ended with status ", final_status, ".")
      }

      .set_audit(
        labels_for_key,
        gbif_status = final_status,
        retrieval_status = retrieval_status,
        message = msg
      )

      if (identical(retrieval_status, "pending_timeout")) {
        cli::cli_alert_warning(
          "[timeout] GBIF download {.val {key}} is still unresolved; key retained for recovery."
        )
      } else {
        cli::cli_alert_danger(
          "[x] GBIF download {.val {key}} ended with status {.val {final_status}}."
        )
      }

      next
    }

    imported_df <- tryCatch(
      {
        zipfile <- rgbif::occ_download_get(
          key,
          overwrite = TRUE,
          path = tempdir()
        )
        rgbif::occ_download_import(zipfile)
      },
      error = function(e) e
    )

    if (inherits(imported_df, "condition")) {
      .set_audit(
        labels_for_key,
        gbif_status = "SUCCEEDED",
        retrieval_status = "import_failed",
        message = conditionMessage(imported_df)
      )
      cli::cli_alert_danger(
        "[x] Failed to import GBIF download {.val {key}}: {conditionMessage(imported_df)}"
      )
      next
    }

    if (!is.data.frame(imported_df)) {
      .set_audit(
        labels_for_key,
        gbif_status = "SUCCEEDED",
        retrieval_status = "import_failed",
        message = "GBIF import did not return a data frame."
      )
      next
    }

    if (!is.null(split_field) && nzchar(as.character(split_field[[1L]]))) {
      split_field <- as.character(split_field[[1L]])

      if (!split_field %in% names(imported_df)) {
        .set_audit(
          labels_for_key,
          gbif_status = "SUCCEEDED",
          retrieval_status = "split_failed",
          message = paste0("Imported GBIF data do not contain split field `", split_field, "`.")
        )
        next
      }

      split_values <- toupper(trimws(as.character(imported_df[[split_field]])))

      for (label in labels_for_key) {
        keep <- !is.na(split_values) & split_values == toupper(trimws(label))
        df_label <- imported_df[keep, , drop = FALSE]
        sf_label <- tryCatch(.to_sf(df_label), error = function(e) NULL)

        if (is.null(sf_label)) {
          .set_audit(
            label,
            gbif_status = "SUCCEEDED",
            retrieval_status = "split_failed",
            message = paste0("Failed to convert country subset `", label, "` to sf.")
          )
          next
        }

        results[[label]] <- sf_label
        n_label <- nrow(sf_label)
        .set_audit(
          label,
          gbif_status = "SUCCEEDED",
          retrieval_status = if (n_label > 0L) "success" else "success_zero",
          n_records = n_label
        )
      }
    } else {
      occ_sf <- tryCatch(.to_sf(imported_df), error = function(e) NULL)

      if (is.null(occ_sf)) {
        .set_audit(
          labels_for_key,
          gbif_status = "SUCCEEDED",
          retrieval_status = "import_failed",
          message = "Failed to convert imported GBIF records to sf."
        )
        next
      }

      for (label in labels_for_key) {
        results[[label]] <- occ_sf
        .set_audit(
          label,
          gbif_status = "SUCCEEDED",
          retrieval_status = if (nrow(occ_sf) > 0L) "success" else "success_zero",
          n_records = nrow(occ_sf)
        )
      }
    }
  }

  attr(results, "gbif_status") <- audit
  results
}
