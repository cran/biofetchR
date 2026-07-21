################################################################################
# summary.R
# -----------------------------------------------------------------------------
# biofetchR: GBIF processing summary helpers
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Provide small, package-level helpers for creating and updating processing
#   summary tables used by the terrestrial/freshwater and marine GBIF pipelines.
#
# MAIN HELPERS
#   - initialize_summary()
#   - append_summary_row()
#
# DESIGN NOTES
#   - These helpers track processing outcomes only; they do not download, import,
#     clean, thin or spatially join GBIF records.
#   - The summary schema is intentionally simple and stable so that multiple
#     pipeline components can append rows consistently.
#   - Extended pipeline-specific metadata, such as failure stage or GBIF download
#     key, can be added by higher-level pipeline functions after row creation.
#
# LICENSING / DATA USE
#   This script does not access, download or redistribute external datasets.
#   Summary rows may contain paths or identifiers linked to provider datasets
#   handled elsewhere in biofetchR. Users remain responsible for following the
#   licence, citation and attribution requirements of those source datasets.
#
################################################################################


#' Initialise a GBIF processing summary table
#'
#' Create an empty, typed summary table for tracking the outcome of GBIF
#' occurrence processing across species and spatial recipient units.
#'
#' `initialize_summary()` is used by biofetchR pipelines before records are
#' downloaded, imported, cleaned, thinned and exported. The returned table has a
#' stable schema that can be extended by higher-level workflows when additional
#' audit fields are needed.
#'
#' @details
#' The summary table is designed around one row per processed species x region
#' combination. The `region_id` and `region_type` fields are intentionally
#' generic so the same helper can be used for terrestrial/freshwater workflows
#' based on GADM units and for marine workflows based on EEZs or other marine
#' overlays.
#'
#' The base columns are:
#'
#' \describe{
#'   \item{`species`}{Scientific name or accepted taxon label used by the pipeline.}
#'   \item{`region_id`}{Spatial recipient identifier, such as a GADM code, EEZ
#'   label or marine-region identifier.}
#'   \item{`region_type`}{Spatial-recipient family, such as `"GADM"`, `"EEZ"` or
#'   another overlay label.}
#'   \item{`n_total`}{Number of records before coordinate cleaning.}
#'   \item{`n_cleaned`}{Number of records remaining after coordinate cleaning.}
#'   \item{`n_thinned`}{Number of records remaining after spatial thinning.}
#'   \item{`output_file`}{Path to the exported CSV file, or `NA` when records are
#'   retained only in memory.}
#'   \item{`status`}{Processing status, such as `"success"`, `"empty"` or
#'   `"failed"`.}
#' }
#'
#' @return A zero-row [tibble::tibble()] with typed columns used to track GBIF
#'   processing outcomes. The returned table contains character columns
#'   `species`, `region_id`, `region_type`, `output_file` and `status`, and
#'   integer columns `n_total`, `n_cleaned` and `n_thinned`. Each future row is
#'   intended to represent one processed species x spatial-recipient
#'   combination, with counts for records before cleaning, after coordinate
#'   cleaning and after spatial thinning. The empty table is used as the starting
#'   summary object for terrestrial/freshwater and marine processing pipelines.
#'
#' @examples
#' summary_tbl <- initialize_summary()
#' summary_tbl
#'
#' @family processing summary helpers
#' @export
#' @md
initialize_summary <- function() {
  tibble::tibble(
    species     = character(),  # Scientific name (e.g., "Rattus norvegicus")
    region_id   = character(),  # GADM country code or EEZ name (e.g., "US", "Gulf of Guinea")
    region_type = character(),  # "GADM" (terrestrial) or "EEZ" (marine)
    n_total     = integer(),    # Number of raw GBIF records before filtering
    n_cleaned   = integer(),    # Number after coordinate cleaning
    n_thinned   = integer(),    # Number after spatial thinning
    output_file = character(),  # Path to exported CSV or NA if in-memory
    status      = character()   # Status of operation ("success", "no_data", "error", etc.)
  )
}

#' Append a row to a GBIF processing summary table
#'
#' Record the processing outcome for one species x region combination and append
#' it to a summary table created by [initialize_summary()].
#'
#' @details
#' This helper performs light type coercion for the count fields so that summary
#' tables remain stable when upstream pipeline steps return numeric, integer or
#' missing values. It does not validate that `n_cleaned <= n_total` or
#' `n_thinned <= n_cleaned`, because some workflows may report values from
#' different stages or use specialised cleaning/thinning backends.
#'
#' When `quiet = FALSE`, the helper prints a short success message using
#' [cli::cli_alert_success()]. Set `quiet = TRUE` in batch tests, package checks
#' or non-interactive workflows where console output should be suppressed.
#'
#' @param summary_tbl A data frame or tibble, usually returned by
#'   [initialize_summary()].
#' @param species Character scalar. Scientific name or accepted taxon label.
#' @param region_id Character scalar. Spatial recipient identifier, such as a
#'   GADM code, EEZ name or marine-region identifier.
#' @param region_type Character scalar. Spatial-recipient family or overlay label,
#'   such as `"GADM"` or `"EEZ"`.
#' @param n_total Integer-compatible value. Number of records before cleaning.
#' @param n_cleaned Integer-compatible value. Number of records after coordinate
#'   cleaning.
#' @param n_thinned Integer-compatible value. Number of records after spatial
#'   thinning.
#' @param output_file Character scalar. Exported CSV path, or `NA` when output is
#'   retained in memory.
#' @param status Character scalar. Processing status label.
#' @param quiet Logical. If `TRUE`, suppress the console status message.
#'
#' @return A tibble containing the existing rows in `summary_tbl` plus one
#'   appended processing-summary row. The appended row records the supplied
#'   `species`, `region_id`, `region_type`, record counts, `output_file` and
#'   `status`. Count fields are coerced to integer so the summary schema remains
#'   stable across pipeline steps. The returned table is used to accumulate one
#'   row per processed species x spatial-recipient combination in
#'   terrestrial/freshwater and marine GBIF workflows.
#'
#' @examples
#' summary_tbl <- initialize_summary()
#'
#' summary_tbl <- append_summary_row(
#'   summary_tbl = summary_tbl,
#'   species = "Rattus norvegicus",
#'   region_id = "GB",
#'   region_type = "GADM",
#'   n_total = 100,
#'   n_cleaned = 90,
#'   n_thinned = 25,
#'   output_file = "Rattus_norvegicus__GB.csv",
#'   status = "success",
#'   quiet = TRUE
#' )
#'
#' summary_tbl
#'
#' @family processing summary helpers
#' @export
#' @md
append_summary_row <- function(summary_tbl,
                               species,
                               region_id,
                               region_type,
                               n_total,
                               n_cleaned,
                               n_thinned,
                               output_file,
                               status,
                               quiet = FALSE) {
  stopifnot(is.data.frame(summary_tbl))

  new_row <- tibble::tibble(
    species     = species,
    region_id   = region_id,
    region_type = region_type,
    n_total     = as.integer(n_total),
    n_cleaned   = as.integer(n_cleaned),
    n_thinned   = as.integer(n_thinned),
    output_file = output_file,
    status      = status
  )

  if (!quiet) {
    cli::cli_alert_success(
      "\U0001F4CA Summary updated: {.val {species}} in {.val {region_id}} ({.val {region_type}}) - {.val {status}} ({.val {n_thinned}} / {.val {n_total}})"
    )
  }

  dplyr::bind_rows(summary_tbl, new_row)
}
