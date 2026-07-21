################################################################################
# utils_pipeline_console.R
# -----------------------------------------------------------------------------
# biofetchR: internal console and progress-message helpers
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Centralise the user-facing console output used by biofetchR's retained GBIF
#   processing workflows. The terrestrial/freshwater and marine pipelines both
#   call these helpers so progress reporting uses the same vocabulary, formatting
#   and quiet-mode behaviour across workflows.
#
# WHAT THIS FILE DOES
#   - formats pipeline start/end messages;
#   - reports species-level progress;
#   - reports GBIF request and import milestones;
#   - summarises coordinate-cleaning and spatial-thinning retention;
#   - summarises spatial-overlay matching success;
#   - formats counts, percentages and elapsed time consistently;
#   - repairs text to UTF-8 before printing, which reduces Windows/locale issues.
#
# DESIGN PRINCIPLES
#   1. Keep console helpers lightweight and dependency-safe.
#   2. Use {cli} when installed, but always provide base R fallbacks.
#   3. Respect `quiet = TRUE` consistently.
#   4. Avoid changing pipeline state. These functions should report progress
#      only; they should not alter occurrence tables, overlay assignments or
#      export outputs.
#   5. Keep all helpers internal. They are documented with `@noRd` so developers
#      can understand the console layer without creating public help pages.
#
# PACKAGE ROLE
#   This is an internal reporting layer, not a data-access layer. It does not
#   download, query, transform or redistribute third-party biodiversity data.
#   Dataset-specific citation and licence responsibilities belong to the GBIF,
#   GRIIS, SInAS, overlay and raster helper scripts that actually access data.
#
################################################################################

# -----------------------------------------------------------------------------
# Internal formatting helpers
# -----------------------------------------------------------------------------

#' Return the current time for pipeline timing
#'
#' Small wrapper around `Sys.time()` used to mark the start of a pipeline, species-level step or GBIF request/import step before elapsed time is reported.
#'
#' @return A `POSIXct` timestamp.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_now <- function() {
  Sys.time()
}


#' Return a scalar console-safe text value
#'
#' Repairs text to UTF-8 and returns the first non-missing scalar value. A fallback is returned when the input is empty or missing.
#'
#' @param x Object to coerce to text.
#' @param fallback Character scalar returned when `x` is empty or missing.
#'
#' @return A character scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_text <- function(x, fallback = "") {
  y <- bf_repair_utf8_chr(x)

  if (is.null(y) || length(y) == 0L) return(fallback)
  if (is.na(y[[1]])) return(fallback)

  y[[1]]
}


#' Return a console-safe text vector
#'
#' Repairs text to UTF-8 and replaces missing values with the literal string `"NA"` so vectors can be collapsed safely in progress messages.
#'
#' @param x Object to coerce to a character vector.
#'
#' @return A character vector.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_text_vec <- function(x) {
  y <- bf_repair_utf8_chr(x)

  if (is.null(y) || length(y) == 0L) {
    return(character(0))
  }

  y[is.na(y)] <- "NA"
  y
}


#' Extract a single numeric value safely
#'
#' Converts the first element of an object to numeric for count, percentage and timing formatting. Returns a fallback when conversion fails.
#'
#' @param x Object containing a numeric-like value.
#' @param fallback Numeric value returned when `x` is empty, missing or non-numeric.
#'
#' @return A numeric scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_num1 <- function(x, fallback = NA_real_) {
  if (is.null(x) || length(x) == 0L) return(fallback)

  y <- suppressWarnings(as.numeric(x[[1]]))

  if (is.na(y)) return(fallback)

  y
}


#' Calculate elapsed seconds from a start time
#'
#' Computes the difference between the current time and a stored start time, returning seconds for downstream formatting.
#'
#' @param start_time Timestamp returned by `bf_console_now()` or `Sys.time()`.
#'
#' @return A numeric scalar giving elapsed seconds, or `NA_real_`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_elapsed <- function(start_time) {
  out <- suppressWarnings(as.numeric(difftime(Sys.time(), start_time, units = "secs")))

  if (length(out) == 0L || is.na(out[[1]])) {
    return(NA_real_)
  }

  out[[1]]
}


#' Format elapsed time for console messages
#'
#' Formats elapsed seconds as a compact human-readable string, using seconds for short operations and minutes plus seconds for longer operations.
#'
#' @param seconds Numeric-like elapsed time in seconds.
#'
#' @return A character scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_fmt_time <- function(seconds) {
  seconds <- bf_console_num1(seconds)

  if (is.na(seconds)) return("NA")

  if (seconds < 1) {
    return(sprintf("%.2f sec", seconds))
  }

  if (seconds < 60) {
    return(sprintf("%.1f sec", seconds))
  }

  mins <- floor(seconds / 60)
  secs <- seconds %% 60

  sprintf("%d min %.1f sec", mins, secs)
}


#' Format integer-like counts for console messages
#'
#' Rounds numeric-like input to an integer and adds thousands separators. Used for species counts, record counts and matched-record counts.
#'
#' @param x Numeric-like value to format.
#'
#' @return A character scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_fmt_int <- function(x) {
  x <- bf_console_num1(x)

  if (is.na(x)) return("NA")

  format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
}


#' Format a percentage for console messages
#'
#' Safely formats `num / den` as a percentage. Returns `"NA"` when either value is missing or the denominator is zero.
#'
#' @param num Numerator.
#' @param den Denominator.
#'
#' @return A character scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_pct <- function(num, den) {
  num <- bf_console_num1(num)
  den <- bf_console_num1(den)

  if (is.na(num) || is.na(den) || den == 0) {
    return("NA")
  }

  sprintf("%.1f%%", 100 * num / den)
}


#' Check whether cli is available
#'
#' Tests whether the optional {cli} package is installed so console helpers can use styled output when possible and base output otherwise.
#'
#' @return Logical scalar.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_has_cli <- function() {
  requireNamespace("cli", quietly = TRUE)
}


#' Print a generic pipeline message
#'
#' Writes a message using {cli} when available, otherwise falls back to `message()`. This is the lowest-level message helper used by other console functions.
#'
#' @param ... Objects pasted together to form the message text.
#' @param quiet Logical; if `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_msg <- function(..., quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  txt <- bf_console_text(paste0(...))

  if (bf_console_has_cli()) {
    tryCatch(
      cli::cli_inform(txt),
      error = function(e) message(txt)
    )
  } else {
    message(txt)
  }

  invisible(NULL)
}


#' Print a section rule
#'
#' Prints a visual section divider used at pipeline starts, species starts and pipeline completion. Falls back to a simple dashed-line divider without {cli}.
#'
#' @param title Optional section title.
#' @param quiet Logical; if `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_rule <- function(title = NULL, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  title <- bf_console_text(title, fallback = "")

  if (bf_console_has_cli()) {
    tryCatch(
      cli::cli_rule(title),
      error = function(e) {
        message("")
        message("-------------------------------------------------------------------------------")
        if (nzchar(title)) message(title)
        message("-------------------------------------------------------------------------------")
      }
    )
  } else {
    message("")
    message("-------------------------------------------------------------------------------")
    if (nzchar(title)) message(title)
    message("-------------------------------------------------------------------------------")
  }

  invisible(NULL)
}


#' Print an informational bullet
#'
#' Prints a compact bullet for routine pipeline status updates such as queued species, selected overlays, cleaning state and thinning state.
#'
#' @param text Message text.
#' @param quiet Logical; if `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_bullet <- function(text, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  text <- bf_console_text(text)

  if (bf_console_has_cli()) {
    tryCatch(
      cli::cli_bullets(c(">" = text)),
      error = function(e) message("  - ", text)
    )
  } else {
    message("  - ", text)
  }

  invisible(NULL)
}


#' Print a successful-step message
#'
#' Prints a success-style message for completed steps such as GBIF imports, overlay matches and finished species.
#'
#' @param text Message text.
#' @param quiet Logical; if `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_tick <- function(text, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  text <- bf_console_text(text)

  if (bf_console_has_cli()) {
    tryCatch(
      cli::cli_bullets(c("v" = text)),
      error = function(e) message("  [OK] ", text)
    )
  } else {
    message("  [OK] ", text)
  }

  invisible(NULL)
}


#' Print a warning-style pipeline message
#'
#' Prints a warning-style message without stopping execution. Used for non-fatal issues such as zero overlay matches.
#'
#' @param text Message text.
#' @param quiet Logical; if `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_warn <- function(text, quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  text <- bf_console_text(text)

  if (bf_console_has_cli()) {
    tryCatch(
      cli::cli_bullets(c("!" = text)),
      error = function(e) warning(text, call. = FALSE)
    )
  } else {
    warning(text, call. = FALSE)
  }

  invisible(NULL)
}

# -----------------------------------------------------------------------------
# Safe record-count helpers
# -----------------------------------------------------------------------------

#' Count records safely across common pipeline objects
#'
#' Returns row counts for data frames and `sf` objects, and recursively sums counts for lists. This avoids repeated object-type checks inside the pipelines.
#'
#' @param x Object to count; commonly a data frame, `sf` object, list of objects or `NULL`.
#'
#' @return Integer count.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_nrow <- function(x) {
  if (is.null(x)) return(0L)

  if (inherits(x, "sf") || is.data.frame(x)) {
    return(nrow(x))
  }

  if (is.list(x)) {
    return(sum(vapply(x, bf_console_nrow, integer(1)), na.rm = TRUE))
  }

  0L
}


#' Count non-missing overlay matches
#'
#' Counts records with a non-missing value in a spatial-overlay identifier column.
#'
#' @param x Data frame containing overlay assignment columns.
#' @param id_col Name of the overlay identifier column to inspect.
#'
#' @return Integer count of non-missing matches.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_count_matched <- function(x, id_col = NULL) {
  if (is.null(x) || !is.data.frame(x) || is.null(id_col) || !id_col %in% names(x)) {
    return(0L)
  }

  sum(!is.na(x[[id_col]]))
}


# -----------------------------------------------------------------------------
# Pipeline-level messages
# -----------------------------------------------------------------------------

#' Print the start message for a GBIF pipeline run
#'
#' Reports the pipeline type, number of queued species, selected terrestrial/freshwater region source or marine overlay, and whether coordinate cleaning and spatial thinning are enabled.
#'
#' @param workflow Pipeline label, such as `"terrestrial"`, `"marine"`, or `"terrestrial/freshwater"`.
#' @param n_species Number of queued species.
#' @param region_source Optional terrestrial/freshwater spatial overlay source.
#' @param overlay Optional marine overlay name.
#' @param gadm_unit Optional GADM level when `region_source = "gadm"`.
#' @param overlay_mode Optional overlay mode label.
#' @param cleaning Optional logical indicating whether coordinate cleaning is enabled.
#' @param thinning Optional logical indicating whether spatial thinning is enabled.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_pipeline_start <- function(workflow,
                                      n_species,
                                      region_source = NULL,
                                      overlay = NULL,
                                      gadm_unit = NULL,
                                      overlay_mode = NULL,
                                      cleaning = NULL,
                                      thinning = NULL,
                                      quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  workflow_txt <- bf_console_text(workflow)

  title <- switch(
    workflow_txt,
    terrestrial = "biofetchR terrestrial GBIF pipeline",
    marine = "biofetchR marine GBIF pipeline",
    paste("biofetchR", workflow_txt, "GBIF pipeline")
  )

  bf_console_rule(title, quiet = quiet)

  bf_console_bullet(
    paste0("Species queued: ", bf_console_fmt_int(n_species)),
    quiet = quiet
  )

  if (!is.null(region_source)) {
    region_source_txt <- bf_console_text(region_source)

    msg <- paste0("Region source: ", region_source_txt)

    if (!is.null(gadm_unit) && identical(region_source_txt, "gadm")) {
      msg <- paste0(msg, " | GADM level: ", bf_console_text(gadm_unit))
    }

    if (!is.null(overlay_mode)) {
      msg <- paste0(msg, " | overlay mode: ", bf_console_text(overlay_mode))
    }

    bf_console_bullet(msg, quiet = quiet)
  }

  if (!is.null(overlay)) {
    bf_console_bullet(
      paste0("Marine overlay: ", bf_console_text(overlay)),
      quiet = quiet
    )
  }

  if (!is.null(cleaning)) {
    bf_console_bullet(
      paste0("Coordinate cleaning: ", if (isTRUE(cleaning)) "enabled" else "skipped"),
      quiet = quiet
    )
  }

  if (!is.null(thinning)) {
    bf_console_bullet(
      paste0("Spatial thinning: ", if (isTRUE(thinning)) "enabled" else "skipped"),
      quiet = quiet
    )
  }

  invisible(NULL)
}


#' Print the completion message for a GBIF pipeline run
#'
#' Reports the number of completed species, number of output records and total elapsed time.
#'
#' @param workflow Pipeline label.
#' @param n_species Number of processed species.
#' @param n_records_out Number of output records generated by the pipeline.
#' @param start_time Start timestamp used to calculate elapsed time.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_pipeline_end <- function(workflow,
                                    n_species,
                                    n_records_out,
                                    start_time,
                                    quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  elapsed <- bf_console_elapsed(start_time)

  bf_console_rule("Pipeline finished", quiet = quiet)

  bf_console_tick(
    paste0(
      "Completed ",
      bf_console_fmt_int(n_species),
      " species with ",
      bf_console_fmt_int(n_records_out),
      " output records in ",
      bf_console_fmt_time(elapsed),
      "."
    ),
    quiet = quiet
  )

  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Species-level messages
# -----------------------------------------------------------------------------

#' Print a species-level start message
#'
#' Creates a rule heading for the species currently being processed, optionally including its index within the run and the selected region/overlay label.
#'
#' @param species Species name displayed in the heading.
#' @param index Optional current species index.
#' @param total Optional total number of species.
#' @param region Optional region or overlay label.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_species_start <- function(species,
                                     index = NULL,
                                     total = NULL,
                                     region = NULL,
                                     quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  prefix <- ""

  if (!is.null(index) && !is.null(total)) {
    prefix <- paste0("[", bf_console_fmt_int(index), "/", bf_console_fmt_int(total), "] ")
  }

  species <- bf_console_text(species)
  region <- bf_console_text(region, fallback = "")

  msg <- paste0(prefix, species)

  if (nzchar(region)) {
    msg <- paste0(msg, " | ", region)
  }

  bf_console_rule(msg, quiet = quiet)

  invisible(NULL)
}


#' Report a submitted GBIF download request
#'
#' Prints the species name, optional GBIF download key and optional elapsed time after a GBIF request has been submitted. Long key strings are truncated for readable console output.
#'
#' @param species Species name.
#' @param elapsed Optional elapsed seconds for the request step.
#' @param gbif_key Optional GBIF download key or vector/list of keys.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_gbif_request <- function(species,
                                    elapsed = NULL,
                                    gbif_key = NULL,
                                    quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  species <- bf_console_text(species)

  msg <- paste0("GBIF request submitted for ", species)

  if (!is.null(gbif_key) && length(gbif_key) > 0L) {
    key_vec <- bf_console_text_vec(gbif_key)
    key_txt <- paste(unique(key_vec), collapse = ", ")
    key_txt <- bf_console_text(key_txt)

    key_n <- suppressWarnings(
      nchar(key_txt, type = "chars", allowNA = FALSE, keepNA = FALSE)
    )

    if (length(key_n) == 0L || is.na(key_n)) key_n <- 0L

    if (key_n > 80) {
      key_txt <- paste0(substr(key_txt, 1, 77), "...")
    }

    msg <- paste0(msg, " | key: ", key_txt)
  }

  if (!is.null(elapsed)) {
    msg <- paste0(msg, " | ", bf_console_fmt_time(elapsed))
  }

  bf_console_tick(msg, quiet = quiet)

  invisible(NULL)
}


#' Report a completed GBIF import
#'
#' Prints the number of imported GBIF records and optional elapsed import time.
#'
#' @param n_records Number of imported records.
#' @param elapsed Optional elapsed seconds for the import step.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_gbif_import <- function(n_records,
                                   elapsed = NULL,
                                   quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  msg <- paste0("GBIF import completed: ", bf_console_fmt_int(n_records), " records")

  if (!is.null(elapsed)) {
    msg <- paste0(msg, " in ", bf_console_fmt_time(elapsed))
  }

  bf_console_tick(msg, quiet = quiet)

  invisible(NULL)
}


#' Report coordinate-cleaning retention
#'
#' Prints before/after row counts and retention percentage for the coordinate-cleaning step.
#'
#' @param n_before Number of records before cleaning.
#' @param n_after Number of records after cleaning.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_cleaning <- function(n_before,
                                n_after,
                                quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  bf_console_bullet(
    paste0(
      "Cleaning: ",
      bf_console_fmt_int(n_before),
      " -> ",
      bf_console_fmt_int(n_after),
      " retained",
      " (",
      bf_console_pct(n_after, n_before),
      ")"
    ),
    quiet = quiet
  )

  invisible(NULL)
}


#' Report spatial-thinning retention
#'
#' Prints before/after row counts and retention percentage for the spatial-thinning step, or a skipped message when thinning was not applied.
#'
#' @param n_before Number of records before thinning.
#' @param n_after Number of records after thinning.
#' @param applied Logical; whether thinning was applied.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_thinning <- function(n_before,
                                n_after,
                                applied = TRUE,
                                quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  if (!isTRUE(applied)) {
    bf_console_bullet("Spatial thinning: skipped", quiet = quiet)
    return(invisible(NULL))
  }

  bf_console_bullet(
    paste0(
      "Spatial thinning: ",
      bf_console_fmt_int(n_before),
      " -> ",
      bf_console_fmt_int(n_after),
      " retained",
      " (",
      bf_console_pct(n_after, n_before),
      ")"
    ),
    quiet = quiet
  )

  invisible(NULL)
}


#' Report the spatial overlay being applied
#'
#' Prints the selected terrestrial/freshwater region source or marine overlay before overlay assignment begins.
#'
#' @param region_source Optional terrestrial/freshwater overlay source.
#' @param overlay Optional marine overlay source.
#' @param gadm_unit Optional GADM level when `region_source = "gadm"`.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_overlay_start <- function(region_source = NULL,
                                     overlay = NULL,
                                     gadm_unit = NULL,
                                     quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  if (!is.null(region_source)) {
    region_source_txt <- bf_console_text(region_source)

    msg <- paste0("Spatial overlay: ", region_source_txt)

    if (identical(region_source_txt, "gadm") && !is.null(gadm_unit)) {
      msg <- paste0(msg, " level ", bf_console_text(gadm_unit))
    }

    bf_console_bullet(msg, quiet = quiet)
  }

  if (!is.null(overlay)) {
    bf_console_bullet(
      paste0("Marine overlay: ", bf_console_text(overlay)),
      quiet = quiet
    )
  }

  invisible(NULL)
}


#' Report spatial-overlay matching success
#'
#' Prints matched/total record counts and a match percentage after spatial-overlay assignment. A warning-style message is used when records were present but no overlay matches were found.
#'
#' @param n_records Total number of records considered for overlay matching.
#' @param n_matched Number of records with a non-missing overlay match.
#' @param id_col Optional overlay identifier column used for matching.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_overlay_result <- function(n_records,
                                      n_matched,
                                      id_col = NULL,
                                      quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  msg <- paste0(
    "Overlay match: ",
    bf_console_fmt_int(n_matched),
    " / ",
    bf_console_fmt_int(n_records),
    " records matched",
    " (",
    bf_console_pct(n_matched, n_records),
    ")"
  )

  if (!is.null(id_col)) {
    msg <- paste0(msg, " using ", bf_console_text(id_col))
  }

  n_records_num <- bf_console_num1(n_records)
  n_matched_num <- bf_console_num1(n_matched)

  if (!is.na(n_records_num) && !is.na(n_matched_num) &&
      n_records_num > 0L && n_matched_num == 0L) {
    bf_console_warn(msg, quiet = quiet)
  } else {
    bf_console_tick(msg, quiet = quiet)
  }

  invisible(NULL)
}


#' Print a species-level completion message
#'
#' Reports the final number of output records for a species and the elapsed time since the species step began.
#'
#' @param species Species name.
#' @param n_final Number of output records for this species.
#' @param start_time Start timestamp for this species step.
#' @param quiet Logical; if `TRUE`, suppress messages.
#'
#' @return Invisibly returns `NULL`.
#'
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_species_end <- function(species,
                                   n_final,
                                   start_time,
                                   quiet = FALSE) {
  if (isTRUE(quiet)) return(invisible(NULL))

  species <- bf_console_text(species)
  elapsed <- bf_console_elapsed(start_time)

  bf_console_tick(
    paste0(
      "Finished ",
      species,
      ": ",
      bf_console_fmt_int(n_final),
      " output records in ",
      bf_console_fmt_time(elapsed)
    ),
    quiet = quiet
  )

  invisible(NULL)
}


#' Return the expected terrestrial overlay identifier column
#'
#' Maps each supported terrestrial/freshwater overlay source to the identifier
#' column used when reporting overlay-match counts. This helper is intentionally
#' limited to reporting: it does not validate that an overlay was actually
#' requested or that a spatial join succeeded.
#'
#' @param region_source Terrestrial/freshwater overlay source name.
#'
#' @return A character scalar giving the expected identifier column. Unknown
#'   overlay names fall back to `"region_id"`.
#'
#' @family internal console helpers
#' @md
#' @keywords internal
#' @noRd
bf_console_expected_id_col_terrestrial <- function(region_source) {
  switch(
    as.character(region_source),
    gadm = "id_gadm",
    teow = "geoname_teow",
    feow = "geoname_feow",
    lakes = "geoname_lakes",
    rivers = "hyriv_id",
    basins = "hybas_id",
    gmba = "gmba_id",
    ne_urban = "geoname_ne_urban",
    resolve2017 = "geoname_resolve2017",
    ne_admin1 = "id_ne_admin1",
    wdpa = "wdpa_id",
    ramsar = "ramsar_id",
    gdw_barriers = "gdw_barrier_id",
    biosphere_reserve = "biosphere_id",
    gloric = "gloric_id",
    hydrowaste = "hydrowaste_id",
    gdw_reservoirs = "gdw_reservoir_id",
    global_mining = "global_mining_id",
    "region_id"
  )
}
