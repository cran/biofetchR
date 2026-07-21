################################################################################
# status_helpers.R
# -----------------------------------------------------------------------------
# biofetchR status and result-binding helpers
# -----------------------------------------------------------------------------
#
# PURPOSE
#   This file contains small package-level helpers used by biofetchR workflows
#   and summary/import utilities.
#
# MAIN HELPERS
#   - list_status_presets(): document native/non-native filtering presets.
#   - filter_by_status(): normalise provider-specific status labels and filter
#     records to one of the supported status presets.
#   - harmonize_column_types(): bind heterogeneous GBIF result chunks while
#     standardising sf geometry handling and avoiding mixed-column type failures.
#
# INTERNAL HELPERS
#   - .bf_normalize_status(): map raw labels to canonical status classes.
#   - .bf_resolve_status_keep(): resolve presets or custom status selections to
#     canonical classes.
#
# LICENSING / DATA USE
#   No external datasets are downloaded or redistributed by this file. Downstream
#   workflows that pass GBIF, GRIIS, Marine Regions, GADM or other provider data
#   into these helpers remain responsible for the licence, citation and
#   attribution requirements of those source datasets.
#
################################################################################

#' List available native/non-native status presets
#'
#' Return the status-filter presets recognised by [filter_by_status()]. These
#' presets define how raw origin, establishment or alien-status labels are grouped
#' into analysis-ready categories such as native, introduced, naturalized,
#' invasive or cryptogenic.
#'
#' @details
#' biofetchR uses these presets when users want consistent filtering across
#' provider-specific status fields. The presets are intentionally conservative:
#' `"non_native"` and its alias `"alien"` exclude cryptogenic records unless the
#' explicit `"non_native_plus_cryptogenic"` preset is requested.
#'
#' @return A base `data.frame` describing the status-filter presets recognised
#'   by [filter_by_status()]. The table contains three character columns:
#'   `preset`, the preset name accepted by `filter_by_status()`; `includes`, the
#'   canonical status classes retained by that preset; and `notes`, a short
#'   interpretation of the preset. Each row represents one supported preset or
#'   alias for filtering native, non-native, invasive, cryptogenic, managed or
#'   unknown status records.
#'
#' @examples
#' list_status_presets()
#'
#' @family status helpers
#' @export
list_status_presets <- function() {
  data.frame(
    preset = c(
      "all",
      "native",
      "non_native",         # alias for "alien"
      "alien",              # alias of non_native
      "introduced",
      "naturalized",
      "invasive",
      "native_plus_cryptogenic",
      "non_native_plus_cryptogenic"
    ),
    includes = c(
      "native, introduced, naturalized, invasive, cryptogenic, managed, unknown",
      "native",
      "introduced, naturalized, invasive",
      "introduced, naturalized, invasive",
      "introduced",
      "naturalized",
      "invasive",
      "native, cryptogenic",
      "introduced, naturalized, invasive, cryptogenic"
    ),
    notes = c(
      "No filtering",
      "Only taxa flagged as native/indigenous/endemic",
      "All alien/non-native categories (excl. cryptogenic)",
      "Alias for non_native",
      "Strictly introduced (often non-established)",
      "Established alien (naturalised)",
      "Alien flagged as invasive",
      "Treat cryptogenic with native (conservative)",
      "Treat cryptogenic with non-native (liberal)"
    ),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Internal mappers
# -----------------------------------------------------------------------------

#' Normalise raw origin/status labels to canonical classes
#'
#' Convert provider-specific native, alien, establishment and uncertainty labels
#' to the canonical classes used by [filter_by_status()]. The mapper uses broad
#' regular-expression matching so it can handle common variants such as
#' `"indigenous"`, `"non-native"`, `"naturalised"` and `"established"`.
#'
#' @param x Character vector of raw status labels.
#'
#' @return Character vector with values in `native`, `introduced`,
#'   `naturalized`, `invasive`, `cryptogenic`, `managed` or `unknown`.
#'
#' @keywords internal
#' @noRd
.bf_normalize_status <- function(x) {
  z <- tolower(trimws(as.character(x)))

  # Treat empty/NA as unknown early
  z[is.na(z) | z %in% c("", "na", "n/a", "none", "unspecified")] <- "unknown"

  # Collapse punctuation/underscores and condense spaces
  z <- gsub("[^a-z]+", " ", z)
  z <- gsub("\\s+", " ", z)

  # canonical outputs: native | introduced | naturalized | invasive | cryptogenic | managed | unknown
  res <- ifelse(grepl("\\b(native|indigenous|endemic|resident)\\b", z), "native",
                ifelse(grepl("\\b(invasive)\\b", z), "invasive",
                       ifelse(grepl("\\b(naturaliz|naturalised|naturalized|establish(ed|ment)?)\\b", z), "naturalized",
                              ifelse(grepl("\\b(alien|non ?native|exotic|introduced|adventive|casual)\\b", z), "introduced",
                                     ifelse(grepl("\\b(cryptogenic)\\b", z), "cryptogenic",
                                            ifelse(grepl("\\b(unknown|uncertain|undetermined|unresolved|unverified|not specified)\\b", z), "unknown",
                                                   ifelse(grepl("\\b(managed|cultivat|captive|controlled)\\b", z), "managed", "unknown")))))))

  res
}


#' Resolve a status preset into canonical classes to keep
#'
#' Internal companion to [filter_by_status()]. When `custom_values` is supplied,
#' the function intersects those values with the recognised canonical status
#' classes. Otherwise it resolves a named preset from [list_status_presets()].
#'
#' @param preset Character scalar preset name.
#' @param custom_values Optional character vector of canonical statuses to keep.
#'
#' @return Character vector of canonical status classes.
#'
#' @keywords internal
#' @noRd
.bf_resolve_status_keep <- function(preset, custom_values = NULL) {
  if (!is.null(custom_values)) {
    canon <- c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown")
    return(intersect(custom_values, canon))
  }
  preset <- tolower(preset)
  if (preset == "alien") preset <- "non_native"
  switch(preset,
         all                          = c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown"),
         native                       = "native",
         non_native                   = c("introduced","naturalized","invasive"),
         introduced                   = "introduced",
         naturalized                  = "naturalized",
         invasive                     = "invasive",
         native_plus_cryptogenic      = c("native","cryptogenic"),
         non_native_plus_cryptogenic  = c("introduced","naturalized","invasive","cryptogenic"),
         c("native","introduced","naturalized","invasive","cryptogenic","managed","unknown")
  )
}

#' Filter or normalise records by native/alien status
#'
#' Normalise a provider-specific origin, establishment or alien-status column to
#' biofetchR's canonical status classes and optionally filter rows to a supported
#' preset.
#'
#' @details
#' The function adds or updates a `status_norm` column using the internal status
#' mapper. When `normalize_only = TRUE`, no rows are removed. Otherwise, rows are
#' retained when `status_norm` belongs to the classes resolved from `preset` or
#' `values`.
#'
#' This helper is deliberately tolerant of missing status columns. If `status_col`
#' is absent, the input data are returned unchanged. A warning is emitted only
#' when a real filtering request was made.
#'
#' @param data A `data.frame` or tibble containing status information.
#' @param status_col Character scalar. Column containing raw status labels.
#'   Defaults to `"alien_status"`.
#' @param preset Character scalar. One of `list_status_presets()$preset`.
#'   Defaults to `"all"`, which retains all canonical classes.
#' @param values Optional character vector of canonical statuses to keep, such as
#'   `c("introduced", "naturalized")`. When supplied, this overrides `preset`.
#' @param normalize_only Logical. If `TRUE`, add `status_norm` but do not filter
#'   rows.
#' @param verbose Logical. If `TRUE`, print a short retained-row summary.
#'
#' @return A data frame or tibble of the same general type as `data`, with a
#'   `status_norm` column added or updated when `status_col` is present.
#'   `status_norm` contains canonical biofetchR status classes such as
#'   `"native"`, `"introduced"`, `"naturalized"`, `"invasive"`,
#'   `"cryptogenic"`, `"managed"` or `"unknown"`. If `normalize_only = TRUE`,
#'   all input rows are returned with the normalised status column. Otherwise,
#'   rows are filtered to the classes selected by `preset` or `values`. If
#'   `status_col` is absent, the input data are returned unchanged.
#'
#' @examples
#' x <- data.frame(
#'   species = c("Species a", "Species b", "Species c"),
#'   alien_status = c("native", "introduced", "cryptogenic")
#' )
#'
#' filter_by_status(x, preset = "non_native", verbose = FALSE)
#' filter_by_status(x, normalize_only = TRUE, verbose = FALSE)
#'
#' @family status helpers
#' @export
#' @md
filter_by_status <- function(
    data,
    status_col = "alien_status",
    preset = "all",
    values = NULL,
    normalize_only = FALSE,
    verbose = TRUE
) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame-like object.")
  if (!status_col %in% names(data)) {
    if (!normalize_only && (!identical(tolower(preset), "all") || !is.null(values))) {
      warning("Status column '", status_col, "' not found; skipping status filter.")
    }
    return(data)
  }

  data$status_norm <- .bf_normalize_status(data[[status_col]])

  if (normalize_only) return(data)

  keep_vals <- .bf_resolve_status_keep(preset, values)
  n0 <- nrow(data)
  data <- data[data$status_norm %in% keep_vals, , drop = FALSE]
  n1 <- nrow(data)

  if (verbose) {
    kept_tab <- sort(table(data$status_norm), decreasing = TRUE)
    msg <- sprintf("Status filter: kept %d / %d rows (preset='%s').", n1, n0, preset)
    message(msg)
    if (length(kept_tab)) utils::capture.output(kept_tab) |> paste(collapse = "; ") |> message()
  }

  if (n1 == 0L) warning("After status filtering, no rows remain.")
  data
}

#' Harmonise and bind GBIF result chunks
#'
#' Bind a list of GBIF result chunks into a single spatial object while
#' standardising geometry handling and reducing mixed-column type conflicts.
#'
#' @details
#' GBIF imports can produce chunks with heterogeneous columns, missing fields or
#' incompatible column classes across species, countries, downloads or publishing
#' datasets. `bf_bind_gbif_chunks()` keeps only non-empty chunks, converts plain
#' data frames with `decimalLongitude` and `decimalLatitude` columns into `sf`
#' point objects, standardises the geometry column name to `"geometry"`, and
#' transforms spatial inputs to WGS84 when needed.
#'
#' The helper then aligns the union of all non-geometry columns across valid
#' chunks and, by default, coerces those attributes to character before binding.
#' This conservative coercion avoids common row-binding failures caused by GBIF
#' fields switching between integer, numeric, logical, factor and character
#' classes across chunks.
#'
#' The function intentionally avoids relying on `sf::st_geometry_name()` so that
#' it remains compatible with older `sf` versions used on some Windows or HPC
#' installations.
#'
#' @param x_list List of GBIF result chunks. Each element should be an `sf` object
#'   or a data frame containing `decimalLongitude` and `decimalLatitude` columns.
#'   Empty, `NULL` or unsupported chunks are silently dropped.
#' @param coerce_all_to_character Logical. If `TRUE`, all non-geometry columns
#'   are coerced to character before binding to avoid type conflicts across GBIF
#'   chunks. If `FALSE`, existing column classes are retained where possible.
#'
#' @return An `sf` point object when at least one valid GBIF chunk is supplied.
#'   The returned object contains the union of non-geometry columns across valid
#'   chunks, a standardised `geometry` column, and WGS84 point geometries derived
#'   from existing `sf` geometries or from `decimalLongitude` and
#'   `decimalLatitude` columns. When `coerce_all_to_character = TRUE`, all
#'   non-geometry columns are coerced to character before binding to avoid mixed
#'   type failures across GBIF chunks. If no valid chunks are available, the
#'   function returns an empty tibble.
#'
#' @section Relationship to older helpers:
#' This function replaces the older, ambiguous `harmonize_column_types()` helper
#' for GBIF result chunks. It is named explicitly to show that it returns one
#' bound object, not a harmonised list.
#'
#' @family GBIF import helpers
#'
#' @examples
#' x1 <- data.frame(
#'   species = "Example species",
#'   decimalLongitude = -6.26,
#'   decimalLatitude = 53.35,
#'   individualCount = 1
#' )
#'
#' x2 <- data.frame(
#'   species = "Example species",
#'   decimalLongitude = -6.30,
#'   decimalLatitude = 53.40,
#'   basisOfRecord = "HUMAN_OBSERVATION"
#' )
#'
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   bf_bind_gbif_chunks(list(x1, x2))
#' }
#'
#' @md
#' @export
bf_bind_gbif_chunks <- function(x_list, coerce_all_to_character = TRUE) {
  # keep only non-empty elements
  x_list <- Filter(function(x) !is.null(x) && NROW(x) > 0, x_list)
  if (!length(x_list)) return(tibble::tibble())

  # helper: standardize to sf with geometry column "geometry"
  .to_sf_std <- function(x) {
    if (!inherits(x, "sf")) {
      if (all(c("decimalLongitude", "decimalLatitude") %in% names(x))) {
        x <- sf::st_as_sf(
          x, coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326, remove = FALSE
        )
      } else {
        # no coords, cannot make sf -> drop this chunk
        return(NULL)
      }
    }
    # ensure WGS84
    if (!is.na(sf::st_crs(x)) && sf::st_crs(x)$epsg != 4326) {
      x <- sf::st_transform(x, 4326)
    }
    # find the sfc column name from attributes/classes, not st_geometry_name()
    gcol <- attr(x, "sf_column", exact = TRUE)
    if (is.null(gcol) || !gcol %in% names(x)) {
      # locate first sfc column by class
      sfc_cols <- names(x)[vapply(x, function(col) inherits(col, "sfc"), logical(1))]
      if (length(sfc_cols)) {
        gcol <- sfc_cols[1]
      } else {
        return(NULL)  # no geometry found, drop
      }
    }
    # rename to "geometry" if needed and set sf attr
    if (!identical(gcol, "geometry")) {
      names(x)[names(x) == gcol] <- "geometry"
      attr(x, "sf_column") <- "geometry"
      sf::st_geometry(x) <- "geometry"
    }
    x
  }

  x_list <- lapply(x_list, .to_sf_std)
  x_list <- Filter(Negate(is.null), x_list)
  if (!length(x_list)) return(tibble::tibble())

  # union of attribute names (excluding geometry)
  name_union <- Reduce(union, lapply(x_list, function(x) setdiff(names(x), "geometry")))

  # align columns & coerce types
  x_list <- lapply(x_list, function(x) {
    missing_cols <- setdiff(name_union, names(x))
    if (length(missing_cols)) x[missing_cols] <- NA_character_
    x <- x[, c(name_union, "geometry"), drop = FALSE]
    if (isTRUE(coerce_all_to_character)) {
      for (nm in name_union) x[[nm]] <- as.character(x[[nm]])
    }
    x
  })

  # bind (sf rbind)
  out <- do.call(rbind, x_list)

  # if geometry types vary, try to cast to POINT (our GBIF outputs are points)
  try({
    gtyp <- unique(sf::st_geometry_type(out))
    if (length(gtyp) > 1) out <- sf::st_cast(out, "POINT", warn = FALSE)
  }, silent = TRUE)

  out
}

