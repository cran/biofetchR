################################################################################
# utils_encoding.R
# -----------------------------------------------------------------------------
# biofetchR: UTF-8 repair and defensive scalar-check helpers
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide small, dependency-free text-encoding utilities used by biofetchR
#   pipelines, overlay loaders, download helpers, summary writers and console
#   reporting functions.
#
# WHY THIS EXISTS
#   Biodiversity, administrative-boundary and spatial-overlay datasets often
#   contain mixed encodings. Non-ASCII names are common in GBIF metadata, GADM
#   regions, Marine Regions attributes, protected-area layers and ecoregion
#   tables. Examples include place names with accents or diacritics such as
#   Cura\\u00e7ao, S\\u00e3o Tom\\u00e9 and R\\u00e9union. If these strings are not repaired before
#   messaging, joining or CSV export, Windows and locale-specific R sessions can
#   fail unexpectedly or write corrupted labels.
#
# DESIGN PRINCIPLES
#   1. Preserve missing values as true NA values.
#   2. Convert character and factor fields to UTF-8 where possible.
#   3. Leave sf geometry/list columns untouched.
#   4. Provide one scalar-safe string helper for filenames, summaries and
#      messages.
#   5. Provide one strict scalar-TRUE helper for defensive condition checks.
#   6. Avoid external package dependencies so these helpers are safe to call from
#      low-level package infrastructure.
#
# PUBLIC API STATUS
#   These functions are intentionally internal package infrastructure. They use
#   non-dot names because multiple package scripts call them directly, but they
#   are documented with @noRd so roxygen records their purpose without generating
#   user-facing help pages.
#
# LICENSING / DATA ACCESS
#   No dataset-specific licence note is required for this file. These helpers do
#   not download, query, redistribute or transform any external provider data by
#   themselves; they only repair text already held in memory by other workflows.
#
################################################################################

#' Repair a vector to UTF-8-safe character text
#'
#' Converts an input vector to character and repairs it to UTF-8 using
#' [base::iconv()]. Missing values are restored after conversion so that `NA`
#' remains missing and is not accidentally written as the literal string
#' `"NA"`.
#'
#' This low-level helper is used before writing CSV outputs, constructing console
#' messages and handling external spatial-resource attributes that may contain
#' mixed encodings. It is intentionally conservative: it does not transliterate
#' characters, drop rows, change vector length or attempt dataset-specific name
#' cleaning.
#'
#' @param x Vector that can be coerced to character. `NULL` is returned unchanged.
#'
#' @return Character vector encoded as UTF-8 where possible, or `NULL` when
#'   `x` is `NULL`.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_repair_utf8_chr <- function(x) {
  if (is.null(x)) return(x)

  was_na <- is.na(x)

  y <- as.character(x)

  y <- iconv(
    y,
    from = "",
    to = "UTF-8",
    sub = "byte"
  )

  y[was_na] <- NA_character_

  ok <- !is.na(y)
  Encoding(y[ok]) <- "UTF-8"

  y
}


#' Repair text columns in a data frame-like object
#'
#' Applies [bf_repair_utf8_chr()] to character columns and to factor columns
#' after converting factor values to character. Geometry columns from `sf`
#' objects are deliberately skipped because they store spatial vectors rather
#' than ordinary text.
#'
#' This helper is useful immediately after reading external overlay resources,
#' after binding heterogeneous GBIF tables, and immediately before writing
#' tabular outputs. Non-data-frame inputs are returned unchanged so the function
#' can be used safely in defensive pipeline code.
#'
#' @param x Data frame, tibble or `sf` object. Non-data-frame inputs are returned
#'   unchanged.
#'
#' @return Object of the same general type as `x`, with character and factor
#'   columns repaired to UTF-8 where applicable.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_repair_utf8_df <- function(x) {
  if (!is.data.frame(x)) return(x)

  for (nm in names(x)) {
    if (inherits(x[[nm]], "sfc")) next

    if (is.character(x[[nm]])) {
      x[[nm]] <- bf_repair_utf8_chr(x[[nm]])
    }

    if (is.factor(x[[nm]])) {
      x[[nm]] <- factor(bf_repair_utf8_chr(as.character(x[[nm]])))
    }
  }

  x
}


#' Return one UTF-8-repaired scalar string
#'
#' Repairs an input vector to UTF-8 and returns the first element as a scalar
#' character value. If the input is empty or the first element is `NA`, the
#' supplied fallback is returned. This prevents console, filename and summary
#' helpers from failing when optional labels are missing.
#'
#' @param x Vector that can be coerced to character.
#' @param fallback Character scalar returned when `x` is empty or the first value
#'   is `NA`. Defaults to `""`.
#'
#' @return Length-one character value.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_repair_utf8_one <- function(x, fallback = "") {
  out <- bf_repair_utf8_chr(x)

  if (length(out) == 0L) return(fallback)
  if (is.na(out[[1]])) return(fallback)

  out[[1]]
}


#' Test for a strict scalar `TRUE` value
#'
#' Returns `TRUE` only when `x` is a length-one, non-missing value equal to
#' `TRUE`. This is a defensive helper for package internals where conditionals
#' must reject vectors, missing values, length-zero inputs and `NULL` without
#' triggering length-related warnings.
#'
#' @param x Object to test.
#'
#' @return Logical scalar: `TRUE` only for strict single-value `TRUE`; otherwise
#'   `FALSE`.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_is_true_scalar <- function(x) {
  # Check length before missingness. Calling is.na() on a vector and then using
  # && can throw "length > 1 in coercion to logical(1)". This helper must be safe
  # for vectors, NULL, length-zero inputs and NA because it is used in defensive
  # pipeline condition checks.
  if (is.null(x) || length(x) != 1L) return(FALSE)
  if (is.na(x)) return(FALSE)

  isTRUE(x)
}

#' Clean simple character text
#'
#' Repairs text to UTF-8, collapses repeated whitespace, trims leading/trailing
#' spaces and converts empty strings to `NA_character_`.
#'
#' This is a generic text-cleaning helper for provider metadata, country names,
#' status labels and API values. It is deliberately simpler than taxonomy-specific
#' name cleaning and should not replace taxonomy helpers such as `.bf_tax_squish()`.
#'
#' @param x Vector coercible to character.
#'
#' @return Character vector with UTF-8 repair, normalised whitespace and blank
#'   strings converted to `NA_character_`.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_clean_text <- function(x) {
  x <- bf_repair_utf8_chr(x)

  if (is.null(x)) {
    return(x)
  }

  x <- trimws(gsub("[[:space:]]+", " ", as.character(x)))
  x[x == ""] <- NA_character_

  x
}


#' Clean one scalar text value
#'
#' Scalar version of [bf_clean_text()]. Returns the first cleaned value, or a
#' fallback when the input is empty, missing or blank.
#'
#' @param x Object coercible to character.
#' @param fallback Value returned when `x` is empty, missing or blank.
#'
#' @return A length-one character value.
#'
#' @family internal encoding helpers
#' @keywords internal
#' @noRd
bf_clean_text_one <- function(x, fallback = NA_character_) {
  x <- bf_clean_text(x)

  if (is.null(x) || length(x) == 0L) {
    return(fallback)
  }

  x <- x[[1L]]

  if (is.na(x) || !nzchar(x)) {
    return(fallback)
  }

  x
}
