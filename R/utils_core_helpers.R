################################################################################
# utils_core_helpers.R
# ------------------------------------------------------------------------------
# biofetchR: minimal shared helpers used across package modules
# ------------------------------------------------------------------------------
#
# PURPOSE
#   This script contains small, general-purpose utilities that are reused across
#   many parts of biofetchR. These helpers are deliberately dependency-light and
#   should remain stable because other scripts rely on them for common tasks such
#   as optional-value handling, quiet-aware messaging, runtime package checks and
#   compact error reporting.
#
# HELPER GROUPS
#   The script is organised into four core helper groups:
#
#     - null-coalescing helpers for optional arguments and list elements;
#     - quiet-aware internal messaging;
#     - runtime package availability checks;
#     - compact error-message formatting for logs and audit tables.
#
################################################################################

# ------------------------------------------------------------------------------
# Null-coalescing helpers
# ------------------------------------------------------------------------------
# These helpers provide a compact way to supply fallbacks for optional arguments
# and optional list fields. They only treat NULL and length-zero objects as
# missing; they intentionally do not treat NA or blank strings as missing.

#' Return a fallback when an object is NULL or length zero
#'
#' Generic internal null-coalescing helper used throughout biofetchR. This is
#' appropriate for optional arguments, optional list fields and default values.
#' It does not treat `NA`, empty strings or blank strings as missing; use a
#' text-cleaning helper for those cases.
#'
#' @param x Primary object.
#' @param y Fallback object returned when `x` is `NULL` or length zero.
#'
#' @return `x` when it is non-`NULL` and has length greater than zero; otherwise
#'   `y`.
#'
#' @keywords internal
#' @noRd
bf_null_coalesce <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}


#' Null-coalescing operator
#'
#' Infix alias for [bf_null_coalesce()]. This is used internally for compact
#' optional-value handling.
#'
#' @param x Primary object.
#' @param y Fallback object.
#'
#' @return `x` unless it is `NULL` or length zero; otherwise `y`.
#'
#' @keywords internal
#' @noRd
`%||%` <- bf_null_coalesce


# ------------------------------------------------------------------------------
# Quiet-aware internal messaging
# ------------------------------------------------------------------------------
# Small helper scripts often need lightweight progress messages without adding
# extra dependencies or formatting assumptions. This wrapper centralises the
# common quiet = TRUE/FALSE behaviour.

#' Emit a quiet-aware internal biofetchR message
#'
#' Lightweight wrapper around [message()] used by helper scripts that need simple
#' progress messages without depending on {cli}. When `quiet = TRUE`, no message
#' is printed.
#'
#' @param ... Message components passed to [message()].
#' @param quiet Logical. If `TRUE`, suppress the message.
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
.bf_msg <- function(..., quiet = FALSE) {
  if (!isTRUE(quiet)) {
    message(...)
  }

  invisible(NULL)
}


# ------------------------------------------------------------------------------
# Runtime dependency checks
# ------------------------------------------------------------------------------
# Some biofetchR features depend on optional packages that should only be
# required when the relevant helper is used. This helper produces consistent,
# readable errors when those optional dependencies are missing.

#' Require optional packages at runtime
#'
#' Checks that one or more optional packages are installed before a helper
#' continues. This centralises repeated `requireNamespace()` checks used across
#' provider-specific and pipeline-specific helper scripts.
#'
#' @param pkgs Character vector of package names.
#' @param context Optional character label included at the start of the error
#'   message, for example `"native-range helper"` or `"marine pipeline"`.
#'
#' @return Invisibly returns `TRUE` when all packages are available. Throws an
#'   error listing missing packages otherwise.
#'
#' @keywords internal
#' @noRd
bf_require_packages <- function(pkgs, context = NULL) {
  pkgs <- unique(as.character(pkgs))
  pkgs <- pkgs[!is.na(pkgs) & nzchar(pkgs)]

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing)) {
    prefix <- if (!is.null(context) && length(context) && nzchar(context[[1]])) {
      paste0(as.character(context[[1]]), ": ")
    } else {
      ""
    }

    stop(
      prefix,
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Error-message compaction
# ------------------------------------------------------------------------------
# Remote services, downloads and spatial operations can return long or formatted
# errors. This helper makes those messages safe to store in summary CSVs,
# diagnostics and console audit outputs.

#' Compact an error message
#'
#' Converts an error condition or object to a single-line message by removing ANSI
#' control sequences, newlines and repeated whitespace. This keeps summary tables,
#' audit outputs and console messages readable when remote downloads, API calls or
#' spatial operations fail.
#'
#' @param e Error condition or object coercible to character.
#'
#' @return A single compact character string.
#'
#' @keywords internal
#' @noRd
bf_compact_error <- function(e) {
  msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

  msg <- gsub("\033\\[[0-9;]*m", "", msg)
  msg <- gsub("\n+", " ", msg)
  msg <- gsub("[[:space:]]+", " ", msg)

  trimws(msg)
}
