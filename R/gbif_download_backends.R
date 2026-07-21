################################################################################
# gbif_download_backends.R
# ------------------------------------------------------------------------------
# biofetchR GBIF download backends
# ------------------------------------------------------------------------------
#
# PURPOSE
#   Package-visible helpers for submitting GBIF occurrence downloads from the
#   terrestrial/freshwater and marine biofetchR pipelines.
#
# SCOPE
#   This file only submits GBIF download requests and returns download keys. It
#   does not wait for downloads to complete, import SIMPLE_CSV archives, clean
#   coordinates, thin occurrences, assign records to spatial overlays, or attach
#   native-range evidence. Those operations are handled by downstream helpers.
#
# CORE RESPONSIBILITIES
#   1. Submit species-level GBIF downloads for marine/global workflows.
#   2. Submit species x ISO2-country GBIF downloads for terrestrial/freshwater
#      workflows.
#   3. Extract download keys robustly across object shapes returned by
#      rgbif::occ_download().
#   4. Provide a compatibility wrapper for import code that expects a
#      wait_and_import_gbif_safe() helper.
#
# DESIGN NOTES
#   - Marine downloads are intentionally not country-filtered because offshore
#     records are assigned to EEZs or other marine overlays downstream.
#   - Terrestrial/freshwater downloads are intentionally country-filtered because
#     pipeline inputs are organised as species x recipient-country combinations.
#   - Download keys should be retained in pipeline outputs because they are the
#     link to GBIF download metadata and citation DOIs.
#
# DATA USE AND ATTRIBUTION
#   GBIF-mediated occurrence data are free to access but remain governed by the
#   GBIF data user agreement and by dataset-level licences from the original
#   publishers. biofetchR submits downloads but does not remove users'
#   responsibility to cite GBIF download DOIs and respect dataset licence terms.
#
################################################################################


#' Extract a GBIF download key from an rgbif return object
#'
#' Recover a GBIF download key from the different object shapes that
#' [rgbif::occ_download()] may return across `rgbif` versions, GBIF API
#' responses, or calling contexts.
#'
#' The helper is deliberately permissive because download submission backends
#' should fail only when no usable key can be recovered. It tries, in order:
#'
#' \enumerate{
#'   \item a non-empty character vector;
#'   \item a list-like object with a `key` element;
#'   \item a data frame with a `key` column;
#'   \item a `key` attribute;
#'   \item a best-effort character coercion.
#' }
#'
#' @param x Object returned by [rgbif::occ_download()], or any compatible object
#'   containing a GBIF download key.
#'
#' @return Character scalar. A GBIF download key when one can be extracted;
#'   otherwise `NA_character_`.
#'
#' @keywords internal
#' @family GBIF download backend internals
#' @noRd
.bf_extract_gbif_download_key <- function(x) {
  if (is.character(x) && length(x) >= 1L) {
    key <- x[!is.na(x) & nzchar(x)]
    if (length(key)) return(as.character(key[[1]]))
  }

  if (is.list(x) && !is.null(x[["key"]])) {
    key <- as.character(x[["key"]])
    key <- key[!is.na(key) & nzchar(key)]
    if (length(key)) return(key[[1]])
  }

  if (is.data.frame(x) && "key" %in% names(x) && nrow(x) >= 1L) {
    key <- as.character(x[["key"]])
    key <- key[!is.na(key) & nzchar(key)]
    if (length(key)) return(key[[1]])
  }

  key_attr <- attr(x, "key", exact = TRUE)
  if (!is.null(key_attr)) {
    key <- as.character(key_attr)
    key <- key[!is.na(key) & nzchar(key)]
    if (length(key)) return(key[[1]])
  }

  x_chr <- suppressWarnings(as.character(x))
  x_chr <- x_chr[!is.na(x_chr) & nzchar(x_chr)]
  if (length(x_chr)) return(x_chr[[1]])

  NA_character_
}


#' Submit global GBIF download jobs for species-level workflows
#'
#' Submit one GBIF occurrence download per species without applying a GBIF
#' country predicate. This backend is primarily used by the marine pipeline,
#' where occurrences are first retrieved globally and then assigned to marine
#' overlays such as EEZs, LMEs, MEOW ecoregions, IHO seas, FAO regions, or other
#' package-managed marine polygons.
#'
#' Each species name is resolved with [rgbif::name_backbone()] and submitted to
#' GBIF using a conservative occurrence predicate requiring:
#'
#' \itemize{
#'   \item the resolved GBIF `taxonKey`;
#'   \item `hasCoordinate = TRUE`;
#'   \item `hasGeospatialIssue = FALSE`;
#'   \item `occurrenceStatus = "PRESENT"`.
#' }
#'
#' Additional GBIF predicates can be appended through `predicate_extra`, allowing
#' callers to add date, basis-of-record, dataset, licence or other GBIF filters
#' without changing the backend.
#'
#' @section GBIF data use and citation:
#' This function submits live GBIF occurrence downloads and returns GBIF download
#' keys. GBIF-mediated data are free to access, but users who use GBIF downloads
#' in research, reporting or policy should cite the DOI generated for each GBIF
#' download, or an appropriate derived dataset DOI where applicable. Individual
#' records may also originate from datasets with their own licences and citation
#' requirements. biofetchR does not redistribute GBIF data or replace GBIF's
#' data-use and citation requirements.
#'
#' @section Side effects:
#' Submits live GBIF occurrence download jobs through [rgbif::occ_download()].
#' The function returns download keys only; it does not wait for, import, clean,
#' thin, or spatially join downloaded occurrence records.
#'
#' @section Errors:
#' The function stops if no valid species names are supplied, if GBIF credentials
#' are missing, if `predicate_extra` is not a list, if a species cannot be
#' resolved to a GBIF usage key, or if GBIF does not return an extractable
#' download key.
#'
#' @param batch Optional data frame containing a `species` column. If supplied,
#'   unique non-empty species names are taken from this column.
#' @param user GBIF username. Required for [rgbif::occ_download()].
#' @param pwd GBIF password. Required for [rgbif::occ_download()].
#' @param email GBIF account email. Required for [rgbif::occ_download()].
#' @param species Optional character vector of species names. Used only when
#'   `batch` is not supplied or does not contain a `species` column.
#' @param predicate_extra Optional list of additional `rgbif` predicates to
#'   append to each download request. Defaults to an empty list.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#' @param ... Additional arguments reserved for compatibility with pipeline code
#'   and future backend extensions. Currently ignored.
#'
#' @return A named character vector of GBIF occurrence-download keys. Each value
#'   is the download key returned by GBIF for one successfully submitted
#'   species-level download, and each name is the corresponding submitted
#'   species name. The returned keys are used by downstream import helpers and
#'   should be retained in summaries or manifests as the link to GBIF download
#'   metadata and citation DOI(s). The function stops if no usable download key
#'   can be extracted for a submitted species.
#'
#' @examples
#' \donttest{
#' gbif_env <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
#'
#' if (interactive() && all(nzchar(Sys.getenv(gbif_env)))) {
#'   keys <- download_gbif_batch(
#'     species = c("Ficopomatus enigmaticus", "Carcinus maenas"),
#'     user = Sys.getenv("GBIF_USER"),
#'     pwd = Sys.getenv("GBIF_PWD"),
#'     email = Sys.getenv("GBIF_EMAIL")
#'   )
#' }
#' }
#'
#' @export
#' @family GBIF download backends
#' @seealso [download_gbif_batch_gadm()], [wait_and_import_gbif_safe()]
#' @importFrom rgbif name_backbone pred pred_and occ_download
#' @importFrom stats na.omit
#' @md
download_gbif_batch <- function(batch = NULL,
                                user,
                                pwd,
                                email,
                                species = NULL,
                                predicate_extra = list(),
                                quiet = FALSE,
                                ...) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package 'rgbif' is required for download_gbif_batch().", call. = FALSE)
  }

  species_vec <- NULL

  if (!is.null(batch) && is.data.frame(batch) && "species" %in% names(batch)) {
    species_vec <- unique(trimws(as.character(batch$species)))
  } else if (!is.null(species)) {
    species_vec <- unique(trimws(as.character(species)))
  }

  species_vec <- species_vec[!is.na(species_vec) & nzchar(species_vec)]

  if (!length(species_vec)) {
    stop("No valid species names supplied to download_gbif_batch().", call. = FALSE)
  }

  if (missing(user) || is.null(user) || !nzchar(as.character(user))) {
    stop("GBIF username `user` is required.", call. = FALSE)
  }

  if (missing(pwd) || is.null(pwd) || !nzchar(as.character(pwd))) {
    stop("GBIF password `pwd` is required.", call. = FALSE)
  }

  if (missing(email) || is.null(email) || !nzchar(as.character(email))) {
    stop("GBIF account email `email` is required.", call. = FALSE)
  }

  if (!is.list(predicate_extra)) {
    stop("`predicate_extra` must be a list of rgbif predicates.", call. = FALSE)
  }

  out <- character(length(species_vec))
  names(out) <- species_vec

  for (sp in species_vec) {
    if (!isTRUE(quiet)) {
      message("Submitting GBIF marine download: ", sp)
    }

    backbone <- tryCatch(
      rgbif::name_backbone(name = sp),
      error = function(e) NULL
    )

    taxon_key <- NULL
    if (!is.null(backbone) && "usageKey" %in% names(backbone)) {
      taxon_key <- backbone$usageKey
    }

    if (is.null(taxon_key) || is.na(taxon_key) || !nzchar(as.character(taxon_key))) {
      stop(
        "Could not resolve a GBIF usageKey for species: ", sp,
        call. = FALSE
      )
    }

    predicates <- c(
      list(
        rgbif::pred("taxonKey", as.integer(taxon_key)),
        rgbif::pred("hasCoordinate", TRUE),
        rgbif::pred("hasGeospatialIssue", FALSE),
        rgbif::pred("occurrenceStatus", "PRESENT")
      ),
      predicate_extra
    )

    dl <- tryCatch(
      do.call(
        rgbif::occ_download,
        c(
          predicates,
          list(
            format = "SIMPLE_CSV",
            user = user,
            pwd = pwd,
            email = email
          )
        )
      ),
      error = function(e) {
        stop(
          "GBIF download submission failed for ", sp, ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    # rgbif::occ_download() has returned different object shapes across
    # versions: sometimes a list/data frame with `key`, sometimes an atomic
    # character vector. Never use `dl$key` before confirming `dl` is list-like.
    key <- .bf_extract_gbif_download_key(dl)

    if (is.na(key) || !nzchar(key)) {
      stop(
        "GBIF download was submitted for ", sp,
        " but no download key could be extracted.",
        call. = FALSE
      )
    }

    out[[sp]] <- key

    if (!isTRUE(quiet)) {
      message("GBIF marine download key for ", sp, ": ", key)
    }
  }

  out
}



#' Import GBIF downloads through the package import helper
#'
#' Compatibility wrapper around `wait_and_import_gbif()`. It exposes a safer,
#' explicit function name expected by some package scripts while preserving the
#' existing import behaviour implemented elsewhere in biofetchR.
#'
#' This wrapper is intentionally minimal: it checks that
#' `wait_and_import_gbif()` is visible, then forwards `download_keys` and any
#' additional arguments to that function.
#'
#' @param download_keys GBIF download keys, usually returned by
#'   [download_gbif_batch()] or [download_gbif_batch_gadm()]. The accepted object
#'   shape is whatever `wait_and_import_gbif()` supports.
#' @param ... Additional arguments passed directly to `wait_and_import_gbif()`,
#'   such as polling or retry settings if supported by that function.
#'
#' @return The object returned by [wait_and_import_gbif()]. In normal biofetchR
#'   workflows this is a named list of imported GBIF occurrence tables or `sf`
#'   point objects, with names corresponding to the supplied download-key labels.
#'   The exact structure depends on the available `wait_and_import_gbif()`
#'   implementation. This wrapper does not alter the imported records; it is
#'   called primarily to provide a stable package-visible import helper name for
#'   pipelines and tests.
#'
#' @section GBIF data use and citation:
#' The imported records should remain traceable to the GBIF download keys and
#' associated GBIF download DOI(s). Downstream workflows should preserve these
#' keys in manifests or summaries so users can cite the exact occurrence
#' downloads used in analysis.
#'
#' @section Errors:
#' Stops if `wait_and_import_gbif()` is not available in the loaded package
#' namespace or search path.
#'
#' @examples
#' \donttest{
#' if (interactive() && exists("keys")) {
#'   imported <- wait_and_import_gbif_safe(download_keys = keys)
#' }
#' }
#'
#' @export
#' @family GBIF download backends
#' @seealso [download_gbif_batch()], [download_gbif_batch_gadm()]
#' @md
wait_and_import_gbif_safe <- function(download_keys, ...) {
  if (!exists("wait_and_import_gbif", mode = "function", inherits = TRUE)) {
    stop(
      "wait_and_import_gbif_safe() requires wait_and_import_gbif(), but it is not available.",
      call. = FALSE
    )
  }

  wait_and_import_gbif(download_keys, ...)
}


#' Submit country-filtered GBIF downloads for terrestrial/freshwater workflows
#'
#' Submit one GBIF occurrence download per ISO2 country for a single species.
#' This backend is used by the terrestrial/freshwater pipeline, where inputs are
#' organised as species x country combinations and GBIF records should be
#' restricted to recipient countries before spatial overlay assignment.
#'
#' The function resolves the species to a GBIF taxon key, then submits one
#' download per country using the predicates:
#'
#' \itemize{
#'   \item `taxonKey = <resolved GBIF taxon key>`;
#'   \item `hasCoordinate = TRUE`;
#'   \item `country = <ISO2 country code>`.
#' }
#'
#' Download submission is retried up to three times per country. Successful
#' submissions are returned as a named list, with names equal to the ISO2 country
#' codes.
#'
#' @section GBIF data use and citation:
#' This function submits live GBIF occurrence downloads and returns GBIF download
#' keys. Users should retain those keys and cite the GBIF download DOI(s) for
#' occurrence records used in research, reports or policy outputs. Because GBIF
#' aggregates records from many data publishers, individual datasets may carry
#' their own licence and attribution requirements. biofetchR does not
#' redistribute GBIF records and does not replace those obligations.
#'
#' @section Side effects:
#' Submits live GBIF occurrence download jobs through [rgbif::occ_download()]
#' and prints submission/retry messages through `cli`.
#'
#' @section Notes:
#' If a package-visible `get_taxon_key()` helper exists, it is used for GBIF
#' name resolution. Otherwise the function falls back to
#' [rgbif::name_backbone()] directly. Download-key extraction is delegated to
#' `.bf_extract_gbif_download_key()` when available.
#'
#' @param species Character scalar. Scientific name of the species.
#' @param iso2_codes Character vector of ISO2 country codes. Codes are trimmed,
#'   upper-cased and deduplicated before submission.
#' @param user GBIF username. Required for [rgbif::occ_download()].
#' @param pwd GBIF password. Required for [rgbif::occ_download()].
#' @param email GBIF account email. Required for [rgbif::occ_download()].
#'
#' @return A named list of GBIF occurrence-download keys, with one element for
#'   each ISO2 country code where download submission succeeds. List names are
#'   the submitted ISO2 country codes and values are the corresponding GBIF
#'   download keys. These keys are used by downstream import helpers and should
#'   be retained in workflow summaries or manifests as the link to GBIF download
#'   metadata and citation DOI(s). Returns `NULL` when no valid species or
#'   country input is supplied, when the species cannot be resolved to a GBIF
#'   taxon key, or when all country-level submissions fail.
#'
#' @examples
#' \donttest{
#' gbif_env <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
#'
#' if (interactive() && all(nzchar(Sys.getenv(gbif_env)))) {
#'   keys <- download_gbif_batch_gadm(
#'     species = "Hydrocotyle ranunculoides",
#'     iso2_codes = c("GB", "IE"),
#'     user = Sys.getenv("GBIF_USER"),
#'     pwd = Sys.getenv("GBIF_PWD"),
#'     email = Sys.getenv("GBIF_EMAIL")
#'   )
#' }
#' }
#'
#' @export
#' @family GBIF download backends
#' @seealso [download_gbif_batch()], [wait_and_import_gbif_safe()]
#' @importFrom cli cli_alert_warning cli_alert_success
#' @importFrom rgbif name_backbone pred pred_and occ_download
#' @md
download_gbif_batch_gadm <- function(species,
                                     iso2_codes,
                                     user,
                                     pwd,
                                     email) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package 'rgbif' is required for GBIF downloads.", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required for GBIF download messages.", call. = FALSE)
  }

  species <- trimws(as.character(species)[1])
  iso2_codes <- unique(toupper(trimws(as.character(iso2_codes))))
  iso2_codes <- iso2_codes[!is.na(iso2_codes) & nzchar(iso2_codes)]

  if (!nzchar(species)) {
    cli::cli_alert_warning("No species name supplied to download_gbif_batch_gadm().")
    return(NULL)
  }

  if (!length(iso2_codes)) {
    cli::cli_alert_warning("No ISO2 country codes supplied for {.val {species}}.")
    return(NULL)
  }

  taxon_key <- tryCatch(
    {
      if (exists("get_taxon_key", mode = "function", inherits = TRUE)) {
        get_taxon_key(species)
      } else {
        bb <- rgbif::name_backbone(name = species)
        if (!is.null(bb$usageKey)) bb$usageKey else NA_integer_
      }
    },
    error = function(e) NA_integer_
  )

  if (is.na(taxon_key)) {
    cli::cli_alert_warning("Could not get GBIF taxonKey for {.val {species}}.")
    return(NULL)
  }

  keys <- list()

  for (cc in iso2_codes) {
    attempt <- 1L
    success <- FALSE

    while (!isTRUE(success) && attempt <= 3L) {
      result <- tryCatch(
        {
          rgbif::occ_download(
            rgbif::pred_and(
              rgbif::pred("taxonKey", taxon_key),
              rgbif::pred("hasCoordinate", TRUE),
              rgbif::pred("country", cc)
            ),
            format = "SIMPLE_CSV",
            user = user,
            pwd = pwd,
            email = email
          )
        },
        error = function(e) e
      )

      key <- NULL

      if (!inherits(result, "condition")) {
        key <- tryCatch(
          {
            if (exists(".bf_extract_gbif_download_key", mode = "function", inherits = TRUE)) {
              .bf_extract_gbif_download_key(result)
            } else if (is.character(result) && length(result) >= 1L) {
              as.character(result[[1]])
            } else if (is.list(result) && !is.null(result$key)) {
              as.character(result$key)
            } else {
              NA_character_
            }
          },
          error = function(e) NA_character_
        )
      }

      if (!is.null(key) && length(key) == 1L && !is.na(key) && nzchar(key)) {
        keys[[cc]] <- key
        cli::cli_alert_success(
          "Submitted terrestrial GBIF download for {.val {species}} in {.val {cc}} | key: {.val {key}}"
        )
        success <- TRUE
      } else {
        msg <- if (inherits(result, "condition")) conditionMessage(result) else "no key returned"
        cli::cli_alert_warning(
          "GBIF download submission failed for {.val {species}} in {.val {cc}} on attempt {.val {attempt}}: {msg}"
        )

        attempt <- attempt + 1L

        if (attempt <= 3L) {
          Sys.sleep(60)
        }
      }
    }
  }

  if (!length(keys)) {
    return(NULL)
  }

  keys
}

