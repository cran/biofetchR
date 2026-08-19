################################################################################
# gbif_download_backends.R
# ------------------------------------------------------------------------------
# biofetchR GBIF download backends
# ------------------------------------------------------------------------------
#
# PURPOSE
#   Package-visible helpers for submitting asynchronous GBIF occurrence downloads
#   from the terrestrial/freshwater and marine biofetchR pipelines.
#
# SCOPE
#   This file is responsible for GBIF download SUBMISSION and submission-level
#   queue management. It does not import completed SIMPLE_CSV archives, clean
#   coordinates, thin occurrences, assign records to spatial overlays, or attach
#   native-range evidence. Download polling and import are handled downstream by
#   wait_and_import_gbif().
#
# CORE RESPONSIBILITIES
#   1. Submit species-level global GBIF downloads for marine workflows.
#   2. Submit one multi-country GBIF download per terrestrial/freshwater species,
#      combining all retained recipient ISO2 countries in one country-IN
#      predicate.
#   3. Inspect the authenticated GBIF account for PREPARING/RUNNING occurrence
#      downloads before a new asynchronous request is submitted.
#   4. Retry the same request when GBIF reports that the account has reached its
#      current incomplete/concurrent-download limit.
#   5. Extract download keys robustly across object shapes returned by
#      rgbif::occ_download().
#   6. Attach submission-level audit information so failed requests do not erase
#      or obscure successful download keys from the same workflow.
#   7. Provide a compatibility wrapper for import code that expects a
#      wait_and_import_gbif_safe() helper.
#
# DESIGN NOTES
#   - Marine downloads are intentionally not country-filtered because offshore
#     records are assigned to EEZs or other marine overlays downstream.
#   - Terrestrial/freshwater species-country combinations remain the logical
#     request and audit units, but they are no longer necessarily separate GBIF
#     asynchronous download jobs. All retained countries for one species are
#     combined into one multi-country download.
#   - The default queue policy is deliberately conservative:
#     max_active_downloads = 1. Before submission, biofetchR checks the
#     authenticated GBIF account and waits while the configured number of
#     PREPARING/RUNNING jobs is already active.
#   - Account-level queue inspection is fail-soft. If the GBIF download list
#     cannot be inspected, submission is allowed to proceed and explicit GBIF
#     concurrency-limit responses are handled by the server-side retry guard.
#   - Submission failures are represented explicitly in submission audits rather
#     than causing successful keys from unrelated species or countries to be
#     silently lost.
#   - GBIF download keys should be retained in pipeline outputs because they are
#     the reproducibility link to GBIF download metadata and citation DOIs.
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


#' Detect a GBIF incomplete-download or concurrency-limit error
#'
#' Identify error messages indicating that GBIF has temporarily refused an
#' asynchronous occurrence-download submission because the authenticated account
#' has reached the number of incomplete or concurrent downloads currently allowed
#' by the service.
#'
#' @details
#' The helper recognises several message forms observed across GBIF/rgbif
#' responses, including references to simultaneous downloads, incomplete
#' downloads, download limitations, concurrent downloads, and GBIF's
#' `"enhance your calm"` rate-limit wording.
#'
#' This helper does not itself wait or retry. It only classifies the error so
#' `.bf_submit_gbif_with_retry()` can distinguish a temporary capacity condition
#' from a non-retryable submission error.
#'
#' @param x A condition object or character value containing a GBIF submission
#'   error message.
#'
#' @return Logical scalar. `TRUE` when the message appears to represent a GBIF
#'   incomplete/concurrent-download capacity limit; otherwise `FALSE`.
#'
#' @keywords internal
#' @family GBIF download backend internals
#' @noRd
.bf_is_gbif_download_limit_error <- function(x) {
  msg <- if (inherits(x, "condition")) conditionMessage(x) else paste(as.character(x), collapse = " ")
  msg <- tolower(msg)

  patterns <- c(
    "too many simultaneous downloads",
    "download limitation is exceeded",
    "too many incomplete downloads",
    "concurrent download",
    "enhance your calm"
  )

  any(vapply(patterns, function(pattern) grepl(pattern, msg, fixed = TRUE), logical(1)))
}


#' Count active GBIF occurrence downloads for the authenticated account
#'
#' Inspect the authenticated user's GBIF occurrence-download list and count
#' downloads whose status is currently `PREPARING` or `RUNNING`.
#'
#' @details
#' This helper supports proactive queue control before a new asynchronous
#' occurrence download is submitted. It is deliberately fail-soft: if the
#' account download list cannot be retrieved, or if no recognisable status column
#' is available, the function returns `NA_integer_` rather than stopping the
#' pipeline. The caller can then fall back to GBIF's server-side capacity response.
#'
#' `list_fun` exists primarily for dependency injection in package tests. Normal
#' production use relies on [rgbif::occ_download_list()].
#'
#' @param user GBIF username used to inspect the authenticated download list.
#' @param pwd GBIF password used to inspect the authenticated download list.
#' @param list_fun Optional replacement for [rgbif::occ_download_list()], used
#'   primarily by unit tests. If `NULL`, `rgbif::occ_download_list` is used.
#'
#' @return Integer scalar giving the number of downloads currently reported as
#'   `PREPARING` or `RUNNING`. Returns `0L` when a valid download list contains no
#'   active jobs, and `NA_integer_` when the active count cannot be determined.
#'
#' @keywords internal
#' @family GBIF download backend internals
#' @noRd
.bf_count_active_gbif_downloads <- function(user,
                                            pwd,
                                            list_fun = NULL) {
  if (is.null(list_fun)) list_fun <- rgbif::occ_download_list

  x <- tryCatch(
    list_fun(user = user, pwd = pwd, limit = 1000),
    error = function(e) NULL
  )

  if (is.null(x)) return(NA_integer_)

  dat <- if (is.list(x) && !is.null(x$results)) x$results else x
  if (!is.data.frame(dat) || !nrow(dat)) return(0L)

  status_col <- intersect(
    c("status", "request.status", "requestStatus"),
    names(dat)
  )

  if (!length(status_col)) return(NA_integer_)

  status <- toupper(trimws(as.character(dat[[status_col[[1L]]]])))
  as.integer(sum(status %in% c("PREPARING", "RUNNING"), na.rm = TRUE))
}


#' Wait until the authenticated GBIF account has submission capacity
#'
#' Poll the authenticated GBIF occurrence-download list until the number of
#' `PREPARING`/`RUNNING` jobs is below `max_active_downloads`, or until the
#' configured polling limit is exhausted.
#'
#' @details
#' biofetchR uses a deliberately conservative default of one active asynchronous
#' download. With `max_active_downloads = 1L`, a new request is submitted only
#' when the account-level active count is zero.
#'
#' If the account-level active count cannot be determined, this helper returns
#' immediately and allows submission to proceed. GBIF's own response is then
#' handled by `.bf_submit_gbif_with_retry()`.
#'
#' Exhausting `max_polls` does not itself throw an error. The helper returns
#' `FALSE`, allowing the caller to attempt submission and use the server-side
#' concurrency-limit retry guard if necessary.
#'
#' `list_fun` and `sleep_fun` are dependency-injection hooks used by tests to
#' simulate queue states and waiting without contacting GBIF or sleeping.
#'
#' @param user GBIF username.
#' @param pwd GBIF password.
#' @param max_active_downloads Integer. Maximum number of active
#'   `PREPARING`/`RUNNING` downloads permitted before another submission is
#'   attempted. Defaults to `1L`.
#' @param poll_seconds Numeric. Seconds between account-level queue checks.
#' @param max_polls Integer. Maximum number of account-level queue checks.
#' @param quiet Logical. If `TRUE`, suppress queue-wait progress messages.
#' @param list_fun Optional replacement for [rgbif::occ_download_list()] used by
#'   tests.
#' @param sleep_fun Optional replacement for [Sys.sleep()] used by tests.
#'
#' @return Invisibly returns `TRUE` when submission capacity appears available or
#'   the active count cannot be determined. Invisibly returns `FALSE` when the
#'   account remained at or above the configured active-download threshold for all
#'   permitted polls.
#'
#' @keywords internal
#' @family GBIF download backend internals
#' @noRd
.bf_wait_for_gbif_download_slot <- function(user,
                                            pwd,
                                            max_active_downloads = 1L,
                                            poll_seconds = 30,
                                            max_polls = 120L,
                                            quiet = FALSE,
                                            list_fun = NULL,
                                            sleep_fun = NULL) {
  if (is.null(sleep_fun)) sleep_fun <- Sys.sleep

  max_active_downloads <- max(1L, as.integer(max_active_downloads[[1L]]))
  poll_seconds <- max(0, as.numeric(poll_seconds[[1L]]))
  max_polls <- max(1L, as.integer(max_polls[[1L]]))

  # To keep at most `max_active_downloads` jobs active, a new submission is safe
  # only when the current active count is strictly below that value.
  for (i in seq_len(max_polls)) {
    n_active <- .bf_count_active_gbif_downloads(
      user = user,
      pwd = pwd,
      list_fun = list_fun
    )

    if (is.na(n_active) || n_active < max_active_downloads) {
      return(invisible(TRUE))
    }

    if (!isTRUE(quiet) && (i == 1L || i %% 10L == 0L)) {
      message(
        "GBIF account currently has ", n_active,
        " active download(s); waiting before the next submission."
      )
    }

    sleep_fun(poll_seconds)
  }

  invisible(FALSE)
}


#' Submit one GBIF download request with queue-aware retry handling
#'
#' Submit one asynchronous GBIF occurrence-download request after applying
#' proactive account-level queue control, and retry the same request when GBIF
#' explicitly reports that download capacity is full.
#'
#' @details
#' Before every submission attempt, the helper calls
#' `.bf_wait_for_gbif_download_slot()`. The request is then submitted as
#' `SIMPLE_CSV` through [rgbif::occ_download()].
#'
#' A successful response is reduced to its GBIF download key. If GBIF reports a
#' recognised incomplete/concurrent-download limit, the helper waits and retries
#' the SAME request so the requested species or countries are not silently
#' dropped.
#'
#' Non-capacity errors are returned immediately as `submit_failed`. A successful
#' GBIF response from which no usable download key can be extracted is returned as
#' `submit_no_key`. Exhaustion of capacity retries is returned as
#' `submit_timeout`.
#'
#' `submit_fun`, `list_fun` and `sleep_fun` are dependency-injection hooks for
#' unit tests and should normally be left `NULL` in production use.
#'
#' @param predicate_args List of `rgbif` predicate objects passed to
#'   [rgbif::occ_download()].
#' @param user GBIF username.
#' @param pwd GBIF password.
#' @param email GBIF account email.
#' @param label Human-readable request label used in progress/retry messages.
#' @param max_active_downloads Integer queue threshold passed to
#'   `.bf_wait_for_gbif_download_slot()`.
#' @param slot_poll_seconds Numeric seconds between queue checks and capacity
#'   retry sleeps.
#' @param slot_max_polls Integer maximum number of proactive queue checks before
#'   submission is attempted.
#' @param submission_max_tries Integer maximum number of attempts when GBIF keeps
#'   returning a recognised capacity-limit error.
#' @param quiet Logical. If `TRUE`, suppress retry progress messages.
#' @param submit_fun Optional replacement for [rgbif::occ_download()] used by
#'   tests.
#' @param list_fun Optional replacement for [rgbif::occ_download_list()] used by
#'   tests.
#' @param sleep_fun Optional replacement for [Sys.sleep()] used by tests.
#'
#' @return A list with four elements:
#'   `key`, containing the GBIF download key or `NA_character_`;
#'   `status`, containing one of `"submitted"`, `"submit_no_key"`,
#'   `"submit_failed"` or `"submit_timeout"`;
#'   `message`, containing an error/recovery message where relevant; and
#'   `attempts`, giving the submission-attempt count.
#'
#' @keywords internal
#' @family GBIF download backend internals
#' @noRd
.bf_submit_gbif_with_retry <- function(predicate_args,
                                       user,
                                       pwd,
                                       email,
                                       label = "GBIF download",
                                       max_active_downloads = 1L,
                                       slot_poll_seconds = 30,
                                       slot_max_polls = 120L,
                                       submission_max_tries = 120L,
                                       quiet = FALSE,
                                       submit_fun = NULL,
                                       list_fun = NULL,
                                       sleep_fun = NULL) {
  if (is.null(submit_fun)) submit_fun <- rgbif::occ_download
  if (is.null(sleep_fun)) sleep_fun <- Sys.sleep

  if (!is.list(predicate_args)) {
    stop("`predicate_args` must be a list of rgbif predicates.", call. = FALSE)
  }

  submission_max_tries <- max(1L, as.integer(submission_max_tries[[1L]]))

  for (attempt in seq_len(submission_max_tries)) {
    .bf_wait_for_gbif_download_slot(
      user = user,
      pwd = pwd,
      max_active_downloads = max_active_downloads,
      poll_seconds = slot_poll_seconds,
      max_polls = slot_max_polls,
      quiet = quiet,
      list_fun = list_fun,
      sleep_fun = sleep_fun
    )

    submitted <- tryCatch(
      do.call(
        submit_fun,
        c(
          predicate_args,
          list(
            format = "SIMPLE_CSV",
            user = user,
            pwd = pwd,
            email = email
          )
        )
      ),
      error = function(e) e
    )

    if (!inherits(submitted, "condition")) {
      key <- .bf_extract_gbif_download_key(submitted)
      if (!is.na(key) && nzchar(key)) {
        return(list(
          key = key,
          status = "submitted",
          message = NA_character_,
          attempts = attempt
        ))
      }

      return(list(
        key = NA_character_,
        status = "submit_no_key",
        message = "GBIF submission returned no extractable download key.",
        attempts = attempt
      ))
    }

    msg <- conditionMessage(submitted)

    if (.bf_is_gbif_download_limit_error(submitted) && attempt < submission_max_tries) {
      if (!isTRUE(quiet)) {
        message(
          label,
          ": GBIF download capacity is currently full; waiting and retrying the same request."
        )
      }
      sleep_fun(slot_poll_seconds)
      next
    }

    return(list(
      key = NA_character_,
      status = if (.bf_is_gbif_download_limit_error(submitted)) "submit_timeout" else "submit_failed",
      message = msg,
      attempts = attempt
    ))
  }

  list(
    key = NA_character_,
    status = "submit_timeout",
    message = "GBIF submission retry limit was exhausted.",
    attempts = submission_max_tries
  )
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
#'
#' @section Queue-aware submission:
#' Each species is submitted through the shared queue-aware submission helper.
#' Before each request, biofetchR checks the authenticated account for
#' `PREPARING`/`RUNNING` GBIF occurrence downloads. The default
#' `max_active_downloads = 1L` therefore avoids intentionally creating multiple
#' unfinished asynchronous downloads from this backend.
#'
#' If GBIF nevertheless reports that download capacity is full, the same species
#' request is retained and retried rather than being discarded.
#'
#' @section Submission audit and partial results:
#' The returned key vector carries a `submission_audit` attribute with one row
#' per requested species and the columns `label`, `gbif_key`,
#' `submission_status`, `message`, and `attempts`.
#'
#' Species-level failures do not remove keys already obtained for other species.
#' Consequently, the returned character vector can contain only the successfully
#' submitted species while the audit records failures such as
#' `taxon_key_failed`, `submit_no_key`, `submit_failed`, and `submit_timeout`.
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
#' @section Errors and failure handling:
#' The function stops for invalid function-level configuration, including no valid
#' species names, missing GBIF credentials, or a non-list `predicate_extra`.
#'
#' Failure to resolve or submit an individual species is handled fail-soft: the
#' species receives an explicit row in the attached `submission_audit`, while
#' processing can continue for other species. This prevents a later failed
#' submission from discarding download keys that were successfully obtained
#' earlier in the call.
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
#' @param max_active_downloads Integer. Conservative maximum number of PREPARING/RUNNING GBIF downloads allowed for the authenticated account before a new request is submitted. Defaults to 1.
#' @param slot_poll_seconds Numeric. Seconds between account-level download-slot checks and retry sleeps.
#' @param slot_max_polls Integer. Maximum number of account-level slot checks before submission is attempted and server-side retry handling takes over.
#' @param submission_max_tries Integer. Maximum number of retries when GBIF explicitly reports that download capacity is full.
#' @param ... Additional arguments reserved for compatibility with pipeline code
#'   and future backend extensions. Currently ignored.
#'
#' @return A named character vector containing GBIF occurrence-download keys for
#'   successfully submitted species. Names are the corresponding species names.
#'   The vector can be empty or contain fewer elements than the number of species
#'   requested when one or more species cannot be resolved or submitted.
#'
#'   The returned vector carries a `submission_audit` attribute with one row per
#'   requested species and the columns `label`, `gbif_key`,
#'   `submission_status`, `message`, and `attempts`. Successful rows use
#'   `submission_status = "submitted"`; failed rows retain the reason and attempt
#'   count so that missing species cannot be silently confused with zero GBIF
#'   occurrences.
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
                                max_active_downloads = 1L,
                                slot_poll_seconds = 30,
                                slot_max_polls = 120L,
                                submission_max_tries = 120L,
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

  out <- character(0)
  audit <- data.frame(
    label = character(),
    gbif_key = character(),
    submission_status = character(),
    message = character(),
    attempts = integer(),
    stringsAsFactors = FALSE
  )

  for (sp in species_vec) {
    if (!isTRUE(quiet)) message("Submitting GBIF marine download: ", sp)

    backbone <- tryCatch(
      rgbif::name_backbone(name = sp),
      error = function(e) NULL
    )

    taxon_key <- NULL
    if (!is.null(backbone) && "usageKey" %in% names(backbone)) {
      taxon_key <- backbone$usageKey
    }

    if (is.null(taxon_key) || is.na(taxon_key) || !nzchar(as.character(taxon_key))) {
      audit <- rbind(
        audit,
        data.frame(
          label = sp,
          gbif_key = NA_character_,
          submission_status = "taxon_key_failed",
          message = paste0("Could not resolve a GBIF usageKey for species: ", sp),
          attempts = 0L,
          stringsAsFactors = FALSE
        )
      )
      next
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

    submission <- .bf_submit_gbif_with_retry(
      predicate_args = predicates,
      user = user,
      pwd = pwd,
      email = email,
      label = paste0("Marine GBIF request for ", sp),
      max_active_downloads = max_active_downloads,
      slot_poll_seconds = slot_poll_seconds,
      slot_max_polls = slot_max_polls,
      submission_max_tries = submission_max_tries,
      quiet = quiet
    )

    audit <- rbind(
      audit,
      data.frame(
        label = sp,
        gbif_key = submission$key,
        submission_status = submission$status,
        message = submission$message,
        attempts = as.integer(submission$attempts),
        stringsAsFactors = FALSE
      )
    )

    if (identical(submission$status, "submitted")) {
      out[[sp]] <- submission$key
      if (!isTRUE(quiet)) {
        message("GBIF marine download key for ", sp, ": ", submission$key)
      }
    } else if (!isTRUE(quiet)) {
      message(
        "GBIF marine submission failed for ", sp,
        " [", submission$status, "]: ", submission$message
      )
    }
  }

  attr(out, "submission_audit") <- audit
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

#' Submit one multi-country GBIF download for a terrestrial/freshwater species
#'
#' Submit a single asynchronous GBIF occurrence download for one species across
#' one or more requested ISO2 countries. This backend is used by the
#' terrestrial/freshwater biofetchR pipeline, where the logical input and audit
#' units are species x recipient-country combinations but the GBIF download unit
#' can contain several countries for the same species.
#'
#' @description
#' `download_gbif_batch_gadm()` provides the queue-aware GBIF submission backend
#' used by `process_gbif_terrestrial_freshwater_pipeline()`.
#'
#' Rather than creating one asynchronous GBIF download job for every
#' species-country combination, all valid requested countries for the supplied
#' species are combined into a single GBIF request. The resulting archive is
#' subsequently imported once and split back to country-level occurrence objects
#' by [wait_and_import_gbif()].
#'
#' This design reduces the number of unfinished asynchronous GBIF jobs created by
#' multi-country workflows while preserving one explicit country-level request
#' label, retrieval outcome and downstream output for every requested country.
#'
#' @details
#' The function performs the following steps:
#'
#' 1. Input standardisation: `species` is trimmed to a single character value.
#'    `iso2_codes` are converted to upper case, trimmed, deduplicated and stripped
#'    of missing or empty values.
#'
#' 2. GBIF taxon resolution: The species is resolved to a GBIF taxon key.
#'    When the package-visible [get_taxon_key()] helper is available it is used;
#'    otherwise the function falls back to [rgbif::name_backbone()].
#'
#' 3. Multi-country predicate construction: All requested ISO2 countries are
#'    combined into one [rgbif::pred_in()] country predicate. The submitted
#'    occurrence-download request therefore contains:
#'
#'    \itemize{
#'      \item the resolved GBIF `taxonKey`;
#'      \item `hasCoordinate = TRUE`;
#'      \item `country IN <all requested ISO2 country codes>`.
#'    }
#'
#' 4. Queue-aware submission: The request is submitted through biofetchR's
#'    queue-aware retry helper. Before each submission attempt, the authenticated
#'    GBIF account is checked for currently `PREPARING` or `RUNNING` occurrence
#'    downloads.
#'
#' 5. Capacity-limit recovery: If GBIF explicitly reports that the account
#'    has reached its current incomplete/concurrent-download capacity, the SAME
#'    multi-country request is retained, paused and retried. Countries are not
#'    dropped merely because a submission attempt encountered a temporary GBIF
#'    capacity limit.
#'
#' 6. Shared-key reconstruction metadata: After a successful submission,
#'    every requested ISO2 country is mapped to the same GBIF download key. The
#'    returned object carries metadata telling [wait_and_import_gbif()] that this
#'    is a shared multi-country download and that the imported archive should be
#'    split using GBIF's `countryCode` column.
#'
#' @section Species-country requests versus GBIF download jobs:
#' A central distinction in this backend is that the biological/request unit and
#' the asynchronous GBIF job unit are not necessarily the same.
#'
#' For example, an input consisting of:
#'
#' \itemize{
#'   \item species A x country GB;
#'   \item species A x country IE;
#'   \item species A x country FR;
#'   \item species A x country ES
#' }
#'
#' remains four species-country retrieval requests for auditing and downstream
#' processing, but is submitted to GBIF as ONE asynchronous occurrence download
#' for species A with:
#'
#' `country IN c("GB", "IE", "FR", "ES")`.
#'
#' The completed GBIF archive is then imported once and split locally back into
#' the four requested country subsets.
#'
#' @section Queue and concurrency protection:
#' GBIF occurrence downloads are asynchronous and the number of incomplete jobs
#' allowed for one account can be limited by the service. The backend therefore
#' applies conservative queue protection before submission.
#'
#' `max_active_downloads` controls the maximum number of account-level
#' `PREPARING`/`RUNNING` downloads permitted before another request is attempted.
#' The default is `1L`, meaning that biofetchR waits while another active download
#' is visible before submitting a new one.
#'
#' Account-level queue checks are repeated every `slot_poll_seconds` seconds for
#' at most `slot_max_polls` checks. If the account-level download list cannot be
#' inspected reliably, submission is allowed to proceed and GBIF's server-side
#' response becomes the final capacity guard.
#'
#' If GBIF then returns a recognised concurrency or incomplete-download limit,
#' the same request is retried up to `submission_max_tries` times. A temporary
#' capacity error therefore does not consume a country request or silently remove
#' that country from the workflow.
#'
#' @section Shared GBIF download keys:
#' On successful submission, the returned object contains one named list element
#' for every requested ISO2 country. All elements intentionally contain the SAME
#' GBIF download key because all countries were included in one asynchronous
#' download.
#'
#' For example, the returned object may conceptually contain:
#'
#' \preformatted{
#' GB = "0001234-260101120000000"
#' IE = "0001234-260101120000000"
#' FR = "0001234-260101120000000"
#' ES = "0001234-260101120000000"
#' }
#'
#' This repeated key is expected behaviour and preserves the pre-0.1.1
#' country-labelled list interface while avoiding multiple redundant GBIF jobs.
#'
#' [wait_and_import_gbif()] deduplicates the shared key, polls and imports the
#' archive once, then reconstructs separate country-labelled results using
#' `countryCode`.
#'
#' @section Returned metadata attributes:
#' The returned object can carry the following attributes:
#'
#' \itemize{
#'   \item `submission_audit`: submission-level audit information;
#'   \item `biofetchR_shared_download = TRUE`: identifies the object as a shared
#'     multi-country GBIF download;
#'   \item `biofetchR_split_field = "countryCode"`: tells the importer which GBIF
#'     field should be used to reconstruct country subsets;
#'   \item `biofetchR_requested_labels`: the original requested ISO2 country
#'     labels.
#' }
#'
#' These attributes are internal workflow metadata and should normally be passed
#' unchanged to [wait_and_import_gbif()] or [wait_and_import_gbif_safe()].
#'
#' @section Submission audit:
#' The `submission_audit` attribute contains one row per requested ISO2 country,
#' even though the countries share one GBIF job. Its fields are:
#'
#' \itemize{
#'   \item `label`: requested ISO2 country code;
#'   \item `gbif_key`: shared GBIF download key when available;
#'   \item `submission_status`: submission outcome;
#'   \item `message`: diagnostic message when relevant;
#'   \item `attempts`: number of GBIF submission attempts used.
#' }
#'
#' Successful countries normally share `submission_status = "submitted"` and the
#' same `gbif_key`. Failed taxon resolution or download submission is represented
#' explicitly in the audit rather than being silently interpreted as a
#' zero-occurrence result.
#'
#' @section Failure handling:
#' Empty or invalid species/country input is rejected before submission.
#'
#' If no GBIF taxon key can be resolved for an otherwise valid species request,
#' the function returns a length-zero object carrying a submission audit with
#' `submission_status = "taxon_key_failed"`.
#'
#' Submission-level failures are returned through statuses produced by the
#' queue-aware submission helper, including:
#'
#' \itemize{
#'   \item `submit_no_key`: GBIF returned a response but no usable download key
#'     could be extracted;
#'   \item `submit_failed`: submission failed for a non-capacity reason;
#'   \item `submit_timeout`: the configured capacity-retry limit was exhausted.
#' }
#'
#' These states remain distinct from downstream polling/import failures handled
#' by [wait_and_import_gbif()]. In particular, a failed submission or unresolved
#' asynchronous download must not be interpreted as evidence that GBIF contained
#' zero occurrence records for the country.
#'
#' @section GBIF data use and citation:
#' This function submits live GBIF occurrence downloads and returns GBIF download
#' keys. Users should retain those keys and cite the corresponding GBIF download
#' DOI(s) for occurrence records used in research, reports or policy outputs.
#'
#' A single multi-country submission normally corresponds to one GBIF download
#' DOI shared by all country subsets reconstructed from that archive. The
#' country-specific outputs should therefore retain the common GBIF key so that
#' all derived records remain traceable to the exact source download.
#'
#' GBIF aggregates occurrence records from many contributing datasets.
#' Individual datasets may carry their own licences and attribution requirements.
#' biofetchR does not redistribute GBIF occurrence data and does not replace
#' users' responsibility to comply with GBIF and publisher data-use terms.
#'
#' @section Side effects:
#' Submits a live asynchronous GBIF occurrence-download request through
#' [rgbif::occ_download()] when a usable taxon key and at least one ISO2 country
#' are available.
#'
#' The function can also query the authenticated user's GBIF download list through
#' [rgbif::occ_download_list()] as part of proactive queue management, and prints
#' submission, waiting and retry messages through `cli`.
#'
#' The function itself does not wait for the submitted occurrence download to
#' finish, download the SIMPLE_CSV archive, clean coordinates, thin occurrences,
#' or assign records to spatial units. Those operations occur downstream.
#'
#' @section Notes:
#' The backend preserves the historical public return shape of
#' `download_gbif_batch_gadm()`: a named list whose names represent requested
#' countries.
#'
#' The important behavioural change is that those country labels can now point to
#' one shared GBIF key. Code should therefore not assume that the number of
#' returned list elements equals the number of distinct GBIF download jobs.
#'
#' @param species Character scalar. Scientific name of the species to resolve and
#'   submit to GBIF. Only the first supplied value is used after trimming.
#' @param iso2_codes Character vector of requested ISO2 country codes. Codes are
#'   trimmed, converted to upper case, deduplicated, and stripped of missing or
#'   empty values before submission.
#' @param user GBIF username used for authenticated occurrence-download
#'   submission and account-level queue inspection.
#' @param pwd GBIF password used for authenticated occurrence-download submission
#'   and account-level queue inspection.
#' @param email GBIF account email passed to [rgbif::occ_download()].
#' @param max_active_downloads Integer. Conservative maximum number of
#'   `PREPARING`/`RUNNING` GBIF occurrence downloads permitted for the
#'   authenticated account before another request is attempted. Defaults to
#'   `1L`.
#' @param slot_poll_seconds Numeric. Seconds between proactive account-level
#'   download-slot checks and between retries following recognised GBIF capacity
#'   errors.
#' @param slot_max_polls Integer. Maximum number of proactive account-level slot
#'   checks performed before submission is attempted and server-side retry
#'   handling becomes the final capacity guard.
#' @param submission_max_tries Integer. Maximum number of attempts when GBIF
#'   repeatedly reports that asynchronous download capacity is full.
#'
#' @return On successful submission, a named list with one element for every
#'   requested ISO2 country. Each element contains the SAME shared GBIF
#'   occurrence-download key because all countries were submitted in one
#'   multi-country request.
#'
#'   The returned object carries a `submission_audit` attribute and, for the
#'   multi-country workflow, the attributes `biofetchR_shared_download`,
#'   `biofetchR_split_field`, and `biofetchR_requested_labels`. These attributes
#'   allow [wait_and_import_gbif()] to poll/import the shared archive once and
#'   reconstruct country-level results.
#'
#'   If taxon resolution or GBIF submission fails for otherwise valid input, a
#'   length-zero object can be returned with submission-audit information
#'   describing the failure. Invalid or empty species/country input can return
#'   `NULL` after a warning.
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
#'
#'   keys
#'   attr(keys, "submission_audit")
#' }
#' }
#'
#' @export
#' @family GBIF download backends
#' @seealso [download_gbif_batch()], [wait_and_import_gbif()],
#'   [wait_and_import_gbif_safe()]
#' @importFrom cli cli_alert_warning cli_alert_success
#' @importFrom rgbif name_backbone pred pred_and pred_in occ_download
#'   occ_download_list
#' @md
download_gbif_batch_gadm <- function(species,
                                     iso2_codes,
                                     user,
                                     pwd,
                                     email,
                                     max_active_downloads = 1L,
                                     slot_poll_seconds = 30,
                                     slot_max_polls = 120L,
                                     submission_max_tries = 120L) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Package 'rgbif' is required for GBIF downloads.", call. = FALSE)
  }

  if (!requireNamespace("cli", quietly = TRUE)) {
    stop("Package 'cli' is required for GBIF download messages.", call. = FALSE)
  }

  species <- trimws(as.character(species)[1])
  iso2_codes <- unique(toupper(trimws(as.character(iso2_codes))))
  iso2_codes <- iso2_codes[!is.na(iso2_codes) & nzchar(iso2_codes)]

  if (is.na(species) || !nzchar(species)) {
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
    keys <- character(0)
    attr(keys, "submission_audit") <- data.frame(
      label = iso2_codes,
      gbif_key = NA_character_,
      submission_status = "taxon_key_failed",
      message = paste0("Could not resolve GBIF taxonKey for ", species, "."),
      attempts = 0L,
      stringsAsFactors = FALSE
    )
    return(keys)
  }

  # One species + all requested countries = one asynchronous GBIF download.
  # The imported SIMPLE_CSV archive is split back to country-level objects in
  # wait_and_import_gbif() using `countryCode`.
  country_predicate <- rgbif::pred_in("country", iso2_codes)

  submission <- .bf_submit_gbif_with_retry(
    predicate_args = list(
      rgbif::pred_and(
        rgbif::pred("taxonKey", taxon_key),
        rgbif::pred("hasCoordinate", TRUE),
        country_predicate
      )
    ),
    user = user,
    pwd = pwd,
    email = email,
    label = paste0(
      "Terrestrial/freshwater GBIF request for ", species,
      " across ", length(iso2_codes), " country/countries"
    ),
    max_active_downloads = max_active_downloads,
    slot_poll_seconds = slot_poll_seconds,
    slot_max_polls = slot_max_polls,
    submission_max_tries = submission_max_tries,
    quiet = FALSE
  )

  audit <- data.frame(
    label = iso2_codes,
    gbif_key = rep(submission$key, length(iso2_codes)),
    submission_status = rep(submission$status, length(iso2_codes)),
    message = rep(submission$message, length(iso2_codes)),
    attempts = rep(as.integer(submission$attempts), length(iso2_codes)),
    stringsAsFactors = FALSE
  )

  if (!identical(submission$status, "submitted")) {
    cli::cli_alert_warning(
      "GBIF download submission failed for {.val {species}} across {.val {length(iso2_codes)}} country/countries: {submission$message}"
    )
    keys <- character(0)
    attr(keys, "submission_audit") <- audit
    attr(keys, "biofetchR_shared_download") <- TRUE
    attr(keys, "biofetchR_split_field") <- "countryCode"
    attr(keys, "biofetchR_requested_labels") <- iso2_codes
    return(keys)
  }

  # Preserve the pre-0.1.1 public return type: one named list element per
  # requested country. The values intentionally share the same GBIF key because
  # the countries were submitted together in one asynchronous download.
  keys <- stats::setNames(
    rep(list(submission$key), length(iso2_codes)),
    iso2_codes
  )

  attr(keys, "submission_audit") <- audit
  attr(keys, "biofetchR_shared_download") <- TRUE
  attr(keys, "biofetchR_split_field") <- "countryCode"
  attr(keys, "biofetchR_requested_labels") <- iso2_codes

  cli::cli_alert_success(
    "Submitted one terrestrial/freshwater GBIF download for {.val {species}} across {.val {length(iso2_codes)}} country/countries | key: {.val {submission$key}}"
  )

  keys
}

