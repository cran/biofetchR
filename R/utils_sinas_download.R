################################################################################
# utils_sinas_download.R
# -----------------------------------------------------------------------------
# biofetchR: package-managed SInAS resource downloader
# -----------------------------------------------------------------------------
#
# Purpose
#   This script provides a small, package-safe downloader for the SInAS 3.1.1
#   resources required by the biofetchR native-evidence workflow. It resolves and
#   caches three required input resources:
#
#     * SInAS_3.1.1.csv
#     * AllLocations.xlsx, recovered from All_Config_Files_SInAS_v3.1.1.zip
#     * SInAS_3.1.1_FullTaxaList.csv, recovered from
#       All_Output_Files_SInAS_v3.1.1.zip for the default Zenodo record
#
# Design notes
#   * The public API is intentionally small:
#       - bf_sinas_default_urls()
#       - bf_download_sinas_resources()
#   * Internal helpers are documented with roxygen2 and @noRd so the source is
#     maintainable without generating extra user-facing help pages.
#   * Downloads are cached locally and can be refreshed with force_refresh = TRUE.
#   * User-supplied local paths take precedence over downloads where supplied.
#
# Important implementation detail
#   The default SInAS 3.1.1 Zenodo record does not expose
#   SInAS_3.1.1_FullTaxaList.csv as a standalone file. The downloader therefore
#   recovers that file from All_Output_Files_SInAS_v3.1.1.zip by default. A
#   standalone FullTaxaList URL is only attempted when the user explicitly passes
#   fulltaxa_url.
#
# Data use and attribution
#   biofetchR does not redistribute SInAS data. This helper only downloads files
#   into a user-managed local cache. Users should check the relevant SInAS/Zenodo
#   record for licence, citation and attribution requirements before publishing,
#   redistributing, or archiving derived outputs.
#
################################################################################


#' Return default SInAS 3.1.1 resource URLs
#'
#' Constructs the direct Zenodo file URLs used by
#' [bf_download_sinas_resources()]. The default record is the SInAS 3.1.1 record
#' used by the biofetchR native-evidence workflow.
#'
#' @param record_id Character scalar. Zenodo record identifier used to construct
#'   the file URLs. Defaults to `"18220953"`.
#'
#' @return A named character vector of direct Zenodo file URLs used by the SInAS
#'   downloader. The vector contains four entries: `main_csv`, the URL for
#'   `SInAS_3.1.1.csv`; `config_zip`, the URL for
#'   `All_Config_Files_SInAS_v3.1.1.zip`, which contains `AllLocations.xlsx`;
#'   `output_zip`, the URL for `All_Output_Files_SInAS_v3.1.1.zip`, which
#'   contains `SInAS_3.1.1_FullTaxaList.csv` for the default record; and
#'   `fulltaxa_csv`, a legacy-style standalone FullTaxaList URL returned for
#'   completeness. The function only constructs URLs and does not download,
#'   cache, read or validate any SInAS files.
#'
#' @section Data access:
#' This function only constructs URLs. It does not download, cache, read or
#' redistribute SInAS data.
#'
#' @family SInAS resource helpers
#'
#' @examples
#' bf_sinas_default_urls()
#'
#' @export
bf_sinas_default_urls <- function(record_id = "18220953") {
  base <- paste0("https://zenodo.org/records/", record_id, "/files/")

  c(
    main_csv = paste0(base, "SInAS_3.1.1.csv?download=1"),
    config_zip = paste0(base, "All_Config_Files_SInAS_v3.1.1.zip?download=1"),
    output_zip = paste0(base, "All_Output_Files_SInAS_v3.1.1.zip?download=1"),
    fulltaxa_csv = paste0(base, "SInAS_3.1.1_FullTaxaList.csv?download=1")
  )
}


#' Test whether a user-supplied local file should be used
#'
#' Checks whether an optional path argument is a single, non-missing,
#' non-empty path that already exists on disk. This is used so manually supplied
#' SInAS files take precedence over remote downloads.
#'
#' @param path Candidate file path.
#'
#' @return Logical scalar. `TRUE` when `path` is a usable existing file;
#'   otherwise `FALSE`.
#'
#' @keywords internal
#' @noRd
.bf_sinas_use_existing <- function(path) {
  !is.null(path) &&
    length(path) == 1L &&
    !is.na(path) &&
    nzchar(path) &&
    file.exists(path)
}


#' Test whether an optional URL has been supplied
#'
#' Checks whether a URL argument is a single, non-missing, non-empty string. This
#' is used for optional resources such as a legacy standalone FullTaxaList CSV.
#'
#' @param url Candidate URL.
#'
#' @return Logical scalar. `TRUE` when a usable URL-like string was supplied;
#'   otherwise `FALSE`.
#'
#' @keywords internal
#' @noRd
.bf_sinas_has_url <- function(url) {
  !is.null(url) &&
    length(url) == 1L &&
    !is.na(url) &&
    nzchar(url)
}

#' Find a required SInAS file inside an extracted archive
#'
#' Recursively searches an extracted SInAS directory for a file matching a regular
#' expression. If multiple matches are found, CSV files are preferred because the
#' workflow expects tabular SInAS resources.
#'
#' @param root Character scalar. Directory to search recursively.
#' @param pattern Character scalar. Regular expression passed to [list.files()].
#' @param required Logical scalar. If `TRUE`, throw an error when no file is
#'   found. If `FALSE`, return `NULL` when no file is found.
#'
#' @return Normalised path to the first matching file, or `NULL` when no match is
#'   found and `required = FALSE`.
#'
#' @keywords internal
#' @noRd
.bf_sinas_find_file <- function(root, pattern, required = TRUE) {
  hits <- list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  hits <- hits[file.exists(hits)]

  if (!length(hits)) {
    if (isTRUE(required)) {
      stop(
        "Could not find required SInAS file matching pattern `",
        pattern,
        "` under: ",
        root,
        call. = FALSE
      )
    }

    return(NULL)
  }

  # Prefer CSV if multiple text-like matches exist.
  csv_hits <- hits[grepl("\\.csv$", hits, ignore.case = TRUE)]
  if (length(csv_hits)) hits <- csv_hits

  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}


#' Download and locate required SInAS resources
#'
#' Downloads, caches and resolves the SInAS 3.1.1 resources required by the
#' biofetchR native-evidence workflow. Existing local files can be supplied to
#' bypass downloading. Otherwise, the function downloads the main SInAS table,
#' extracts `AllLocations.xlsx` from the configuration archive, and recovers
#' `SInAS_3.1.1_FullTaxaList.csv` from the output archive for the default record.
#'
#' @param cache_dir Character scalar. Local cache directory for downloaded and
#'   extracted SInAS resources. Must be supplied explicitly when any required
#'   SInAS resource needs to be downloaded or extracted. If `main_path`,
#'   `allloc_path` and `fulltaxa_path` all point to existing local files,
#'   `cache_dir` is not required. In examples, tests and vignettes, use a path
#'   under `tempdir()`.
#' @param force_refresh Logical scalar. If `TRUE`, re-download and re-extract
#'   cached resources even when existing files are present.
#' @param quiet Logical scalar. If `TRUE`, suppress progress messages.
#' @param main_path Optional character scalar. Existing local path to
#'   `SInAS_3.1.1.csv`. When supplied and valid, this path is used instead of
#'   downloading `main_url`.
#' @param allloc_path Optional character scalar. Existing local path to
#'   `AllLocations.xlsx`. When supplied and valid, this path is used instead of
#'   searching the extracted configuration archive.
#' @param fulltaxa_path Optional character scalar. Existing local path to
#'   `SInAS_3.1.1_FullTaxaList.csv`. When supplied and valid, this path is used
#'   instead of downloading or extracting FullTaxaList.
#' @param record_id Character scalar. Zenodo record identifier used by
#'   [bf_sinas_default_urls()] to construct default URLs.
#' @param main_url Character scalar. Direct URL for the main SInAS CSV.
#' @param fulltaxa_url Optional character scalar. Direct URL for a standalone
#'   FullTaxaList CSV. Defaults to `NULL` because the default SInAS 3.1.1 record
#'   stores FullTaxaList inside the output ZIP rather than as a standalone file.
#' @param config_zip_url Character scalar. Direct URL for the configuration ZIP
#'   containing `AllLocations.xlsx`.
#' @param output_zip_url Character scalar. Direct URL for the output ZIP
#'   containing `SInAS_3.1.1_FullTaxaList.csv` for the default record.
#'
#' @return A named list of normalised local file paths to the SInAS resources
#'   required by biofetchR native-origin evidence workflows. The list always
#'   contains `main_csv`, the local path to `SInAS_3.1.1.csv`;
#'   `alllocations_xlsx`, the local path to `AllLocations.xlsx` or an equivalent
#'   location table found in the extracted configuration archive; `fulltaxa_csv`,
#'   the local path to `SInAS_3.1.1_FullTaxaList.csv` or a matching FullTaxaList
#'   file recovered from the output archive; `config_zip`, the cached
#'   configuration archive path; and `config_dir`, the extracted configuration
#'   directory path. When FullTaxaList is recovered from
#'   `All_Output_Files_SInAS_v3.1.1.zip`, the list also contains `output_zip` and
#'   `output_dir`. User-supplied local paths are returned in normalised form when
#'   provided and valid; otherwise paths point to files downloaded or extracted
#'   inside `cache_dir`.
#'
#' @section Cache behaviour:
#' Cached files are reused when present and above the minimum size threshold.
#' Set `force_refresh = TRUE` to re-download archives and rebuild extracted
#' directories. User-supplied local paths always take precedence over downloads.
#'
#' @section Data access and attribution:
#' biofetchR does not ship or redistribute SInAS resources. This function only
#' downloads files into a user-managed local cache. Users should check the
#' relevant SInAS/Zenodo record for licence, citation and attribution
#' requirements before publishing, redistributing or archiving derived outputs.
#'
#' @family SInAS resource helpers
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   sinas_paths <- bf_download_sinas_resources(
#'     cache_dir = file.path(tempdir(), "biofetchR_sinas"),
#'     quiet = FALSE
#'   )
#'
#'   names(sinas_paths)
#'   sinas_paths$main_csv
#' }
#' }
#'
#' @export
bf_download_sinas_resources <- function(cache_dir = NULL,
                                        force_refresh = FALSE,
                                        quiet = FALSE,
                                        main_path = NULL,
                                        allloc_path = NULL,
                                        fulltaxa_path = NULL,
                                        record_id = "18220953",
                                        main_url = bf_sinas_default_urls(record_id)[["main_csv"]],
                                        fulltaxa_url = NULL,
                                        config_zip_url = bf_sinas_default_urls(record_id)[["config_zip"]],
                                        output_zip_url = bf_sinas_default_urls(record_id)[["output_zip"]]) {
  needs_cache <- !.bf_sinas_use_existing(main_path) ||
    !.bf_sinas_use_existing(allloc_path) ||
    !.bf_sinas_use_existing(fulltaxa_path)

  if (isTRUE(needs_cache)) {
    if (is.null(cache_dir) || length(cache_dir) == 0L ||
        !nzchar(trimws(as.character(cache_dir[[1L]])))) {
      stop(
        "`cache_dir` must be supplied explicitly when SInAS resources need to be downloaded or extracted. ",
        "Alternatively, supply existing local paths for `main_path`, `allloc_path` and `fulltaxa_path`.",
        call. = FALSE
      )
    }

    cache_dir <- normalizePath(
      as.character(cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  main_csv <- if (.bf_sinas_use_existing(main_path)) {
    normalizePath(main_path, winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(
      bf_download_cached(
        url = main_url,
        dest = file.path(cache_dir, "SInAS_3.1.1.csv"),
        force_refresh = force_refresh,
        quiet = quiet,
        min_bytes = 1000L,
        validate_not_html = TRUE
      ),
      winslash = "/",
      mustWork = TRUE
    )
  }

  config_zip <- NA_character_
  config_dir <- NA_character_

  alllocations_xlsx <- if (.bf_sinas_use_existing(allloc_path)) {
    normalizePath(allloc_path, winslash = "/", mustWork = TRUE)
  } else {
    config_zip <- normalizePath(
      bf_download_cached(
        url = config_zip_url,
        dest = file.path(cache_dir, "All_Config_Files_SInAS_v3.1.1.zip"),
        force_refresh = force_refresh,
        quiet = quiet,
        min_bytes = 1000L,
        validate_not_html = TRUE
      ),
      winslash = "/",
      mustWork = TRUE
    )

    config_dir <- normalizePath(
      bf_unzip_cached(
        zipfile = config_zip,
        exdir = file.path(cache_dir, "All_Config_Files_SInAS_v3.1.1"),
        force_refresh = force_refresh,
        quiet = quiet
      ),
      winslash = "/",
      mustWork = TRUE
    )

    .bf_sinas_find_file(config_dir, "^AllLocations\\.(xlsx|xls|csv|tsv)$", required = TRUE)
  }

  output_zip <- NA_character_
  output_dir <- NA_character_

  fulltaxa_csv <- if (.bf_sinas_use_existing(fulltaxa_path)) {
    normalizePath(fulltaxa_path, winslash = "/", mustWork = TRUE)
  } else {
    standalone <- NULL

    # Only attempt the legacy direct FullTaxaList URL if the user explicitly
    # supplies one. The default SInAS 3.1.1 Zenodo record does not have this file
    # as a standalone download.
    if (.bf_sinas_has_url(fulltaxa_url)) {
      standalone <- tryCatch(
        normalizePath(
          bf_download_cached(
            url = fulltaxa_url,
            dest = file.path(cache_dir, "SInAS_3.1.1_FullTaxaList.csv"),
            force_refresh = force_refresh,
            quiet = quiet,
            min_bytes = 1000L,
            validate_not_html = TRUE
          ),
          winslash = "/",
          mustWork = TRUE
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(standalone) && file.exists(standalone)) {
      standalone
    } else {
      if (!.bf_sinas_has_url(output_zip_url)) {
        stop(
          "FullTaxaList was not supplied and no `output_zip_url` was available. ",
          "For SInAS record 18220953, FullTaxaList must be recovered from ",
          "All_Output_Files_SInAS_v3.1.1.zip.",
          call. = FALSE
        )
      }

      output_zip <- normalizePath(
        bf_download_cached(
          url = output_zip_url,
          dest = file.path(cache_dir, "All_Output_Files_SInAS_v3.1.1.zip"),
          force_refresh = force_refresh,
          quiet = quiet,
          min_bytes = 1000L,
          validate_not_html = TRUE
        ),
        winslash = "/",
        mustWork = TRUE
      )

      output_dir <- normalizePath(
        bf_unzip_cached(
          zipfile = output_zip,
          exdir = file.path(cache_dir, "All_Output_Files_SInAS_v3.1.1"),
          force_refresh = force_refresh,
          quiet = quiet
        ),
        winslash = "/",
        mustWork = TRUE
      )

      .bf_sinas_find_file(
        output_dir,
        "(FullTaxaList|Full_Taxa_List).*\\.(csv|tsv|txt)$|SInAS.*FullTaxaList.*\\.(csv|tsv|txt)$",
        required = TRUE
      )
    }
  }

  out <- list(
    main_csv = main_csv,
    alllocations_xlsx = alllocations_xlsx,
    fulltaxa_csv = fulltaxa_csv
  )

  if (!is.na(config_zip) && nzchar(config_zip)) out$config_zip <- config_zip
  if (!is.na(config_dir) && nzchar(config_dir)) out$config_dir <- config_dir
  if (!is.na(output_zip) && nzchar(output_zip)) out$output_zip <- output_zip
  if (!is.na(output_dir) && nzchar(output_dir)) out$output_dir <- output_dir

  out
}
