################################################################################
# utils_sinas_download.R
# -----------------------------------------------------------------------------
# biofetchR: package-managed SInAS resource downloader
# -----------------------------------------------------------------------------
#
# Default release
#   SInAS dataset 3.2
#   Zenodo record 21933976
#   DOI 10.5281/zenodo.21933976
#   Generated with SInAS workflow v2.0
#
# The current SInAS 3.2 record provides the main database as SInAS_3.2.csv and
# packages configuration and workflow outputs in All_Config_Files_v3.2.zip and
# All_Output_Files_v3.2.zip. The FullTaxaList used for alias reconciliation is
# recovered from the output archive.
#
# SInAS 3.1.1 / Zenodo 18220953 remains supported when explicitly requested.
################################################################################


#' Resolve metadata for a supported SInAS release
#'
#' @param record_id Zenodo record identifier.
#'
#' @return Named list describing filenames and release provenance.
#'
#' @keywords internal
#' @noRd
.bf_sinas_release <- function(record_id = "21933976") {
  record_id <- trimws(as.character(record_id)[1L])

  if (identical(record_id, "21933976")) {
    return(list(
      record_id = "21933976",
      dataset_version = "3.2",
      workflow_version = "2.0",
      doi = "10.5281/zenodo.21933976",
      main_file = "SInAS_3.2.csv",
      config_file = "All_Config_Files_v3.2.zip",
      output_file = "All_Output_Files_v3.2.zip",
      fulltaxa_file = "SInAS_3.2_FullTaxaList.csv"
    ))
  }

  if (identical(record_id, "18220953")) {
    return(list(
      record_id = "18220953",
      dataset_version = "3.1.1",
      workflow_version = "2.0",
      doi = "10.5281/zenodo.18220953",
      main_file = "SInAS_3.1.1.csv",
      config_file = "All_Config_Files_SInAS_v3.1.1.zip",
      output_file = "All_Output_Files_SInAS_v3.1.1.zip",
      fulltaxa_file = "SInAS_3.1.1_FullTaxaList.csv"
    ))
  }

  stop(
    "Unsupported package-managed SInAS Zenodo record: ", record_id, ". ",
    "biofetchR currently knows the resource layout for SInAS 3.2 ",
    "(record 21933976) and SInAS 3.1.1 (record 18220953). ",
    "For another release, supply explicit local resource paths or update the ",
    "SInAS release manifest.",
    call. = FALSE
  )
}


#' Return default SInAS resource URLs
#'
#' Constructs direct Zenodo file URLs for a supported SInAS release. The default
#' is SInAS dataset 3.2 (Zenodo record 21933976), generated with SInAS workflow
#' v2.0.
#'
#' @param record_id Character scalar. Supported Zenodo record identifier.
#'   Defaults to `"21933976"` (SInAS 3.2). Record `"18220953"` resolves the
#'   previous SInAS 3.1.1 layout.
#'
#' @return A named character vector with `main_csv`, `config_zip`,
#'   `output_zip`, and `fulltaxa_csv`. `fulltaxa_csv` is `NA` because the
#'   package-managed releases recover FullTaxaList from `output_zip`.
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
bf_sinas_default_urls <- function(record_id = "21933976") {
  release <- .bf_sinas_release(record_id)
  base <- paste0(
    "https://zenodo.org/records/",
    release$record_id,
    "/files/"
  )

  c(
    main_csv = paste0(base, release$main_file, "?download=1"),
    config_zip = paste0(base, release$config_file, "?download=1"),
    output_zip = paste0(base, release$output_file, "?download=1"),
    fulltaxa_csv = NA_character_
  )
}


#' Test whether a user-supplied local file should be used
#'
#' @param path Candidate file path.
#'
#' @return Logical scalar.
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
#' @param url Candidate URL.
#'
#' @return Logical scalar.
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
#' @param root Directory to search recursively.
#' @param pattern Regular expression passed to [list.files()].
#' @param required Throw an error when no file is found.
#'
#' @return Normalised path or `NULL`.
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

  csv_hits <- hits[grepl("\\.csv$", hits, ignore.case = TRUE)]
  if (length(csv_hits)) hits <- csv_hits

  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}


#' Download a SInAS resource through biofetchR's cache helper
#'
#' @keywords internal
#' @noRd
.bf_sinas_download <- function(url,
                               dest,
                               force_refresh = FALSE,
                               quiet = FALSE,
                               min_bytes = 1000L) {
  if (!exists("bf_download_cached", mode = "function", inherits = TRUE)) {
    stop(
      "`bf_download_cached()` is required by the SInAS downloader but was not ",
      "available in the package namespace.",
      call. = FALSE
    )
  }

  f <- get("bf_download_cached", mode = "function", inherits = TRUE)
  args <- list(
    url = url,
    dest = dest,
    force_refresh = force_refresh,
    quiet = quiet,
    min_bytes = min_bytes,
    validate_not_html = TRUE
  )
  keep <- names(args) %in% names(formals(f)) | "..." %in% names(formals(f))
  ans <- do.call(f, args[keep])

  normalizePath(ans, winslash = "/", mustWork = TRUE)
}


#' Extract a cached SInAS archive
#'
#' This local extractor deliberately avoids depending on one of several historical
#' `bf_unzip_cached()` implementations that may coexist in older biofetchR trees.
#'
#' @keywords internal
#' @noRd
.bf_sinas_unzip <- function(zipfile,
                            exdir,
                            force_refresh = FALSE,
                            quiet = FALSE) {
  if (!file.exists(zipfile)) {
    stop("SInAS archive does not exist: ", zipfile, call. = FALSE)
  }

  marker <- file.path(exdir, ".biofetchR_sinas_unzip_complete")

  if (isTRUE(force_refresh) && dir.exists(exdir)) {
    unlink(exdir, recursive = TRUE, force = TRUE)
  }

  if (dir.exists(exdir) &&
      file.exists(marker) &&
      length(list.files(exdir, recursive = TRUE, all.files = TRUE)) > 1L) {
    if (!isTRUE(quiet)) {
      message("Using extracted SInAS directory: ", exdir)
    }
    return(normalizePath(exdir, winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(exdir)) {
    unlink(exdir, recursive = TRUE, force = TRUE)
  }
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  if (!isTRUE(quiet)) {
    message("Extracting SInAS archive: ", zipfile)
  }

  extracted <- tryCatch(
    utils::unzip(zipfile, exdir = exdir),
    error = function(e) e
  )

  if (inherits(extracted, "condition")) {
    stop(
      "Could not extract SInAS archive: ",
      zipfile,
      ". ",
      conditionMessage(extracted),
      call. = FALSE
    )
  }

  files_after <- list.files(
    exdir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  files_after <- files_after[file.exists(files_after)]

  if (!length(files_after)) {
    stop(
      "SInAS archive extracted but produced no files: ",
      zipfile,
      call. = FALSE
    )
  }

  writeLines(
    c(
      paste0("archive=", normalizePath(zipfile, winslash = "/", mustWork = TRUE)),
      paste0("extracted=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    ),
    marker,
    useBytes = TRUE
  )

  normalizePath(exdir, winslash = "/", mustWork = TRUE)
}


#' Download and locate required SInAS resources
#'
#' Downloads, caches and resolves the SInAS resources required by biofetchR's
#' native-origin evidence workflow. The default is SInAS dataset 3.2
#' (Zenodo record 21933976), generated with SInAS workflow v2.0.
#'
#' The main SInAS table is downloaded directly. `AllLocations` is recovered from
#' the configuration archive, and `FullTaxaList` is recovered from the workflow
#' output archive unless an explicit local file or standalone URL is supplied.
#'
#' @param cache_dir Local cache directory.
#' @param force_refresh Re-download and re-extract cached resources.
#' @param quiet Suppress progress messages.
#' @param main_path Optional local main SInAS table.
#' @param allloc_path Optional local `AllLocations` table.
#' @param fulltaxa_path Optional local FullTaxaList.
#' @param record_id Zenodo record identifier. Defaults to `"21933976"`.
#' @param main_url Optional direct main-table URL. By default it is derived from
#'   [bf_sinas_default_urls()].
#' @param fulltaxa_url Optional standalone FullTaxaList URL. Defaults to `NULL`;
#'   the current SInAS release stores FullTaxaList in `output_zip`.
#' @param config_zip_url Optional configuration-archive URL.
#' @param output_zip_url Optional workflow-output archive URL.
#'
#' @return Named list containing local paths `main_csv`,
#'   `alllocations_xlsx`, `fulltaxa_csv`, `config_zip`, `config_dir`,
#'   `output_zip`, and `output_dir`.
#'
#' @section Data access and attribution:
#' biofetchR does not redistribute SInAS resources. Users should cite the exact
#' SInAS release used. The default release is SInAS 3.2,
#' DOI `10.5281/zenodo.21933976`.
#'
#' @family SInAS resource helpers
#'
#' @examples
#' \donttest{
#' sinas_paths <- bf_download_sinas_resources(
#'   cache_dir = file.path(tempdir(), "biofetchR_sinas"),
#'   quiet = TRUE
#' )
#' names(sinas_paths)
#' }
#'
#' @export
bf_download_sinas_resources <- function(
    cache_dir = tools::R_user_dir("biofetchR", "data"),
    force_refresh = FALSE,
    quiet = FALSE,
    main_path = NULL,
    allloc_path = NULL,
    fulltaxa_path = NULL,
    record_id = "21933976",
    main_url = bf_sinas_default_urls(record_id)[["main_csv"]],
    fulltaxa_url = NULL,
    config_zip_url = bf_sinas_default_urls(record_id)[["config_zip"]],
    output_zip_url = bf_sinas_default_urls(record_id)[["output_zip"]]) {

  release <- .bf_sinas_release(record_id)

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  main_csv <- if (.bf_sinas_use_existing(main_path)) {
    normalizePath(main_path, winslash = "/", mustWork = TRUE)
  } else {
    .bf_sinas_download(
      url = main_url,
      dest = file.path(cache_dir, release$main_file),
      force_refresh = force_refresh,
      quiet = quiet
    )
  }

  config_zip <- .bf_sinas_download(
    url = config_zip_url,
    dest = file.path(cache_dir, release$config_file),
    force_refresh = force_refresh,
    quiet = quiet
  )

  config_dir_name <- sub("\\.zip$", "", release$config_file, ignore.case = TRUE)
  config_dir <- .bf_sinas_unzip(
    zipfile = config_zip,
    exdir = file.path(cache_dir, config_dir_name),
    force_refresh = force_refresh,
    quiet = quiet
  )

  alllocations_xlsx <- if (.bf_sinas_use_existing(allloc_path)) {
    normalizePath(allloc_path, winslash = "/", mustWork = TRUE)
  } else {
    .bf_sinas_find_file(
      config_dir,
      "^AllLocations\\.(xlsx|xls|csv|tsv)$",
      required = TRUE
    )
  }

  output_zip <- NA_character_
  output_dir <- NA_character_

  fulltaxa_csv <- if (.bf_sinas_use_existing(fulltaxa_path)) {
    normalizePath(fulltaxa_path, winslash = "/", mustWork = TRUE)
  } else {
    standalone <- NULL

    if (.bf_sinas_has_url(fulltaxa_url)) {
      standalone <- tryCatch(
        .bf_sinas_download(
          url = fulltaxa_url,
          dest = file.path(cache_dir, release$fulltaxa_file),
          force_refresh = force_refresh,
          quiet = quiet
        ),
        error = function(e) NULL
      )
    }

    if (!is.null(standalone) && file.exists(standalone)) {
      standalone
    } else {
      if (!.bf_sinas_has_url(output_zip_url)) {
        stop(
          "FullTaxaList was not supplied and `output_zip_url` is unavailable. ",
          "For SInAS ", release$dataset_version,
          ", FullTaxaList is expected inside ",
          release$output_file,
          ".",
          call. = FALSE
        )
      }

      output_zip <- .bf_sinas_download(
        url = output_zip_url,
        dest = file.path(cache_dir, release$output_file),
        force_refresh = force_refresh,
        quiet = quiet
      )

      output_dir_name <- sub(
        "\\.zip$",
        "",
        release$output_file,
        ignore.case = TRUE
      )

      output_dir <- .bf_sinas_unzip(
        zipfile = output_zip,
        exdir = file.path(cache_dir, output_dir_name),
        force_refresh = force_refresh,
        quiet = quiet
      )

      exact_hits <- list.files(
        output_dir,
        recursive = TRUE,
        full.names = TRUE
      )
      exact_hits <- exact_hits[
        basename(exact_hits) == release$fulltaxa_file &
          file.exists(exact_hits)
      ]

      exact <- if (length(exact_hits)) {
        normalizePath(exact_hits[[1L]], winslash = "/", mustWork = TRUE)
      } else {
        NULL
      }

      if (!is.null(exact)) {
        exact
      } else {
        .bf_sinas_find_file(
          output_dir,
          "(FullTaxaList|Full_Taxa_List).*\\.(csv|tsv|txt)$|SInAS.*FullTaxaList.*\\.(csv|tsv|txt)$",
          required = TRUE
        )
      }
    }
  }

  # `output_zip` / `output_dir` remain NA only when an explicitly supplied
  # standalone/local FullTaxaList made the output archive unnecessary.
  list(
    main_csv = main_csv,
    alllocations_xlsx = alllocations_xlsx,
    fulltaxa_csv = fulltaxa_csv,
    config_zip = config_zip,
    config_dir = config_dir,
    output_zip = output_zip,
    output_dir = output_dir
  )
}
