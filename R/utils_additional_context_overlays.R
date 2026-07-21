################################################################################
# utils_additional_context_overlays.R
# -----------------------------------------------------------------------------
# biofetchR: additional terrestrial/freshwater context-overlay loaders
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Provide package-managed loaders for additional terrestrial and freshwater
#   context layers used by process_gbif_terrestrial_freshwater_pipeline(). These
#   loaders convert external spatial resources into predictable WGS84 `sf`
#   objects with stable identifier/name fields so GBIF occurrence points can be
#   assigned to contextual features during the processing workflow.
#
# RETAINED OVERLAYS IN THIS FILE
#   - HydroWASTE wastewater-treatment plant points.
#   - Global Dam Watch reservoir polygons.
#   - Global-scale mining polygons.
#
# DESIGN PRINCIPLES
#   1. Prefer package-managed downloads from stable public endpoints where a
#      defensible single-file or metadata-selected source is available.
#   2. Allow `source_path` and `source_url` developer overrides so tests,
#      institutional mirrors and manually cached files can be used without
#      changing package code.
#   3. Download only the selected file needed from metadata-driven archives where
#      possible, rather than fetching entire multi-file repositories.
#   4. Fail clearly when a requested overlay is not package-managed yet, instead
#      of silently proceeding with missing context data.
#   5. Return standardised WGS84 `sf` objects with stable ID/name columns used by
#      the terrestrial/freshwater pipeline.
#
# DATA ACCESS, LICENSING AND ATTRIBUTION
#   This file can download and cache third-party spatial datasets. biofetchR does
#   not claim ownership of those datasets and does not alter their provider
#   licence terms. Users should cite and attribute the exact dataset, version and
#   provider used in their analyses. The exported loader documentation below
#   includes source-specific reminders where this is relevant.
#
# DEPENDENCIES PROVIDED ELSEWHERE IN THE PACKAGE
#   This file assumes the core overlay utilities are available, including:
#     bf_require_packages()
#     bf_download_cached()
#     bf_unpack_archive()
#     bf_read_vector_any()
#     bf_sf_wgs84()
#     bf_sf_make_valid()
#     bf_clip_to_countries()
#
################################################################################

# -----------------------------------------------------------------------------
# Small internal utilities
# -----------------------------------------------------------------------------

#' Add a standardised field from candidate source columns
#'
#' Many external overlay datasets use different attribute names for the same
#' concept, such as `NAME`, `name`, `OBJECTID`, or dataset-specific ID fields.
#' This helper searches for the first matching candidate column, copies it to a
#' stable output column, and fills a default value when no candidate is present.
#'
#' @param x Data frame or `sf` object.
#' @param out_name Name of the standardised output column to create.
#' @param candidates Character vector of possible source column names.
#' @param default Fallback value used when no candidate source column is found.
#'
#' @return `x` with `out_name` added or overwritten.
#'
#' @keywords internal
#' @noRd
.bf_add_field_from_candidates <- function(x, out_name, candidates, default = NULL) {
  nms <- names(x)
  ln  <- tolower(nms)
  pick <- NULL

  for (cand in candidates) {
    hit <- which(ln == tolower(cand))
    if (length(hit)) {
      pick <- nms[hit[1]]
      break
    }
  }

  if (!is.null(pick)) {
    x[[out_name]] <- x[[pick]]
  } else {
    x[[out_name]] <- default
  }

  x
}

#' Public alias for candidate-field standardisation
#'
#' This alias is retained for compatibility with scripts that already call the
#' helper directly after sourcing the file. It is not exported because it is an
#' implementation detail of the context-overlay loaders.
#'
#' @keywords internal
#' @noRd
# Public alias used by the loaders below. The dot-prefixed version is retained
# for compatibility with scripts that may already source/call it internally.
bf_add_field_from_candidates <- .bf_add_field_from_candidates

#' Return the first non-NULL value from a set of candidates
#'
#' Used when remote metadata records provide alternative download-link fields,
#' for example `self`, `content`, or `download` links.
#'
#' @param ... Candidate values.
#'
#' @return The first non-`NULL` value, or `NULL` if all candidates are `NULL`.
#'
#' @keywords internal
#' @noRd
.bf_first_existing <- function(...) {
  x <- list(...)
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) return(NULL)
  x[[1]]
}

#' Infer a file extension from a URL
#'
#' Query strings and fragments are removed before extracting the extension. If
#' no extension is present, a supplied default is returned.
#'
#' @param url Character URL.
#' @param default Extension to use when the URL has no visible extension.
#'
#' @return Character file extension without a leading dot.
#'
#' @keywords internal
#' @noRd
.bf_url_ext <- function(url, default = "zip") {
  u <- sub("[?#].*$", "", as.character(url)[1])
  ext <- tools::file_ext(u)
  if (!nzchar(ext)) ext <- default
  ext
}

#' Test whether a function is available on the search path
#'
#' @param nm Function name.
#'
#' @return Logical scalar.
#'
#' @keywords internal
#' @noRd
.bf_has_function <- function(nm) exists(nm, mode = "function", inherits = TRUE)

#' Stop with a clear message for unsupported package-managed overlays
#'
#' This helper is used when an overlay loader is retained as an explicit failure
#' point but no stable public download route has been configured yet.
#'
#' @param cache_stem Dataset/cache identifier used in the error message.
#' @param reason Human-readable explanation of why the overlay is unavailable.
#'
#' @return This function always errors.
#'
#' @keywords internal
#' @noRd
.bf_unsupported_package_managed <- function(cache_stem, reason) {
  stop(
    cache_stem, " is not package-managed yet. ", reason,
    " This loader is retained so the pipeline can fail explicitly when the overlay is requested, ",
    "rather than silently asking users for local/manual files.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# Metadata-driven download helpers
# ------------------------------------------------------------------------------
# These helpers download selected files from metadata APIs, instead of downloading
# an entire multi-gigabyte record. This is essential for low-memory integration
# tests and normal package use.
# ------------------------------------------------------------------------------

#' Download selected files from a Zenodo record
#'
#' Reads Zenodo record metadata, selects files whose keys match optional regular
#' expression patterns, and downloads only the selected file(s) to the cache.
#' This avoids downloading large multi-file records when a loader only needs one
#' vector archive or table.
#'
#' @param record_id Zenodo record identifier.
#' @param cache_dir Directory used for metadata and downloaded files.
#' @param force_refresh Logical; re-download metadata/files even if cached.
#' @param quiet Logical; suppress download messages where supported.
#' @param file_patterns Optional character vector of regular expressions used to
#'   select files by Zenodo file key.
#' @param max_files Maximum number of matching files to download.
#'
#' @return Character vector of downloaded file paths.
#'
#' @family additional context overlay helpers
#' @keywords internal
#' @noRd
bf_download_zenodo_record_selected <- function(record_id,
                                               cache_dir,
                                               force_refresh = FALSE,
                                               quiet = TRUE,
                                               file_patterns = NULL,
                                               max_files = 1L) {
  bf_require_packages(c("jsonlite"), context = "additional context overlay")

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  meta_path <- file.path(cache_dir, paste0("zenodo_record_", record_id, ".json"))

  bf_download_cached(
    sprintf("https://zenodo.org/api/records/%s", record_id),
    meta_path,
    force_refresh = force_refresh,
    quiet = quiet,
    mode = "wb",
    min_bytes = 100,
    validate_not_html = TRUE
  )

  meta <- jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
  files <- meta$files

  if (!length(files)) {
    stop("Zenodo record ", record_id, " returned no files.", call. = FALSE)
  }

  keys <- vapply(files, function(z) as.character(z$key %||% NA_character_), character(1))

  if (!is.null(file_patterns) && length(file_patterns)) {
    keep <- rep(FALSE, length(keys))
    for (pat in file_patterns) {
      keep <- keep | grepl(pat, keys, ignore.case = TRUE)
    }
    files <- files[keep]
    keys  <- keys[keep]
  }

  if (!length(files)) {
    stop(
      "Zenodo record ", record_id,
      " contained no files matching pattern(s): ",
      paste(file_patterns, collapse = ", "),
      call. = FALSE
    )
  }

  max_files <- suppressWarnings(as.integer(max_files))
  if (!is.finite(max_files) || max_files < 1L) max_files <- 1L
  files <- files[seq_len(min(length(files), max_files))]

  out <- character()

  for (z in files) {
    key <- as.character(z$key)
    link <- .bf_first_existing(z$links$self, z$links$content, z$links$download)

    if (is.null(link) || !nzchar(link)) {
      stop("Zenodo file has no usable download link: ", key, call. = FALSE)
    }

    dest <- file.path(cache_dir, key)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

    bf_download_cached(link, dest, force_refresh = force_refresh, quiet = quiet)
    out <- c(out, dest)
  }

  out
}

#' Download selected files from a Figshare article
#'
#' Reads Figshare article metadata, selects files whose names match optional
#' regular expression patterns, and downloads only the selected file(s) to the
#' cache. This keeps overlay loading predictable and avoids unnecessary large
#' downloads.
#'
#' @param article_id Figshare article identifier.
#' @param cache_dir Directory used for metadata and downloaded files.
#' @param force_refresh Logical; re-download metadata/files even if cached.
#' @param quiet Logical; suppress download messages where supported.
#' @param file_patterns Optional character vector of regular expressions used to
#'   select files by Figshare file name.
#' @param max_files Maximum number of matching files to download.
#'
#' @return Character vector of downloaded file paths.
#'
#' @family additional context overlay helpers
#' @keywords internal
#' @noRd
bf_download_figshare_article_selected <- function(article_id,
                                                  cache_dir,
                                                  force_refresh = FALSE,
                                                  quiet = TRUE,
                                                  file_patterns = NULL,
                                                  max_files = 1L) {
  bf_require_packages(c("jsonlite"), context = "additional context overlay")

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  meta_path <- file.path(cache_dir, paste0("figshare_article_", article_id, ".json"))

  bf_download_cached(
    sprintf("https://api.figshare.com/v2/articles/%s", article_id),
    meta_path,
    force_refresh = force_refresh,
    quiet = quiet,
    mode = "wb",
    min_bytes = 100,
    validate_not_html = TRUE
  )

  meta <- jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
  files <- meta$files

  if (!length(files)) {
    stop("Figshare article ", article_id, " returned no files.", call. = FALSE)
  }

  nms <- vapply(files, function(z) as.character(z$name %||% NA_character_), character(1))

  if (!is.null(file_patterns) && length(file_patterns)) {
    keep <- rep(FALSE, length(nms))
    for (pat in file_patterns) {
      keep <- keep | grepl(pat, nms, ignore.case = TRUE)
    }
    files <- files[keep]
    nms   <- nms[keep]
  }

  if (!length(files)) {
    stop(
      "Figshare article ", article_id,
      " contained no files matching pattern(s): ",
      paste(file_patterns, collapse = ", "),
      call. = FALSE
    )
  }

  max_files <- suppressWarnings(as.integer(max_files))
  if (!is.finite(max_files) || max_files < 1L) max_files <- 1L
  files <- files[seq_len(min(length(files), max_files))]

  out <- character()

  for (z in files) {
    nm <- as.character(z$name)
    link <- as.character(z$download_url)

    if (!nzchar(link)) {
      stop("Figshare file has no usable download_url: ", nm, call. = FALSE)
    }

    dest <- file.path(cache_dir, nm)
    bf_download_cached(link, dest, force_refresh = force_refresh, quiet = quiet)
    out <- c(out, dest)
  }

  out
}

#' Discover and download selected files from a PANGAEA dataset
#'
#' Performs a conservative HTML-based discovery of binary/vector download links
#' from a PANGAEA DOI landing page, then downloads only links matching optional
#' file patterns. The helper fails clearly if no suitable links are found.
#'
#' @param doi PANGAEA DOI suffix, for example `"10.1594/PANGAEA.942325"`.
#' @param cache_dir Directory used for metadata and downloaded files.
#' @param force_refresh Logical; re-download files even if cached.
#' @param quiet Logical; suppress download messages where supported.
#' @param file_patterns Optional character vector of regular expressions used to
#'   select discovered links.
#' @param max_files Maximum number of matching files to download.
#'
#' @return Character vector of downloaded file paths.
#'
#' @family additional context overlay helpers
#' @keywords internal
#' @noRd
bf_download_pangaea_dataset_selected <- function(doi,
                                                 cache_dir,
                                                 force_refresh = FALSE,
                                                 quiet = TRUE,
                                                 file_patterns = NULL,
                                                 max_files = 1L) {
  # Best-effort resolver for PANGAEA datasets with binary-object links.
  # PANGAEA delivery can change, so this helper fails clearly when binary links
  # cannot be discovered.
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  html_path <- file.path(cache_dir, paste0("pangaea_", gsub("[^A-Za-z0-9]+", "_", doi), ".html"))
  page_url <- sprintf("https://doi.pangaea.de/%s", doi)

  bf_download_cached(
    page_url,
    html_path,
    force_refresh = force_refresh,
    quiet = quiet,
    mode = "wb",
    min_bytes = 100,
    validate_not_html = FALSE
  )

  txt <- paste(readLines(html_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

  # Extract href attributes. This is deliberately conservative.
  m <- gregexpr("href=[\"']([^\"']+)[\"']", txt, perl = TRUE)[[1]]
  if (identical(m[1], -1L)) {
    stop("Could not discover downloadable files from PANGAEA dataset: ", doi, call. = FALSE)
  }

  hrefs <- regmatches(txt, list(m))[[1]]
  hrefs <- sub("^href=[\"']", "", hrefs)
  hrefs <- sub("[\"']$", "", hrefs)
  hrefs <- unique(hrefs)

  # Keep plausible binary/vector archive links.
  hrefs <- hrefs[grepl("(download|Binary|binary|\\.gpkg|\\.zip|\\.gz|\\.csv)", hrefs, ignore.case = TRUE)]

  if (!length(hrefs)) {
    stop("PANGAEA dataset ", doi, " exposed no discoverable binary/vector download links.", call. = FALSE)
  }

  hrefs <- vapply(hrefs, function(h) {
    if (grepl("^https?://", h, ignore.case = TRUE)) return(h)
    if (startsWith(h, "/")) return(paste0("https://doi.pangaea.de", h))
    paste0("https://doi.pangaea.de/", h)
  }, character(1))

  if (!is.null(file_patterns) && length(file_patterns)) {
    keep <- rep(FALSE, length(hrefs))
    for (pat in file_patterns) {
      keep <- keep | grepl(pat, hrefs, ignore.case = TRUE)
    }
    hrefs <- hrefs[keep]
  }

  if (!length(hrefs)) {
    stop(
      "PANGAEA dataset ", doi,
      " contained no links matching pattern(s): ",
      paste(file_patterns, collapse = ", "),
      call. = FALSE
    )
  }

  max_files <- suppressWarnings(as.integer(max_files))
  if (!is.finite(max_files) || max_files < 1L) max_files <- 1L
  hrefs <- hrefs[seq_len(min(length(hrefs), max_files))]

  out <- character()

  for (u in hrefs) {
    stem <- basename(sub("[?#].*$", "", u))
    if (!nzchar(stem)) stem <- paste0("pangaea_", length(out) + 1L, ".zip")
    dest <- file.path(cache_dir, stem)
    bf_download_cached(u, dest, force_refresh = force_refresh, quiet = quiet)
    out <- c(out, dest)
  }

  out
}

# ------------------------------------------------------------------------------
# Vector source preparation
# ------------------------------------------------------------------------------
#' Collect candidate vector files from paths or extracted archives
#'
#' Expands directories and archives into a candidate set of supported vector
#' sources (`.gpkg`, `.shp`, GeoJSON/JSON, CSV, and file geodatabases). Optional
#' preference patterns are used to rank likely primary layers first.
#'
#' @param paths Character vector of files or directories.
#' @param force_refresh Logical; passed to archive unpacking.
#' @param quiet Logical; suppress unpacking messages where supported.
#' @param prefer_patterns Optional regular expressions used to rank candidate
#'   filenames.
#' @param cache_dir Optional cache directory used for archive extraction.
#' @param cache_stem Dataset/cache identifier used in extraction folder names.
#'
#' @return Character vector of candidate vector-file paths, ordered by preference.
#'
#' @keywords internal
#' @noRd
.bf_collect_vector_candidates <- function(paths,
                                          force_refresh = FALSE,
                                          quiet = TRUE,
                                          prefer_patterns = NULL,
                                          cache_dir = NULL,
                                          cache_stem = "vector") {
  candidates <- character()

  for (p in paths) {
    if (!file.exists(p)) next

    if (dir.exists(p)) {
      candidates <- c(
        candidates,
        list.files(p, pattern = "\\.(gpkg|shp|geojson|json|csv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
        list.dirs(p, recursive = TRUE, full.names = TRUE)[grepl("\\.gdb$", list.dirs(p, recursive = TRUE, full.names = TRUE), ignore.case = TRUE)]
      )
      next
    }

    if (grepl("\\.(gpkg|shp|geojson|json|csv)$", p, ignore.case = TRUE)) {
      candidates <- c(candidates, p)
      next
    }

    if (grepl("\\.(zip|tar\\.gz|tgz|gz)$", p, ignore.case = TRUE)) {
      if (is.null(cache_dir)) cache_dir <- dirname(p)

      exdir <- file.path(
        cache_dir,
        paste0(cache_stem, "_", tools::file_path_sans_ext(basename(p)), "_extracted")
      )

      bf_unpack_archive(p, exdir, force_refresh = force_refresh, quiet = quiet)

      candidates <- c(
        candidates,
        list.files(exdir, pattern = "\\.(gpkg|shp|geojson|json|csv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
        list.dirs(exdir, recursive = TRUE, full.names = TRUE)[grepl("\\.gdb$", list.dirs(exdir, recursive = TRUE, full.names = TRUE), ignore.case = TRUE)]
      )
    }
  }

  candidates <- unique(candidates)

  if (!length(candidates)) return(character(0))

  if (!is.null(prefer_patterns) && length(prefer_patterns)) {
    bn <- basename(candidates)
    score <- rep(0L, length(candidates))

    for (i in seq_along(prefer_patterns)) {
      score <- score +
        as.integer(grepl(prefer_patterns[i], bn, ignore.case = TRUE)) *
        (length(prefer_patterns) - i + 1L)
    }

    candidates <- candidates[order(score, decreasing = TRUE)]
  }

  candidates
}

#' Prepare a local vector source for an overlay loader
#'
#' Resolves a vector source in priority order: local source path, explicit source
#' URL, package-managed direct URL, Zenodo record, Figshare article, or PANGAEA
#' dataset. The selected file is downloaded/cached if needed, archives are
#' unpacked, and the best candidate vector file is returned.
#'
#' @param cache_dir Directory used for downloads and extracted files.
#' @param cache_stem Dataset/cache identifier used in generated filenames.
#' @param source_path Optional local file or directory override.
#' @param source_url Optional remote URL override.
#' @param force_refresh Logical; re-download/re-extract cached files.
#' @param quiet Logical; suppress messages where supported.
#' @param prefer_patterns Optional filename patterns used to rank vector files.
#' @param default_source_url Optional package-managed direct download URL.
#' @param default_zenodo_record Optional Zenodo record identifier.
#' @param default_zenodo_file_patterns Optional file-selection patterns for
#'   Zenodo metadata.
#' @param default_figshare_article Optional Figshare article identifier.
#' @param default_figshare_file_patterns Optional file-selection patterns for
#'   Figshare metadata.
#' @param default_pangaea_doi Optional PANGAEA DOI suffix.
#' @param default_pangaea_file_patterns Optional link-selection patterns for
#'   PANGAEA discovery.
#' @param max_default_files Maximum number of package-managed files to download.
#'
#' @return Path to a selected local vector file.
#'
#' @keywords internal
#' @noRd
.bf_prepare_vector_source <- function(cache_dir,
                                      cache_stem,
                                      source_path = NULL,
                                      source_url = NULL,
                                      force_refresh = FALSE,
                                      quiet = TRUE,
                                      prefer_patterns = NULL,
                                      default_source_url = NULL,
                                      default_zenodo_record = NULL,
                                      default_zenodo_file_patterns = NULL,
                                      default_figshare_article = NULL,
                                      default_figshare_file_patterns = NULL,
                                      default_pangaea_doi = NULL,
                                      default_pangaea_file_patterns = NULL,
                                      max_default_files = 1L) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- character()

  # Developer/local mirror override.
  if (!is.null(source_path) && file.exists(source_path)) {
    paths <- c(paths, source_path)
  }

  # Developer/source mirror URL override.
  if (!length(paths) && !is.null(source_url) && nzchar(as.character(source_url)[1])) {
    ext <- .bf_url_ext(source_url, default = "zip")
    dest <- file.path(cache_dir, paste0(cache_stem, "_download.", ext))
    bf_download_cached(source_url, dest, force_refresh = force_refresh, quiet = quiet)
    paths <- c(paths, dest)
  }

  # Package-managed direct URL.
  if (!length(paths) && !is.null(default_source_url) && nzchar(as.character(default_source_url)[1])) {
    ext <- .bf_url_ext(default_source_url, default = "zip")
    dest <- file.path(cache_dir, paste0(cache_stem, "_download.", ext))
    bf_download_cached(default_source_url, dest, force_refresh = force_refresh, quiet = quiet)
    paths <- c(paths, dest)
  }

  # Package-managed Zenodo metadata selection.
  if (!length(paths) && !is.null(default_zenodo_record)) {
    paths <- bf_download_zenodo_record_selected(
      record_id = default_zenodo_record,
      cache_dir = cache_dir,
      force_refresh = force_refresh,
      quiet = quiet,
      file_patterns = default_zenodo_file_patterns,
      max_files = max_default_files
    )
  }

  # Package-managed Figshare metadata selection.
  if (!length(paths) && !is.null(default_figshare_article)) {
    paths <- bf_download_figshare_article_selected(
      article_id = default_figshare_article,
      cache_dir = cache_dir,
      force_refresh = force_refresh,
      quiet = quiet,
      file_patterns = default_figshare_file_patterns,
      max_files = max_default_files
    )
  }

  # Package-managed PANGAEA metadata selection.
  if (!length(paths) && !is.null(default_pangaea_doi)) {
    paths <- bf_download_pangaea_dataset_selected(
      doi = default_pangaea_doi,
      cache_dir = cache_dir,
      force_refresh = force_refresh,
      quiet = quiet,
      file_patterns = default_pangaea_file_patterns,
      max_files = max_default_files
    )
  }

  if (!length(paths)) {
    .bf_unsupported_package_managed(
      cache_stem,
      "No stable single-file public download has been configured for this dataset."
    )
  }

  candidates <- .bf_collect_vector_candidates(
    paths = paths,
    force_refresh = force_refresh,
    quiet = quiet,
    prefer_patterns = prefer_patterns,
    cache_dir = cache_dir,
    cache_stem = cache_stem
  )

  if (!length(candidates)) {
    stop(cache_stem, " download completed, but no supported vector file was found.", call. = FALSE)
  }

  candidates[[1]]
}

#' Read a vector file, optionally selecting a preferred layer
#'
#' For multi-layer sources such as GeoPackage or file geodatabase inputs, this
#' helper selects the layer whose name best matches `layer_patterns`. Otherwise
#' it delegates to the package's generic vector reader.
#'
#' @param path Path to a vector source.
#' @param layer_patterns Optional regular expressions used to rank layer names.
#' @param quiet Logical; suppress read messages where supported.
#'
#' @return An `sf` object read from `path`.
#'
#' @keywords internal
#' @noRd
.bf_read_vector_layer_any <- function(path, layer_patterns = NULL, quiet = TRUE) {
  bf_require_packages(c("sf"), context = "additional context overlay")

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    return(bf_read_vector_any(path, quiet = quiet))
  }

  if (!is.null(layer_patterns) && length(layer_patterns) &&
      (grepl("\\.gpkg$", path, ignore.case = TRUE) || grepl("\\.gdb$", path, ignore.case = TRUE))) {
    lyr <- tryCatch(sf::st_layers(path)$name, error = function(e) character())

    if (length(lyr)) {
      score <- rep(0L, length(lyr))
      for (i in seq_along(layer_patterns)) {
        score <- score +
          as.integer(grepl(layer_patterns[i], lyr, ignore.case = TRUE)) *
          (length(layer_patterns) - i + 1L)
      }

      if (any(score > 0)) {
        return(sf::read_sf(path, layer = lyr[which.max(score)], quiet = isTRUE(quiet)))
      }
    }
  }

  bf_read_vector_any(path, quiet = quiet)
}

# Read CSV tables that contain coordinates under dataset-specific names, then
# promote them to sf point objects. This is deliberately separate from
# bf_read_vector_any(), because several biofetchR overlays use recognised but
# non-generic coordinate names such as HydroWASTE's LON_WWTP/LAT_WWTP.
#' Pick the first matching column from candidate names
#'
#' Matching is case-insensitive and returns the original column name from the
#' input object.
#'
#' @param x Data frame or `sf` object.
#' @param candidates Character vector of possible column names.
#'
#' @return Matching column name, or `NULL` when no candidate is present.
#'
#' @keywords internal
#' @noRd
.bf_pick_col <- function(x, candidates) {
  nms <- names(x)
  low <- tolower(nms)
  for (cand in candidates) {
    hit <- which(low == tolower(cand))
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}

#' Read a coordinate table as an `sf` point object
#'
#' Reads a CSV file, identifies longitude and latitude columns from candidate
#' names, removes rows with invalid coordinates, and converts the result to WGS84
#' points. This supports point-based overlay datasets whose coordinate columns
#' are not named `decimalLongitude` and `decimalLatitude`.
#'
#' @param path Path to a CSV file.
#' @param lon_candidates Candidate longitude column names.
#' @param lat_candidates Candidate latitude column names.
#' @param quiet Logical; suppress read progress where supported.
#'
#' @return WGS84 `sf` point object.
#'
#' @keywords internal
#' @noRd
.bf_read_csv_points_as_sf <- function(path,
                                      lon_candidates,
                                      lat_candidates,
                                      quiet = TRUE) {
  bf_require_packages(c("sf"), context = "additional context overlay")

  dat <- if (requireNamespace("readr", quietly = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE, progress = !isTRUE(quiet))
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (!nrow(dat)) {
    stop("CSV contains no rows: ", path, call. = FALSE)
  }

  lon_col <- .bf_pick_col(dat, lon_candidates)
  lat_col <- .bf_pick_col(dat, lat_candidates)

  if (is.null(lon_col) || is.null(lat_col)) {
    stop(
      "CSV lacks recognised lon/lat columns: ", path,
      ". Available columns: ", paste(names(dat), collapse = ", "),
      call. = FALSE
    )
  }

  dat[[lon_col]] <- suppressWarnings(as.numeric(dat[[lon_col]]))
  dat[[lat_col]] <- suppressWarnings(as.numeric(dat[[lat_col]]))

  dat <- dat[is.finite(dat[[lon_col]]) & is.finite(dat[[lat_col]]), , drop = FALSE]
  dat <- dat[dat[[lon_col]] >= -180 & dat[[lon_col]] <= 180 & dat[[lat_col]] >= -90 & dat[[lat_col]] <= 90, , drop = FALSE]

  if (!nrow(dat)) {
    stop("CSV has recognised coordinate columns but no valid coordinate rows: ", path, call. = FALSE)
  }

  sf::st_as_sf(dat, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
}

#' Read HydroWASTE CSV coordinates as points
#'
#' HydroWASTE mirrors may provide wastewater treatment plant coordinates or
#' outfall coordinates under several possible field names. This helper prefers
#' treatment-plant coordinates and falls back to outfall coordinates.
#'
#' @param path Path to a HydroWASTE CSV file.
#' @param quiet Logical; suppress read progress where supported.
#'
#' @return WGS84 `sf` point object.
#'
#' @keywords internal
#' @noRd
.bf_read_hydrowaste_csv_as_sf <- function(path, quiet = TRUE) {
  # Prefer the WWTP coordinates because this loader represents treatment plants.
  # Fall back to outfall coordinates if a mirror only provides those fields.
  .bf_read_csv_points_as_sf(
    path = path,
    lon_candidates = c("LON_WWTP", "lon_wwtp", "LONGITUDE", "longitude", "decimalLongitude", "decimal_longitude", "LON_OUT", "lon_out"),
    lat_candidates = c("LAT_WWTP", "lat_wwtp", "LATITUDE", "latitude", "decimalLatitude", "decimal_latitude", "LAT_OUT", "lat_out"),
    quiet = quiet
  )
}

# ------------------------------------------------------------------------------
# Public overlay loaders
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# HydroWASTE wastewater treatment plants
# ------------------------------------------------------------------------------
#' Load HydroWASTE wastewater treatment plant points
#'
#' Downloads, caches, reads, and standardises the HydroWASTE wastewater treatment
#' plant layer for use as a terrestrial/freshwater context overlay. The returned
#' object contains stable columns used by the GBIF processing pipeline:
#' `hydrowaste_id`, `hydrowaste_name`, `hydrowaste_level`, and geometry.
#'
#' HydroWASTE is treated as a point overlay. If the source is not already point
#' geometry, centroids are used so occurrence points can be assigned to the
#' nearest treatment-plant feature by the main pipeline.
#'
#' @param cache_dir Directory used for downloads, extracted files, and cached
#'   source data. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param force_refresh Logical; re-download and re-read cached files.
#' @param iso2c Optional ISO2 country code vector used to crop the layer.
#' @param clip_to_countries Logical; if `TRUE` and `iso2c` is supplied, crop the
#'   returned features to those countries.
#' @param source_path Optional local source-file or directory override.
#' @param source_url Optional remote source URL override.
#' @param quiet Logical; suppress messages where supported.
#'
#' @section Data access and attribution:
#' This loader can download HydroWASTE v1.0 from a package-managed public mirror
#' or read a user-supplied local/source URL. HydroWASTE should be cited using the
#' original dataset publication and any provider guidance associated with the
#' exact file/version used. biofetchR caches the source locally for reproducible
#' processing but does not redistribute HydroWASTE.
#'
#' @return An `sf` point object in WGS84 longitude/latitude coordinates
#'   containing HydroWASTE wastewater-treatment plant features. The returned
#'   object includes `hydrowaste_id`, a stable treatment-plant identifier;
#'   `hydrowaste_name`, a human-readable plant or feature label;
#'   `hydrowaste_level`, the treatment-level attribute where available; and the
#'   active geometry column. If the source layer is not already point geometry,
#'   centroids are returned so occurrence points can be assigned to the nearest
#'   HydroWASTE feature by downstream workflows. If `clip_to_countries = TRUE`
#'   and `iso2c` is supplied, the returned object is cropped to the requested
#'   countries where possible.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   hydrowaste <- bf_load_hydrowaste(
#'     iso2c = c("GB", "IE"),
#'     cache_dir = file.path(tempdir(), "biofetchR_hydrowaste"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   hydrowaste
#' }
#' }
#'
#' @family additional context overlay loaders
#' @md
#' @export
bf_load_hydrowaste <- function(cache_dir = NULL,
                               force_refresh = FALSE,
                               iso2c = NULL,
                               clip_to_countries = TRUE,
                               source_path = NULL,
                               source_url = NULL,
                               quiet = TRUE) {
  bf_require_packages(c("sf"), context = "additional context overlay")

  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly. In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  cache_dir <- normalizePath(
    as.character(cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  path <- .bf_prepare_vector_source(
    cache_dir = cache_dir,
    cache_stem = "hydrowaste",
    source_path = source_path,
    source_url = source_url,
    force_refresh = force_refresh,
    quiet = quiet,
    prefer_patterns = c("HydroWASTE_v10", "hydrowaste", "wastewater", "wwtp", "outfall"),
    default_source_url = "https://ckan.rdas.live/dataset/0fb578e3-aeb6-4acb-997e-b5ee6c6127e1/resource/2a3f32eb-1ea4-426f-931c-555941c3448b/download/hydrowaste_v10.zip"
  )

  x <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    .bf_read_hydrowaste_csv_as_sf(path, quiet = quiet)
  } else {
    .bf_read_vector_layer_any(path, quiet = quiet)
  }

  x <- bf_add_field_from_candidates(x, "hydrowaste_id", c("hydrowaste_id", "WASTE_ID", "waste_id", "WWTP_ID", "wwtp_id", "OBJECTID", "id", "ID"), default = as.character(seq_len(nrow(x))))
  x <- bf_add_field_from_candidates(x, "hydrowaste_name", c("hydrowaste_name", "WWTP_NAME", "wwtp_name", "NAME", "name", "plant_name"), default = x$hydrowaste_id)
  x <- bf_add_field_from_candidates(x, "hydrowaste_level", c("treatment", "TREATMENT", "treatment_level", "LEVEL"), default = NA_character_)

  x <- bf_sf_make_valid(bf_sf_wgs84(x), quiet = quiet)

  if (!all(grepl("POINT", as.character(unique(sf::st_geometry_type(x))), ignore.case = TRUE))) {
    x <- suppressWarnings(sf::st_centroid(x))
  }

  gcol <- attr(x, "sf_column")
  keep <- intersect(c("hydrowaste_id", "hydrowaste_name", "hydrowaste_level", gcol), names(x))
  x <- x[, keep, drop = FALSE]

  if (isTRUE(clip_to_countries) && length(iso2c)) x <- bf_clip_to_countries(x, iso2c, quiet = quiet)
  x
}

# ------------------------------------------------------------------------------
# GDW reservoir polygons
# ------------------------------------------------------------------------------
#' Load Global Dam Watch reservoir polygons
#'
#' Downloads, caches, reads, and standardises Global Dam Watch reservoir polygon
#' data for use as a terrestrial/freshwater context overlay. The returned object
#' contains stable columns used by the GBIF processing pipeline:
#' `gdw_reservoir_id`, `gdw_reservoir_name`, and geometry.
#'
#' @param cache_dir Directory used for downloads, extracted files, and cached
#'   source data. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param force_refresh Logical; re-download and re-read cached files.
#' @param iso2c Optional ISO2 country code vector used to crop the layer.
#' @param clip_to_countries Logical; if `TRUE` and `iso2c` is supplied, crop the
#'   returned features to those countries.
#' @param source_path Optional local source-file or directory override.
#' @param source_url Optional remote source URL override.
#' @param quiet Logical; suppress messages where supported.
#'
#' @section Data access and attribution:
#' This loader can download Global Dam Watch reservoir data from provider-hosted
#' metadata/download services or read a user-supplied local/source URL. Users
#' should cite the Global Dam Watch data release and associated publication for
#' the exact version used. biofetchR caches the source locally for reproducible
#' processing but does not redistribute Global Dam Watch data.
#'
#' @return An `sf` polygon object in WGS84 longitude/latitude coordinates
#'   containing Global Dam Watch reservoir features. The returned object includes
#'   `gdw_reservoir_id`, a stable reservoir identifier;
#'   `gdw_reservoir_name`, a human-readable reservoir or dam label where
#'   available; and the active geometry column. Each row represents one reservoir
#'   polygon or multipolygon used for assigning occurrence records to reservoir
#'   context. If `clip_to_countries = TRUE` and `iso2c` is supplied, the returned
#'   object is cropped to the requested countries where possible.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   reservoirs <- bf_load_gdw_reservoirs(
#'     iso2c = "GB",
#'     cache_dir = file.path(tempdir(), "biofetchR_gdw_reservoirs"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   reservoirs
#' }
#' }
#'
#' @family additional context overlay loaders
#' @md
#' @export
bf_load_gdw_reservoirs <- function(cache_dir = NULL,
                                   force_refresh = FALSE,
                                   iso2c = NULL,
                                   clip_to_countries = TRUE,
                                   source_path = NULL,
                                   source_url = NULL,
                                   quiet = TRUE) {
  bf_require_packages(c("sf"), context = "additional context overlay")

  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly. In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  cache_dir <- normalizePath(
    as.character(cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  path <- .bf_prepare_vector_source(
    cache_dir = cache_dir,
    cache_stem = "gdw_reservoirs",
    source_path = source_path,
    source_url = source_url,
    force_refresh = force_refresh,
    quiet = quiet,
    prefer_patterns = c("reservoir", "Reservoir", "GDW_v1_0_shp", "gdw", "waterbody"),
    default_figshare_article = "25988293",
    default_figshare_file_patterns = c("GDW_v1_0_shp\\.zip$"),
    max_default_files = 1L
  )

  x <- .bf_read_vector_layer_any(path, layer_patterns = c("reservoir", "waterbody"), quiet = quiet)
  x <- bf_add_field_from_candidates(x, "gdw_reservoir_id", c("gdw_reservoir_id", "RESERVOIR_ID", "reservoir_id", "GRAND_ID", "GDW_ID", "OBJECTID", "id"), default = as.character(seq_len(nrow(x))))
  x <- bf_add_field_from_candidates(x, "gdw_reservoir_name", c("gdw_reservoir_name", "RES_NAME", "res_name", "DAM_NAME", "NAME", "name"), default = x$gdw_reservoir_id)

  x <- bf_sf_make_valid(bf_sf_wgs84(x), quiet = quiet)
  gcol <- attr(x, "sf_column")
  keep <- intersect(c("gdw_reservoir_id", "gdw_reservoir_name", gcol), names(x))
  x <- x[, keep, drop = FALSE]

  if (isTRUE(clip_to_countries) && length(iso2c)) x <- bf_clip_to_countries(x, iso2c, quiet = quiet)
  x
}

# ------------------------------------------------------------------------------
# Global mining polygons
# ------------------------------------------------------------------------------
#' Load global mining polygons
#'
#' Downloads, caches, reads, and standardises global mining polygon data for use
#' as a terrestrial context overlay. The returned object contains stable columns
#' used by the GBIF processing pipeline: `global_mining_id`,
#' `global_mining_name`, and geometry.
#'
#' @param cache_dir Directory used for downloads, extracted files, and cached
#'   source data. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param force_refresh Logical; re-download and re-read cached files.
#' @param iso2c Optional ISO2 country code vector used to crop the layer.
#' @param clip_to_countries Logical; if `TRUE` and `iso2c` is supplied, crop the
#'   returned features to those countries.
#' @param source_path Optional local source-file or directory override.
#' @param source_url Optional remote source URL override.
#' @param quiet Logical; suppress messages where supported.
#'
#' @section Data access and attribution:
#' This loader can download the global-scale mining polygons from PANGAEA or read
#' a user-supplied local/source URL. Users should cite the PANGAEA dataset DOI,
#' version and associated publication and follow the dataset licence attached to
#' the downloaded file. biofetchR caches the source locally for reproducible
#' processing but does not redistribute the mining dataset.
#'
#' @return An `sf` polygon object in WGS84 longitude/latitude coordinates
#'   containing global mining-area features. The returned object includes
#'   `global_mining_id`, a stable mining-feature identifier;
#'   `global_mining_name`, a human-readable mine, country or feature label where
#'   available; and the active geometry column. Each row represents one mining
#'   polygon or multipolygon used for assigning occurrence records to mining-area
#'   context. If `clip_to_countries = TRUE` and `iso2c` is supplied, the returned
#'   object is cropped to the requested countries where possible.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   mines <- bf_load_global_mining(
#'     iso2c = "ZA",
#'     cache_dir = file.path(tempdir(), "biofetchR_global_mining"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   mines
#' }
#' }
#'
#' @family additional context overlay loaders
#' @md
#' @export
bf_load_global_mining <- function(cache_dir = NULL,
                                  force_refresh = FALSE,
                                  iso2c = NULL,
                                  clip_to_countries = TRUE,
                                  source_path = NULL,
                                  source_url = NULL,
                                  quiet = TRUE) {
  bf_require_packages(c("sf"), context = "additional context overlay")

  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly. In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  cache_dir <- normalizePath(
    as.character(cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  path <- .bf_prepare_vector_source(
    cache_dir = cache_dir,
    cache_stem = "global_mining",
    source_path = source_path,
    source_url = source_url,
    force_refresh = force_refresh,
    quiet = quiet,
    prefer_patterns = c("global_mining_polygons_v2", "global_mining_polygons", "mining_polygons", "mining", "mine", "global", "maus"),
    default_source_url = "https://download.pangaea.de/dataset/942325/files/global_mining_polygons_v2.gpkg",
    max_default_files = 1L
  )

  x <- .bf_read_vector_layer_any(path, layer_patterns = c("mining", "polygon", "global"), quiet = quiet)
  x <- bf_add_field_from_candidates(x, "global_mining_id", c("global_mining_id", "MINE_ID", "mine_id", "FID", "OBJECTID", "id"), default = as.character(seq_len(nrow(x))))
  x <- bf_add_field_from_candidates(x, "global_mining_name", c("global_mining_name", "NAME", "name", "mine_name", "COUNTRY_NAME", "ISO3_CODE"), default = x$global_mining_id)

  x <- bf_sf_make_valid(bf_sf_wgs84(x), quiet = quiet)
  gcol <- attr(x, "sf_column")
  keep <- intersect(c("global_mining_id", "global_mining_name", gcol), names(x))
  x <- x[, keep, drop = FALSE]

  if (isTRUE(clip_to_countries) && length(iso2c)) x <- bf_clip_to_countries(x, iso2c, quiet = quiet)
  x
}
