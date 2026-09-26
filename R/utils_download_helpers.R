################################################################################
# utils_download_helpers.R
# ------------------------------------------------------------------------------
# biofetchR download and archive utility helpers
# ------------------------------------------------------------------------------
#
# This script contains low-level file/download helpers used by multiple biofetchR
# overlay and context-loading functions. These helpers are deliberately generic:
# they do not know anything about GBIF, GADM, SInAS, TEOW, FEOW, Marine Regions,
# WDPA, Ramsar, or any other specific data provider. Instead, higher-level
# loader functions call these utilities whenever they need to:
#
#   1. check whether a cached file already exists and is large enough to reuse;
#   2. inspect the first bytes of a file to confirm whether it is a real archive;
#   3. reject downloaded HTML landing pages, login pages, or error pages;
#   4. download large files robustly, with retries and optional resume support;
#   5. unpack ZIP, tar.gz, tgz, and single-file gzip archives;
#   6. recursively unpack nested ZIP archives inside provider downloads.
#
# Why this matters:
#   Spatial biodiversity data providers often serve large archives, redirected
#   URLs, HTML error pages, or nested zip files. Without validation, a failed
#   download can silently leave an HTML file or tiny partial file in the cache,
#   which later spatial loaders would try to read as a shapefile/geopackage.
#   These helpers centralise those checks so provider-specific loaders stay
#   simpler and safer.
#
# Main helper groups:
#   File validation:
#     bf_file_size()
#     bf_read_magic()
#     bf_is_probably_html()
#     bf_is_zip_file()
#     bf_is_gzip_file()
#
#   Downloading:
#     bf_download_cached()
#     bf_download_file()
#
#   Archive extraction:
#     bf_unzip_cached()
#     bf_unpack_archive()
#     bf_unzip_nested()
#
################################################################################

# -----------------------------------------------------------------------------
# File-inspection helpers
# -----------------------------------------------------------------------------
# The first group of helpers performs cheap checks on local files. These are
# intentionally small and dependency-light because they are used by the download
# and archive-unpacking functions below before any heavier GIS-reading steps are
# attempted.

#' Return file size safely
#'
#' Internal helper that returns zero for missing files or unreadable file-size
#' metadata.
#'
#' @param path Character path to a file.
#'
#' @return Numeric file size in bytes, or `0` when unavailable.
#' @keywords internal
#' @noRd
bf_file_size <- function(path) {
  if (!file.exists(path)) return(0)
  sz <- suppressWarnings(file.info(path)$size)
  if (is.na(sz)) 0 else as.numeric(sz)
}

#' Read file signature bytes
#'
#' Reads the first bytes of a file for lightweight archive and HTML validation.
#'
#' @param path Character path to a file.
#' @param n Integer number of bytes to read.
#'
#' @return Raw vector of up to `n` bytes.
#' @keywords internal
#' @noRd
bf_read_magic <- function(path, n = 512L) {
  if (!file.exists(path) || bf_file_size(path) <= 0) return(raw())
  con <- file(path, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  readBin(con, what = raw(), n = n)
}

#' Detect likely HTML error pages
#'
#' Used after downloads to distinguish real GIS archives from provider landing,
#' login or error pages saved to disk.
#'
#' @param path Character path to a downloaded file.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
bf_is_probably_html <- function(path) {
  magic <- bf_read_magic(path, n = 512L)
  if (!length(magic)) return(FALSE)

  ints <- as.integer(magic)
  ints <- ints[ints >= 9L & ints <= 126L]
  if (!length(ints)) return(FALSE)

  txt <- paste(intToUtf8(ints, multiple = TRUE), collapse = "")
  txt <- tolower(trimws(txt))

  grepl("^<!doctype html|^<html|<head|<body", txt)
}

#' Detect ZIP archives from file signature bytes
#'
#' @param path Character path to a file.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
bf_is_zip_file <- function(path) {
  magic <- bf_read_magic(path, n = 4L)
  if (length(magic) < 4L) return(FALSE)
  sig <- paste(as.character(magic[1:4]), collapse = " ")
  sig %in% c("50 4b 03 04", "50 4b 05 06", "50 4b 07 08")
}

#' Detect gzip archives from file signature bytes
#'
#' @param path Character path to a file.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
bf_is_gzip_file <- function(path) {
  magic <- bf_read_magic(path, n = 2L)
  length(magic) >= 2L && identical(as.integer(magic[1:2]), c(31L, 139L))
}

# -----------------------------------------------------------------------------
# Robust cached downloader
# -----------------------------------------------------------------------------
# This is the central download helper used by most provider loaders. It handles
# cache reuse, timeout management, retry attempts, validation of downloaded files,
# and optional partial-file resume for large archives.

#' Download a remote file into a local cache
#'
#' Robust internal downloader with cache reuse, optional resume support, retry
#' handling and basic validation against tiny files or HTML error pages.
#'
#' @param url Character URL to download.
#' @param dest Character destination path in the cache.
#' @param force_refresh Logical. If `TRUE`, ignore any completed cached file.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#' @param mode Download mode passed to download backends.
#' @param timeout Numeric timeout in seconds.
#' @param retries Integer number of retries after the first attempt.
#' @param min_bytes Minimum acceptable downloaded file size.
#' @param validate_not_html Logical. If `TRUE`, reject likely HTML files.
#' @param resume Logical. If `TRUE`, try to resume partial downloads.
#' @param keep_partial Logical. If `TRUE`, retain valid partial downloads for a
#'   future resume attempt.
#'
#' @return Character path to the completed cached file.
#' @keywords internal
#' @noRd
bf_download_cached <- function(url, dest, force_refresh = FALSE, quiet = TRUE, mode = "wb",
                               timeout = as.numeric(Sys.getenv("BF_DOWNLOAD_TIMEOUT", unset = "7200")),
                               retries = as.integer(Sys.getenv("BF_DOWNLOAD_RETRIES", unset = "2")),
                               min_bytes = 1024,
                               validate_not_html = TRUE,
                               resume = TRUE,
                               keep_partial = TRUE) {
  stopifnot(is.character(url), length(url) == 1, nzchar(url))
  stopifnot(is.character(dest), length(dest) == 1, nzchar(dest))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

  timeout <- suppressWarnings(as.numeric(timeout))
  if (!is.finite(timeout) || timeout < 60) timeout <- 7200

  retries <- suppressWarnings(as.integer(retries))
  if (!is.finite(retries) || retries < 0) retries <- 2L

  min_bytes <- suppressWarnings(as.numeric(min_bytes))
  if (!is.finite(min_bytes) || min_bytes < 1) min_bytes <- 1

  # A completed cached file is only reused if it exists, is not too small, and
  # does not look like an HTML error/landing page. This prevents broken provider
  # responses being treated as valid GIS data in later functions.
  cache_ok <- file.exists(dest) &&
    bf_file_size(dest) >= min_bytes &&
    (!isTRUE(validate_not_html) || !bf_is_probably_html(dest))

  if (cache_ok && !isTRUE(force_refresh)) return(dest)

  if (file.exists(dest) && (!cache_ok || isTRUE(force_refresh))) {
    try(unlink(dest, force = TRUE), silent = TRUE)
  }

  tmp <- paste0(dest, ".part")

  # If force_refresh is TRUE we remove the completed cache, but we do not
  # automatically remove a .part file. For very large archives, keeping the
  # partial download allows resume across repeated calls.
  if (file.exists(tmp) && !isTRUE(resume)) {
    try(unlink(tmp, force = TRUE), silent = TRUE)
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(as.numeric(old_timeout), timeout, na.rm = TRUE))

  last_msg <- NULL
  attempts <- seq_len(retries + 1L)

  .system2_ok <- function(status) {
    exit_status <- attr(
      status,
      "status",
      exact = TRUE
    )

    # With stdout/stderr capture, system2() returns character output.
    # A non-zero exit code is stored in attr(x, "status"); successful
    # captured commands have no status attribute.
    if (!is.null(exit_status)) {
      return(
        identical(
          as.integer(exit_status),
          0L
        )
      )
    }

    # Without output capture, system2() returns the integer exit code.
    if (is.numeric(status) && length(status) == 1L) {
      return(
        identical(
          as.integer(status),
          0L
        )
      )
    }

    # Character output with no status attribute means command success.
    is.character(status)
  }

  # Move the temporary .part file into the final cache location. On Windows,
  # file.rename() can fail across devices or locked paths, so we fall back to
  # file.copy() and then remove the temporary file.
  .move_tmp_to_dest <- function() {
    if (file.exists(dest)) try(unlink(dest, force = TRUE), silent = TRUE)
    renamed <- file.rename(tmp, dest)
    if (!isTRUE(renamed)) {
      ok_copy <- file.copy(tmp, dest, overwrite = TRUE)
      if (isTRUE(ok_copy)) try(unlink(tmp, force = TRUE), silent = TRUE)
      if (!isTRUE(ok_copy)) stop("Downloaded file could not be moved into cache: ", dest, call. = FALSE)
    }
    dest
  }

  # Validate the temporary file before accepting it as a completed download.
  # This repeats the same size/HTML checks used for completed cached files.
  .tmp_valid_enough <- function() {
    file.exists(tmp) &&
      bf_file_size(tmp) >= min_bytes &&
      (!isTRUE(validate_not_html) || !bf_is_probably_html(tmp))
  }

  # Prefer curl-backed downloading when available. For resumed downloads, system
  # curl is used because `curl -C -` is generally more reliable for very large
  # interrupted files than base R's download.file().
  .download_with_rcurl_resume <- function() {
    if (!requireNamespace("curl", quietly = TRUE)) return(FALSE)

    bytes_existing <- bf_file_size(tmp)

    # For resumed downloads, use system curl because it handles "-C -" robustly
    # for very large files.
    if (isTRUE(resume) && bytes_existing > 0) {
      curl_bin <- Sys.which("curl")
      if (!nzchar(curl_bin)) return(FALSE)

      args <- c(
        "-L",
        "--fail",
        "--retry", "2",
        "--retry-delay", "5",
        "--connect-timeout", "60",
        "--max-time", as.character(timeout),
        "-A", "biofetchR",
        "-C", as.character(bytes_existing),
        "-o", shQuote(tmp),
        shQuote(url)
      )

      status <- system2(
        curl_bin,
        args = args,
        stdout = if (isTRUE(quiet)) TRUE else "",
        stderr = if (isTRUE(quiet)) TRUE else ""
      )

      return(.system2_ok(status))
    }

    # Fresh download: curl package is fine here.
    h <- curl::new_handle()
    curl::handle_setopt(
      h,
      connecttimeout = max(60, min(timeout, 600)),
      timeout = timeout,
      followlocation = TRUE,
      failonerror = TRUE,
      useragent = "biofetchR"
    )

    curl::curl_download(
      url = url,
      destfile = tmp,
      quiet = isTRUE(quiet),
      mode = mode,
      handle = h
    )

    TRUE
  }

  # Fallback resume-capable downloader using the system curl binary. This path is
  # used when the curl R package is unavailable but command-line curl exists.
  .download_with_system_curl_resume <- function() {
    curl_bin <- Sys.which("curl")
    if (!nzchar(curl_bin)) return(FALSE)

    args <- c(
      "-L",
      "--fail",
      "--retry", "2",
      "--retry-delay", "5",
      "--connect-timeout", "60",
      "--max-time", as.character(timeout),
      "-A", "biofetchR",
      "-o", shQuote(tmp)
    )

    if (isTRUE(resume) && file.exists(tmp) && bf_file_size(tmp) > 0) {
      args <- c(args, "-C", as.character(bf_file_size(tmp)))
    }

    args <- c(args, shQuote(url))

    status <- system2(
      curl_bin,
      args = args,
      stdout = if (isTRUE(quiet)) TRUE else "",
      stderr = if (isTRUE(quiet)) TRUE else ""
    )

    .system2_ok(status)
  }

  # Final dependency-free fallback. This does not reliably resume large files, but
  # keeps the package usable on systems without curl.
  .download_with_utils <- function() {
    # utils::download.file cannot reliably resume large interrupted downloads,
    # but remains a dependency-free fallback.
    if (file.exists(tmp) && (!isTRUE(resume) || bf_file_size(tmp) == 0)) {
      try(unlink(tmp, force = TRUE), silent = TRUE)
    }

    status <- utils::download.file(
      url,
      destfile = tmp,
      mode = mode,
      quiet = isTRUE(quiet),
      method = "libcurl"
    )

    is.null(status) || identical(status, 0L)
  }

  # Download loop: try the best available backend, validate the temporary file,
  # then either promote it to the final cache path or keep/remove the partial file
  # depending on whether it looks resumable.
  for (attempt in attempts) {
    part_size <- bf_file_size(tmp)

    if (!isTRUE(quiet)) {
      message("Downloading: ", url)
      message("  -> ", dest)
      message("  attempt ", attempt, " of ", length(attempts), "; timeout = ", getOption("timeout"), " sec")
      if (part_size > 0) {
        message("  resuming from partial file: ", round(part_size / 1024^2, 2), " MB")
      }
    }

    ok <- tryCatch({
      # Prefer the system curl executable when available. On Windows, the R {curl}
      # backend can occasionally fail after a completed download when it tries to
      # rename an internal ".curltmp" file to the package ".part" file. System curl
      # writes directly to `tmp` through "-o", avoiding that rename failure.
      if (nzchar(Sys.which("curl"))) {
        .download_with_system_curl_resume()
      } else if (isTRUE(resume) && requireNamespace("curl", quietly = TRUE)) {
        .download_with_rcurl_resume()
      } else {
        .download_with_utils()
      }
    }, warning = function(w) {
      last_msg <<- conditionMessage(w)
      FALSE
    }, error = function(e) {
      last_msg <<- conditionMessage(e)
      FALSE
    })

    if (ok && .tmp_valid_enough()) {
      .move_tmp_to_dest()
      return(dest)
    }

    # If the partial file is HTML or too tiny, it is probably an error page rather
    # than a resumable archive. Remove it. Otherwise keep it for resume.
    if (file.exists(tmp)) {
      tiny_or_html <- bf_file_size(tmp) < min_bytes ||
        (isTRUE(validate_not_html) && bf_is_probably_html(tmp))
      if (tiny_or_html || !isTRUE(keep_partial)) {
        try(unlink(tmp, force = TRUE), silent = TRUE)
      }
    }

    if (!isTRUE(quiet)) {
      message("Download attempt failed", if (!is.null(last_msg)) paste0(": ", last_msg) else ".")
      if (file.exists(tmp) && bf_file_size(tmp) >= min_bytes) {
        message("Keeping partial file for resume: ", tmp,
                " (", round(bf_file_size(tmp) / 1024^2, 2), " MB)")
      }
    }
  }

  stop("Download failed for ", url,
       if (!is.null(last_msg)) paste0(": ", last_msg) else "",
       ". A partial file may have been kept at: ", tmp,
       ". Re-run the same command to resume, or delete the cache and try again.",
       call. = FALSE)
}

# -----------------------------------------------------------------------------
# One-off download wrapper
# -----------------------------------------------------------------------------

#' Download a file without reusing an existing completed cache
#'
#' Thin wrapper around `bf_download_cached()` for temporary or one-off downloads.
#' Unlike `bf_download_cached()`, this helper always refreshes the destination
#' file and disables resume/partial-file retention by default. It is useful where
#' older code expected a simple `.bf_download_file()` helper but we still want to
#' use the package's central download validation logic.
#'
#' @param url Character. Remote URL.
#' @param destfile Character. Local destination path.
#' @param quiet Logical. If `TRUE`, suppress download messages where possible.
#' @param mode Download mode passed to `bf_download_cached()`. Defaults to `"wb"`.
#' @param min_bytes Minimum acceptable file size in bytes.
#' @param validate_not_html Logical. If `TRUE`, reject likely HTML error pages.
#'
#' @return Invisibly returns `destfile`.
#'
#' @keywords internal
#' @noRd
bf_download_file <- function(url,
                             destfile,
                             quiet = TRUE,
                             mode = "wb",
                             min_bytes = 1,
                             validate_not_html = TRUE) {
  bf_download_cached(
    url = url,
    dest = destfile,
    force_refresh = TRUE,
    quiet = quiet,
    mode = mode,
    min_bytes = min_bytes,
    validate_not_html = validate_not_html,
    resume = FALSE,
    keep_partial = FALSE
  )

  invisible(destfile)
}

# -----------------------------------------------------------------------------
# ZIP extraction helpers
# -----------------------------------------------------------------------------
# These helpers unpack provider archives into stable cache folders. Reusing an
# existing extraction saves time during repeated pipeline runs, while
# force_refresh allows users/tests to rebuild cached files when needed.

#' Unzip a cached archive into a stable extraction directory
#'
#' Extracts a zip archive into a stable cache directory and reuses an existing
#' extraction when at least one shapefile is already present. If
#' `force_refresh = TRUE`, the extraction directory is deleted and rebuilt.
#'
#' This helper is primarily used for cached vector-overlay downloads distributed
#' as zipped shapefiles.
#'
#' @param zipfile Character. Path to a `.zip` archive.
#' @param exdir Character. Directory where the archive should be extracted.
#' @param force_refresh Logical. If `TRUE`, delete `exdir` and re-extract the
#'   archive even when extracted shapefiles already exist.
#' @param quiet Logical. Retained for interface consistency with download
#'   helpers. Currently only used to suppress non-critical extraction warnings
#'   where possible.
#'
#' @return Character path to the extraction directory, `exdir`.
#'
#' @keywords internal
#' @noRd
bf_unzip_cached <- function(zipfile, exdir, force_refresh = FALSE, quiet = TRUE) {
  stopifnot(is.character(zipfile), length(zipfile) == 1L, nzchar(zipfile))
  stopifnot(is.character(exdir), length(exdir) == 1L, nzchar(exdir))
  stopifnot(file.exists(zipfile))

  if (dir.exists(exdir) && isTRUE(force_refresh)) {
    try(unlink(exdir, recursive = TRUE, force = TRUE), silent = TRUE)
  }

  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  shp <- list.files(
    exdir,
    pattern = "\\.shp$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(shp) > 0L && !isTRUE(force_refresh)) {
    return(exdir)
  }

  if (isTRUE(quiet)) {
    suppressWarnings(utils::unzip(zipfile, exdir = exdir, overwrite = TRUE))
  } else {
    utils::unzip(zipfile, exdir = exdir, overwrite = TRUE)
  }

  exdir
}

# -----------------------------------------------------------------------------
# General archive unpacking
# -----------------------------------------------------------------------------
# This helper handles common compressed formats beyond simple ZIP files and runs
# file-signature checks before extraction so corrupted or HTML cache files fail
# early with useful errors.

#' Unpack a cached archive
#'
#' Expands ZIP, tar.gz, tgz or gzip files after validating their file signature.
#' Existing extraction folders are reused unless `force_refresh = TRUE`.
#'
#' @param path Character path to an archive.
#' @param exdir Character extraction directory.
#' @param force_refresh Logical. If `TRUE`, rebuild an existing extraction folder.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Character path to `exdir`.
#'
#' @keywords internal
#' @noRd
bf_unpack_archive <- function(path, exdir, force_refresh = FALSE, quiet = TRUE) {
  if (!file.exists(path)) {
    stop("Archive does not exist: ", path, call. = FALSE)
  }

  if (isTRUE(force_refresh) && dir.exists(exdir)) {
    try(unlink(exdir, recursive = TRUE, force = TRUE), silent = TRUE)
  }

  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  existing <- list.files(
    exdir,
    recursive = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )

  if (length(existing) > 0L && !isTRUE(force_refresh)) {
    return(exdir)
  }

  # Reject HTML before inspecting archive extensions. Provider error pages often
  # arrive with the requested filename but contain HTML rather than spatial data.
  if (bf_is_probably_html(path)) {
    stop(
      "Downloaded archive appears to be HTML rather than GIS data: ",
      path,
      ". The remote provider may have returned a landing, login or error page.",
      call. = FALSE
    )
  }

  # ZIP/KMZ branch: validate the ZIP signature and list the archive before
  # extraction. Listing first catches empty/corrupt files with a clearer message.
  if (grepl("\\.(zip|kmz)$", path, ignore.case = TRUE)) {
    if (!bf_is_zip_file(path)) {
      stop(
        "Downloaded file is not a valid ZIP archive: ",
        path,
        ". Delete the cache or rerun with force_refresh = TRUE.",
        call. = FALSE
      )
    }

    listed <- tryCatch(
      utils::unzip(path, list = TRUE),
      error = function(e) e
    )

    if (inherits(listed, "error") || !NROW(listed)) {
      stop("ZIP archive could not be listed or is empty: ", path, call. = FALSE)
    }

    utils::unzip(path, exdir = exdir)
  } else if (grepl("\\.(tar\\.gz|tgz)$", path, ignore.case = TRUE)) {
    # tar.gz/tgz branch: require a gzip signature before untarring.
    if (!bf_is_gzip_file(path)) {
      stop("Downloaded file is not a valid gzip/tar archive: ", path, call. = FALSE)
    }

    utils::untar(path, exdir = exdir)
  } else if (grepl("\\.gz$", path, ignore.case = TRUE) &&
             !grepl("\\.tar\\.gz$", path, ignore.case = TRUE)) {
    # Single-file gzip branch: stream-decompress the archive into `exdir` without
    # requiring additional packages.
    if (!bf_is_gzip_file(path)) {
      stop("Downloaded file is not a valid gzip file: ", path, call. = FALSE)
    }

    out <- file.path(
      exdir,
      sub("\\.gz$", "", basename(path), ignore.case = TRUE)
    )

    con_in <- gzfile(path, "rb")
    con_out <- file(out, "wb")

    on.exit(try(close(con_in), silent = TRUE), add = TRUE)
    on.exit(try(close(con_out), silent = TRUE), add = TRUE)

    repeat {
      buf <- readBin(con_in, what = raw(), n = 1024 * 1024)
      if (!length(buf)) break
      writeBin(buf, con_out)
    }
  } else {
    stop("Unsupported archive type: ", path, call. = FALSE)
  }

  if (!isTRUE(quiet)) {
    message("Archive unpacked to: ", exdir)
  }

  exdir
}

# -----------------------------------------------------------------------------
# Nested archive extraction
# -----------------------------------------------------------------------------
# Some data providers distribute one top-level archive containing further ZIP
# files. This helper expands any nested ZIPs into same-name folders so later
# vector-file discovery can search a normal directory tree.

#' Recursively unpack nested ZIP archives
#'
#' Searches a directory for ZIP files and extracts each one into a same-name
#' subdirectory. This is useful for provider downloads that contain a top-level
#' archive with one or more nested data archives.
#'
#' @param root_dir Directory to search recursively.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Invisibly returns `root_dir`.
#'
#' @keywords internal
#' @noRd
bf_unzip_nested <- function(root_dir, quiet = TRUE) {
  if (!dir.exists(root_dir)) {
    stop("Directory does not exist: ", root_dir, call. = FALSE)
  }

  inner_zips <- list.files(
    root_dir,
    pattern = "\\.zip$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (!length(inner_zips)) {
    return(invisible(root_dir))
  }

  for (z in inner_zips) {
    target <- file.path(
      dirname(z),
      tools::file_path_sans_ext(basename(z))
    )

    dir.create(target, showWarnings = FALSE, recursive = TRUE)
    utils::unzip(z, exdir = target)

    if (!isTRUE(quiet)) {
      message("Nested ZIP unpacked: ", z)
    }
  }

  invisible(root_dir)
}
