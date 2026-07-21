################################################################################
# utils_griis.R
# ------------------------------------------------------------------------------
# biofetchR: GRIIS download, standardisation and invasive-status utilities
# ------------------------------------------------------------------------------
#
# PURPOSE
#   This file contains the retained Global Register of Introduced and Invasive
#   Species (GRIIS) helpers used by biofetchR's pre-GBIF origin-evidence gates.
#
#   The processing workflow is:
#     1. Download or reuse the cached country-level GRIIS compendium.
#     2. Read the country-level table from a CSV/TSV-like source.
#     3. Standardise taxon names, country identifiers and invasive-status fields.
#     4. Build lookup tables either by species x country or by species only.
#     5. Attach GRIIS evidence to terrestrial/freshwater and marine inputs.
#     6. Optionally filter retained rows to GRIIS-listed invasive records.
#
# DESIGN NOTES
#   - Terrestrial/freshwater workflows normally use species x recipient-country
#     matching because the GRIIS country compendium is country-specific.
#   - Marine/global workflows can use species-only matching with
#     require_country = FALSE. In that mode, country parsing is deliberately
#     bypassed so species-only marine inputs do not trigger zero-length ISO
#     assignment errors.
#   - The downloader writes to a temporary file before promoting to the final
#     cache path. This protects first-run and Windows workflows from partial
#     cache files and brittle file.rename() behaviour.
#   - The helpers reject the GRIIS metadata/DwC-A archive because it does not
#     contain the species x country fields required by biofetchR's GRIIS gate.
#
# OUTPUT SCOPE
#   These helpers create cached GRIIS source files, standardised GRIIS tables,
#   species/country lookup tables and input tables with appended GRIIS evidence.
#   They do not submit GBIF downloads, write occurrence outputs or create plots.
#
# DATA SOURCE AND ATTRIBUTION
#   GRIIS is an external data product. biofetchR can download, cache and
#   standardise the country-level compendium for analysis, but it does not own or
#   redistribute the underlying source. Users should cite the GRIIS source, its
#   version/record used in their workflow, and any associated data publication or
#   provider guidance in downstream analyses.
#
################################################################################

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

#' Normalise column names for flexible GRIIS schema matching
#'
#' Converts column names to lower-case, underscore-separated strings. The GRIIS
#' country compendium and related tables may use slightly different column naming
#' conventions, so helpers compare cleaned names rather than relying on exact raw
#' column names.
#'
#' @param x Character vector of column names.
#'
#' @return Character vector of cleaned column names.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_clean_colname <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

#' Pick the first matching GRIIS column from candidate names
#'
#' Searches a data frame for the first column whose cleaned name matches one of the
#' candidate names. This allows the standardisation step to support small schema
#' variations while still failing clearly when required fields are absent.
#'
#' @param dat Data frame to search.
#' @param candidates Character vector of possible column names.
#'
#' @return Name of the first matching column, or `NA_character_` if no candidate is
#' present.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_pick_col <- function(dat, candidates) {
  if (is.null(dat) || !is.data.frame(dat)) return(NA_character_)

  nms <- names(dat)
  nms_clean <- .bf_griis_clean_colname(nms)
  cand_clean <- .bf_griis_clean_colname(candidates)

  hit <- which(nms_clean %in% cand_clean)

  if (length(hit)) {
    return(nms[hit[[1]]])
  }

  NA_character_
}

#' Convert GRIIS invasive-status labels to logical values
#'
#' Maps common yes/no and invasive/non-invasive labels to `TRUE`, `FALSE` or `NA`.
#' The helper is intentionally conservative: unrecognised values remain `NA` rather
#' than being forced into a binary class.
#'
#' @param x Vector of invasive-status labels.
#'
#' @return Logical vector with `TRUE` for invasive/yes labels, `FALSE` for
#' non-invasive/no labels and `NA` for ambiguous or missing values.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_yes_no <- function(x) {
  x <- tolower(trimws(as.character(x)))

  out <- rep(NA, length(x))

  true_vals <- c(
    "yes", "y", "true", "t", "1",
    "invasive",
    "invasive alien",
    "alien invasive",
    "invasive alien species",
    "ias"
  )

  false_vals <- c(
    "no", "n", "false", "f", "0",
    "not invasive",
    "non-invasive",
    "introduced",
    "alien",
    "naturalized",
    "naturalised",
    "established",
    "native"
  )

  out[x %in% true_vals] <- TRUE
  out[x %in% false_vals] <- FALSE

  out
}

#' Recycle optional vectors to a target length
#'
#' Used when country, ISO2 and ISO3 inputs may be supplied as scalars, vectors or
#' missing values. Length-zero inputs become an `NA_character_` vector of the target
#' length; scalars are recycled; incompatible vector lengths raise an error.
#'
#' @param x Input vector.
#' @param n Target length.
#' @param arg Name of the argument, used in error messages.
#'
#' @return Character vector of length `n`.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_recycle_to <- function(x, n, arg = "input") {
  if (is.null(x) || length(x) == 0L) {
    return(rep(NA_character_, n))
  }

  x <- as.character(x)

  if (length(x) == n) {
    return(x)
  }

  if (length(x) == 1L) {
    return(rep(x, n))
  }

  stop(
    "`", arg, "` has length ", length(x),
    ", but expected length 1 or ", n, ".",
    call. = FALSE
  )
}

#' Create a normalised species-name key for GRIIS joins
#'
#' Builds the species-name key used to match user input against GRIIS accepted and
#' scientific names. When `bf_clean_taxon_names()` is available, it is applied first
#' so GRIIS matching follows the same taxon-name cleaning used elsewhere in
#' biofetchR.
#'
#' @param x Character vector of taxon names.
#'
#' @return Lower-case, whitespace-normalised character vector suitable for joins.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_name_key <- function(x) {
  x <- bf_clean_text(x)

  if (exists("bf_clean_taxon_names", mode = "function")) {
    cleaned <- suppressWarnings(bf_clean_taxon_names(x))
    x <- cleaned
  }

  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x[x == "" | is.na(x)] <- NA_character_

  x
}

#' Resolve country evidence to ISO3 codes
#'
#' Combines optional country names, ISO2 codes and ISO3 codes into a single ISO3
#' vector. Existing ISO3 values take priority, followed by ISO2 conversion, then
#' country-name conversion. The helper includes fallbacks for common aliases such
#' as `UK`, `NA` for Namibia and `XK` for Kosovo.
#'
#' @param country Optional country-name vector.
#' @param iso2c Optional ISO2 country-code vector.
#' @param iso3c Optional ISO3 country-code vector.
#'
#' @return Character vector of ISO3 codes. Returns `character(0)` when all inputs
#' are absent.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_country_to_iso3 <- function(country = NULL,
                                      iso2c = NULL,
                                      iso3c = NULL) {
  n <- max(
    length(country %||% character(0)),
    length(iso2c %||% character(0)),
    length(iso3c %||% character(0)),
    0L
  )

  if (n == 0L) {
    return(character(0))
  }

  country <- .bf_griis_recycle_to(country, n, "country")
  iso2c   <- .bf_griis_recycle_to(iso2c, n, "iso2c")
  iso3c   <- .bf_griis_recycle_to(iso3c, n, "iso3c")

  out <- rep(NA_character_, n)

  # ---------------------------------------------------------------------------
  # ISO3 values, if supplied, have highest priority.
  # ---------------------------------------------------------------------------

  x3 <- toupper(trimws(as.character(iso3c)))
  x3[x3 == ""] <- NA_character_

  fill3 <- !is.na(x3)
  out[fill3] <- x3[fill3]

  # ---------------------------------------------------------------------------
  # ISO2 values.
  # ---------------------------------------------------------------------------

  x2 <- toupper(trimws(as.character(iso2c)))
  x2[x2 == ""] <- NA_character_
  x2[x2 == "UK"] <- "GB"

  y2 <- rep(NA_character_, n)

  if (requireNamespace("countrycode", quietly = TRUE)) {
    y2 <- suppressWarnings(countrycode::countrycode(
      x2,
      origin = "iso2c",
      destination = "iso3c",
      custom_match = c(
        "GB" = "GBR",
        "UK" = "GBR",
        "NA" = "NAM",
        "XK" = "XKX"
      ),
      warn = FALSE
    ))
  } else {
    y2[x2 == "GB"] <- "GBR"
    y2[x2 == "UK"] <- "GBR"
    y2[x2 == "NA"] <- "NAM"
    y2[x2 == "US"] <- "USA"
    y2[x2 == "NZ"] <- "NZL"
    y2[x2 == "ZA"] <- "ZAF"
    y2[x2 == "AU"] <- "AUS"
    y2[x2 == "CA"] <- "CAN"
  }

  fill2 <- is.na(out) & !is.na(y2)
  out[fill2] <- y2[fill2]

  # ---------------------------------------------------------------------------
  # Country names.
  # ---------------------------------------------------------------------------

  xc <- bf_clean_text(country)

  yc <- rep(NA_character_, n)

  if (requireNamespace("countrycode", quietly = TRUE)) {
    yc <- suppressWarnings(countrycode::countrycode(
      xc,
      origin = "country.name",
      destination = "iso3c",
      custom_match = c(
        "Namibia" = "NAM",
        "United Kingdom" = "GBR",
        "Great Britain" = "GBR",
        "Britain" = "GBR",
        "United States" = "USA",
        "United States of America" = "USA",
        "New Zealand" = "NZL",
        "Kosovo" = "XKX"
      ),
      warn = FALSE
    ))
  } else {
    yc[grepl("^namibia$", xc, ignore.case = TRUE)] <- "NAM"
    yc[grepl("^united kingdom$|^great britain$|^britain$", xc, ignore.case = TRUE)] <- "GBR"
    yc[grepl("^united states$|^united states of america$", xc, ignore.case = TRUE)] <- "USA"
    yc[grepl("^new zealand$", xc, ignore.case = TRUE)] <- "NZL"
  }

  fillc <- is.na(out) & !is.na(yc)
  out[fillc] <- yc[fillc]

  out
}

#' Ensure a GRIIS table has biofetchR standard columns
#'
#' Accepts either a raw GRIIS-like table or a table that has already been
#' standardised by [bf_standardise_griis()]. If `griis` is `NULL`, the default GRIIS
#' reader is used.
#'
#' @param griis Raw, standardised or `NULL` GRIIS table.
#'
#' @return A standardised GRIIS data frame with the columns required by lookup and
#' join helpers.
#'
#'
#' @md
#' @keywords internal
#' @noRd
.bf_griis_prepare_standardised <- function(griis) {
  if (is.null(griis)) {
    griis <- bf_read_griis()
  }

  if (!is.data.frame(griis)) {
    stop("`griis` must be a data frame.", call. = FALSE)
  }

  required <- c(
    "griis_scientific_name",
    "griis_accepted_name",
    "griis_name_key",
    "griis_scientific_name_key",
    "griis_iso3c",
    "griis_is_invasive"
  )

  if (!all(required %in% names(griis))) {
    griis <- bf_standardise_griis(griis)
  }

  griis
}

# Internal default GRIIS country-compendium URL.
.bf_griis_default_source_url <- function() {
  "https://zenodo.org/records/6348164/files/GRIIS%20-%20Country%20Compendium%20V1_0.csv?download=1"
}

#' Download the GRIIS country compendium
#'
#' Downloads the country-level Global Register of Introduced and Invasive Species
#' (GRIIS) compendium and stores it in a local cache. The downloader is designed to
#' be robust on Windows: files are first downloaded to a `.tmp` path and then
#' promoted to the final cache path by rename or copy.
#'
#' @details
#' biofetchR requires the country-level GRIIS table because the terrestrial and
#' freshwater workflows use species x recipient-country evidence. The helper
#' therefore rejects the `griis_meta_recoded` metadata/DwC-A archive, which lacks
#' the required country-level fields used by the GRIIS gate.
#'
#' If a previous run left a valid temporary file in the cache, the helper attempts
#' to finalise that file before downloading again. This avoids unnecessary repeated
#' downloads when a prior session succeeded but file finalisation failed.
#'
#' @section Data source and attribution:
#' This function downloads the external GRIIS country compendium into a local
#' user cache. biofetchR does not redistribute or claim ownership of the source
#' data. Users should cite the GRIIS data source, version/record and associated
#' provider guidance used in their analysis.
#'
#' @param cache_dir Directory used for the cached GRIIS file. Must be supplied
#'   explicitly. In examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, re-download even if a cached file is
#' already present.
#' @param source_url Character URL for the country-level GRIIS compendium. The
#' default points to the Zenodo-hosted country compendium CSV used by biofetchR.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return A character vector of length one giving the normalised path to the
#'   cached GRIIS country compendium file, using forward slashes. The path points
#'   to either an existing valid cached file or a newly downloaded file. The
#'   function is also called for the side effect of downloading, validating and
#'   caching the external GRIIS source file.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   path <- bf_download_griis(
#'     cache_dir = file.path(tempdir(), "biofetchR_griis"),
#'     quiet = FALSE
#'   )
#'
#'   file.exists(path)
#' }
#' }
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_download_griis <- function(
    cache_dir = NULL,
    force_refresh = FALSE,
    source_url = .bf_griis_default_source_url(),
    quiet = FALSE
) {
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

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # The Zenodo 14758504 `griis_meta_recoded.tar.gz` object is a DwC-A /
  # metadata-style GBIF export. It does not provide the country-level fields
  # required by biofetchR's GRIIS gate (`species`, `countryCode_alpha2`,
  # `countryCode_alpha3`, `isInvasive`). Do not silently accept it, because it
  # produces a species-only lookup with empty ISO and invasive fields.
  if (
    !is.null(source_url) &&
    length(source_url) == 1L &&
    !is.na(source_url) &&
    grepl("griis_meta_recoded|14758504", source_url, ignore.case = TRUE)
  ) {
    stop(
      "The supplied GRIIS source appears to be the metadata/DwC-A archive ",
      "(`griis_meta_recoded`). biofetchR requires the country-level GRIIS ",
      "Country Compendium CSV with species, countryCode_alpha2, ",
      "countryCode_alpha3 and isInvasive columns. Use the default source_url ",
      "or provide a local/download URL for that country-level table.",
      call. = FALSE
    )
  }

  source_no_query <- sub("\\?.*$", "", source_url)
  dest_name <- basename(utils::URLdecode(source_no_query))

  if (!nzchar(dest_name) || dest_name %in% c(".", "/")) {
    dest_name <- "GRIIS_Country_Compendium.csv"
  }

  # Keep a stable, filesystem-safe filename while preserving the extension.
  dest_name <- gsub("[^A-Za-z0-9._-]+", "_", dest_name)
  dest_path <- file.path(cache_dir, dest_name)
  tmp_path  <- paste0(dest_path, ".tmp")

  # If a valid cached final file exists, use it unless explicitly refreshing.
  if (file.exists(dest_path) && !isTRUE(force_refresh)) {
    if (file.info(dest_path)$size >= 1000000) {
      .bf_msg("GRIIS country compendium already cached: ", dest_path, quiet = quiet)
      return(normalizePath(dest_path, winslash = "/", mustWork = TRUE))
    } else {
      unlink(dest_path, force = TRUE)
    }
  }

  # If a previous run left a valid .tmp file, promote it before downloading again.
  if (file.exists(tmp_path) && !isTRUE(force_refresh)) {
    tmp_size <- suppressWarnings(file.info(tmp_path)$size)

    if (!is.na(tmp_size) && tmp_size >= 1000000) {
      .bf_msg("Found existing GRIIS temporary file; finalising cache...", quiet = quiet)

      ok_copy <- file.copy(tmp_path, dest_path, overwrite = TRUE)

      if (isTRUE(ok_copy) && file.exists(dest_path) && file.info(dest_path)$size >= 1000000) {
        unlink(tmp_path, force = TRUE)
        .bf_msg("GRIIS country compendium cached: ", dest_path, quiet = quiet)
        return(normalizePath(dest_path, winslash = "/", mustWork = TRUE))
      }
    }
  }

  # Clean stale/incomplete temporary file before a fresh download.
  if (file.exists(tmp_path)) {
    unlink(tmp_path, force = TRUE)
  }

  if (file.exists(dest_path) && isTRUE(force_refresh)) {
    unlink(dest_path, force = TRUE)
  }

  .bf_msg("Downloading GRIIS country compendium...", quiet = quiet)
  .bf_msg("Source: ", source_url, quiet = quiet)

  tryCatch(
    {
      bf_download_cached(
        url = source_url,
        dest = tmp_path,
        force_refresh = TRUE,
        quiet = quiet,
        mode = "wb",
        min_bytes = 1000000,
        validate_not_html = TRUE,
        resume = TRUE,
        keep_partial = TRUE
      )

      if (!file.exists(tmp_path) || file.info(tmp_path)$size < 1000000) {
        stop("Downloaded GRIIS file is missing or unexpectedly small.", call. = FALSE)
      }

      if (file.exists(dest_path)) {
        unlink(dest_path, force = TRUE)
      }

      # First try rename. If Windows blocks it, fall back to copy + delete.
      ok_rename <- suppressWarnings(file.rename(tmp_path, dest_path))

      if (!isTRUE(ok_rename)) {
        ok_copy <- file.copy(tmp_path, dest_path, overwrite = TRUE)

        if (isTRUE(ok_copy)) {
          unlink(tmp_path, force = TRUE)
        } else {
          stop(
            "Could not move or copy temporary GRIIS file into cache: ",
            tmp_path,
            " -> ",
            dest_path,
            call. = FALSE
          )
        }
      }

      if (!file.exists(dest_path) || file.info(dest_path)$size < 1000000) {
        stop("Final cached GRIIS file is missing or unexpectedly small.", call. = FALSE)
      }
    },
    error = function(e) {
      # Do NOT delete tmp_path here. Keeping it allows recovery if the download
      # succeeded but finalisation failed on Windows.
      stop(
        "Failed to download/finalise GRIIS country compendium: ",
        conditionMessage(e),
        "\nTemporary file retained, if present: ",
        tmp_path,
        call. = FALSE
      )
    }
  )

  .bf_msg("GRIIS country compendium cached: ", dest_path, quiet = quiet)

  normalizePath(dest_path, winslash = "/", mustWork = TRUE)
}

#' Unpack a cached GRIIS archive
#'
#' Unpacks a downloaded GRIIS archive into a cache directory. This helper is kept
#' for archive-style GRIIS sources, although the default biofetchR source is now a
#' country-level CSV that can be read directly without unpacking.
#'
#' @param archive_path Path returned by [bf_download_griis()] or another downloaded
#' GRIIS archive file.
#' @param exdir Directory to unpack into. Must be supplied explicitly. In
#'   examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, remove and recreate `exdir` before
#' unpacking.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return A character vector of length one giving the normalised path to the
#'   directory containing the unpacked GRIIS archive, using forward slashes. If
#'   the archive has already been unpacked and `force_refresh = FALSE`, the
#'   existing directory path is returned. The function is also called for the
#'   side effect of creating or refreshing the unpacked cache directory.
#'
#' @examples
#' archive_root <- file.path(tempdir(), "biofetchR_griis_unpack_example")
#' unlink(archive_root, recursive = TRUE, force = TRUE)
#' dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
#'
#' archive_src <- file.path(archive_root, "src")
#' dir.create(archive_src, recursive = TRUE, showWarnings = FALSE)
#'
#' griis_file <- file.path(archive_src, "griis_country_compendium.csv")
#'
#' write.csv(
#'   data.frame(species = "Example species"),
#'   griis_file,
#'   row.names = FALSE
#' )
#'
#' archive_path <- file.path(archive_root, "griis_example.tar")
#'
#' utils::tar(
#'   tarfile = archive_path,
#'   files = griis_file,
#'   tar = "internal"
#' )
#'
#' unpacked <- bf_unpack_griis(
#'   archive_path = archive_path,
#'   exdir = file.path(archive_root, "unpacked"),
#'   quiet = TRUE
#' )
#'
#' dir.exists(unpacked)
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_unpack_griis <- function(
    archive_path,
    exdir = NULL,
    force_refresh = FALSE,
    quiet = FALSE
) {
  if (!file.exists(archive_path)) {
    stop("GRIIS archive does not exist: ", archive_path, call. = FALSE)
  }

  if (is.null(exdir) || length(exdir) == 0L ||
      !nzchar(trimws(as.character(exdir[[1L]])))) {
    stop(
      "`exdir` must be supplied explicitly. In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  exdir <- normalizePath(
    as.character(exdir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (dir.exists(exdir) && !isTRUE(force_refresh)) {
    .bf_msg("GRIIS archive already unpacked: ", exdir, quiet = quiet)
    return(normalizePath(exdir, winslash = "/", mustWork = TRUE))
  }

  if (dir.exists(exdir) && isTRUE(force_refresh)) {
    unlink(exdir, recursive = TRUE, force = TRUE)
  }

  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  .bf_msg("Unpacking GRIIS archive...", quiet = quiet)

  tryCatch(
    utils::untar(archive_path, exdir = exdir),
    error = function(e) {
      stop(
        "Failed to unpack GRIIS archive: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  normalizePath(exdir, winslash = "/", mustWork = TRUE)
}

#' Find the most likely GRIIS table in an unpacked archive
#'
#' Searches an unpacked GRIIS directory for CSV, TSV, TXT or TAB files and selects
#' the most likely data table using filename and file-size heuristics. Metadata,
#' readme and citation files are down-weighted.
#'
#' @param exdir Directory returned by [bf_unpack_griis()].
#'
#' @return A character vector of length one giving the path to the selected
#'   GRIIS table file within `exdir`. The selected file is chosen from CSV, TSV,
#'   TXT or TAB files using filename and file-size heuristics that favour likely
#'   species or compendium data tables and down-weight metadata, readme and
#'   citation files.
#'
#' @examples
#' exdir <- tempfile("griis_example_")
#' dir.create(exdir)
#'
#' write.csv(
#'   data.frame(species = "Example species"),
#'   file.path(exdir, "griis_country_compendium.csv"),
#'   row.names = FALSE
#' )
#'
#' table_path <- bf_find_griis_table(exdir)
#' basename(table_path)
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_find_griis_table <- function(exdir) {
  if (!dir.exists(exdir)) {
    stop("GRIIS unpack directory does not exist: ", exdir, call. = FALSE)
  }

  files <- list.files(
    exdir,
    pattern = "\\.(csv|tsv|txt|tab)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (!length(files)) {
    stop(
      "No CSV/TSV/TXT/TAB files were found inside the GRIIS archive: ",
      exdir,
      call. = FALSE
    )
  }

  base <- basename(files)
  score <- rep(0, length(files))

  score <- score + ifelse(grepl("griis", base, ignore.case = TRUE), 5, 0)
  score <- score + ifelse(grepl("taxon|species|record|compendium|data", base, ignore.case = TRUE), 4, 0)
  score <- score - ifelse(grepl("metadata|meta|eml|readme|citation", base, ignore.case = TRUE), 6, 0)

  sizes <- file.info(files)$size
  score <- score + rank(sizes, ties.method = "first") / length(sizes)

  files[order(score, decreasing = TRUE)][[1]]
}

#' Read and standardise the GRIIS country compendium
#'
#' Downloads or reuses the cached GRIIS country compendium, reads the selected table
#' and converts it to biofetchR's standard GRIIS schema. This is the main entry
#' point used by the pipeline origin-evidence gates when a pre-read GRIIS table is
#' not supplied.
#'
#' @section Data source and attribution:
#' This helper reads the external GRIIS country compendium after download/caching
#' and standardises it for biofetchR joins. Users remain responsible for citing
#' the GRIIS source, version/record and any associated data publication in
#' downstream outputs.
#'
#' @param cache_dir Directory used for the cached GRIIS file.
#' @param force_refresh Logical; if `TRUE`, re-download and re-read the GRIIS
#' source.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return A tibble returned by [bf_standardise_griis()] containing the
#'   standardised GRIIS country-compendium fields used by biofetchR. The output
#'   includes stable species-name keys, country identifiers, ISO2 and ISO3 codes,
#'   invasive-status evidence and optional taxonomic/status metadata. Each row
#'   represents a standardised GRIIS species-country evidence record.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   griis <- bf_read_griis(
#'     cache_dir = file.path(tempdir(), "biofetchR_griis"),
#'     quiet = FALSE
#'   )
#'
#'   names(griis)
#' }
#' }
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_read_griis <- function(
    cache_dir = NULL,
    force_refresh = FALSE,
    quiet = FALSE
) {
  bf_require_packages(c("readr", "dplyr", "tibble"), context = "GRIIS helper")

  source_path <- bf_download_griis(
    cache_dir = cache_dir,
    force_refresh = force_refresh,
    quiet = quiet
  )

  ext <- tolower(tools::file_ext(source_path))

  if (ext %in% c("csv", "tsv", "txt", "tab")) {
    table_path <- source_path
  } else {
    exdir <- bf_unpack_griis(
      archive_path = source_path,
      exdir = file.path(dirname(source_path), "griis_unpacked"),
      force_refresh = force_refresh,
      quiet = quiet
    )

    table_path <- bf_find_griis_table(exdir)
    ext <- tolower(tools::file_ext(table_path))
  }

  .bf_msg("Reading GRIIS country-level table: ", table_path, quiet = quiet)

  raw <- tryCatch(
    {
      if (ext == "csv") {
        readr::read_csv(table_path, show_col_types = FALSE, progress = !quiet)
      } else {
        readr::read_tsv(table_path, show_col_types = FALSE, progress = !quiet)
      }
    },
    error = function(e) {
      stop(
        "Failed to read GRIIS table: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  bf_standardise_griis(raw)
}

# -----------------------------------------------------------------------------
# Standardisation
# -----------------------------------------------------------------------------

#' Standardise a raw GRIIS table
#'
#' Converts a raw country-level GRIIS table into the stable column names used by
#' biofetchR. The function identifies accepted/scientific name fields, country/ISO
#' fields, invasive-status fields and optional taxonomic metadata using flexible
#' column-name matching.
#'
#' @details
#' The returned table includes both accepted-name and scientific-name join keys.
#' Downstream lookup construction uses both keys so that records can match either
#' the accepted GRIIS name or a reported scientific-name field.
#'
#' This helper deliberately checks that the table resembles the country-level GRIIS
#' compendium. If required country-level fields are missing, it stops rather than
#' silently producing a species-only table with empty country or invasive-status
#' evidence.
#'
#' @section Interpretation:
#' Standardised columns are evidence fields derived from GRIIS records. They
#' should be interpreted as GRIIS-listed status evidence, not as proof that a
#' species is absent from countries where no match is found.
#'
#' @param griis_raw Raw data frame read from a GRIIS source.
#'
#' @return A tibble with stable `griis_*` columns used by biofetchR lookup and
#'   origin-evidence workflows. The output includes `griis_scientific_name`,
#'   `griis_accepted_name`, `griis_name_key`,
#'   `griis_scientific_name_key`, `griis_country`, `griis_iso2c`,
#'   `griis_iso3c` and `griis_is_invasive`, plus optional occurrence,
#'   establishment and taxonomic metadata where available. Each row represents a
#'   standardised GRIIS species-country evidence record. The invasive-status
#'   fields should be interpreted as GRIIS evidence, not as proof of absence
#'   where no match exists.
#'
#' @examples
#' raw <- data.frame(
#'   species = "Example species",
#'   scientificName = "Example species",
#'   countryCode_alpha2 = "GB",
#'   countryCode_alpha3 = "GBR",
#'   isInvasive = "yes",
#'   stringsAsFactors = FALSE
#' )
#' std <- bf_standardise_griis(raw)
#' names(std)
#'
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_standardise_griis <- function(griis_raw) {
  bf_require_packages(c("dplyr", "tibble"), context = "GRIIS helper")

  if (!is.data.frame(griis_raw)) {
    stop("`griis_raw` must be a data frame.", call. = FALSE)
  }

  dat <- tibble::as_tibble(griis_raw)

  # In the country-level GRIIS compendium, `species` is the matching name.
  # `scientificName` may contain a more verbatim/reported taxon string and is
  # retained as supporting information, not used in preference to `species`.
  accepted_col <- .bf_griis_pick_col(
    dat,
    c(
      "species",
      "acceptedNameUsage",
      "accepted_name_usage",
      "acceptedName",
      "accepted_name",
      "taxonConcept",
      "taxon_concept",
      "canonicalName",
      "canonical_name"
    )
  )

  scientific_col <- .bf_griis_pick_col(
    dat,
    c(
      "scientificName",
      "scientific_name",
      "reportedTaxon",
      "reported_taxon",
      "taxonName",
      "taxon_name",
      "verbatimScientificName"
    )
  )

  country_col <- .bf_griis_pick_col(
    dat,
    c(
      "country",
      "location",
      "locationName",
      "location_name",
      "area",
      "territory"
    )
  )

  iso2_col <- .bf_griis_pick_col(
    dat,
    c(
      "countryCode_alpha2",
      "country_code_alpha2",
      "countryCode",
      "country_code",
      "iso2",
      "iso2c",
      "country_iso2"
    )
  )

  iso3_col <- .bf_griis_pick_col(
    dat,
    c(
      "countryCode_alpha3",
      "country_code_alpha3",
      "iso3",
      "iso3c",
      "country_iso3"
    )
  )

  invasive_col <- .bf_griis_pick_col(
    dat,
    c(
      "isInvasive",
      "is_invasive",
      "invasive",
      "invasiveAlien",
      "invasive_alien",
      "impact",
      "hasImpact",
      "has_impact"
    )
  )

  occurrence_col <- .bf_griis_pick_col(
    dat,
    c(
      "occurrenceStatus",
      "occurrence_status",
      "presence",
      "status"
    )
  )

  establishment_col <- .bf_griis_pick_col(
    dat,
    c(
      "establishmentMeans",
      "establishment_means",
      "degreeOfEstablishment",
      "degree_of_establishment",
      "origin"
    )
  )

  rank_col <- .bf_griis_pick_col(
    dat,
    c(
      "taxonRank",
      "taxon_rank",
      "rank"
    )
  )

  kingdom_col <- .bf_griis_pick_col(dat, c("kingdom"))
  phylum_col  <- .bf_griis_pick_col(dat, c("phylum"))
  class_col   <- .bf_griis_pick_col(dat, c("class"))
  order_col   <- .bf_griis_pick_col(dat, c("order"))
  family_col  <- .bf_griis_pick_col(dat, c("family"))
  genus_col   <- .bf_griis_pick_col(dat, c("genus"))

  if (is.na(accepted_col) && is.na(scientific_col)) {
    stop(
      "Could not identify a scientific-name or species column in the GRIIS table. ",
      "For the country-level GRIIS compendium the required matching column is ",
      "usually `species`.",
      call. = FALSE
    )
  }

  required_country_cols <- c("countryCode_alpha2", "countryCode_alpha3", "isInvasive")
  missing_required_country_cols <- required_country_cols[
    !.bf_griis_clean_colname(required_country_cols) %in% .bf_griis_clean_colname(names(dat))
  ]

  if (length(missing_required_country_cols)) {
    stop(
      "The selected GRIIS table does not look like the country-level ",
      "compendium. Missing expected column(s): ",
      paste(missing_required_country_cols, collapse = ", "),
      ". Check that bf_download_griis() is using the Country Compendium CSV, ",
      "not the `griis_meta_recoded` metadata/DwC-A archive.",
      call. = FALSE
    )
  }

  get_col <- function(col, default = NA_character_) {
    if (is.na(col) || !col %in% names(dat)) {
      rep(default, nrow(dat))
    } else {
      dat[[col]]
    }
  }

  scientific_name <- bf_clean_text(get_col(scientific_col))
  accepted_name   <- bf_clean_text(get_col(accepted_col))

  accepted_name[is.na(accepted_name)] <- scientific_name[is.na(accepted_name)]
  scientific_name[is.na(scientific_name)] <- accepted_name[is.na(scientific_name)]

  country <- bf_clean_text(get_col(country_col))
  iso2c   <- toupper(bf_clean_text(get_col(iso2_col)))
  iso3c   <- toupper(bf_clean_text(get_col(iso3_col)))

  iso3_from_any <- .bf_griis_country_to_iso3(
    country = country,
    iso2c = iso2c,
    iso3c = iso3c
  )

  is_invasive <- .bf_griis_yes_no(get_col(invasive_col))

  out <- tibble::tibble(
    griis_scientific_name = scientific_name,
    griis_accepted_name = accepted_name,
    griis_name_key = .bf_griis_name_key(accepted_name),
    griis_scientific_name_key = .bf_griis_name_key(scientific_name),
    griis_country = country,
    griis_iso2c = iso2c,
    griis_iso3c = iso3_from_any,
    griis_is_invasive = is_invasive,
    griis_occurrence_status = bf_clean_text(get_col(occurrence_col)),
    griis_establishment_means = bf_clean_text(get_col(establishment_col)),
    griis_taxon_rank = toupper(bf_clean_text(get_col(rank_col))),
    griis_kingdom = bf_clean_text(get_col(kingdom_col)),
    griis_phylum = bf_clean_text(get_col(phylum_col)),
    griis_class = bf_clean_text(get_col(class_col)),
    griis_order = bf_clean_text(get_col(order_col)),
    griis_family = bf_clean_text(get_col(family_col)),
    griis_genus = bf_clean_text(get_col(genus_col))
  )

  out |>
    dplyr::filter(
      !is.na(.data$griis_name_key) |
        !is.na(.data$griis_scientific_name_key)
    )
}

# -----------------------------------------------------------------------------
# GRIIS lookup and join
# -----------------------------------------------------------------------------

#' Build a compact GRIIS lookup table
#'
#' Creates a deduplicated lookup table from a raw or standardised GRIIS table. The
#' lookup can be built either at species x country resolution or at species-only
#' resolution.
#'
#' @details
#' Use `require_country = TRUE` for terrestrial and freshwater species-country
#' workflows, where invasive status should be matched against the recipient country.
#' Use `require_country = FALSE` for marine or global workflows that only need
#' species-level GRIIS evidence before GBIF download submission.
#'
#' The lookup combines accepted-name and scientific-name keys, then summarises all
#' matching GRIIS rows into a compact evidence table.
#'
#' @section Interpretation:
#' The lookup reports whether a supplied name/country combination matches GRIIS
#' evidence and whether any matched GRIIS row is flagged invasive. Non-matches
#' should be treated as `not_listed_or_no_match`, not as confirmed absence of
#' introduction or impact.
#'
#' @param griis Raw or standardised GRIIS table.
#' @param require_country Logical. If `TRUE`, build a species x ISO3 lookup. If
#' `FALSE`, build a species-only lookup.
#'
#' @return A deduplicated tibble used as a compact GRIIS lookup table. When
#'   `require_country = TRUE`, each row represents a species-name key by ISO3
#'   country-code combination. When `require_country = FALSE`, each row
#'   represents species-level GRIIS evidence without country matching. The output
#'   includes `griis_listed`, `griis_invasive`,
#'   `griis_any_invasive_field_present`, `griis_n_matches` and collapsed
#'   supporting fields such as matched names, countries, establishment means,
#'   occurrence status and taxon rank. These fields summarise matching GRIIS
#'   evidence for downstream joins and filters.
#'
#' @examples
#' raw <- data.frame(
#'   species = "Example species",
#'   scientificName = "Example species",
#'   countryCode_alpha2 = "GB",
#'   countryCode_alpha3 = "GBR",
#'   isInvasive = "yes",
#'   stringsAsFactors = FALSE
#' )
#' lookup <- bf_griis_lookup(raw, require_country = TRUE)
#' lookup
#'
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_griis_lookup <- function(griis,
                            require_country = TRUE) {
  bf_require_packages(c("dplyr", "tibble"), context = "GRIIS helper")

  griis <- .bf_griis_prepare_standardised(griis)

  required <- c(
    "griis_name_key",
    "griis_scientific_name_key",
    "griis_iso3c",
    "griis_is_invasive"
  )

  missing <- setdiff(required, names(griis))

  if (length(missing)) {
    stop(
      "`griis` is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  optional_cols <- c(
    "griis_country",
    "griis_accepted_name",
    "griis_scientific_name",
    "griis_establishment_means",
    "griis_occurrence_status",
    "griis_taxon_rank"
  )

  for (nm in optional_cols) {
    if (!nm %in% names(griis)) griis[[nm]] <- NA_character_
  }

  by_accepted <- griis |>
    dplyr::transmute(
      griis_join_name_key = .data$griis_name_key,
      griis_iso3c = .data$griis_iso3c,
      griis_country = .data$griis_country,
      griis_matched_name = .data$griis_accepted_name,
      griis_is_invasive = .data$griis_is_invasive,
      griis_establishment_means = .data$griis_establishment_means,
      griis_occurrence_status = .data$griis_occurrence_status,
      griis_taxon_rank = .data$griis_taxon_rank
    )

  by_scientific <- griis |>
    dplyr::transmute(
      griis_join_name_key = .data$griis_scientific_name_key,
      griis_iso3c = .data$griis_iso3c,
      griis_country = .data$griis_country,
      griis_matched_name = .data$griis_scientific_name,
      griis_is_invasive = .data$griis_is_invasive,
      griis_establishment_means = .data$griis_establishment_means,
      griis_occurrence_status = .data$griis_occurrence_status,
      griis_taxon_rank = .data$griis_taxon_rank
    )

  long <- dplyr::bind_rows(by_accepted, by_scientific) |>
    dplyr::filter(!is.na(.data$griis_join_name_key))

  group_cols <- if (isTRUE(require_country)) {
    c("griis_join_name_key", "griis_iso3c")
  } else {
    "griis_join_name_key"
  }

  long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      griis_listed = TRUE,
      griis_invasive = any(.data$griis_is_invasive %in% TRUE, na.rm = TRUE),
      griis_any_invasive_field_present = any(!is.na(.data$griis_is_invasive)),
      griis_n_matches = dplyr::n(),
      griis_matched_names = paste(sort(unique(stats::na.omit(.data$griis_matched_name))), collapse = ";"),
      griis_countries = paste(sort(unique(stats::na.omit(.data$griis_country))), collapse = ";"),
      griis_establishment_means = paste(sort(unique(stats::na.omit(.data$griis_establishment_means))), collapse = ";"),
      griis_occurrence_status = paste(sort(unique(stats::na.omit(.data$griis_occurrence_status))), collapse = ";"),
      griis_taxon_rank = paste(sort(unique(stats::na.omit(.data$griis_taxon_rank))), collapse = ";"),
      .groups = "drop"
    )
}

#' Attach GRIIS status to a species or species-country table
#'
#' Adds GRIIS listing and invasive-status evidence to an input table. The function
#' can join by species + country for terrestrial/freshwater workflows or by species
#' only for marine/global workflows.
#'
#' @details
#' When `require_country = TRUE`, at least one of `iso2c_col`, `iso3c_col` or
#' `country_col` must identify the recipient country. These fields are resolved to
#' ISO3 and joined against a species x country GRIIS lookup.
#'
#' When `require_country = FALSE`, no country parsing is attempted and the join is
#' species-only. This is important for marine pipelines where the input table may
#' contain only species names before marine spatial overlays are assigned.
#'
#' The function appends:
#'
#' - `griis_listed`: whether the species/country or species-only record appears in
#'   GRIIS;
#' - `griis_invasive`: whether any matching GRIIS row is flagged invasive;
#' - `griis_status`: a compact categorical summary of the match.
#'
#' @section Interpretation:
#' GRIIS columns appended by this helper are evidence fields for filtering and
#' audit. A row with `griis_listed = FALSE` means no matching GRIIS evidence was
#' found under the selected matching mode; it should not be interpreted as proof
#' that the species is native, absent or harmless in that recipient region.
#'
#' @param df Input data frame.
#' @param species_col Name of the column containing species names.
#' @param iso2c_col Optional name of an ISO2 recipient-country column.
#' @param iso3c_col Optional name of an ISO3 recipient-country column.
#' @param country_col Optional name of a recipient country-name column.
#' @param griis Optional raw or standardised GRIIS table. If `NULL`, GRIIS is read
#' with [bf_read_griis()].
#' @param cache_dir Cache directory used when `griis = NULL`. Must be supplied
#'   explicitly when the helper needs to download/read GRIIS internally. In
#'   examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, re-download/re-read GRIIS when the
#' helper loads it.
#' @param require_country Logical. If `TRUE`, join by species + recipient country.
#' If `FALSE`, join by species only.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return A tibble containing the input rows with GRIIS evidence columns
#'   appended. The returned object preserves the input columns and adds fields
#'   such as `griis_listed`, `griis_invasive`, `griis_status`,
#'   `griis_any_invasive_field_present`, `griis_n_matches` and supporting
#'   matched-name/status metadata where available. `griis_status` is a character
#'   summary with values such as `"listed_invasive"`,
#'   `"listed_introduced_or_alien"` and `"not_listed_or_no_match"`. These
#'   columns provide GRIIS evidence for filtering and auditing; non-matches
#'   should not be interpreted as confirmed native status, absence or lack of
#'   impact.
#'
#' @examples
#' raw <- data.frame(
#'   species = "Example species",
#'   scientificName = "Example species",
#'   countryCode_alpha2 = "GB",
#'   countryCode_alpha3 = "GBR",
#'   isInvasive = "yes",
#'   stringsAsFactors = FALSE
#' )
#' input <- data.frame(species = "Example species", iso2c = "GB")
#' bf_attach_griis_status(input, iso2c_col = "iso2c", griis = raw)
#'
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_attach_griis_status <- function(
    df,
    species_col = "species",
    iso2c_col = NULL,
    iso3c_col = NULL,
    country_col = NULL,
    griis = NULL,
    cache_dir = NULL,
    force_refresh = FALSE,
    require_country = TRUE,
    quiet = FALSE
) {
  bf_require_packages(c("dplyr", "tibble"), context = "GRIIS helper")

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  if (!species_col %in% names(df)) {
    stop("`species_col` was not found in `df`: ", species_col, call. = FALSE)
  }

  if (isTRUE(require_country)) {
    has_country <- any(c(iso2c_col, iso3c_col, country_col) %in% names(df))

    if (!has_country) {
      stop(
        "`require_country = TRUE`, but no valid `iso2c_col`, `iso3c_col`, or `country_col` was supplied.",
        call. = FALSE
      )
    }
  }

  if (is.null(griis)) {
    griis <- bf_read_griis(
      cache_dir = cache_dir,
      force_refresh = force_refresh,
      quiet = quiet
    )
  }

  input <- tibble::as_tibble(df)

  input$.bf_griis_rowid <- seq_len(nrow(input))
  input$.bf_griis_join_name_key <- .bf_griis_name_key(input[[species_col]])

  if (isTRUE(require_country)) {
    lookup <- bf_griis_lookup(griis, require_country = TRUE)

    input_iso2 <- if (!is.null(iso2c_col) && iso2c_col %in% names(input)) input[[iso2c_col]] else NULL
    input_iso3 <- if (!is.null(iso3c_col) && iso3c_col %in% names(input)) input[[iso3c_col]] else NULL
    input_country <- if (!is.null(country_col) && country_col %in% names(input)) input[[country_col]] else NULL

    input$.bf_griis_iso3c <- .bf_griis_country_to_iso3(
      country = input_country,
      iso2c = input_iso2,
      iso3c = input_iso3
    )

    joined <- input |>
      dplyr::left_join(
        lookup,
        by = c(
          ".bf_griis_join_name_key" = "griis_join_name_key",
          ".bf_griis_iso3c" = "griis_iso3c"
        )
      )
  } else {
    # Critical marine/global fix:
    # In species-only mode, do not create .bf_griis_iso3c and do not parse
    # country columns. This avoids zero-length country-vector assignment when
    # marine inputs contain species names only.
    lookup <- bf_griis_lookup(griis, require_country = FALSE)

    joined <- input |>
      dplyr::left_join(
        lookup,
        by = c(".bf_griis_join_name_key" = "griis_join_name_key")
      )
  }

  joined <- joined |>
    dplyr::mutate(
      griis_listed = dplyr::coalesce(.data$griis_listed, FALSE),
      griis_invasive = dplyr::coalesce(.data$griis_invasive, FALSE),
      griis_status = dplyr::case_when(
        .data$griis_invasive ~ "listed_invasive",
        .data$griis_listed ~ "listed_introduced_or_alien",
        TRUE ~ "not_listed_or_no_match"
      )
    ) |>
    dplyr::select(
      -dplyr::any_of(c(
        ".bf_griis_rowid",
        ".bf_griis_join_name_key",
        ".bf_griis_iso3c"
      ))
    )

  .bf_msg(
    "GRIIS join complete: ",
    sum(joined$griis_listed, na.rm = TRUE),
    " / ",
    nrow(joined),
    " row(s) listed; ",
    sum(joined$griis_invasive, na.rm = TRUE),
    " flagged invasive.",
    quiet = quiet
  )

  joined
}

#' Keep rows flagged as invasive in GRIIS
#'
#' Convenience wrapper around [bf_attach_griis_status()] that attaches GRIIS
#' evidence and then retains only rows where `griis_invasive` is `TRUE`.]
#'
#' @section Interpretation:
#' This is a conservative convenience filter based only on GRIIS invasive-status
#' evidence available under the selected matching mode. Rows removed by the
#' filter are not necessarily native or impact-free; they simply lack a matching
#' invasive GRIIS flag under the supplied criteria.
#'
#' @param df Input data frame.
#' @param species_col Name of the column containing species names.
#' @param iso2c_col Optional name of an ISO2 recipient-country column.
#' @param iso3c_col Optional name of an ISO3 recipient-country column.
#' @param country_col Optional name of a recipient country-name column.
#' @param griis Optional raw or standardised GRIIS table. If `NULL`, GRIIS is read
#' with [bf_read_griis()].
#' @param cache_dir Cache directory used when `griis = NULL`. Must be supplied
#'   explicitly when the helper needs to download/read GRIIS internally. In
#'   examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, re-download/re-read GRIIS when the
#' helper loads it.
#' @param require_country Logical. If `TRUE`, filter by species + recipient-country
#' GRIIS evidence. If `FALSE`, filter by species-level GRIIS evidence.
#' @param quiet Logical; if `TRUE`, suppress progress messages.
#'
#' @return A tibble containing the subset of input rows for which attached GRIIS
#'   evidence has `griis_invasive == TRUE`. The returned object includes the
#'   original input columns plus the GRIIS evidence columns added by
#'   [bf_attach_griis_status()], including `griis_listed`, `griis_invasive`,
#'   `griis_status` and supporting matched-name/status metadata where available.
#'   Rows excluded by this filter are not necessarily native or impact-free; they
#'   simply lack a matching invasive GRIIS flag under the selected matching mode.
#'
#' @examples
#' raw <- data.frame(
#'   species = "Example species",
#'   scientificName = "Example species",
#'   countryCode_alpha2 = "GB",
#'   countryCode_alpha3 = "GBR",
#'   isInvasive = "yes",
#'   stringsAsFactors = FALSE
#' )
#' input <- data.frame(species = "Example species", iso2c = "GB")
#' bf_filter_griis_invasive(input, iso2c_col = "iso2c", griis = raw)
#'
#'
#' @family GRIIS origin-evidence helpers
#'
#' @md
#' @export
bf_filter_griis_invasive <- function(
    df,
    species_col = "species",
    iso2c_col = NULL,
    iso3c_col = NULL,
    country_col = NULL,
    griis = NULL,
    cache_dir = NULL,
    force_refresh = FALSE,
    require_country = TRUE,
    quiet = FALSE
) {
  out <- bf_attach_griis_status(
    df = df,
    species_col = species_col,
    iso2c_col = iso2c_col,
    iso3c_col = iso3c_col,
    country_col = country_col,
    griis = griis,
    cache_dir = cache_dir,
    force_refresh = force_refresh,
    require_country = require_country,
    quiet = quiet
  )

  out[out$griis_invasive %in% TRUE, , drop = FALSE]
}
