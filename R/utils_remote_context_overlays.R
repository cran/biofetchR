###############################################################################
# utils_remote_context_overlays.R
# -----------------------------------------------------------------------------
# biofetchR remote context overlay utilities
# -----------------------------------------------------------------------------
# Remote-first loaders and enrichment helpers for terrestrial/freshwater context
# layers used by biofetchR occurrence-processing workflows.
#
# Main exported helpers:
#   * bf_load_gdw_barriers()       - Global Dam Watch barrier points.
#   * bf_load_biosphere_reserve()  - UNESCO biosphere reserve point locations.
#   * bf_load_gloric()             - GloRiC river reach lines.
#   * bf_enrich_raster_context()   - Raster-derived point covariates.
#
# Design principles:
#   * Prefer official/public remote sources plus local caching.
#   * Allow local paths or user-supplied URLs as explicit overrides.
#   * Use conservative validation for downloaded archives and error pages.
#   * Avoid shipping third-party data inside the package.
#   * Keep provider-specific licensing/citation responsibilities with the user.
#
###############################################################################

# -----------------------------------------------------------------------------
# Internal utilities
# -----------------------------------------------------------------------------

#' Validate an explicit cache directory
#'
#' @param cache_dir Candidate cache directory.
#' @param context Character label used in the error message.
#'
#' @return Normalised cache directory path.
#'
#' @keywords internal
#' @noRd
bf_require_explicit_cache_dir <- function(cache_dir, context = "this helper") {
  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly for ", context,
      ". In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  normalizePath(
    as.character(cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )
}

#' Guess the most appropriate vector file in a directory
#'
#' Searches extracted archives for common GIS vector formats and optionally ranks
#' candidates using filename patterns.
#'
#' @param path_or_dir Character path to a file or directory.
#' @param prefer_patterns Optional character vector of regular expressions used
#'   to prioritise candidate filenames.
#'
#' @return Character path to the selected vector file.
#' @keywords internal
#' @noRd
bf_guess_vector_file <- function(path_or_dir, prefer_patterns = NULL) {
  if (file.exists(path_or_dir) && !dir.exists(path_or_dir)) return(path_or_dir)
  stopifnot(dir.exists(path_or_dir))

  cand <- c(
    list.files(path_or_dir, pattern = "\\.gpkg$", recursive = TRUE, full.names = TRUE),
    list.files(path_or_dir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE),
    list.files(path_or_dir, pattern = "\\.geojson$", recursive = TRUE, full.names = TRUE),
    list.files(path_or_dir, pattern = "\\.json$", recursive = TRUE, full.names = TRUE),
    list.files(path_or_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  )

  if (!length(cand)) stop("No supported vector file found in: ", path_or_dir, call. = FALSE)

  if (!is.null(prefer_patterns) && length(prefer_patterns)) {
    bn <- basename(cand)
    score <- rep(0L, length(cand))
    for (i in seq_along(prefer_patterns)) {
      score <- score + as.integer(grepl(prefer_patterns[i], bn, ignore.case = TRUE)) * (length(prefer_patterns) - i + 1L)
    }
    cand <- cand[order(score, decreasing = TRUE)]
  }

  cand[1]
}



#' Read a vector layer from common file formats
#'
#' Reads spatial vector files with `sf`, or converts CSV files with recognisable
#' longitude/latitude columns to point `sf` objects.
#'
#' @param path Character path to a vector file.
#' @param quiet Logical. If `TRUE`, suppress read messages.
#'
#' @return An `sf` object.
#' @keywords internal
#' @noRd
bf_read_vector_any <- function(path, quiet = TRUE) {
  bf_require_packages(c("sf"), context = "remote context overlay")
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    x <- utils::read.csv(path, stringsAsFactors = FALSE)
    lon_col <- intersect(tolower(names(x)), c("longitude", "lon", "x", "decimallongitude"))
    lat_col <- intersect(tolower(names(x)), c("latitude", "lat", "y", "decimallatitude"))
    if (!length(lon_col) || !length(lat_col)) stop("CSV lacks recognisable lon/lat columns: ", path, call. = FALSE)
    lon_nm <- names(x)[match(lon_col[1], tolower(names(x)))]
    lat_nm <- names(x)[match(lat_col[1], tolower(names(x)))]
    sf::st_as_sf(x, coords = c(lon_nm, lat_nm), crs = 4326, remove = FALSE)
  } else {
    sf::read_sf(path, quiet = isTRUE(quiet))
  }
}

#' Clip or filter an overlay to requested countries
#'
#' Uses Natural Earth country boundaries as a lightweight optional spatial filter.
#' Point layers are filtered by intersection; non-point layers are cropped to the
#' target country bounding box.
#'
#' @param x An `sf` object.
#' @param iso2c Character vector of ISO2 country codes.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Filtered or cropped `sf` object.
#' @keywords internal
#' @noRd
bf_clip_to_countries <- function(x, iso2c, quiet = TRUE) {
  if (is.null(x) || !inherits(x, "sf") || !nrow(x) || !length(iso2c)) return(x)
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(x)

  world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
  keep <- world[world$iso_a2 %in% toupper(iso2c), , drop = FALSE]
  if (!nrow(keep)) return(x)

  keep <- sf::st_transform(keep, 4326)
  x <- bf_sf_wgs84(x)

  pre_n <- nrow(x)

  if (all(grepl("POINT", as.character(unique(sf::st_geometry_type(x))), ignore.case = TRUE))) {
    inside <- lengths(sf::st_intersects(x, keep)) > 0
    x <- x[inside, , drop = FALSE]
  } else {
    bb <- sf::st_as_sfc(sf::st_bbox(keep))
    x <- suppressWarnings(sf::st_crop(x, sf::st_bbox(bb)))
  }

  .bf_msg("Clip kept ", nrow(x), " / ", pre_n, " features for countries: ", paste(iso2c, collapse = ", "), quiet = quiet)
  x
}

#' Standardise identifier and label columns
#'
#' Finds the first matching candidate identifier/name columns and writes them to
#' standard output names expected by downstream overlay joins.
#'
#' @param x Data frame or `sf` object.
#' @param id_candidates Candidate input column names for identifiers.
#' @param name_candidates Candidate input column names for names/classes.
#' @param id_out Output identifier column name.
#' @param name_out Optional output name/class column name.
#'
#' @return `x` with standardised identifier and optional name columns added.
#' @keywords internal
#' @noRd
bf_normalize_columns <- function(x, id_candidates, name_candidates, id_out, name_out = NULL) {
  nms <- names(x)
  ln <- tolower(nms)
  pick <- function(cands) {
    for (cand in cands) {
      hit <- which(ln == tolower(cand))
      if (length(hit)) return(nms[hit[1]])
    }
    NULL
  }
  id_col <- pick(id_candidates)
  name_col <- pick(name_candidates)
  x[[id_out]] <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))
  if (!is.null(name_out)) x[[name_out]] <- if (!is.null(name_col)) as.character(x[[name_col]]) else x[[id_out]]
  x
}

# -----------------------------------------------------------------------------
# Official/default remote sources
# -----------------------------------------------------------------------------

#' Default remote source URLs for context overlays
#'
#' Internal constants used by the remote-first loaders when users do not provide
#' `source_path` or `source_url` overrides. These are not exported because source
#' URLs may change independently of the public user interface.
#'
#' @keywords internal
#' @noRd
BF_GLORIC_URL     <- "https://data.hydrosheds.org/file/hydrosheds-associated/gloric/GloRiC_v10_shapefile.zip"
BF_BIOSPHERE_CSV  <- "https://ihp-wins.unesco.org/dataset/98c49238-85e9-4a52-9a4a-d3e8c33993cc/resource/c52902df-09a3-4ece-b2ce-e581a17a7d25/download/mab_biospheres.csv"
BF_HII_URL_2020   <- "https://storage.googleapis.com/hii-export/2020-01-01/hii_2020-01-01.tif"
BF_HII_URL_2019   <- "https://storage.googleapis.com/hii-export/2019-01-01/hii_2019-01-01.tif"


#' Download a selected file from a Figshare article
#'
#' Fetches Figshare article metadata, selects the first file whose name matches a
#' requested pattern and downloads it into the cache.
#'
#' @param article_id Character or numeric Figshare article identifier.
#' @param cache_dir Character cache directory.
#' @param file_patterns Character vector of regular expressions matched against
#'   Figshare file names.
#' @param force_refresh Logical. If `TRUE`, refresh metadata and file cache.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Character path to the downloaded file.
#' @keywords internal
#' @noRd
bf_download_figshare_article_file <- function(article_id,
                                              cache_dir,
                                              file_patterns,
                                              force_refresh = FALSE,
                                              quiet = TRUE) {
  bf_require_packages(c("jsonlite"), context = "remote context overlay")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  meta_path <- file.path(cache_dir, paste0("figshare_article_", article_id, ".json"))
  if (!file.exists(meta_path) || isTRUE(force_refresh)) {
    bf_download_cached(
      sprintf("https://api.figshare.com/v2/articles/%s", article_id),
      meta_path,
      force_refresh = force_refresh,
      quiet = quiet,
      mode = "wb",
      min_bytes = 100
    )
  }

  meta <- jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
  files <- meta$files
  if (!length(files)) stop("Figshare article ", article_id, " returned no files.", call. = FALSE)

  file_names <- vapply(files, function(f) as.character(f$name %||% ""), character(1))
  keep <- rep(FALSE, length(file_names))
  for (pat in file_patterns) {
    keep <- keep | grepl(pat, file_names, ignore.case = TRUE)
  }

  if (!any(keep)) {
    stop(
      "Figshare article ", article_id,
      " contained no files matching pattern(s): ",
      paste(file_patterns, collapse = ", "),
      call. = FALSE
    )
  }

  f <- files[[which(keep)[1]]]
  nm <- as.character(f$name)
  url <- as.character(f$download_url %||% "")
  if (!nzchar(url)) stop("Figshare file ", nm, " has no download_url.", call. = FALSE)

  dest <- file.path(cache_dir, nm)
  bf_download_cached(url, dest, force_refresh = force_refresh, quiet = quiet)
  dest
}

# -----------------------------------------------------------------------------
# Vector overlay loaders
# -----------------------------------------------------------------------------

#' Load Global Dam Watch river barrier points
#'
#' Downloads, caches and reads the Global Dam Watch barrier layer, then returns a
#' standardised point `sf` object. By default, the function uses the public
#' Figshare article metadata to locate the Global Dam Watch v1.0 shapefile
#' archive. Local files or alternate URLs can be supplied with `source_path` or
#' `source_url` for testing, mirrors or manually downloaded data.
#'
#' @param cache_dir Character. Directory used for downloaded archives, extracted
#'   files and metadata caches. Must be supplied explicitly. In examples, tests
#'   and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, re-download remote files and rebuild
#'   extracted caches where possible.
#' @param iso2c Optional character vector of ISO2 country codes used to clip the
#'   returned layer when `clip_to_countries = TRUE`.
#' @param clip_to_countries Logical. If `TRUE` and `iso2c` is supplied, keep only
#'   features intersecting the requested countries.
#' @param source_path Optional local vector file or archive path. Used instead of
#'   downloading when the path exists.
#' @param source_url Optional remote URL to a vector file or archive. Used instead
#'   of the package default source.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` point object in WGS84 longitude/latitude coordinates
#'   containing Global Dam Watch barrier features. The returned object includes
#'   `gdw_barrier_id`, a stable barrier, dam or feature identifier;
#'   `gdw_barrier_name`, a human-readable barrier or dam label where available;
#'   and the active geometry column. If the selected source layer is not already
#'   point geometry, centroids are returned so occurrence points can be assigned
#'   to nearby barrier features by downstream workflows. If
#'   `clip_to_countries = TRUE` and `iso2c` is supplied, the returned object is
#'   filtered to the requested countries where possible. The function stops if no
#'   rows remain after processing or clipping.
#'
#' @section Data licensing and attribution:
#' biofetchR does not redistribute Global Dam Watch data. Users are responsible
#' for checking and following the licence, citation and attribution requirements
#' of the specific Global Dam Watch release used in their workflow.
#'
#' @family remote overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   gdw <- bf_load_gdw_barriers(
#'     iso2c = c("GB", "IE"),
#'     cache_dir = file.path(tempdir(), "biofetchR_gdw_barriers"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   gdw
#' }
#' }
#'
#' @export
bf_load_gdw_barriers <- function(cache_dir = NULL,
                                 force_refresh = FALSE,
                                 iso2c = NULL,
                                 clip_to_countries = TRUE,
                                 source_path = NULL,
                                 source_url = NULL,
                                 quiet = TRUE) {
  bf_require_packages(c("sf", "jsonlite"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_load_gdw_barriers()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  path <- NULL

  if (!is.null(source_path) && file.exists(source_path)) {
    path <- source_path
  } else {
    if (!is.null(source_url)) {
      ext <- if (grepl("\\.zip$", source_url, ignore.case = TRUE)) ".zip" else tools::file_ext(source_url)
      if (!nzchar(ext)) ext <- ".zip"
      arc <- file.path(cache_dir, paste0("gdw_barriers_download", if (startsWith(ext, ".")) ext else paste0(".", ext)))
      bf_download_cached(source_url, arc, force_refresh = force_refresh, quiet = quiet)
    } else {
      arc <- bf_download_figshare_article_file(
        article_id = "25988293",
        cache_dir = cache_dir,
        file_patterns = c("^GDW_v1_0_shp\\.zip$", "GDW.*shp.*\\.zip$"),
        force_refresh = force_refresh,
        quiet = quiet
      )
    }

    if (grepl("\\.(zip|tar\\.gz|tgz|gz)$", arc, ignore.case = TRUE)) {
      exdir <- file.path(cache_dir, "gdw_v1_0_shp_extracted")
      bf_unpack_archive(arc, exdir, force_refresh = force_refresh, quiet = quiet)
      path <- bf_guess_vector_file(
        exdir,
        prefer_patterns = c("barrier", "barriers", "dam", "point", "GDW")
      )
    } else {
      path <- arc
    }
  }

  x <- bf_read_vector_any(path, quiet = quiet)

  # If the guessed file was not a point layer, search the same extracted folder
  # for a point layer whose filename looks like the GDW barrier layer.
  if (!all(grepl("POINT", as.character(unique(sf::st_geometry_type(x))), ignore.case = TRUE))) {
    parent <- dirname(path)
    candidates <- c(
      list.files(parent, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE),
      list.files(parent, pattern = "\\.gpkg$", recursive = TRUE, full.names = TRUE)
    )
    candidates <- candidates[grepl("barrier|dam|point|gdw", basename(candidates), ignore.case = TRUE)]
    for (cand in candidates) {
      tmp <- tryCatch(bf_read_vector_any(cand, quiet = TRUE), error = function(e) NULL)
      if (inherits(tmp, "sf") && nrow(tmp) && all(grepl("POINT", as.character(unique(sf::st_geometry_type(tmp))), ignore.case = TRUE))) {
        x <- tmp
        path <- cand
        break
      }
    }
  }

  if (!all(grepl("POINT", as.character(unique(sf::st_geometry_type(x))), ignore.case = TRUE))) {
    x <- suppressWarnings(sf::st_centroid(x))
  }

  x <- bf_normalize_columns(
    x,
    id_candidates = c("gdw_barrier_id", "GDW_ID", "barrier_id", "BAR_ID", "DAM_ID", "dam_id", "GRAND_ID", "OBJECTID", "id", "ID"),
    name_candidates = c("gdw_barrier_name", "DAM_NAME", "dam_name", "barrier_name", "NAME", "name", "RES_NAME", "res_name"),
    id_out = "gdw_barrier_id",
    name_out = "gdw_barrier_name"
  )

  x <- bf_sf_make_valid(bf_sf_wgs84(x), quiet = quiet)
  gcol <- attr(x, "sf_column")
  keep <- intersect(c("gdw_barrier_id", "gdw_barrier_name", gcol), names(x))
  x <- x[, keep, drop = FALSE]

  if (isTRUE(clip_to_countries) && length(iso2c)) x <- bf_clip_to_countries(x, iso2c, quiet = quiet)
  if (!nrow(x)) stop("GDW barrier loader returned 0 rows after processing/clipping.", call. = FALSE)
  x
}

#' Load UNESCO biosphere reserve locations
#'
#' Downloads, caches and reads a UNESCO biosphere reserve CSV, then converts
#' recognised longitude/latitude columns to an `sf` point object. The function is
#' remote-first but accepts local or alternate remote sources for reproducible
#' workflows and tests.
#'
#' @param cache_dir Character. Directory used for the downloaded CSV cache. Must
#'   be supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, re-download the CSV even when a cache
#'   file is already present.
#' @param iso2c Optional character vector of ISO2 country codes used to clip the
#'   returned layer when `clip_to_countries = TRUE`.
#' @param clip_to_countries Logical. If `TRUE` and `iso2c` is supplied, keep only
#'   points intersecting the requested countries.
#' @param source_path Optional local CSV path. Used instead of downloading when
#'   the path exists.
#' @param source_url Optional remote CSV URL. Used instead of the package default
#'   source.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` point object in WGS84 longitude/latitude coordinates
#'   containing UNESCO biosphere reserve locations. The returned object includes
#'   `biosphere_id`, a stable reserve or row identifier; `biosphere_name`, a
#'   human-readable reserve label where available; the original longitude and
#'   latitude columns used to create point geometry; and the active geometry
#'   column. If `clip_to_countries = TRUE` and `iso2c` is supplied, the returned
#'   object is filtered to biosphere reserve points intersecting the requested
#'   countries where possible.
#'
#' @section Data licensing and attribution:
#' biofetchR does not redistribute UNESCO biosphere reserve data. Users are
#' responsible for checking and following the licence, citation and attribution
#' requirements of the source used in their workflow.
#'
#' @family remote overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   biosphere <- bf_load_biosphere_reserve(
#'     iso2c = "FR",
#'     cache_dir = file.path(tempdir(), "biofetchR_biosphere"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   biosphere
#' }
#' }
#'
#' @export
bf_load_biosphere_reserve <- function(cache_dir = NULL,
                                      force_refresh = FALSE,
                                      iso2c = NULL,
                                      clip_to_countries = TRUE,
                                      source_path = NULL,
                                      source_url = NULL,
                                      quiet = TRUE) {
  bf_require_packages(c("sf"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_load_biosphere_reserve()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- NULL
  if (!is.null(source_path) && file.exists(source_path)) {
    path <- source_path
  } else {
    src <- if (is.null(source_url)) BF_BIOSPHERE_CSV else source_url
    csv_path <- file.path(cache_dir, "biosphere_reserves.csv")
    bf_download_cached(src, csv_path, force_refresh = force_refresh, quiet = quiet, mode = "wb")
    path <- csv_path
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  nms <- names(x); ln <- tolower(nms)
  pick <- function(cands) {
    for (cand in cands) {
      hit <- which(ln == tolower(cand))
      if (length(hit)) return(nms[hit[1]])
    }
    NULL
  }
  lon_col <- pick(c("longitude", "lon", "x", "decimalLongitude", "long"))
  lat_col <- pick(c("latitude", "lat", "y", "decimalLatitude"))
  if (is.null(lon_col) || is.null(lat_col)) stop("UNESCO biosphere CSV did not contain recognisable longitude/latitude columns.", call. = FALSE)
  id_col <- pick(c("id", "ID", "site_id", "identifier", "OBJECTID"))
  name_col <- pick(c("name", "Name", "site_name", "biosphere_reserve", "title"))
  x$biosphere_id <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))
  x$biosphere_name <- if (!is.null(name_col)) as.character(x[[name_col]]) else x$biosphere_id
  sf_x <- sf::st_as_sf(x, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
  sf_x <- bf_sf_make_valid(bf_sf_wgs84(sf_x), quiet = quiet)
  gcol <- attr(sf_x, "sf_column")
  keep <- intersect(c("biosphere_id", "biosphere_name", lon_col, lat_col, gcol), names(sf_x))
  sf_x <- sf_x[, keep, drop = FALSE]
  if (isTRUE(clip_to_countries) && length(iso2c)) sf_x <- bf_clip_to_countries(sf_x, iso2c, quiet = quiet)
  sf_x
}

#' Load GloRiC river reach lines
#'
#' Downloads, caches and reads the Global River Classification (GloRiC) line
#' layer, then standardises core identifier/class columns for use as a
#' terrestrial/freshwater context overlay.
#'
#' @param cache_dir Character. Directory used for downloaded archives and
#'   extracted files. Must be supplied explicitly. In examples, tests and
#'   vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, re-download remote files and rebuild
#'   extracted caches where possible.
#' @param iso2c Optional character vector of ISO2 country codes used to clip the
#'   returned layer when `clip_to_countries = TRUE`.
#' @param clip_to_countries Logical. If `TRUE` and `iso2c` is supplied, crop the
#'   river layer to the requested country boundaries.
#' @param source_path Optional local vector file or archive path. Used instead of
#'   downloading when the path exists.
#' @param source_url Optional remote URL to a vector file or archive. Used instead
#'   of the package default source.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` line object in WGS84 longitude/latitude coordinates
#'   containing GloRiC river-reach features. The returned object includes
#'   `gloric_id`, a stable reach identifier; `gloric_class`, the selected
#'   GloRiC reach-type or classification label where available; and the active
#'   geometry column. Each row represents one river reach line or multiline
#'   feature used as terrestrial/freshwater context. If
#'   `clip_to_countries = TRUE` and `iso2c` is supplied, the returned object is
#'   cropped to the requested countries where possible. The function stops if no
#'   rows remain after processing or clipping.
#'
#' @section Data licensing and attribution:
#' biofetchR does not redistribute GloRiC data. Users are responsible for checking
#' and following the licence, citation and attribution requirements of the GloRiC
#' version used in their workflow.
#'
#' @family remote overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   rivers <- bf_load_gloric(
#'     iso2c = "GB",
#'     cache_dir = file.path(tempdir(), "biofetchR_gloric"),
#'     clip_to_countries = TRUE
#'   )
#'
#'   rivers
#' }
#' }
#'
#' @export
bf_load_gloric <- function(cache_dir = NULL,
                           force_refresh = FALSE,
                           iso2c = NULL,
                           clip_to_countries = TRUE,
                           source_path = NULL,
                           source_url = NULL,
                           quiet = TRUE) {
  bf_require_packages(c("sf"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_load_gloric()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  path <- NULL

  if (!is.null(source_path) && file.exists(source_path)) {
    path <- source_path
  } else {
    src <- if (is.null(source_url)) BF_GLORIC_URL else source_url
    arc <- file.path(cache_dir, basename(src))
    if (!grepl("\\.(zip|shp|gpkg|geojson|json)$", arc, ignore.case = TRUE)) {
      arc <- paste0(arc, ".zip")
    }

    bf_download_cached(src, arc, force_refresh = force_refresh, quiet = quiet)

    if (grepl("\\.(zip|tar\\.gz|tgz|gz)$", arc, ignore.case = TRUE)) {
      exdir <- file.path(cache_dir, "gloric_extracted")
      bf_unpack_archive(arc, exdir, force_refresh = force_refresh, quiet = quiet)
      path <- bf_guess_vector_file(
        exdir,
        prefer_patterns = c("gloric", "reach", "river", "hyriv")
      )
    } else {
      path <- arc
    }
  }

  x <- bf_read_vector_any(path, quiet = quiet)

  x <- bf_normalize_columns(
    x,
    id_candidates = c("gloric_id", "GLORIC_ID", "HYRIV_ID", "REACH_ID", "id", "ID"),
    name_candidates = c("Reach_type", "reach_type", "gloric_class", "Class_hydr", "Kmeans_30"),
    id_out = "gloric_id",
    name_out = "gloric_class"
  )

  x <- bf_sf_make_valid(bf_sf_wgs84(x), quiet = quiet)

  gcol <- attr(x, "sf_column")
  keep <- intersect(c("gloric_id", "gloric_class", gcol), names(x))
  x <- x[, keep, drop = FALSE]

  .bf_msg("GLORIC source path: ", path, quiet = quiet)
  .bf_msg("GLORIC rows after processing: ", nrow(x), quiet = quiet)
  .bf_msg(
    "GLORIC geometry types: ",
    paste(unique(as.character(sf::st_geometry_type(x))), collapse = ", ", quiet = quiet)
  )
  .bf_msg(
    "GLORIC columns present: ",
    paste(intersect(c("gloric_id", "gloric_class"), names(x)), collapse = ", ", quiet = quiet)
  )

  if (isTRUE(clip_to_countries) && length(iso2c)) {
    pre_n <- nrow(x)
    x <- bf_clip_to_countries(x, iso2c, quiet = quiet)
    .bf_msg("GLORIC rows after clip: ", nrow(x), " / ", pre_n, quiet = quiet)
  }

  if (!nrow(x)) {
    stop("GLORIC loader returned 0 rows after clipping.", call. = FALSE)
  }

  x
}

# -----------------------------------------------------------------------------
# Raster helpers
# -----------------------------------------------------------------------------

#' Download all files from a Zenodo record
#'
#' Fetches Zenodo record metadata and downloads all linked files into a cache
#' directory.
#'
#' @param record_id Character or numeric Zenodo record identifier.
#' @param cache_dir Character cache directory.
#' @param force_refresh Logical. If `TRUE`, refresh metadata and files.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Character vector of downloaded file paths.
#' @keywords internal
#' @noRd
bf_download_zenodo_record <- function(record_id, cache_dir, force_refresh = FALSE, quiet = TRUE) {
  bf_require_packages(c("jsonlite"), context = "remote context overlay")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  meta_path <- file.path(cache_dir, paste0("zenodo_record_", record_id, ".json"))
  if (!file.exists(meta_path) || isTRUE(force_refresh)) {
    bf_download_cached(sprintf("https://zenodo.org/api/records/%s", record_id), meta_path, force_refresh = force_refresh, quiet = quiet, mode = "wb")
  }
  meta <- jsonlite::fromJSON(meta_path, simplifyVector = FALSE)
  files <- meta$files
  if (!length(files)) stop("Zenodo record ", record_id, " returned no files.", call. = FALSE)
  out <- character()
  for (f in files) {
    key <- f$key
    link <- f$links$self
    dest <- file.path(cache_dir, key)
    if (!file.exists(dest) || isTRUE(force_refresh)) bf_download_cached(link, dest, force_refresh = force_refresh, quiet = quiet)
    out <- c(out, dest)
  }
  out
}

#' Prepare a WorldCover raster source
#'
#' Resolves a local, user-supplied remote or default Zenodo WorldCover source to
#' one raster path. Multiple downloaded tiles are combined into a VRT.
#'
#' @param cache_dir Character cache directory.
#' @param year Numeric or character WorldCover year.
#' @param source Optional local path or URL.
#' @param force_refresh Logical. If `TRUE`, refresh downloads and extractions.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Character path to a GeoTIFF or VRT.
#' @keywords internal
#' @noRd
bf_prepare_worldcover <- function(cache_dir = NULL,
                                  year = 2020,
                                  source = NULL,
                                  force_refresh = FALSE,
                                  quiet = TRUE) {
  bf_require_packages(c("terra"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_prepare_worldcover()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  if (!is.null(source)) {
    if (file.exists(source)) {
      files <- source
    } else {
      dest <- file.path(cache_dir, basename(source))
      bf_download_cached(source, dest, force_refresh = force_refresh, quiet = quiet)
      files <- dest
    }
  } else {
    record_id <- if (as.integer(year) == 2020) 5571936 else 7254221
    wc_dir <- file.path(cache_dir, paste0("worldcover_", year))
    files <- bf_download_zenodo_record(record_id, wc_dir, force_refresh = force_refresh, quiet = quiet)
  }
  tif_files <- character()
  for (f in files) {
    if (grepl("\\.(tif|tiff)$", f, ignore.case = TRUE)) {
      tif_files <- c(tif_files, f)
    } else if (grepl("\\.(zip|tar\\.gz|tgz|gz)$", f, ignore.case = TRUE)) {
      exdir <- file.path(dirname(f), paste0(tools::file_path_sans_ext(basename(f)), "_extracted"))
      bf_unpack_archive(f, exdir, force_refresh = force_refresh, quiet = quiet)
      tif_files <- c(tif_files, list.files(exdir, pattern = "\\.(tif|tiff)$", recursive = TRUE, full.names = TRUE))
    }
  }
  tif_files <- unique(tif_files)
  if (!length(tif_files)) stop("WorldCover download completed but no GeoTIFF files were found.", call. = FALSE)
  if (length(tif_files) == 1) return(tif_files[1])
  vrt_path <- file.path(cache_dir, paste0("worldcover_", year, ".vrt"))
  terra::vrt(tif_files, filename = vrt_path, overwrite = TRUE)
  vrt_path
}

#' Build default SoilGrids raster URLs
#'
#' Constructs default 0--5 cm mean, 5 km SoilGrids GeoTIFF URLs for requested
#' variable names.
#'
#' @param vars Character vector of SoilGrids variable names.
#'
#' @return Named character vector of URLs.
#' @keywords internal
#' @noRd
bf_default_soilgrids_urls <- function(vars) {
  vars <- unique(as.character(vars))
  out <- setNames(vector("character", length(vars)), vars)
  for (v in vars) {
    out[[v]] <- sprintf("https://files.isric.org/soilgrids/latest/data_aggregated/5000m/%s/%s_0-5cm_mean_5000.tif", v, v)
  }
  out
}

#' Prepare SoilGrids raster sources
#'
#' Resolves named local or remote SoilGrids sources to cached raster paths, or
#' builds default SoilGrids URLs from requested variable names.
#'
#' @param cache_dir Character cache directory.
#' @param sources Optional named vector/list of local paths or URLs.
#' @param vars Optional character vector of SoilGrids variable names.
#' @param force_refresh Logical. If `TRUE`, refresh downloads.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Named list of raster file paths.
#' @keywords internal
#' @noRd
bf_prepare_soilgrids <- function(cache_dir = NULL,
                                 sources = NULL,
                                 vars = NULL,
                                 force_refresh = FALSE,
                                 quiet = TRUE) {
  bf_require_packages(c("terra"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_prepare_soilgrids()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(sources)) {
    if (is.null(vars) || !length(vars)) stop("Provide `soilgrids_vars` (and optionally `soilgrids_sources`) to use SoilGrids enrichment.", call. = FALSE)
    sources <- bf_default_soilgrids_urls(vars)
  }
  if (is.null(names(sources)) || any(!nzchar(names(sources)))) stop("`soilgrids_sources` must be a named character vector/list keyed by variable.", call. = FALSE)
  out <- list()
  for (nm in names(sources)) {
    src <- sources[[nm]]
    if (file.exists(src)) {
      out[[nm]] <- src
    } else {
      dest <- file.path(cache_dir, paste0(nm, ".tif"))
      bf_download_cached(src, dest, force_refresh = force_refresh, quiet = quiet)
      out[[nm]] <- dest
    }
  }
  out
}

#' Prepare a Human Footprint raster source
#'
#' Resolves a local, user-supplied remote or default Human Impact/Human Footprint
#' raster source to one cached raster path.
#'
#' @param cache_dir Character cache directory.
#' @param source Optional local path or URL.
#' @param force_refresh Logical. If `TRUE`, refresh downloads.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return Character path to a raster file.
#' @keywords internal
#' @noRd
bf_prepare_human_footprint <- function(cache_dir = NULL,
                                       source = NULL,
                                       force_refresh = FALSE,
                                       quiet = TRUE) {
  bf_require_packages(c("terra"), context = "remote context overlay")
  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_prepare_human_footprint()")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  src <- if (is.null(source)) BF_HII_URL_2020 else source
  if (file.exists(src)) return(src)
  dest <- file.path(cache_dir, basename(src))
  if (!grepl("\\.(tif|tiff)$", dest, ignore.case = TRUE)) dest <- file.path(cache_dir, "human_footprint_2020.tif")
  bf_download_cached(src, dest, force_refresh = force_refresh, quiet = quiet)
  dest
}

#' Extract raster values for point features
#'
#' Internal wrapper around `terra::extract()` that adds one output column to an
#' input point `sf` object. Optional buffers can be used for local summaries.
#'
#' @param sf_points Point `sf` object.
#' @param raster_path Character raster path readable by `terra::rast()`.
#' @param column_name Output column name.
#' @param buffer_m Numeric buffer radius in raster units/metres, depending on CRS.
#' @param fun Optional summary function passed to `terra::extract()`.
#'
#' @return Input `sf` object with one additional column.
#' @keywords internal
#' @noRd
bf_extract_raster_values <- function(sf_points, raster_path, column_name, buffer_m = 0, fun = NULL) {
  bf_require_packages(c("terra", "sf"), context = "remote context overlay")
  stopifnot(inherits(sf_points, "sf"))
  pts <- sf_points
  if (!".bf_rowid" %in% names(pts)) pts$.bf_rowid <- seq_len(nrow(pts))
  r <- terra::rast(raster_path)
  if (is.null(fun)) fun <- if (buffer_m > 0) mean else NULL
  pts_v <- terra::vect(pts)
  if (buffer_m > 0) pts_v <- terra::buffer(pts_v, width = buffer_m)
  vals <- terra::extract(r, pts_v, fun = fun, ID = FALSE)
  vals <- as.data.frame(vals)
  if (!nrow(vals)) {
    pts[[column_name]] <- NA_real_
    return(pts)
  }
  if (ncol(vals) > 1) vals <- vals[, 1, drop = FALSE]
  names(vals)[1] <- column_name
  pts[[column_name]] <- vals[[column_name]]
  pts
}

#' Enrich GBIF points with raster context layers
#'
#' Adds raster-derived covariates to an `sf` point object. Supported contexts are
#' currently WorldCover, SoilGrids and Human Footprint. Raster files are prepared
#' through remote-first helper functions unless local sources are supplied.
#'
#' @param sf_points Point `sf` object. Empty inputs, non-`sf` inputs or calls with
#'   no requested `raster_context` are returned unchanged where applicable.
#' @param raster_context Character vector of raster contexts to extract. Supported
#'   values are `"worldcover"`, `"soilgrids"` and `"human_footprint"`.
#' @param cache_dir Character. Directory used for raster download and preparation
#'   caches. Must be supplied explicitly when raster context layers are requested.
#'   In examples, tests and vignettes, use a path under `tempdir()`.
#' @param worldcover_source Optional local path or URL for a WorldCover raster.
#'   When `NULL`, the default WorldCover Zenodo record for `worldcover_year` is
#'   used.
#' @param worldcover_year Numeric or character WorldCover year. Used to choose the
#'   default source and to name the output column.
#' @param worldcover_buffer_m Numeric. Buffer radius, in metres, for modal
#'   WorldCover extraction. A value of `0` performs point extraction.
#' @param soilgrids_sources Optional named vector/list of SoilGrids raster paths
#'   or URLs. Names are used in output column names.
#' @param soilgrids_vars Optional character vector of SoilGrids variable names.
#'   Used to build default SoilGrids URLs when `soilgrids_sources = NULL`.
#' @param soilgrids_buffer_m Numeric. Buffer radius, in metres, for mean
#'   SoilGrids extraction. A value of `0` performs point extraction.
#' @param human_footprint_source Optional local path or URL for a Human Footprint
#'   raster. When `NULL`, the package default source is used.
#' @param human_footprint_buffer_m Numeric. Buffer radius, in metres, for mean
#'   Human Footprint extraction. A value of `0` performs point extraction.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` object with the same rows and geometry as `sf_points`, with
#'   additional raster-context columns appended for the requested contexts. When
#'   `"worldcover"` is requested, the returned object includes a
#'   `worldcover_<year>` column containing extracted WorldCover class values for
#'   the selected `worldcover_year`. When `"soilgrids"` is requested, one column
#'   is added for each prepared SoilGrids variable, using names of the form
#'   `soilgrids_<variable>`. When `"human_footprint"` is requested, the returned
#'   object includes a `human_footprint` column. Raster values are extracted
#'   directly at point locations when the relevant buffer distance is `0`;
#'   otherwise they are summarised within buffers using modal extraction for
#'   WorldCover and mean extraction for SoilGrids and Human Footprint. Empty
#'   inputs, non-`sf` inputs or calls with no requested raster contexts are
#'   returned unchanged.
#'
#' @section Data licensing and attribution:
#' This function can download or extract values from third-party raster products,
#' including WorldCover, SoilGrids and Human Footprint layers. biofetchR does not
#' ship, redistribute or modify these raster datasets. Users are responsible for
#' checking the licence, citation and attribution requirements of each raster
#' source before use in analysis, publication or redistribution. Always cite the
#' specific product, version, year and provider used in the workflow.
#'
#' @family raster context helpers
#'
#' @examples
#' if (
#'   requireNamespace("sf", quietly = TRUE) &&
#'   requireNamespace("terra", quietly = TRUE)
#' ) {
#'   gbif_points <- sf::st_as_sf(
#'     data.frame(
#'       species = "Example species",
#'       decimalLongitude = -5,
#'       decimalLatitude = 55
#'     ),
#'     coords = c("decimalLongitude", "decimalLatitude"),
#'     crs = 4326,
#'     remove = FALSE
#'   )
#'
#'   r <- terra::rast(
#'     nrows = 2,
#'     ncols = 2,
#'     xmin = -10,
#'     xmax = 0,
#'     ymin = 50,
#'     ymax = 60,
#'     crs = "EPSG:4326"
#'   )
#'   terra::values(r) <- c(10, 20, 30, 40)
#'
#'   raster_path <- file.path(tempdir(), "worldcover_example.tif")
#'   terra::writeRaster(r, raster_path, overwrite = TRUE)
#'
#'   enriched <- bf_enrich_raster_context(
#'     sf_points = gbif_points,
#'     raster_context = "worldcover",
#'     cache_dir = tempdir(),
#'     worldcover_source = raster_path,
#'     worldcover_year = 2020
#'   )
#'
#'   names(enriched)
#' }
#'
#' @export
bf_enrich_raster_context <- function(sf_points,
                                     raster_context = character(0),
                                     cache_dir = NULL,
                                     worldcover_source = NULL,
                                     worldcover_year = 2020,
                                     worldcover_buffer_m = 0,
                                     soilgrids_sources = NULL,
                                     soilgrids_vars = NULL,
                                     soilgrids_buffer_m = 0,
                                     human_footprint_source = NULL,
                                     human_footprint_buffer_m = 0,
                                     quiet = TRUE) {
  bf_require_packages(c("sf", "terra"), context = "remote context overlay")

  if (!inherits(sf_points, "sf") || !nrow(sf_points) || !length(raster_context)) {
    return(sf_points)
  }

  cache_dir <- bf_require_explicit_cache_dir(cache_dir, context = "bf_enrich_raster_context()")

  out <- sf_points
  if ("worldcover" %in% raster_context) {
    wc <- bf_prepare_worldcover(cache_dir = file.path(cache_dir, "worldcover"), year = worldcover_year, source = worldcover_source, quiet = quiet)
    out <- bf_extract_raster_values(out, wc, column_name = paste0("worldcover_", worldcover_year), buffer_m = worldcover_buffer_m, fun = if (worldcover_buffer_m > 0) terra::modal else NULL)
  }
  if ("soilgrids" %in% raster_context) {
    sg <- bf_prepare_soilgrids(cache_dir = file.path(cache_dir, "soilgrids"), sources = soilgrids_sources, vars = soilgrids_vars, quiet = quiet)
    for (nm in names(sg)) {
      out <- bf_extract_raster_values(out, sg[[nm]], column_name = paste0("soilgrids_", nm), buffer_m = soilgrids_buffer_m, fun = if (soilgrids_buffer_m > 0) mean else NULL)
    }
  }
  if ("human_footprint" %in% raster_context) {
    hf <- bf_prepare_human_footprint(cache_dir = file.path(cache_dir, "human_footprint"), source = human_footprint_source, quiet = quiet)
    out <- bf_extract_raster_values(out, hf, column_name = "human_footprint", buffer_m = human_footprint_buffer_m, fun = if (human_footprint_buffer_m > 0) mean else NULL)
  }
  out
}
