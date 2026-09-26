###############################################################################
# utils_raster_context_layers.R
# -----------------------------------------------------------------------------
# biofetchR raster context enrichment helpers
# -----------------------------------------------------------------------------
# Adds raster-derived environmental and anthropogenic context variables to point
# `sf` occurrence data. These helpers are designed for enrichment rather than
# spatial regionalisation: they add columns to the input points but do not split
# records into separate terrestrial, freshwater, marine, or administrative units.
#
# Supported raster contexts
#   * ESA WorldCover-style categorical land-cover classes.
#   * SoilGrids-style continuous soil variables.
#   * Human Footprint-style continuous anthropogenic-pressure rasters.
#
# Design notes
#   * Inputs are intentionally flexible: local file path, URL, or pre-loaded
#     `terra::SpatRaster` object where possible.
#   * Raster downloads are cached locally so repeated pipeline runs do not
#     repeatedly download the same large files.
#   * Point extraction is used when buffer size is zero.
#   * Buffered extraction is used when a positive buffer radius is supplied:
#       - modal extraction for categorical WorldCover;
#       - mean extraction for continuous SoilGrids and Human Footprint layers.
#   * This file intentionally keeps all helper functions internal except the
#     public enrichment function `bf_enrich_raster_context()`.
#
###############################################################################

#' Load or download a raster source
#'
#' Internal raster-source resolver. Accepts an already loaded
#' `terra::SpatRaster`, a local file path, or a URL. URL inputs are downloaded
#' into a context-specific cache directory before being opened with
#' [terra::rast()].
#'
#' @param x Raster source. May be a `terra::SpatRaster`, a local path, or a URL.
#' @param cache_dir Directory used to cache downloaded raster files.
#' @param stem File-name stem used when saving a downloaded raster.
#' @param force_refresh Logical. If `TRUE`, re-download URL sources even when a
#'   cached copy already exists.
#' @param quiet Logical. Passed to [utils::download.file()] to suppress download
#'   messages.
#'
#' @return A `terra::SpatRaster`.
#'
#' @keywords internal
#' @noRd
.bf_rast_get <- function(x, cache_dir, stem, force_refresh = FALSE, quiet = TRUE) {
  if (inherits(x, "SpatRaster")) return(x)
  if (is.null(x) || !length(x) || !nzchar(as.character(x)[1])) stop("A raster source must be provided.", call. = FALSE)
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required for raster enrichment.", call. = FALSE)
  src <- as.character(x)[1]
  if (file.exists(src)) return(terra::rast(src))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  ext <- tools::file_ext(src)
  if (!nzchar(ext)) ext <- "tif"
  dest <- file.path(cache_dir, paste0(stem, ".", ext))
  if (!file.exists(dest) || isTRUE(force_refresh)) {
    bf_download_cached(
      url = src,
      dest = dest,
      force_refresh = force_refresh,
      quiet = quiet,
      min_bytes = 1024,
      validate_not_html = TRUE
    )
  }
  terra::rast(dest)
}

#' Extract raster values at points or within buffers
#'
#' Internal extraction helper used by all raster-context layers. With
#' `buffer_m = 0`, values are extracted directly at point locations. With a
#' positive buffer, the helper extracts a summary statistic over buffered
#' geometries using the supplied `fun` argument.
#'
#' @param r A `terra::SpatRaster`.
#' @param pts_sf Point `sf` object.
#' @param buffer_m Numeric buffer radius in metres. Use `0` for point extraction.
#' @param fun Optional summary function passed to [terra::extract()] for
#'   buffered extraction, for example [base::mean()] or [terra::modal()].
#' @param quiet Logical. Included for interface consistency with the other
#'   raster helpers.
#'
#' @return A data frame returned by [terra::extract()], including the extraction
#'   `ID` column and one or more raster-value columns.
#'
#' @keywords internal
#' @noRd
.bf_rast_extract_point_or_buffer <- function(r, pts_sf, buffer_m = 0, fun = NULL, quiet = TRUE) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required.", call. = FALSE)
  pts <- pts_sf
  if (!inherits(pts, "sf")) stop("sf_points must be an sf object.", call. = FALSE)
  if (is.na(sf::st_crs(pts))) pts <- sf::st_set_crs(pts, 4326)
  if (sf::st_crs(pts) != sf::st_crs(terra::crs(r, proj = TRUE))) {
    pts <- tryCatch(sf::st_transform(pts, terra::crs(r, proj = TRUE)), error = function(e) pts)
  }
  v <- terra::vect(pts)
  if (isTRUE(buffer_m > 0)) {
    if (is.na(sf::st_crs(pts)) || sf::st_is_longlat(pts)) {
      pts <- sf::st_transform(pts, 3857)
      v <- terra::vect(sf::st_buffer(pts, buffer_m))
    } else {
      v <- terra::buffer(v, width = buffer_m)
    }
  }
  out <- terra::extract(r, v, fun = fun, na.rm = TRUE)
  out
}

#' ESA WorldCover class-code lookup table
#'
#' Named character vector used to convert integer WorldCover class codes to
#' human-readable land-cover labels after raster extraction.
#'
#' @return Named character vector, stored internally.
#'
#' @keywords internal
#' @noRd
.bf_worldcover_labels <- c(
  `10` = "Tree cover",
  `20` = "Shrubland",
  `30` = "Grassland",
  `40` = "Cropland",
  `50` = "Built-up",
  `60` = "Bare / sparse vegetation",
  `70` = "Snow and ice",
  `80` = "Permanent water bodies",
  `90` = "Herbaceous wetland",
  `95` = "Mangroves",
  `100` = "Moss and lichen"
)

#' Enrich point occurrences from explicitly supplied raster context layers
#'
#' Adds one or more raster-derived contextual variables to an input point
#' `sf` object using raster sources supplied directly by the user. This helper is
#' intended for biofetchR occurrence-enrichment workflows where environmental or
#' anthropogenic raster covariates need to be attached to GBIF-derived point
#' records before downstream modelling, mapping, filtering or summary steps.
#'
#' @details
#' `bf_enrich_raster_context_from_sources()` is the source-explicit companion to
#' [bf_enrich_raster_context()]. Use this function when the raster products are
#' already available locally, when you want to provide custom raster URLs, or when
#' you want to control the exact raster versions used in the analysis.
#'
#' For package-managed/default raster downloads, use [bf_enrich_raster_context()]
#' instead.
#'
#' This function currently supports three optional raster-context groups:
#'
#' \describe{
#'   \item{`"worldcover"`}{Extracts an ESA WorldCover-style categorical raster,
#'   returning `worldcover_class`, `worldcover_label` and `worldcover_year`.}
#'   \item{`"soilgrids"`}{Extracts one or more named continuous SoilGrids-style
#'   rasters and stores them as `soilgrids_<variable>` columns.}
#'   \item{`"human_footprint"`}{Extracts a continuous Human Footprint-style
#'   raster and stores it as `human_footprint`.}
#' }
#'
#' Raster sources may be local file paths, URLs or pre-loaded
#' `terra::SpatRaster` objects. URL sources are downloaded into `cache_dir` and
#' reused on later runs unless `force_refresh = TRUE`.
#'
#' Extraction is performed directly at point coordinates when the relevant buffer
#' argument is `0`. If a positive buffer radius is supplied, raster values are
#' summarised within buffered geometries: modal class for WorldCover and mean
#' value for SoilGrids and Human Footprint.
#'
#' @param sf_points Point `sf` object in any coordinate reference system.
#' @param raster_context Character vector selecting one or more raster contexts.
#'   Supported values are `"worldcover"`, `"soilgrids"` and
#'   `"human_footprint"`.
#' @param cache_dir Directory used to cache downloaded raster files when any
#'   supplied raster source is a URL. Must be supplied explicitly if URL raster
#'   sources are used. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param worldcover_source Path, URL or `terra::SpatRaster` for the WorldCover
#'   raster. Required when `"worldcover"` is requested.
#' @param worldcover_year Numeric or character label stored in the output
#'   `worldcover_year` column. This is metadata only and does not automatically
#'   select or download a raster.
#' @param worldcover_buffer_m Numeric buffer radius in metres for modal
#'   WorldCover extraction. Use `0` for direct point extraction.
#' @param soilgrids_sources Named list or named vector of SoilGrids raster
#'   sources. Names are used to construct output columns of the form
#'   `soilgrids_<variable>`. Required when `"soilgrids"` is requested.
#' @param soilgrids_vars Optional character vector of SoilGrids variable names to
#'   extract. Defaults to all names in `soilgrids_sources`.
#' @param soilgrids_buffer_m Numeric buffer radius in metres for mean SoilGrids
#'   extraction. Use `0` for direct point extraction.
#' @param human_footprint_source Path, URL or `terra::SpatRaster` for the Human
#'   Footprint raster. Required when `"human_footprint"` is requested.
#' @param human_footprint_buffer_m Numeric buffer radius in metres for mean
#'   Human Footprint extraction. Use `0` for direct point extraction.
#' @param force_refresh Logical. If `TRUE`, re-download URL raster sources even
#'   when cached files already exist.
#' @param quiet Logical. If `TRUE`, suppress routine messages and download output
#'   where possible.
#'
#' @return An `sf` object with the same rows and geometry as `sf_points`, with
#'   additional raster-context columns appended for the requested contexts. When
#'   `"worldcover"` is requested, the returned object includes
#'   `worldcover_class`, `worldcover_label` and `worldcover_year`. When
#'   `"soilgrids"` is requested, one numeric column is added for each requested
#'   SoilGrids variable, using names of the form `soilgrids_<variable>`. When
#'   `"human_footprint"` is requested, the returned object includes a numeric
#'   `human_footprint` column. Raster values are extracted directly at point
#'   locations when the relevant buffer distance is `0`; otherwise they are
#'   summarised within buffers using modal extraction for categorical WorldCover
#'   data and mean extraction for continuous SoilGrids and Human Footprint data.
#'   No rows are intentionally added or removed by this enrichment helper.
#'
#' @section Output columns:
#' Depending on the requested raster contexts, the returned object may include:
#'
#' \describe{
#'   \item{`worldcover_class`}{Integer categorical class extracted from the
#'   supplied WorldCover-style raster.}
#'   \item{`worldcover_label`}{Human-readable land-cover label corresponding to
#'   `worldcover_class`, where the class is recognised.}
#'   \item{`worldcover_year`}{Year or label supplied through `worldcover_year`.}
#'   \item{`soilgrids_<variable>`}{Continuous raster values extracted from each
#'   named SoilGrids-style raster.}
#'   \item{`human_footprint`}{Continuous Human Footprint-style raster values.}
#' }
#'
#' @section Data licensing and attribution:
#' This function extracts values from user-supplied raster products, including
#' ESA WorldCover, SoilGrids and Human Footprint-style layers. biofetchR does not
#' ship, redistribute or modify these raster datasets. Users are responsible for
#' checking the licence, citation and attribution requirements of each raster
#' source before use in analysis, publication or redistribution.
#'
#' Commonly used sources include ESA WorldCover and SoilGrids, which are commonly
#' distributed under Creative Commons Attribution licences, and Human Footprint
#' products, whose terms vary by dataset and provider. Always cite the specific
#' raster product, version, year and provider used in your workflow.
#'
#' @family raster context helpers
#'
#' @examples
#' if (
#'   requireNamespace("sf", quietly = TRUE) &&
#'   requireNamespace("terra", quietly = TRUE)
#' ) {
#'   occ_sf <- sf::st_as_sf(
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
#'   worldcover_path <- file.path(tempdir(), "worldcover_example.tif")
#'   footprint_path <- file.path(tempdir(), "human_footprint_example.tif")
#'
#'   terra::writeRaster(r, worldcover_path, overwrite = TRUE)
#'   terra::writeRaster(r, footprint_path, overwrite = TRUE)
#'
#'   enriched <- bf_enrich_raster_context_from_sources(
#'     sf_points = occ_sf,
#'     raster_context = c("worldcover", "human_footprint"),
#'     cache_dir = tempdir(),
#'     worldcover_source = worldcover_path,
#'     human_footprint_source = footprint_path,
#'     worldcover_year = 2021
#'   )
#'
#'   names(enriched)
#' }
#'
#' @md
#' @export
bf_enrich_raster_context_from_sources <- function(
    sf_points,
    raster_context = c("worldcover", "soilgrids", "human_footprint"),
    cache_dir = NULL,
    worldcover_source = NULL,
    worldcover_year = 2021,
    worldcover_buffer_m = 0,
    soilgrids_sources = NULL,
    soilgrids_vars = NULL,
    soilgrids_buffer_m = 0,
    human_footprint_source = NULL,
    human_footprint_buffer_m = 0,
    force_refresh = FALSE,
    quiet = TRUE
) {
  if (!inherits(sf_points, "sf")) stop("`sf_points` must be an sf object.", call. = FALSE)
  if (!length(raster_context)) return(sf_points)
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is required for raster enrichment.", call. = FALSE)

  out <- sf_points
  raster_context <- unique(match.arg(raster_context, c("worldcover","soilgrids","human_footprint"), several.ok = TRUE))

  .bf_source_needs_raster_cache <- function(x) {
    if (inherits(x, "SpatRaster")) {
      return(FALSE)
    }

    if (is.null(x) || length(x) == 0L || !nzchar(as.character(x[[1L]]))) {
      return(FALSE)
    }

    src <- as.character(x[[1L]])

    !file.exists(src)
  }

  needs_cache <- FALSE

  if ("worldcover" %in% raster_context) {
    needs_cache <- needs_cache || .bf_source_needs_raster_cache(worldcover_source)
  }

  if ("soilgrids" %in% raster_context && !is.null(soilgrids_sources) && length(soilgrids_sources)) {
    needs_cache <- needs_cache || any(vapply(
      soilgrids_sources,
      .bf_source_needs_raster_cache,
      logical(1)
    ))
  }

  if ("human_footprint" %in% raster_context) {
    needs_cache <- needs_cache || .bf_source_needs_raster_cache(human_footprint_source)
  }

  if (isTRUE(needs_cache)) {
    if (is.null(cache_dir) || length(cache_dir) == 0L ||
        !nzchar(trimws(as.character(cache_dir[[1L]])))) {
      stop(
        "`cache_dir` must be supplied explicitly when URL raster sources are used. ",
        "In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
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

  if ("worldcover" %in% raster_context) {
    r_wc <- .bf_rast_get(worldcover_source, file.path(cache_dir, "worldcover"), stem = paste0("worldcover_", worldcover_year), force_refresh = force_refresh, quiet = quiet)
    vals <- .bf_rast_extract_point_or_buffer(r_wc, out, buffer_m = worldcover_buffer_m, fun = if (isTRUE(worldcover_buffer_m > 0)) terra::modal else NULL, quiet = quiet)
    val_col <- names(vals)[names(vals) != "ID"][1]
    out$worldcover_class <- suppressWarnings(as.integer(vals[[val_col]]))
    out$worldcover_label <- unname(.bf_worldcover_labels[as.character(out$worldcover_class)])
    out$worldcover_year <- as.character(worldcover_year)
  }

  if ("soilgrids" %in% raster_context) {
    if (is.null(soilgrids_sources) || !length(soilgrids_sources)) stop("For SoilGrids enrichment, provide `soilgrids_sources` as a named list/vector of raster paths or URLs.", call. = FALSE)
    if (is.null(soilgrids_vars)) soilgrids_vars <- names(soilgrids_sources)
    if (is.null(soilgrids_vars) || !length(soilgrids_vars)) stop("`soilgrids_vars` could not be inferred. Name the entries in `soilgrids_sources`.", call. = FALSE)
    for (nm in soilgrids_vars) {
      if (is.null(soilgrids_sources[[nm]])) stop("Missing SoilGrids source for variable: ", nm, call. = FALSE)
      r_sg <- .bf_rast_get(soilgrids_sources[[nm]], file.path(cache_dir, "soilgrids"), stem = paste0("soilgrids_", nm), force_refresh = force_refresh, quiet = quiet)
      vals <- .bf_rast_extract_point_or_buffer(r_sg, out, buffer_m = soilgrids_buffer_m, fun = if (isTRUE(soilgrids_buffer_m > 0)) mean else NULL, quiet = quiet)
      val_col <- names(vals)[names(vals) != "ID"][1]
      out[[paste0("soilgrids_", nm)]] <- suppressWarnings(as.numeric(vals[[val_col]]))
    }
  }

  if ("human_footprint" %in% raster_context) {
    r_hf <- .bf_rast_get(human_footprint_source, file.path(cache_dir, "human_footprint"), stem = "human_footprint", force_refresh = force_refresh, quiet = quiet)
    vals <- .bf_rast_extract_point_or_buffer(r_hf, out, buffer_m = human_footprint_buffer_m, fun = if (isTRUE(human_footprint_buffer_m > 0)) mean else NULL, quiet = quiet)
    val_col <- names(vals)[names(vals) != "ID"][1]
    out$human_footprint <- suppressWarnings(as.numeric(vals[[val_col]]))
  }

  out
}
