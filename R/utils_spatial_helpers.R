################################################################################
# utils_spatial_helpers.R
# ------------------------------------------------------------------------------
# biofetchR: small spatial utility helpers used across overlay and join workflows
# ------------------------------------------------------------------------------
#
# PURPOSE
#   This script collects lightweight `sf` helper functions that are reused by
#   the package's terrestrial, freshwater, marine and contextual-overlay
#   workflows.
#
#   The helpers here intentionally stay small and dependency-light. They handle
#   common spatial housekeeping tasks that otherwise become repeated across
#   loaders and pipelines:
#
#     - standardising spatial layers to WGS84 / EPSG:4326;
#     - repairing invalid geometries before spatial joins;
#     - finding and standardising the active geometry column in `sf` objects;
#     - converting simple numeric bounding boxes into `sf` polygons.
#
# DESIGN NOTES
#   These functions are defensive by design. Several optional overlay loaders can
#   return `NULL`, empty objects, or non-`sf` objects when a provider download is
#   unavailable or when an overlay is skipped. Helpers therefore return such
#   objects unchanged where possible, rather than failing early and breaking a
#   larger pipeline.
#
################################################################################


#' Transform or assign an `sf` object to WGS84
#'
#' Standardises spatial objects to EPSG:4326 for downstream spatial joins and
#' overlay operations. If the input has no coordinate reference system, WGS84 is
#' assigned rather than transformed. If the input already uses EPSG:4326, it is
#' returned unchanged. Empty, `NULL`, or non-`sf` inputs are returned unchanged so
#' optional overlay loaders can fail softly.
#'
#' @param x An `sf` object, or an object that may be `NULL`/empty during optional
#'   overlay workflows.
#'
#' @return The input object with CRS EPSG:4326 where possible.
#'
#' @keywords internal
#' @noRd
bf_sf_wgs84 <- function(x) {
  if (is.null(x) || !inherits(x, "sf") || !nrow(x)) {
    return(x)
  }

  cr <- suppressWarnings(sf::st_crs(x))

  if (is.na(cr)) {
    x <- sf::st_set_crs(x, 4326)
  } else if (!identical(sf::st_crs(x), sf::st_crs(4326))) {
    x <- sf::st_transform(x, 4326)
  }

  x
}

#' Repair invalid `sf` geometries using the safest available method
#'
#' Applies `sf::st_make_valid()` when available. On older `sf` installations,
#' falls back to the common `sf::st_buffer(x, 0)` repair approach. Empty,
#' `NULL`, or non-`sf` inputs are returned unchanged so optional overlay loaders
#' can fail softly.
#'
#' @param x An `sf` object, or an object that may be `NULL`/empty during optional
#'   overlay workflows.
#' @param quiet Logical. If `FALSE`, report when the buffer-based fallback is
#'   used.
#'
#' @return The input object with best-effort valid geometries where possible.
#'
#' @keywords internal
#' @noRd
bf_sf_make_valid <- function(x, quiet = TRUE) {
  if (is.null(x) || !inherits(x, "sf") || !nrow(x)) {
    return(x)
  }

  if ("st_make_valid" %in% getNamespaceExports("sf")) {
    return(suppressWarnings(sf::st_make_valid(x)))
  }

  if (!isTRUE(quiet)) {
    message("sf::st_make_valid not available; using sf::st_buffer(., 0) fallback")
  }

  suppressWarnings(sf::st_buffer(x, 0))
}

#' Return the active geometry column name for an sf object
#'
#' Uses the `sf_column` attribute where possible and falls back to locating the
#' first `sfc` column by class. This avoids relying on `sf::st_geometry_name()`,
#' which improves compatibility with older sf versions.
#'
#' @param x An `sf` object.
#'
#' @return Character scalar geometry-column name, or `NULL` when no geometry
#'   column can be identified.
#'
#' @keywords internal
#' @noRd
bf_sf_geometry_col <- function(x) {
  if (!inherits(x, "sf")) {
    return(NULL)
  }

  gcol <- attr(x, "sf_column", exact = TRUE)

  if (!is.null(gcol) && gcol %in% names(x)) {
    return(gcol)
  }

  sfc_cols <- names(x)[vapply(x, function(col) inherits(col, "sfc"), logical(1))]

  if (length(sfc_cols)) {
    return(sfc_cols[[1L]])
  }

  NULL
}


#' Standardise an sf geometry column name
#'
#' Renames the active geometry column to a predictable name, usually
#' `"geometry"`, and resets the `sf_column` attribute. This is useful after
#' reading shapefiles, GeoPackages or service outputs that use provider-specific
#' geometry-column names.
#'
#' @param x An `sf` object.
#' @param geometry_name Character. Desired geometry column name.
#'
#' @return `x` with a standardised active geometry column.
#'
#' @keywords internal
#' @noRd
bf_sf_standardise_geometry <- function(x, geometry_name = "geometry") {
  if (!inherits(x, "sf")) {
    return(x)
  }

  gcol <- bf_sf_geometry_col(x)

  if (is.null(gcol)) {
    stop("No geometry column found in sf object.", call. = FALSE)
  }

  if (!identical(gcol, geometry_name)) {
    names(x)[names(x) == gcol] <- geometry_name
  }

  attr(x, "sf_column") <- geometry_name
  sf::st_geometry(x) <- geometry_name

  x
}


#' Convert a numeric bounding box to an sf polygon
#'
#' Builds a WGS84 polygon from a numeric bounding box supplied as
#' `c(xmin, ymin, xmax, ymax)`.
#'
#' @param bbox Numeric vector of length four: `xmin`, `ymin`, `xmax`, `ymax`.
#' @param crs Coordinate reference system for the output polygon. Defaults to
#'   EPSG:4326.
#'
#' @return An `sf` polygon object.
#'
#' @keywords internal
#' @noRd
bf_bbox_to_sf <- function(bbox, crs = 4326) {
  stopifnot(is.numeric(bbox), length(bbox) == 4L)

  coords <- matrix(
    c(
      bbox[[1]], bbox[[2]],
      bbox[[3]], bbox[[2]],
      bbox[[3]], bbox[[4]],
      bbox[[1]], bbox[[4]],
      bbox[[1]], bbox[[2]]
    ),
    ncol = 2,
    byrow = TRUE
  )

  sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = crs)
  )
}
