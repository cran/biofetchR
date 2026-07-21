################################################################################
# clean_and_thin_occurrences.R
# -----------------------------------------------------------------------------
# biofetchR: coordinate cleaning and spatial thinning helpers
# -----------------------------------------------------------------------------
#
# Purpose
#   Provide the package-level helper used by terrestrial/freshwater and marine
#   pipelines to remove low-quality GBIF-style coordinates and, when requested,
#   reduce spatial clustering among occurrence records.
#
# Scope
#   This script only operates on occurrence records already present in memory.
#   It does not submit GBIF downloads, resolve taxonomy, infer native range,
#   query GRIIS, or attach environmental overlays.
#
# Input design
#   The main function is intentionally input-flexible:
#     - package pipelines usually pass sf point objects;
#     - users and tests may pass ordinary data frames/tibbles with
#       decimalLongitude and decimalLatitude columns.
#
# Output design
#   The returned object is always an sf object, with longitude/latitude columns
#   retained for downstream CSV export, diagnostics, joins and tests.
#
################################################################################


#' Clean and spatially thin GBIF-style occurrence points
#'
#' `thin_spatial_points()` performs conservative coordinate screening and
#' optional spatial thinning for GBIF-style occurrence records. It accepts either
#' an [`sf`][sf::sf] point object or a plain data frame/tibble containing
#' `decimalLongitude` and `decimalLatitude` columns.
#'
#' This function is used internally by the main biofetchR pipelines, but it is
#' also exported so users can apply the same cleaning and thinning logic to
#' occurrence tables before, after or outside a full pipeline run.
#'
#' @details
#' The function applies its checks in a fixed order:
#'
#' 1. **Input standardisation**. A plain data frame is screened for usable
#'    longitude/latitude values and converted to an `sf` object in EPSG:4326.
#'    Existing `sf` inputs are retained as `sf` objects.
#' 2. **Coordinate screening**. Records with missing, non-finite, impossible
#'    coordinates, or exact `0,0` coordinates, are removed.
#' 3. **Coordinate-uncertainty handling**. If
#'    `coordinateUncertaintyInMeters` is present, it is coerced to numeric. When
#'    `filter_uncertain = TRUE`, records above `max_uncertainty_m` are removed.
#'    When `filter_uncertain = FALSE`, those records are retained but uncertainty
#'    still acts as a small priority penalty during thinning.
#' 4. **Priority ordering**. Records are sorted before thinning so that higher
#'    values in `priority_col` are preferentially retained. By default,
#'    `individualCount` is used when present.
#' 5. **Spatial thinning**. If `max_density` is supplied, density-based thinning
#'    is applied. Otherwise, if `dist_km > 0`, greedy distance-based thinning is
#'    applied using great-circle distances from [geosphere::distHaversine()].
#'
#' Set `dist_km = 0` and leave `max_density = NULL` to use the function as a
#' cleaning-only helper. In that mode no distance thinning is applied, and all
#' cleaned records are returned.
#'
#' @section Input requirements:
#' `sf_obj` may be either:
#'
#' - an `sf` point object; or
#' - a data frame/tibble with `decimalLongitude` and `decimalLatitude` columns.
#'
#' Non-`sf` inputs are converted to `sf` with EPSG:4326 and with longitude and
#' latitude columns retained (`remove = FALSE`). Existing `sf` inputs are
#' assumed to represent point geometries with longitude/latitude-like
#' coordinates.
#'
#' @section Cleaning rules:
#' The coordinate-cleaning step removes records with:
#'
#' - missing longitude or latitude;
#' - non-finite longitude or latitude;
#' - longitude outside `[-180, 180]`;
#' - latitude outside `[-90, 90]`; or
#' - exact `0,0` coordinates.
#'
#' If `coordinateUncertaintyInMeters` is present and
#' `filter_uncertain = TRUE`, records with uncertainty greater than
#' `max_uncertainty_m` are also removed.
#'
#' @section Thinning behaviour:
#' Distance-based thinning is greedy and order-dependent. The function therefore
#' sorts records before thinning so that higher-priority records are considered
#' first. If `coordinateUncertaintyInMeters` is available, higher uncertainty
#' reduces the thinning priority slightly, even when high-uncertainty records are
#' not filtered out.
#'
#' Density-based thinning is activated by `max_density`. In that mode, the
#' function estimates the input bounding-box area and retains approximately
#' `max_density` records per km^2, after priority sorting.
#'
#' @section Output columns:
#' The returned object is always an `sf` object. The function ensures that
#' `decimalLongitude` and `decimalLatitude` are available in the output. If
#' `return_all = TRUE`, a logical `.thinned` column indicates which rows were
#' retained by the thinning step. If `coordinateUncertaintyInMeters` is present,
#' a logical `.high_uncertainty` column flags retained records above
#' `warn_uncertainty_m`.
#'
#' @section Workflow note:
#' This function is intentionally conservative and local-only. It does not
#' perform taxonomic cleaning, native-range filtering, GRIIS filtering, GBIF
#' download submission or environmental overlay extraction. Those steps belong
#' to higher-level biofetchR pipeline functions.
#'
#' @param sf_obj An `sf` point object, or a data frame/tibble containing
#'   `decimalLongitude` and `decimalLatitude` columns.
#' @param dist_km Numeric. Minimum distance, in kilometres, between retained
#'   points during distance-based thinning. Use `dist_km = 0` for cleaning-only
#'   behaviour.
#' @param priority_col Character. Name of a column used to prioritise retained
#'   records during thinning. Higher values are retained first. If the column is
#'   absent, all points are treated as equal priority. Defaults to
#'   `"individualCount"`.
#' @param return_all Logical. If `TRUE`, return all post-cleaning records and add
#'   `.thinned` to indicate whether each record was retained. If `FALSE`, return
#'   only retained records. Defaults to `FALSE`.
#' @param plot Logical. If `TRUE` and **ggplot2** is installed, print a diagnostic
#'   retained-versus-removed point plot. Defaults to `FALSE`.
#' @param max_density Optional numeric. If supplied, apply density-based thinning
#'   by retaining at most this many points per km^2 within the input bounding box.
#'   If `NULL`, distance-based thinning is used. Defaults to `NULL`.
#' @param max_uncertainty_m Numeric. Maximum accepted coordinate uncertainty, in
#'   metres, when `filter_uncertain = TRUE`. Defaults to `10000`.
#' @param warn_uncertainty_m Numeric. Uncertainty threshold, in metres, above
#'   which retained records are flagged in `.high_uncertainty`. Defaults to
#'   `5000`.
#' @param filter_uncertain Logical. If `TRUE`, remove records with
#'   `coordinateUncertaintyInMeters > max_uncertainty_m`. If `FALSE`, retain
#'   those records but still coerce uncertainty to numeric and use it as a small
#'   priority penalty during thinning. Defaults to `TRUE`.
#' @param quiet Logical. If `TRUE`, suppress console messages. Defaults to
#'   `FALSE`.
#'
#' @return An `sf` point object containing occurrence records that passed the
#'   coordinate-cleaning rules and, when thinning is requested, the spatial
#'   thinning rules. The returned object retains `decimalLongitude` and
#'   `decimalLatitude` columns for downstream CSV export, diagnostics and spatial
#'   joins. If `return_all = TRUE`, the output includes a logical `.thinned`
#'   column, where `TRUE` marks records retained by the thinning step and
#'   `FALSE` marks records removed by thinning. If
#'   `coordinateUncertaintyInMeters` is present, the output may also include a
#'   logical `.high_uncertainty` column indicating records with uncertainty above
#'   `warn_uncertainty_m`. When no valid records remain, an empty `sf` object is
#'   returned with the same output conventions.
#'
#' @examples
#' pts <- data.frame(
#'   species = "Example species",
#'   decimalLongitude = c(-6.26, -6.2601, 0, 181),
#'   decimalLatitude = c(53.35, 53.3501, 0, 53),
#'   coordinateUncertaintyInMeters = c(20, 30, 10, 10),
#'   individualCount = c(5, 1, 1, 1)
#' )
#'
#' cleaned <- thin_spatial_points(
#'   pts,
#'   dist_km = 0,
#'   filter_uncertain = TRUE,
#'   quiet = TRUE
#' )
#'
#' thinned <- thin_spatial_points(
#'   pts,
#'   dist_km = 5,
#'   priority_col = "individualCount",
#'   filter_uncertain = TRUE,
#'   quiet = TRUE
#' )
#'
#' @seealso [sf::st_as_sf()], [sf::st_coordinates()],
#'   [geosphere::distHaversine()]
#'
#' @family occurrence cleaning helpers
#'
#' @export
#' @md

thin_spatial_points <- function(sf_obj,
                                dist_km = 5,
                                priority_col = "individualCount",
                                return_all = FALSE,
                                plot = FALSE,
                                max_density = NULL,
                                max_uncertainty_m = 10000,
                                warn_uncertainty_m = 5000,
                                filter_uncertain = TRUE,
                                quiet = FALSE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    cli::cli_abort("[x] The {.pkg sf} package is required.")
  }

  if (!requireNamespace("geosphere", quietly = TRUE)) {
    cli::cli_abort("[x] The {.pkg geosphere} package is required.")
  }

  # ---------------------------------------------------------------------------
  # Small internal helpers
  # ---------------------------------------------------------------------------

  # Internal nested helper: add a high-uncertainty diagnostic flag when the
  # source uncertainty column is available. This is deliberately scoped inside
  # thin_spatial_points() because it depends on warn_uncertainty_m.
  add_uncertainty_flag <- function(x) {
    if ("coordinateUncertaintyInMeters" %in% names(x) && isTRUE(warn_uncertainty_m > 0)) {
      unc_vals <- suppressWarnings(as.numeric(x$coordinateUncertaintyInMeters))
      x$.high_uncertainty <- !is.na(unc_vals) & unc_vals > warn_uncertainty_m
    }
    x
  }

  # Internal nested helper: add `.thinned` only when the caller has requested
  # all cleaned records rather than only retained records.
  add_thinning_flag <- function(x, value = TRUE) {
    if (isTRUE(return_all)) {
      x$.thinned <- rep(isTRUE(value), nrow(x))
    }
    x
  }

  # Internal nested helper: return an empty sf object with the same diagnostic
  # columns that non-empty outputs would receive.
  return_empty <- function(x, msg = "No input points provided; returning empty result.") {
    if (!isTRUE(quiet)) cli::cli_alert_warning(msg)
    x <- add_thinning_flag(x, value = FALSE)
    x <- add_uncertainty_flag(x)
    x
  }

  # ---------------------------------------------------------------------------
  # Validate scalar settings
  # ---------------------------------------------------------------------------

  dist_km <- suppressWarnings(as.numeric(dist_km)[1])
  if (is.na(dist_km) || dist_km < 0) {
    cli::cli_abort("{.arg dist_km} must be a non-negative numeric value.")
  }

  if (!is.null(max_density)) {
    max_density <- suppressWarnings(as.numeric(max_density)[1])
    if (is.na(max_density) || max_density <= 0) {
      cli::cli_abort("{.arg max_density} must be a positive numeric value when supplied.")
    }
  }

  max_uncertainty_m <- suppressWarnings(as.numeric(max_uncertainty_m)[1])
  if (is.na(max_uncertainty_m) || max_uncertainty_m < 0) {
    cli::cli_abort("{.arg max_uncertainty_m} must be a non-negative numeric value.")
  }

  warn_uncertainty_m <- suppressWarnings(as.numeric(warn_uncertainty_m)[1])
  if (is.na(warn_uncertainty_m) || warn_uncertainty_m < 0) {
    cli::cli_abort("{.arg warn_uncertainty_m} must be a non-negative numeric value.")
  }

  # ---------------------------------------------------------------------------
  # Accept both sf objects and plain GBIF-style data frames
  # ---------------------------------------------------------------------------

  if (is.null(sf_obj)) {
    cli::cli_abort("{.arg sf_obj} is NULL. Supply an sf object or a GBIF-style data frame.")
  }

  if (!inherits(sf_obj, "sf")) {
    if (!is.data.frame(sf_obj)) {
      cli::cli_abort(
        "thin_spatial_points() expects either an sf object or a data frame with {.field decimalLongitude} and {.field decimalLatitude} columns."
      )
    }

    if (!all(c("decimalLongitude", "decimalLatitude") %in% names(sf_obj))) {
      cli::cli_abort(
        "thin_spatial_points() received a non-sf data frame but could not find {.field decimalLongitude} and {.field decimalLatitude} columns."
      )
    }

    sf_obj$decimalLongitude <- suppressWarnings(as.numeric(sf_obj$decimalLongitude))
    sf_obj$decimalLatitude  <- suppressWarnings(as.numeric(sf_obj$decimalLatitude))

    keep_valid_coords <- is.finite(sf_obj$decimalLongitude) &
      is.finite(sf_obj$decimalLatitude) &
      sf_obj$decimalLongitude >= -180 &
      sf_obj$decimalLongitude <= 180 &
      sf_obj$decimalLatitude >= -90 &
      sf_obj$decimalLatitude <= 90 &
      !(sf_obj$decimalLongitude == 0 & sf_obj$decimalLatitude == 0)

    n_removed_coords <- sum(!keep_valid_coords, na.rm = TRUE)

    if (!isTRUE(quiet) && n_removed_coords > 0) {
      cli::cli_alert_info(
        "Filtered {.val {n_removed_coords}} records with missing, invalid, impossible or zero coordinates."
      )
    }

    sf_obj <- sf_obj[keep_valid_coords, , drop = FALSE]

    # `st_as_sf()` can convert a zero-row data frame as long as the coordinate
    # columns exist. Returning sf here keeps the output type stable.
    sf_obj <- sf::st_as_sf(
      sf_obj,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326,
      remove = FALSE
    )
  }

  if (!inherits(sf_obj, "sf")) {
    cli::cli_abort("Internal error: input could not be converted to an sf object.")
  }

  # ---------------------------------------------------------------------------
  # Clean invalid sf geometries/coordinates as well
  # ---------------------------------------------------------------------------

  if (nrow(sf_obj) == 0) {
    return(return_empty(sf_obj))
  }

  coords <- sf::st_coordinates(sf_obj)

  if (!all(c("X", "Y") %in% colnames(coords))) {
    cli::cli_abort("thin_spatial_points() currently expects point geometries with X/Y coordinates.")
  }

  keep_valid_geom <- is.finite(coords[, "X"]) &
    is.finite(coords[, "Y"]) &
    coords[, "X"] >= -180 &
    coords[, "X"] <= 180 &
    coords[, "Y"] >= -90 &
    coords[, "Y"] <= 90 &
    !(coords[, "X"] == 0 & coords[, "Y"] == 0)

  n_removed_geom <- sum(!keep_valid_geom, na.rm = TRUE)

  if (!isTRUE(quiet) && n_removed_geom > 0) {
    cli::cli_alert_info(
      "Filtered {.val {n_removed_geom}} sf record(s) with missing, invalid, impossible or zero coordinates."
    )
  }

  sf_obj <- sf_obj[keep_valid_geom, , drop = FALSE]

  if (nrow(sf_obj) == 0) {
    return(return_empty(sf_obj, "No valid coordinates remain after coordinate cleaning."))
  }

  coords <- sf::st_coordinates(sf_obj)

  # Keep coordinate columns available for downstream CSV exports and diagnostics.
  sf_obj$decimalLongitude <- coords[, "X"]
  sf_obj$decimalLatitude  <- coords[, "Y"]

  # ---------------------------------------------------------------------------
  # Normalise coordinate uncertainty once, regardless of filtering mode
  # ---------------------------------------------------------------------------
  # Even when `filter_uncertain = FALSE`, uncertainty is used below as a small
  # priority penalty. Therefore it must always be numeric before priority sorting.

  if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
    sf_obj$coordinateUncertaintyInMeters <- suppressWarnings(
      as.numeric(sf_obj$coordinateUncertaintyInMeters)
    )
  }

  # ---------------------------------------------------------------------------
  # Remove high-uncertainty records, when requested
  # ---------------------------------------------------------------------------

  if (isTRUE(filter_uncertain) && "coordinateUncertaintyInMeters" %in% names(sf_obj)) {
    unc <- sf_obj$coordinateUncertaintyInMeters
    keep_uncertainty <- is.na(unc) | unc <= max_uncertainty_m
    n_removed_unc <- sum(!keep_uncertainty, na.rm = TRUE)

    if (!isTRUE(quiet) && n_removed_unc > 0) {
      cli::cli_alert_info(
        "Filtered {.val {n_removed_unc}} records with uncertainty > {max_uncertainty_m} m."
      )
    }

    sf_obj <- sf_obj[keep_uncertainty, , drop = FALSE]
  }

  if (nrow(sf_obj) == 0) {
    return(return_empty(sf_obj, "No records remain after coordinate-uncertainty filtering."))
  }

  coords <- sf::st_coordinates(sf_obj)

  # ---------------------------------------------------------------------------
  # If dist_km = 0 and no density limit is requested, this call is being used for
  # cleaning only. Return all cleaned records without accidental duplicate-point
  # thinning.
  # ---------------------------------------------------------------------------

  if (is.null(max_density) && dist_km <= 0) {
    if (!isTRUE(quiet)) {
      cli::cli_alert_info("Distance-based thinning skipped because {.arg dist_km} is 0; returning cleaned records.")
    }

    sf_obj <- add_thinning_flag(sf_obj, value = TRUE)
    sf_obj <- add_uncertainty_flag(sf_obj)
    return(sf_obj)
  }

  n_pts <- nrow(coords)

  if (n_pts < 2) {
    if (!isTRUE(quiet)) cli::cli_alert_info("Only one point provided; skipping thinning.")
    sf_obj <- add_thinning_flag(sf_obj, value = TRUE)
    sf_obj <- add_uncertainty_flag(sf_obj)
    return(sf_obj)
  }

  # ---------------------------------------------------------------------------
  # Sort records by retention priority before thinning
  # ---------------------------------------------------------------------------
  # Distance-based thinning is greedy: the first point retained in a local
  # cluster suppresses nearby points within `dist_km`. Therefore, row order
  # matters. We sort the records before thinning so that the most useful or
  # reliable occurrence in each spatial cluster is considered first.
  #
  # By default, `individualCount` is used as the priority column because records
  # representing more individuals are often more informative than singletons.
  # Missing priority values are treated as the lowest possible priority, so they
  # are only retained when no higher-priority nearby record is available.

  count_vals <- if (!is.null(priority_col) && priority_col %in% names(sf_obj)) {
    val <- suppressWarnings(as.numeric(sf_obj[[priority_col]]))
    ifelse(is.na(val), -Inf, val)
  } else {
    rep(0, nrow(sf_obj))
  }

  # Penalise coordinate uncertainty, if present. This does not remove records;
  # it simply makes high-uncertainty points less likely to be chosen over nearby
  # lower-uncertainty points during thinning. This remains active even when
  # `filter_uncertain = FALSE`, because low-uncertainty records are usually
  # better representatives of a local cluster.
  if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
    unc_vals <- suppressWarnings(as.numeric(sf_obj$coordinateUncertaintyInMeters))
    penalty <- ifelse(is.na(unc_vals), 0, unc_vals / 1000)
    count_vals <- count_vals - penalty
  }

  order_priority <- order(-count_vals)
  coords <- coords[order_priority, , drop = FALSE]
  sf_obj <- sf_obj[order_priority, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # Optional density-based thinning
  # ---------------------------------------------------------------------------

  if (!is.null(max_density)) {
    bbox <- sf::st_bbox(sf_obj)

    area_km2 <- suppressWarnings(
      geosphere::areaPolygon(matrix(c(
        bbox["xmin"], bbox["ymin"],
        bbox["xmin"], bbox["ymax"],
        bbox["xmax"], bbox["ymax"],
        bbox["xmax"], bbox["ymin"],
        bbox["xmin"], bbox["ymin"]
      ), ncol = 2, byrow = TRUE)) / 1e6
    )

    if (!is.finite(area_km2) || area_km2 <= 0) {
      area_km2 <- 0
    }

    target_n <- max(1L, ceiling(area_km2 * max_density))

    if (target_n < nrow(sf_obj)) {
      if (!isTRUE(quiet)) {
        cli::cli_alert_info(
          "Thinning to approximately {.val {target_n}} points based on {.val {max_density}} pts/km^2 over approximately {.val {round(area_km2)}} km^2."
        )
      }

      keep_idx <- seq_len(target_n)

      if (isTRUE(return_all)) {
        sf_obj$.thinned <- FALSE
        sf_obj$.thinned[keep_idx] <- TRUE
        sf_obj <- add_uncertainty_flag(sf_obj)
        return(sf_obj)
      }

      out <- sf_obj[keep_idx, , drop = FALSE]
      out <- add_uncertainty_flag(out)
      return(out)
    }

    if (!isTRUE(quiet)) {
      cli::cli_alert_info("Max density not exceeded; returning all {.val {nrow(sf_obj)}} points.")
    }

    sf_obj <- add_thinning_flag(sf_obj, value = TRUE)
    sf_obj <- add_uncertainty_flag(sf_obj)
    return(sf_obj)
  }

  # ---------------------------------------------------------------------------
  # Distance-based thinning
  # ---------------------------------------------------------------------------

  keep <- logical(nrow(coords))
  remaining <- seq_len(nrow(coords))

  while (length(remaining) > 0) {
    idx <- remaining[1]
    keep[idx] <- TRUE

    dists <- geosphere::distHaversine(
      coords[idx, , drop = FALSE],
      coords[remaining, , drop = FALSE]
    ) / 1000

    remaining <- remaining[dists > dist_km]
  }

  n_retained <- sum(keep)

  if (!isTRUE(quiet)) {
    cli::cli_alert_info(
      "Distance-based thinning retained {.val {n_retained}} of {.val {length(keep)}} points."
    )

    if ("coordinateUncertaintyInMeters" %in% names(sf_obj)) {
      unc_vals <- suppressWarnings(as.numeric(sf_obj$coordinateUncertaintyInMeters))
      high_unc <- !is.na(unc_vals) & unc_vals > warn_uncertainty_m
      n_flagged <- sum(high_unc & keep, na.rm = TRUE)

      if (n_flagged > 0) {
        cli::cli_alert_warning(
          "{.val {n_flagged}} retained points have coordinate uncertainty > {warn_uncertainty_m} m."
        )
      }
    }
  }

  out <- if (isTRUE(return_all)) {
    sf_obj$.thinned <- keep
    sf_obj
  } else {
    sf_obj[keep, , drop = FALSE]
  }

  out <- add_uncertainty_flag(out)

  # ---------------------------------------------------------------------------
  # Optional visualisation
  # ---------------------------------------------------------------------------

  if (isTRUE(plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    coords_df <- as.data.frame(sf::st_coordinates(sf_obj))
    coords_df$.thinned <- if (isTRUE(return_all)) out$.thinned else keep
    coords_df$category <- ifelse(coords_df$.thinned, "Retained", "Removed")

    p <- ggplot2::ggplot(coords_df, ggplot2::aes(X, Y, colour = category)) +
      ggplot2::geom_point(alpha = 0.7, size = 1.5) +
      ggplot2::coord_fixed() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "Spatial thinning results",
        x = "Longitude",
        y = "Latitude"
      )

    print(p)
  }

  out
}
