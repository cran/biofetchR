################################################################################
# utils_ne_resolve_overlays.R
# -----------------------------------------------------------------------------
# biofetchR: Natural Earth and RESOLVE context-overlay loaders
# -----------------------------------------------------------------------------
#
# PURPOSE
#   Provide lightweight, package-managed overlay loaders for terrestrial context
#   layers used by the retained GBIF processing pipelines:
#
#     - Natural Earth urban-area polygons (`ne_urban`);
#     - Natural Earth Admin-1 state/province polygons (`ne_admin1`);
#     - RESOLVE 2017 terrestrial ecoregions (`resolve2017`).
#
# ROLE IN THE PACKAGE
#   These helpers download, cache, read and standardise external spatial layers
#   into small, predictable `sf` objects in EPSG:4326. They are processing
#   helpers for point-in-polygon context assignment, not plotting functions.
#
# DESIGN PRINCIPLES
#   1. Users should not need to pre-download shapefiles for supported overlays.
#   2. Downloaded archives, extracted shapefiles and GeoPackage subsets are
#      cached on disk to make repeated runs resumable.
#   3. Returned objects expose stable ID/name columns used by downstream joins.
#   4. Country clipping is deliberately lightweight and bbox-based to keep
#      first-run downloads and release-gate tests manageable.
#   5. External data ownership remains with the original providers. biofetchR
#      standardises and caches these layers but does not redistribute them.
#
# DATA SOURCES AND ATTRIBUTION
#   Natural Earth layers are downloaded from Natural Earth public S3 endpoints.
#   RESOLVE 2017 ecoregions are queried from the UNEP-WCMC ArcGIS service. Users
#   should cite the exact source, version/service and access date used in their
#   workflow, and should check provider licence and attribution requirements
#   before publication or redistribution.
#
################################################################################

#' Validate an explicit overlay cache directory
#'
#' @param cache_dir Candidate cache directory.
#' @param context Character label used in error messages.
#'
#' @return Normalised cache directory path.
#'
#' @keywords internal
#' @noRd
bf_require_overlay_cache_dir <- function(cache_dir, context = "this overlay loader") {
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

#' Read the first shapefile found in a directory
#'
#' @param exdir Character. Directory containing extracted shapefile(s).
#' @param quiet Logical. If TRUE, reduce GDAL messages.
#'
#' @return An sf object.
#' @keywords internal
#' @noRd
bf_read_first_shp <- function(exdir, quiet = TRUE) {
  shp <- list.files(exdir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  if (!length(shp)) stop("No .shp found in: ", exdir)
  sf::st_read(shp[1], quiet = isTRUE(quiet))
}

#' Crop an sf overlay to the bounding box of a set of ISO2 countries
#'
#' Uses rnaturalearth to get country polygons. If rnaturalearth is missing, returns overlay unmodified.
#' Cropping is bbox-based and fail-soft: if cropping yields 0 features, returns original overlay.
#'
#' @param overlay sf. Overlay polygons.
#' @param iso2c Character vector. ISO2 country codes.
#' @param tag Character. Label for console messages.
#' @param quiet Logical. If TRUE, reduce messages.
#' @param buffer_deg Numeric. Degrees to expand bbox on each side (default 1).
#'
#' @return sf. Cropped overlay (or original if crop fails).
#' @keywords internal
#' @noRd
bf_crop_to_countries_bbox <- function(overlay, iso2c, tag = "overlay", quiet = TRUE, buffer_deg = 1) {
  if (is.null(overlay) || !inherits(overlay, "sf") || !nrow(overlay)) return(overlay)
  iso2c <- toupper(trimws(as.character(iso2c)))
  iso2c <- iso2c[nzchar(iso2c)]
  if (!length(iso2c)) return(overlay)

  if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(overlay)

  tryCatch({
    world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
    sel <- world[world$iso_a2 %in% iso2c, , drop = FALSE]
    if (!nrow(sel)) return(overlay)

    sel <- sf::st_transform(sel, 4326)
    bb <- sf::st_bbox(sel)
    bb["xmin"] <- bb["xmin"] - buffer_deg
    bb["xmax"] <- bb["xmax"] + buffer_deg
    bb["ymin"] <- bb["ymin"] - buffer_deg
    bb["ymax"] <- bb["ymax"] + buffer_deg

    pre <- nrow(overlay)
    cropped <- suppressWarnings(sf::st_crop(sf::st_transform(overlay, 4326), bb))

    if (!isTRUE(quiet)) message(tag, " crop: ", nrow(cropped), " / ", pre, " features.")
    if (nrow(cropped) > 0) cropped else overlay
  }, error = function(e) overlay)
}

#' Filter Natural Earth Admin-1 polygons to requested ISO2 countries
#'
#' Natural Earth Admin-1 layers contain country-identifying attributes such as
#' `iso_a2`, `iso_3166_2`, `adm0_a3`, `gu_a3`, `admin` or `geonunit`. For
#' Admin-1 overlays, an attribute filter is safer than bbox-only cropping because
#' bbox cropping can silently return the full global layer if geometries are
#' invalid or if a crop operation fails.
#'
#' @param x sf. Natural Earth Admin-1 layer before column reduction.
#' @param iso2c Character vector of ISO2 country codes.
#' @param quiet Logical. If TRUE, reduce messages.
#'
#' @return sf. Filtered Admin-1 layer when matching attributes are found;
#'   otherwise the original layer.
#' @keywords internal
#' @noRd
bf_filter_ne_admin1_to_iso2 <- function(x, iso2c, quiet = TRUE) {
  if (is.null(x) || !inherits(x, "sf") || !nrow(x)) return(x)

  iso2c <- toupper(trimws(as.character(iso2c)))
  iso2c <- unique(iso2c[!is.na(iso2c) & nzchar(iso2c)])
  if (!length(iso2c)) return(x)

  nms <- names(x)
  ln <- tolower(nms)

  keep <- rep(FALSE, nrow(x))

  add_exact_match <- function(candidates, values) {
    values <- unique(toupper(trimws(as.character(values))))
    values <- values[!is.na(values) & nzchar(values)]

    if (!length(values)) return(invisible(NULL))

    for (cand in candidates) {
      hit <- which(ln == tolower(cand))
      if (!length(hit)) next

      vals <- toupper(trimws(as.character(x[[nms[hit[1]]]])))
      keep <<- keep | (!is.na(vals) & vals %in% values)
    }

    invisible(NULL)
  }

  # Direct ISO2 columns, if present.
  add_exact_match(
    candidates = c("iso_a2", "adm0_a2"),
    values = iso2c
  )

  # ISO 3166-2 codes usually look like "GB-ENG", "IE-M", etc.
  hit_iso3166 <- which(ln == "iso_3166_2")
  if (length(hit_iso3166)) {
    vals <- toupper(trimws(as.character(x[[nms[hit_iso3166[1]]]])))
    vals_iso2 <- sub("-.*$", "", vals)
    keep <- keep | (!is.na(vals_iso2) & vals_iso2 %in% iso2c)
  }

  # Build ISO3 and country-name lookup from Natural Earth countries when
  # available. This avoids adding a hard dependency on countrycode here.
  iso3_values <- character()
  country_names <- character()

  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    world <- tryCatch(
      rnaturalearth::ne_countries(scale = "small", returnclass = "sf"),
      error = function(e) NULL
    )

    if (!is.null(world) && inherits(world, "sf") && nrow(world)) {
      w_nms <- names(world)
      w_ln <- tolower(w_nms)

      iso2_col <- w_nms[match("iso_a2", w_ln)]
      if (!is.na(iso2_col)) {
        sel <- world[toupper(as.character(world[[iso2_col]])) %in% iso2c, , drop = FALSE]

        if (nrow(sel)) {
          iso3_cols <- intersect(
            tolower(names(sel)),
            c("adm0_a3", "iso_a3", "sov_a3")
          )

          if (length(iso3_cols)) {
            for (cc in iso3_cols) {
              orig <- names(sel)[match(cc, tolower(names(sel)))]
              iso3_values <- c(iso3_values, as.character(sel[[orig]]))
            }
          }

          name_cols <- intersect(
            tolower(names(sel)),
            c("admin", "name", "name_long", "formal_en", "sovereignt", "geounit")
          )

          if (length(name_cols)) {
            for (cc in name_cols) {
              orig <- names(sel)[match(cc, tolower(names(sel)))]
              country_names <- c(country_names, as.character(sel[[orig]]))
            }
          }
        }
      }
    }
  }

  iso3_values <- unique(toupper(trimws(iso3_values)))
  iso3_values <- iso3_values[!is.na(iso3_values) & nzchar(iso3_values)]

  country_names <- unique(toupper(trimws(country_names)))
  country_names <- country_names[!is.na(country_names) & nzchar(country_names)]

  # Common Natural Earth Admin-1 ISO3/country-name fields.
  add_exact_match(
    candidates = c("adm0_a3", "gu_a3", "sov_a3", "sr_adm0_a3", "iso_a3"),
    values = iso3_values
  )

  add_exact_match(
    candidates = c("admin", "geonunit", "geounit", "sovereignt", "adm0_name", "name_0"),
    values = country_names
  )

  if (!any(keep, na.rm = TRUE)) {
    if (!isTRUE(quiet)) {
      message("NE_ADMIN1 attribute filter found no country matches; returning original layer.")
    }
    return(x)
  }

  out <- x[keep, , drop = FALSE]

  if (!isTRUE(quiet)) {
    message("NE_ADMIN1 attribute filter: ", nrow(out), " / ", nrow(x), " features.")
  }

  out
}

#' Load Natural Earth Urban Areas (auto-download + cache)
#'
#' Creates two columns:
#' - `id_ne_urban`: stable ID per polygon
#' - `geoname_ne_urban`: urban area name (if available)
#'
#' @param cache_dir Character. Folder to store downloaded zips and extracted
#'   shapefiles. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param scale Character. Natural Earth scale, typically "10m".
#' @param force_refresh Logical. If TRUE, re-download and re-extract.
#' @param clip_to_countries Logical. If TRUE, crop to bbox of `iso2c` (if provided).
#' @param iso2c Optional character vector of ISO2 country codes used for cropping.
#' @param quiet Logical. If TRUE, reduce messages.
#' @param url Optional override URL.
#'
#' @return An `sf` polygon object containing Natural Earth urban-area features
#'   in WGS84 longitude/latitude coordinates. The returned object contains
#'   `id_ne_urban`, a character identifier assigned by biofetchR;
#'   `geoname_ne_urban`, a character urban-area name where available from the
#'   source layer or a generated fallback label otherwise; and the active
#'   geometry column. Each row represents one urban-area polygon or multipolygon
#'   used for assigning occurrence records to Natural Earth urban-context
#'   overlays. If `clip_to_countries = TRUE` and `iso2c` is supplied, the object
#'   is cropped to a buffered bounding box around the requested countries.
#'
#' @section Data source and attribution:
#' Downloads Natural Earth urban-area polygons from the Natural Earth public
#' distribution endpoint unless `url` is overridden. biofetchR caches and
#' standardises the layer but does not redistribute or own the data. Cite Natural
#' Earth and the scale used in any downstream analysis.
#'
#' @family Natural Earth and RESOLVE overlay loaders
#' @md
#' @export
bf_load_ne_urban <- function(cache_dir,
                             scale = c("10m"),
                             force_refresh = FALSE,
                             clip_to_countries = TRUE,
                             iso2c = NULL,
                             quiet = TRUE,
                             url = NULL) {

  cache_dir <- bf_require_overlay_cache_dir(
    cache_dir,
    context = "bf_load_ne_urban()"
  )

  scale <- match.arg(scale)
  if (is.null(url)) {
    url <- "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_urban_areas.zip"
  }

  zip_path <- file.path(cache_dir, basename(url))
  exdir <- file.path(cache_dir, sub("\\.zip$", "", basename(url)))

  bf_download_cached(url, zip_path, force_refresh = force_refresh, quiet = quiet)
  bf_unzip_cached(zip_path, exdir, force_refresh = force_refresh, quiet = quiet)

  x <- bf_read_first_shp(exdir, quiet = quiet)
  x <- sf::st_transform(x, 4326)

  # Normalise name field (Natural Earth uses NAME in many layers)
  nm <- intersect(tolower(names(x)), c("name", "name_en", "nam"))
  if (length(nm)) {
    orig <- names(x)[match(nm[1], tolower(names(x)))]
    x$geoname_ne_urban <- as.character(x[[orig]])
  } else {
    x$geoname_ne_urban <- paste0("URBAN_", seq_len(nrow(x)))
  }
  x$id_ne_urban <- paste0("NE_URBAN_", seq_len(nrow(x)))

  x <- x[, c("id_ne_urban", "geoname_ne_urban", attr(x, "sf_column")), drop = FALSE]

  if (isTRUE(clip_to_countries) && !is.null(iso2c)) {
    x <- bf_crop_to_countries_bbox(x, iso2c = iso2c, tag = "NE_URBAN", quiet = quiet, buffer_deg = 1)
  }
  x
}

#' Load Natural Earth Admin-1 (States/Provinces) (auto-download + cache)
#'
#' Creates:
#' - `id_ne_admin1`: stable polygon identifier
#' - `geoname_ne_admin1`: admin-1 name
#'
#' @param cache_dir Character. Folder to store Natural Earth downloads. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param scale Character. "50m" (smaller) or "10m" (more detailed).
#' @param force_refresh Logical. If TRUE, re-download and re-extract.
#' @param clip_to_countries Logical. If TRUE, crop to bbox of `iso2c` (if provided).
#' @param iso2c Optional ISO2 code vector for cropping.
#' @param quiet Logical. If TRUE, reduce messages.
#' @param url Optional override URL.
#'
#' @return An `sf` polygon object containing Natural Earth Admin-1
#'   state/province features in WGS84 longitude/latitude coordinates. The
#'   returned object contains `id_ne_admin1`, a character identifier for each
#'   administrative unit; `geoname_ne_admin1`, a character state/province name
#'   where available from the source layer or a generated fallback label
#'   otherwise; and the active geometry column. Each row represents one
#'   first-level administrative polygon or multipolygon used for assigning
#'   occurrence records to subnational Natural Earth administrative units. If
#'   `clip_to_countries = TRUE` and `iso2c` is supplied, the object is filtered
#'   to the requested countries where country-identifying attributes are
#'   available, with a bounding-box crop used as a fallback.
#'
#' @section Data source and attribution:
#' Downloads Natural Earth Admin-1 state/province polygons from the Natural Earth
#' public distribution endpoint unless `url` is overridden. biofetchR caches and
#' standardises the layer but does not redistribute or own the data. Cite Natural
#' Earth and the scale used in any downstream analysis.
#'
#' @family Natural Earth and RESOLVE overlay loaders
#' @md
#' @export
bf_load_ne_admin1 <- function(cache_dir,
                              scale = c("50m","10m"),
                              force_refresh = FALSE,
                              clip_to_countries = TRUE,
                              iso2c = NULL,
                              quiet = TRUE,
                              url = NULL) {

  cache_dir <- bf_require_overlay_cache_dir(
    cache_dir,
    context = "bf_load_ne_admin1()"
  )

  scale <- match.arg(scale)
  if (is.null(url)) {
    url <- if (scale == "10m") {
      "https://naturalearth.s3.amazonaws.com/10m_cultural/ne_10m_admin_1_states_provinces.zip"
    } else {
      "https://naturalearth.s3.amazonaws.com/50m_cultural/ne_50m_admin_1_states_provinces.zip"
    }
  }

  zip_path <- file.path(cache_dir, basename(url))
  exdir <- file.path(cache_dir, sub("\\.zip$", "", basename(url)))

  bf_download_cached(url, zip_path, force_refresh = force_refresh, quiet = quiet)
  bf_unzip_cached(zip_path, exdir, force_refresh = force_refresh, quiet = quiet)

  x <- bf_read_first_shp(exdir, quiet = quiet)
  x <- sf::st_transform(x, 4326)

  # Admin-1 country filtering must happen before reducing columns. Natural
  # Earth Admin-1 contains useful country attributes such as iso_a2, adm0_a3,
  # iso_3166_2, admin and geonunit, but these would be lost after we keep only
  # id/name/geometry columns.
  if (isTRUE(clip_to_countries) && !is.null(iso2c)) {
    x_pre_filter <- x

    x <- bf_filter_ne_admin1_to_iso2(
      x = x,
      iso2c = iso2c,
      quiet = quiet
    )

    # Fall back to bbox crop only if the attribute filter did not reduce the
    # global Admin-1 layer. This keeps the previous behaviour as a fallback but
    # prevents silent full-layer returns when country attributes are available.
    if (nrow(x) == nrow(x_pre_filter)) {
      x <- bf_crop_to_countries_bbox(
        x,
        iso2c = iso2c,
        tag = "NE_ADMIN1",
        quiet = quiet,
        buffer_deg = 1
      )
    }
  }

  # Name field: many NE layers have `name`
  nm <- intersect(tolower(names(x)), c("name", "name_en", "name_alt", "gn_name"))

  if (length(nm)) {
    orig <- names(x)[match(nm[1], tolower(names(x)))]
    x$geoname_ne_admin1 <- as.character(x[[orig]])
  } else {
    x$geoname_ne_admin1 <- paste0("ADMIN1_", seq_len(nrow(x)))
  }

  if ("adm1_code" %in% tolower(names(x))) {
    orig <- names(x)[match("adm1_code", tolower(names(x)))]
    x$id_ne_admin1 <- as.character(x[[orig]])
  } else if ("ne_id" %in% tolower(names(x))) {
    orig <- names(x)[match("ne_id", tolower(names(x)))]
    x$id_ne_admin1 <- paste0("NEA1_", as.character(x[[orig]]))
  } else {
    x$id_ne_admin1 <- paste0("NEA1_", seq_len(nrow(x)))
  }

  x <- x[, c("id_ne_admin1", "geoname_ne_admin1", attr(x, "sf_column")), drop = FALSE]

  x
}

#' Query an ArcGIS FeatureServer layer and return GeoJSON as sf (bbox-filtered)
#'
#' This is used for RESOLVE 2017. We query only the bbox that covers the requested countries
#' (when available) to keep payload sizes small and avoid GeoJSON size limits.
#'
#' @param query_url Character. Full query endpoint (ending in `/query`).
#' @param bbox Numeric vector length 4: c(xmin, ymin, xmax, ymax) in EPSG:4326.
#' @param out_fields Character. Comma-separated list of fields to return.
#' @param max_allowable_offset Numeric. Optional ArcGIS simplification in degrees (outSR=4326).
#' @param geometry_precision Integer. Optional geometry precision (number of decimal places).
#' @param chunk_size Integer. Number of records per page.
#' @param quiet Logical. Reduce messages.
#'
#' @return sf. Combined features.
#' @keywords internal
#' @noRd
bf_arcgis_query_geojson_bbox <- function(query_url,
                                         bbox,
                                         out_fields,
                                         max_allowable_offset = 0.01,
                                         geometry_precision = 5,
                                         chunk_size = 200,
                                         quiet = TRUE) {

  stopifnot(length(bbox) == 4)
  stopifnot(is.character(query_url), length(query_url) == 1)
  stopifnot(is.character(out_fields), length(out_fields) == 1)

  # Ensure GeoJSON reads won't fail due to GDAL object size limits
  old_lim <- Sys.getenv("OGR_GEOJSON_MAX_OBJ_SIZE", unset = NA_character_)
  on.exit({
    if (is.na(old_lim)) Sys.unsetenv("OGR_GEOJSON_MAX_OBJ_SIZE") else Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = old_lim)
  }, add = TRUE)
  Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = "0")

  # Basic pagination using resultOffset/resultRecordCount
  out <- list()
  offset <- 0L

  repeat {
    qs <- list(
      where = "1=1",
      outFields = out_fields,
      returnGeometry = "true",
      f = "geojson",
      outSR = 4326,
      inSR = 4326,
      geometry = paste(bbox, collapse = ","),
      geometryType = "esriGeometryEnvelope",
      spatialRel = "esriSpatialRelIntersects",
      resultOffset = offset,
      resultRecordCount = chunk_size
    )
    if (!is.null(max_allowable_offset) && is.finite(max_allowable_offset)) qs$maxAllowableOffset <- max_allowable_offset
    if (!is.null(geometry_precision) && is.finite(geometry_precision)) qs$geometryPrecision <- as.integer(geometry_precision)

    resp <- httr::GET(query_url, query = qs)
    httr::stop_for_status(resp)
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")

    # Write to temp file for sf to read
    tf <- tempfile(fileext = ".geojson")
    writeLines(txt, tf, useBytes = TRUE)

    # If no features, break
    js <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
    nfeat <- 0L
    if (!is.null(js) && !is.null(js$features)) nfeat <- nrow(js$features)

    if (!nfeat) break

    sf_part <- suppressWarnings(sf::st_read(tf, quiet = TRUE))
    out[[length(out) + 1L]] <- sf_part

    if (!isTRUE(quiet)) message("RESOLVE2017: fetched ", nfeat, " features (offset ", offset, ").")
    if (nfeat < chunk_size) break
    offset <- offset + chunk_size
  }

  if (!length(out)) return(NULL)

  # Combine safely (columns should match, but be defensive)
  all_cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(x) {
    miss <- setdiff(all_cols, names(x))
    if (length(miss)) for (m in miss) x[[m]] <- NA
    x[, all_cols, drop = FALSE]
  })

  do.call(rbind, out)
}

#' Load RESOLVE Ecoregions (2017) (auto-download + cache)
#'
#' RESOLVE ecoregions provide a modern, consistent terrestrial biogeographic layer that is
#' frequently used in biodiversity and macroecology workflows. We access the UNEP-WCMC
#' ArcGIS service, download a country-bbox subset (by default), cache it as a GeoPackage,
#' and return it as an sf object ready for point-in-polygon tagging.
#'
#' Output columns:
#' - `id_resolve2017` (stable ID)
#' - `geoname_resolve2017` (ecoregion name)
#' - `biome_resolve2017` (biome name if available)
#' - `realm_resolve2017` (realm if available)
#'
#' @param cache_dir Character. Cache directory. Must be supplied explicitly. In
#'   examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical. Re-download even if cached.
#' @param clip_to_countries Logical. If TRUE and iso2c provided, query only bbox around those countries.
#' @param iso2c Optional ISO2 codes for bbox query when clip_to_countries=TRUE.
#' @param quiet Logical. Reduce messages.
#' @param service_url Optional character URL. If `NULL`, the default RESOLVE provider URL is used.
#' @param max_allowable_offset Numeric. ArcGIS geometry simplification in degrees (outSR=4326).
#' @param geometry_precision Integer. ArcGIS geometry precision (decimal places).
#' @param chunk_size Integer. Pagination chunk size.
#'
#' @return An `sf` polygon object containing RESOLVE 2017 terrestrial ecoregion
#'   features in WGS84 longitude/latitude coordinates, or `NULL` when the
#'   provider query returns no readable features for the requested extent. When
#'   features are returned, the object contains `id_resolve2017`, a character
#'   ecoregion identifier; `geoname_resolve2017`, a character ecoregion name;
#'   `biome_resolve2017`, the biome name where available; `realm_resolve2017`,
#'   the realm name where available; and the active geometry column. Each row
#'   represents one ecoregion polygon or multipolygon used for assigning
#'   occurrence records to RESOLVE terrestrial ecoregion context. The returned
#'   subset is also cached locally as a GeoPackage for reuse.
#'
#' @section Data source and attribution:
#' Queries the UNEP-WCMC ArcGIS FeatureServer for the RESOLVE 2017 terrestrial
#' ecoregions layer and caches the returned subset as a GeoPackage. biofetchR
#' does not redistribute or own the service data. Cite the RESOLVE ecoregions
#' data source/service and record the access date and query settings used.
#'
#' @family Natural Earth and RESOLVE overlay loaders
#' @md
#' @export
bf_load_resolve2017 <- function(cache_dir,
                                force_refresh = FALSE,
                                clip_to_countries = TRUE,
                                iso2c = NULL,
                                quiet = TRUE,
                                service_url = NULL,
                                max_allowable_offset = 0.01,
                                geometry_precision = 5,
                                chunk_size = 200) {

  if (is.null(service_url)) {
    service_url <- paste0(
      "https://data-gis.unep-wcmc.org/server/rest/services/",
      "Bio-geographicalRegions/Resolve_Ecoregions/FeatureServer/0"
    )
  }

  if (!requireNamespace("httr", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("bf_load_resolve2017 requires packages: httr, jsonlite")
  }

  cache_dir <- bf_require_overlay_cache_dir(
    cache_dir,
    context = "bf_load_resolve2017()"
  )

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  iso2c_clean <- toupper(trimws(as.character(iso2c)))
  iso2c_clean <- iso2c_clean[nzchar(iso2c_clean)]
  key <- if (isTRUE(clip_to_countries) && length(iso2c_clean)) {
    paste(sort(unique(iso2c_clean)), collapse = "_")
  } else {
    "GLOBAL"
  }

  gpkg_path <- file.path(cache_dir, paste0("resolve2017_", key, ".gpkg"))
  layer_name <- "resolve2017"

  if (file.exists(gpkg_path) && !isTRUE(force_refresh)) {
    x <- suppressWarnings(sf::st_read(gpkg_path, layer = layer_name, quiet = TRUE))
    return(sf::st_transform(x, 4326))
  }

  # Determine bbox (default global)
  bbox <- c(-180, -90, 180, 90)
  if (isTRUE(clip_to_countries) && length(iso2c_clean) && requireNamespace("rnaturalearth", quietly = TRUE)) {
    world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
    sel <- world[world$iso_a2 %in% iso2c_clean, , drop = FALSE]
    if (nrow(sel)) {
      sel <- sf::st_transform(sel, 4326)
      bb <- sf::st_bbox(sel)
      # Light buffer so near-boundary points still match
      bbox <- c(bb["xmin"] - 1, bb["ymin"] - 1, bb["xmax"] + 1, bb["ymax"] + 1)
    }
  }

  query_url <- paste0(sub("/+$", "", service_url), "/query")

  # Request only the bbox subset (or global)
  sf_res <- bf_arcgis_query_geojson_bbox(
    query_url = query_url,
    bbox = bbox,
    out_fields = "eco_name,biome_name,realm,eco_id",
    max_allowable_offset = max_allowable_offset,
    geometry_precision = geometry_precision,
    chunk_size = chunk_size,
    quiet = quiet
  )

  if (is.null(sf_res) || !inherits(sf_res, "sf") || !nrow(sf_res)) {
    if (!isTRUE(quiet)) message("RESOLVE2017: no features returned for bbox; skipping.")
    return(NULL)
  }

  sf_res <- sf::st_transform(sf_res, 4326)

  # Standardise column names (service is lower-case)
  nms <- names(sf_res)
  ln <- tolower(nms)

  get_col <- function(cand) {
    hit <- which(ln == tolower(cand))
    if (length(hit)) nms[hit[1]] else NULL
  }

  c_name <- get_col("eco_name")
  c_id   <- get_col("eco_id")
  c_bio  <- get_col("biome_name")
  c_rlm  <- get_col("realm")

  sf_res$geoname_resolve2017 <- if (!is.null(c_name)) as.character(sf_res[[c_name]]) else paste0("RESOLVE_", seq_len(nrow(sf_res)))
  sf_res$id_resolve2017 <- if (!is.null(c_id)) paste0("RESOLVE_", as.integer(sf_res[[c_id]])) else paste0("RESOLVE_", seq_len(nrow(sf_res)))
  sf_res$biome_resolve2017 <- if (!is.null(c_bio)) as.character(sf_res[[c_bio]]) else NA_character_
  sf_res$realm_resolve2017 <- if (!is.null(c_rlm)) as.character(sf_res[[c_rlm]]) else NA_character_

  gcol <- attr(sf_res, "sf_column")
  keep <- c("id_resolve2017", "geoname_resolve2017", "biome_resolve2017", "realm_resolve2017", gcol)
  keep <- keep[keep %in% names(sf_res)]
  sf_res <- sf_res[, keep, drop = FALSE]

  # Cache to gpkg
  try({
    if (file.exists(gpkg_path)) unlink(gpkg_path, force = TRUE)
    sf::st_write(sf_res, gpkg_path, layer = layer_name, quiet = TRUE)
  }, silent = TRUE)

  if (!isTRUE(quiet)) message("RESOLVE2017 cached: ", gpkg_path)
  sf_res
}
