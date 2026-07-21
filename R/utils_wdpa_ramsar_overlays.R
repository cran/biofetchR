###############################################################################
# utils_wdpa_ramsar_overlays.R
# -----------------------------------------------------------------------------
# WDPA / Ramsar overlay utilities
# -----------------------------------------------------------------------------
# Package-safe helpers for downloading, caching, standardising and loading
# protected-area overlays used by biofetchR spatial workflows.
#
# This script intentionally does not ship WDPA or Ramsar geometries. Instead,
# users opt in to local downloads from the WDPA/Protected Planet ArcGIS service
# and cache country-level extracts on their own machine.
#
# Main exported helpers:
#   bf_load_wdpa()
#   bf_load_ramsar()
#
# Internal helpers:
#   - standardise CRS and geometry validity;
#   - convert ISO2/ISO3 country codes;
#   - safely bind mixed-schema sf objects;
#   - standardise WDPA/Ramsar cache schemas;
#   - query paged ArcGIS REST layers as sf objects.
#
###############################################################################

#' Convert ISO2 country codes to ISO3
#'
#' Thin internal wrapper around `countrycode::countrycode()` used to harmonise
#' user-supplied country codes before querying WDPA/Ramsar country extracts.
#' The helper also handles common special cases used in biodiversity data
#' workflows, including `"UK"` as `"GB"` and Kosovo as `"XKX"`.
#'
#' @param iso2c Character vector of ISO 3166-1 alpha-2 country codes.
#'
#' @return Character vector of ISO 3166-1 alpha-3 country codes.
#'
#' @keywords internal
#' @noRd
bf_iso2_to_iso3 <- function(iso2c) {
  iso2c <- toupper(trimws(as.character(iso2c)))

  if (!requireNamespace("countrycode", quietly = TRUE)) {
    stop("Package 'countrycode' required for ISO conversion.", call. = FALSE)
  }

  iso2c[iso2c == "UK"] <- "GB"

  suppressWarnings(countrycode::countrycode(
    iso2c,
    "iso2c",
    "iso3c",
    custom_match = c("XK" = "XKX", "NA" = "NAM"),
    warn = FALSE
  ))
}

#' Convert ISO3 country codes to ISO2
#'
#' Internal country-code conversion helper used when Ramsar downloads are routed
#' through the WDPA loader, which expects ISO2 country codes at the public
#' function boundary.
#'
#' @param iso3 Character vector of ISO 3166-1 alpha-3 country codes.
#'
#' @return Character vector of ISO 3166-1 alpha-2 country codes.
#'
#' @keywords internal
#' @noRd
bf_iso3_to_iso2 <- function(iso3) {
  iso3 <- toupper(trimws(as.character(iso3)))

  if (!requireNamespace("countrycode", quietly = TRUE)) {
    stop("Package 'countrycode' required for ISO conversion.", call. = FALSE)
  }

  suppressWarnings(countrycode::countrycode(
    iso3,
    "iso3c",
    "iso2c",
    warn = FALSE
  ))
}

#' Bind a list of mixed-schema `sf` objects safely
#'
#' Combines country-level spatial extracts that may contain different attribute
#' columns because they came from different cache versions or service responses.
#' The helper validates geometries, enforces WGS84, standardises the geometry
#' column name to `geometry`, adds missing attribute columns as `NA`, coerces
#' factor columns to character, and then row-binds the objects.
#'
#' @param x_list A list of `sf` objects.
#' @param quiet Logical. If `FALSE`, allow helper-level progress or fallback
#'   messages.
#'
#' @return A single bound `sf` object, or `NULL` when no non-empty `sf` inputs
#'   are available.
#'
#' @keywords internal
#' @noRd
bf_bind_sf_safe <- function(x_list, quiet = TRUE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }

  x_list <- Filter(function(x) inherits(x, "sf") && nrow(x) > 0, x_list)

  if (!length(x_list)) return(NULL)

  x_list <- lapply(x_list, function(x) {
    x <- bf_sf_make_valid(x, quiet = quiet)
    x <- bf_sf_wgs84(x)

    gcol <- attr(x, "sf_column")
    if (is.null(gcol) || !gcol %in% names(x)) {
      stop("sf object is missing a valid geometry column.", call. = FALSE)
    }

    if (gcol != "geometry") {
      names(x)[names(x) == gcol] <- "geometry"
      sf::st_geometry(x) <- "geometry"
    }

    x
  })

  all_cols <- unique(unlist(lapply(x_list, function(x) setdiff(names(x), "geometry"))))

  x_list <- lapply(x_list, function(x) {
    miss <- setdiff(all_cols, setdiff(names(x), "geometry"))
    if (length(miss)) {
      for (nm in miss) x[[nm]] <- NA
    }

    x <- x[, c(all_cols, "geometry"), drop = FALSE]

    for (nm in all_cols) {
      if (is.factor(x[[nm]])) x[[nm]] <- as.character(x[[nm]])
    }

    x
  })

  out <- dplyr::bind_rows(x_list)
  out <- sf::st_as_sf(out, sf_column_name = "geometry")
  sf::st_geometry(out) <- "geometry"
  out
}

#' Create an empty WDPA polygon template
#'
#' Returns a zero-row `sf` object with the standard WDPA columns used throughout
#' biofetchR. This provides a stable return type when a country has no matching
#' WDPA features or a query returns no records.
#'
#' @return A zero-row `sf` object in EPSG:4326 with columns `wdpa_id`,
#'   `wdpa_name`, `wdpa_desig_eng`, `wdpa_iucn_cat`, `wdpa_iso3`,
#'   `wdpa_area_km2`, and `geometry`.
#'
#' @keywords internal
#' @noRd
bf_empty_wdpa_sf <- function() {
  sf::st_sf(
    wdpa_id = character(),
    wdpa_name = character(),
    wdpa_desig_eng = character(),
    wdpa_iucn_cat = character(),
    wdpa_iso3 = character(),
    wdpa_area_km2 = numeric(),
    geometry = sf::st_sfc(crs = 4326)
  )
}

#' Create an empty Ramsar polygon template
#'
#' Returns a zero-row `sf` object with the standard Ramsar columns used
#' throughout biofetchR. This provides a stable return type when a country has
#' no matching Ramsar features or a query returns no records.
#'
#' @return A zero-row `sf` object in EPSG:4326 with columns `ramsar_id`,
#'   `ramsar_name`, `ramsar_iso3`, `ramsar_desig_eng`, and `geometry`.
#'
#' @keywords internal
#' @noRd
bf_empty_ramsar_sf <- function() {
  sf::st_sf(
    ramsar_id = character(),
    ramsar_name = character(),
    ramsar_iso3 = character(),
    ramsar_desig_eng = character(),
    geometry = sf::st_sfc(crs = 4326)
  )
}

#' Standardise a WDPA `sf` object to the biofetchR schema
#'
#' Harmonises WDPA attributes from fresh ArcGIS service responses and older
#' cached files that may use different column names. The function also applies
#' geometry repair and WGS84 standardisation before returning only the columns
#' required by downstream overlay workflows.
#'
#' @param x An `sf` object containing WDPA polygon features.
#' @param quiet Logical. If `FALSE`, allow geometry-repair fallback messages.
#'
#' @return An `sf` object with standard WDPA columns: `wdpa_id`, `wdpa_name`,
#'   `wdpa_desig_eng`, `wdpa_iucn_cat`, `wdpa_iso3`, `wdpa_area_km2`, and
#'   geometry.
#'
#' @keywords internal
#' @noRd
bf_standardize_wdpa_sf <- function(x, quiet = TRUE) {
  if (is.null(x) || !inherits(x, "sf")) {
    stop("Expected an sf object for WDPA standardisation.", call. = FALSE)
  }

  x <- bf_sf_make_valid(x, quiet = quiet)
  x <- bf_sf_wgs84(x)

  nms <- names(x)
  ln  <- tolower(nms)

  pick <- function(cands) {
    for (cand in cands) {
      hit <- which(ln == tolower(cand))
      if (length(hit)) return(nms[hit[1]])
    }
    NULL
  }

  id_col    <- pick(c("wdpa_id", "wdpaid", "site_pid", "site_id", "id", "objectid"))
  name_col  <- pick(c("wdpa_name", "name", "name_eng", "site_name"))
  desig_col <- pick(c("wdpa_desig_eng", "desig_eng"))
  iucn_col  <- pick(c("wdpa_iucn_cat", "iucn_cat"))
  iso3_col  <- pick(c("wdpa_iso3", "iso3", "prnt_iso3"))
  area_col  <- pick(c("wdpa_area_km2", "gis_area", "area_km2"))

  x$wdpa_id        <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))
  x$wdpa_name      <- if (!is.null(name_col)) as.character(x[[name_col]]) else x$wdpa_id
  x$wdpa_desig_eng <- if (!is.null(desig_col)) as.character(x[[desig_col]]) else NA_character_
  x$wdpa_iucn_cat  <- if (!is.null(iucn_col)) as.character(x[[iucn_col]]) else NA_character_
  x$wdpa_iso3      <- if (!is.null(iso3_col)) as.character(x[[iso3_col]]) else NA_character_
  x$wdpa_area_km2  <- if (!is.null(area_col)) suppressWarnings(as.numeric(x[[area_col]])) else NA_real_

  gcol <- attr(x, "sf_column")
  x <- x[, c("wdpa_id", "wdpa_name", "wdpa_desig_eng", "wdpa_iucn_cat", "wdpa_iso3", "wdpa_area_km2", gcol), drop = FALSE]

  if ("wdpa_area_km2" %in% names(x)) {
    x <- x[order(x$wdpa_area_km2, decreasing = FALSE), , drop = FALSE]
  }

  x
}

#' Standardise a Ramsar `sf` object to the biofetchR schema
#'
#' Harmonises Ramsar attributes from fresh extracts and older cached files that
#' may use WDPA-style names. The function also applies geometry repair and WGS84
#' standardisation before returning only the columns required by downstream
#' overlay workflows.
#'
#' @param x An `sf` object containing Ramsar polygon features.
#' @param quiet Logical. If `FALSE`, allow geometry-repair fallback messages.
#'
#' @return An `sf` object with standard Ramsar columns: `ramsar_id`,
#'   `ramsar_name`, `ramsar_iso3`, `ramsar_desig_eng`, and geometry.
#'
#' @keywords internal
#' @noRd
bf_standardize_ramsar_sf <- function(x, quiet = TRUE) {
  if (is.null(x) || !inherits(x, "sf")) {
    stop("Expected an sf object for RAMSAR standardisation.", call. = FALSE)
  }

  x <- bf_sf_make_valid(x, quiet = quiet)
  x <- bf_sf_wgs84(x)

  nms <- names(x)
  ln  <- tolower(nms)

  pick <- function(cands) {
    for (cand in cands) {
      hit <- which(ln == tolower(cand))
      if (length(hit)) return(nms[hit[1]])
    }
    NULL
  }

  id_col    <- pick(c("ramsar_id", "wdpa_id", "site_pid", "id", "objectid"))
  name_col  <- pick(c("ramsar_name", "wdpa_name", "name", "name_eng", "site_name"))
  iso3_col  <- pick(c("ramsar_iso3", "wdpa_iso3", "iso3", "prnt_iso3"))
  desig_col <- pick(c("ramsar_desig_eng", "wdpa_desig_eng", "desig_eng"))

  x$ramsar_id        <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))
  x$ramsar_name      <- if (!is.null(name_col)) as.character(x[[name_col]]) else x$ramsar_id
  x$ramsar_iso3      <- if (!is.null(iso3_col)) as.character(x[[iso3_col]]) else NA_character_
  x$ramsar_desig_eng <- if (!is.null(desig_col)) as.character(x[[desig_col]]) else NA_character_

  gcol <- attr(x, "sf_column")
  x <- x[, c("ramsar_id", "ramsar_name", "ramsar_iso3", "ramsar_desig_eng", gcol), drop = FALSE]

  x
}

#' Query a paged ArcGIS REST layer and return features as `sf`
#'
#' Internal downloader for ArcGIS FeatureServer or MapServer query endpoints. It
#' pages through results using `resultOffset` and `resultRecordCount`, downloads
#' each page as GeoJSON, reads the GeoJSON with `sf`, and safely binds all pages
#' into one WGS84 `sf` object.
#'
#' @param url Character. ArcGIS layer query endpoint URL, usually ending in
#'   `/query`.
#' @param where Character. ArcGIS SQL `WHERE` clause.
#' @param out_fields Character vector of attribute fields to request.
#' @param out_sr Integer or character. Output spatial reference EPSG code.
#'   Defaults to `4326`.
#' @param page_size Integer. Number of features to request per page. ArcGIS
#'   services commonly cap this around 2,000 records.
#' @param quiet Logical. If `FALSE`, report page-level download progress.
#'
#' @return An `sf` object in EPSG:4326. Returns a zero-row `sf` object if no
#'   features are downloaded or readable.
#'
#' @keywords internal
#' @noRd
bf_arcgis_query_sf <- function(
    url,
    where,
    out_fields,
    out_sr = 4326,
    page_size = 2000,
    max_allowable_offset = NULL,
    geometry_precision = NULL,
    quiet = TRUE
) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf required", call. = FALSE)
  }

  old_env <- Sys.getenv("OGR_GEOJSON_MAX_OBJ_SIZE", unset = NA_character_)
  on.exit({
    if (is.na(old_env)) {
      Sys.unsetenv("OGR_GEOJSON_MAX_OBJ_SIZE")
    } else {
      Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = old_env)
    }
  }, add = TRUE)
  Sys.setenv(OGR_GEOJSON_MAX_OBJ_SIZE = "0")

  offset <- 0L
  out <- list()

  repeat {
    q <- list(
      where = where,
      outFields = paste(out_fields, collapse = ","),
      returnGeometry = "true",
      f = "geojson",
      outSR = as.character(out_sr),
      resultOffset = as.character(offset),
      resultRecordCount = as.character(page_size)
    )

    if (!is.null(max_allowable_offset)) {
      q$maxAllowableOffset <- as.character(max_allowable_offset)
    }

    if (!is.null(geometry_precision)) {
      q$geometryPrecision <- as.character(geometry_precision)
    }

    q_str <- paste0(
      names(q), "=",
      vapply(q, utils::URLencode, character(1), reserved = TRUE),
      collapse = "&"
    )

    full <- paste0(url, "?", q_str)
    tmp <- tempfile(fileext = ".geojson")

    ok <- tryCatch({
      if (requireNamespace("curl", quietly = TRUE)) {
        h <- curl::new_handle()

        curl::handle_setopt(
          h,
          post = TRUE,
          postfields = q_str,
          useragent = "biofetchR",
          connecttimeout = 60,
          timeout = max(600, as.numeric(getOption("timeout")))
        )

        curl::handle_setheaders(
          h,
          "Content-Type" = "application/x-www-form-urlencoded"
        )

        res <- curl::curl_fetch_memory(url, handle = h)

        if (!is.null(res$status_code) && res$status_code >= 400) {
          stop("ArcGIS HTTP status ", res$status_code, call. = FALSE)
        }

        writeBin(res$content, tmp)
      } else {
        bf_download_file(
          url = full,
          destfile = tmp,
          quiet = TRUE,
          min_bytes = 1,
          validate_not_html = FALSE
        )
      }

      if (!file.exists(tmp) || bf_file_size(tmp) < 1) {
        stop("ArcGIS response file is empty.", call. = FALSE)
      }

      TRUE
    }, error = function(e) {
      if (!isTRUE(quiet)) {
        message("ArcGIS request failed: ", conditionMessage(e))
        message("ArcGIS GET fallback URL: ", full)
      }
      FALSE
    })

    if (!ok) break

    raw_txt <- tryCatch(
      paste(readLines(tmp, warn = FALSE, n = 40), collapse = "\n"),
      error = function(e) ""
    )

    if (!isTRUE(quiet)) {
      message("ArcGIS response file: ", tmp)
      message("ArcGIS response preview: ", substr(raw_txt, 1, 800))
    }

    # ArcGIS often returns JSON error objects even when f=geojson was requested.
    # Those are not readable by sf, so detect and expose them before returning empty.
    if (grepl('"error"\\s*:', raw_txt) || grepl('"features"\\s*:\\s*\\[\\s*\\]', raw_txt)) {
      if (!isTRUE(quiet)) {
        message("ArcGIS returned an error or empty feature collection.")
        message("ArcGIS URL: ", full)
      }
      try(unlink(tmp), silent = TRUE)
      break
    }

    x <- tryCatch(
      suppressWarnings(sf::st_read(tmp, quiet = TRUE)),
      error = function(e) {
        if (!isTRUE(quiet)) {
          message("sf::st_read() failed for ArcGIS response: ", conditionMessage(e))
          message("ArcGIS URL: ", full)
          message("ArcGIS response preview: ", substr(raw_txt, 1, 800))
        }
        NULL
      }
    )

    try(unlink(tmp), silent = TRUE)

    if (is.null(x) || !nrow(x)) {
      if (!isTRUE(quiet)) {
        message("ArcGIS page produced no readable rows.")
        message("ArcGIS URL: ", full)
      }
      break
    }

    out[[length(out) + 1L]] <- x

    n <- nrow(x)
    offset <- offset + n

    if (!isTRUE(quiet)) {
      message("ArcGIS page fetched: ", n, " (offset now ", offset, ")")
    }

    if (n < page_size) break
  }

  if (!length(out)) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = 4326))[0, ])
  }

  x <- bf_bind_sf_safe(out, quiet = quiet)
  if (is.null(x)) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = 4326))[0, ])
  }

  x <- bf_sf_wgs84(x)
  x
}

#' Load WDPA protected-area polygons for selected countries
#'
#' Downloads, standardises and caches country-level WDPA protected-area polygons
#' from the UNEP-WCMC/Protected Planet licensed ArcGIS service. biofetchR does
#' not redistribute WDPA data; this function creates local user-side cache files
#' that can be reused in later overlay runs.
#'
#' @details
#' The function is deliberately opt-in because WDPA data are licensed. By
#' default, `require_opt_in = TRUE` and users must explicitly set `opt_in = TRUE`
#' after confirming that their intended use complies with WDPA/Protected Planet
#' terms. Cached files are written as one GeoPackage per ISO3 country code.
#'
#' When `exclude_marine = TRUE`, features whose WDPA `realm` is exactly
#' `"Marine"` are excluded from the service query. Mixed or terrestrial records
#' are retained.
#'
#' @param iso2c Character vector of ISO 3166-1 alpha-2 country codes, for
#'   example `"GB"` or `c("GB", "IE")`.
#' @param cache_dir Character. Directory used to read/write cached country
#'   extracts. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, ignore existing cache files and
#'   re-download country extracts.
#' @param quiet Logical. If `FALSE`, print cache/download progress messages.
#' @param exclude_marine Logical. If `TRUE`, exclude WDPA features where
#'   `realm == "Marine"` in the ArcGIS query.
#' @param extra_where Optional character string containing an additional ArcGIS
#'   SQL filter to append to the query with `AND`. Intended for advanced use.
#' @param require_opt_in Logical. If `TRUE`, require `opt_in = TRUE` before any
#'   download is attempted.
#' @param opt_in Logical. Must be `TRUE` when `require_opt_in = TRUE`.
#'
#' @return An `sf` polygon object in WGS84 longitude/latitude coordinates
#'   containing WDPA protected-area features for the requested countries. The
#'   returned object uses the standard biofetchR WDPA schema: `wdpa_id`, a stable
#'   protected-area identifier; `wdpa_name`, the protected-area name;
#'   `wdpa_desig_eng`, the English designation label where available;
#'   `wdpa_iucn_cat`, the IUCN protected-area category where available;
#'   `wdpa_iso3`, the ISO3 country code associated with the feature;
#'   `wdpa_area_km2`, the reported GIS area in square kilometres where
#'   available; and the active geometry column. If no valid country codes are
#'   supplied, or if no features are returned for the requested countries, a
#'   zero-row `sf` object with the same standard columns is returned.
#'
#' @section Data licensing:
#' WDPA/Protected Planet data are subject to external licensing and attribution
#' requirements. Do not redistribute downloaded WDPA geometries through
#' biofetchR, package examples, tests, or derived cache folders.
#'
#' @family overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   wdpa_gb <- bf_load_wdpa(
#'     iso2c = "GB",
#'     cache_dir = file.path(tempdir(), "biofetchR_wdpa"),
#'     opt_in = TRUE,
#'     quiet = FALSE
#'   )
#'
#'   wdpa_gb
#' }
#' }
#'
#' @export
bf_load_wdpa <- function(
    iso2c,
    cache_dir = NULL,
    force_refresh = FALSE,
    quiet = TRUE,
    exclude_marine = TRUE,
    extra_where = NULL,
    require_opt_in = TRUE,
    opt_in = FALSE
) {
  if (isTRUE(require_opt_in) && !isTRUE(opt_in)) {
    stop("WDPA download is opt-in. Set opt_in=TRUE after confirming WDPA terms apply to your use.", call. = FALSE)
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.", call. = FALSE)
  }

  iso2c <- toupper(trimws(as.character(iso2c)))
  iso2c <- iso2c[nzchar(iso2c)]

  if (!length(iso2c)) return(bf_empty_wdpa_sf()[0, ])

  iso3 <- bf_iso2_to_iso3(iso2c)
  bad <- is.na(iso3) | !nzchar(iso3)
  if (any(bad)) {
    stop("Unrecognised ISO2 code(s): ", paste(unique(iso2c[bad]), collapse = ", "), call. = FALSE)
  }
  iso3 <- unique(iso3)

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

  out_list <- vector("list", length(iso3))
  names(out_list) <- iso3

  for (cc in iso3) {
    f_cache <- file.path(cache_dir, paste0("wdpa_poly_", cc, ".gpkg"))

    if (file.exists(f_cache) && !isTRUE(force_refresh)) {
      if (!isTRUE(quiet)) {
        message("WDPA: using cache for ", cc, " -> ", f_cache)
      }

      x <- tryCatch(
        suppressWarnings(sf::st_read(f_cache, quiet = TRUE)),
        error = function(e) NULL
      )

      x <- tryCatch(
        bf_standardize_wdpa_sf(x, quiet = quiet),
        error = function(e) NULL
      )

      if (!is.null(x) && inherits(x, "sf") && nrow(x) > 0) {
        out_list[[cc]] <- x
        next
      } else {
        if (!isTRUE(quiet)) {
          message(
            "WDPA: cache for ", cc,
            " is empty, outdated, unreadable, or corrupt; rebuilding automatically."
          )
        }
        try(unlink(f_cache, force = TRUE), silent = TRUE)
      }
    }

    if (!isTRUE(quiet)) {
      message("WDPA: downloading for ", cc, " (licensed service)")
    }

    where <- paste0(
      "(",
      "prnt_iso3 = '", cc, "'",
      " OR iso3 = '", cc, "'",
      ")"
    )

    if (isTRUE(exclude_marine)) {
      where <- paste0(where, " AND realm IN ('Terrestrial', 'Coastal')")
    }

    if (!is.null(extra_where) && nzchar(extra_where)) {
      where <- paste0(where, " AND (", extra_where, ")")
    }

    poly_url <- "https://data-gis.unep-wcmc.org/server/rest/services/ProtectedSites/wdpa_licensed/MapServer/1/query"

    if (!isTRUE(quiet)) {
      message("WDPA where clause: ", where)
    }

    x <- bf_arcgis_query_sf(
      url = poly_url,
      where = where,
      out_fields = c(
        "site_pid",
        "name",
        "name_eng",
        "desig_eng",
        "iucn_cat",
        "realm",
        "iso3",
        "prnt_iso3",
        "gis_area"
      ),
      out_sr = 4326,
      page_size = 25,
      max_allowable_offset = 0.01,
      geometry_precision = 5,
      quiet = quiet
    )

    if (!nrow(x)) {
      out_list[[cc]] <- bf_empty_wdpa_sf()[0, ]
      try(suppressWarnings(sf::st_write(out_list[[cc]], f_cache, delete_dsn = TRUE, quiet = TRUE)), silent = TRUE)
      next
    }

    x <- bf_standardize_wdpa_sf(x, quiet = quiet)

    try(suppressWarnings(sf::st_write(x, f_cache, delete_dsn = TRUE, quiet = TRUE)), silent = TRUE)
    out_list[[cc]] <- x
  }

  out <- bf_bind_sf_safe(out_list, quiet = quiet)

  if (is.null(out)) {
    out <- bf_empty_wdpa_sf()
  }

  out <- bf_sf_make_valid(out, quiet = quiet)
  out <- bf_sf_wgs84(out)
  out
}

#' Load Ramsar wetland polygons for selected countries
#'
#' Downloads, standardises and caches Ramsar-designated polygon features for one
#' or more countries. Ramsar features are obtained as a designation subset of
#' WDPA, which keeps provenance, licensing, geometry handling and update cadence
#' consistent with the WDPA overlay workflow.
#'
#' @details
#' The function is deliberately opt-in because Ramsar features are queried
#' through the WDPA/Protected Planet service. By default, `require_opt_in = TRUE`
#' and users must explicitly set `opt_in = TRUE` after confirming that their
#' intended use complies with the relevant terms. Cached files are written as one
#' GeoPackage per ISO3 country code.
#'
#' @param iso2c Character vector of ISO 3166-1 alpha-2 country codes.
#' @param cache_dir Character. Directory used to read/write cached country
#'   extracts. Must be supplied explicitly. In examples, tests and vignettes,
#'   use a path under `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, ignore existing cache files and
#'   re-download country extracts.
#' @param quiet Logical. If `FALSE`, print cache/download progress messages.
#' @param require_opt_in Logical. If `TRUE`, require `opt_in = TRUE` before any
#'   download is attempted.
#' @param opt_in Logical. Must be `TRUE` when `require_opt_in = TRUE`.
#'
#' @return An `sf` polygon object in WGS84 longitude/latitude coordinates
#'   containing Ramsar wetland features for the requested countries. The returned
#'   object uses the standard biofetchR Ramsar schema: `ramsar_id`, a stable
#'   Ramsar or WDPA-derived site identifier; `ramsar_name`, the Ramsar site name;
#'   `ramsar_iso3`, the ISO3 country code associated with the feature;
#'   `ramsar_desig_eng`, the English designation label where available; and the
#'   active geometry column. If no valid country codes are supplied, or if no
#'   Ramsar features are returned for the requested countries, a zero-row `sf`
#'   object with the same standard columns is returned.
#'
#' @section Data licensing:
#' Ramsar records returned here are derived from the WDPA/Protected Planet
#' service and remain subject to the relevant external licensing and attribution
#' requirements. Do not redistribute downloaded geometries through biofetchR,
#' package examples, tests, or derived cache folders.
#'
#' @family overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   ramsar_gb <- bf_load_ramsar(
#'     iso2c = "GB",
#'     cache_dir = file.path(tempdir(), "biofetchR_ramsar"),
#'     opt_in = TRUE,
#'     quiet = FALSE
#'   )
#'
#'   ramsar_gb
#' }
#' }
#'
#' @export
bf_load_ramsar <- function(
    iso2c,
    cache_dir = NULL,
    force_refresh = FALSE,
    quiet = TRUE,
    require_opt_in = TRUE,
    opt_in = FALSE
) {
  if (isTRUE(require_opt_in) && !isTRUE(opt_in)) {
    stop("RAMSAR download is opt-in. Set opt_in=TRUE after confirming WDPA/Protected Planet terms apply to your use.", call. = FALSE)
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.", call. = FALSE)
  }

  iso2c <- toupper(trimws(as.character(iso2c)))
  iso2c <- iso2c[nzchar(iso2c)]

  if (!length(iso2c)) return(bf_empty_ramsar_sf()[0, ])

  iso3 <- bf_iso2_to_iso3(iso2c)
  bad <- is.na(iso3) | !nzchar(iso3)
  if (any(bad)) {
    stop("Unrecognised ISO2 code(s): ", paste(unique(iso2c[bad]), collapse = ","), call. = FALSE)
  }
  iso3 <- unique(iso3)

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

  out_list <- vector("list", length(iso3))
  names(out_list) <- iso3

  for (cc in iso3) {
    f_cache <- file.path(cache_dir, paste0("ramsar_poly_", cc, ".gpkg"))

    if (file.exists(f_cache) && !isTRUE(force_refresh)) {
      if (!isTRUE(quiet)) {
        message("RAMSAR: using cache for ", cc, " -> ", f_cache)
      }

      x <- tryCatch(
        suppressWarnings(sf::st_read(f_cache, quiet = TRUE)),
        error = function(e) NULL
      )

      x <- tryCatch(
        bf_standardize_ramsar_sf(x, quiet = quiet),
        error = function(e) NULL
      )

      if (!is.null(x) && inherits(x, "sf") && nrow(x) > 0) {
        out_list[[cc]] <- x
        next
      } else {
        if (!isTRUE(quiet)) {
          message(
            "RAMSAR: cache for ", cc,
            " is empty, outdated, unreadable, or corrupt; rebuilding automatically."
          )
        }
        try(unlink(f_cache, force = TRUE), silent = TRUE)
      }
    }

    if (!isTRUE(quiet)) {
      message("RAMSAR: downloading for ", cc, " (WDPA designation subset)")
    }

    extra <- paste0(
      "(",
      "desig_eng LIKE '%Ramsar%'",
      " OR desig_eng LIKE '%RAMSAR%'",
      " OR name LIKE '%Ramsar%'",
      " OR name_eng LIKE '%Ramsar%'",
      ")"
    )

    wdpa <- bf_load_wdpa(
      iso2c = bf_iso3_to_iso2(cc),
      cache_dir = tempdir(),
      force_refresh = TRUE,
      quiet = quiet,
      exclude_marine = FALSE,
      extra_where = extra,
      require_opt_in = require_opt_in,
      opt_in = opt_in
    )

    if (!nrow(wdpa)) {
      out_list[[cc]] <- bf_empty_ramsar_sf()[0, ]
      try(suppressWarnings(sf::st_write(out_list[[cc]], f_cache, delete_dsn = TRUE, quiet = TRUE)), silent = TRUE)
      next
    }

    wdpa$ramsar_id        <- wdpa$wdpa_id
    wdpa$ramsar_name      <- wdpa$wdpa_name
    wdpa$ramsar_iso3      <- wdpa$wdpa_iso3
    wdpa$ramsar_desig_eng <- wdpa$wdpa_desig_eng

    gcol <- attr(wdpa, "sf_column")
    x <- wdpa[, c("ramsar_id", "ramsar_name", "ramsar_iso3", "ramsar_desig_eng", gcol), drop = FALSE]
    x <- bf_standardize_ramsar_sf(x, quiet = quiet)

    try(suppressWarnings(sf::st_write(x, f_cache, delete_dsn = TRUE, quiet = TRUE)), silent = TRUE)
    out_list[[cc]] <- x
  }

  out <- bf_bind_sf_safe(out_list, quiet = quiet)

  if (is.null(out)) {
    out <- bf_empty_ramsar_sf()
  }

  out <- bf_sf_make_valid(out, quiet = quiet)
  out <- bf_sf_wgs84(out)
  out
}
