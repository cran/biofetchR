################################################################################
# util_gadm_joins.R
# -----------------------------------------------------------------------------
# biofetchR: GADM loading and point-to-administrative-unit joins
# -----------------------------------------------------------------------------
#
# PURPOSE
#   This script provides the GADM-specific spatial helper layer used by the
#   terrestrial/freshwater biofetchR workflow. It handles first-run GADM
#   downloads through geodata, normalises schema differences across cached or
#   newly downloaded GADM objects, and attaches stable GADM identifiers and names
#   to GBIF-style point records.
#
# CORE RESPONSIBILITIES
#   1. Convert ISO2 country codes to ISO3 codes used by GADM/geodata.
#   2. Download, cache and reload GADM administrative polygons by level.
#   3. Standardise GADM ID/name columns across full and packed schemas.
#   4. Bind country-level GADM layers into a single sf object where needed.
#   5. Join occurrence points to GADM polygons without duplicating point rows.
#
# DESIGN NOTES
#   - `load_gadm()` returns a named list of per-country sf objects.
#   - `load_all_gadm()` is a convenience wrapper that returns one bound sf layer.
#   - `gadm_join()` performs a one-to-one point assignment by taking the first
#     polygon hit for each point, preventing row multiplication during export.
#   - The implementation is intentionally robust to GADM/geodata schema drift,
#     including packed objects containing only generic `GID` and `GADM_name`
#     fields.
#
# DATA ACCESS AND ATTRIBUTION
#   These helpers can download GADM boundary data through geodata, but this file
#   does not redistribute GADM data with biofetchR. Users should check and follow
#   the GADM and geodata terms, citation requirements and redistribution rules
#   for the exact boundary version used in their workflow.
#
################################################################################

#' Coalesce two character vectors
#'
#' Fill missing or empty values in `x` with corresponding values from `y`.
#'
#' @param x Primary vector.
#' @param y Fallback vector.
#'
#' @return Character vector with missing or blank values replaced where possible.
#'
#' @keywords internal
#' @noRd
.coalesce_chr <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  out <- x
  bad <- is.na(out) | !nzchar(trimws(out))
  out[bad] <- y[bad]
  out
}

#' Pick the first case-insensitive column-name match
#'
#' Searches existing names for the first candidate, ignoring case. This tolerates
#' schema differences across GADM/geodata versions and cache formats.
#'
#' @param nms Character vector of existing column names.
#' @param candidates Character vector of candidate names in priority order.
#'
#' @return Matching column name from `nms`, or `NULL` if no candidate is present.
#' @keywords internal
#' @noRd
.bf_pick_ci <- function(nms, candidates) {
  ln <- tolower(nms)
  for (cand in candidates) {
    hit <- which(ln == tolower(cand))
    if (length(hit)) return(nms[hit[1]])
  }
  NULL
}

#' Find the active geometry column
#'
#' Returns the sf geometry column name using the `sf_column` attribute where
#' possible, with a class-based fallback for older or repaired sf objects.
#'
#' @param x An sf object or data frame-like object.
#'
#' @return Character scalar naming the geometry column, or `"geometry"` as a
#'   fallback.
#' @keywords internal
#' @noRd
.bf_geom_col <- function(x) {
  g <- tryCatch(attr(x, "sf_column", exact = TRUE), error = function(e) NULL)
  if (!is.null(g) && g %in% names(x)) return(g)
  i <- which(vapply(x, function(v) inherits(v, "sfc"), logical(1)))
  if (length(i)) return(names(x)[i[1]])
  "geometry"
}

#' Test whether values look like GADM identifiers
#'
#' Checks whether non-missing values follow the typical GADM ID pattern, such as
#' `GBR.1_1` or `NAM.3.2_1`.
#'
#' @param x Vector of candidate identifiers.
#'
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.bf_gid_like <- function(x) {
  # typical GADM IDs look like "GBR.1_1" / "NAM.3.2_1" etc.
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(FALSE)
  all(grepl("^[A-Z]{3}\\.[^.]+(?:\\.[^.]+)*_\\d+$", x))
}

#' Repair known missing GADM names from stable GADM IDs
#'
#' Some GADM releases/cache formats contain valid administrative-unit IDs but
#' missing human-readable names for particular units. This helper applies a
#' narrow, conservative lookup using stable GADM identifiers. It only fills names
#' that are missing or blank; it does not overwrite names already supplied by
#' GADM.
#'
#' The immediate case is GADM level 1 for Great Britain / United Kingdom, where
#' `GBR.1_1` can appear with a missing name even though it represents England.
#'
#' @param gid Character vector of GADM identifiers.
#' @param name Character vector of GADM names.
#'
#' @return Character vector of repaired names.
#'
#' @keywords internal
#' @noRd
.bf_repair_gadm_names_from_gid <- function(gid, name) {
  gid <- as.character(gid)
  name <- as.character(name)

  manual <- c(
    # GADM level 1, United Kingdom / Great Britain.
    # These are only used when the source name is missing or blank.
    "GBR.1_1" = "England",
    "GBR.2_1" = "Northern Ireland",
    "GBR.3_1" = "Scotland",
    "GBR.4_1" = "Wales"
  )

  bad <- is.na(name) | !nzchar(trimws(name))
  hit <- gid %in% names(manual)

  name[bad & hit] <- unname(manual[gid[bad & hit]])

  name
}

# -----------------------------------------------------------------------------
# Direct GADM country-shapefile fallback
# -----------------------------------------------------------------------------
# This fallback keeps the source as GADM 4.1, but avoids relying only on
# geodata::gadm(). It downloads the official per-country GADM shapefile archive
# and reads the requested administrative level from that archive.
#
# The country archive is level-complete: the same ZIP can contain level 0, 1, 2,
# etc. The helper then selects the requested level, for example:
#   gadm41_GBR_0.shp
#   gadm41_GBR_1.shp
#   gadm41_GBR_2.shp
#
# The returned object is not standardised here. Standardisation remains in
# load_gadm(), so geodata and direct-download objects pass through the same
# schema-cleaning code.

#' Build a direct GADM 4.1 country shapefile ZIP URL
#'
#' @param cc3 ISO3 country code.
#'
#' @return Character URL.
#'
#' @keywords internal
#' @noRd
.bf_gadm_direct_shp_url <- function(cc3) {
  sprintf(
    "https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_%s_shp.zip",
    toupper(as.character(cc3)[[1L]])
  )
}


#' Read one GADM level from the direct per-country shapefile archive
#'
#' @param cc3 ISO3 country code.
#' @param level Integer GADM level.
#' @param cache_dir Cache directory.
#' @param force_refresh Logical. Re-download and re-unpack when TRUE.
#' @param quiet Logical. Suppress messages when TRUE.
#'
#' @return An sf object, or NULL if the fallback fails.
#'
#' @keywords internal
#' @noRd
.bf_gadm_read_direct_country_shp <- function(cc3,
                                             level,
                                             cache_dir,
                                             force_refresh = FALSE,
                                             quiet = FALSE) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required to read GADM shapefiles.", call. = FALSE)
  }

  if (!exists("bf_download_cached", mode = "function", inherits = TRUE)) {
    .bf_msg(
      "[!]  bf_download_cached() is not available; direct GADM fallback skipped.",
      quiet = quiet
    )
    return(NULL)
  }

  cc3 <- toupper(trimws(as.character(cc3)[[1L]]))
  level <- as.integer(level[[1L]])

  if (!nzchar(cc3) || is.na(level)) {
    return(NULL)
  }

  raw_dir <- file.path(cache_dir, "gadm_direct_raw")
  exdir <- file.path(cache_dir, sprintf("gadm41_%s_shp_extracted", cc3))

  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

  zip_url <- .bf_gadm_direct_shp_url(cc3)
  zip_path <- file.path(raw_dir, sprintf("gadm41_%s_shp.zip", cc3))

  .bf_msg(
    "-> trying direct GADM shapefile fallback for ",
    cc3,
    " L",
    level,
    quiet = quiet
  )

  dl <- tryCatch(
    bf_download_cached(
      url = zip_url,
      dest = zip_path,
      force_refresh = force_refresh,
      quiet = quiet,
      min_bytes = 1024,
      validate_not_html = TRUE
    ),
    error = function(e) e
  )

  if (inherits(dl, "condition") || !file.exists(zip_path)) {
    .bf_msg(
      "[!]  direct GADM download failed for ",
      cc3,
      ": ",
      if (inherits(dl, "condition")) conditionMessage(dl) else "file not found after download",
      quiet = quiet
    )
    return(NULL)
  }

  unpacked <- tryCatch(
    {
      if (exists("bf_unpack_archive", mode = "function", inherits = TRUE)) {
        bf_unpack_archive(
          path = zip_path,
          exdir = exdir,
          force_refresh = force_refresh,
          quiet = quiet
        )
      } else {
        utils::unzip(zip_path, exdir = exdir)
      }
      TRUE
    },
    error = function(e) e
  )

  if (inherits(unpacked, "condition")) {
    .bf_msg(
      "[!]  direct GADM unzip failed for ",
      cc3,
      ": ",
      conditionMessage(unpacked),
      quiet = quiet
    )
    return(NULL)
  }

  # Prefer the exact level file.
  exact_pattern <- sprintf("^gadm41_%s_%d\\.shp$", cc3, level)

  shp <- list.files(
    exdir,
    pattern = exact_pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Fallback for small naming variations.
  if (!length(shp)) {
    level_pattern <- sprintf("_%d\\.shp$", level)

    shp <- list.files(
      exdir,
      pattern = level_pattern,
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    # Keep only files that also contain the ISO3 code where possible.
    shp_iso <- shp[grepl(cc3, basename(shp), ignore.case = TRUE)]
    if (length(shp_iso)) {
      shp <- shp_iso
    }
  }

  if (!length(shp)) {
    .bf_msg(
      "[!]  direct GADM archive did not contain a level ",
      level,
      " shapefile for ",
      cc3,
      quiet = quiet
    )
    return(NULL)
  }

  # If several candidates are present, use the shortest basename first. This
  # usually favours the canonical gadm41_ISO_LEVEL.shp file over side-products.
  shp <- shp[order(nchar(basename(shp)), basename(shp))]

  gadm_sf <- tryCatch(
    sf::st_read(shp[[1L]], quiet = isTRUE(quiet)),
    error = function(e) e
  )

  if (inherits(gadm_sf, "condition") || !inherits(gadm_sf, "sf") || !nrow(gadm_sf)) {
    .bf_msg(
      "[!]  direct GADM shapefile read failed for ",
      cc3,
      " L",
      level,
      ": ",
      if (inherits(gadm_sf, "condition")) conditionMessage(gadm_sf) else "no sf rows",
      quiet = quiet
    )
    return(NULL)
  }

  attr(gadm_sf, "biofetchr_gadm_source") <- "direct_gadm_shp"

  gadm_sf
}

# -----------------------------------------------------------------------------
# Harmonise columns across a list of sf/data.frame objects
# -----------------------------------------------------------------------------

#' Harmonise GADM result columns before row binding
#'
#' Prepare a list of GADM-derived `sf` objects or data frames so they can be
#' combined safely by downstream row-binding code.
#'
#' @details
#' GADM loading and joining workflows can return objects with slightly different
#' column sets or incompatible attribute classes across countries, administrative
#' levels or fallback paths. `bf_harmonize_gadm_columns()` drops `NULL` elements,
#' repairs `sf` geometries where possible, transforms spatial inputs to WGS84
#' when needed, adds missing columns to each object, and coerces conflicting
#' attribute classes to safer common representations before binding.
#'
#' The function returns a harmonised **list**, not a single bound object. This is
#' intentional: GADM loaders may still need to call `dplyr::bind_rows()` with an
#' `.id` argument, for example to retain the source ISO2 country code.
#'
#' @param lst List of `sf` objects or data frames produced by GADM loading or
#'   spatial-join helper steps. `NULL` elements are dropped.
#'
#' @return A list of non-`NULL` objects with aligned column sets and safer
#' attribute types for downstream row binding. If `lst` is empty, or contains
#' only `NULL` elements, the empty list is returned.
#'
#' @section Relationship to older helpers:
#' This function replaces the older, ambiguous `harmonize_column_types()` helper
#' for GADM workflows. It is named explicitly to show that it harmonises
#' GADM-derived objects and returns a list, unlike GBIF chunk-binding helpers
#' that return one bound spatial object.
#'
#' @keywords internal
#' @noRd
bf_harmonize_gadm_columns <- function(lst) {
  if (!length(lst)) return(lst)

  # drop NULLs
  lst <- lst[!vapply(lst, is.null, logical(1))]
  if (!length(lst)) return(lst)

  # standardize sf CRS/geometry where present
  lst <- lapply(lst, function(x) {
    if (inherits(x, "sf")) {
      x <- suppressWarnings(sf::st_make_valid(x))
      crs <- tryCatch(sf::st_crs(x), error = function(e) NA)
      if (!is.na(crs) && !is.null(crs) && sf::st_crs(x) != sf::st_crs(4326)) {
        x <- tryCatch(sf::st_transform(x, 4326), error = function(e) x)
      }
    }
    x
  })

  all_cols <- unique(unlist(lapply(lst, names)))
  gcols <- lapply(lst, function(x) if (inherits(x, "sf")) .bf_geom_col(x) else NULL)

  # add missing columns
  lst <- lapply(seq_along(lst), function(i) {
    x <- lst[[i]]
    miss <- setdiff(all_cols, names(x))
    if (length(miss)) {
      for (m in miss) x[[m]] <- NA
    }
    x
  })

  # coerce column types to avoid bind_rows explosions
  # Rule: if any character/factor -> character; else if any numeric/integer -> numeric; else keep as-is.
  for (col in all_cols) {
    classes <- unique(unlist(lapply(lst, function(x) class(x[[col]])[1])))
    if (any(classes %in% c("character", "factor"))) {
      lst <- lapply(lst, function(x) { x[[col]] <- as.character(x[[col]]); x })
    } else if (any(classes %in% c("numeric", "integer", "double"))) {
      lst <- lapply(lst, function(x) { suppressWarnings(x[[col]] <- as.numeric(x[[col]])); x })
    } else if (any(classes %in% c("logical"))) {
      lst <- lapply(lst, function(x) { x[[col]] <- as.logical(x[[col]]); x })
    }
  }

  lst
}

# -----------------------------------------------------------------------------
# load_gadm(): robust loader for a vector of ISO2 codes at a given level
# - Creates cache_dir if missing (first-time users)
# - Downloads via geodata::gadm() when needed
# - Falls back to direct GADM 4.1 per-country shapefile archives when geodata fails
# - Accepts packed schemas containing only GID + GADM_name
# - ALWAYS emits gid_{level} / name_{level} + keeps GID / GADM_name
# -----------------------------------------------------------------------------

#' Load GADM administrative geometries
#'
#' Download, cache and standardise GADM administrative boundary polygons for one
#' or more ISO2 country codes. The function returns a named list of `sf` objects,
#' keyed by ISO2 code, with stable GADM identifier and name columns that can be
#' used by downstream joins.
#'
#' @details
#' `load_gadm()` is designed to be robust across first-run downloads, cached RDS
#' objects and schema differences in GADM/geodata outputs. It converts ISO2 codes
#' to ISO3 codes, downloads boundaries through `geodata::gadm()` when no valid
#' cache exists, and falls back to the official GADM 4.1 per-country shapefile
#' archive when `geodata::gadm()` fails. The resulting object is converted to
#' `sf`, repaired where possible, transformed to WGS84, and normalised to the key
#' columns expected by the terrestrial/freshwater pipeline.
#'
#' The returned layers always include generic `GID` and `GADM_name` columns, plus
#' level-specific `gid_<level>` and `name_<level>` columns. This makes downstream
#' code less sensitive to whether the source object contains full GADM schemas or
#' a compact packed schema.
#'
#' @section Data access and attribution:
#' This function can download GADM boundary data through `geodata`. biofetchR does
#' not ship or redistribute GADM data from this script. Users are responsible for
#' checking the GADM licence, citation and redistribution terms for the boundary
#' version used in their analyses.
#'
#' @param iso2c Character vector of ISO2 country codes, for example `"GB"` or
#'   `"IE"`. The common alias `"UK"` is converted to `"GB"`, and the literal
#'   ISO2 code `"NA"` is treated as Namibia.
#' @param level Integer GADM administrative level to load, usually `0`, `1` or
#'   `2` in biofetchR workflows.
#' @param cache_dir Character. Directory used for GADM downloads and cached
#'   objects. Must be supplied explicitly. In examples, tests and vignettes,
#'   use `file.path(tempdir(), ...)`.
#' @param quiet Logical; suppress progress messages when `TRUE`.
#' @param force_refresh Logical; re-download/rebuild the cleaned cache even when
#'   a cached RDS file exists.
#'
#' @return A named list of `sf` polygon objects, keyed by ISO2 country code.
#'   Each list element contains the successfully loaded GADM polygons for one
#'   country at the requested administrative `level`, transformed to WGS84 where
#'   possible. Returned layers include generic `GID` and `GADM_name` columns plus
#'   level-specific `gid_<level>` and `name_<level>` columns used by downstream
#'   joins. Invalid ISO2 codes or countries that cannot be loaded are skipped,
#'   so the returned list may be shorter than the input `iso2c` vector or empty.
#'   The function is also called for the side effect of downloading, repairing
#'   and caching GADM boundary data locally.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   gadm <- load_gadm(
#'     iso2c = c("GB", "IE"),
#'     level = 1,
#'     cache_dir = file.path(tempdir(), "biofetchR_gadm_cache")
#'   )
#'
#'   names(gadm)
#' }
#' }
#'
#' @family GADM spatial helpers
#' @md
#' @export
load_gadm <- function(iso2c,
                      level = 1,
                      cache_dir = NULL,
                      quiet = FALSE,
                      force_refresh = FALSE) {

  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(requireNamespace("geodata", quietly = TRUE))
  stopifnot(requireNamespace("countrycode", quietly = TRUE))

  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly for GADM downloads/cache use. ",
      "In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
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

  iso2c <- toupper(trimws(as.character(iso2c)))
  iso2c <- iso2c[!is.na(iso2c) & nzchar(iso2c)]
  iso2c[iso2c == "UK"] <- "GB"  # common alias

  # ISO2 -> ISO3 (ensure literal "NA" maps to Namibia)
  iso3c <- countrycode::countrycode(
    iso2c, "iso2c", "iso3c",
    custom_match = c("XK" = "XKX", "NA" = "NAM")
  )

  out <- list()

  for (i in seq_along(iso2c)) {
    cc2 <- iso2c[i]
    cc3 <- iso3c[i]

    if (is.na(cc3) || !nzchar(cc3)) {
      .bf_msg("[!]  ISO2->ISO3 failed for ", cc2, " - skipping", quiet = quiet)
      next
    }

    # We store a cleaned sf cache as RDS to avoid schema drift across runs
    rds <- file.path(cache_dir, sprintf("gadm41_%s_%d_pk.rds", cc3, level))

    gadm_sf <- NULL

    if (file.exists(rds) && !isTRUE(force_refresh)) {
      gadm_sf <- tryCatch(readRDS(rds), error = function(e) NULL)
      if (!is.null(gadm_sf)) .bf_msg("[OK] cached GADM L", level, " for ", cc2, " (", cc3, ")", quiet = quiet)
    }

    if (is.null(gadm_sf)) {
      .bf_msg(
        "[download] downloading GADM L",
        level,
        " for ",
        cc2,
        " (",
        cc3,
        ") via geodata::gadm()",
        quiet = quiet
      )

      gadm_obj <- tryCatch(
        geodata::gadm(country = cc3, level = level, path = cache_dir),
        error = function(e) e
      )

      if (inherits(gadm_obj, "error")) {
        .bf_msg(
          "[!]  geodata::gadm() failed for ",
          cc3,
          " L",
          level,
          ": ",
          conditionMessage(gadm_obj),
          quiet = quiet
        )
        gadm_sf <- NULL
      } else {
        gadm_sf <- tryCatch(
          sf::st_as_sf(gadm_obj),
          error = function(e) NULL
        )

        if (is.null(gadm_sf) || !inherits(gadm_sf, "sf") || nrow(gadm_sf) == 0) {
          .bf_msg(
            "[!]  st_as_sf() produced 0 rows for ",
            cc3,
            " L",
            level,
            " after geodata::gadm(); trying direct GADM fallback.",
            quiet = quiet
          )
          gadm_sf <- NULL
        } else {
          attr(gadm_sf, "biofetchr_gadm_source") <- "geodata"
        }
      }
    }

    if (is.null(gadm_sf)) {
      gadm_sf <- .bf_gadm_read_direct_country_shp(
        cc3 = cc3,
        level = level,
        cache_dir = cache_dir,
        force_refresh = force_refresh,
        quiet = quiet
      )

      if (is.null(gadm_sf) || !inherits(gadm_sf, "sf") || nrow(gadm_sf) == 0) {
        .bf_msg(
          "[!]  both geodata::gadm() and direct GADM fallback failed for ",
          cc3,
          " L",
          level,
          " - skipping",
          quiet = quiet
        )
        next
      }
    }

    # validate + CRS
    gadm_sf <- suppressWarnings(sf::st_make_valid(gadm_sf))
    crs_now <- tryCatch(sf::st_crs(gadm_sf), error = function(e) NA)
    if (!is.na(crs_now) && !is.null(crs_now) && sf::st_crs(gadm_sf) != sf::st_crs(4326)) {
      gadm_sf <- tryCatch(sf::st_transform(gadm_sf, 4326), error = function(e) gadm_sf)
    }

    nms <- names(gadm_sf)

    # Name column selection
    # -------------------------------------------------------------------------
    # Prefer the requested level-specific GADM name column, e.g. NAME_1 for
    # level 1 and NAME_2 for level 2. For direct GADM level-0 shapefiles, the
    # human-readable country name is often stored in COUNTRY rather than NAME_0,
    # while GID_0 stores the ISO-like identifier. Therefore, for level 0, check
    # COUNTRY-style fields before falling back to compact aliases such as
    # GADM_name.
    if (identical(as.integer(level), 0L)) {
      name_candidates <- c(
        "NAME_0",
        "name_0",
        "COUNTRY",
        "country",
        "CNTRY_NAME",
        "cntry_name",
        "ADMIN",
        "admin",
        "NAME_EN",
        "name_en",
        "GADM_name",
        "NAME",
        "name"
      )
    } else {
      name_candidates <- c(
        "GADM_name",
        sprintf("NAME_%d", level),
        sprintf("name_%d", level),
        sprintf("NL_NAME_%d", level),
        sprintf("nl_name_%d", level),
        "NAME",
        "name"
      )
    }

    name_src <- .bf_pick_ci(nms, unique(name_candidates))

    # If still missing, try to find a non-ID character-ish attribute. Avoid
    # geometry and avoid GID-like fields so level-0 direct GADM files do not use
    # GID_0 as the human-readable name when COUNTRY is absent.
    if (is.null(name_src)) {
      gcol <- .bf_geom_col(gadm_sf)

      cand <- setdiff(names(gadm_sf), gcol)

      cand <- cand[
        vapply(
          cand,
          function(z) is.character(gadm_sf[[z]]) || is.factor(gadm_sf[[z]]),
          logical(1)
        )
      ]

      cand <- cand[
        !grepl(
          "^(GID|gid|ID|id)($|_|\\.)",
          cand,
          ignore.case = TRUE
        )
      ]

      if (length(cand)) {
        name_src <- cand[1]
      }
    }

    if (is.null(name_src)) {
      .bf_msg("[!]  No usable name column for ", cc3, " L", level, " - skipping. Cols: ", paste(nms, collapse = ", "), quiet = quiet)
      next
    }

    nm_vals <- as.character(gadm_sf[[name_src]])
    if (!length(nm_vals) || all(is.na(nm_vals) | !nzchar(trimws(nm_vals)))) {
      .bf_msg("[!]  Name column empty for ", cc3, " L", level, " - skipping", quiet = quiet)
      next
    }

    # gid column selection
    gid_src <- .bf_pick_ci(nms, c(
      "GID",
      sprintf("GID_%d", level),
      sprintf("gid_%d", level),
      "gid"
    ))

    gid_vals <- if (!is.null(gid_src)) as.character(gadm_sf[[gid_src]]) else NA_character_
    if (all(is.na(gid_vals) | !nzchar(trimws(gid_vals)))) {
      # fallback: generate stable IDs if absent
      gid_vals <- sprintf("%s_L%d_%05d", cc3, level, seq_len(nrow(gadm_sf)))
    }

    # Repair known missing names using stable GADM IDs before standardising.
    # This is deliberately conservative: it only fills blank/missing names and
    # never overwrites a valid source-supplied GADM name.
    nm_vals <- .bf_repair_gadm_names_from_gid(
      gid = gid_vals,
      name = nm_vals
    )

    # standard names
    gadm_sf$GADM_name <- nm_vals
    gadm_sf$GID       <- gid_vals

    # Ensure level-specific columns exist (what the terrestrial pipeline expects to discover)
    gadm_sf[[sprintf("name_%d", level)]] <- gadm_sf$GADM_name
    gadm_sf[[sprintf("gid_%d",  level)]] <- gadm_sf$GID

    # Clean minimal columns (keep both generic + level-specific)
    gcol <- .bf_geom_col(gadm_sf)
    keep <- c(sprintf("gid_%d", level),
              sprintf("name_%d", level),
              "GID", "GADM_name", gcol)
    gadm_sf <- gadm_sf[, keep, drop = FALSE]

    # cache cleaned object (best effort)
    tryCatch(saveRDS(gadm_sf, rds), error = function(e) NULL)

    out[[cc2]] <- gadm_sf
    .bf_msg("[OK] GADM ready: ", cc2, " L", level, " n=", nrow(gadm_sf),quiet = quiet)
  }

  out <- out[vapply(out, function(x) inherits(x, "sf") && nrow(x) > 0, logical(1))]
  out
}

# -----------------------------------------------------------------------------
# load_all_gadm(): convenience wrapper returning one sf for many countries
# -----------------------------------------------------------------------------

#' Load and bind GADM geometries for multiple countries
#'
#' Convenience wrapper around [load_gadm()] that loads one or more countries and
#' binds the returned list into a single `sf` object. A new `iso2c` column records
#' the source country for each bound polygon.
#'
#' @param iso_codes Character vector of ISO2 country codes.
#' @param level Integer GADM administrative level.
#' @param path Cache directory passed to [load_gadm()] as `cache_dir`. If
#'   `NULL`, [load_gadm()] uses a package cache directory. In examples, tests and
#'   vignettes, supply a path under `tempdir()`.
#' @param quiet Logical; suppress progress messages when `TRUE`.
#' @param force_refresh Logical; re-download/rebuild cached GADM objects when
#'   `TRUE`.
#'
#' @return An `sf` polygon object containing all successfully loaded GADM
#'   countries bound into one layer. The returned object includes an `iso2c`
#'   column identifying the source country for each polygon, along with the GADM
#'   identifier and name columns standardised by [load_gadm()], including `GID`,
#'   `GADM_name`, `gid_<level>` and `name_<level>`. If no country can be loaded,
#'   the function returns an empty `sf` object.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   gadm_l1 <- load_all_gadm(
#'     iso_codes = c("GB", "IE"),
#'     level = 1,
#'     path = file.path(tempdir(), "biofetchR_gadm_cache")
#'   )
#'
#'   names(gadm_l1)
#' }
#' }
#'
#' @family GADM spatial helpers
#' @md
#' @export
load_all_gadm <- function(iso_codes, level = 1, path = NULL, quiet = FALSE, force_refresh = FALSE) {
  lst <- load_gadm(iso2c = iso_codes, level = level, cache_dir = path, quiet = quiet, force_refresh = force_refresh)
  if (!length(lst)) return(sf::st_as_sf(data.frame()))
  lst <- bf_harmonize_gadm_columns(lst)
  dplyr::bind_rows(lst, .id = "iso2c")
}

# -----------------------------------------------------------------------------
# gadm_join(): 1-to-1 join (never duplicates point rows)
# -----------------------------------------------------------------------------

#' Attach GADM labels to occurrence points
#'
#' Join point records to GADM polygons and append one set of GADM identifier/name
#' columns to each point. The join is intentionally one-to-one: if a point
#' intersects multiple polygons, only the first polygon hit is used, which avoids
#' accidental row multiplication during GBIF export.
#'
#' @details
#' `gadm_join()` accepts either the named list returned by [load_gadm()] or a
#' single pre-bound GADM `sf` polygon layer. Both points and polygons are
#' transformed to WGS84 where possible before the intersection test. When no
#' valid polygon layer is available, the output columns are still created and
#' filled with `NA_character_`, keeping the downstream schema stable.
#'
#' @param pts_sf `sf` point object. EPSG:4326 is recommended; other CRSs are
#'   transformed to WGS84 where possible.
#' @param gadm Either a named list from [load_gadm()] or a single `sf` polygon
#'   layer containing GADM geometries.
#' @param level Integer GADM level used to select `gid_<level>` and
#'   `name_<level>` columns.
#' @param id_out Name of the output GADM identifier column. Defaults to
#'   `"id_gadm"`.
#' @param name_out Name of the output GADM name column. Defaults to
#'   `"name_gadm"`.
#' @param geoname_out Name of the output human-readable GADM label column.
#'   Defaults to `"geoname_gadm"`.
#' @param quiet Logical; suppress warning/progress messages when `TRUE`.
#'
#' @return An `sf` point object with the same rows as the input `pts_sf` and with
#'   three GADM assignment columns appended. The column named by `id_out`
#'   contains the matched GADM identifier, the column named by `name_out`
#'   contains the matched administrative-unit name, and the column named by
#'   `geoname_out` contains the human-readable GADM label used by downstream
#'   exports. Points with no polygon match receive `NA_character_` in the added
#'   columns. The function performs a one-to-one assignment and does not
#'   duplicate point rows when polygons overlap.
#'
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   pts <- sf::st_as_sf(
#'     data.frame(
#'       species = "Example species",
#'       decimalLongitude = -5.93,
#'       decimalLatitude = 54.60
#'     ),
#'     coords = c("decimalLongitude", "decimalLatitude"),
#'     crs = 4326,
#'     remove = FALSE
#'   )
#'
#'   poly <- sf::st_sf(
#'     gid_1 = "GBR.1_1",
#'     name_1 = "Example region",
#'     GID = "GBR.1_1",
#'     GADM_name = "Example region",
#'     geometry = sf::st_sfc(
#'       sf::st_polygon(list(rbind(
#'         c(-7, 53),
#'         c(-4, 53),
#'         c(-4, 56),
#'         c(-7, 56),
#'         c(-7, 53)
#'       ))),
#'       crs = 4326
#'     )
#'   )
#'
#'   joined <- gadm_join(pts, poly, level = 1)
#'   names(joined)
#' }
#'
#' @family GADM spatial helpers
#' @md
#' @export
gadm_join <- function(pts_sf, gadm, level = 1,
                      id_out = "id_gadm",
                      name_out = "name_gadm",
                      geoname_out = "geoname_gadm",
                      quiet = TRUE) {

  stopifnot(inherits(pts_sf, "sf"))
  if (!nrow(pts_sf)) {
    pts_sf[[id_out]] <- character()
    pts_sf[[name_out]] <- character()
    pts_sf[[geoname_out]] <- character()
    return(pts_sf)
  }

  if (inherits(gadm, "sf")) {
    poly <- gadm
  } else if (is.list(gadm)) {
    poly <- gadm[!vapply(gadm, is.null, logical(1))]
    if (!length(poly)) {
      pts_sf[[id_out]] <- NA_character_
      pts_sf[[name_out]] <- NA_character_
      pts_sf[[geoname_out]] <- NA_character_
      return(pts_sf)
    }
    poly <- dplyr::bind_rows(bf_harmonize_gadm_columns(poly), .id = "iso2c")
  } else {
    poly <- gadm
  }

  if (is.null(poly) || !inherits(poly, "sf") || !nrow(poly)) {
    pts_sf[[id_out]] <- NA_character_
    pts_sf[[name_out]] <- NA_character_
    pts_sf[[geoname_out]] <- NA_character_
    return(pts_sf)
  }

  # Ensure CRS match
  if (!is.na(sf::st_crs(pts_sf)) && sf::st_crs(pts_sf) != sf::st_crs(4326)) {
    pts_sf <- tryCatch(sf::st_transform(pts_sf, 4326), error = function(e) pts_sf)
  }
  if (!is.na(sf::st_crs(poly)) && sf::st_crs(poly) != sf::st_crs(4326)) {
    poly <- tryCatch(sf::st_transform(poly, 4326), error = function(e) poly)
  }

  id_col   <- sprintf("gid_%d", level)
  name_col <- sprintf("name_%d", level)

  if (!id_col %in% names(poly))   id_col   <- if ("GID" %in% names(poly)) "GID" else id_col
  if (!name_col %in% names(poly)) name_col <- if ("GADM_name" %in% names(poly)) "GADM_name" else name_col

  if (!id_col %in% names(poly) || !name_col %in% names(poly)) {
    .bf_msg("[!] gadm_join(): polygon layer missing expected id/name columns.", quiet = quiet)
    pts_sf[[id_out]] <- NA_character_
    pts_sf[[name_out]] <- NA_character_
    pts_sf[[geoname_out]] <- NA_character_
    return(pts_sf)
  }

  vals_id   <- as.character(poly[[id_col]])
  vals_name <- as.character(poly[[name_col]])

  # 1-to-1 assignment: pick first hit only
  idx <- tryCatch(suppressWarnings(sf::st_intersects(pts_sf, poly)), error = function(e) NULL)
  if (is.null(idx)) {
    pts_sf[[id_out]] <- NA_character_
    pts_sf[[name_out]] <- NA_character_
    pts_sf[[geoname_out]] <- NA_character_
    return(pts_sf)
  }

  out_id <- rep(NA_character_, length(idx))
  out_nm <- rep(NA_character_, length(idx))

  for (i in seq_along(idx)) {
    h <- idx[[i]]
    if (length(h)) {
      out_id[i] <- vals_id[h[1]]
      out_nm[i] <- vals_name[h[1]]
    }
  }

  pts_sf[[id_out]]      <- out_id
  pts_sf[[name_out]]    <- out_nm
  pts_sf[[geoname_out]] <- out_nm

  pts_sf
}
