###############################################################################
# utils_spatial_joins.R
# -----------------------------------------------------------------------------
# biofetchR spatial join and context-overlay utilities
# -----------------------------------------------------------------------------
# Helpers for joining GBIF occurrence points to marine, freshwater, terrestrial
# and administrative spatial layers used by biofetchR workflows.
#
# Main exported helpers include:
#   - eez_join()
#   - overlay_join()
#   - bf_load_teow()
#   - bf_teow_cache_info()
#   - bf_teow_clear_cache()
#   - bf_load_feow()
#   - bf_load_lakes()
#   - bf_load_rivers()
#   - bf_name_rivers_osm()
#   - bf_load_basins()
#   - bf_cache_dir()
#   - bf_load_ne_urban_areas()
#   - bf_load_ne_admin1()
#   - bf_load_resolve_ecoregions2017()
#
#   biofetchR downloads and caches third-party spatial layers for user-side
#   analysis only. Users remain responsible for checking each data provider's
#   licence, citation and redistribution requirements.
#
###############################################################################

#' Validate an explicitly supplied spatial cache directory
#'
#' @param cache_dir Candidate cache directory.
#' @param context Character label used in error messages.
#' @param create Logical; if `TRUE`, create the directory if missing.
#'
#' @return Normalised cache directory path.
#'
#' @keywords internal
#' @noRd
.bf_require_explicit_spatial_cache <- function(cache_dir,
                                               context = "this helper",
                                               create = TRUE) {
  if (is.null(cache_dir) || length(cache_dir) == 0L ||
      !nzchar(trimws(as.character(cache_dir[[1L]])))) {
    stop(
      "`cache_dir` must be supplied explicitly for ", context,
      ". In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  cache_dir <- normalizePath(
    as.character(cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (isTRUE(create) && !dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  cache_dir
}

#' Join GBIF points to Exclusive Economic Zone polygons
#'
#' Loads or caches Exclusive Economic Zone (EEZ) polygons via \pkg{mregions2},
#' repairs invalid geometries where possible, harmonises common EEZ attribute
#' names across Marine Regions releases, and joins each occurrence-point
#' \code{sf} object in \code{results_list} to an EEZ polygon.
#'
#' The function returns the same named list structure as the input and appends
#' standardised EEZ metadata used by downstream native-status and marine-overlay
#' workflows. When a point does not fall strictly within an EEZ polygon, the
#' function attempts a nearest-feature fallback so that points lying close to
#' polygon boundaries can still be assigned cautiously.
#'
#' @param results_list Named list of point \code{sf} objects. Each list element is
#'   usually one species or one species-country download result. A single
#'   \code{sf} object is also accepted and will be wrapped into a length-one list.
#' @param use_planar Logical. If \code{TRUE}, disable spherical geometry
#'   operations with \code{sf::sf_use_s2(FALSE)} and run joins in a projected CRS.
#'   If \code{FALSE}, the function first tries the layer CRS and falls back to
#'   planar joins when spherical operations fail.
#' @param eez_cache_file Path to a GeoPackage used to cache the EEZ layer. Must
#'   be supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#'
#' @return A named list with the same names as `results_list`. Each element is
#'   an `sf` point object containing the original occurrence records plus
#'   standardised Exclusive Economic Zone attributes where a spatial match could
#'   be made. Added columns include `geoname`, `mrgid`, `territory1`,
#'   `sovereign1`, `iso3_best` and `iso2_best`. These columns identify the
#'   matched marine region and provide best-effort country codes for downstream
#'   filtering, summarisation and export. Records that cannot be assigned retain
#'   the input geometry and receive missing values in the added fields.
#'
#' @details
#' EEZ source attributes differ across Marine Regions releases. This helper uses
#' case-insensitive matching to find common territory, sovereign and ISO-code
#' fields, then derives best-effort ISO3 and ISO2 codes for downstream matching.
#' If \pkg{countrycode} is installed, missing ISO codes are inferred from cleaned
#' territory, sovereign or EEZ names where possible.
#'
#' @section Data access and attribution:
#' EEZ polygons are retrieved from Marine Regions through \pkg{mregions2} and
#' cached locally. biofetchR does not redistribute Marine Regions data. Users
#' should cite Marine Regions and check the applicable data licence and
#' attribution requirements for the version used in their analysis.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   pts <- data.frame(
#'     decimalLongitude = c(-5, 10),
#'     decimalLatitude = c(50, 35)
#'   )
#'
#'   pts_sf <- sf::st_as_sf(
#'     pts,
#'     coords = c("decimalLongitude", "decimalLatitude"),
#'     crs = 4326,
#'     remove = FALSE
#'   )
#'
#'   joined <- eez_join(
#'     list("Carcinus maenas" = pts_sf),
#'     eez_cache_file = file.path(tempdir(), "biofetchR_eez_cache.gpkg")
#'   )
#'
#'   names(joined[["Carcinus maenas"]])
#' }
#' }
#'
#' @family marine overlay joins
#' @export
eez_join <- function(results_list, use_planar = TRUE, eez_cache_file = NULL) {
  if (!requireNamespace("mregions2", quietly = TRUE)) {
    stop("Package 'mregions2' is required for this operation. Please install it before using this function.", call. = FALSE)
  }

  if (is.null(eez_cache_file) || length(eez_cache_file) == 0L ||
      !nzchar(trimws(as.character(eez_cache_file[[1L]])))) {
    stop(
      "`eez_cache_file` must be supplied explicitly. In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
      call. = FALSE
    )
  }

  eez_cache_file <- normalizePath(
    as.character(eez_cache_file[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (!dir.exists(dirname(eez_cache_file))) {
    dir.create(dirname(eez_cache_file), recursive = TRUE, showWarnings = FALSE)
  }

  # load or download
  if (file.exists(eez_cache_file)) {
    message("Using cached EEZ polygons: ", eez_cache_file)
    eez_sf <- sf::read_sf(eez_cache_file)
  } else {
    message("Downloading EEZ polygons via mregions2::mrp_get('eez')...")
    eez_sf <- mregions2::mrp_get("eez")
    sf::write_sf(eez_sf, eez_cache_file)
  }

  # fix invalid geoms (prefer lwgeom if present; else buffer(0) trick)
  eez_sf <- suppressWarnings(sf::st_make_valid(eez_sf))
  bad <- !sf::st_is_valid(eez_sf)
  if (any(bad)) {
    if (requireNamespace("lwgeom", quietly = TRUE)) {
      eez_sf[bad, ] <- sf::st_make_valid(eez_sf[bad, ])
    } else {
      eez_sf <- sf::st_buffer(eez_sf, 0)
    }
  }

  # helper: case-insensitive picker + safe column getter (no warnings)
  pick_ci <- function(x, candidates) {
    nx <- names(x); for (cand in candidates) {
      hit <- which(tolower(nx) == tolower(cand)); if (length(hit)) return(nx[hit[1]])
    }; NULL
  }
  gcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL

  # normalize handles
  geoname_col   <- pick_ci(eez_sf, c("geoname","name","geoname_lng","geoname_long"))
  mrgid_col     <- pick_ci(eez_sf, c("mrgid","MRGID"))
  terr1_col     <- pick_ci(eez_sf, c("TERRITORY1","territory1","territory"))
  sov1_col      <- pick_ci(eez_sf, c("SOVEREIGN1","sovereign1","sovereign"))
  iso_ter1_col  <- pick_ci(eez_sf, c("ISO_Ter1","iso_ter1"))
  iso_sov1_col  <- pick_ci(eez_sf, c("ISO_SOV1","iso_sov1"))
  iso2_ter1_col <- pick_ci(eez_sf, c("ISO2_Ter1","iso2_ter1"))
  iso2_sov1_col <- pick_ci(eez_sf, c("ISO2_SOV1","iso2_sov1"))

  keep_cols <- unique(na.omit(c(geoname_col, mrgid_col,
                                terr1_col,  sov1_col,
                                iso_ter1_col, iso_sov1_col,
                                iso2_ter1_col, iso2_sov1_col)))
  eez_keep <- eez_sf[, keep_cols, drop = FALSE]

  # standardize names
  ren <- function(old, new) if (!is.null(old) && old %in% names(eez_keep)) names(eez_keep)[names(eez_keep)==old] <<- new
  ren(geoname_col, "geoname"); ren(mrgid_col, "mrgid"); ren(terr1_col, "territory1")
  ren(sov1_col, "sovereign1"); ren(iso_ter1_col, "iso_ter1"); ren(iso_sov1_col, "iso_sov1")
  ren(iso2_ter1_col, "iso2_ter1"); ren(iso2_sov1_col, "iso2_sov1")
  if (!"geoname" %in% names(eez_keep) && "name" %in% names(eez_keep)) eez_keep$geoname <- eez_keep$name

  # compute iso3_best / iso2_best without warning on missing cols
  coalesce_char <- function(...) {
    args <- list(...); n <- nrow(eez_keep); out <- rep(NA_character_, n)
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v); fill <- is.na(out) | !nzchar(out); out[fill] <- v[fill]
    }; out
  }
  iso3_direct <- coalesce_char(gcol(eez_keep, "iso_ter1"), gcol(eez_keep, "iso_sov1"))
  if (requireNamespace("countrycode", quietly = TRUE)) {
    need <- is.na(iso3_direct) | !nzchar(iso3_direct)
    if (any(need)) {
      nm_guess <- coalesce_char(gcol(eez_keep, "territory1"), gcol(eez_keep, "sovereign1"), gcol(eez_keep, "geoname"))
      nm_guess <- gsub("Exclusive Economic Zone|Joint regime area|Overlapping claim area|\\(.*?\\)|-{1,2}\\s*EEZ|\\s+EEZ", "", nm_guess, ignore.case = TRUE)
      nm_guess <- trimws(nm_guess)
      iso3_from_name <- suppressWarnings(countrycode::countrycode(nm_guess, "country.name", "iso3c"))
      iso3_direct[need] <- toupper(iso3_from_name[need])
    }
  }
  eez_keep$iso3_best <- toupper(iso3_direct)

  iso2_direct <- coalesce_char(gcol(eez_keep, "iso2_ter1"), gcol(eez_keep, "iso2_sov1"))
  if ((is.null(iso2_direct) || all(is.na(iso2_direct))) && requireNamespace("countrycode", quietly = TRUE)) {
    iso2_direct <- suppressWarnings(countrycode::countrycode(eez_keep$iso3_best, "iso3c", "iso2c"))
  }
  eez_keep$iso2_best <- toupper(iso2_direct)

  # s2 mode handling + auto-fallback on error
  old_s2 <- sf::sf_use_s2()
  on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)
  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE)

  joined_list <- list()
  for (key in names(results_list)) {
    gbif_sf <- results_list[[key]]
    if (is.null(gbif_sf) || !nrow(gbif_sf)) {
      gbif_sf$geoname <- NA_character_
      joined_list[[key]] <- gbif_sf
      next
    }

    # both layers in the same CRS; try s2 if enabled, else planar; fallback if needed
    try_join <- function(force_planar = FALSE) {
      if (force_planar) sf::sf_use_s2(FALSE)
      crs_target <- if (!force_planar && !isTRUE(use_planar)) sf::st_crs(eez_keep) else 3857
      gb  <- try(suppressWarnings(sf::st_transform(gbif_sf, crs_target)), silent = TRUE)
      ply <- try(suppressWarnings(sf::st_transform(eez_keep, crs_target)), silent = TRUE)
      if (inherits(gb, "try-error") || inherits(ply, "try-error")) stop("transform_failed")
      out <- try(suppressWarnings(
        sf::st_join(gb, ply[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
                    join = sf::st_within, left = TRUE)
      ), silent = TRUE)
      if (inherits(out, "try-error") || (!("geoname" %in% names(out)) || all(is.na(out$geoname)))) {
        out <- try(suppressWarnings(
          sf::st_join(gb, ply[, c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best")],
                      join = sf::st_nearest_feature, left = TRUE)
        ), silent = TRUE)
      }
      if (inherits(out, "try-error")) stop(attr(out, "condition")$message)
      out
    }

    res <- try(try_join(force_planar = FALSE), silent = TRUE)
    if (inherits(res, "try-error")) {
      # any s2/wk failure -> hard fallback to planar
      res <- try(try_join(force_planar = TRUE), silent = TRUE)
      if (inherits(res, "try-error")) {
        # as a last resort, return input with empty columns
        gbif_sf$geoname <- NA_character_
        gbif_sf$mrgid <- NA_integer_
        gbif_sf$iso3_best <- NA_character_
        gbif_sf$iso2_best <- NA_character_
        res <- gbif_sf
      }
    }
    joined_list[[key]] <- res
  }

  joined_list
}

#' Join GBIF points to Marine Regions overlays (robust, with s2 fallback)
#'
#' Reads/caches a Marine Regions polygon layer (EEZ/LME/IHO/MEOW/FAO, etc.),
#' fixes invalid geometries, harmonizes name/id columns, and joins each sf in
#' `results_list` to the polygons. If any s2 error occurs, the join
#' transparently falls back to planar (sf_use_s2(FALSE)).
#'
#' @param results_list Named list of sf point data (e.g., species -> sf).
#' @param overlay Character; which overlay family to use. One of
#'   c("eez","lme","iho","meow","fao","custom"). Only affects defaults if
#'   `overlay_layer` is not supplied.
#' @param overlay_layer Character; exact layer id for Marine Regions. If given,
#'   it is used directly. If NULL, a canonical id is chosen for the `overlay`.
#'   Examples: "eez", "lme", "iho", "meow", "fao".
#' @param overlay_sf Optional `sf` polygon object (preloaded overlay). If
#'   provided, no download/caching is attempted.
#' @param use_planar Logical; if TRUE, run joins with s2 disabled from the
#'   start (planar ops). Regardless, the function auto-falls back to planar
#'   if a spherical join fails.
#' @param cache_dir Directory to cache the overlay as a GeoPackage. Must be
#'   supplied explicitly when `overlay_sf = NULL`. In examples, tests and
#'   vignettes, use a path under `tempdir()`.
#'
#' @return A named list with the same names as `results_list`. Each element is
#'   an `sf` point object containing the original occurrence records plus
#'   standardised attributes from the matched Marine Regions overlay. The
#'   returned columns may include `geoname`, `mrgid`, `territory1`,
#'   `sovereign1`, `iso3_best` and `iso2_best`, depending on the overlay used.
#'   `geoname` is the matched marine-region label, `mrgid` is the Marine
#'   Regions identifier when available, and `iso3_best`/`iso2_best` are
#'   best-effort ISO country codes mainly available for Exclusive Economic Zone
#'   overlays. Records that cannot be assigned retain their input geometry and
#'   receive missing values in the added fields.
#'
#' @details
#' - Uses **mregions2** through `mrp_list`/`mrp_get`.
#'   Legacy Marine Regions fallback support has been removed.
#' - For non-EEZ overlays, `iso3_best`/`iso2_best` will likely be `NA`; your
#'   pipeline's `ensure_iso3_on_points()` will still fill ISO by overlaying
#'   a world-country layer later.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   pts <- data.frame(
#'     decimalLongitude = c(-5, 10),
#'     decimalLatitude = c(50, 35)
#'   )
#'
#'   pts_sf <- sf::st_as_sf(
#'     pts,
#'     coords = c("decimalLongitude", "decimalLatitude"),
#'     crs = 4326,
#'     remove = FALSE
#'   )
#'
#'   joined <- overlay_join(
#'     list("Carcinus maenas" = pts_sf),
#'     overlay = "eez",
#'     cache_dir = file.path(tempdir(), "biofetchR_overlay_cache")
#'   )
#'
#'   names(joined[["Carcinus maenas"]])
#' }
#' }
#'
#' @section Data access and attribution:
#' Marine Regions overlays are accessed through \pkg{mregions2}
#' and cached locally. biofetchR does not redistribute these layers. Users should
#' cite Marine Regions and check the data provider's licence and attribution
#' requirements for the specific overlay used.
#'
#' @family marine overlay joins
#' @export
overlay_join <- function(results_list,
                         overlay = c("eez","lme","iho","meow","fao","custom"),
                         overlay_layer = NULL,
                         overlay_sf = NULL,
                         use_planar = TRUE,
                         cache_dir = NULL) {
  overlay <- match.arg(overlay)

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # helpers to discover & fetch overlays
  list_overlays <- function(pattern = NULL) {
    has_mr2 <- requireNamespace("mregions2", quietly = TRUE)

    if (has_mr2 && "mrp_list" %in% getNamespaceExports("mregions2")) {
      L <- mregions2::mrp_list()
      L <- as.data.frame(L)
      if (!"layername" %in% names(L) && "layer" %in% names(L)) L$layername <- L$layer
    } else if (FALSE) {
      stop("Legacy Marine Regions fallback has been removed; use mregions2.", call. = FALSE)
      if (!"layername" %in% names(L) && "name" %in% names(L)) L$layername <- L$name
    } else {
      stop("Need mregions2 with mrp_list/mrp_get.")
    }

    if (!is.null(pattern)) {
      keep <- grepl(pattern, L$title, ignore.case = TRUE) |
        ("layername" %in% names(L) && grepl(pattern, L$layername, ignore.case = TRUE))
      L <- L[keep, , drop = FALSE]
    }
    L
  }

  get_overlay_sf <- function(id_or_family, which = 1, cache_dir = cache_dir, quiet = FALSE) {
    # If a direct id works, use it; else search catalog by family/id.
    L <- try(list_overlays(), silent = TRUE)
    layer <- NA_character_

    if (!inherits(L, "try-error")) {
      # choose canonical id for family if none given
      canon <- switch(tolower(id_or_family),
                      eez  = "eez",
                      lme  = "lme",
                      iho  = "iho",
                      meow = "meow",
                      fao  = "fao",
                      id_or_family)
      if (canon %in% L$layername) {
        layer <- canon
      } else {
        hits <- list_overlays(id_or_family)
        if (!nrow(hits)) stop("No layers match: ", id_or_family)
        if (which < 1 || which > nrow(hits)) stop("`which` out of range (", nrow(hits), " hits).")
        layer <- hits$layername[which]
      }
    } else {
      layer <- id_or_family
    }

    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    safe <- gsub("[^A-Za-z0-9_-]+", "_", layer)
    cache_path <- file.path(cache_dir, paste0("overlay_", safe, ".gpkg"))

    if (file.exists(cache_path)) {
      if (!quiet) message("Using cached overlay: ", cache_path)
      ov <- try(sf::read_sf(cache_path), silent = TRUE)
      if (!inherits(ov, "try-error")) return(ov)
    }

    if (requireNamespace("mregions2", quietly = TRUE) &&
        "mrp_get" %in% getNamespaceExports("mregions2")) {
      ov <- mregions2::mrp_get(layer)
        } else if (FALSE) {
      stop("Legacy Marine Regions fallback has been removed; use mregions2.", call. = FALSE)
    } else {
      stop("Need mregions2::mrp_get to fetch overlay.")
    }

    ov <- suppressWarnings(sf::st_make_valid(ov))
    # cache (best effort)
    try(sf::write_sf(ov, cache_path), silent = TRUE)
    ov
  }

  pick_ci <- function(x, candidates) {
    nx <- names(x)
    for (cand in candidates) {
      hit <- which(tolower(nx) == tolower(cand))
      if (length(hit)) return(nx[hit[1]])
    }
    NULL
  }

  # fetch overlay (preloaded, or download+cache)
  ov <- overlay_sf
  if (is.null(ov)) {
    cache_dir <- .bf_require_explicit_spatial_cache(
      cache_dir,
      context = "overlay_join()"
    )

    layer_id <- overlay_layer %||% overlay
    ov <- get_overlay_sf(layer_id, which = 1, cache_dir = cache_dir, quiet = FALSE)
  }

  # try to standardize a label and id column across overlays
  # label candidates by family (first hit wins)
  name_col <- pick_ci(ov, c(
    # EEZ
    "geoname","GEONAME","GeoName","name",
    # LME
    "LME_NAME","lme_name","LME_NAME_LME66",
    # IHO
    "NAME","name","IHO_SEA",
    # MEOW
    "ECOREGION","ecoregion","ECO_NAME","eco_name",
    # FAO
    "F_AREA","F_LEVEL","OCEAN","title"
  ))

  # id (MRGID is only guaranteed on some layers, e.g., EEZ)
  mrgid_col <- pick_ci(ov, c("mrgid","MRGID","Mrgid","MRGID_EEZ","MRGID_POLY"))

  # Try to carry useful country-ish fields for EEZ ISO derivation
  terr1_col     <- pick_ci(ov, c("TERRITORY1","territory1","territory"))
  sov1_col      <- pick_ci(ov, c("SOVEREIGN1","sovereign1","sovereign"))
  iso_ter1_col  <- pick_ci(ov, c("ISO_Ter1","iso_ter1","ISO_TER1"))
  iso_sov1_col  <- pick_ci(ov, c("ISO_SOV1","iso_sov1","ISO_SOV"))
  iso2_ter1_col <- pick_ci(ov, c("ISO2_Ter1","iso2_ter1","ISO2_TER1"))
  iso2_sov1_col <- pick_ci(ov, c("ISO2_SOV1","iso2_sov1","ISO2_SOV"))

  keep_cols <- unique(na.omit(c(name_col, mrgid_col,
                                terr1_col,  sov1_col,
                                iso_ter1_col, iso_sov1_col,
                                iso2_ter1_col, iso2_sov1_col)))
  ov_keep <- ov[, keep_cols, drop = FALSE]

  # standardize names used downstream
  if (!is.null(name_col)  && name_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==name_col]   <- "geoname"
  if (!is.null(mrgid_col) && mrgid_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==mrgid_col]  <- "mrgid"
  if (!is.null(terr1_col) && terr1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==terr1_col]  <- "territory1"
  if (!is.null(sov1_col)  && sov1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==sov1_col]   <- "sovereign1"
  if (!is.null(iso_ter1_col)  && iso_ter1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso_ter1_col]  <- "iso_ter1"
  if (!is.null(iso_sov1_col)  && iso_sov1_col  %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso_sov1_col]  <- "iso_sov1"
  if (!is.null(iso2_ter1_col) && iso2_ter1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso2_ter1_col] <- "iso2_ter1"
  if (!is.null(iso2_sov1_col) && iso2_sov1_col %in% names(ov_keep)) names(ov_keep)[names(ov_keep)==iso2_sov1_col] <- "iso2_sov1"

  if (!"geoname" %in% names(ov_keep)) {
    # last resort: take any non-geometry column
    non_geom <- setdiff(names(ov_keep), attr(ov_keep, "sf_column"))
    if (length(non_geom)) ov_keep$geoname <- as.character(ov_keep[[non_geom[1]]]) else ov_keep$geoname <- NA_character_
  }

  # compute iso3_best / iso2_best (robust for EEZ; NA for other overlays)
  gcol <- function(df, nm) if (nm %in% names(df)) df[[nm]] else NULL
  coalesce_char <- function(...) {
    args <- list(...); n <- nrow(ov_keep); out <- rep(NA_character_, n)
    for (v in args) {
      if (is.null(v)) next
      v <- as.character(v); fill <- is.na(out) | !nzchar(out); out[fill] <- v[fill]
    }
    out
  }

  iso3_best <- iso2_best <- rep(NA_character_, nrow(ov_keep))
  if (tolower(overlay) == "eez" || (!is.null(overlay_layer) && grepl("^eez$", overlay_layer, ignore.case = TRUE))) {
    iso3_direct <- coalesce_char(gcol(ov_keep, "iso_ter1"), gcol(ov_keep, "iso_sov1"))
    if (requireNamespace("countrycode", quietly = TRUE)) {
      need <- is.na(iso3_direct) | !nzchar(iso3_direct)
      if (any(need)) {
        nm_guess <- coalesce_char(gcol(ov_keep, "territory1"), gcol(ov_keep, "sovereign1"), gcol(ov_keep, "geoname"))
        nm_guess <- gsub("Exclusive Economic Zone|Joint regime area|Overlapping claim area|\\(.*?\\)|-{1,2}\\s*EEZ|\\s+EEZ",
                         "", nm_guess, ignore.case = TRUE)
        nm_guess <- trimws(nm_guess)
        iso3_from_name <- suppressWarnings(countrycode::countrycode(nm_guess, "country.name", "iso3c"))
        iso3_direct[need] <- toupper(iso3_from_name[need])
      }
    }
    iso3_best <- toupper(iso3_direct)

    iso2_direct <- coalesce_char(gcol(ov_keep, "iso2_ter1"), gcol(ov_keep, "iso2_sov1"))
    if ((is.null(iso2_direct) || all(is.na(iso2_direct))) && requireNamespace("countrycode", quietly = TRUE)) {
      iso2_direct <- suppressWarnings(countrycode::countrycode(iso3_best, "iso3c", "iso2c"))
    }
    iso2_best <- toupper(iso2_direct)
  }

  ov_keep$iso3_best <- iso3_best
  ov_keep$iso2_best <- iso2_best

  # s2 handling and robust fallback
  old_s2 <- sf::sf_use_s2()
  on.exit(try(sf::sf_use_s2(old_s2), silent = TRUE), add = TRUE)
  if (isTRUE(use_planar)) sf::sf_use_s2(FALSE)

  out_list <- list()
  cols_to_join <- intersect(c("geoname","mrgid","territory1","sovereign1","iso3_best","iso2_best"), names(ov_keep))

  for (key in names(results_list)) {
    pts <- results_list[[key]]
    if (is.null(pts) || !inherits(pts, "sf") || !nrow(pts)) {
      if (inherits(pts, "sf")) {
        if (!"geoname"   %in% names(pts)) pts$geoname   <- NA_character_
        if (!"mrgid"     %in% names(pts)) pts$mrgid     <- NA_integer_
        if (!"iso3_best" %in% names(pts)) pts$iso3_best <- NA_character_
        if (!"iso2_best" %in% names(pts)) pts$iso2_best <- NA_character_
      }
      out_list[[key]] <- pts
      next
    }

    try_join <- function(force_planar = FALSE) {
      if (force_planar) sf::sf_use_s2(FALSE)
      crs_target <- if (!force_planar && !isTRUE(use_planar)) sf::st_crs(ov_keep) else 3857
      pts_tr <- try(suppressWarnings(sf::st_transform(pts, crs_target)), silent = TRUE)
      ply_tr <- try(suppressWarnings(sf::st_transform(ov_keep, crs_target)), silent = TRUE)
      if (inherits(pts_tr, "try-error") || inherits(ply_tr, "try-error")) stop("transform_failed")

      out <- try(suppressWarnings(
        sf::st_join(pts_tr, ply_tr[, cols_to_join, drop = FALSE], join = sf::st_within, left = TRUE)
      ), silent = TRUE)

      if (inherits(out, "try-error") || (!("geoname" %in% names(out)) || all(is.na(out$geoname)))) {
        # nearest as a last-ditch effort (kept because it salvages geometries that barely miss boundaries)
        out <- try(suppressWarnings(
          sf::st_join(pts_tr, ply_tr[, cols_to_join, drop = FALSE], join = sf::st_nearest_feature, left = TRUE)
        ), silent = TRUE)
      }
      if (inherits(out, "try-error")) stop(attr(out, "condition")$message)
      out
    }

    res <- try(try_join(FALSE), silent = TRUE)
    if (inherits(res, "try-error")) {
      res <- try(try_join(TRUE),  silent = TRUE)
      if (inherits(res, "try-error")) {
        # absolute fallback: return input plus empty columns
        if (!"geoname"   %in% names(pts)) pts$geoname   <- NA_character_
        if (!"mrgid"     %in% names(pts)) pts$mrgid     <- NA_integer_
        if (!"iso3_best" %in% names(pts)) pts$iso3_best <- NA_character_
        if (!"iso2_best" %in% names(pts)) pts$iso2_best <- NA_character_
        res <- pts
      }
    }
    out_list[[key]] <- res
  }

  out_list
}

#' Load & cache WWF Terrestrial Ecoregions (TEOW)
#'
#' Downloads (once) and caches the WWF terrestrial ecoregions as an sf object.
#' Prefer `method = "mapme"` if the \pkg{mapme.biodiversity} package is installed;
#' otherwise, use `method = "direct"` with a `source_url` to the official TEOW ZIP.
#'
#' The cached, normalized layer is written to `<cache_dir>/teow/teow.gpkg`
#' with the ecoregion name exposed as `geoname` and the geometry column
#' standardized to `"geometry"`.
#'
#' @param cache_dir Cache root used for the TEOW cache. Must be supplied
#'   explicitly. In examples, tests and vignettes, use a path under `tempdir()`.
#' @param method Character. One of \code{"auto"}, \code{"mapme"} or \code{"direct"}.
#' @param source_url Optional direct URL to a TEOW ZIP file when \code{method = "direct"}.
#' @param aoi Optional \code{sf} or \code{sfc} object used to clip the TEOW layer.
#' @param bbox Optional bounding box used to clip the TEOW layer.
#' @param force_refresh Logical; if TRUE, rebuild the TEOW cache even if it already exists.
#' @param simplify Logical; if TRUE, simplify TEOW geometries after loading.
#' @param simplify_tolerance Numeric simplification tolerance used when \code{simplify = TRUE}.
#' @param quiet Logical; if TRUE, suppress progress messages.
#'
#' @return An `sf` polygon object containing Terrestrial Ecoregions of the World
#'   features in WGS84 longitude/latitude coordinates. The returned object
#'   includes a standard `geoname` column and, where available, original
#'   ecoregion, biome and realm attributes. If `aoi` or `bbox` is supplied, the
#'   returned layer is clipped to that area; otherwise it contains the cached or
#'   freshly loaded global layer. The normalised layer is also cached locally as
#'   a GeoPackage for reuse.
#'
#' @section Data access and attribution:
#' TEOW is a third-party terrestrial ecoregion product. biofetchR downloads and
#' caches it locally for user-side analysis only. Users should cite the original
#' provider and check licence and redistribution terms before sharing derived
#' products.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   teow <- bf_load_teow(
#'     cache_dir = file.path(tempdir(), "biofetchR_teow"),
#'     quiet = TRUE
#'   )
#'
#'   names(teow)
#' }
#' }
#'
#' @family terrestrial overlay loaders
#' @export
bf_load_teow <- function(
    cache_dir = NULL,
    method = c("auto", "mapme", "direct"),
    source_url = NULL,
    aoi = NULL,
    bbox = NULL,
    force_refresh = FALSE,
    simplify = FALSE,
    simplify_tolerance = 0.01,
    quiet = FALSE
) {
  method <- match.arg(method)

  cache_dir <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_load_teow()"
  )

  teow_dir  <- file.path(cache_dir, "teow")
  teow_gpkg <- file.path(teow_dir, "teow.gpkg")
  if (!dir.exists(teow_dir)) dir.create(teow_dir, recursive = TRUE, showWarnings = FALSE)

  # helpers
  .global_bbox_sf <- function() {
    coords <- matrix(
      c(-179.999, -89.999,
        179.999, -89.999,
        179.999,  89.999,
        -179.999,  89.999,
        -179.999, -89.999),
      ncol = 2, byrow = TRUE
    )
    sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(coords)), crs = 4326))
  }
  .as_sf <- function(g) {
    if (inherits(g, "sf"))  return(g)
    if (inherits(g, "sfc")) return(sf::st_sf(geometry = g))
    stop("AOI must be 'sf' or 'sfc'.")
  }

  .read_best_teow <- function(path, quiet = FALSE) {
    if (grepl("\\.gpkg$", path, ignore.case = TRUE)) {
      lay <- tryCatch(sf::st_layers(path), error = function(e) NULL)
      if (!is.null(lay)) {
        best <- NULL; best_n <- -Inf
        for (ln in lay$name) {
          x <- tryCatch(sf::read_sf(path, layer = ln, quiet = quiet), error = function(e) NULL)
          if (is.null(x)) next
          has_cols <- any(c("ECO_NAME","eco_name","ECO_ID","ECOREGION","ecoregion","ECO_CODE") %in% names(x))
          if (has_cols && nrow(x) > best_n) { best <- x; best_n <- nrow(x) }
        }
        if (!is.null(best)) return(best)
      }
    }
    sf::read_sf(path, quiet = quiet)
  }

  aoi_clip  <- aoi
  bbox_clip <- bbox

  # fast path: cached
  if (!force_refresh && file.exists(teow_gpkg)) {
    out <- .read_best_teow(teow_gpkg, quiet = quiet)
    out <- bf_sf_wgs84(out)
    out <- bf_sf_standardise_geometry(out)

    if (nrow(out) < 100 && !force_refresh) {
      if (!quiet) {
        message("Cached TEOW has only ", nrow(out), " features; refetching via 'direct'.")
      }
      force_refresh <- TRUE
    } else {
      return(.bf_teow_clip(out, aoi_clip, bbox_clip, quiet))
    }
  }

  # choose method
  if (method == "auto") {
    method <- if (requireNamespace("mapme.biodiversity", quietly = TRUE)) "mapme" else "direct"
  }

  # fetch
  teow <- NULL

  if (identical(method, "mapme")) {
    if (!requireNamespace("mapme.biodiversity", quietly = TRUE)) {
      stop("mapme.biodiversity not installed; use method='direct' with a source_url.")
    }
    aoi_dl <- if (is.null(aoi) && is.null(bbox)) {
      .global_bbox_sf()
    } else if (!is.null(bbox)) {
      bf_bbox_to_sf(bbox)
    } else {
      bf_sf_wgs84(.as_sf(aoi))
    }

    mapme.biodiversity::mapme_options(outdir = teow_dir, verbose = FALSE)
    invisible(mapme.biodiversity::get_resources(aoi_dl, mapme.biodiversity::get_teow()))

    cand <- c(
      list.files(teow_dir, pattern = "wwf_terr_ecos\\.gpkg$", full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "\\.gpkg$",              full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "wwf_terr_ecos\\.shp$",  full.names = TRUE, recursive = TRUE),
      list.files(teow_dir, pattern = "\\.shp$",               full.names = TRUE, recursive = TRUE)
    )
    cand <- unique(cand)
    if (length(cand)) {
      info <- file.info(cand)
      path <- cand[order(info$mtime, decreasing = TRUE)][1]
      teow <- .read_best_teow(path, quiet = quiet)
    } else {
      if (is.null(source_url)) {
        source_url <- "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
      }
      if (!quiet) message("Mapme did not materialize TEOW; falling back to direct.")
      method <- "direct"
    }
  }

  if (identical(method, "direct")) {
    if (is.null(source_url)) stop("For method='direct', provide source_url to the official TEOW ZIP.")
    zip_path  <- file.path(teow_dir, "teow.zip")
    unzip_dir <- file.path(teow_dir, "unzipped")

    bf_download_cached(
      url = source_url,
      dest = zip_path,
      force_refresh = force_refresh,
      quiet = quiet,
      min_bytes = 1024,
      validate_not_html = TRUE
    )

    bf_unpack_archive(
      path = zip_path,
      exdir = unzip_dir,
      force_refresh = TRUE,
      quiet = quiet
    )

    shp <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (!length(shp)) stop("No .shp found inside the TEOW ZIP.")
    teow <- sf::read_sf(shp[1], quiet = quiet)
  }

  if (!inherits(teow, "sf")) stop("TEOW layer could not be read as an sf object.")

  # Ensure WGS84 and standardise geometry column name
  teow <- bf_sf_wgs84(teow)
  teow <- bf_sf_standardise_geometry(teow)

  # Sanity check: if suspiciously small and not already 'direct', refetch via 'direct'
  if (nrow(teow) < 100 && !identical(method, "direct")) {
    if (!quiet) message("TEOW layer has only ", nrow(teow), " features; refetching via 'direct'.")
    zip_path  <- file.path(teow_dir, "teow.zip")
    unzip_dir <- file.path(teow_dir, "unzipped")
    if (is.null(source_url)) {
      source_url <- "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
    }

    bf_download_cached(
      url = source_url,
      dest = zip_path,
      force_refresh = TRUE,
      quiet = quiet,
      min_bytes = 1024,
      validate_not_html = TRUE
    )

    bf_unpack_archive(
      path = zip_path,
      exdir = unzip_dir,
      force_refresh = TRUE,
      quiet = quiet
    )

    shp <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (!length(shp)) stop("No .shp found inside the TEOW ZIP (refetch).")
    teow <- sf::read_sf(shp[1], quiet = quiet)
    teow <- bf_sf_wgs84(teow)
    teow <- bf_sf_standardise_geometry(teow)
  }

  # Normalize fields
  name_candidates <- c("ECO_NAME","eco_name","ECOREGION","ecoregion","ECO_NAME_E","ECO_ID","ECO_CODE")
  nm <- intersect(name_candidates, names(teow))
  teow$geoname <- if (length(nm)) as.character(teow[[nm[1]]]) else paste0("ECO_", seq_len(nrow(teow)))

  # Keep useful attrs + geometry
  keep <- intersect(c("geoname","ECO_NAME","ECO_ID","ECO_CODE","BIOME_NAME","BIOME","REALM"), names(teow))
  teow <- teow[, unique(c(keep, "geometry")), drop = FALSE]

  # Validity + optional simplify
  teow <- bf_sf_make_valid(teow, quiet = quiet)
  if (isTRUE(simplify)) {
    teow <- sf::st_simplify(teow, dTolerance = simplify_tolerance, preserveTopology = TRUE)
  }

  # Optional clip (use original user AOI/bbox)
  teow <- .bf_teow_clip(teow, aoi_clip, bbox_clip, quiet)

  # Write normalized cache
  if (file.exists(teow_gpkg)) unlink(teow_gpkg)
  sf::write_sf(teow, teow_gpkg, quiet = quiet)
  if (!quiet) message("TEOW cached at: ", teow_gpkg)

  teow
}

#' Clip TEOW polygons to an area of interest or bounding box
#'
#' Internal helper used by [bf_load_teow()] to clip terrestrial ecoregions after
#' loading or rebuilding the local cache.
#'
#' @param x An `sf` polygon object.
#' @param aoi Optional `sf` or `sfc` clipping geometry.
#' @param bbox Optional numeric bounding box, supplied as `c(xmin, ymin, xmax, ymax)`.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` object clipped to `aoi` or `bbox`.
#'
#' @keywords internal
#' @noRd
.bf_teow_clip <- function(x, aoi = NULL, bbox = NULL, quiet = TRUE) {
  stopifnot(inherits(x, "sf"))

  # Nothing to do when no clipping geometry was supplied.
  if (is.null(aoi) && is.null(bbox)) {
    return(x)
  }

  # Local coercion is retained because this helper accepts either sf or sfc AOIs.
  .as_sf <- function(g) {
    if (inherits(g, "sf")) {
      return(g)
    }

    if (inherits(g, "sfc")) {
      return(sf::st_sf(geometry = g))
    }

    stop("AOI must be sf/sfc.", call. = FALSE)
    }

  # Build clip geometry in WGS84.
  clipper <- if (!is.null(bbox)) {
    bf_bbox_to_sf(bbox)
  } else {
    bf_sf_wgs84(.as_sf(aoi))
  }

  # Repair geometries using the shared spatial helper.
  x <- bf_sf_make_valid(x, quiet = quiet)
  clipper <- bf_sf_make_valid(clipper, quiet = quiet)

  # Fast prefilter by intersection. This reduces geometry complexity before the
  # more expensive precise intersection step.
  x_pref <- tryCatch(
    sf::st_filter(x, clipper, .predicate = sf::st_intersects, drop = FALSE),
    error = function(e) x
  )

  # Precise clip; if it fails, return the prefiltered object.
  out <- tryCatch(
    suppressWarnings(sf::st_intersection(x_pref, clipper)),
    error = function(e) x_pref
  )

  if (!quiet) {
    message("TEOW clipped to AOI/BBOX (", nrow(out), " features).")
  }

  out
}

#' Return cache path/info for TEOW
#'
#' @param cache_dir Cache root containing the TEOW cache. Must be supplied
#'   explicitly. In examples, tests and vignettes, use a path under `tempdir()`.
#'
#' @return A named list describing the Terrestrial Ecoregions of the World cache.
#'   The list contains `path`, a character string giving the expected GeoPackage
#'   cache path; `exists`, a logical value indicating whether the cache file is
#'   present; `size_bytes`, a numeric file size in bytes or `NA_real_` when the
#'   file is absent; and `modified`, a `POSIXct` modification time or `NA` when
#'   the file is absent. The values are intended for cache diagnostics and do
#'   not load the spatial layer itself.
#'
#' @family terrestrial overlay loaders
#' @export
bf_teow_cache_info <- function(cache_dir = NULL) {
  cache_dir <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_teow_cache_info()",
    create = FALSE
  )
  teow_gpkg <- file.path(cache_dir, "teow", "teow.gpkg")
  info <- if (file.exists(teow_gpkg)) file.info(teow_gpkg) else NULL
  list(
    path = teow_gpkg,
    exists = file.exists(teow_gpkg),
    size_bytes = if (!is.null(info)) unname(info$size) else NA_real_,
    modified = if (!is.null(info)) info$mtime else NA
  )
}

#' Clear the cached TEOW dataset
#'
#' @param cache_dir Cache root containing the TEOW cache to remove. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param ask If TRUE, prompt before deleting. Set FALSE for non-interactive use.
#'
#' @return Invisibly returns `TRUE` when the Terrestrial Ecoregions of the World
#'   cache directory does not exist or is successfully removed. Invisibly
#'   returns `FALSE` when deletion is cancelled after an interactive prompt.
#'   The function is called primarily for its side effect of removing cached
#'   TEOW files.
#'
#' @family terrestrial overlay loaders
#' @export
bf_teow_clear_cache <- function(cache_dir = NULL, ask = interactive()) {
  cache_dir <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_teow_clear_cache()",
    create = FALSE
  )
  teow_dir <- file.path(cache_dir, "teow")
  if (!dir.exists(teow_dir)) return(invisible(TRUE))
  if (ask) {
    ans <- utils::menu(c("No", "Yes"), title = sprintf("Delete TEOW cache at %s ?", teow_dir))
    if (ans != 2) return(invisible(FALSE))
  }
  unlink(teow_dir, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}


#' Load Freshwater Ecoregions of the World polygons
#'
#' Loads the Freshwater Ecoregions of the World (FEOW) polygon layer, standardises
#' key fields, and caches the result as a single GeoPackage. The default
#' \code{method = "auto"} first attempts a direct ArcGIS FeatureServer query,
#' then falls back to a downloadable archive if
#' available.
#'
#' @param force_refresh Logical. If \code{TRUE}, ignore any existing cached FEOW
#'   GeoPackage and rebuild the cache.
#' @param quiet Logical. If \code{TRUE}, suppress progress messages.
#' @param cache_dir Cache root used for the FEOW cache. Must be supplied
#'   explicitly. In examples, tests and vignettes, use a path under `tempdir()`.
#' @param method Character. One of \code{"auto"}, \code{"arcgis"},
#'   \code{"download"} or \code{"download"}. \code{"auto"} tries the methods in
#'   that order.
#' @param source_url Optional source URL used when \code{method = "download"}.
#'   If omitted, the FEOW downloads page is queried for a ZIP URL.
#' @param ... Ignored. Retained for backward compatibility with older calls.
#'
#' @return An `sf` polygon object containing Freshwater Ecoregions of the World
#'   features in WGS84 longitude/latitude coordinates. The returned object
#'   includes standardised columns such as `feow_id`, `ecoregion`, `biome` and
#'   `realm`, where these fields are available from the source layer. Each row
#'   represents a freshwater ecoregion polygon used for assigning occurrence
#'   records to freshwater spatial units.
#'
#' @section Data access and attribution:
#' FEOW is a third-party freshwater ecoregion product. biofetchR only downloads
#' and caches the layer locally for the user. Users should cite FEOW and check the
#' provider's licence and redistribution requirements before using outputs in
#' publications or shared data products.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   feow <- bf_load_feow(
#'     cache_dir = file.path(tempdir(), "biofetchR_feow"),
#'     quiet = TRUE
#'   )
#'
#'   names(feow)
#' }
#' }
#'
#' @family freshwater overlay loaders
#' @export
bf_load_feow <- function(force_refresh = FALSE,
                         quiet = TRUE,
                         cache_dir = NULL,
                         method = "auto",
                         source_url = NULL,
                         ...) {
  FEOW_ARCGIS_LAYER <- "https://services.arcgis.com/uUvqNMGPm7axC2dD/arcgis/rest/services/FEOWv1_TNC/FeatureServer/0"

  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_feow(). Please install it.")
  }

  method <- tolower(method)
  if (!method %in% c("auto","arcgis","download","feowr")) {
    stop("bf_load_feow(): unknown method: ", method)
  }

  cache_root <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_load_feow()"
  )
  feow_dir   <- file.path(cache_root, "feow")
  if (!dir.exists(feow_dir)) dir.create(feow_dir, recursive = TRUE, showWarnings = FALSE)
  feow_gpkg  <- file.path(feow_dir, "feow.gpkg")

  # Use cache if present
  if (file.exists(feow_gpkg) && !isTRUE(force_refresh)) {
    if (!quiet) message("biofetchR: FEOW cache found: ", feow_gpkg)
    return(.bf_read_and_standardise_feow(feow_gpkg, quiet = quiet))
  }

  # Fetchers
  fetch_arcgis <- function() {
    if (!quiet) message("biofetchR: Fetching FEOW from ArcGIS FeatureServer...")
    x <- .bf_arcgis_fetch_all(FEOW_ARCGIS_LAYER, quiet = quiet)
    if (!is.null(x) && is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    x
  }

  fetch_download <- function() {
    url <- source_url
    if (is.null(url)) {
      if (!quiet) message("biofetchR: Discovering FEOW download URL...")
      url <- .bf_feow_find_zip_url("https://www.feow.org/download")
    }
    if (is.null(url)) return(NULL)
    if (!quiet) message("biofetchR: Downloading FEOW archive...")

    tf <- tempfile(fileext = if (grepl("\\.zip($|\\?)", url, ignore.case = TRUE)) ".zip" else "")
    bf_download_file(url, tf, quiet = quiet)

    root <- if (grepl("\\.zip($|\\?)", tolower(tf))) {
      td <- tempfile("feow_unzip_")
      bf_unpack_archive(
        path = tf,
        exdir = td,
        force_refresh = TRUE,
        quiet = quiet
      )
      td
    } else {
      dirname(tf)
    }

    .bf_read_any_vector(root, quiet = quiet)
  }

  fetch_feowr <- function() {
    return(NULL)
    if (!is.null(x) && is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    x
  }

  # Orchestrate per method
  feow <- switch(method,
                 arcgis   = fetch_arcgis(),
                 download = fetch_download(),
                 feowr    = fetch_feowr(),
                 auto     = { x <- fetch_arcgis(); if (is.null(x)) x <- fetch_download(); if (is.null(x)) x <- fetch_feowr(); x }
  )

  if (is.null(feow)) {
    stop("bf_load_feow(): failed to obtain FEOW via method='", method,
         "'. Tried ArcGIS, downloadable archive, (if installed).")
  }

  suppressWarnings(sf::st_write(feow, feow_gpkg, delete_dsn = TRUE, quiet = TRUE))
  .bf_read_and_standardise_feow(feow_gpkg, quiet = quiet)
}

# -----------------------------------------------------------------------------
# Helpers (internal)
# -----------------------------------------------------------------------------

#' Read and standardise a cached FEOW GeoPackage
#'
#' @param gpkg_path Path to a cached FEOW GeoPackage.
#' @param quiet Logical. If `TRUE`, suppress read messages.
#'
#' @return An `sf` object with standardised FEOW columns.
#'
#' @keywords internal
#' @noRd
.bf_read_and_standardise_feow <- function(gpkg_path, quiet = TRUE) {
  x  <- sf::st_read(gpkg_path, quiet = quiet)
  nm <- names(x)
  x$feow_id   <- if ("ECO_ID_U" %in% nm) x$ECO_ID_U else if ("ECO_ID" %in% nm) x$ECO_ID else if ("FEOW_ID" %in% nm) x$FEOW_ID else if ("ECO_CODE" %in% nm) x$ECO_CODE else NA
  x$ecoregion <- if ("ECOREGION" %in% nm) x$ECOREGION else if ("ECO_NAME" %in% nm) x$ECO_NAME else NA
  x$biome     <- if ("MHT_TXT"  %in% nm) x$MHT_TXT  else if ("BIOME" %in% nm) x$BIOME else NA
  x$realm     <- if ("REALM"    %in% nm) x$REALM    else NA
  x <- sf::st_make_valid(x)
  sf::st_set_agr(x, "constant")
  x
}

#' Find a FEOW ZIP download URL from the FEOW downloads page
#'
#' @param feow_downloads_page URL of the FEOW downloads page.
#'
#' @return Character URL, or `NULL` when no suitable ZIP link can be found.
#'
#' @keywords internal
#' @noRd
.bf_feow_find_zip_url <- function(feow_downloads_page) {
  raw <- tryCatch(
    suppressWarnings(readLines(feow_downloads_page, warn = FALSE)),
    error = function(e) NULL
  )
  if (is.null(raw)) return(NULL)
  m     <- regmatches(raw, gregexpr('href\\s*=\\s*"([^"]+)"', raw, perl = TRUE))
  hrefs <- unique(unlist(m))
  if (!length(hrefs)) return(NULL)
  urls <- sub('^href\\s*=\\s*"', "", hrefs)
  urls <- sub('"$', "", urls)
  urls <- urls[grepl("\\.zip(\\?.*)?$", urls, ignore.case = TRUE)]
  if (!length(urls)) return(NULL)
  urls <- ifelse(
    grepl("^https?://", urls, ignore.case = TRUE),
    urls,
    paste0("https://www.feow.org/", sub("^/", "", urls))
  )
  urls[[1]]
}

#' Read the first supported vector layer under a directory
#'
#' @param dirpath Directory to search recursively for GeoPackage, shapefile or
#'   GeoJSON vector data.
#' @param quiet Logical. If `TRUE`, suppress read messages.
#'
#' @return An `sf` object.
#'
#' @keywords internal
#' @noRd
.bf_read_any_vector <- function(dirpath, quiet = TRUE) {
  gpkg <- list.files(dirpath, pattern = "\\.gpkg$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(gpkg)) return(sf::st_read(gpkg[1], quiet = quiet))

  gdbs <- list.dirs(dirpath, full.names = TRUE, recursive = TRUE)
  gdbs <- gdbs[grepl("\\.gdb$", gdbs, ignore.case = TRUE)]
  if (length(gdbs)) {
    lay <- try(sf::st_layers(gdbs[1]), silent = TRUE)
    if (!inherits(lay, "try-error")) {
      nms <- lay$name
      gtypes <- NULL
      if (!is.null(lay$geomtype)) gtypes <- lay$geomtype else if (!is.null(lay$geometry_type)) gtypes <- lay$geometry_type
      score <- tolower(paste(nms, gtypes))
      idx <- grep("poly", score)
      lyr <- if (length(idx)) nms[idx[1]] else nms[1]
      return(sf::st_read(gdbs[1], layer = lyr, quiet = quiet))
    }
  }

  shp <- list.files(dirpath, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(shp)) return(sf::st_read(shp[1], quiet = quiet))

  jsn <- list.files(dirpath, pattern = "\\.(geo)?json$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(jsn)) return(sf::st_read(jsn[1], quiet = quiet))

  kml <- list.files(dirpath, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(kml)) return(sf::st_read(kml[1], quiet = quiet))

  kmz <- list.files(dirpath, pattern = "\\.kmz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  if (length(kmz)) {
    td <- tempfile("kmz_")

    bf_unpack_archive(
      path = kmz[1],
      exdir = td,
      force_refresh = TRUE,
      quiet = quiet
    )

    kml2 <- list.files(td, pattern = "\\.kml$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (length(kml2)) return(sf::st_read(kml2[1], quiet = quiet))
  }

  NULL
}

#' Query one paged ArcGIS FeatureServer result as GeoJSON
#'
#' @param layer_url ArcGIS FeatureServer layer URL.
#' @param offset Result offset.
#' @param n Maximum records requested for the page.
#' @param quiet Logical. If `TRUE`, suppress download and read messages.
#'
#' @return An `sf` object for the page, or `NULL` when the page is empty or unreadable.
#'
#' @keywords internal
#' @noRd
.bf_arcgis_get_page <- function(layer_url, offset = 0L, n = 2000L, quiet = TRUE) {
  q <- paste0(
    layer_url,
    "/query?",
    "where=1%3D1",
    "&outFields=*",
    "&returnGeometry=true",
    "&outSR=4326",
    "&f=geojson",
    "&resultRecordCount=", n,
    "&resultOffset=", offset,
    "&orderByFields=FID"
  )
  tf <- tempfile(fileext = ".geojson")
  bf_download_file(q, tf, quiet = quiet)
  g <- try(suppressWarnings(sf::st_read(tf, quiet = quiet)), silent = TRUE)
  if (inherits(g, "try-error") || !inherits(g, "sf") || nrow(g) == 0L) return(NULL)
  g
}

#' Fetch all pages from an ArcGIS FeatureServer layer
#'
#' @param layer_url ArcGIS FeatureServer layer URL.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return A combined `sf` object, or `NULL` if no pages can be read.
#'
#' @keywords internal
#' @noRd
.bf_arcgis_fetch_all <- function(layer_url, quiet = TRUE) {
  pgsize <- 2000L
  out <- list()
  off <- 0L
  i <- 1L
  repeat {
    g <- .bf_arcgis_get_page(layer_url, offset = off, n = pgsize, quiet = quiet)
    if (is.null(g)) break
    out[[i]] <- g
    if (nrow(g) < pgsize) break
    off <- off + pgsize
    i <- i + 1L
  }
  if (!length(out)) return(NULL)
  x <- do.call(rbind, out)
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  x
}

#' Load global lake polygons (HydroLAKES)
#'
#' Downloads HydroLAKES once (RAW), standardizes key fields, then applies
#' optional filters in-memory on each call. Cached as a single GeoPackage.
#'
#' @param force_refresh logical; re-download even if cache exists.
#' @param quiet logical; suppress messages.
#' @param cache_dir Cache root used for the HydroLAKES cache. Must be supplied
#'   explicitly. In examples, tests and vignettes, use a path under `tempdir()`.
#' @param lakes_only logical; keep only natural lakes (exclude reservoirs/controls).
#' @param min_area_km2 optional numeric; drop lakes smaller than this area.
#'
#' @return An `sf` polygon object containing HydroLAKES features in WGS84
#'   longitude/latitude coordinates. The returned object includes standardised
#'   columns `lake_id`, `name`, `is_reservoir` and `area_km2`, alongside any
#'   retained source attributes. Each row represents a lake polygon. If
#'   `lakes_only` or `min_area_km2` is used, the returned object contains only
#'   the lakes retained by those filters.
#'
#' @family freshwater overlay loaders
#' @export
bf_load_lakes <- function(force_refresh = FALSE,
                          quiet = TRUE,
                          cache_dir = NULL,
                          lakes_only = TRUE,
                          min_area_km2 = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_lakes(). Please install it.")
  }

  HYDROL_SHP_ZIP <- "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip"

  cache_root <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_load_lakes()"
  )
  lakes_dir  <- file.path(cache_root, "lakes")
  if (!dir.exists(lakes_dir)) dir.create(lakes_dir, recursive = TRUE, showWarnings = FALSE)

  gpkg_raw <- file.path(lakes_dir, "hydrolakes_raw.gpkg")

  # ensure RAW cache
  if (!file.exists(gpkg_raw) || isTRUE(force_refresh)) {
    if (!quiet) message("biofetchR: Downloading HydroLAKES (first run; large file)...")
    tmp_zip <- tempfile(fileext = ".zip")
    bf_download_file(HYDROL_SHP_ZIP, tmp_zip, quiet = quiet)

    tmp_dir <- tempfile("hydrolakes_")

    bf_unpack_archive(
      path = tmp_zip,
      exdir = tmp_dir,
      force_refresh = TRUE,
      quiet = quiet
    )

    shp <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (!length(shp)) stop("HydroLAKES archive did not contain a shapefile.")
    x <- sf::st_read(shp[1], quiet = quiet)

    # Standardize key fields (per HydroLAKES schema)
    nm <- names(x)
    x$lake_id      <- if ("Hylak_id"  %in% nm) x$Hylak_id  else NA_integer_
    x$name         <- if ("Lake_name" %in% nm) x$Lake_name else NA_character_
    x$is_reservoir <- if ("Lake_type" %in% nm) x$Lake_type != 1L else NA
    x$area_km2     <- if ("Lake_area" %in% nm) x$Lake_area else NA_real_

    x <- sf::st_make_valid(x)
    if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326

    suppressWarnings(sf::st_write(x, gpkg_raw, delete_dsn = TRUE, quiet = TRUE))
    if (!quiet) message("biofetchR: HydroLAKES cached at: ", gpkg_raw)
  }

  # read RAW cache, then filter in-memory
  x <- sf::st_read(gpkg_raw, quiet = quiet)

  # apply filters without altering the RAW cache
  nm <- names(x)
  if (isTRUE(lakes_only) && "Lake_type" %in% nm) {
    x <- x[is.na(x$Lake_type) | x$Lake_type == 1L, , drop = FALSE]
  }
  if (!is.null(min_area_km2) && "Lake_area" %in% nm) {
    x <- x[is.na(x$Lake_area) | x$Lake_area >= min_area_km2, , drop = FALSE]
  }

  # return standardized + original attrs
  std_cols <- c("lake_id","name","is_reservoir","area_km2")
  std_cols <- std_cols[std_cols %in% names(x)]
  x <- x[, unique(c(std_cols, setdiff(names(x), std_cols))), drop = FALSE]
  sf::st_set_agr(x, "constant")
  x
}

#' Load global river reaches (HydroRIVERS) with optional caching
#'
#' Downloads + standardizes HydroRIVERS. If `cache="disk"`, writes an RDS cache
#' (and optionally a GPKG). If `cache="memory"`, returns in-memory only.
#' Filters (`min_strahler`, `min_discharge_cms`) are applied on read, not baked into cache.
#'
#' @param force_refresh logical; re-download even if cache exists (disk mode).
#' @param quiet logical; suppress messages.
#' @param cache_dir Cache root used when `cache = "disk"`. Must be supplied
#'   explicitly for disk caching. Ignored when `cache = "memory"`. In examples,
#'   tests and vignettes, use a path under `tempdir()`.
#' @param cache "disk" or "memory". Disk persists a cache; memory does not write files.
#' @param cache_format "rds","gpkg","both" (used only when cache="disk"; default "rds").
#' @param regions Character vector. Supported values include \code{"global"}, \code{"af"}, \code{"ar"}, \code{"as"}, \code{"au"}, \code{"eu"}, \code{"gr"}, \code{"na"}, \code{"sa"}, and \code{"si"}. If \code{NULL}, \code{"global"} is used.
#' @param min_strahler optional integer; keep ORD_STRA >= this (applied on read).
#' @param min_discharge_cms optional numeric; keep DIS_AV_CMS >= this (applied on read).
#'
#' @return An `sf` line object containing HydroRIVERS reaches in WGS84
#'   longitude/latitude coordinates. The returned object includes canonical
#'   reach attributes `hyriv_id`, `ord_stra`, `ord_clas`, `ord_flow`,
#'   `dis_av_cms`, `length_km` and `hybas_l12`. Each row represents a river
#'   reach. If `min_strahler` or `min_discharge_cms` is supplied, the returned
#'   object contains only reaches retained by those filters.
#'
#' @family freshwater overlay loaders
#' @export
bf_load_rivers <- function(force_refresh = FALSE,
                           quiet = TRUE,
                           cache_dir = NULL,
                           cache = c("disk","memory"),
                           cache_format = c("rds","gpkg","both"),
                           regions = NULL,
                           min_strahler = NULL,
                           min_discharge_cms = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for bf_load_rivers(). Please install it.")
  }
  cache <- match.arg(cache)
  cache_format <- match.arg(cache_format)

  base <- "https://data.hydrosheds.org/file/hydrorivers"
  pick_urls <- function(regs) {
    if (is.null(regs) || identical(regs, "global"))
      return(file.path(base, "HydroRIVERS_v10_shp.zip"))
    regs <- match.arg(regs, c("af","ar","as","au","eu","gr","na","sa","si"), several.ok = TRUE)
    file.path(base, paste0("HydroRIVERS_v10_", regs, "_shp.zip"))
  }

  has_gpkg <- function() {
    drv <- try(sf::st_drivers(), silent = TRUE)
    !inherits(drv, "try-error") && any(drv$name == "GPKG" & drv$write)
  }

  # build one big sf in memory
  build_in_memory <- function(regs) {
    urls <- pick_urls(regs)
    if (!quiet) message("biofetchR: Downloading HydroRIVERS (", length(urls), " file(s)) ...")
    layers <- vector("list", length(urls))
    for (i in seq_along(urls)) {
      u <- urls[i]
      tmp_zip <- tempfile(fileext = ".zip")
      bf_download_file(u, tmp_zip, quiet = quiet)

      tmp_dir <- tempfile("hydrorivers_")

      bf_unpack_archive(
        path = tmp_zip,
        exdir = tmp_dir,
        force_refresh = TRUE,
        quiet = quiet
      )

      shp <- list.files(tmp_dir, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
      if (!length(shp)) stop("HydroRIVERS archive did not contain a shapefile: ", u)
      layers[[i]] <- sf::st_read(shp[1], quiet = quiet)
    }
    x <- do.call(rbind, layers)

    # pure/local renamer
    ren <- function(obj, old, new) {
      if (old %in% names(obj)) {
        nm <- names(obj); nm[nm == old] <- new; names(obj) <- nm
      }
      obj
    }
    x <- ren(x, "HYRIV_ID","hyriv_id")
    x <- ren(x, "ORD_STRA","ord_stra")
    x <- ren(x, "ORD_CLAS","ord_clas")
    x <- ren(x, "ORD_FLOW","ord_flow")
    x <- ren(x, "DIS_AV_CMS","dis_av_cms")
    x <- ren(x, "LENGTH_KM","length_km")
    x <- ren(x, "HYBAS_L12","hybas_l12")

    for (nm in c("hyriv_id","ord_stra","ord_clas","ord_flow","dis_av_cms","length_km","hybas_l12")) {
      if (!nm %in% names(x)) x[[nm]] <- NA
    }
    x <- x[, !duplicated(tolower(names(x))), drop = FALSE]

    x <- suppressWarnings(sf::st_make_valid(x))
    if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
    sf::st_set_agr(x, "constant")
    x
  }

  if (cache == "memory") {
    x <- build_in_memory(regions)
    if (!is.null(min_strahler) && "ord_stra" %in% names(x))
      x <- x[is.na(x$ord_stra) | x$ord_stra >= min_strahler, , drop = FALSE]
    if (!is.null(min_discharge_cms) && "dis_av_cms" %in% names(x))
      x <- x[is.na(x$dis_av_cms) | x$dis_av_cms >= min_discharge_cms, , drop = FALSE]
    return(x)
  }


  # disk caching
  cache_root <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_load_rivers(cache = 'disk')"
  )
  rivers_dir <- file.path(cache_root, "rivers"); if (!dir.exists(rivers_dir)) dir.create(rivers_dir, TRUE, FALSE)
  rds_raw    <- file.path(rivers_dir, "hydrorivers_raw.rds")
  gpkg_raw   <- file.path(rivers_dir, "hydrorivers_raw.gpkg")

  write_cache <- function(x) {
    if (cache_format %in% c("rds","both")) {
      tmp_rds <- tempfile(fileext = ".rds"); saveRDS(x, tmp_rds)
      if (file.exists(rds_raw)) unlink(rds_raw)
      file.rename(tmp_rds, rds_raw)
    }
    if (cache_format %in% c("gpkg","both") && has_gpkg()) {
      tmp_gpkg <- tempfile(fileext = ".gpkg")
      ok <- try({
        suppressWarnings(sf::st_write(x, tmp_gpkg, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE))
        TRUE
      }, silent = TRUE)
      if (isTRUE(ok)) {
        if (file.exists(gpkg_raw)) unlink(gpkg_raw)
        file.rename(tmp_gpkg, gpkg_raw)
      } else {
        if (!quiet) message("biofetchR: GPKG write failed; continuing with RDS cache.")
        if (file.exists(tmp_gpkg)) unlink(tmp_gpkg)
      }
    }
  }

  read_cached <- function() {
    if (file.exists(rds_raw)) {
      x <- try(readRDS(rds_raw), silent = TRUE)
      if (!inherits(x, "try-error") && inherits(x, "sf")) return(x)
    }
    if (file.exists(gpkg_raw)) {
      x <- try(sf::st_read(gpkg_raw, quiet = quiet), silent = TRUE)
      if (!inherits(x, "try-error") && inherits(x, "sf")) return(x)
    }
    NULL
  }

  corrupt_file <- function(p) file.exists(p) && isTRUE(tryCatch(file.info(p)$size < 1024, error = function(e) FALSE))
  need_build <- isTRUE(force_refresh) ||
    (is.null(read_cached()) && !(file.exists(rds_raw) || file.exists(gpkg_raw))) ||
    corrupt_file(gpkg_raw)

  if (need_build) {
    x <- build_in_memory(regions)
    write_cache(x)
  } else {
    x <- read_cached()
    if (is.null(x)) {
      if (!quiet) message("biofetchR: Cache unreadable; rebuilding ...")
      x <- build_in_memory(regions); write_cache(x)
    }
  }

  if (!is.null(min_strahler) && "ord_stra" %in% names(x))
    x <- x[is.na(x$ord_stra) | x$ord_stra >= min_strahler, , drop = FALSE]
  if (!is.null(min_discharge_cms) && "dis_av_cms" %in% names(x))
    x <- x[is.na(x$dis_av_cms) | x$dis_av_cms >= min_discharge_cms, , drop = FALSE]

  x
}

#' Annotate HydroRIVERS reaches with OSM river names
#' @param rivers sf LINESTRING with field `hyriv_id` (from bf_load_rivers()).
#' @param iso2c character ISO2 country codes used to build a bbox query.
#' @param max_snap_m numeric, max distance (m) to accept a name match.
#' @param quiet logical.
#'
#' @return An `sf` line object with the same river reaches as `rivers`, returned
#'   in WGS84 longitude/latitude coordinates, with an added character column
#'   `river_name`. The `river_name` column contains the nearest accepted
#'   OpenStreetMap waterway name within `max_snap_m`, or `NA_character_` when no
#'   suitable named waterway is found. The output is used to add human-readable
#'   river labels to HydroRIVERS features.
#'
#' @family freshwater overlay loaders
#' @export
bf_name_rivers_osm <- function(rivers, iso2c, max_snap_m = 500, quiet = TRUE) {
  stopifnot(inherits(rivers, "sf"), "hyriv_id" %in% names(rivers))
  if (!requireNamespace("osmdata", quietly = TRUE) ||
      !requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Packages 'osmdata' and 'rnaturalearth' are required.")
  }

  # 1) Country bbox (WGS84)
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  sel   <- world[world$iso_a2 %in% toupper(iso2c), , drop = FALSE]
  if (!nrow(sel)) stop("No countries found for: ", paste(iso2c, collapse=", "))
  bb    <- sf::st_bbox(sf::st_transform(sel, 4326))

  # 2) OSM query for named waterways
  q <- osmdata::opq(bbox = bb) |>
    osmdata::add_osm_feature(key = "waterway",
                             value = c("river","stream","canal","drain")) |>
    osmdata::add_osm_feature(key = "name")
  osm <- osmdata::osmdata_sf(q)
  lines <- osm$osm_lines
  if (is.null(lines) || !nrow(lines)) {
    if (!quiet) message("No named OSM waterways returned for bbox.")
    rivers$river_name <- NA_character_
    return(rivers)
  }
  lines <- lines[!is.na(lines$name) & nzchar(lines$name), c("name","geometry")]
  lines <- sf::st_transform(lines, 4326)

  # 3) Snap by nearest feature with distance screen (meters)
  riv4326 <- sf::st_transform(rivers, 4326)
  idx <- sf::st_nearest_feature(riv4326, lines)

  # distance check in a metric CRS
  riv_m  <- tryCatch(sf::st_transform(riv4326, 3857), error = function(e) riv4326)
  line_m <- tryCatch(sf::st_transform(lines[idx, ], 3857), error = function(e) lines[idx, ])
  d <- suppressWarnings(sf::st_distance(sf::st_geometry(riv_m),
                                        sf::st_geometry(line_m),
                                        by_element = TRUE))
  nm <- lines$name[idx]
  nm[as.numeric(d) > max_snap_m] <- NA_character_

  riv4326$river_name <- nm
  riv4326
}


#' Load nested basins (HydroBASINS; Pfafstetter)
#'
#' Downloads selected regional tiles at a chosen Pfafstetter level (L1-L12),
#' caches once, and returns polygons with canonical fields:
#' hybas_id, next_down, next_sink, main_bas, sub_area, up_area, endo, order.
#'
#' @param level Integer from 1 to 12 giving the Pfafstetter basin level.
#' @param with_lakes Logical; if TRUE, use the customised HydroBASINS "with lakes" variant.
#' @param cache_dir Cache root used for downloaded HydroBASINS files. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param quiet Logical; if TRUE, suppress download and status messages.
#'
#' @return An `sf` polygon object containing HydroBASINS catchments for the
#'   requested Pfafstetter `level`, in WGS84 longitude/latitude coordinates.
#'   The returned object includes canonical basin attributes such as `hybas_id`,
#'   `next_down`, `next_sink`, `main_bas`, `sub_area`, `up_area`, `endo` and
#'   `order`, where available. Each row represents a basin polygon used for
#'   assigning occurrence records to freshwater catchment units.
#'
#' @family freshwater overlay loaders
#' @export
bf_load_basins <- function(level = 12, with_lakes = FALSE,
                           cache_dir = NULL,
                           quiet = TRUE) {
  stopifnot(requireNamespace("sf", quietly = TRUE))

  cache_dir <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_load_basins()"
  )

  hb_dir <- file.path(cache_dir, "hydrobasins")
  dir.create(hb_dir, recursive = TRUE, showWarnings = FALSE)

  regions <- getOption("biofetchR.basins.regions",
                       c("eu","na","sa","af","as","au","ar","an"))

  url_roots <- c(
    "https://data.hydrosheds.org/file/hydrobasins",
    "https://data.hydrosheds.org/file/hydrobasins/standard"
  )

  # candidate filenames (some mirrors use gdb.zip)
  fname_patterns <- c(
    sprintf("hybas_%%s_lev%02d_v1c.zip", level),
    sprintf("hybas_%%s_lev%02d_v1c.gdb.zip", level),
    sprintf("hybas_%%s_lev%02d_v1c_gdb.zip", level)
  )

  # downloader (quietly tries several URL variants)
  fetch_zip <- function(region) {
    dest <- file.path(hb_dir, sprintf("hybas_%s_lev%02d_v1c.zip", region, level))
    if (file.exists(dest)) return(dest)

    for (root in url_roots) {
      for (pat in fname_patterns) {
        fname <- sprintf(pat, region)
        url <- file.path(root, fname)
        tmp <- tempfile(fileext = ".zip")
        ok <- tryCatch({
          curl::curl_download(url, tmp, quiet = quiet)
          TRUE
        }, error = function(e) FALSE)
        if (ok) {
          file.rename(tmp, dest)
          if (!quiet) message("HydroBASINS: downloaded ", url)
          return(dest)
        }
      }
    }
    if (!quiet) message("HydroBASINS: FAILED to fetch region '", region,
                        "' for level ", level, " (tried ", paste(url_roots, collapse=" | "), ")")
    return(NA_character_)
  }

  zips <- vapply(regions, fetch_zip, character(1))
  zips <- zips[!is.na(zips) & file.exists(zips)]
  if (!length(zips)) stop("HydroBASINS download failed for all regions.")

  # unpack + read the target level shapefile (or gdb)
  read_one <- function(zip) {
    exdir <- tempfile("hybas_")

    bf_unpack_archive(
      path = zip,
      exdir = exdir,
      force_refresh = TRUE,
      quiet = quiet
    )

    # prefer levXX shp
    shp <- list.files(exdir, pattern = sprintf("lev%02d_v1c\\.shp$", level),
                      recursive = TRUE, full.names = TRUE)
    if (!length(shp)) {
      # any shapefile as fallback
      shp <- list.files(exdir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
    }
    if (length(shp)) {
      sf <- suppressMessages(sf::st_read(shp[1], quiet = TRUE))
    } else {
      # try geodatabase
      gdb <- list.files(exdir, pattern = "\\.gdb$", recursive = TRUE, full.names = TRUE)
      if (!length(gdb)) return(NULL)
      # layer name usually includes level; read first layer if unknown
      layers <- suppressWarnings(sf::st_layers(gdb[1])$name)
      lyr <- if (length(grep(sprintf("lev%02d", level), layers))) {
        layers[grep(sprintf("lev%02d", level), layers)][1]
      } else layers[1]
      sf <- suppressMessages(sf::st_read(gdb[1], layer = lyr, quiet = TRUE))
    }
    # normalize id column
    if (!"hybas_id" %in% names(sf)) {
      if ("HYBAS_ID" %in% names(sf)) sf$hybas_id <- sf$HYBAS_ID
      else if ("HYBAS_ID" %in% toupper(names(sf))) {
        nm <- names(sf)[toupper(names(sf)) == "HYBAS_ID"][1]
        sf$hybas_id <- sf[[nm]]
      } else {
        sf$hybas_id <- NA_integer_
      }
    }
    # keep only id + geometry to stay light
    geom_col <- attr(sf, "sf_column") %||% "geometry"
    sf[, c("hybas_id", geom_col)]
  }

  sfs <- Filter(Negate(is.null), lapply(zips, read_one))
  if (!length(sfs)) stop("HydroBASINS: unzip/read produced no layers.")
  out <- do.call(rbind, sfs)
  # make valid just in case
  out <- tryCatch(suppressWarnings(sf::st_make_valid(out)), error = function(e) out)
  out
}

#' GMBA hierarchy-aware point-in-polygon join (0-360 degrees safe)
#'
#' For each input point, finds all intersecting GMBA Mountains v2 polygons and
#' returns the **single most specific** polygon label per point. "Specificity"
#' is chosen by either the deepest hierarchy level (preferred) or the smallest
#' polygon area (fallback).
#'
#' @param pts_sf  An `sf` POINT layer (any CRS).
#' @param gmba_sf An `sf` POLYGON/MULTIPOLYGON layer for GMBA v2 (any CRS). Must
#'   contain a label column named in `name_col` (default `"geoname_gmba"`).
#' @param name_col Character scalar; column in `gmba_sf` to return as the label
#'   (e.g., `"geoname_gmba"`, or a native GMBA name/ID column you prepared).
#' @param level_cols Character vector of candidate columns encoding hierarchical
#'   depth (e.g., `c("LEVEL","level","HIERARCHY","hierarchy")`). The first
#'   present column that can be coerced to numeric is used. Higher numbers are
#'   treated as deeper (more specific) levels.
#' @param prefer One of `c("deepest_level","smallest_area")`. If
#'   `"deepest_level"` is requested but no usable level column exists, the
#'   function falls back to `"smallest_area"`.
#' @param quiet Logical; suppress progress messages.
#'
#' @return A character vector of length `nrow(pts_sf)` with the chosen GMBA
#'   label for each point, or `NA` when a point intersects no GMBA polygon.
#'
#' @details
#' - **0-360 degrees safety:** If `gmba_sf` spans 0-360 degrees longitudes (xmin >= 0 & xmax > 180),
#'   points with negative longitudes are shifted by +360 degrees for the spatial join,
#'   mirroring the behavior used elsewhere in the pipeline.
#' - **Tie-breaking:** When multiple polygons intersect a point:
#'   1) Use the polygon with the **deepest hierarchy level** (largest numeric
#'      value in the chosen `level_cols`).
#'   2) If levels are unavailable or tied, pick the **smallest area** polygon
#'      (computed in a projected CRS when possible).
#'
#' @keywords internal
#' @noRd
.join_gmba_labels <- function(
    pts_sf, gmba_sf,
    name_col   = "geoname_gmba",
    level_cols = c("LEVEL","level","HIERARCHY","hierarchy"),
    prefer     = c("deepest_level","smallest_area"),
    quiet      = FALSE
) {
  stopifnot(inherits(pts_sf, "sf"), inherits(gmba_sf, "sf"), name_col %in% names(gmba_sf))
  prefer <- match.arg(prefer)

  # Work in WGS84
  pts4326  <- if (is.na(sf::st_crs(pts_sf))  || sf::st_crs(pts_sf)  != sf::st_crs(4326)) sf::st_transform(pts_sf,  4326) else pts_sf
  gmba4326 <- if (is.na(sf::st_crs(gmba_sf)) || sf::st_crs(gmba_sf) != sf::st_crs(4326)) sf::st_transform(gmba_sf, 4326) else gmba_sf

  .is_0360 <- function(x) {
    bb <- sf::st_bbox(x)
    isTRUE(!is.na(bb["xmin"]) && !is.na(bb["xmax"]) && bb["xmin"] >= 0 && bb["xmax"] > 180)
  }

  if (!.is_0360(gmba4326)) {
    pts_for_join <- pts4326
  } else {
    cc <- sf::st_coordinates(pts4326)
    lon360 <- ifelse(cc[,1] < 0, cc[,1] + 360, cc[,1])
    lat    <- cc[,2]
    tmp <- cbind(sf::st_drop_geometry(pts4326), ..lon360 = lon360, ..lat = lat)
    pts_for_join <- sf::st_as_sf(tmp, coords = c("..lon360","..lat"), crs = 4326, remove = TRUE)
    if (!quiet) message("GMBA join using 0-360 degrees point shift")
  }

  # Specificity metrics
  lvl_col <- intersect(level_cols, names(gmba4326))
  lvl_vec <- if (length(lvl_col)) suppressWarnings(as.numeric(gmba4326[[lvl_col[1]]])) else rep(NA_real_, nrow(gmba4326))

  gmba_m <- tryCatch(sf::st_transform(gmba4326, 3857), error = function(e) gmba4326)
  areas  <- tryCatch(as.numeric(sf::st_area(gmba_m)),    error = function(e) rep(NA_real_, nrow(gmba4326)))

  hits <- sf::st_intersects(pts_for_join, gmba4326)

  pick_one <- function(cands) {
    if (!length(cands)) return(NA_integer_)
    if (prefer == "deepest_level" && any(!is.na(lvl_vec[cands]))) {
      cands <- cands[order(-lvl_vec[cands], areas[cands], na.last = TRUE)]
      return(cands[1])
    }
    cands[order(areas[cands], na.last = TRUE)][1]
  }

  sel_idx <- vapply(hits, pick_one, integer(1))
  out <- rep(NA_character_, nrow(pts_sf))
  ok  <- !is.na(sel_idx)
  out[ok] <- as.character(gmba4326[[name_col]][sel_idx[ok]])

  if (!quiet) message("GMBA join (most-specific): ", sum(ok), " / ", length(out), " labeled")
  out
}

# -----------------------------------------------------------------------------
# Natural Earth + RESOLVE Ecoregions 2017 loaders (auto-download + cache)
# -----------------------------------------------------------------------------

#' Get a writable cache directory for biofetchR spatial layers
#'
#' This helper returns a cache directory (creating it if needed) where biofetchR
#' can store downloaded spatial layers (zips + extracted shapefiles). Users are
#' not expected to download or manage raw data files themselves.
#'
#' @param cache_dir Character. Root cache directory. Must be supplied explicitly.
#'   In examples, tests and vignettes, use a path under `tempdir()`.
#' @param subdir Optional character. Subdirectory to create/use inside cache_dir.
#'
#' @return A character vector of length one giving the normalised path to a
#'   writable cache directory. The directory is created if it does not already
#'   exist, so the function is also called for the side effect of preparing a
#'   cache location for downloaded spatial layers. If `subdir` is supplied, the
#'   returned path points to that subdirectory inside the cache root.
#'
#' @family spatial cache helpers
#' @export
bf_cache_dir <- function(cache_dir = NULL, subdir = NULL) {
  cache_dir <- .bf_require_explicit_spatial_cache(
    cache_dir,
    context = "bf_cache_dir()"
  )

  if (!is.null(subdir) && nzchar(subdir)) {
    cache_dir <- file.path(cache_dir, subdir)

    if (!dir.exists(cache_dir)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
}

#' Find a shapefile within a directory
#'
#' @param xdir Character. Directory to search.
#' @param prefer Optional character vector. If provided, shapefiles whose basename
#'   matches any of these patterns are preferred.
#'
#' @return A single shapefile path (character length 1).
#' @keywords internal
#' @noRd
bf_find_shp <- function(xdir, prefer = NULL) {
  shp <- list.files(xdir, pattern = "\\.shp$", recursive = TRUE, full.names = TRUE)
  if (!length(shp)) stop("No .shp found in: ", xdir)

  if (!is.null(prefer) && length(prefer)) {
    b <- basename(shp)
    for (p in prefer) {
      hit <- grep(p, b, ignore.case = TRUE)
      if (length(hit)) return(shp[hit[1]])
    }
  }
  shp[1]
}

#' Load Natural Earth Urban Areas (polygons), auto-downloaded and cached
#'
#' This downloads the Natural Earth "Urban Areas" polygons (10m cultural vectors),
#' caches the zip, extracts it, and returns an sf object in EPSG:4326 with a
#' consistent ID + label column for downstream point-tagging.
#'
#' Users do not need to pre-download shapefiles; biofetchR handles retrieval and caching.
#'
#' @param cache_dir Character. Cache directory root. Must be supplied explicitly.
#'   In examples, tests and vignettes, use a path under `tempdir()`.
#' @param url Character. Optional override download URL.
#' @param force_refresh Logical. Re-download and re-extract even if cached.
#' @param quiet Logical. Reduce messages.
#'
#' @return An `sf` polygon object containing Natural Earth urban-area features
#'   in WGS84 longitude/latitude coordinates. The returned object contains
#'   `urban_id`, a character identifier for each feature; `urban_name`, a
#'   character label for the urban area; and the active geometry column. Each
#'   row represents one urban-area polygon or multipolygon used for assigning
#'   occurrence records to urban contextual overlays.
#'
#' @family terrestrial overlay loaders
#' @export
bf_load_ne_urban_areas <- function(cache_dir = NULL,
                                   url = "http://www.naturalearthdata.com/download/10m/cultural/ne_10m_urban_areas.zip",
                                   force_refresh = FALSE,
                                   quiet = TRUE) {

  cache_dir <- bf_cache_dir(cache_dir, "naturalearth_urban_10m")
  zipfile   <- file.path(cache_dir, "ne_10m_urban_areas.zip")
  exdir     <- file.path(cache_dir, "unzipped")

  bf_download_cached(
    url = url,
    dest = zipfile,
    force_refresh = force_refresh,
    quiet = quiet
  )

  bf_unzip_cached(zipfile, exdir, force_refresh = force_refresh)

  shp <- bf_find_shp(exdir, prefer = c("ne_10m_urban_areas"))
  x <- sf::st_read(shp, quiet = quiet)
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  if (sf::st_crs(x) != sf::st_crs(4326)) x <- sf::st_transform(x, 4326)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]

  nms <- names(x); ln <- tolower(nms)
  name_col <- NULL
  for (cand in c("name", "name_en", "namealt", "name_alt")) {
    hit <- which(ln == cand)
    if (length(hit)) { name_col <- nms[hit[1]]; break }
  }
  if (is.null(name_col)) {
    x$urban_name <- paste0("URBAN_", seq_len(nrow(x)))
  } else {
    x$urban_name <- as.character(x[[name_col]])
    x$urban_name[!nzchar(x$urban_name)] <- paste0("URBAN_", which(!nzchar(x$urban_name)))
  }

  id_col <- NULL
  for (cand in c("ne_id", "fid", "id")) {
    hit <- which(ln == cand)
    if (length(hit)) { id_col <- nms[hit[1]]; break }
  }
  x$urban_id <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))

  gcol <- attr(x, "sf_column"); if (is.null(gcol) || !gcol %in% names(x)) gcol <- "geometry"
  x <- x[, c("urban_id", "urban_name", gcol), drop = FALSE]
  x
}

#' Load RESOLVE Ecoregions 2017 (polygons), auto-downloaded and cached
#'
#' Downloads the RESOLVE terrestrial ecoregions dataset (2017 update), caches
#' and extracts it, and returns an sf object in EPSG:4326 with consistent ID +
#' label columns.
#'
#' Users do not need to pre-download shapefiles; biofetchR handles retrieval and caching.
#'
#' @param cache_dir Character. Cache directory root. Must be supplied explicitly.
#'   In examples, tests and vignettes, use a path under `tempdir()`.
#' @param url Character. Optional override download URL.
#' @param force_refresh Logical. Re-download and re-extract even if cached.
#' @param quiet Logical. Reduce messages.
#'
#' @return An `sf` polygon object containing RESOLVE 2017 terrestrial ecoregion
#'   features in WGS84 longitude/latitude coordinates. The returned object
#'   contains `resolve_eco_id`, a character ecoregion identifier;
#'   `resolve_eco_name`, a character ecoregion label; `resolve_biome`, the
#'   biome name where available; `resolve_realm`, the realm name where
#'   available; and the active geometry column. Each row represents one
#'   ecoregion polygon or multipolygon used for terrestrial contextual
#'   assignment.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   resolve <- bf_load_resolve_ecoregions2017(
#'     cache_dir = file.path(tempdir(), "biofetchR_resolve_ecoregions")
#'   )
#'
#'   names(resolve)
#' }
#' }
#'
#' @family terrestrial overlay loaders
#' @export
bf_load_resolve_ecoregions2017 <- function(cache_dir = NULL,
                                           url = "https://storage.googleapis.com/teow2016/Ecoregions2017.zip",
                                           force_refresh = FALSE,
                                           quiet = TRUE) {

  cache_dir <- bf_cache_dir(cache_dir, "resolve_ecoregions_2017")
  zipfile   <- file.path(cache_dir, "Ecoregions2017.zip")
  exdir     <- file.path(cache_dir, "unzipped")

  bf_download_cached(
    url = url,
    dest = zipfile,
    force_refresh = force_refresh,
    quiet = quiet
  )

  bf_unzip_cached(zipfile, exdir, force_refresh = force_refresh)

  shp <- bf_find_shp(exdir, prefer = c("Ecoregions2017", "ecoregions2017"))
  x <- sf::st_read(shp, quiet = quiet)
  if (is.na(sf::st_crs(x))) sf::st_crs(x) <- 4326
  if (sf::st_crs(x) != sf::st_crs(4326)) x <- sf::st_transform(x, 4326)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]

  nms <- names(x); ln <- tolower(nms)

  # ECO_ID-like
  id_col <- NULL
  for (cand in c("eco_id", "eco_id_u", "ecoid", "objectid", "fid")) {
    hit <- which(ln == cand)
    if (length(hit)) { id_col <- nms[hit[1]]; break }
  }
  x$resolve_eco_id <- if (!is.null(id_col)) as.character(x[[id_col]]) else as.character(seq_len(nrow(x)))

  # ECO_NAME-like
  name_col <- NULL
  for (cand in c("eco_name", "ecoregion", "ecoregionname", "name")) {
    hit <- which(ln == cand)
    if (length(hit)) { name_col <- nms[hit[1]]; break }
  }
  x$resolve_eco_name <- if (!is.null(name_col)) as.character(x[[name_col]]) else paste0("ECO_", x$resolve_eco_id)

  # BIOME / REALM (best effort)
  biome_col <- NULL
  for (cand in c("biome_name", "biome", "biome_nm")) {
    hit <- which(ln == cand)
    if (length(hit)) { biome_col <- nms[hit[1]]; break }
  }
  realm_col <- NULL
  for (cand in c("realm", "realm_name")) {
    hit <- which(ln == cand)
    if (length(hit)) { realm_col <- nms[hit[1]]; break }
  }

  x$resolve_biome <- if (!is.null(biome_col)) as.character(x[[biome_col]]) else NA_character_
  x$resolve_realm <- if (!is.null(realm_col)) as.character(x[[realm_col]]) else NA_character_

  gcol <- attr(x, "sf_column"); if (is.null(gcol) || !gcol %in% names(x)) gcol <- "geometry"
  x <- x[, c("resolve_eco_id", "resolve_eco_name", "resolve_biome", "resolve_realm", gcol), drop = FALSE]
  x
}
