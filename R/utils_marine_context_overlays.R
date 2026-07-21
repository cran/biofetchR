################################################################################
# utils_marine_context_overlays.R
# -----------------------------------------------------------------------------
# biofetchR - package-managed marine overlay loaders
# -----------------------------------------------------------------------------
#
# PURPOSE
#   This file contains the retained marine overlay-loading utilities used by the
#   marine GBIF processing workflow. These helpers resolve, download/cache, read
#   and standardise polygon layers from Marine Regions through mregions2.
#
#   The resulting overlay objects are used by `process_gbif_marine_pipeline()` to
#   assign GBIF occurrence points to marine spatial units such as EEZs, MEOW
#   ecoregions, LME regions, IHO seas, Longhurst provinces and related Marine
#   Regions products.
#
# DESIGN PRINCIPLES
#   1. Keep marine overlays package-managed wherever possible. Users should not
#      need to manually download shapefiles for the supported overlays.
#   2. Resolve overlay aliases to one canonical biofetchR name. This keeps the
#      public API stable even when Marine Regions product identifiers differ
#      between mregions2 or Marine Regions releases.
#   3. Cache downloaded products as RDS files after first load. This avoids
#      repeated network calls during normal use and during release-gate checks.
#   4. Standardise every overlay to canonical columns used downstream:
#        - `marine_region_id`
#        - `marine_region_name`
#        - `marine_region_source`
#        - `marine_regions_layer`
#        - `geoname`
#   5. Fail clearly when an overlay cannot be resolved or loaded. Silent fallbacks
#      are avoided because they can otherwise produce misleading zero-row outputs.
#
# SUPPORTED PACKAGE-MANAGED OVERLAYS
#   Core overlays:
#      eez, lme, meow, iho
#
#   Extended overlays:
#      high_seas, territorial_seas, contiguous_zone, internal_waters,
#      archipelagic_waters, ecs, goas, eez_iho, eez_land, longhurst,
#      worldheritagemarine
#
# SCOPE
#   This file is about spatial *assignment*. The functions return
#   standardised `sf` objects for use in processing and CSV export workflows.
#   Boundary-only or line-only products are intentionally excluded because they
#   require nearest-line logic rather than polygon containment.
#
# DEPENDENCIES
#   The exported loader requires sf and mregions2. Internal utilities avoid
#   heavy dependencies where possible and repair UTF-8 text to reduce
#   Windows/locale-related failures.
#
# DATA USE AND ATTRIBUTION
#   Marine overlay products loaded through this file are external Marine Regions
#   data products accessed through mregions2. biofetchR standardises and caches
#   these layers for reproducible processing, but it does not own the underlying
#   datasets or replace their citation/licence requirements. Users should cite
#   Marine Regions and, where relevant, the specific Marine Regions data product
#   reported in `mregions2::mrp_list` / the `marine_regions_layer` output field.
#
################################################################################

#' Normalise marine overlay names and aliases
#'
#' Produces lower-case, underscore-separated keys for robust matching of overlay
#' names, aliases and Marine Regions product identifiers.
#'
#' @param x Character vector.
#'
#' @return Normalised character vector.
#'
#' @keywords internal
#' @noRd
.bf_marine_norm <- function(x) {
  x <- bf_repair_utf8_chr(x)
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "_", x)
}

#' Pick the first matching column from candidate names
#'
#' Searches an object for a column matching any supplied candidate, ignoring case
#' and punctuation differences. This supports Marine Regions products whose
#' attribute names differ across products or releases.
#'
#' @param x Data frame or `sf` object.
#' @param candidates Character vector of possible column names.
#'
#' @return Matching column name, or `NULL` if no candidate is present.
#'
#' @keywords internal
#' @noRd
.bf_marine_pick_col <- function(x, candidates) {
  nms <- names(x)
  if (!length(nms)) return(NULL)
  nms_norm <- .bf_marine_norm(nms)
  cand_norm <- .bf_marine_norm(candidates)
  hit <- match(cand_norm, nms_norm, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) nms[hit[[1]]] else NULL
}

#' Convert marine overlay attributes to safe character values
#'
#' Converts common non-character attribute classes to character and repairs UTF-8
#' encoding. This is used when standardising IDs and names from source layers.
#'
#' @param x Vector to convert.
#'
#' @return UTF-8 character vector.
#'
#' @keywords internal
#' @noRd
.bf_marine_as_chr <- function(x) {
  if (is.factor(x)) return(bf_repair_utf8_chr(as.character(x)))
  if (inherits(x, c("integer64", "Date", "POSIXct", "POSIXt"))) {
    return(bf_repair_utf8_chr(as.character(x)))
  }
  bf_repair_utf8_chr(as.character(x))
}

# -----------------------------------------------------------------------------
# Overlay registry
# -----------------------------------------------------------------------------
# The first candidate is the preferred mregions2 product identifier. Additional
# candidates are kept because Marine Regions product names occasionally differ
# between legacy and current package releases.
# -----------------------------------------------------------------------------
#' Marine overlay registry
#'
#' Defines the canonical marine overlay names recognised by biofetchR, accepted
#' aliases for each overlay and a short description of each spatial unit type.
#'
#' @return Data frame with columns `overlay`, `candidates` and `description`.
#'
#' @keywords internal
#' @noRd
bf_marine_overlay_registry <- function() {
  data.frame(
    overlay = c(
      "eez",
      "lme",
      "meow",
      "iho",
      "high_seas",
      "territorial_seas",
      "contiguous_zone",
      "internal_waters",
      "archipelagic_waters",
      "ecs",
      "goas",
      "eez_iho",
      "eez_land",
      "longhurst",
      "worldheritagemarine"
    ),
    candidates = I(list(
      c("eez"),
      c("lme"),
      c("meow", "ecoregions", "marine_ecoregions_of_the_world"),
      c("iho", "iho_seas"),
      c("high_seas", "highseas"),
      c("territorial_seas", "territorial_sea", "eez_12nm", "12nm"),
      c("contiguous_zone", "contiguous_zones", "eez_24nm", "24nm"),
      c("internal_waters", "internalwaters"),
      c("archipelagic_waters", "archipelagicwaters"),
      c("ecs", "extended_continental_shelves", "extended_continental_shelf"),
      c("goas", "global_oceans_and_seas", "global_ocean_and_seas"),
      c("eez_iho", "eez_iho_intersection", "eez_iho_seas"),
      c("eez_land", "eez_land_union", "eez_land"),
      c("longhurst", "longhurst_provinces", "longhurst_biogeochemical_provinces"),
      c("worldheritagemarine", "world_heritage_marine", "unesco_world_heritage_marine")
    )),
    description = c(
      "Exclusive Economic Zones",
      "Large Marine Ecosystems",
      "Marine Ecoregions of the World",
      "IHO Sea Areas",
      "High seas",
      "Territorial seas",
      "Contiguous zones",
      "Internal waters",
      "Archipelagic waters",
      "Extended continental shelves",
      "Global Oceans and Seas",
      "EEZ by IHO intersection",
      "EEZ-land union / marine-land zones",
      "Longhurst biogeochemical provinces",
      "UNESCO World Heritage marine sites"
    ),
    stringsAsFactors = FALSE
  )
}

#' Return supported marine overlay names
#'
#' Convenience helper used internally to retrieve canonical marine overlay names,
#' optionally including aliases that can be resolved to those canonical names.
#'
#' @param include_aliases Logical. If `TRUE`, return canonical names and aliases.
#'
#' @return Character vector of overlay names.
#'
#' @keywords internal
#' @noRd
bf_marine_overlay_sources <- function(include_aliases = FALSE) {
  registry <- bf_marine_overlay_registry()

  if (!isTRUE(include_aliases)) {
    return(registry$overlay)
  }

  unique(unlist(Map(c, registry$overlay, registry$candidates), use.names = FALSE))
}

#' List package-supported marine overlays
#'
#' Returns the marine spatial overlays that biofetchR can load through its
#' package-managed Marine Regions workflow.
#'
#' @details
#' By default, this helper returns the core overlays that are expected to be
#' resolvable in most mregions2/Marine Regions installations: `eez`, `meow`,
#' `lme` and `iho`. Set `include_extended = TRUE` to list the full registry of
#' supported extended overlays used by the marine pipeline and release-gate
#' checks.
#'
#' The returned overlay names can be supplied to `process_gbif_marine_pipeline()`
#' through its marine overlay argument. Aliases can also be inspected with
#' `include_aliases = TRUE`.
#'
#' @param include_aliases Logical. If `TRUE`, include a list-column containing
#'   accepted aliases and candidate Marine Regions product identifiers.
#' @param include_extended Logical. If `TRUE`, include all registered overlays.
#'   If `FALSE`, return only the core overlays.
#'
#' @return A data frame listing marine overlays registered by biofetchR. The
#'   returned object contains `overlay`, a character canonical overlay name;
#'   `description`, a character description of the spatial unit type; and
#'   `workflow`, a character label identifying the overlay as part of the marine
#'   workflow. When `include_aliases = TRUE`, the returned data frame also
#'   contains `aliases`, a list-column of accepted aliases and candidate Marine
#'   Regions product identifiers. The table is used to inspect valid overlay
#'   names before loading or passing them to marine processing workflows.
#'
#' @section Scope:
#' This helper only lists overlay names registered by biofetchR. It does not
#' download or cache Marine Regions products. Data access happens in
#' [bf_load_marine_regions_overlay()].
#'
#' @family marine overlay loaders
#'
#'
#' @examples
#' bf_available_marine_overlays()
#' bf_available_marine_overlays(include_extended = TRUE)
#'
#' @md
#' @export
bf_available_marine_overlays <- function(include_aliases = FALSE,
                                         include_extended = FALSE) {
  registry <- bf_marine_overlay_registry()

  # The live release gate should default to overlays that are normally resolvable
  # from Marine Regions without extra external datasets. Extended aliases/layers
  # can still be requested explicitly with BIOFETCHR_MARINE_OVERLAYS or by
  # calling this helper with `include_extended = TRUE`.
  core <- c("eez", "meow", "lme", "iho")

  if (!isTRUE(include_extended)) {
    registry <- registry[registry$overlay %in% core, , drop = FALSE]
  }

  out <- registry[, c("overlay", "description"), drop = FALSE]
  out$workflow <- "marine"

  if (isTRUE(include_aliases)) {
    out$aliases <- registry$candidates
  }

  out
}

#' Resolve a marine overlay alias to its canonical name
#'
#' Converts a user-supplied overlay name or accepted alias to the canonical
#' biofetchR overlay name used internally.
#'
#' @details
#' Marine Regions product names can vary across data releases and package
#' versions. biofetchR therefore keeps a registry of accepted aliases. This
#' function performs the alias lookup and returns the canonical name used by the
#' marine pipeline and overlay loader.
#'
#' @param region_source Character. Marine overlay name or alias, for example
#'   `"eez"`, `"iho_seas"` or `"longhurst_provinces"`.
#'
#' @return A character vector of length one containing the canonical biofetchR
#'   marine overlay name corresponding to `region_source`. The returned value is
#'   one of the registered overlay identifiers used internally by the marine
#'   pipeline and by [bf_load_marine_regions_overlay()]. The function errors if
#'   `region_source` cannot be matched to a supported overlay or alias.
#'
#' @family marine overlay loaders
#'
#' @examples
#' bf_marine_overlay_canonical("eez")
#' bf_marine_overlay_canonical("iho_seas")
#'
#' @md
#' @export
bf_marine_overlay_canonical <- function(region_source) {
  registry <- bf_marine_overlay_registry()
  src <- .bf_marine_norm(region_source)

  for (i in seq_len(nrow(registry))) {
    keys <- unique(c(registry$overlay[[i]], registry$candidates[[i]]))
    if (src %in% .bf_marine_norm(keys)) {
      return(registry$overlay[[i]])
    }
  }

  stop(
    "Unsupported marine overlay: '", region_source, "'. Supported overlays are: ",
    paste(registry$overlay, collapse = ", "),
    call. = FALSE
  )
}

#' Return the registry row for one marine overlay
#'
#' Resolves a requested overlay or alias and returns the corresponding row from
#' the marine overlay registry.
#'
#' @param region_source Character. Overlay name or alias.
#'
#' @return One-row data frame from `bf_marine_overlay_registry()`.
#'
#' @keywords internal
#' @noRd
.bf_marine_registry_row <- function(region_source) {
  registry <- bf_marine_overlay_registry()
  canonical <- bf_marine_overlay_canonical(region_source)
  hit <- which(.bf_marine_norm(registry$overlay) == .bf_marine_norm(canonical))
  registry[hit[[1]], , drop = FALSE]
}

#' Retrieve the Marine Regions product list
#'
#' Accesses `mregions2::mrp_list` and converts it to a data frame for product
#' resolution.
#'
#' @return Data frame of Marine Regions products exposed by mregions2.
#'
#' @keywords internal
#' @noRd
.bf_marine_mrp_list <- function() {
  bf_require_packages("mregions2", context = "marine overlay helper")
  x <- tryCatch(mregions2::mrp_list, error = function(e) NULL)
  if (is.null(x)) {
    stop("Could not access mregions2::mrp_list.", call. = FALSE)
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

#' Extract possible product identifiers from a product table
#'
#' Pulls identifier-like columns from a Marine Regions product table so overlay
#' candidate names can be matched robustly.
#'
#' @param product_table Data frame returned by `.bf_marine_mrp_list()`.
#'
#' @return Character vector of possible product identifiers.
#'
#' @keywords internal
#' @noRd
.bf_marine_product_ids <- function(product_table) {
  id_cols <- intersect(
    c("layer", "name", "id", "product", "identifier", "Name", "Layer", "title", "Title"),
    names(product_table)
  )
  if (!length(id_cols)) {
    id_cols <- names(product_table)[vapply(product_table, is.character, logical(1))]
  }
  ids <- unique(unlist(product_table[id_cols], use.names = FALSE))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  as.character(ids)
}

#' Resolve a canonical overlay to a Marine Regions product
#'
#' Matches the requested overlay against product identifiers and text fields in
#' the available mregions2 product list. Falls back to the preferred candidate
#' when the product list cannot be queried.
#'
#' @param region_source Character. Overlay name or alias.
#' @param quiet Logical. If `FALSE`, print resolution messages.
#'
#' @return Character string giving the selected Marine Regions product/layer id.
#'
#' @keywords internal
#' @noRd
.bf_resolve_marine_product <- function(region_source, quiet = TRUE) {
  row <- .bf_marine_registry_row(region_source)
  candidates <- row$candidates[[1]]

  product_table <- tryCatch(.bf_marine_mrp_list(), error = function(e) NULL)
  if (is.null(product_table) || !nrow(product_table)) {
    if (!isTRUE(quiet)) {
      message("Marine Regions product list unavailable; trying preferred product id: ", candidates[[1]])
    }
    return(candidates[[1]])
  }

  ids <- .bf_marine_product_ids(product_table)
  ids_norm <- .bf_marine_norm(ids)

  # Exact candidate match against product identifiers.
  for (cand in candidates) {
    hit <- which(ids_norm == .bf_marine_norm(cand))
    if (length(hit)) return(ids[[hit[[1]]]])
  }

  # Flexible keyword search across all character columns in mrp_list.
  text_cols <- names(product_table)[vapply(product_table, is.character, logical(1))]
  if (length(text_cols)) {
    hay <- apply(product_table[text_cols], 1, function(z) paste(z, collapse = " "))
    hay_norm <- .bf_marine_norm(hay)

    keywords <- unique(c(row$overlay, candidates))
    keywords <- keywords[nzchar(keywords)]
    keywords <- .bf_marine_norm(keywords)

    for (kw in keywords) {
      hit <- grep(kw, hay_norm, fixed = TRUE)
      if (length(hit)) {
        # Prefer an identifier-like column if available.
        id_col <- intersect(c("layer", "name", "id", "product", "identifier"), names(product_table))
        if (length(id_col)) {
          val <- product_table[[id_col[[1]]]][hit[[1]]]
          if (!is.na(val) && nzchar(as.character(val))) return(as.character(val))
        }
        return(ids[[hit[[1]]]])
      }
    }
  }

  stop(
    "Could not resolve marine overlay '", region_source,
    "' to a Marine Regions data product. Checked candidate product ids: ",
    paste(candidates, collapse = ", "),
    ". Run `mregions2::mrp_list` to inspect available products in your installed mregions2/Marine Regions setup.",
    call. = FALSE
  )
}

#' Load a Marine Regions product with an RDS cache
#'
#' Retrieves a Marine Regions product through mregions2, caches the loaded `sf`
#' object as an RDS file, and reuses the cache on later calls unless refresh is
#' requested.
#'
#' @param layer Character. Marine Regions product/layer identifier.
#' @param cache_dir Directory used for the overlay cache.
#' @param force_refresh Logical. If `TRUE`, ignore any existing RDS cache.
#' @param quiet Logical. If `FALSE`, print cache messages.
#'
#' @return An `sf` object returned by `mregions2::mrp_get()`.
#'
#' @keywords internal
#' @noRd
.bf_mrp_get_cached <- function(layer, cache_dir, force_refresh = FALSE, quiet = TRUE) {
  bf_require_packages(c("sf", "mregions2"), context = "marine overlay helper")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  safe_layer <- gsub("[^A-Za-z0-9_]+", "_", layer)
  rds_cache <- file.path(cache_dir, paste0("marine_regions_", safe_layer, ".rds"))

  if (file.exists(rds_cache) && !isTRUE(force_refresh)) {
    x <- tryCatch(readRDS(rds_cache), error = function(e) NULL)
    if (inherits(x, "sf") && nrow(x)) {
      if (!isTRUE(quiet)) message("Using cached Marine Regions overlay: ", rds_cache)
      return(x)
    }
  }

  old_path <- getOption("mregions2.download_path", default = NULL)
  options(mregions2.download_path = cache_dir)
  on.exit(options(mregions2.download_path = old_path), add = TRUE)

  x <- NULL
  last_err <- NULL

  # mregions2 1.1.x documents `mrp_get(layer, path = ...)`. Use positional
  # matching first because some installed builds do not accept a named `layer`
  # argument. Do not try invalid names such as `name = layer` because that masks
  # the real failure with a misleading final error.
  x <- tryCatch(
    mregions2::mrp_get(layer, path = cache_dir),
    error = function(e) {
      last_err <<- e
      NULL
    }
  )

  if (!inherits(x, "sf") || !nrow(x)) {
    x <- tryCatch(
      mregions2::mrp_get(layer),
      error = function(e) {
        last_err <<- e
        NULL
      }
    )
  }

  if (!inherits(x, "sf") || !nrow(x)) {
    msg <- if (inherits(last_err, "condition")) conditionMessage(last_err) else "unknown error"
    stop("Failed to load Marine Regions product '", layer, "': ", msg, call. = FALSE)
  }

  try(saveRDS(x, rds_cache), silent = TRUE)
  x
}

#' Standardise a marine overlay object
#'
#' Repairs text encodings, ensures WGS84 coordinates, validates geometries and
#' adds canonical marine-region columns required by downstream biofetchR
#' workflows.
#'
#' @param x `sf` object loaded from Marine Regions.
#' @param region_source Canonical overlay name.
#' @param marine_regions_layer Marine Regions layer/product identifier used to
#'   load `x`.
#' @param quiet Logical. Reserved for consistency with other helpers.
#'
#' @return Standardised `sf` object in EPSG:4326.
#'
#' @keywords internal
#' @noRd
.bf_standardise_marine_overlay <- function(x, region_source, marine_regions_layer, quiet = TRUE) {
  bf_require_packages("sf", context = "marine overlay helper")

  x <- bf_repair_utf8_df(x)

  if (!inherits(x, "sf")) {
    stop("Marine overlay '", region_source, "' did not load as an sf object.", call. = FALSE)
  }
  if (!nrow(x)) {
    stop("Marine overlay '", region_source, "' loaded 0 features.", call. = FALSE)
  }

  x <- bf_sf_wgs84(x)

  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  if (!nrow(x)) {
    stop("Marine overlay '", region_source, "' has only empty geometries.", call. = FALSE)
  }

  x <- bf_sf_make_valid(x, quiet = quiet)

  name_col <- .bf_marine_pick_col(
    x,
    c(
      "geoname", "GeoName", "GEONAME", "name", "Name", "NAME",
      "area_name", "AREANAME", "region_name", "marine_region_name",
      "LME_NAME", "ECOREGION", "ECO_NAME", "IHO_NAME", "FAO_NAME",
      "province", "Province", "provdescr", "POL_NAME", "POL_TYPE"
    )
  )

  id_col <- .bf_marine_pick_col(
    x,
    c(
      "marine_region_id", "mrgid", "MRGID", "mrgid_eez", "MRGID_EEZ",
      "id", "ID", "fid", "FID", "objectid", "OBJECTID", "gid", "GID",
      "LME_NUMBER", "ECO_CODE", "FAO", "POL_ID"
    )
  )

  if (is.null(name_col)) {
    x$geoname <- paste0(toupper(region_source), "_", seq_len(nrow(x)))
  } else {
    x$geoname <- .bf_marine_as_chr(x[[name_col]])
    x$geoname[is.na(x$geoname) | !nzchar(x$geoname)] <- paste0(toupper(region_source), "_", which(is.na(x$geoname) | !nzchar(x$geoname)))
  }

  if (is.null(id_col)) {
    x$marine_region_id <- paste0(.bf_marine_norm(region_source), "_", seq_len(nrow(x)))
  } else {
    x$marine_region_id <- .bf_marine_as_chr(x[[id_col]])
    x$marine_region_id[is.na(x$marine_region_id) | !nzchar(x$marine_region_id)] <- paste0(.bf_marine_norm(region_source), "_", which(is.na(x$marine_region_id) | !nzchar(x$marine_region_id)))
  }

  x$marine_region_name <- bf_repair_utf8_chr(x$geoname)
  x$marine_region_source <- .bf_marine_norm(region_source)
  x$marine_regions_layer <- bf_repair_utf8_one(marine_regions_layer)
  x <- bf_repair_utf8_df(x)

  # Preserve useful source columns but place canonical columns first.
  gcol <- attr(x, "sf_column")
  canonical <- c(
    "marine_region_id", "marine_region_name", "marine_region_source",
    "marine_regions_layer", "geoname"
  )
  keep <- unique(c(canonical, setdiff(names(x), c(canonical, gcol)), gcol))
  out <- x[, keep, drop = FALSE]
  out <- bf_repair_utf8_df(out)
  out
}

#' Load a package-managed Marine Regions overlay
#'
#' Loads one supported marine spatial overlay, caches the source product locally,
#' and returns a standardised `sf` object for use in the marine GBIF processing
#' pipeline.
#'
#' @details
#' The loader first resolves `region_source` to a canonical biofetchR overlay
#' name, then attempts to load the corresponding Marine Regions product through
#' mregions2. Candidate product identifiers are tried in order because product
#' names may differ between Marine Regions releases.
#'
#' The returned object is transformed to EPSG:4326, repaired where possible and
#' standardised to include:
#'
#' - `marine_region_id`
#' - `marine_region_name`
#' - `marine_region_source`
#' - `marine_regions_layer`
#' - `geoname`
#'
#' These columns are used by `process_gbif_marine_pipeline()` when assigning
#' occurrence records to marine spatial units and exporting grouped CSV outputs.
#'
#' @section Caching:
#' Loaded overlays are cached as RDS files under `cache_dir` using one subfolder
#' per canonical overlay. Set `force_refresh = TRUE` to ignore the cached object
#' and request a fresh copy through mregions2.
#'
#' @section Data source, licensing and attribution:
#' Marine overlays are retrieved from Marine Regions data products through
#' mregions2. biofetchR only standardises and caches those products for local
#' processing; it does not redistribute them as package data. Users should cite
#' Marine Regions and the specific product loaded, using the
#' `marine_regions_layer` output column and the `license` / `citation` fields in
#' `mregions2::mrp_list` where available.
#'
#' @param region_source Character. Marine overlay name or alias. See
#'   `bf_available_marine_overlays(include_extended = TRUE)` for registered
#'   options.
#' @param cache_dir Directory used for Marine Regions product caches. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param force_refresh Logical. If `TRUE`, ignore the local RDS cache and
#'   re-download/reload the product through mregions2.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return An `sf` polygon object in WGS84 longitude/latitude coordinates
#'   containing the requested Marine Regions overlay. The returned object
#'   includes canonical columns used by downstream marine workflows:
#'   `marine_region_id`, a character identifier for each spatial unit;
#'   `marine_region_name`, a character region label; `marine_region_source`,
#'   the canonical biofetchR overlay name; `marine_regions_layer`, the Marine
#'   Regions product identifier used to load the data; and `geoname`, a
#'   standardised human-readable region name. Source attributes may also be
#'   retained after these canonical columns. Each row represents one marine
#'   polygon or multipolygon used for assigning occurrence records to marine
#'   contextual regions. The source overlay is cached locally as an RDS file for
#'   reuse.
#'
#' @family marine overlay loaders
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   meow <- bf_load_marine_regions_overlay(
#'     region_source = "meow",
#'     cache_dir = file.path(tempdir(), "biofetchR_marine_regions")
#'   )
#'
#'   meow
#' }
#' }
#'
#' @md
#' @export
bf_load_marine_regions_overlay <- function(region_source = "eez",
                                           cache_dir = NULL,
                                           force_refresh = FALSE,
                                           quiet = TRUE) {
  bf_require_packages(c("sf", "mregions2"), context = "marine overlay helper")

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

  region_source <- .bf_marine_norm(region_source)
  .bf_marine_registry_row(region_source)

  row <- .bf_marine_registry_row(region_source)
  candidates <- unique(as.character(row$candidates[[1]]))
  resolved <- tryCatch(.bf_resolve_marine_product(region_source, quiet = quiet), error = function(e) NA_character_)
  layers_to_try <- unique(c(resolved[!is.na(resolved) & nzchar(resolved)], candidates))

  x <- NULL
  layer <- NA_character_
  errors <- character()

  for (candidate_layer in layers_to_try) {
    if (!isTRUE(quiet)) {
      message("Trying package-managed marine overlay '", toupper(region_source), "' from Marine Regions layer '", candidate_layer, "'.")
    }

    candidate_x <- tryCatch(
      .bf_mrp_get_cached(
        layer = candidate_layer,
        cache_dir = file.path(cache_dir, region_source),
        force_refresh = force_refresh,
        quiet = quiet
      ),
      error = function(e) {
        errors <<- c(errors, paste0(candidate_layer, ": ", conditionMessage(e)))
        NULL
      }
    )

    if (inherits(candidate_x, "sf") && nrow(candidate_x)) {
      x <- candidate_x
      layer <- candidate_layer
      break
    }
  }

  if (!inherits(x, "sf") || !nrow(x)) {
    stop(
      "Could not load marine overlay '", region_source, "'. Tried Marine Regions layer candidate(s): ",
      paste(layers_to_try, collapse = ", "), ". Last errors: ", paste(utils::tail(errors, 5), collapse = " | "),
      call. = FALSE
    )
  }

  x <- .bf_standardise_marine_overlay(
    x = x,
    region_source = region_source,
    marine_regions_layer = layer,
    quiet = quiet
  )

  x <- bf_repair_utf8_df(x)

  if (!isTRUE(quiet)) {
    message("Loaded marine overlay '", toupper(region_source), "' (", nrow(x), " features).")
  }

  x
}
