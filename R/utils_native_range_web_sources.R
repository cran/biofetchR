################################################################################
# utils_native_range_web_sources.R
# ------------------------------------------------------------------------------
# biofetchR: optional web/API native-origin evidence aggregator
# ------------------------------------------------------------------------------
#
# PURPOSE
#   Compile species-level native-origin evidence from optional web/API sources
#   when users do not already have a complete curated native-range table. The
#   outputs are deliberately shaped to feed directly into the core native-range
#   classifier in utils_native_range.R, especially bf_attach_native_status() and
#   bf_standardise_native_ranges().
#
# RELATIONSHIP TO OTHER NATIVE-RANGE FILES
#   This is the canonical web-source aggregator. It is not a replacement for:
#     - utils_native_range.R, which classifies recipient records as native,
#       non-native or origin-unknown; or
#     - utils_native_range_sinas.R, which contains the provider-specific SInAS
#       parser used when sources includes "sinas".
#
#   Older filenames such as utils_native_range_web_sources_with_sinas.R were
#   transitional copies. The standard utils_native_range_web_sources.R file now
#   already includes SInAS support and should be the one retained in R/.
#
# CURRENT SOURCES
#   - SInAS 3.1.1 native-location records, accessed through the package-managed
#     SInAS parser/download helpers.
#   - GBIF Species API distribution records with native-like status language.
#   - WoRMS REST Aphia distribution records with native-like status language.
#
# DESIGN PRINCIPLES
#   1. Keep web evidence optional and modular.
#   2. Use cache-first request handling to avoid repeated API calls.
#   3. Preserve unmapped source strings for audit and later manual correction.
#   4. Use conservative native-status parsing: accept native/endemic/indigenous
#      signals and reject explicitly alien/introduced/non-native signals.
#   5. Treat web-derived evidence as support for screening/auditing, not as a
#      perfect replacement for curated native-range databases.
#
# OUTPUT CONTRACT
#   bf_fetch_native_ranges_web(return = "list") returns:
#     - long: one source-evidence row per retained distribution record;
#     - species: collapsed species-level native_origin_iso3 evidence;
#     - unmapped: native-like rows whose locality could not be mapped to ISO3;
#     - summary: counts and sources used.
#
# DATA USE AND ATTRIBUTION
#   This file can query or cache third-party web/API responses. Users remain
#   responsible for citing the exact sources used in their workflow, including
#   SInAS/Zenodo records, GBIF API-derived evidence and WoRMS records where
#   applicable. biofetchR standardises and caches evidence; it does not own or
#   redistribute the underlying source databases.
#
################################################################################


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Build a filesystem-safe cache key
#'
#' @param x Text to convert into a cache-safe slug.
#'
#' @return Lowercase character slug.
#'
#' @keywords internal
#' @noRd
.bf_nrw_slug <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  if (!nzchar(x)) x <- "unknown"

  x
}


#' Detect native-like source language
#'
#' @param x Character vector containing establishment, distribution, status, or
#'   nativity text from a source.
#'
#' @return Logical vector. `TRUE` indicates a native-like signal after excluding
#'   explicitly non-native wording.
#'
#' @keywords internal
#' @noRd
.bf_nrw_is_native_like <- function(x) {
  x <- toupper(bf_clean_text(x))
  x[is.na(x)] <- ""

  native_hit <- grepl(
    "\\b(NATIVE|ENDEMIC|INDIGENOUS|AUTOCHTHONOUS|ORIGIN|NATURAL RANGE)\\b",
    x
  )

  alien_hit <- grepl(
    "\\b(NON[- ]?NATIVE|ALIEN|INTRODUCED|INVASIVE|EXOTIC|NATURALI[ZS]ED|ADVENTIVE|CULTIVATED|DOMESTICATED)\\b",
    x
  )

  native_hit & !alien_hit
}


#' Convert country, ISO2, or ISO3 evidence to ISO3
#'
#' Uses biofetchR's native-range country parser when available. Otherwise falls
#' back to `countrycode` plus a small set of common manual aliases.
#'
#' @param country Optional country/place names.
#' @param iso2c Optional ISO2 codes.
#' @param iso3c Optional ISO3 codes.
#' @param n Expected output length.
#'
#' @return Character vector of ISO3 codes, with unresolved values as `NA`.
#'
#' @keywords internal
#' @noRd
.bf_nrw_country_to_iso3 <- function(country = NULL,
                                    iso2c = NULL,
                                    iso3c = NULL,
                                    n = NULL) {
  if (exists(".bf_native_country_to_iso3", mode = "function", inherits = TRUE)) {
    return(.bf_native_country_to_iso3(
      country = country,
      iso2c = iso2c,
      iso3c = iso3c,
      n = n
    ))
  }

  if (!requireNamespace("countrycode", quietly = TRUE)) {
    if (is.null(n)) {
      n <- max(
        length(bf_null_coalesce(country, character())),
        length(bf_null_coalesce(iso2c, character())),
        length(bf_null_coalesce(iso3c, character())),
        0L
      )
    }
    return(rep(NA_character_, n))
  }

  lengths <- c(
    length(bf_null_coalesce(country, character(0))),
    length(bf_null_coalesce(iso2c, character(0))),
    length(bf_null_coalesce(iso3c, character(0)))
  )

  if (is.null(n)) {
    n <- max(lengths, 0L)
  }

  if (!length(n) || is.na(n) || n < 1L) {
    return(character(0))
  }

  recycle_to <- function(x) {
    if (is.null(x) || length(x) == 0L) return(rep(NA_character_, n))
    x <- as.character(x)
    if (length(x) == 1L && n > 1L) return(rep(x, n))
    if (length(x) != n) {
      stop("Country/ISO vector length mismatch while parsing web-native evidence.", call. = FALSE)
    }
    x
  }

  country <- recycle_to(country)
  iso2c <- recycle_to(iso2c)
  iso3c <- recycle_to(iso3c)

  out <- rep(NA_character_, n)

  x3 <- toupper(trimws(iso3c))
  x3[!grepl("^[A-Z]{3}$|^XKX$", x3)] <- NA_character_
  out[!is.na(x3)] <- x3[!is.na(x3)]

  x2 <- toupper(trimws(iso2c))
  x2[x2 == "UK"] <- "GB"
  y2 <- suppressWarnings(countrycode::countrycode(
    x2,
    origin = "iso2c",
    destination = "iso3c",
    custom_match = c("XK" = "XKX", "NA" = "NAM"),
    warn = FALSE
  ))
  fill2 <- is.na(out) & !is.na(y2)
  out[fill2] <- y2[fill2]

  # GBIF/WoRMS/native-web sources may place ISO2-like area strings in the
  # country/place field rather than in iso2c. Resolve these before attempting
  # country-name matching.
  xc_iso2 <- toupper(trimws(bf_clean_text(country)))
  xc_iso2[xc_iso2 == "UK"] <- "GB"
  xc_iso2[!grepl("^[A-Z]{2}$|^XK$", xc_iso2)] <- NA_character_

  yc_iso2 <- suppressWarnings(countrycode::countrycode(
    xc_iso2,
    origin = "iso2c",
    destination = "iso3c",
    custom_match = c("XK" = "XKX", "NA" = "NAM"),
    warn = FALSE
  ))

  fill_country_iso2 <- is.na(out) & !is.na(yc_iso2)
  out[fill_country_iso2] <- yc_iso2[fill_country_iso2]

  # Manual aliases for country, territory and narrowly bounded native-area
  # labels returned by web/API native-range sources.
  #
  # Keep this conservative. Broad regions such as "Europe", "Siberia",
  # "Middle Asia" and "Malesia" are intentionally not mapped here because they
  # would over-expand native-origin evidence across many countries.
  manual <- c(
    # Existing country and territory aliases.
    "United Kingdom" = "GBR",
    "Great Britain" = "GBR",
    "Britain" = "GBR",
    "England" = "GBR",
    "Scotland" = "GBR",
    "Wales" = "GBR",
    "Northern Ireland" = "GBR",
    "United States" = "USA",
    "United States of America" = "USA",
    "USA" = "USA",
    "Russia" = "RUS",
    "Russian Federation" = "RUS",
    "Czechia" = "CZE",
    "Czech Republic" = "CZE",
    "R\\u00e9union" = "REU",
    "Reunion" = "REU",
    "Cura\\u00e7ao" = "CUW",
    "Curacao" = "CUW",
    "Kosovo" = "XKX",

    # Conservative aliases discovered by the live native-alias harvester.
    # These are bounded enough for ISO3-level native-origin evidence.
    "British Isles" = "GBR;IRL",
    "Baltic States" = "EST;LVA;LTU",
    "Transcaucasus" = "ARM;AZE;GEO",
    "Trans Caucasus" = "ARM;AZE;GEO",
    "South Caucasus" = "ARM;AZE;GEO",
    "North Caucasus" = "RUS",
    "Gruziya" = "GEO",
    "Lebanon-Syria" = "LBN;SYR",
    "Lebanon Syria" = "LBN;SYR",
    "Sinai" = "EGY",

    # Greek island/archipelago labels.
    "East Aegean Is." = "GRC",
    "East Aegean Is" = "GRC",
    "East Aegean Islands" = "GRC",
    "Kriti" = "GRC",
    "Crete" = "GRC",

    # Spanish island/archipelago labels.
    "Baleares" = "ESP",
    "Balearic Islands" = "ESP",
    "Balearic Is." = "ESP",
    "Balearic Is" = "ESP",
    "Canary Is." = "ESP",
    "Canary Is" = "ESP",
    "Canary Islands" = "ESP",

    # Portuguese island labels.
    "Madeira" = "PRT",
    "Madeira Islands" = "PRT",

    # Italian and French island labels.
    "Sardegna" = "ITA",
    "Sardinia" = "ITA",
    "Sicilia" = "ITA",
    "Sicily" = "ITA",
    "Corse" = "FRA",
    "Corsica" = "FRA",

    # Indonesian island-region label.
    "Maluku" = "IDN",
    "Moluccas" = "IDN",

    # Marine country-region aliases discovered by test 02c.
    # These mostly come from WoRMS native-distribution strings that explicitly
    # identify a country or exclusive economic zone. These aliases only map the
    # geography string to ISO3; they do not override native/uncertain status.
    "Philippine Exclusive Economic Zone" = "PHL",
    "Spanish part of the Balearic Sea" = "ESP",
    "Belgian part of the North Sea" = "BEL",
    "Danish part of the North Sea" = "DNK",
    "Dutch part of the North Sea" = "NLD",
    "French part of Celtic Seas" = "FRA",
    "Greek part of the Adriatic Sea" = "GRC",
    "Hawaiian part of the North Pacific Ocean" = "USA",
    "Maltese part of the Mediterranean Sea - Eastern Basin" = "MLT",
    "Swedish part of the Skagerrak" = "SWE"
  )

  xc <- bf_clean_text(country)
  yc <- unname(manual[xc])
  still <- is.na(yc) & !is.na(xc)

  if (any(still)) {
    yc[still] <- suppressWarnings(countrycode::countrycode(
      xc[still],
      origin = "country.name",
      destination = "iso3c",
      custom_match = manual,
      warn = FALSE
    ))
  }

  still <- is.na(yc) & !is.na(xc)
  if (any(still)) {
    xc2 <- gsub("\\s*\\(.*\\)\\s*$", "", xc[still])
    yc[still] <- suppressWarnings(countrycode::countrycode(
      xc2,
      origin = "country.name",
      destination = "iso3c",
      custom_match = manual,
      warn = FALSE
    ))
  }

  fillc <- is.na(out) & !is.na(yc)
  out[fillc] <- yc[fillc]

  out <- toupper(out)
  out[out == ""] <- NA_character_
  out
}


#' Fetch a URL with retries and local caching
#'
#' @param url URL to fetch.
#' @param cache_file File used to cache the raw response body.
#' @param force_refresh Logical; if `TRUE`, ignore an existing cache file.
#' @param sleep_sec Delay between requests/retries.
#' @param tries Number of attempts.
#' @param user_agent HTTP user-agent string.
#' @param quiet Logical; suppress messages.
#'
#' @return Response text, or `NULL` if all attempts fail.
#'
#' @keywords internal
#' @noRd
.bf_nrw_get_cached_text <- function(url,
                                    cache_file,
                                    force_refresh = FALSE,
                                    sleep_sec = 0.25,
                                    tries = 4L,
                                    user_agent = "biofetchR native-range web helper (academic use)",
                                    quiet = TRUE) {
  bf_require_packages("httr", context = "native-range web helper")

  if (!is.null(cache_file) && file.exists(cache_file) && !isTRUE(force_refresh)) {
    txt <- tryCatch(
      paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n"),
      error = function(e) NULL
    )

    if (!is.null(txt) && nzchar(txt)) {
      return(txt)
    }

    suppressWarnings(unlink(cache_file, force = TRUE))
  }

  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)

  tries <- suppressWarnings(as.integer(tries)[1])
  if (!is.finite(tries) || tries < 1L) tries <- 1L

  for (i in seq_len(tries)) {
    resp <- tryCatch(
      httr::GET(
        url,
        httr::user_agent(user_agent),
        httr::timeout(40),
        httr::config(connecttimeout = 15)
      ),
      error = function(e) NULL
    )

    if (!is.null(resp)) {
      status <- httr::status_code(resp)

      if (status < 400L) {
        txt <- tryCatch(
          httr::content(resp, as = "text", encoding = "UTF-8"),
          error = function(e) NULL
        )

        if (!is.null(txt) && nzchar(txt)) {
          writeLines(txt, cache_file, useBytes = TRUE)
          Sys.sleep(sleep_sec)
          return(txt)
        }
      }

      if (status == 429L || status >= 500L) {
        Sys.sleep(sleep_sec * (2 ^ (i - 1L)))
      }
    }

    Sys.sleep(sleep_sec * (2 ^ (i - 1L)))
  }

  .bf_msg("Failed to fetch URL after retries: ", url, quiet = quiet)
  NULL
}


#' Fetch cached JSON from a URL
#'
#' @inheritParams .bf_nrw_get_cached_text
#'
#' @return Parsed JSON object, or `NULL` if the request or parsing fails.
#'
#' @keywords internal
#' @noRd
.bf_nrw_get_cached_json <- function(url,
                                    cache_file,
                                    force_refresh = FALSE,
                                    sleep_sec = 0.25,
                                    tries = 4L,
                                    user_agent = "biofetchR native-range web helper (academic use)",
                                    quiet = TRUE) {
  bf_require_packages("jsonlite", context = "native-range web helper")

  txt <- .bf_nrw_get_cached_text(
    url = url,
    cache_file = cache_file,
    force_refresh = force_refresh,
    sleep_sec = sleep_sec,
    tries = tries,
    user_agent = user_agent,
    quiet = quiet
  )

  if (is.null(txt) || !nzchar(txt)) {
    return(NULL)
  }

  tryCatch(
    jsonlite::fromJSON(txt, simplifyVector = FALSE),
    error = function(e) NULL
  )
}


#' Return an empty native-range web evidence table
#'
#' @return Empty tibble with the standard long-evidence columns.
#'
#' @keywords internal
#' @noRd
.bf_nrw_empty_long <- function() {
  tibble::tibble(
    species = character(),
    source = character(),
    accepted_name = character(),
    source_taxon_id = character(),
    raw_native_area = character(),
    raw_status = character(),
    origin_iso3 = character(),
    evidence_type = character(),
    source_url = character()
  )
}


#' Collapse long web evidence to species-level native-origin rows
#'
#' @param long Long evidence table returned by `bf_fetch_native_ranges_web()`.
#'
#' @return Species-level tibble compatible with `bf_standardise_native_ranges()`.
#'
#' @keywords internal
#' @noRd
.bf_nrw_collapse_species <- function(long) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range web helper")

  if (is.null(long) || !nrow(long)) {
    return(tibble::tibble(
      species = character(),
      native_origin_iso3 = character(),
      native_sources_used = character(),
      native_web_unmapped_strings = character(),
      native_web_n_records = integer(),
      native_has_origin = logical()
    ))
  }

  long <- tibble::as_tibble(long)

  split_by_species <- split(long, long$species)

  rows <- lapply(split_by_species, function(x) {
    iso <- unique(toupper(as.character(x$origin_iso3)))
    iso <- iso[!is.na(iso) & nzchar(iso)]

    src <- unique(as.character(x$source[!is.na(x$source) & nzchar(x$source)]))

    unmapped <- unique(as.character(x$raw_native_area[is.na(x$origin_iso3) | !nzchar(x$origin_iso3)]))
    unmapped <- unmapped[!is.na(unmapped) & nzchar(unmapped)]

    tibble::tibble(
      species = as.character(x$species[[1]]),
      native_origin_iso3 = if (length(iso)) paste(sort(iso), collapse = ";") else NA_character_,
      native_sources_used = if (length(src)) paste(sort(src), collapse = ";") else NA_character_,
      native_web_unmapped_strings = if (length(unmapped)) paste(sort(unmapped), collapse = ";") else NA_character_,
      native_web_n_records = nrow(x),
      native_has_origin = length(iso) > 0L
    )
  })

  dplyr::bind_rows(rows)
}


# -----------------------------------------------------------------------------
# Public helpers
# -----------------------------------------------------------------------------

#' List web sources supported for native-range evidence
#'
#' Returns the web/API sources currently implemented by
#' `bf_fetch_native_ranges_web()`. These are optional enrichment sources used to
#' compile species-level native-origin evidence when a complete local native
#' range table is not available.
#'
#' @return A data frame listing native-origin web/API evidence sources supported
#'   by biofetchR. The returned table contains `source`, the source identifier
#'   accepted by [bf_fetch_native_ranges_web()]; `description`, a short
#'   explanation of the evidence provider; `access_type`, the way the source is
#'   accessed or cached; and `default`, a logical flag indicating whether the
#'   source is included in the default source set. The table is intended for
#'   inspecting valid source names before requesting web-derived native-origin
#'   evidence.
#'
#' @examples
#' bf_available_native_web_sources()
#'
#' @section Relationship to other helpers:
#' This function lists optional evidence providers only. The returned sources
#' are consumed by [bf_fetch_native_ranges_web()] and then passed to the core
#' native-status classifier in [bf_attach_native_status()].
#'
#' @family native-range web-source helpers
#' @md
#' @export
bf_available_native_web_sources <- function() {
  data.frame(
    source = c("sinas", "gbif", "worms"),
    description = c(
      "SInAS 3.1.1 native-location records resolved through AllLocations",
      "GBIF Species API distributions with native-like establishment/status fields",
      "WoRMS REST Aphia distributions for marine and aquatic taxa"
    ),
    access_type = c("zenodo_cache", "json_api", "json_api"),
    default = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}


#' Fetch native-range evidence from the GBIF Species API
#'
#' Matches species names to GBIF taxon keys, retrieves checklist distribution
#' records, keeps distribution rows with native-like establishment/status
#' wording, and attempts to convert reported areas or country codes to ISO3.
#'
#' @param species Character vector of species names.
#' @param cache_dir Directory used to cache GBIF JSON responses. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, ignore cached responses and fetch
#'   from the API again.
#' @param sleep_sec Delay between requests, in seconds.
#' @param user_agent HTTP user-agent string.
#' @param quiet Logical; suppress progress messages.
#'
#' @return A tibble with one row per retained native-like GBIF Species API
#'   distribution record. The output uses the standard long native-evidence
#'   schema: `species`, the requested species name; `source`, set to `"GBIF"`;
#'   `accepted_name`, the matched GBIF scientific or canonical name;
#'   `source_taxon_id`, the GBIF taxon key; `raw_native_area`, the source
#'   country, area, locality or location string; `raw_status`, the establishment
#'   or status text used for native-like screening; `origin_iso3`, the resolved
#'   ISO3 country code where available; `evidence_type`, usually
#'   `"native_distribution"`; and `source_url`, the GBIF API endpoint used. If
#'   no native-like records are retained, an empty tibble with the same columns
#'   is returned.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   gbif_native <- bf_web_native_gbif(
#'     species = "Carcinus maenas",
#'     cache_dir = file.path(tempdir(), "biofetchR_native_web"),
#'     quiet = FALSE
#'   )
#'
#'   gbif_native
#' }
#' }
#'
#' @section Data source and interpretation:
#' This helper uses GBIF Species API checklist-distribution records, not GBIF
#' occurrence downloads. Retained rows are evidence of native-like distribution
#' language in the API response and should be audited before being treated as a
#' definitive native range.
#'
#' @family native-range web-source helpers
#' @md
#' @export
bf_web_native_gbif <- function(species,
                               cache_dir = NULL,
                               force_refresh = FALSE,
                               sleep_sec = 0.25,
                               user_agent = "biofetchR native-range web helper (academic use)",
                               quiet = FALSE) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range web helper")

  species <- unique(bf_clean_text(species))
  species <- species[!is.na(species) & nzchar(species)]

  if (!length(species)) {
    return(.bf_nrw_empty_long())
  }

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

  cache_dir <- file.path(cache_dir, "native_range_web", "gbif")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  out <- list()
  out_i <- 0L

  for (sp in species) {
    .bf_msg("GBIF native-range query: ", sp, quiet = quiet)

    sp_slug <- .bf_nrw_slug(sp)
    match_url <- paste0(
      "https://api.gbif.org/v1/species/match?name=",
      utils::URLencode(sp, reserved = TRUE)
    )

    match_json <- .bf_nrw_get_cached_json(
      url = match_url,
      cache_file = file.path(cache_dir, paste0(sp_slug, "_match.json")),
      force_refresh = force_refresh,
      sleep_sec = sleep_sec,
      user_agent = user_agent,
      quiet = quiet
    )

    taxon_key <- suppressWarnings(as.integer(bf_clean_text_one(bf_null_coalesce(match_json$usageKey, match_json$speciesKey))))
    accepted <- bf_clean_text_one(bf_null_coalesce(match_json$scientificName, match_json$canonicalName), fallback = sp)

    if (is.na(taxon_key)) {
      next
    }

    dist_url <- paste0(
      "https://api.gbif.org/v1/species/",
      taxon_key,
      "/distributions?limit=1000"
    )

    dist_json <- .bf_nrw_get_cached_json(
      url = dist_url,
      cache_file = file.path(cache_dir, paste0("key_", taxon_key, "_distributions.json")),
      force_refresh = force_refresh,
      sleep_sec = sleep_sec,
      user_agent = user_agent,
      quiet = quiet
    )

    records <- dist_json$results
    if (is.null(records) || !length(records)) next

    for (rec in records) {
      status_txt <- paste(
        bf_clean_text_one(rec$establishmentMeans, fallback = ""),
        bf_clean_text_one(rec$status, fallback = ""),
        bf_clean_text_one(rec$occurrenceStatus, fallback = ""),
        sep = " "
      )

      if (!.bf_nrw_is_native_like(status_txt)) {
        next
      }

      country_code <- bf_clean_text_one(rec$countryCode)
      area <- bf_clean_text_one(bf_null_coalesce(rec$country, bf_null_coalesce(rec$area, bf_null_coalesce(rec$locality, rec$locationID))))

      iso3 <- NA_character_

      if (!is.na(country_code) && nzchar(country_code)) {
        if (grepl("^[A-Za-z]{2}$", country_code)) {
          iso3 <- .bf_nrw_country_to_iso3(iso2c = country_code, n = 1L)
        } else if (grepl("^[A-Za-z]{3}$|^XKX$", country_code)) {
          iso3 <- toupper(country_code)
        }
      }

      if ((is.na(iso3) || !nzchar(iso3)) && !is.na(area) && nzchar(area)) {
        iso3 <- .bf_nrw_country_to_iso3(country = area, n = 1L)
      }

      out_i <- out_i + 1L
      out[[out_i]] <- tibble::tibble(
        species = sp,
        source = "GBIF",
        accepted_name = accepted,
        source_taxon_id = as.character(taxon_key),
        raw_native_area = area,
        raw_status = bf_clean_text(status_txt),
        origin_iso3 = iso3,
        evidence_type = "native_distribution",
        source_url = dist_url
      )
    }
  }

  if (!length(out)) {
    return(.bf_nrw_empty_long())
  }

  dplyr::bind_rows(out)
}


#' Fetch native-range evidence from WoRMS REST distributions
#'
#' Resolves species names to valid AphiaIDs through the WoRMS REST API, retrieves
#' Aphia distribution records, keeps native-like distribution/status rows, and
#' attempts to map reported localities to ISO3 countries.
#'
#' @param species Character vector of species names.
#' @param cache_dir Directory used to cache WoRMS JSON responses. Must be
#'   supplied explicitly. In examples, tests and vignettes, use a path under
#'   `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, ignore cached responses and fetch
#'   from the API again.
#' @param sleep_sec Delay between requests, in seconds.
#' @param marine_only Logical passed to the WoRMS name-resolution endpoint.
#' @param user_agent HTTP user-agent string.
#' @param quiet Logical; suppress progress messages.
#'
#' @return A tibble with one row per retained native-like WoRMS Aphia
#'   distribution record. The output uses the standard long native-evidence
#'   schema: `species`, the requested species name; `source`, set to `"WoRMS"`;
#'   `accepted_name`, the valid or accepted WoRMS name; `source_taxon_id`, the
#'   AphiaID; `raw_native_area`, the source locality, area or geounit string;
#'   `raw_status`, the status, origin or establishment text used for native-like
#'   screening; `origin_iso3`, the resolved ISO3 country code where available;
#'   `evidence_type`, usually `"native_distribution"`; and `source_url`, the
#'   WoRMS REST endpoint used. If no native-like records are retained, an empty
#'   tibble with the same columns is returned.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   worms_native <- bf_web_native_worms(
#'     species = "Carcinus maenas",
#'     cache_dir = file.path(tempdir(), "biofetchR_native_web"),
#'     quiet = FALSE
#'   )
#'
#'   worms_native
#' }
#' }
#'
#' @section Data source and interpretation:
#' This helper uses WoRMS REST Aphia name and distribution endpoints. Retained
#' rows are native-like distribution evidence and may remain unmapped when the
#' source locality cannot be resolved to an ISO3 country code.
#'
#' @family native-range web-source helpers
#' @md
#' @export
bf_web_native_worms <- function(species,
                                cache_dir = NULL,
                                force_refresh = FALSE,
                                sleep_sec = 0.25,
                                marine_only = FALSE,
                                user_agent = "biofetchR native-range web helper (academic use)",
                                quiet = FALSE) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range web helper")

  species <- unique(bf_clean_text(species))
  species <- species[!is.na(species) & nzchar(species)]

  if (!length(species)) {
    return(.bf_nrw_empty_long())
  }

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

  cache_dir <- file.path(cache_dir, "native_range_web", "worms")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  out <- list()
  out_i <- 0L

  for (sp in species) {
    .bf_msg("WoRMS native-range query: ", sp, quiet = quiet)

    sp_slug <- .bf_nrw_slug(sp)
    name_url <- paste0(
      "https://www.marinespecies.org/rest/AphiaRecordsByName/",
      utils::URLencode(sp, reserved = TRUE),
      "?like=false&marine_only=",
      tolower(as.character(isTRUE(marine_only)))
    )

    name_json <- .bf_nrw_get_cached_json(
      url = name_url,
      cache_file = file.path(cache_dir, paste0(sp_slug, "_name.json")),
      force_refresh = force_refresh,
      sleep_sec = sleep_sec,
      user_agent = user_agent,
      quiet = quiet
    )

    if (is.null(name_json) || !length(name_json)) next

    recs <- name_json
    if (is.data.frame(recs)) recs <- split(recs, seq_len(nrow(recs)))

    score <- vapply(recs, function(z) {
      status <- tolower(bf_clean_text_one(z$status, fallback = ""))
      valid <- tolower(bf_clean_text_one(z$valid_name, fallback = ""))
      aphia <- suppressWarnings(as.integer(bf_clean_text_one(z$AphiaID)))
      as.integer(!is.na(aphia)) +
        as.integer(grepl("accepted|valid", status)) * 2L +
        as.integer(nzchar(valid))
    }, integer(1))

    rec <- recs[[which.max(score)]]
    aphia_id <- suppressWarnings(as.integer(bf_clean_text_one(rec$AphiaID)))
    accepted <- bf_clean_text_one(bf_null_coalesce(rec$valid_name, bf_null_coalesce(rec$scientificname, rec$scientificName)), fallback = sp)

    if (is.na(aphia_id)) next

    dist_url <- paste0(
      "https://www.marinespecies.org/rest/AphiaDistributionsByAphiaID/",
      aphia_id
    )

    dist_json <- .bf_nrw_get_cached_json(
      url = dist_url,
      cache_file = file.path(cache_dir, paste0("aphia_", aphia_id, "_distributions.json")),
      force_refresh = force_refresh,
      sleep_sec = sleep_sec,
      user_agent = user_agent,
      quiet = quiet
    )

    if (is.null(dist_json) || !length(dist_json)) next

    dists <- dist_json
    if (is.data.frame(dists)) dists <- split(dists, seq_len(nrow(dists)))

    for (d in dists) {
      locality <- bf_clean_text_one(bf_null_coalesce(d$locality, bf_null_coalesce(d$location, bf_null_coalesce(d$area, d$geounit))))
      status_txt <- paste(
        bf_clean_text_one(d$status, fallback = ""),
        bf_clean_text_one(d$origin, fallback = ""),
        bf_clean_text_one(d$establishmentMeans, fallback = ""),
        sep = " "
      )

      if (!.bf_nrw_is_native_like(status_txt)) {
        next
      }

      iso3 <- if (!is.na(locality) && nzchar(locality)) {
        .bf_nrw_country_to_iso3(country = locality, n = 1L)
      } else {
        NA_character_
      }

      out_i <- out_i + 1L
      out[[out_i]] <- tibble::tibble(
        species = sp,
        source = "WoRMS",
        accepted_name = accepted,
        source_taxon_id = as.character(aphia_id),
        raw_native_area = locality,
        raw_status = bf_clean_text(status_txt),
        origin_iso3 = iso3,
        evidence_type = "native_distribution",
        source_url = dist_url
      )
    }
  }

  if (!length(out)) {
    return(.bf_nrw_empty_long())
  }

  dplyr::bind_rows(out)
}


#' Compile species-level native-origin evidence from web sources
#'
#' Queries one or more web/API sources for species-level native-origin evidence,
#' caches raw responses, converts source locality/country strings to ISO3 where
#' possible, and returns both long evidence and species-level collapsed outputs.
#'
#' This helper is designed as an optional complement to curated native-range
#' tables. Its collapsed `species` output can be passed to
#' `bf_attach_native_status()` as `native_ranges`, because it includes a
#' `species` column and a semicolon-delimited `native_origin_iso3` column.
#'
#' @param species Character vector of species names, or a data frame containing
#'   a species column.
#' @param species_col Species column name when `species` is a data frame.
#' @param sources Character vector of native evidence sources. Currently supports
#'   `"sinas"`, `"gbif"` and `"worms"`.
#' @param cache_dir Directory used for all web-response caches. Must be supplied
#'   explicitly when web/API sources need to be queried or cached. In examples,
#'   tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; if `TRUE`, ignore cached responses and fetch
#'   from web sources again.
#' @param sleep_sec Delay between requests, in seconds.
#' @param quiet Logical; suppress progress messages.
#' @param sinas_main_path Optional local path to `SInAS_3.1.1.csv`.
#' @param sinas_alllocations_path Optional local path to `AllLocations.xlsx`,
#'   `.csv` or `.tsv`.
#' @param sinas_fulltaxa_path Optional local path to
#'   `SInAS_3.1.1_FullTaxaList.csv`.
#' @param sinas_record_id Zenodo record identifier used for package-managed
#'   SInAS downloads.
#' @param sinas_main_url Optional direct URL for the main SInAS CSV.
#' @param sinas_fulltaxa_url Optional direct URL for the FullTaxaList CSV.
#' @param sinas_config_zip_url Optional direct URL for the SInAS config archive.
#' @param return One of `"list"`, `"species"` or `"long"`.
#'
#' @return The returned object depends on `return`. If `return = "list"`, the
#'   function returns a named list with four elements: `long`, a tibble with one
#'   row per retained source evidence record using the standard columns
#'   `species`, `source`, `accepted_name`, `source_taxon_id`,
#'   `raw_native_area`, `raw_status`, `origin_iso3`, `evidence_type` and
#'   `source_url`; `species`, a collapsed species-level tibble containing
#'   `species`, `native_origin_iso3`, `native_sources_used`,
#'   `native_web_unmapped_strings`, `native_web_n_records` and
#'   `native_has_origin`; `unmapped`, the subset of long-format native-like
#'   evidence rows whose locality or area string could not be resolved to ISO3;
#'   and `summary`, a one-row tibble reporting input species counts, matched web
#'   records, ISO3-resolved species, unmapped records and sources used. If
#'   `return = "species"`, only the collapsed species-level tibble is returned.
#'   If `return = "long"`, only the long evidence tibble is returned. These
#'   outputs provide species-level native-origin evidence for downstream use in
#'   [bf_attach_native_status()] and do not themselves classify recipient records
#'   as native or non-native.
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   native_web <- bf_fetch_native_ranges_web(
#'     species = c("Carcinus maenas", "Ficopomatus enigmaticus"),
#'     sources = c("gbif", "worms"),
#'     cache_dir = file.path(tempdir(), "biofetchR_native_web"),
#'     return = "species",
#'     quiet = FALSE
#'   )
#'
#'   native_web
#' }
#' }
#'
#' @section Relationship to SInAS:
#' SInAS support is now part of this canonical web-source aggregator. When
#' `sources` includes `"sinas"`, this function delegates provider-specific
#' parsing to [bf_fetch_native_ranges_sinas()]. Keep that SInAS provider file in
#' the package, but do not also keep a duplicate `*_with_sinas` aggregator file.
#'
#' @section Data source and interpretation:
#' The returned evidence is web/API-derived and source-dependent. Treat it as an
#' auditable native-origin evidence layer rather than proof that a species is
#' native, non-native or absent from any place.
#'
#' @family native-range web-source helpers
#' @md
#' @export
bf_fetch_native_ranges_web <- function(species,
                                       species_col = "species",
                                       sources = c("sinas", "gbif", "worms"),
                                       cache_dir = NULL,
                                       force_refresh = FALSE,
                                       sleep_sec = 0.25,
                                       quiet = FALSE,
                                       sinas_main_path = NULL,
                                       sinas_alllocations_path = NULL,
                                       sinas_fulltaxa_path = NULL,
                                       sinas_record_id = "18220953",
                                       sinas_main_url = NULL,
                                       sinas_fulltaxa_url = NULL,
                                       sinas_config_zip_url = NULL,
                                       return = c("list", "species", "long")) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range web helper")

  return <- match.arg(return)

  if (is.data.frame(species)) {
    if (!species_col %in% names(species)) {
      stop("`species_col` was not found in the supplied data frame: ", species_col, call. = FALSE)
    }

    species_vec <- species[[species_col]]
  } else {
    species_vec <- species
  }

  species_vec <- unique(bf_clean_text(species_vec))
  species_vec <- species_vec[!is.na(species_vec) & nzchar(species_vec)]

  if (!length(species_vec)) {
    long <- .bf_nrw_empty_long()
    species_out <- .bf_nrw_collapse_species(long)
  } else {
    sources <- unique(tolower(trimws(as.character(sources))))
    supported <- bf_available_native_web_sources()$source
    bad <- setdiff(sources, supported)

    if (length(bad)) {
      stop(
        "Unsupported native web source(s): ",
        paste(bad, collapse = ", "),
        ". Supported sources are: ",
        paste(supported, collapse = ", "),
        call. = FALSE
      )
    }

    if (is.null(cache_dir) || length(cache_dir) == 0L ||
        !nzchar(trimws(as.character(cache_dir[[1L]])))) {
      stop(
        "`cache_dir` must be supplied explicitly when native web/API evidence is requested. ",
        "In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
        call. = FALSE
      )
    }

    cache_dir <- normalizePath(
      as.character(cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )

    parts <- list()

    if ("sinas" %in% sources) {
      if (!exists("bf_fetch_native_ranges_sinas", mode = "function", inherits = TRUE)) {
        stop(
          "`sources` includes 'sinas', but `bf_fetch_native_ranges_sinas()` is not available. ",
          "Add `R/utils_native_range_sinas.R` and `R/utils_sinas_download.R` to the package.",
          call. = FALSE
        )
      }

      parts$sinas <- bf_fetch_native_ranges_sinas(
        species = species_vec,
        cache_dir = cache_dir,
        force_refresh = force_refresh,
        quiet = quiet,
        main_path = sinas_main_path,
        alllocations_path = sinas_alllocations_path,
        fulltaxa_path = sinas_fulltaxa_path,
        record_id = sinas_record_id,
        main_url = sinas_main_url,
        fulltaxa_url = sinas_fulltaxa_url,
        config_zip_url = sinas_config_zip_url,
        return = "long"
      )
    }

    if ("gbif" %in% sources) {
      parts$gbif <- bf_web_native_gbif(
        species = species_vec,
        cache_dir = cache_dir,
        force_refresh = force_refresh,
        sleep_sec = sleep_sec,
        quiet = quiet
      )
    }

    if ("worms" %in% sources) {
      parts$worms <- bf_web_native_worms(
        species = species_vec,
        cache_dir = cache_dir,
        force_refresh = force_refresh,
        sleep_sec = sleep_sec,
        quiet = quiet
      )
    }

    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0L, parts)

    long <- if (length(parts)) {
      dplyr::bind_rows(parts)
    } else {
      .bf_nrw_empty_long()
    }

    species_out <- .bf_nrw_collapse_species(long)
  }

  unmapped <- long[is.na(long$origin_iso3) | !nzchar(long$origin_iso3), , drop = FALSE]

  sources_used <- paste(sort(unique(long$source)), collapse = ";")

  summary <- tibble::tibble(
    n_species_input = length(species_vec),
    n_species_with_web_records = length(unique(long$species)),
    n_species_with_iso3 = sum(species_out$native_has_origin %in% TRUE),
    n_long_records = nrow(long),
    n_unmapped_records = nrow(unmapped),
    sources = sources_used,
    sources_used = sources_used
  )

  out <- list(
    long = tibble::as_tibble(long),
    species = tibble::as_tibble(species_out),
    unmapped = tibble::as_tibble(unmapped),
    summary = summary
  )

  if (return == "species") {
    return(out$species)
  }

  if (return == "long") {
    return(out$long)
  }

  out
}


#' Write native-range web evidence outputs
#'
#' Writes the list returned by `bf_fetch_native_ranges_web()` to CSV files so
#' the web-derived native-origin evidence can be audited and reused.
#'
#' @param x Result from `bf_fetch_native_ranges_web(return = "list")`.
#' @param output_dir Directory for CSV outputs.
#' @param prefix Filename prefix. Defaults to `"native_web"`.
#'
#' @return Invisibly returns a named character vector of file paths written to
#'   `output_dir`. The vector contains paths named `long`, `species`,
#'   `unmapped` and `summary`, corresponding to the CSV files written from the
#'   matching elements of `x`. The function is called primarily for its side
#'   effect of writing auditable native-web evidence tables to disk.
#'
#' @examples
#' native_web <- list(
#'   long = data.frame(
#'     species = "Example species",
#'     source = "example",
#'     accepted_name = "Example species",
#'     source_taxon_id = "example_id",
#'     raw_native_area = "Exampleland",
#'     raw_status = "native",
#'     origin_iso3 = "GBR",
#'     evidence_type = "native_distribution",
#'     source_url = NA_character_
#'   ),
#'   species = data.frame(
#'     species = "Example species",
#'     native_origin_iso3 = "GBR",
#'     native_sources_used = "example",
#'     native_web_unmapped_strings = NA_character_,
#'     native_web_n_records = 1L,
#'     native_has_origin = TRUE
#'   ),
#'   unmapped = data.frame(),
#'   summary = data.frame(
#'     n_species_input = 1L,
#'     n_species_with_origin = 1L,
#'     n_records_long = 1L,
#'     n_records_unmapped = 0L
#'   )
#' )
#'
#' out_dir <- tempfile("native_web_outputs_")
#' paths <- bf_write_native_web_outputs(native_web, output_dir = out_dir)
#' names(paths)
#'
#' @section Audit role:
#' Writing these outputs preserves the exact web-derived evidence used by a
#' pipeline run, including unmapped source strings. This is useful for later
#' manual checking and for reporting which sources contributed native-origin
#' evidence.
#'
#' @family native-range web-source helpers
#' @md
#' @export
bf_write_native_web_outputs <- function(x,
                                        output_dir,
                                        prefix = "native_web") {
  bf_require_packages("readr", context = "native-range web helper")

  if (!is.list(x) || !all(c("long", "species", "unmapped", "summary") %in% names(x))) {
    stop("`x` must be the list output from bf_fetch_native_ranges_web(return = 'list').", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- c(
    long = file.path(output_dir, paste0(prefix, "_long.csv")),
    species = file.path(output_dir, paste0(prefix, "_species.csv")),
    unmapped = file.path(output_dir, paste0(prefix, "_unmapped.csv")),
    summary = file.path(output_dir, paste0(prefix, "_summary.csv"))
  )

  readr::write_csv(x$long, paths[["long"]])
  readr::write_csv(x$species, paths[["species"]])
  readr::write_csv(x$unmapped, paths[["unmapped"]])
  readr::write_csv(x$summary, paths[["summary"]])

  invisible(paths)
}
