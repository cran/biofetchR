################################################################################
# utils_native_range.R
# ------------------------------------------------------------------------------
# biofetchR: core native-origin status classification utilities
# ------------------------------------------------------------------------------
#
# PURPOSE
#   Core native-range classifier used by the terrestrial/freshwater and marine
#   pipelines. This file accepts species-level native-origin evidence from
#   provider-specific helpers, such as SInAS, GBIF and WoRMS, or from user-supplied
#   native-range tables. It then classifies whether a recipient country is inside
#   or outside each species' native-origin set.
#
# RELATIONSHIP TO utils_native_range_sinas.R
#   This file is not a duplicate of utils_native_range_sinas.R. The SInAS script
#   is a provider parser/evidence source. This script is the shared classifier
#   and filter layer. It standardises native-origin evidence, attaches recipient
#   status, filters confirmed native/non-native rows, summarises status columns
#   and reconciles native-origin evidence with GRIIS evidence.
#
# MAIN PUBLIC HELPERS
#   bf_standardise_native_ranges()
#   bf_native_range_lookup()
#   bf_attach_native_status()
#   bf_filter_non_native()
#   bf_filter_native()
#   bf_native_status_summary()
#   bf_reconcile_griis_native_status()
#
# INTERNAL CONTRACT
#   Provider evidence tables can be passed directly via native_ranges when they
#   contain the following species-level fields:
#     species
#     native_origin_iso3
#     native_sources_used
#     native_has_origin
#
# DATA USE AND ATTRIBUTION
#   This classifier does not own provider datasets. When native_ranges are derived
#   from SInAS, GBIF, WoRMS or another source, users should cite the exact data
#   products/releases used to build that evidence and follow the corresponding
#   licence or attribution requirements.
#
################################################################################


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Create a conservative species-name key
#'
#' @param x Character vector of species names.
#'
#' @return Lowercase name key.
#'
#' @keywords internal
#' @noRd
.bf_native_name_key <- function(x) {
  x <- bf_clean_text(x)

  if (exists("bf_clean_taxon_names", mode = "function", inherits = TRUE)) {
    cleaned <- tryCatch(
      bf_clean_taxon_names(x),
      error = function(e) x
    )
    cleaned[is.na(cleaned)] <- x[is.na(cleaned)]
    x <- cleaned
  } else {
    x <- gsub("x", "x", x, fixed = TRUE)
    x <- gsub("\\s*\\([^)]*\\)", " ", x)
    x <- gsub("\\b(Linnaeus|Pallas|Fauvel|Darwin|Thunberg|Molina|Michx\\.|L\\.|Gmelin)\\b.*$", " ", x)
    x <- gsub("[[:space:]]+", " ", x)
    x <- trimws(x)
  }

  x <- tolower(bf_clean_text(x))
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(gsub("[[:space:]]+", " ", x))
  x[x == ""] <- NA_character_
  x
}


#' Split semicolon/comma/pipe native ISO3 strings
#'
#' @param x Character vector containing ISO3 values.
#' @param max_origins Maximum number of origin codes to keep per element.
#'
#' @return List of character vectors.
#'
#' @keywords internal
#' @noRd
.bf_native_split_iso3 <- function(x, max_origins = 100L) {
  x <- bf_clean_text(x)

  lapply(x, function(z) {
    if (is.na(z) || !nzchar(z)) return(character(0))

    z <- gsub("[,|/]+", ";", z)
    bits <- unlist(strsplit(z, ";", fixed = TRUE), use.names = FALSE)
    bits <- toupper(trimws(bits))
    bits <- bits[grepl("^[A-Z]{3}$|^XKX$", bits)]
    bits <- unique(bits)

    if (length(bits) > max_origins) {
      bits <- bits[seq_len(max_origins)]
    }

    bits
  })
}


#' Collapse ISO3-code vectors to semicolon-delimited strings
#'
#' @param x List of character vectors.
#'
#' @return Character vector.
#'
#' @keywords internal
#' @noRd
.bf_native_collapse_iso3 <- function(x) {
  vapply(x, function(z) {
    z <- unique(toupper(trimws(as.character(z))))
    z <- z[grepl("^[A-Z]{3}$|^XKX$", z)]
    if (!length(z)) return(NA_character_)
    paste(sort(z), collapse = ";")
  }, character(1))
}


#' Convert country, ISO2 or ISO3 evidence to ISO3
#'
#' This helper is also used by `utils_native_range_web_sources.R` when available.
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
.bf_native_country_to_iso3 <- function(country = NULL,
                                       iso2c = NULL,
                                       iso3c = NULL,
                                       n = NULL) {
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

  if (is.null(n)) n <- max(lengths, 0L)
  if (!length(n) || is.na(n) || n < 1L) return(character(0))

  recycle_to <- function(x) {
    if (is.null(x) || length(x) == 0L) return(rep(NA_character_, n))
    x <- as.character(x)
    if (length(x) == 1L && n > 1L) return(rep(x, n))
    if (length(x) != n) {
      stop("Country/ISO vector length mismatch while parsing native-origin evidence.", call. = FALSE)
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

  # Some web/API sources, especially GBIF distribution records, can report
  # native areas as short ISO2-like strings in the country/place field rather
  # than in a dedicated iso2c column. Treat two-letter country strings such as
  # BR, VE and ZA as ISO2 codes before attempting country-name matching.
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

    "Sardegna" = "ITA",
    "Sardinia" = "ITA",
    "Sicilia" = "ITA",
    "Sicily" = "ITA",
    "Corse" = "FRA",
    "Corsica" = "FRA",
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
    "Swedish part of the Skagerrak" = "SWE",

    # Indonesian island-region label.
    "Maluku" = "IDN",
    "Moluccas" = "IDN"
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


#' Pick native-origin columns from a native-range table
#'
#' @param x Native-range data frame.
#' @param origin_iso3_col Optional explicit origin column.
#'
#' @return Character vector of origin columns.
#'
#' @keywords internal
#' @noRd
.bf_native_origin_cols <- function(x, origin_iso3_col = NULL) {
  if (!is.null(origin_iso3_col) && origin_iso3_col %in% names(x)) {
    return(origin_iso3_col)
  }

  priority <- c(
    "native_origin_iso3",
    "origin_iso3",
    "origin_iso3_sinas",
    "origin_iso3_gbif",
    "origin_iso3_cabi",
    "origin_iso3_griis",
    "native_iso3",
    "native_range_iso3"
  )

  hit <- priority[priority %in% names(x)]

  if (length(hit)) return(hit)

  grep("origin.*iso3|native.*iso3", names(x), ignore.case = TRUE, value = TRUE)
}


# -----------------------------------------------------------------------------
# Public helpers
# -----------------------------------------------------------------------------

#' Standardise native-range evidence to a species-level lookup
#'
#' Converts one or more native-origin ISO3 columns into a compact, one-row-per-
#' species lookup table. This is the normal entry point for provider outputs from
#' SInAS/GBIF/WoRMS helpers or for user-supplied native-range tables.
#'
#' @details
#' The output is species-level evidence only. Recipient-level native/non-native
#' classification is performed by [bf_attach_native_status()], which compares the
#' recipient ISO3 code against the native-origin ISO3 set.
#'
#' @param native_raw Data frame containing native-origin evidence.
#' @param species_col Species column in `native_raw`.
#' @param origin_iso3_col Optional explicit ISO3-origin column.
#' @param strategy Origin-column strategy. Currently `"tiered"` and `"union"`
#'   both collapse all available ISO3 origin columns conservatively.
#' @param max_origins Maximum number of ISO3 origins to retain per species.
#' @param quiet Logical; suppress messages.
#'
#' @return A tibble with one row per standardised species-name key. The output
#'   contains `native_name`, the cleaned species label; `native_name_key`, the
#'   normalised name used for joins; `native_origin_iso3`, a semicolon-delimited
#'   character string of ISO3 native-origin country codes; `native_origin_flag`,
#'   a character flag indicating whether usable origin evidence was found;
#'   `native_has_origin`, a logical indicator of origin-evidence availability;
#'   and `native_sources_used`, a semicolon-delimited provenance field where
#'   source information is available. The table represents species-level
#'   native-origin evidence only; recipient-level native/non-native
#'   classification is performed by [bf_attach_native_status()].
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_standardise_native_ranges <- function(native_raw,
                                         species_col = "species",
                                         origin_iso3_col = NULL,
                                         strategy = c("tiered", "union"),
                                         max_origins = 100L,
                                         quiet = FALSE) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range helper")

  strategy <- match.arg(strategy)

  if (is.null(native_raw) || !is.data.frame(native_raw)) {
    stop("`native_raw` must be a data frame.", call. = FALSE)
  }

  if (!species_col %in% names(native_raw)) {
    stop("`species_col` was not found in `native_raw`: ", species_col, call. = FALSE)
  }

  x <- tibble::as_tibble(native_raw)

  origin_cols <- .bf_native_origin_cols(x, origin_iso3_col = origin_iso3_col)

  if (!length(origin_cols)) {
    stop(
      "Could not find a native-origin ISO3 column in `native_raw`. ",
      "Expected one of: native_origin_iso3, origin_iso3, origin_iso3_sinas, ",
      "origin_iso3_gbif, origin_iso3_cabi, native_iso3.",
      call. = FALSE
    )
  }

  name <- bf_clean_text(x[[species_col]])
  key <- .bf_native_name_key(name)

  rows <- lapply(seq_len(nrow(x)), function(i) {
    iso_list <- list()

    for (col in origin_cols) {
      iso_list[[col]] <- .bf_native_split_iso3(x[[col]][[i]], max_origins = max_origins)[[1]]
    }

    iso <- unique(unlist(iso_list, use.names = FALSE))
    iso <- iso[grepl("^[A-Z]{3}$|^XKX$", iso)]

    sources <- character(0)

    if ("native_sources_used" %in% names(x)) {
      src <- bf_clean_text(x$native_sources_used[[i]])
      if (!is.na(src)) sources <- unlist(strsplit(src, ";", fixed = TRUE), use.names = FALSE)
    }

    # Add origin-column provenance where source strings are not supplied.
    if (!length(sources)) {
      non_empty_cols <- origin_cols[vapply(iso_list, length, integer(1)) > 0L]
      sources <- sub("^origin_iso3_", "", non_empty_cols)
      sources <- sub("^native_origin_iso3$", "native", sources)
      sources <- toupper(sources)
    }

    native_has_origin <- length(iso) > 0L

    tibble::tibble(
      native_name = name[[i]],
      native_name_key = key[[i]],
      native_origin_iso3 = if (native_has_origin) paste(sort(iso), collapse = ";") else NA_character_,
      native_origin_flag = if (native_has_origin) "HAS_ORIGIN" else "NO_ORIGIN",
      native_has_origin = native_has_origin,
      native_sources_used = if (length(sources)) paste(sort(unique(bf_clean_text(sources))), collapse = ";") else NA_character_
    )
  })

  out <- dplyr::bind_rows(rows)
  out <- out[!is.na(out$native_name_key) & nzchar(out$native_name_key), , drop = FALSE]

  # Collapse duplicate species rows by unioning their ISO3 origins.
  split_rows <- split(out, out$native_name_key)

  out2 <- lapply(split_rows, function(z) {
    iso <- unique(unlist(.bf_native_split_iso3(z$native_origin_iso3, max_origins = max_origins), use.names = FALSE))
    iso <- iso[grepl("^[A-Z]{3}$|^XKX$", iso)]

    src <- unique(unlist(strsplit(
      paste(z$native_sources_used[!is.na(z$native_sources_used)], collapse = ";"),
      ";",
      fixed = TRUE
    ), use.names = FALSE))
    src <- bf_clean_text(src)
    src <- src[!is.na(src) & nzchar(src)]

    tibble::tibble(
      native_name = z$native_name[[which(!is.na(z$native_name))[1] %||% 1L]],
      native_name_key = z$native_name_key[[1]],
      native_origin_iso3 = if (length(iso)) paste(sort(iso), collapse = ";") else NA_character_,
      native_origin_flag = if (length(iso)) "HAS_ORIGIN" else "NO_ORIGIN",
      native_has_origin = length(iso) > 0L,
      native_sources_used = if (length(src)) paste(sort(unique(src)), collapse = ";") else NA_character_
    )
  })

  dplyr::bind_rows(out2)
}


#' Build a native-range lookup table
#'
#' Thin wrapper around [bf_standardise_native_ranges()] retained for readability
#' in pipelines and tests.
#'
#' @inheritParams bf_standardise_native_ranges
#' @param native_ranges Data frame containing native-origin evidence to standardise.
#'
#' @return A tibble with the same structure as [bf_standardise_native_ranges()]:
#'   one row per standardised species-name key and columns describing
#'   species-level native-origin evidence, including `native_name`,
#'   `native_name_key`, `native_origin_iso3`, `native_origin_flag`,
#'   `native_has_origin` and `native_sources_used`. This function is a wrapper
#'   around [bf_standardise_native_ranges()] and returns the lookup table used by
#'   downstream native-status joins.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_native_range_lookup <- function(native_ranges,
                                   species_col = "species",
                                   origin_iso3_col = NULL,
                                   strategy = c("tiered", "union"),
                                   max_origins = 100L,
                                   quiet = FALSE) {
  bf_standardise_native_ranges(
    native_raw = native_ranges,
    species_col = species_col,
    origin_iso3_col = origin_iso3_col,
    strategy = strategy,
    max_origins = max_origins,
    quiet = quiet
  )
}


#' Attach native/non-native recipient status to a data frame
#'
#' Joins species-level native-origin evidence to an input table and, when a
#' recipient country is available, classifies each row as native, non-native or
#' origin-unknown for that recipient.
#'
#' @details
#' This is the core classifier used by the terrestrial/freshwater and marine
#' pipelines. With `require_country = TRUE`, the helper resolves recipient
#' country information to ISO3 and compares it with each species' native-origin
#' ISO3 set. With `require_country = FALSE`, it attaches species-level origin
#' availability without making recipient-level native/non-native claims.
#'
#' If `native_ranges = NULL` and `use_native_web = TRUE`, provider evidence is
#' compiled through [bf_fetch_native_ranges_web()] before classification. This
#' allows SInAS/GBIF/WoRMS evidence to be generated inside the helper, while still
#' keeping the classification logic centralised here.
#'
#' @param df Input data frame.
#' @param species_col Species column in `df`.
#' @param iso2c_col Optional ISO2 recipient column in `df`.
#' @param iso3c_col Optional ISO3 recipient column in `df`.
#' @param country_col Optional recipient country-name column in `df`. Used only
#'   when `require_country = TRUE` and ISO2/ISO3 columns are unavailable.
#' @param native_ranges Native-origin evidence table.
#' @param native_species_col Species column in `native_ranges`.
#' @param origin_iso3_col Optional origin ISO3 column in `native_ranges`.
#' @param require_country Logical; when `TRUE`, classify recipient status using
#'   recipient ISO2/ISO3. When `FALSE`, attach species-level origin availability.
#' @param strategy Native-origin column strategy passed to
#'   [bf_native_range_lookup()].
#' @param native_filter_mode One of `"audit_only"`, `"non_native_only"`,
#'   `"keep_non_native"`, `"non_native_or_unknown"`,
#'   `"keep_non_native_or_unknown"`, `"native_only"` or `"keep_native"`.
#' @param max_origins Maximum origin count per species.
#' @param use_native_web Logical. If `TRUE` and `native_ranges = NULL`, compile
#'   native-origin evidence from provider helpers before classification.
#' @param native_web_sources Character vector of provider names passed to
#'   `bf_fetch_native_ranges_web()`, commonly `c("sinas", "gbif", "worms")`.
#' @param native_web_cache_dir Cache directory for web/API native-origin evidence.
#'   Must be supplied explicitly when `use_native_web = TRUE`. In examples,
#'   tests and vignettes, use a path under `tempdir()`.
#' @param native_web_force_refresh Logical; refresh cached provider evidence where
#'   supported.
#' @param native_web_sleep_sec Delay in seconds between provider requests where
#'   supported.
#' @param native_web_quiet Logical; suppress provider-specific progress messages.
#' @param native_web_sinas_main_path,native_web_sinas_alllocations_path,native_web_sinas_fulltaxa_path
#'   Optional local SInAS resource paths passed through to the SInAS provider.
#' @param export_native_web_audit Logical; write native-web evidence audit files
#'   when provider compilation is used and the writer helper is available.
#' @param native_web_audit_dir Directory for optional native-web audit files. If
#'   `NULL`, audit files are written to `native_web_cache_dir` when
#'   `export_native_web_audit = TRUE`.
#' @param native_web_audit_prefix File prefix for optional native-web audit files.
#' @param quiet Logical; suppress messages.
#'
#' @return A tibble containing the input rows with native-origin and
#'   recipient-status columns attached. The returned object preserves the input
#'   columns and adds fields such as `native_name_key`, `native_name`,
#'   `native_origin_iso3`, `native_has_origin`, `native_sources_used`,
#'   `native_recipient_iso3`, `native_is_native_recipient`,
#'   `native_is_non_native_recipient`, `native_status`,
#'   `native_filter_keep` and `native_filter_rejection_reason`. When
#'   `require_country = TRUE`, `native_status` classifies each row as
#'   `"native_recipient"`, `"non_native_recipient"` or `"origin_unknown"` by
#'   comparing the recipient ISO3 code with the species-level native-origin ISO3
#'   set. When `require_country = FALSE`, the function attaches species-level
#'   origin availability without making recipient-level native/non-native
#'   claims. Strict `*_only` filter modes return only retained rows; `keep_*`
#'   modes annotate retention decisions without dropping rows.
#'
#' @section Interpretation:
#' `origin_unknown` means that the available native-origin evidence did not
#' confirm a usable origin set for that species. It should not be interpreted as
#' proof that the species is native, non-native, absent or data-deficient in an
#' ecological sense.
#'
#' @section Data source and attribution:
#' This function classifies evidence supplied by other sources. If evidence is
#' generated from SInAS, GBIF, WoRMS or another provider, cite the exact source,
#' version/release and access date used to create `native_ranges`.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_attach_native_status <- function(df,
                                    species_col = "species",
                                    iso2c_col = NULL,
                                    iso3c_col = NULL,
                                    country_col = NULL,
                                    native_ranges = NULL,
                                    native_species_col = "species",
                                    origin_iso3_col = NULL,
                                    require_country = TRUE,
                                    strategy = c("tiered", "union"),
                                    native_filter_mode = c(
                                      "audit_only",
                                      "non_native_only",
                                      "keep_non_native",
                                      "non_native_or_unknown",
                                      "keep_non_native_or_unknown",
                                      "native_only",
                                      "keep_native"
                                    ),
                                    max_origins = 100L,
                                    use_native_web = FALSE,
                                    native_web_sources = c("sinas", "gbif", "worms"),
                                    native_web_cache_dir = NULL,
                                    native_web_force_refresh = FALSE,
                                    native_web_sleep_sec = 0.25,
                                    native_web_quiet = quiet,
                                    native_web_sinas_main_path = NULL,
                                    native_web_sinas_alllocations_path = NULL,
                                    native_web_sinas_fulltaxa_path = NULL,
                                    export_native_web_audit = FALSE,
                                    native_web_audit_dir = NULL,
                                    native_web_audit_prefix = "native_web",
                                    quiet = FALSE) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range helper")

  strategy <- match.arg(strategy)
  native_filter_mode <- match.arg(native_filter_mode)

  if (is.null(df) || !is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  if (!species_col %in% names(df)) {
    stop("`species_col` was not found in `df`: ", species_col, call. = FALSE)
  }

  x <- tibble::as_tibble(df)

  if (is.null(native_ranges)) {
    if (!isTRUE(use_native_web)) {
      stop(
        "`native_ranges` is NULL. Supply a native-origin evidence table, or set ",
        "`use_native_web = TRUE` so biofetchR can compile native evidence from ",
        "the configured web/API sources.",
        call. = FALSE
      )
    }

    if (!exists("bf_fetch_native_ranges_web", mode = "function", inherits = TRUE)) {
      stop(
        "`use_native_web = TRUE`, but `bf_fetch_native_ranges_web()` is not visible. ",
        "Make sure R/utils_native_range_web_sources.R has been sourced/loaded.",
        call. = FALSE
      )
    }

    species_vec <- unique(bf_clean_text(x[[species_col]]))
    species_vec <- species_vec[!is.na(species_vec) & nzchar(species_vec)]

    if (is.null(native_web_cache_dir) || length(native_web_cache_dir) == 0L ||
        !nzchar(trimws(as.character(native_web_cache_dir[[1L]])))) {
      stop(
        "`native_web_cache_dir` must be supplied explicitly when `use_native_web = TRUE`. ",
        "In examples, tests and vignettes, use `file.path(tempdir(), ...)`.",
        call. = FALSE
      )
    }

    native_web_cache_dir <- normalizePath(
      as.character(native_web_cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )

    native_web_args <- list(
      species = species_vec,
      sources = native_web_sources,
      cache_dir = native_web_cache_dir,
      force_refresh = native_web_force_refresh,
      sleep_sec = native_web_sleep_sec,
      quiet = native_web_quiet,
      return = "list"
    )

    # These SInAS path overrides are optional. They are used in tests when the
    # live Zenodo download has already happened once and the parser should reuse
    # those package-managed cached files rather than manually supplied datasets.
    if (!is.null(native_web_sinas_main_path)) {
      native_web_args$sinas_main_path <- native_web_sinas_main_path
    }
    if (!is.null(native_web_sinas_alllocations_path)) {
      native_web_args$sinas_alllocations_path <- native_web_sinas_alllocations_path
    }
    if (!is.null(native_web_sinas_fulltaxa_path)) {
      native_web_args$sinas_fulltaxa_path <- native_web_sinas_fulltaxa_path
    }

    web_evidence <- tryCatch(
      do.call(bf_fetch_native_ranges_web, native_web_args),
      error = function(e) e
    )

    if (inherits(web_evidence, "condition")) {
      stop(
        "Native web evidence compilation failed inside `bf_attach_native_status()`: ",
        conditionMessage(web_evidence),
        call. = FALSE
      )
    }

    if (!is.list(web_evidence) || !"species" %in% names(web_evidence) ||
        !is.data.frame(web_evidence$species)) {
      stop(
        "`bf_fetch_native_ranges_web()` did not return a list with a data-frame ",
        "`species` element.",
        call. = FALSE
      )
    }

    if (isTRUE(export_native_web_audit)) {
      audit_dir <- native_web_audit_dir

      if (is.null(audit_dir) || length(audit_dir) == 0L ||
          !nzchar(trimws(as.character(audit_dir[[1L]])))) {
        audit_dir <- native_web_cache_dir
      } else {
        audit_dir <- normalizePath(
          as.character(audit_dir[[1L]]),
          winslash = "/",
          mustWork = FALSE
        )
      }

      if (exists("bf_write_native_web_outputs", mode = "function", inherits = TRUE)) {
        bf_write_native_web_outputs(
          web_evidence,
          output_dir = audit_dir,
          prefix = native_web_audit_prefix
        )
      }
    }

    native_ranges <- web_evidence$species
    native_species_col <- "species"
    origin_iso3_col <- "native_origin_iso3"
  }

  lookup <- bf_native_range_lookup(
    native_ranges = native_ranges,
    species_col = native_species_col,
    origin_iso3_col = origin_iso3_col,
    strategy = strategy,
    max_origins = max_origins,
    quiet = quiet
  )

  x$native_name_key <- .bf_native_name_key(x[[species_col]])

  out <- dplyr::left_join(
    x,
    lookup,
    by = "native_name_key"
  )

  out$native_has_origin <- out$native_has_origin %in% TRUE

  if (isTRUE(require_country)) {
    if (!is.null(iso3c_col) && iso3c_col %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(iso3c = out[[iso3c_col]], n = nrow(out))
    } else if (!is.null(iso2c_col) && iso2c_col %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(iso2c = out[[iso2c_col]], n = nrow(out))
    } else if (!is.null(country_col) && country_col %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(country = out[[country_col]], n = nrow(out))
    } else if ("iso3c" %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(iso3c = out$iso3c, n = nrow(out))
    } else if ("iso2c" %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(iso2c = out$iso2c, n = nrow(out))
    } else if ("country" %in% names(out)) {
      recipient_iso3 <- .bf_native_country_to_iso3(country = out$country, n = nrow(out))
    } else {
      stop(
        "`require_country = TRUE`, but no recipient ISO2/ISO3/country column was found. ",
        "Supply `iso2c_col`, `iso3c_col` or `country_col`, or use `require_country = FALSE`.",
        call. = FALSE
      )
    }

    origin_list <- .bf_native_split_iso3(out$native_origin_iso3, max_origins = max_origins)

    is_native <- vapply(seq_along(origin_list), function(i) {
      !is.na(recipient_iso3[[i]]) &&
        nzchar(recipient_iso3[[i]]) &&
        recipient_iso3[[i]] %in% origin_list[[i]]
    }, logical(1))

    has_origin <- out$native_has_origin %in% TRUE

    out$native_recipient_iso3 <- recipient_iso3
    out$native_is_native_recipient <- has_origin & is_native
    out$native_is_non_native_recipient <- has_origin & !is_native & !is.na(recipient_iso3) & nzchar(recipient_iso3)

    out$native_status <- dplyr::case_when(
      !has_origin ~ "origin_unknown",
      out$native_is_native_recipient ~ "native_recipient",
      out$native_is_non_native_recipient ~ "non_native_recipient",
      TRUE ~ "origin_unknown"
    )
  } else {
    out$native_recipient_iso3 <- NA_character_
    out$native_is_native_recipient <- FALSE
    out$native_is_non_native_recipient <- FALSE
    out$native_status <- ifelse(
      out$native_has_origin %in% TRUE,
      "origin_available_species_only",
      "origin_unknown"
    )
  }

  native_unknown <- is.na(out$native_status) |
    out$native_status %in% c(
      "origin_unknown",
      "recipient_country_missing",
      "cosmopolitan_or_uncertain",
      "origin_available_species_only"
    ) |
    !(out$native_has_origin %in% TRUE)

  # `native_filter_keep` is an annotation column used by the main pipelines and
  # the shared origin-gate logic. The `keep_*` modes intentionally annotate rows
  # rather than dropping them, so pipeline audit tables can still record accepted
  # and rejected rows before GBIF submission.
  out$native_filter_keep <- switch(
    native_filter_mode,
    audit_only = rep(TRUE, nrow(out)),
    non_native_only = out$native_is_non_native_recipient %in% TRUE,
    keep_non_native = out$native_is_non_native_recipient %in% TRUE,
    non_native_or_unknown = (out$native_is_non_native_recipient %in% TRUE) | native_unknown,
    keep_non_native_or_unknown = (out$native_is_non_native_recipient %in% TRUE) | native_unknown,
    native_only = out$native_is_native_recipient %in% TRUE,
    keep_native = out$native_is_native_recipient %in% TRUE
  )

  out$native_filter_rejection_reason <- ifelse(
    out$native_filter_keep %in% TRUE,
    NA_character_,
    switch(
      native_filter_mode,
      audit_only = NA_character_,
      non_native_only = "not_confirmed_non_native",
      keep_non_native = "not_confirmed_non_native",
      non_native_or_unknown = "not_non_native_or_unknown",
      keep_non_native_or_unknown = "not_non_native_or_unknown",
      native_only = "not_confirmed_native",
      keep_native = "not_confirmed_native"
    )
  )

  # Strict `*_only` modes are public convenience filters. `keep_*` modes are
  # pipeline-safe annotation modes and should not drop rows.
  if (native_filter_mode == "non_native_only") {
    out <- out[out$native_filter_keep %in% TRUE, , drop = FALSE]
  }

  if (native_filter_mode == "non_native_or_unknown") {
    out <- out[out$native_filter_keep %in% TRUE, , drop = FALSE]
  }

  if (native_filter_mode == "native_only") {
    out <- out[out$native_filter_keep %in% TRUE, , drop = FALSE]
  }

  tibble::as_tibble(out)
}


#' Filter a data frame to confirmed non-native recipient records
#'
#' @inheritParams bf_attach_native_status
#'
#' @return A tibble containing only rows classified by
#'   [bf_attach_native_status()] as confirmed non-native recipients under the
#'   supplied native-origin evidence and recipient-country settings. The returned
#'   object includes the original input columns plus the native-origin/status
#'   columns added by [bf_attach_native_status()], including
#'   `native_is_non_native_recipient`, `native_status`,
#'   `native_filter_keep` and `native_filter_rejection_reason`. Rows excluded by
#'   this filter are not necessarily native; they simply lack confirmed
#'   non-native recipient evidence under the selected criteria.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_filter_non_native <- function(df,
                                 species_col = "species",
                                 iso2c_col = NULL,
                                 iso3c_col = NULL,
                                 country_col = NULL,
                                 native_ranges,
                                 native_species_col = "species",
                                 origin_iso3_col = NULL,
                                 require_country = TRUE,
                                 strategy = c("tiered", "union"),
                                 max_origins = 100L,
                                 quiet = FALSE) {
  bf_attach_native_status(
    df = df,
    species_col = species_col,
    iso2c_col = iso2c_col,
    iso3c_col = iso3c_col,
    country_col = country_col,
    native_ranges = native_ranges,
    native_species_col = native_species_col,
    origin_iso3_col = origin_iso3_col,
    require_country = require_country,
    strategy = strategy,
    native_filter_mode = "non_native_only",
    max_origins = max_origins,
    quiet = quiet
  )
}


#' Filter a data frame to confirmed native recipient records
#'
#' @inheritParams bf_attach_native_status
#'
#' @return A tibble containing only rows classified by
#'   [bf_attach_native_status()] as confirmed native recipients under the
#'   supplied native-origin evidence and recipient-country settings. The returned
#'   object includes the original input columns plus the native-origin/status
#'   columns added by [bf_attach_native_status()], including
#'   `native_is_native_recipient`, `native_status`, `native_filter_keep` and
#'   `native_filter_rejection_reason`. Rows excluded by this filter are not
#'   necessarily non-native; they simply lack confirmed native recipient evidence
#'   under the selected criteria.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_filter_native <- function(df,
                             species_col = "species",
                             iso2c_col = NULL,
                             iso3c_col = NULL,
                             country_col = NULL,
                             native_ranges,
                             native_species_col = "species",
                             origin_iso3_col = NULL,
                             require_country = TRUE,
                             strategy = c("tiered", "union"),
                             max_origins = 100L,
                             quiet = FALSE) {
  bf_attach_native_status(
    df = df,
    species_col = species_col,
    iso2c_col = iso2c_col,
    iso3c_col = iso3c_col,
    country_col = country_col,
    native_ranges = native_ranges,
    native_species_col = native_species_col,
    origin_iso3_col = origin_iso3_col,
    require_country = require_country,
    strategy = strategy,
    native_filter_mode = "native_only",
    max_origins = max_origins,
    quiet = quiet
  )
}


#' Summarise native-origin status columns
#'
#' @param x Data frame returned by [bf_attach_native_status()].
#'
#' @return A tibble summarising the `native_status` column in `x`. The output
#'   contains one row per native-status class and two columns: `native_status`,
#'   the character status label, and `n`, the number of rows assigned to that
#'   status. Rows are ordered from the most frequent status to the least
#'   frequent, then alphabetically by status label. The summary is intended for
#'   quick auditing of native-origin classification outcomes.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_native_status_summary <- function(x) {
  bf_require_packages(c("tibble", "dplyr"), context = "native-range helper")

  if (is.null(x) || !is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }

  if (!"native_status" %in% names(x)) {
    stop("`x` must contain a `native_status` column.", call. = FALSE)
  }

  tibble::as_tibble(x) |>
    dplyr::count(.data$native_status, name = "n") |>
    dplyr::arrange(dplyr::desc(.data$n), .data$native_status)
}


#' Reconcile GRIIS and native-origin evidence
#'
#' @param x Data frame containing GRIIS columns and native-status columns.
#'
#' @return A tibble containing the input rows with one additional character
#'   column, `biofetchr_origin_evidence_status`. The added column summarises
#'   agreement, conflict or missingness between available GRIIS evidence and
#'   native-origin recipient evidence, with labels such as
#'   `"griis_invasive_native_origin_supported"`,
#'   `"griis_invasive_native_origin_conflict"`,
#'   `"native_origin_non_native_only"` and `"origin_unknown"`. The returned
#'   table preserves the input columns and should be interpreted as an audit
#'   table for evidence consistency, not as an independent biological
#'   classification.
#'
#' @section Interpretation:
#' Reconciliation labels describe agreement or conflict between available GRIIS
#' and native-origin evidence. They are audit labels, not independent biological
#' proof of invasion status.
#'
#' @family native-origin status helpers
#' @md
#' @export
bf_reconcile_griis_native_status <- function(x) {
  bf_require_packages("tibble", context = "native-range helper")

  if (is.null(x) || !is.data.frame(x)) {
    stop("`x` must be a data frame.", call. = FALSE)
  }

  out <- tibble::as_tibble(x)

  griis_listed <- if ("griis_listed" %in% names(out)) {
    out$griis_listed %in% TRUE
  } else if ("griis_status" %in% names(out)) {
    !is.na(out$griis_status) & out$griis_status != "not_listed"
  } else {
    rep(FALSE, nrow(out))
  }

  griis_invasive <- if ("griis_invasive" %in% names(out)) {
    out$griis_invasive %in% TRUE
  } else if ("griis_is_invasive" %in% names(out)) {
    out$griis_is_invasive %in% TRUE
  } else {
    rep(FALSE, nrow(out))
  }

  native_has_origin <- if ("native_has_origin" %in% names(out)) {
    out$native_has_origin %in% TRUE
  } else {
    rep(FALSE, nrow(out))
  }

  native_is_native <- if ("native_is_native_recipient" %in% names(out)) {
    out$native_is_native_recipient %in% TRUE
  } else {
    rep(FALSE, nrow(out))
  }

  native_is_non_native <- if ("native_is_non_native_recipient" %in% names(out)) {
    out$native_is_non_native_recipient %in% TRUE
  } else {
    rep(FALSE, nrow(out))
  }

  out$biofetchr_origin_evidence_status <- dplyr::case_when(
    griis_invasive & native_is_non_native ~ "griis_invasive_native_origin_supported",
    griis_invasive & native_is_native ~ "griis_invasive_native_origin_conflict",
    griis_invasive & !native_has_origin ~ "griis_invasive_origin_unknown",
    griis_listed & native_is_non_native ~ "griis_listed_native_origin_supported",
    griis_listed & native_is_native ~ "griis_listed_native_origin_conflict",
    griis_listed & !native_has_origin ~ "griis_listed_origin_unknown",
    !griis_listed & native_is_non_native ~ "native_origin_non_native_only",
    !griis_listed & native_is_native ~ "native_origin_native_recipient_only",
    !griis_listed & native_has_origin ~ "native_origin_species_only",
    TRUE ~ "origin_unknown"
  )

  out
}
