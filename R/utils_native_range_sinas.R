################################################################################
# utils_native_range_sinas.R
# ------------------------------------------------------------------------------
# biofetchR: SInAS native-range evidence provider
# ------------------------------------------------------------------------------
#
# PURPOSE
#   Add SInAS as a first-class native-origin evidence source for biofetchR. This
#   provider reads local or package-resolved SInAS 3.1.1 resources, converts
#   SInAS native-location evidence into ISO3 country codes, and returns the same
#   evidence contract used by the wider native-origin workflow.
#
# RELATIONSHIP TO utils_native_range.R
#   This file is not a replacement for utils_native_range.R and should not be
#   treated as a duplicate of it. This file is provider-specific: it knows how to
#   parse SInAS tables, detect native-like SInAS establishment records and build
#   SInAS-derived origin evidence. utils_native_range.R is the core classifier:
#   it accepts provider evidence from SInAS/GBIF/WoRMS or user tables and decides
#   whether a recipient country is native, non-native or unresolved.
#
#   In normal package use the flow is:
#     bf_fetch_native_ranges_sinas()  -> provider evidence table
#     bf_attach_native_status()       -> recipient-level native/non-native status
#
# SCOPE
#   This file parses SInAS resources that are already local or have been resolved
#   by bf_download_sinas_resources(). It does not submit GBIF occurrence
#   downloads, does not perform spatial overlay assignment, and does not classify
#   recipient-country native/non-native status by itself.
#
# OUTPUT CONTRACT
#   The exported bf_fetch_native_ranges_sinas() returns evidence that can feed
#   directly into bf_attach_native_status() via the native_ranges argument:
#     species
#     native_origin_iso3
#     native_sources_used
#     native_has_origin
#
# DATA USE AND ATTRIBUTION
#   SInAS records remain third-party data. biofetchR resolves, caches and
#   standardises those resources but does not own or redistribute them. Users
#   should cite the exact SInAS release/Zenodo record used in their workflow and
#   follow any licence or attribution requirements attached to that release.
#
################################################################################


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' Extract a simple canonical binomial
#'
#' @param x Character vector of scientific names.
#'
#' @return Character vector containing `Genus species` where detectable.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_binomial <- function(x) {
  x <- bf_clean_text(x)
  x <- gsub("x", "x", x, fixed = TRUE)
  x <- gsub("\\s*\\([^)]*\\)", " ", x)
  x <- gsub("[[:space:]]+", " ", x)
  x <- trimws(x)

  m <- regexec("([A-Z][A-Za-z.-]+)\\s+([a-z][A-Za-z.-]+)", x)
  parts <- regmatches(x, m)

  out <- vapply(parts, function(z) {
    if (length(z) >= 3L) paste(z[[2]], z[[3]]) else NA_character_
  }, character(1))

  out[!nzchar(out)] <- NA_character_
  out
}


#' Return a normalised key for species-name matching
#'
#' @param x Species names.
#'
#' @return Lowercase canonical-binomial keys.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_name_key <- function(x) {
  if (exists(".bf_native_name_key", mode = "function", inherits = TRUE)) {
    return(.bf_native_name_key(x))
  }

  x <- .bf_sinas_native_binomial(x)
  x <- tolower(trimws(as.character(x)))
  x[x == "" | is.na(x)] <- NA_character_
  x
}


#' Detect native-like SInAS establishment text
#'
#' SInAS establishment fields can contain native-like labels as well as explicit
#' non-native labels. This helper keeps native/endemic/indigenous signals while
#' excluding introduced, alien, invasive, cultivated and similar signals.
#'
#' @param x Establishment or status text.
#'
#' @return Logical vector.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_is_native_like <- function(x) {
  x <- toupper(bf_clean_text(x))
  x[is.na(x)] <- ""

  native_hit <- grepl(
    "\\b(NATIVE|ENDEMIC|INDIGENOUS|AUTOCHTHONOUS)\\b",
    x
  )

  alien_hit <- grepl(
    "\\b(NON[- ]?NATIVE|ALIEN|INTRODUCED|INVASIVE|EXOTIC|NATURALI[ZS]ED|ADVENTIVE|CULTIVATED|DOMESTICATED)\\b",
    x
  )

  native_hit & !alien_hit
}


#' Clean column names for schema matching
#'
#' @param x Column names.
#'
#' @return Lowercase, underscore-separated names.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_clean_colname <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}


#' Pick the first available column from candidate names
#'
#' @param dat Data frame.
#' @param candidates Candidate column names.
#'
#' @return Matched raw column name or `NA_character_`.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_pick_col <- function(dat, candidates) {
  if (is.null(dat) || !is.data.frame(dat)) return(NA_character_)

  nms <- names(dat)
  nms_clean <- .bf_sinas_native_clean_colname(nms)
  cand_clean <- .bf_sinas_native_clean_colname(candidates)

  # Important: respect the order of `candidates`, not the order of columns in
  # the data. SInAS has `occurrenceStatus` before `establishmentMeans`, and both
  # can be plausible status-like fields. If we return the first matching data
  # column, the parser incorrectly chooses occurrenceStatus = "present" instead
  # of establishmentMeans = "native"/"introduced", which removes all native rows.
  for (cand in cand_clean) {
    hit <- which(nms_clean == cand)
    if (length(hit)) return(nms[hit[[1]]])
  }

  NA_character_
}


#' Read a SInAS-related table from CSV/TSV/XLSX
#'
#' @param path File path.
#' @param combine_sheets Logical; for XLSX files, read all sheets and bind rows.
#'
#' @return A tibble.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_read_table <- function(path, combine_sheets = TRUE) {
  if (!file.exists(path)) {
    stop("SInAS table does not exist: ", path, call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    bf_require_packages(c("readxl", "tibble", "dplyr"), context = "SInAS native-range helper")

    sheets <- readxl::excel_sheets(path)
    if (!length(sheets)) {
      stop("No sheets found in SInAS workbook: ", path, call. = FALSE)
    }

    if (!isTRUE(combine_sheets)) {
      return(tibble::as_tibble(readxl::read_excel(path, sheet = sheets[[1]], col_types = "text")))
    }

    parts <- lapply(sheets, function(sh) {
      x <- tibble::as_tibble(readxl::read_excel(path, sheet = sh, col_types = "text"))
      x$sheet <- sh
      x
    })

    return(dplyr::bind_rows(parts))
  }

  bf_require_packages("tibble", context = "SInAS native-range helper")

  # SInAS 3.1.1 uses a quoted, whitespace-separated text format despite the
  # .csv extension, e.g.
  #
  #   "location" "locationID" "taxon" "taxonID" ...
  #
  # readr::read_csv() and readr::read_table() are not reliable for this file:
  # read_csv() reads the whole header as one column, while read_table() can
  # generate many row-width parsing failures. Base R read.table(sep = "") is the
  # most reliable here because it treats whitespace as the delimiter while
  # respecting quoted multi-word fields.
  read_base_whitespace <- function() {
    tryCatch(
      tibble::as_tibble(utils::read.table(
        file = path,
        sep = "",
        header = TRUE,
        stringsAsFactors = FALSE,
        fill = FALSE,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        na.strings = c("NA"),
        blank.lines.skip = TRUE,
        allowEscapes = FALSE
      )),
      error = function(e) e
    )
  }

  clean_names <- function(x) .bf_sinas_native_clean_colname(names(x))

  has_expected_sinas_shape <- function(x) {
    if (inherits(x, "condition") || is.null(x) || !is.data.frame(x) || ncol(x) <= 1L) {
      return(FALSE)
    }

    nms <- clean_names(x)

    has_taxon <- any(nms %in% .bf_sinas_native_clean_colname(c(
      "taxon", "scientificName", "scientific_name", "species", "taxonName",
      "acceptedNameUsage", "canonicalName", "verbatimScientificName"
    )))

    has_loc <- any(nms %in% .bf_sinas_native_clean_colname(c(
      "locationID", "location_id", "locationid", "areaID", "area_id",
      "regionID", "region_id", "locationCode"
    ))) || any(nms %in% .bf_sinas_native_clean_colname(c(
      "location", "location_name", "country", "area", "locality",
      "countryCode", "country_code", "region"
    )))

    has_status_or_taxalist <- any(nms %in% .bf_sinas_native_clean_colname(c(
      "establishmentMeans", "establishment_means", "establishment", "status",
      "nativeStatus", "native_status", "nativity", "occurrenceStatus",
      "occurrence_status", "degreeOfEstablishment",
      "degree_of_establishment"
    ))) || any(nms %in% .bf_sinas_native_clean_colname(c(
      "GBIFstatus", "GBIFstatus_Synonym", "GBIFmatchtype", "GBIFtaxonRank",
      "taxonID", "taxaGroup"
    )))

    has_taxon && has_loc && has_status_or_taxalist
  }

  base_ws <- read_base_whitespace()
  if (has_expected_sinas_shape(base_ws)) {
    return(base_ws)
  }

  # Fallbacks for future SInAS releases that may use conventional delimiters.
  parsed <- list(base_whitespace = base_ws)

  if (requireNamespace("readr", quietly = TRUE)) {
    for (delim in c(",", ";", "\t", "|")) {
      nm <- paste0("readr_delim_", if (delim == "\t") "tab" else delim)
      parsed[[nm]] <- tryCatch(
        tibble::as_tibble(suppressWarnings(readr::read_delim(
          file = path,
          delim = delim,
          col_types = readr::cols(.default = readr::col_character()),
          show_col_types = FALSE,
          progress = FALSE,
          locale = readr::locale(encoding = "UTF-8"),
          guess_max = 10000
        ))),
        error = function(e) e
      )
    }

    parsed$readr_csv <- tryCatch(
      tibble::as_tibble(suppressWarnings(readr::read_csv(
        file = path,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress = FALSE,
        locale = readr::locale(encoding = "UTF-8"),
        guess_max = 10000
      ))),
      error = function(e) e
    )
  }

  for (sep in c(",", ";", "\t", "|")) {
    nm <- paste0("base_sep_", if (sep == "\t") "tab" else sep)
    parsed[[nm]] <- tryCatch(
      tibble::as_tibble(utils::read.table(
        file = path,
        sep = sep,
        header = TRUE,
        stringsAsFactors = FALSE,
        fill = TRUE,
        quote = "\"",
        comment.char = "",
        check.names = FALSE,
        na.strings = c("NA"),
        blank.lines.skip = TRUE
      )),
      error = function(e) e
    )
  }

  score_table <- function(x) {
    if (!has_expected_sinas_shape(x)) return(-Inf)

    nms <- clean_names(x)

    target_hits <- sum(c(
      any(nms %in% .bf_sinas_native_clean_colname(c("taxon", "scientificName", "species"))),
      any(nms %in% .bf_sinas_native_clean_colname(c("locationID", "location", "country"))),
      any(nms %in% .bf_sinas_native_clean_colname(c("establishmentMeans", "status", "occurrenceStatus", "GBIFstatus")))
    ))

    (1000 * target_hits) + ncol(x)
  }

  scores <- vapply(parsed, score_table, numeric(1))
  if (any(is.finite(scores))) {
    return(parsed[[which.max(scores)]])
  }

  err <- if (inherits(base_ws, "condition")) conditionMessage(base_ws) else "no valid parse"
  stop(
    "Could not read SInAS table with base whitespace, CSV, semicolon, TSV or pipe parsers: ",
    path,
    ". Base whitespace parse detail: ",
    err,
    call. = FALSE
  )
}


#' Convert country/place information to ISO3
#'
#' @param country Optional country/place names.
#' @param iso2c Optional ISO2 codes.
#' @param iso3c Optional ISO3 codes.
#' @param n Expected output length.
#'
#' @return Character vector of ISO3 codes.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_country_to_iso3 <- function(country = NULL,
                                             iso2c = NULL,
                                             iso3c = NULL,
                                             n = NULL) {
  if (exists(".bf_native_country_to_iso3", mode = "function", inherits = TRUE)) {
    return(.bf_native_country_to_iso3(country = country, iso2c = iso2c, iso3c = iso3c, n = n))
  }

  if (exists(".bf_nrw_country_to_iso3", mode = "function", inherits = TRUE)) {
    return(.bf_nrw_country_to_iso3(country = country, iso2c = iso2c, iso3c = iso3c, n = n))
  }

  if (is.null(n)) {
    n <- max(
      length(bf_null_coalesce(country, character(0))),
      length(bf_null_coalesce(iso2c, character(0))),
      length(bf_null_coalesce(iso3c, character(0))),
      0L
    )
  }

  if (n == 0L) return(character(0))

  recycle_to_n <- function(x) {
    if (is.null(x)) return(rep(NA_character_, n))
    x <- as.character(x)
    if (length(x) == 0L) return(rep(NA_character_, n))
    if (length(x) == 1L && n > 1L) return(rep(x, n))
    if (length(x) != n) {
      stop("Country/ISO vector length mismatch in SInAS helper.", call. = FALSE)
    }
    x
  }

  out <- rep(NA_character_, n)

  if (!is.null(iso3c)) {
    x <- toupper(trimws(recycle_to_n(iso3c)))
    x[x == ""] <- NA_character_
    out[grepl("^[A-Z]{3}$|^XKX$", x)] <- x[grepl("^[A-Z]{3}$|^XKX$", x)]
  }

  if (!is.null(iso2c) && requireNamespace("countrycode", quietly = TRUE)) {
    x <- toupper(trimws(recycle_to_n(iso2c)))
    y <- suppressWarnings(countrycode::countrycode(
      x,
      origin = "iso2c",
      destination = "iso3c",
      warn = FALSE
    ))
    fill <- is.na(out) & !is.na(y)
    out[fill] <- toupper(y[fill])
  }

  if (!is.null(country) && requireNamespace("countrycode", quietly = TRUE)) {
    x <- bf_clean_text(recycle_to_n(country))
    y <- suppressWarnings(countrycode::countrycode(
      x,
      origin = "country.name",
      destination = "iso3c",
      warn = FALSE
    ))
    fill <- is.na(out) & !is.na(y)
    out[fill] <- toupper(y[fill])
  }

  out
}


#' Build a locationID-to-ISO3 crosswalk from AllLocations
#'
#' @param alllocations Data frame read from AllLocations.
#'
#' @return Tibble with `locationID`, `iso3c` and `raw_location`.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_build_location_crosswalk <- function(alllocations) {
  bf_require_packages(c("tibble", "dplyr"), context = "SInAS native-range helper")

  dat <- tibble::as_tibble(alllocations)

  loc_col <- .bf_sinas_native_pick_col(dat, c("locationID", "location_id", "locationid"))
  if (is.na(loc_col)) {
    stop("AllLocations is missing a `locationID` column.", call. = FALSE)
  }

  iso3_col <- .bf_sinas_native_pick_col(dat, c("ISO3", "iso3", "iso3c", "country_iso3"))
  iso2_col <- .bf_sinas_native_pick_col(dat, c("ISO_1", "ISO2", "iso2", "iso2c", "country_code"))
  name_col <- .bf_sinas_native_pick_col(
    dat,
    c(
      "gadm0_name", "glonaf_country", "country", "Country", "country_name",
      "location", "Location", "admin0", "admin0_name", "name"
    )
  )

  n <- nrow(dat)

  iso3 <- .bf_sinas_native_country_to_iso3(
    country = if (!is.na(name_col)) dat[[name_col]] else NULL,
    iso2c = if (!is.na(iso2_col)) dat[[iso2_col]] else NULL,
    iso3c = if (!is.na(iso3_col)) dat[[iso3_col]] else NULL,
    n = n
  )

  raw_location <- if (!is.na(name_col)) {
    bf_clean_text(dat[[name_col]])
  } else {
    bf_clean_text(dat[[loc_col]])
  }

  tibble::tibble(
    locationID = bf_clean_text(dat[[loc_col]]),
    iso3c = toupper(iso3),
    raw_location = raw_location
  ) |>
    dplyr::filter(!is.na(.data$locationID), nzchar(.data$locationID)) |>
    dplyr::distinct(.data$locationID, .data$iso3c, .data$raw_location)
}


#' Return an empty SInAS long-evidence table
#'
#' @return Empty tibble with native-evidence columns.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_empty_long <- function() {
  bf_require_packages("tibble", context = "SInAS native-range helper")

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


#' Collapse SInAS long evidence to species-level output
#'
#' @param long Long SInAS evidence table.
#'
#' @return Species-level tibble.
#'
#' @keywords internal
#' @noRd
.bf_sinas_native_collapse_species <- function(long) {
  if (exists(".bf_nrw_collapse_species", mode = "function", inherits = TRUE)) {
    return(.bf_nrw_collapse_species(long))
  }

  bf_require_packages(c("tibble", "dplyr"), context = "SInAS native-range helper")

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

  rows <- lapply(split(long, long$species), function(x) {
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
# Public SInAS evidence helper
# -----------------------------------------------------------------------------

#' Fetch native-origin evidence from SInAS
#'
#' Reads SInAS 3.1.1 native-location records, converts SInAS `locationID` values
#' to ISO3 countries using `AllLocations`, and returns native-origin evidence in
#' the same structure as the optional native web/API helpers. When local paths
#' are not supplied, the function uses `bf_download_sinas_resources()` to resolve
#' package-managed cached copies from Zenodo.
#'
#' @details
#' This is a **provider-specific** helper. It extracts native-origin evidence
#' from SInAS and returns species-level origin countries, but it does not decide
#' whether a recipient country is native or non-native. Pass the returned
#' species-level table to [bf_attach_native_status()] for recipient-level
#' classification.
#'
#' SInAS 3.1.1 can use a quoted whitespace-delimited text format despite the
#' `.csv` extension. The internal reader therefore tries a SInAS-aware parser
#' before falling back to more conventional CSV, TSV and pipe-delimited readers.
#'
#' @section Relationship to `utils_native_range.R`:
#' Keep this script alongside `utils_native_range.R`. This file provides the
#' SInAS evidence source; `utils_native_range.R` provides the general
#' standardisation, attachment, filtering and reconciliation logic used by the
#' terrestrial/freshwater and marine pipelines.
#'
#' @section Data source and attribution:
#' biofetchR resolves and standardises SInAS resources but does not own or
#' redistribute the underlying data. Users should cite the exact SInAS release,
#' Zenodo record, version and access date used in their workflow, and should
#' follow the licence terms attached to that release.
#'
#' @param species Character vector of species names, or a data frame containing
#'   a species column.
#' @param species_col Species column name when `species` is a data frame.
#' @param cache_dir Cache directory used when SInAS resources need to be
#'   downloaded. Must be supplied explicitly unless all three local resource
#'   paths are supplied via `main_path`, `alllocations_path` and `fulltaxa_path`.
#'   In examples, tests and vignettes, use a path under `tempdir()`.
#' @param force_refresh Logical; re-download/re-read package-managed resources.
#' @param quiet Logical; suppress progress messages.
#' @param main_path Optional local path to `SInAS_3.1.1.csv`.
#' @param alllocations_path Optional local path to `AllLocations.xlsx`, `.csv`,
#'   or `.tsv`.
#' @param fulltaxa_path Optional local path to
#'   `SInAS_3.1.1_FullTaxaList.csv`.
#' @param record_id Zenodo record identifier used by
#'   `bf_download_sinas_resources()`.
#' @param main_url Optional direct URL for the main SInAS CSV.
#' @param fulltaxa_url Optional direct URL for the FullTaxaList CSV.
#' @param config_zip_url Optional direct URL for the config archive containing
#'   `AllLocations.xlsx`.
#' @param return One of `"long"`, `"species"` or `"list"`.
#'
#' @return The returned object depends on `return`. If `return = "long"`, the
#'   function returns a tibble with one row per retained SInAS native-origin
#'   evidence record, including `species`, `source`, `accepted_name`,
#'   `source_taxon_id`, `raw_native_area`, `raw_status`, `origin_iso3`,
#'   `evidence_type` and `source_url`. If `return = "species"`, the function
#'   returns a species-level tibble with collapsed native-origin evidence,
#'   including `species`, `native_origin_iso3`, `native_sources_used`,
#'   `native_web_unmapped_strings`, `native_web_n_records` and
#'   `native_has_origin`. If `return = "list"`, the function returns a named
#'   list with `long`, `species`, `unmapped` and `summary` elements. The
#'   `unmapped` element contains long-format evidence rows without resolved ISO3
#'   origin codes, and `summary` is a one-row tibble reporting input species
#'   counts, matched SInAS records, ISO3-resolved species and unmapped records.
#'   These outputs provide species-level native-origin evidence only; recipient-
#'   level native/non-native classification is performed by
#'   [bf_attach_native_status()].
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   sinas_native <- bf_fetch_native_ranges_sinas(
#'     species = c("Carcinus maenas", "Ficopomatus enigmaticus"),
#'     cache_dir = file.path(tempdir(), "biofetchR_sinas"),
#'     return = "species",
#'     quiet = FALSE
#'   )
#'
#'   sinas_native
#' }
#' }
#'
#' @seealso [bf_attach_native_status()], [bf_download_sinas_resources()]
#' @family native-origin evidence providers
#' @md
#' @export
bf_fetch_native_ranges_sinas <- function(species,
                                         species_col = "species",
                                         cache_dir = NULL,
                                         force_refresh = FALSE,
                                         quiet = FALSE,
                                         main_path = NULL,
                                         alllocations_path = NULL,
                                         fulltaxa_path = NULL,
                                         record_id = "18220953",
                                         main_url = NULL,
                                         fulltaxa_url = NULL,
                                         config_zip_url = NULL,
                                         return = c("long", "species", "list")) {
  bf_require_packages(c("tibble", "dplyr"), context = "SInAS native-range helper")

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
    long <- .bf_sinas_native_empty_long()
    species_out <- .bf_sinas_native_collapse_species(long)
    unmapped <- long
    summary <- tibble::tibble(
      n_species_input = 0L,
      n_species_with_sinas_records = 0L,
      n_species_with_iso3 = 0L,
      n_long_records = 0L,
      n_unmapped_records = 0L,
      sources = "SInAS",
      sources_used = "SInAS"
    )
    out <- list(long = long, species = species_out, unmapped = unmapped, summary = summary)
    if (return == "species") return(out$species)
    if (return == "long") return(out$long)
    return(out)
  }

  urls <- if (exists("bf_sinas_default_urls", mode = "function", inherits = TRUE)) {
    bf_sinas_default_urls(record_id)
  } else {
    base <- paste0("https://zenodo.org/records/", record_id, "/files/")
    c(
      main_csv = paste0(base, "SInAS_3.1.1.csv?download=1"),
      fulltaxa_csv = paste0(base, "SInAS_3.1.1_FullTaxaList.csv?download=1"),
      config_zip = paste0(base, "All_Config_Files_SInAS_v3.1.1.zip?download=1")
    )
  }

  main_url <- bf_null_coalesce(main_url, urls[["main_csv"]])

  # Keep NULL by default.
  # The default SInAS 3.1.1 Zenodo record does not expose FullTaxaList as a
  # standalone CSV. Leaving this as NULL makes bf_download_sinas_resources()
  # recover FullTaxaList from All_Output_Files_SInAS_v3.1.1.zip instead.
  fulltaxa_url <- fulltaxa_url

  config_zip_url <- bf_null_coalesce(config_zip_url, urls[["config_zip"]])

  if (is.null(main_path) || is.null(alllocations_path) || is.null(fulltaxa_path)) {
    if (!exists("bf_download_sinas_resources", mode = "function", inherits = TRUE)) {
      stop(
        "SInAS resources were not supplied and `bf_download_sinas_resources()` is not available. ",
        "Add `R/utils_sinas_download.R` to the package or supply `main_path`, ",
        "`alllocations_path` and `fulltaxa_path`.",
        call. = FALSE
      )
    }

    if (is.null(cache_dir) || length(cache_dir) == 0L ||
        !nzchar(trimws(as.character(cache_dir[[1L]])))) {
      stop(
        "`cache_dir` must be supplied explicitly when SInAS resources need to be downloaded. ",
        "Alternatively, supply `main_path`, `alllocations_path` and `fulltaxa_path`.",
        call. = FALSE
      )
    }

    cache_dir <- normalizePath(
      as.character(cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )

    resolved <- bf_download_sinas_resources(
      cache_dir = file.path(cache_dir, "native_range_sinas"),
      force_refresh = force_refresh,
      quiet = quiet,
      main_path = main_path,
      allloc_path = alllocations_path,
      fulltaxa_path = fulltaxa_path,
      record_id = record_id,
      main_url = main_url,
      fulltaxa_url = fulltaxa_url,
      config_zip_url = config_zip_url
    )

    main_path <- resolved$main_csv
    alllocations_path <- resolved$alllocations_xlsx
    fulltaxa_path <- resolved$fulltaxa_csv
  }

  .bf_msg("Reading SInAS native-origin evidence.", quiet = quiet)

  sinas <- .bf_sinas_native_read_table(main_path, combine_sheets = FALSE)
  allloc <- .bf_sinas_native_read_table(alllocations_path, combine_sheets = TRUE)

  taxon_col <- .bf_sinas_native_pick_col(sinas, c("taxon", "scientificName", "scientific_name", "species", "taxonName", "acceptedNameUsage", "canonicalName", "originalNameUsage", "verbatimScientificName"))
  taxon_id_col <- .bf_sinas_native_pick_col(sinas, c("taxonID", "taxon_id", "id"))
  loc_id_col <- .bf_sinas_native_pick_col(sinas, c("locationID", "location_id", "locationid", "areaID", "area_id", "regionID", "region_id", "locationCode"))
  loc_col <- .bf_sinas_native_pick_col(sinas, c("location", "location_name", "country", "area", "locality", "countryCode", "country_code", "region"))
  status_col <- .bf_sinas_native_pick_col(sinas, c("establishmentMeans", "establishment_means", "establishment", "status", "nativeStatus", "native_status", "nativity", "occurrenceStatus", "occurrence_status", "degreeOfEstablishment", "degree_of_establishment"))

  missing_cols <- c(
    if (is.na(taxon_col)) "taxon",
    if (is.na(loc_id_col)) "locationID",
    if (is.na(status_col)) "establishmentMeans/status"
  )

  if (length(missing_cols)) {
    stop(
      "SInAS table is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ". Parsed columns were: ",
      paste(utils::head(names(sinas), 80), collapse = "; "),
      ". The SInAS file should be read as quoted whitespace-separated text.",
      call. = FALSE
    )
  }

  crosswalk <- .bf_sinas_native_build_location_crosswalk(allloc)

  sinas_tbl <- tibble::as_tibble(sinas) |>
    dplyr::mutate(
      .bf_sinas_taxon = bf_clean_text(.data[[taxon_col]]),
      .bf_sinas_taxon_key = .bf_sinas_native_name_key(.data[[taxon_col]]),
      .bf_sinas_location_id = bf_clean_text(.data[[loc_id_col]]),
      .bf_sinas_location = if (!is.na(loc_col)) bf_clean_text(.data[[loc_col]]) else .bf_sinas_location_id,
      .bf_sinas_status = bf_clean_text(.data[[status_col]]),
      .bf_sinas_taxon_id = if (!is.na(taxon_id_col)) as.character(.data[[taxon_id_col]]) else NA_character_
    ) |>
    dplyr::filter(.bf_sinas_native_is_native_like(.data$.bf_sinas_status)) |>
    dplyr::left_join(crosswalk, by = c(".bf_sinas_location_id" = "locationID"))

  requested <- tibble::tibble(
    species = species_vec,
    .bf_requested_key = .bf_sinas_native_name_key(species_vec)
  ) |>
    dplyr::filter(!is.na(.data$.bf_requested_key), nzchar(.data$.bf_requested_key)) |>
    dplyr::distinct(.data$species, .data$.bf_requested_key)

  matched <- requested |>
    dplyr::left_join(
      sinas_tbl,
      by = c(".bf_requested_key" = ".bf_sinas_taxon_key")
    )

  # Optional FullTaxaList expansion. This lets aliases inherit the accepted
  # SInAS native-origin evidence where SInAS provides an accepted-name link.
  if (!is.null(fulltaxa_path) && file.exists(fulltaxa_path)) {
    ft <- tryCatch(.bf_sinas_native_read_table(fulltaxa_path, combine_sheets = FALSE), error = function(e) NULL)

    if (!is.null(ft) && nrow(ft)) {
      alias_col <- .bf_sinas_native_pick_col(ft, c("scientificName", "scientific_name", "alias", "name"))
      accepted_col <- .bf_sinas_native_pick_col(ft, c("taxon", "acceptedName", "accepted_name", "acceptedNameUsage"))

      if (!is.na(alias_col) && !is.na(accepted_col)) {
        alias_tbl <- tibble::as_tibble(ft) |>
          dplyr::mutate(
            .bf_alias_key = .bf_sinas_native_name_key(.data[[alias_col]]),
            .bf_accepted_key = .bf_sinas_native_name_key(.data[[accepted_col]])
          ) |>
          dplyr::filter(
            !is.na(.data$.bf_alias_key), nzchar(.data$.bf_alias_key),
            !is.na(.data$.bf_accepted_key), nzchar(.data$.bf_accepted_key)
          ) |>
          dplyr::distinct(.data$.bf_alias_key, .data$.bf_accepted_key)

        alias_requested <- requested |>
          dplyr::inner_join(alias_tbl, by = c(".bf_requested_key" = ".bf_alias_key"))

        if (nrow(alias_requested)) {
          alias_matched <- alias_requested |>
            dplyr::left_join(
              sinas_tbl,
              by = c(".bf_accepted_key" = ".bf_sinas_taxon_key")
            )

          matched <- dplyr::bind_rows(matched, alias_matched)
        }
      }
    }
  }

  # Defensive exact-name fallback.
  #
  # SInAS 3.1.1 taxon names are already canonical enough for exact comparison.
  # The normal key-based join above should work, but exact lower-case matching is
  # kept as a fallback because user environments may load different taxonomy
  # cleaning helpers before this file, changing .bf_sinas_native_name_key()
  # behaviour. This fallback prevents valid SInAS rows from being lost merely
  # because the normalised key function was too aggressive or inconsistent.
  if (!any(!is.na(matched$.bf_sinas_taxon))) {
    requested_exact <- requested |>
      dplyr::mutate(
        .bf_requested_exact = tolower(trimws(as.character(.data$species)))
      ) |>
      dplyr::filter(!is.na(.data$.bf_requested_exact), nzchar(.data$.bf_requested_exact))

    sinas_exact <- sinas_tbl |>
      dplyr::mutate(
        .bf_sinas_exact = tolower(trimws(as.character(.data$.bf_sinas_taxon)))
      ) |>
      dplyr::filter(!is.na(.data$.bf_sinas_exact), nzchar(.data$.bf_sinas_exact))

    matched_exact <- requested_exact |>
      dplyr::left_join(
        sinas_exact,
        by = c(".bf_requested_exact" = ".bf_sinas_exact")
      )

    matched <- dplyr::bind_rows(matched, matched_exact)
  }

  matched <- matched |>
    dplyr::filter(!is.na(.data$.bf_sinas_taxon))

  if (!nrow(matched)) {
    long <- .bf_sinas_native_empty_long()
  } else {
    raw_area <- dplyr::coalesce(matched$raw_location, matched$.bf_sinas_location)

    long <- tibble::tibble(
      species = matched$species,
      source = "SInAS",
      accepted_name = matched$.bf_sinas_taxon,
      source_taxon_id = matched$.bf_sinas_taxon_id,
      raw_native_area = raw_area,
      raw_status = matched$.bf_sinas_status,
      origin_iso3 = toupper(matched$iso3c),
      evidence_type = "native_distribution",
      source_url = paste0("https://zenodo.org/records/", record_id)
    ) |>
      dplyr::distinct()
  }

  species_out <- .bf_sinas_native_collapse_species(long)
  unmapped <- long[is.na(long$origin_iso3) | !nzchar(long$origin_iso3), , drop = FALSE]

  summary <- tibble::tibble(
    n_species_input = length(species_vec),
    n_species_with_sinas_records = length(unique(long$species)),
    n_species_with_iso3 = sum(species_out$native_has_origin %in% TRUE),
    n_long_records = nrow(long),
    n_unmapped_records = nrow(unmapped),
    sources = "SInAS",
    sources_used = "SInAS"
  )

  out <- list(
    long = tibble::as_tibble(long),
    species = tibble::as_tibble(species_out),
    unmapped = tibble::as_tibble(unmapped),
    summary = summary
  )

  if (return == "species") return(out$species)
  if (return == "long") return(out$long)

  out
}
