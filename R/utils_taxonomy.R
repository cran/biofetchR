###############################################################################
# utils_taxonomy.R
# -----------------------------------------------------------------------------
# biofetchR taxonomy utilities
# -----------------------------------------------------------------------------
# Package-safe taxonomy helpers for preparing, cleaning and resolving taxon
# names used by the biofetchR GBIF occurrence pipelines.
#
# Design principles:
#   * keep raw names, cleaned names and resolved names auditable;
#   * filter obvious non-taxonomic strings before remote API calls;
#   * use GBIF as the primary resolver, with optional WoRMS support for marine
#     taxa or missing rank fields;
#   * keep cache paths package-safe and user-local;
#   * return tidy tibbles that can be joined directly to pipeline input data.
#
# Main exported helpers:
#   bf_clean_taxon_names()
#   bf_prepare_taxa_for_gbif()
#   bf_resolve_gbif_taxonomy()
#   bf_resolve_gbif_taxonomy_batch()
#   bf_attach_gbif_taxonomy()
#   bf_check_taxonomic_resolution()
#   bf_taxonomic_summary()
#   bf_apply_manual_taxonomy_fixes()
#
# Backward-compatible wrappers retained:
#   check_entrez_key()
#   resolve_species_names()
#   is_marine_species()
#
# Internal helpers use the `.bf_tax_*()` prefix and are documented with
# Roxygen2 `@noRd` blocks so their intent remains clear without generating
# exported help pages.
#
###############################################################################

# -----------------------------------------------------------------------------
# Small internal utilities
# -----------------------------------------------------------------------------

#' Squish taxonomy strings and standardise common encoding artefacts
#'
#' Converts input to character, replaces non-breaking spaces, normalises common
#' mojibake and Unicode dash variants, collapses repeated whitespace and trims
#' leading/trailing spaces.
#'
#' @param x Vector coercible to character.
#' @return Character vector with normalised spacing and dash-like characters.
#' @keywords internal
#' @noRd
.bf_tax_squish <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00a0", " ", x, fixed = TRUE)

  # Common mojibake / dash variants from copied PDFs, spreadsheets and Windows encodings.
  x <- gsub(
    "\u00e2[\u0080\u20ac][\u0093\u0094\u0098\u0099\u009c\u009d\u02dc\u0153\u2122\"-]",
    "-",
    x,
    perl = TRUE
  )

  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x)

  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

#' Remove simple authority strings from taxon names
#'
#' Removes bracketed comments and common trailing author/year fragments before
#' taxonomy matching. This is intentionally conservative and does not attempt
#' full nomenclatural parsing.
#'
#' @param x Character vector of raw taxon names.
#' @return Character vector with simple authority fragments removed.
#' @keywords internal
#' @noRd
.bf_tax_strip_authority <- function(x) {
  x <- .bf_tax_squish(x)

  # Remove bracketed authorities/comments.
  x <- gsub("\\s*\\([^)]*\\)", "", x)

  # Remove trailing author/year fragments after a comma.
  x <- gsub("\\s*,\\s*[A-Z][A-Za-z .'-]+\\s*,?\\s*[0-9]{3,4}\\s*$", "", x)
  x <- gsub("\\s*,\\s*[0-9]{3,4}\\s*$", "", x)

  .bf_tax_squish(x)
}

#' Detect obvious mojibake in taxonomy strings
#'
#' @param x Character vector.
#' @return Logical vector indicating likely broken text encoding.
#' @keywords internal
#' @noRd
.bf_tax_has_bad_encoding <- function(x) {
  x <- as.character(x)
  grepl("\\u00e2|\\ufffd|\\u20ac", x)
}

#' Detect measurement-like strings that are unlikely to be taxa
#'
#' Flags common environmental variables, measurement terms and simple chemical
#' ratio expressions that should be excluded before taxonomy resolution.
#'
#' @param x Character vector.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.bf_tax_has_measurement_pattern <- function(x) {
  x0 <- tolower(.bf_tax_squish(x))

  measurement_terms <- paste(
    c(
      "temperature", "conductivity", "salinity", "turbidity", "ph",
      "oxygen", "nitrogen", "phosphorus", "carbon", "ratio",
      "concentration", "biomass", "density", "cover", "abundance",
      "richness", "diversity", "area", "volume", "depth", "height",
      "weight", "mass", "length", "width", "diameter"
    ),
    collapse = "|"
  )

  has_term <- grepl(paste0("\\b(", measurement_terms, ")\\b"), x0)

  # # Examples: C-N ratio, C:N ratio, N/P ratio, C-N ratio.
  has_ratio_formula <- grepl(
    "^[a-z]{1,3}\\s*[-:/]+\\s*[a-z]{1,3}\\s+ratio$",
    x0
  )

  has_term | has_ratio_formula
}

#' Detect open nomenclature markers
#'
#' Flags strings containing markers such as `sp.`, `spp.`, `cf.`, `aff.`,
#' `complex`, `group`, `aggregate` or similar unresolved-name indicators.
#'
#' @param x Character vector.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.bf_tax_has_open_nomenclature <- function(x) {
  x0 <- tolower(.bf_tax_squish(x))

  grepl(
    "\\b(sp|spp|cf|aff|nr|indet|undet|complex|group|aggregate|agg)\\.?\\b",
    x0
  )
}

#' Detect abbreviated epithets within taxon names
#'
#' Flags patterns such as `Oncorhynchus c. henshawl`, which are ambiguous for
#' automated backbone resolution.
#'
#' @param x Character vector.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.bf_tax_has_abbreviated_epithet <- function(x) {
  x0 <- .bf_tax_squish(x)

  # Examples: Oncorhynchus c. henshawl, Genus s. str.
  grepl("^[A-Z][A-Za-z.-]+\\s+[a-z]\\.?\\s+[a-z][A-Za-z.-]+", x0)
}

#' Convert blank strings to missing values
#'
#' @param x Vector coercible to character.
#' @return Character vector with trimmed blanks converted to `NA_character_`.
#' @keywords internal
#' @noRd
.bf_tax_clean_blank <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

#' Transliterate taxonomy strings to ASCII
#'
#' Uses `stringi::stri_trans_general()` when available and falls back to
#' `iconv()`. This is used for stable cache and join keys, not for displayed
#' taxon names.
#'
#' @param x Character vector or `NULL`.
#' @return ASCII-transliterated character vector, or `NULL` if `x` is `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_to_ascii <- function(x) {
  if (is.null(x)) return(x)
  x <- as.character(x)
  if (requireNamespace("stringi", quietly = TRUE)) {
    return(stringi::stri_trans_general(x, "Latin-ASCII"))
  }
  iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
}


#' Create normalised taxonomy join/cache keys
#'
#' Converts strings to ASCII, lower case and underscore-separated tokens.
#'
#' @param x Character vector.
#' @return Character vector of stable normalised keys.
#' @keywords internal
#' @noRd
.bf_tax_norm_key <- function(x) {
  x <- .bf_tax_to_ascii(x)
  x <- tolower(.bf_tax_squish(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

#' Run a taxonomy query with simple retries
#'
#' Evaluates a zero-argument function, retries after errors and optionally
#' reports the final compact error. Failed attempts return `NULL` rather than
#' stopping the whole pipeline.
#'
#' @param expr_fun Zero-argument function that performs the remote call.
#' @param retries Number of attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return The successful result, or `NULL` if all attempts fail.
#' @keywords internal
#' @noRd
.bf_tax_with_retries <- function(expr_fun,
                                 retries = 3,
                                 sleep_sec = 0.15,
                                 quiet = TRUE) {
  last_error <- NULL

  for (i in seq_len(max(1L, as.integer(retries)))) {
    out <- tryCatch(
      expr_fun(),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )

    if (!is.null(out)) return(out)

    if (i < retries && is.finite(sleep_sec) && sleep_sec > 0) {
      Sys.sleep(sleep_sec * i)
    }
  }

  if (!isTRUE(quiet) && !is.null(last_error)) {
    message("Taxonomy lookup failed after retries: ", bf_compact_error(last_error))
  }

  NULL
}

#' Pick the first matching column name
#'
#' Performs case-insensitive, underscore-insensitive matching against a set of
#' candidate column names.
#'
#' @param df Data frame to inspect.
#' @param choices Character vector of acceptable column names.
#' @return Matching column name, or `NULL` if none is found.
#' @keywords internal
#' @noRd
.bf_tax_pick_col <- function(df, choices) {
  if (is.null(df) || !is.data.frame(df) || !length(choices)) return(NULL)
  n <- names(df)
  nn <- gsub("_", "", tolower(n))
  cc <- gsub("_", "", tolower(choices))
  hit <- match(cc, nn)
  hit <- hit[!is.na(hit)][1]
  if (length(hit) == 0L || is.na(hit)) return(NULL)
  n[[hit]]
}

#' Extract the first non-missing cleaned value
#'
#' @param x Vector coercible to character.
#' @return First non-missing value, or `NA_character_`.
#' @keywords internal
#' @noRd
.bf_tax_first_nonmissing <- function(x) {
  x <- .bf_tax_clean_blank(x)
  hit <- which(!is.na(x))[1]
  if (is.na(hit)) NA_character_ else x[[hit]]
}

#' Create an empty taxonomy-resolution record
#'
#' Returns a zero-row template or a one-row unresolved record with the same
#' schema used by GBIF/WoRMS resolver outputs.
#'
#' @param input_name Original input name. Use length-zero input to request a
#'   zero-row template.
#' @param clean_name Cleaned name associated with the input.
#' @param resolved_by Resolver label to record.
#' @param warning Optional warning/status message.
#' @return Tibble with the standard taxonomy-resolution columns.
#' @keywords internal
#' @noRd
.bf_tax_empty_record <- function(input_name = NA_character_,
                                 clean_name = NA_character_,
                                 resolved_by = NA_character_,
                                 warning = NA_character_) {
  if (length(input_name) == 0L) {
    return(tibble::tibble(
      input_name = character(),
      clean_name = character(),
      query_name = character(),
      usageKey = integer(),
      acceptedUsageKey = integer(),
      taxonKey_for_download = integer(),
      scientificName = character(),
      canonicalName = character(),
      rank = character(),
      status = character(),
      matchType = character(),
      confidence = numeric(),
      kingdom = character(),
      phylum = character(),
      class = character(),
      order = character(),
      family = character(),
      genus = character(),
      species = character(),
      resolved_by = character(),
      taxon_match_status = character(),
      taxon_match_warning = character(),
      is_genus_level = logical(),
      is_non_taxon = logical()
    ))
  }

  tibble::tibble(
    input_name = as.character(input_name),
    clean_name = as.character(clean_name),
    query_name = as.character(clean_name),
    usageKey = NA_integer_,
    acceptedUsageKey = NA_integer_,
    taxonKey_for_download = NA_integer_,
    scientificName = NA_character_,
    canonicalName = NA_character_,
    rank = NA_character_,
    status = NA_character_,
    matchType = NA_character_,
    confidence = NA_real_,
    kingdom = NA_character_,
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    species = NA_character_,
    resolved_by = as.character(resolved_by),
    taxon_match_status = "unresolved",
    taxon_match_warning = as.character(warning),
    is_genus_level = FALSE,
    is_non_taxon = FALSE
  )
}

# -----------------------------------------------------------------------------
# Name cleaning and diagnostics
# -----------------------------------------------------------------------------

#' Detect likely non-taxonomic strings
#'
#' Flags common environmental variables, habitat descriptors, measurement names,
#' and other strings that should not be submitted to taxonomy APIs.
#'
#' @param x Character vector.
#'
#' @return A logical vector with the same length as `x`. Values are `TRUE` for
#'   missing values, placeholders, measurement-like strings, environmental
#'   variables or other inputs that are unlikely to be biological taxon names.
#'   Values are `FALSE` for strings that remain plausible taxonomic names after
#'   the conservative screening rules.
#'
#' Detect strings that are unlikely to be biological taxon names
#' @export
bf_tax_looks_non_taxon <- function(x) {
  x_raw <- as.character(x)
  x0 <- tolower(.bf_tax_squish(x_raw))

  out <- is.na(x_raw) | !nzchar(x0)

  out <- out |
    x0 %in% c(
      "na", "n/a", "none", "null", "unknown", "unidentified",
      "not identified", "not recorded", "missing", "blank"
    )

  out <- out |
    grepl("^(unknown|unidentified|indeterminate)(\\s|$)", x0)

  out <- out |
    .bf_tax_has_measurement_pattern(x_raw)

  # Catch obvious data/variable fields rather than taxa.
  out <- out |
    grepl(
      "\\b(sample|plot|site|station|transect|quadrant|replicate|treatment|control|survey|record|observation|environment|substrate|sediment|soil|water)\\b",
      x0
    ) & !grepl("^[A-Z][a-z]+\\s+[a-z-]+", .bf_tax_squish(x_raw))

  # Strings with broken encoding plus measurement-like content should be removed.
  out <- out |
    (.bf_tax_has_bad_encoding(x_raw) & .bf_tax_has_measurement_pattern(x_raw))

  out[is.na(out)] <- TRUE
  out
}

#' Detect genus-level or open-nomenclature names
#'
#' @param x Character vector.
#'
#' @return A logical vector with the same length as `x`. Values are `TRUE` for
#'   names that appear to be genus-level, higher-level or open-nomenclature
#'   strings, including single-token genus names, `sp.`/`spp.` names,
#'   `cf.`/`aff.` names, groups, complexes or abbreviated epithets. Values are
#'   `FALSE` for names that are not flagged by these rules.
#'
#' Flag genus-level, higher-level, or open-name taxonomic strings
#' @export
bf_tax_is_genus_level <- function(x) {
  x0 <- bf_clean_taxon_names(x)

  out <- is.na(x0)

  # Single-token genus names.
  out <- out | grepl("^[A-Z][A-Za-z.-]+$", x0)

  # Open nomenclature: sp., spp., cf., aff., nr., group, complex, aggregate.
  out <- out | .bf_tax_has_open_nomenclature(x)

  # Abbreviated species/subspecies epithet, e.g. Oncorhynchus c. henshawl.
  out <- out | .bf_tax_has_abbreviated_epithet(x0)

  out[is.na(out)] <- FALSE
  out
}

#' Extract the genus component from a taxon name
#'
#' @param x Character vector.
#'
#' @return A character vector with the same length as `x`. Each element contains
#'   the extracted genus name from the corresponding cleaned taxon string, or
#'   `NA_character_` when no plausible genus component can be identified. The
#'   output is used for genus-level diagnostics and fallback taxonomy matching.
#'
#' @export
bf_tax_extract_genus <- function(x) {
  x0 <- bf_clean_taxon_names(x)

  out <- rep(NA_character_, length(x0))
  ok <- !is.na(x0) & grepl("^[A-Z][A-Za-z.-]+", x0)

  out[ok] <- sub("^([A-Z][A-Za-z.-]+).*", "\\1", x0[ok])
  out
}

#' Clean taxon names before taxonomy resolution
#'
#' Applies conservative string cleaning before GBIF or WoRMS lookup. The helper
#' removes simple authorship fragments, parenthetical notes, selected descriptor
#' words and implausible non-taxonomic strings while retaining plausible genus,
#' binomial and trinomial names for downstream auditing.
#'
#' @param x Character vector of raw taxon names.
#' @param max_tokens Integer maximum number of whitespace-separated tokens to
#'   retain after cleaning. The default, `3`, preserves trinomial names.
#' @param keep_genus_for_open_names Logical. If `TRUE`, strings such as
#'   `Genus sp.` or `Genus spp.` are converted to `Genus`; if `FALSE`, they are
#'   converted to `NA_character_`.
#'
#' @return A character vector with the same length as `x`. Each element contains
#'   a cleaned taxon-name candidate after whitespace, authority, punctuation,
#'   descriptor and open-name handling. Blank, non-taxonomic or implausible
#'   entries are returned as `NA_character_`. The returned names are intended for
#'   GBIF or WoRMS taxonomy resolution and for audit columns in downstream
#'   pipelines.
#'
#' @family taxonomy cleaning helpers
#' @export
bf_clean_taxon_names <- function(x,
                                 max_tokens = 3,
                                 keep_genus_for_open_names = TRUE) {
  x_raw <- as.character(x)
  x1 <- .bf_tax_strip_authority(x_raw)

  # Remove clear non-taxa early.
  non_taxon <- bf_tax_looks_non_taxon(x1)
  x1[non_taxon] <- NA_character_

  # Standardise whitespace and punctuation.
  x1 <- .bf_tax_squish(x1)
  x1 <- gsub("[;:,]+$", "", x1)
  x1 <- .bf_tax_squish(x1)

  # Remove common qualifier words that are not part of taxon names.
  x1 <- gsub("\\b(domestic|feral|wild|cultivated|introduced|invasive|alien|native)\\b", "", x1, ignore.case = TRUE)
  x1 <- .bf_tax_squish(x1)

  # Convert open names according to the requested policy. Keeping the genus
  # allows genus-level auditing/fallbacks; setting this to FALSE removes them
  # from downstream species-level candidate sets.
  open_genus <- grepl("^([A-Z][A-Za-z.-]+)\\s+spp?\\.?$", x1, ignore.case = FALSE)
  open_genus[is.na(open_genus)] <- FALSE

  if (isTRUE(keep_genus_for_open_names)) {
    x1 <- gsub("^([A-Z][A-Za-z.-]+)\\s+spp?\\.?$", "\\1", x1, ignore.case = FALSE)
  } else {
    x1[open_genus] <- NA_character_
  }

  # Preserve cf./aff./nr. forms for review rather than collapsing to genus.
  x1 <- gsub("\\b(cf|aff|nr)\\s+", "\\1. ", x1, ignore.case = TRUE)
  x1 <- .bf_tax_squish(x1)

  # Limit very long names/comments to the requested number of taxonomic
  # tokens. The default of three preserves trinomial names while removing
  # common trailing notes that survived earlier cleaning.
  if (!is.null(max_tokens) && is.finite(max_tokens) && max_tokens > 0) {
    max_tokens <- as.integer(max_tokens)
    x1 <- vapply(
      strsplit(x1, "[[:space:]]+"),
      function(z) {
        if (length(z) == 1L && is.na(z)) return(NA_character_)
        paste(utils::head(z, max_tokens), collapse = " ")
      },
      character(1)
    )
  }

  # Keep only plausible taxonomic strings.
  plausible <- !is.na(x1) & (
    grepl("^[A-Z][A-Za-z.-]+$", x1) |
      grepl("^[A-Z][A-Za-z.-]+\\s+[a-z][A-Za-z.-]+", x1) |
      grepl("^[A-Z][A-Za-z.-]+\\s+(cf|aff|nr)\\.\\s+[a-z][A-Za-z.-]+", x1, ignore.case = TRUE)
  )

  x1[!plausible] <- NA_character_
  x1[!nzchar(x1)] <- NA_character_

  x1
}

# -----------------------------------------------------------------------------
# GBIF and optional resolver helpers
# -----------------------------------------------------------------------------

#' Flatten a GBIF backbone result into biofetchR taxonomy columns
#'
#' Converts the list returned by `rgbif::name_backbone()` into the standard
#' one-row tibble schema used throughout this file.
#'
#' @param res GBIF backbone result list, or `NULL`.
#' @param input_name Original input name.
#' @param clean_name Cleaned input name.
#' @param query_name Name actually submitted to GBIF.
#' @param resolved_by Resolver label to record.
#' @return One-row taxonomy-resolution tibble.
#' @keywords internal
#' @noRd
.bf_tax_flatten_gbif <- function(res,
                                 input_name = NA_character_,
                                 clean_name = NA_character_,
                                 query_name = clean_name,
                                 resolved_by = "gbif_backbone") {
  if (is.null(res) || !length(res)) {
    return(.bf_tax_empty_record(
      input_name = input_name,
      clean_name = clean_name,
      resolved_by = resolved_by,
      warning = "No GBIF match returned."
    ))
  }

  usage_key <- suppressWarnings(as.integer(res$usageKey %||% NA_integer_))
  accepted_key <- suppressWarnings(as.integer(res$acceptedUsageKey %||% NA_integer_))
  key_for_download <- if (!is.na(accepted_key)) accepted_key else usage_key

  status <- as.character(res$status %||% res$taxonomicStatus %||% NA_character_)
  match_type <- as.character(res$matchType %||% NA_character_)
  confidence <- suppressWarnings(as.numeric(res$confidence %||% NA_real_))

  match_status <- if (!is.na(usage_key)) {
    if (!is.na(match_type) && match_type %in% c("EXACT", "HIGHERRANK")) "matched" else "matched_non_exact"
  } else {
    "unresolved"
  }

  warning <- NA_character_
  if (identical(match_status, "matched_non_exact")) {
    warning <- paste0("GBIF matchType is ", match_type, ". Check before using in strict analyses.")
  }

  tibble::tibble(
    input_name = as.character(input_name),
    clean_name = as.character(clean_name),
    query_name = as.character(query_name),
    usageKey = usage_key,
    acceptedUsageKey = accepted_key,
    taxonKey_for_download = key_for_download,
    scientificName = as.character(res$scientificName %||% NA_character_),
    canonicalName = as.character(res$canonicalName %||% NA_character_),
    rank = as.character(res$rank %||% NA_character_),
    status = status,
    matchType = match_type,
    confidence = confidence,
    kingdom = as.character(res$kingdom %||% NA_character_),
    phylum = as.character(res$phylum %||% NA_character_),
    class = as.character(res$class %||% NA_character_),
    order = as.character(res$order %||% NA_character_),
    family = as.character(res$family %||% NA_character_),
    genus = as.character(res$genus %||% NA_character_),
    species = as.character(res$species %||% NA_character_),
    resolved_by = as.character(resolved_by),
    taxon_match_status = match_status,
    taxon_match_warning = warning,
    is_genus_level = identical(as.character(res$rank %||% NA_character_), "GENUS"),
    is_non_taxon = FALSE
  )
}

#' Query the GBIF backbone with retries
#'
#' Thin wrapper around `rgbif::name_backbone()` that returns `NULL` when
#' `rgbif` is unavailable, the name is blank or all retry attempts fail.
#'
#' @param name Character taxon name.
#' @param strict Passed to `rgbif::name_backbone()`.
#' @param rank Optional GBIF rank constraint.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return GBIF backbone result list, or `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_gbif_backbone <- function(name,
                                  strict = FALSE,
                                  rank = NULL,
                                  retries = 3,
                                  sleep_sec = 0.15,
                                  quiet = TRUE) {
  if (is.na(name) || !nzchar(name)) return(NULL)
  if (!requireNamespace("rgbif", quietly = TRUE)) return(NULL)

  .bf_tax_with_retries(
    function() {
      args <- list(name = name, strict = strict, verbose = FALSE)
      if (!is.null(rank)) args$rank <- rank
      do.call(rgbif::name_backbone, args)
    },
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )
}

#' Parse GBIF canonical names in batches
#'
#' Uses `rgbif::name_parse()` to extract canonical names while preserving input
#' order and returning `NA` for failed chunks.
#'
#' @param x Character vector of taxon names.
#' @param batch_size Number of names per GBIF parse request.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return Character vector of canonical names.
#' @keywords internal
#' @noRd
.bf_tax_gbif_name_parse_canonical <- function(x,
                                              batch_size = 500,
                                              retries = 2,
                                              sleep_sec = 0.1,
                                              quiet = TRUE) {
  out <- rep(NA_character_, length(x))
  if (!requireNamespace("rgbif", quietly = TRUE)) return(out)
  if (!length(x)) return(out)

  idx <- seq_along(x)
  batches <- split(idx, ceiling(seq_along(idx) / max(1L, as.integer(batch_size))))

  for (b in batches) {
    chunk <- x[b]
    res <- .bf_tax_with_retries(
      function() rgbif::name_parse(scientificname = chunk),
      retries = retries,
      sleep_sec = sleep_sec,
      quiet = quiet
    )

    if (is.data.frame(res) && nrow(res) == length(chunk)) {
      canon <- NULL
      if ("canonicalName" %in% names(res)) canon <- res$canonicalName
      if (is.null(canon) && "canonical" %in% names(res)) canon <- res$canonical
      if (is.null(canon) && "scientificName" %in% names(res)) canon <- res$scientificName

      if (!is.null(canon)) {
        canon <- .bf_tax_clean_blank(.bf_tax_squish(canon))
        out[b] <- canon
      }
    }
  }

  out
}

#' Choose the best GBIF name suggestion
#'
#' Queries `rgbif::name_suggest()` and prioritises accepted species/subspecies
#' matches before genus-level suggestions.
#'
#' @param q Query string.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return Suggested canonical/scientific name, or `NA_character_`.
#' @keywords internal
#' @noRd
.bf_tax_gbif_suggest_best <- function(q,
                                      retries = 2,
                                      sleep_sec = 0.15,
                                      quiet = TRUE) {
  if (is.na(q) || !nzchar(q)) return(NA_character_)
  if (!requireNamespace("rgbif", quietly = TRUE)) return(NA_character_)

  res <- .bf_tax_with_retries(
    function() rgbif::name_suggest(q = q, limit = 10),
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )

  dat <- res$data %||% NULL
  if (is.null(dat) || !is.data.frame(dat) || !nrow(dat)) return(NA_character_)

  status <- if ("status" %in% names(dat)) as.character(dat$status) else NA_character_
  rank <- if ("rank" %in% names(dat)) as.character(dat$rank) else NA_character_
  nm <- if ("canonicalName" %in% names(dat)) {
    as.character(dat$canonicalName)
  } else if ("scientificName" %in% names(dat)) {
    as.character(dat$scientificName)
  } else {
    rep(NA_character_, nrow(dat))
  }

  score <- rep(0L, nrow(dat))
  score[status %in% c("ACCEPTED", "ACCEPTED_NAME")] <- score[status %in% c("ACCEPTED", "ACCEPTED_NAME")] + 100L
  score[rank %in% c("SPECIES", "SUBSPECIES", "VARIETY")] <- score[rank %in% c("SPECIES", "SUBSPECIES", "VARIETY")] + 10L
  score[rank %in% c("GENUS")] <- score[rank %in% c("GENUS")] + 1L

  ord <- order(score, decreasing = TRUE, na.last = TRUE)
  nm <- .bf_tax_clean_blank(.bf_tax_squish(nm[ord]))
  nm <- nm[!is.na(nm)]
  if (!length(nm)) NA_character_ else nm[[1]]
}

#' Retrieve a GBIF name-usage record
#'
#' @param key GBIF usage key.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return GBIF usage list, or `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_gbif_usage <- function(key,
                               retries = 3,
                               sleep_sec = 0.15,
                               quiet = TRUE) {
  if (is.na(key) || !is.finite(key)) return(NULL)
  if (!requireNamespace("rgbif", quietly = TRUE)) return(NULL)
  key <- as.integer(key)

  .bf_tax_with_retries(
    function() rgbif::name_usage(key = key),
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )
}

#' Retrieve a GBIF classification table
#'
#' @param key GBIF usage key.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return GBIF classification result, or `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_gbif_classification <- function(key,
                                        retries = 3,
                                        sleep_sec = 0.15,
                                        quiet = TRUE) {
  if (is.na(key) || !is.finite(key)) return(NULL)
  if (!requireNamespace("rgbif", quietly = TRUE)) return(NULL)
  key <- as.integer(key)

  .bf_tax_with_retries(
    function() rgbif::name_usage(key = key, data = "classification"),
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )
}

#' Extract standard ranks from a classification table
#'
#' Reads a GBIF/WoRMS-like classification table and returns one row containing
#' kingdom, phylum, class, order, family, genus and species.
#'
#' @param class_df Classification data frame or list containing one.
#' @param name_choices Candidate columns containing taxon names.
#' @return One-row tibble of standard rank columns.
#' @keywords internal
#' @noRd
.bf_tax_extract_ranks_from_classification <- function(class_df,
                                                      name_choices = c("name", "scientificName", "scientificname", "canonicalName", "canonicalname")) {
  empty <- tibble::tibble(
    kingdom = NA_character_,
    phylum = NA_character_,
    class = NA_character_,
    order = NA_character_,
    family = NA_character_,
    genus = NA_character_,
    species = NA_character_
  )

  if (is.null(class_df)) return(empty)
  if (is.list(class_df) && length(class_df) == 1L && is.data.frame(class_df[[1]])) class_df <- class_df[[1]]
  if (!is.data.frame(class_df) || !nrow(class_df)) return(empty)

  rank_col <- .bf_tax_pick_col(class_df, c("rank", "taxonRank", "taxonomicRank"))
  name_col <- .bf_tax_pick_col(class_df, name_choices)
  if (is.null(rank_col) || is.null(name_col)) return(empty)

  rk <- tolower(.bf_tax_squish(class_df[[rank_col]]))
  nm <- .bf_tax_clean_blank(.bf_tax_squish(class_df[[name_col]]))

  pick <- function(r) {
    hit <- which(rk == r & !is.na(nm))
    if (!length(hit)) return(NA_character_)
    nm[tail(hit, 1)]
  }

  tibble::tibble(
    kingdom = pick("kingdom"),
    phylum = pick("phylum"),
    class = pick("class"),
    order = pick("order"),
    family = pick("family"),
    genus = pick("genus"),
    species = pick("species")
  )
}

#' Fill missing rank fields from a GBIF key
#'
#' Uses GBIF usage and classification endpoints to complete missing standard rank
#' fields in a one-row taxonomy result.
#'
#' @param row One-row taxonomy-resolution tibble.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return Updated one-row taxonomy-resolution tibble.
#' @keywords internal
#' @noRd
.bf_tax_complete_ranks_by_key <- function(row,
                                          retries = 3,
                                          sleep_sec = 0.15,
                                          quiet = TRUE) {
  rank_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  key <- row$taxonKey_for_download[[1]]

  if (is.na(key) || !is.finite(key)) return(row)

  missing_ranks <- rank_cols[rank_cols %in% names(row)][vapply(row[rank_cols], function(z) is.na(z[[1]]) || !nzchar(as.character(z[[1]])), logical(1))]
  if (!length(missing_ranks)) return(row)

  usage <- .bf_tax_gbif_usage(key, retries = retries, sleep_sec = sleep_sec, quiet = quiet)
  if (!is.null(usage)) {
    for (rk in rank_cols) {
      val <- usage[[rk]] %||% NA_character_
      if (rk %in% names(row) && (is.na(row[[rk]][[1]]) || !nzchar(as.character(row[[rk]][[1]]))) && !is.na(val)) {
        row[[rk]] <- as.character(val)
      }
    }
  }

  missing_ranks <- rank_cols[rank_cols %in% names(row)][vapply(row[rank_cols], function(z) is.na(z[[1]]) || !nzchar(as.character(z[[1]])), logical(1))]
  if (!length(missing_ranks)) return(row)

  cls <- .bf_tax_gbif_classification(key, retries = retries, sleep_sec = sleep_sec, quiet = quiet)
  dat <- cls$data %||% cls$classification %||% cls %||% NULL
  ranks <- .bf_tax_extract_ranks_from_classification(dat)

  for (rk in rank_cols) {
    val <- ranks[[rk]][[1]]
    if (rk %in% names(row) && (is.na(row[[rk]][[1]]) || !nzchar(as.character(row[[rk]][[1]]))) && !is.na(val)) {
      row[[rk]] <- val
    }
  }

  row
}

#' Query WoRMS records for a taxon name
#'
#' Uses `{worrms}` when installed and tries both broad and marine-only searches.
#'
#' @param name Character taxon name.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return WoRMS records data frame, or `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_worms_record <- function(name,
                                 retries = 2,
                                 sleep_sec = 0.2,
                                 quiet = TRUE) {
  if (!requireNamespace("worrms", quietly = TRUE)) return(NULL)
  if (is.na(name) || !nzchar(name)) return(NULL)

  tries <- list(
    list(name = name, marine_only = FALSE),
    list(name = name, marine_only = TRUE),
    list(name = name)
  )

  for (args in tries) {
    res <- .bf_tax_with_retries(
      function() do.call(worrms::wm_records_name, args),
      retries = retries,
      sleep_sec = sleep_sec,
      quiet = quiet
    )
    if (is.data.frame(res) && nrow(res)) return(res)
  }

  NULL
}

#' Select an AphiaID from WoRMS records
#'
#' Prioritises accepted records and valid AphiaID columns when available.
#'
#' @param records_df Data frame returned by `{worrms}`.
#' @return Integer AphiaID, or `NA_integer_`.
#' @keywords internal
#' @noRd
.bf_tax_worms_pick_aphia <- function(records_df) {
  if (is.null(records_df) || !is.data.frame(records_df) || !nrow(records_df)) return(NA_integer_)
  nms <- names(records_df)
  id_col <- nms[grepl("aphia", tolower(nms))][1]
  if (is.na(id_col)) return(NA_integer_)

  rec <- records_df
  status_col <- nms[tolower(nms) == "status"][1]
  valid_col <- nms[grepl("valid", tolower(nms)) & grepl("aphia", tolower(nms))][1]

  if (!is.na(status_col)) {
    accepted <- rec[tolower(as.character(rec[[status_col]])) == "accepted", , drop = FALSE]
    if (nrow(accepted)) rec <- accepted
  }

  if (!is.na(valid_col) && any(!is.na(rec[[valid_col]]))) {
    return(suppressWarnings(as.integer(rec[[valid_col]][which(!is.na(rec[[valid_col]]))[1]])))
  }

  suppressWarnings(as.integer(rec[[id_col]][1]))
}

#' Retrieve WoRMS classification for an AphiaID
#'
#' @param aphia_id WoRMS AphiaID.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return WoRMS classification object, or `NULL`.
#' @keywords internal
#' @noRd
.bf_tax_worms_classification <- function(aphia_id,
                                         retries = 2,
                                         sleep_sec = 0.2,
                                         quiet = TRUE) {
  if (!requireNamespace("worrms", quietly = TRUE)) return(NULL)
  if (is.na(aphia_id)) return(NULL)

  .bf_tax_with_retries(
    function() worrms::wm_classification(aphia_id),
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )
}

#' Resolve a taxon name with WoRMS as an external fallback
#'
#' Returns the standard taxonomy-resolution schema, but does not provide a GBIF
#' taxon key. Used only when marine/rank fallback is explicitly requested.
#'
#' @param name Name submitted to WoRMS.
#' @param input_name Original input name.
#' @param clean_name Cleaned input name.
#' @param retries Number of remote-call attempts.
#' @param sleep_sec Base sleep interval between failed attempts.
#' @param quiet Suppress retry-failure messages.
#' @return One-row taxonomy-resolution tibble.
#' @keywords internal
#' @noRd
.bf_tax_worms_resolve <- function(name,
                                  input_name = name,
                                  clean_name = name,
                                  retries = 2,
                                  sleep_sec = 0.2,
                                  quiet = TRUE) {
  rec <- .bf_tax_worms_record(name, retries = retries, sleep_sec = sleep_sec, quiet = quiet)
  aphia <- .bf_tax_worms_pick_aphia(rec)
  if (is.na(aphia)) {
    return(.bf_tax_empty_record(
      input_name = input_name,
      clean_name = clean_name,
      resolved_by = "worms",
      warning = "No WoRMS AphiaID match."
    ))
  }

  cls <- .bf_tax_worms_classification(aphia, retries = retries, sleep_sec = sleep_sec, quiet = quiet)
  ranks <- .bf_tax_extract_ranks_from_classification(cls, name_choices = c("scientificname", "scientificName", "name"))

  scientific <- NA_character_
  canonical <- NA_character_
  status <- NA_character_
  if (is.data.frame(rec) && nrow(rec)) {
    sci_col <- .bf_tax_pick_col(rec, c("scientificname", "scientificName", "valid_name"))
    status_col <- .bf_tax_pick_col(rec, c("status"))
    if (!is.null(sci_col)) scientific <- as.character(rec[[sci_col]][1])
    if (!is.null(status_col)) status <- as.character(rec[[status_col]][1])
    canonical <- scientific
  }

  tibble::tibble(
    input_name = as.character(input_name),
    clean_name = as.character(clean_name),
    query_name = as.character(name),
    usageKey = NA_integer_,
    acceptedUsageKey = NA_integer_,
    taxonKey_for_download = NA_integer_,
    scientificName = scientific,
    canonicalName = canonical,
    rank = NA_character_,
    status = status,
    matchType = "WORMS",
    confidence = NA_real_,
    kingdom = ranks$kingdom,
    phylum = ranks$phylum,
    class = ranks$class,
    order = ranks$order,
    family = ranks$family,
    genus = ranks$genus,
    species = ranks$species,
    resolved_by = "worms",
    taxon_match_status = if (any(!is.na(ranks))) "matched_external" else "unresolved",
    taxon_match_warning = if (any(!is.na(ranks))) NA_character_ else "WoRMS matched but returned no usable ranks.",
    is_genus_level = FALSE,
    is_non_taxon = FALSE
  )
}

# -----------------------------------------------------------------------------
# Manual override support
# -----------------------------------------------------------------------------

#' Standardise manual taxonomy override inputs
#'
#' Converts `NULL`, named vectors or data frames into a two-column override table
#' with `input_name` and `replacement_name`.
#'
#' @param manual_overrides `NULL`, named character vector or data frame.
#' @return Tibble with `input_name` and `replacement_name`.
#' @keywords internal
#' @noRd
.bf_tax_standardise_manual_overrides <- function(manual_overrides = NULL) {
  if (is.null(manual_overrides)) {
    return(tibble::tibble(input_name = character(), replacement_name = character()))
  }

  if (is.character(manual_overrides) && !is.null(names(manual_overrides))) {
    return(tibble::tibble(
      input_name = names(manual_overrides),
      replacement_name = unname(as.character(manual_overrides))
    ))
  }

  if (is.data.frame(manual_overrides)) {
    nms <- names(manual_overrides)
    key_col <- .bf_tax_pick_col(manual_overrides, c("input_name", "original_name", "native_species_name", "species", "name"))
    fix_col <- .bf_tax_pick_col(manual_overrides, c("replacement_name", "manual_fix", "resolved_name", "accepted_name", "taxon_name"))

    if (is.null(key_col) || is.null(fix_col)) {
      stop(
        "manual_overrides data frame must contain a key column ",
        "(input_name/original_name/native_species_name/species/name) and a fix column ",
        "(replacement_name/manual_fix/resolved_name/accepted_name/taxon_name).",
        call. = FALSE
      )
    }

    return(tibble::tibble(
      input_name = as.character(manual_overrides[[key_col]]),
      replacement_name = as.character(manual_overrides[[fix_col]])
    ))
  }

  stop("manual_overrides must be NULL, a named character vector, or a data frame.", call. = FALSE)
}

#' Apply manual taxonomy name and rank fixes to a data frame
#'
#' @param x Data frame.
#' @param manual_fixes Named character vector or data frame with manual fixes.
#' @param name_col Name column in `x`.
#' @param overwrite_names Replace names when a manual fix is available.
#' @param overwrite_ranks Overwrite existing rank columns when manual rank values are supplied.
#' @param rank_cols Rank columns to fill if they exist in both tables.
#'
#' @return A data frame with the same rows as `x` and with manual taxonomy fixes
#'   applied where matching keys are found. If `overwrite_names = TRUE`, values
#'   in `name_col` are replaced by matched manual replacement names. If
#'   `overwrite_ranks = TRUE`, matching rank columns are also updated from
#'   `manual_fixes`; otherwise rank values are filled only where the input value
#'   is missing or blank. If no usable manual fixes are supplied, the original
#'   input data frame is returned unchanged.
#'
#' @export
bf_apply_manual_taxonomy_fixes <- function(x,
                                           manual_fixes,
                                           name_col = "species",
                                           overwrite_names = TRUE,
                                           overwrite_ranks = FALSE,
                                           rank_cols = c("kingdom", "phylum", "class", "order", "family", "genus", "species_rank")) {
  bf_require_packages(c("dplyr", "tibble"), context = "taxonomy helper")

  if (!is.data.frame(x)) stop("`x` must be a data frame.", call. = FALSE)
  if (!name_col %in% names(x)) stop("Name column not found: ", name_col, call. = FALSE)

  mf <- .bf_tax_standardise_manual_overrides(manual_fixes)
  if (!nrow(mf)) return(x)

  mf <- mf[!is.na(mf$input_name) & nzchar(mf$input_name), , drop = FALSE]
  mf$.key <- .bf_tax_norm_key(mf$input_name)
  mf <- mf[!duplicated(mf$.key, fromLast = TRUE), , drop = FALSE]

  out <- x
  out$.bf_tax_key <- .bf_tax_norm_key(out[[name_col]])
  idx <- match(out$.bf_tax_key, mf$.key)

  if (isTRUE(overwrite_names)) {
    replace <- !is.na(idx) & !is.na(mf$replacement_name[idx]) & nzchar(mf$replacement_name[idx])
    out[[name_col]][replace] <- mf$replacement_name[idx[replace]]
  }

  rank_cols <- intersect(rank_cols, intersect(names(out), names(manual_fixes)))
  if (length(rank_cols) && is.data.frame(manual_fixes)) {
    for (rk in rank_cols) {
      manual_vals <- manual_fixes[[rk]][idx]
      manual_vals <- .bf_tax_clean_blank(manual_vals)
      if (isTRUE(overwrite_ranks)) {
        fill <- !is.na(manual_vals)
      } else {
        fill <- (is.na(out[[rk]]) | !nzchar(as.character(out[[rk]]))) & !is.na(manual_vals)
      }
      out[[rk]][fill] <- manual_vals[fill]
    }
  }

  out$.bf_tax_key <- NULL
  out
}

#' Prepare taxon names for GBIF download pipelines
#'
#' Applies manual fixes first, then cleans names, removes obvious non-taxa,
#' blocks genus-level/open/abbreviated names, and returns accepted, rejected,
#' audit and summary tables.
#'
#' @param df Input data frame.
#' @param name_col Column containing raw taxon names.
#' @param iso2_col Optional ISO2 country column. Use `"iso2c"` for terrestrial/freshwater workflows.
#' @param manual_fixes Optional named character vector or data frame passed to
#'   `bf_apply_manual_taxonomy_fixes()`.
#' @param require_species_level Logical. If TRUE, only species-level candidates are accepted.
#' @param deduplicate Logical. If TRUE, accepted names are deduplicated.
#'
#' @return A named list with four elements: `accepted`, a tibble/data frame of
#'   cleaned taxon records retained for GBIF download workflows; `rejected`, a
#'   tibble/data frame of records removed before GBIF submission with rejection
#'   diagnostics; `audit`, a tibble recording raw names, fixed names, cleaned
#'   names, genus values, acceptance flags and rejection reasons for every input
#'   row; and `summary`, a one-row tibble with counts and percentages describing
#'   the preparation outcome. The output separates records suitable for download
#'   from records rejected as non-taxonomic, open-name or non-species-level
#'   candidates.
#'
#' @export
bf_prepare_taxa_for_gbif <- function(df,
                                     name_col = "species",
                                     iso2_col = NULL,
                                     manual_fixes = NULL,
                                     require_species_level = TRUE,
                                     deduplicate = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  if (!name_col %in% names(df)) {
    stop("`name_col` was not found in `df`: ", name_col, call. = FALSE)
  }

  if (!is.null(iso2_col) && !iso2_col %in% names(df)) {
    stop("`iso2_col` was not found in `df`: ", iso2_col, call. = FALSE)
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required.", call. = FALSE)
  }

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required.", call. = FALSE)
  }

  raw_df <- tibble::as_tibble(df)
  raw_names <- as.character(raw_df[[name_col]])

  fixed_df <- raw_df

  if (!is.null(manual_fixes)) {
    fixed_df <- bf_apply_manual_taxonomy_fixes(
      fixed_df,
      manual_fixes = manual_fixes,
      name_col = name_col,
      overwrite_names = TRUE
    )
  }

  fixed_names <- as.character(fixed_df[[name_col]])
  cleaned_names <- bf_clean_taxon_names(fixed_names)

  is_non_taxon <- bf_tax_looks_non_taxon(fixed_names) | is.na(cleaned_names)

  is_open_or_higher <- rep(FALSE, length(cleaned_names))
  taxon_rows <- !is_non_taxon & !is.na(cleaned_names)

  if (any(taxon_rows)) {
    is_open_or_higher[taxon_rows] <- bf_tax_is_genus_level(cleaned_names[taxon_rows])
  }

  genus <- bf_tax_extract_genus(cleaned_names)

  looks_species_level <- !is_non_taxon &
    !is_open_or_higher &
    !is.na(cleaned_names) &
    grepl(
      "^[A-Z][A-Za-z.-]+[[:space:]]+[a-z][A-Za-z.-]+",
      cleaned_names
    )

  accepted_flag <- if (isTRUE(require_species_level)) {
    looks_species_level
  } else {
    !is_non_taxon & !is.na(cleaned_names)
  }

  rejection_reason <- dplyr::case_when(
    accepted_flag ~ NA_character_,
    is_non_taxon ~ "non_taxon_or_placeholder",
    is_open_or_higher ~ "open_name_or_higher_rank",
    is.na(cleaned_names) ~ "missing_cleaned_name",
    TRUE ~ "not_species_level_candidate"
  )

  audit <- tibble::tibble(
    raw_name = raw_names,
    fixed_name = fixed_names,
    cleaned_name = cleaned_names,
    was_manual_fixed = raw_names != fixed_names,
    is_non_taxon = is_non_taxon,
    is_open_or_higher = is_open_or_higher,
    genus = genus,
    species_level_candidate = looks_species_level,
    accepted = accepted_flag,
    rejection_reason = rejection_reason
  )

  if (!is.null(iso2_col)) {
    audit[[iso2_col]] <- as.character(raw_df[[iso2_col]])
  }

  accepted <- fixed_df[accepted_flag, , drop = FALSE]
  accepted[[name_col]] <- cleaned_names[accepted_flag]

  # Pipelines expect a column called `species`.
  accepted$species <- cleaned_names[accepted_flag]

  if (!is.null(iso2_col)) {
    accepted[[iso2_col]] <- as.character(accepted[[iso2_col]])

    # Terrestrial/freshwater pipeline expects `iso2c`.
    if (!"iso2c" %in% names(accepted)) {
      accepted$iso2c <- accepted[[iso2_col]]
    }
  }

  rejected <- raw_df[!accepted_flag, , drop = FALSE]
  rejected$raw_name <- raw_names[!accepted_flag]
  rejected$fixed_name <- fixed_names[!accepted_flag]
  rejected$cleaned_name <- cleaned_names[!accepted_flag]
  rejected$rejection_reason <- rejection_reason[!accepted_flag]
  rejected$is_non_taxon <- is_non_taxon[!accepted_flag]
  rejected$is_open_or_higher <- is_open_or_higher[!accepted_flag]
  rejected$species_level_candidate <- looks_species_level[!accepted_flag]

  if (isTRUE(deduplicate) && nrow(accepted)) {
    if (!is.null(iso2_col) && "iso2c" %in% names(accepted)) {
      accepted <- accepted |>
        dplyr::distinct(.data$species, .data$iso2c, .keep_all = TRUE)
    } else {
      accepted <- accepted |>
        dplyr::distinct(.data$species, .keep_all = TRUE)
    }
  }

  summary <- tibble::tibble(
    n_input = length(raw_names),
    n_accepted = sum(accepted_flag, na.rm = TRUE),
    n_rejected = sum(!accepted_flag, na.rm = TRUE),
    n_manual_fixed = sum(raw_names != fixed_names, na.rm = TRUE),
    n_non_taxon = sum(is_non_taxon, na.rm = TRUE),
    n_open_or_higher = sum(is_open_or_higher, na.rm = TRUE),
    n_species_level_candidates = sum(looks_species_level, na.rm = TRUE),
    pct_accepted = round(100 * sum(accepted_flag, na.rm = TRUE) / length(raw_names), 1)
  )

  list(
    accepted = accepted,
    rejected = rejected,
    audit = audit,
    summary = summary
  )
}

# -----------------------------------------------------------------------------
# Main taxonomy resolution API
# -----------------------------------------------------------------------------

#' Resolve one taxon name against GBIF
#'
#' @param name Taxon name.
#' @param strict Passed to `rgbif::name_backbone()`.
#' @param use_gbif_parse Use `rgbif::name_parse()` for canonicalisation.
#' @param use_gbif_suggest Try `rgbif::name_suggest()` if backbone fails.
#' @param use_genus_fallback If TRUE, retry a genus-level match for open names.
#' @param use_worms_fallback If TRUE and `{worrms}` is installed, use WoRMS as an
#'   external fallback for unresolved names/ranks. This does not provide a GBIF
#'   taxonKey but can fill rank fields.
#' @param manual_overrides Optional named vector or data frame of manual name fixes.
#' @param retries Number of retries for remote calls.
#' @param sleep_sec Delay between failed retries.
#' @param quiet Reduce messages.
#'
#' @return A one-row tibble describing the taxonomy-resolution outcome for
#'   `name`. The output contains the original and cleaned names, the query name,
#'   GBIF usage keys, the taxon key selected for downloads, resolved scientific
#'   and canonical names, rank, status, match type, confidence, standard
#'   classification columns from kingdom to species, resolver provenance,
#'   `taxon_match_status`, any warning text, and logical flags for genus-level
#'   or non-taxonomic inputs. Unresolved or skipped names are returned in the
#'   same schema with missing key fields and diagnostic status values.
#'
#' @export
bf_resolve_gbif_taxonomy <- function(name,
                                     strict = FALSE,
                                     use_gbif_parse = TRUE,
                                     use_gbif_suggest = TRUE,
                                     use_genus_fallback = TRUE,
                                     use_worms_fallback = FALSE,
                                     manual_overrides = NULL,
                                     retries = 3,
                                     sleep_sec = 0.15,
                                     quiet = TRUE) {
  bf_require_packages(c("tibble"), context = "taxonomy helper")

  input_name <- as.character(name)[1]
  if (is.na(input_name) || !nzchar(trimws(input_name))) {
    return(.bf_tax_empty_record(input_name, NA_character_, warning = "Input name is empty."))
  }

  is_non_taxon <- bf_tax_looks_non_taxon(input_name)
  if (isTRUE(is_non_taxon)) {
    out <- .bf_tax_empty_record(input_name, NA_character_, warning = "Input looks non-taxonomic; skipped.")
    out$is_non_taxon <- TRUE
    return(out)
  }

  overrides <- .bf_tax_standardise_manual_overrides(manual_overrides)
  override_name <- NA_character_
  if (nrow(overrides)) {
    idx <- match(.bf_tax_norm_key(input_name), .bf_tax_norm_key(overrides$input_name))
    if (!is.na(idx)) override_name <- overrides$replacement_name[[idx]]
  }

  query_seed <- if (!is.na(override_name) && nzchar(override_name)) override_name else input_name
  clean_name <- bf_clean_taxon_names(query_seed, max_tokens = 3, keep_genus_for_open_names = TRUE)[[1]]
  genus_level <- bf_tax_is_genus_level(query_seed)

  if (is.na(clean_name) || !nzchar(clean_name)) {
    out <- .bf_tax_empty_record(input_name, clean_name, warning = "No usable cleaned taxon name.")
    out$is_genus_level <- genus_level
    return(out)
  }

  query_name <- clean_name
  if (isTRUE(use_gbif_parse) && requireNamespace("rgbif", quietly = TRUE)) {
    parsed <- .bf_tax_gbif_name_parse_canonical(
      query_name,
      batch_size = 1,
      retries = retries,
      sleep_sec = sleep_sec,
      quiet = quiet
    )[[1]]
    if (!is.na(parsed) && nzchar(parsed)) query_name <- parsed
  }

  res <- .bf_tax_gbif_backbone(
    query_name,
    strict = strict,
    rank = if (isTRUE(genus_level)) "GENUS" else NULL,
    retries = retries,
    sleep_sec = sleep_sec,
    quiet = quiet
  )

  out <- .bf_tax_flatten_gbif(
    res,
    input_name = input_name,
    clean_name = clean_name,
    query_name = query_name,
    resolved_by = if (!is.na(override_name)) "manual_override+gbif_backbone" else "gbif_backbone"
  )

  # Rescue 1: GBIF suggest -> backbone.
  if (is.na(out$usageKey[[1]]) && isTRUE(use_gbif_suggest)) {
    suggested <- .bf_tax_gbif_suggest_best(
      query_name,
      retries = retries,
      sleep_sec = sleep_sec,
      quiet = quiet
    )

    if (!is.na(suggested) && nzchar(suggested)) {
      res2 <- .bf_tax_gbif_backbone(
        suggested,
        strict = strict,
        rank = if (isTRUE(genus_level)) "GENUS" else NULL,
        retries = retries,
        sleep_sec = sleep_sec,
        quiet = quiet
      )
      out <- .bf_tax_flatten_gbif(
        res2,
        input_name = input_name,
        clean_name = clean_name,
        query_name = suggested,
        resolved_by = "gbif_suggest+backbone"
      )
    }
  }

  # Rescue 2: explicit genus fallback.
  if (is.na(out$usageKey[[1]]) && isTRUE(use_genus_fallback)) {
    genus <- bf_tax_extract_genus(query_seed)[[1]]
    if (!is.na(genus) && nzchar(genus)) {
      res3 <- .bf_tax_gbif_backbone(
        genus,
        strict = strict,
        rank = "GENUS",
        retries = retries,
        sleep_sec = sleep_sec,
        quiet = quiet
      )
      out <- .bf_tax_flatten_gbif(
        res3,
        input_name = input_name,
        clean_name = clean_name,
        query_name = genus,
        resolved_by = "gbif_genus_fallback"
      )
    }
  }

  # Complete ranks from usageKey/classification when name_backbone omits fields.
  if (!is.na(out$taxonKey_for_download[[1]])) {
    out <- .bf_tax_complete_ranks_by_key(
      out,
      retries = retries,
      sleep_sec = sleep_sec,
      quiet = quiet
    )
  }

  # Optional external fallback for ranks only.
  if (isTRUE(use_worms_fallback)) {
    rank_cols <- c("kingdom", "phylum", "class", "order", "family")
    missing_rank <- any(vapply(out[rank_cols], function(z) is.na(z[[1]]) || !nzchar(as.character(z[[1]])), logical(1)))

    if (is.na(out$usageKey[[1]]) || missing_rank) {
      wrm <- .bf_tax_worms_resolve(
        query_name,
        input_name = input_name,
        clean_name = clean_name,
        retries = retries,
        sleep_sec = sleep_sec,
        quiet = quiet
      )

      if (!is.na(wrm$taxon_match_status[[1]]) && wrm$taxon_match_status[[1]] != "unresolved") {
        for (rk in c("kingdom", "phylum", "class", "order", "family", "genus", "species")) {
          if (rk %in% names(out) && (is.na(out[[rk]][[1]]) || !nzchar(as.character(out[[rk]][[1]]))) && !is.na(wrm[[rk]][[1]])) {
            out[[rk]] <- wrm[[rk]][[1]]
          }
        }
        if (is.na(out$usageKey[[1]])) {
          out$scientificName <- wrm$scientificName
          out$canonicalName <- wrm$canonicalName
          out$status <- wrm$status
          out$matchType <- wrm$matchType
          out$resolved_by <- "worms_fallback"
          out$taxon_match_status <- wrm$taxon_match_status
          out$taxon_match_warning <- "Resolved by WoRMS only; no GBIF taxonKey available for download."
        } else if (any(!is.na(wrm[rank_cols]))) {
          out$resolved_by <- paste0(out$resolved_by, "+worms_rank_fill")
        }
      }
    }
  }

  out$is_genus_level <- isTRUE(genus_level) || identical(out$rank[[1]], "GENUS")
  out$is_non_taxon <- FALSE

  if (!is.na(override_name) && nzchar(override_name)) {
    existing <- out$taxon_match_warning[[1]]
    add <- paste0("Manual override applied: ", input_name, " -> ", override_name, ".")
    out$taxon_match_warning <- if (is.na(existing) || !nzchar(existing)) add else paste(existing, add)
  }

  out
}

#' Resolve a vector of taxon names against GBIF
#'
#' @param names Character vector of taxon names.
#' @param cache_dir Optional cache directory. If supplied, an RDS cache is read
#'   and written there. If `NULL`, no taxonomy cache is read or written. In
#'   examples, tests and vignettes, use a path under `tempdir()` when caching is
#'   needed.
#' @param cache_file Cache file name.
#' @param refresh_cache Ignore existing cached rows.
#' @param strict Passed to `rgbif::name_backbone()`.
#' @param strict_taxonomy If TRUE, stop if any non-empty/non-junk name is unresolved.
#' @param quiet Reduce messages.
#' @param ... Passed to `bf_resolve_gbif_taxonomy()`.
#'
#' @return A tibble with one row per input name, preserving input order and using
#'   the same taxonomy-resolution columns returned by
#'   [bf_resolve_gbif_taxonomy()]. The output includes cleaned names, GBIF keys,
#'   resolved names, rank and classification fields, resolver provenance,
#'   diagnostic status columns and flags for genus-level or non-taxonomic inputs.
#'   When a cache is used, cached and newly resolved rows are combined into the
#'   same output schema.
#'
#' @export
bf_resolve_gbif_taxonomy_batch <- function(names,
                                           cache_dir = NULL,
                                           cache_file = "bf_gbif_taxonomy_cache.rds",
                                           refresh_cache = FALSE,
                                           strict = FALSE,
                                           strict_taxonomy = FALSE,
                                           quiet = TRUE,
                                           ...) {
  bf_require_packages(c("dplyr", "tibble"), context = "taxonomy helper")

  input <- as.character(names)
  if (!length(input)) return(.bf_tax_empty_record(character(0), character(0)))

  clean <- bf_clean_taxon_names(input, max_tokens = 3, keep_genus_for_open_names = TRUE)
  keys <- .bf_tax_norm_key(paste(input, clean, sep = "||"))

  cache_path <- NULL
  cache <- tibble::tibble(.cache_key = character())

  if (!is.null(cache_dir) && length(cache_dir) > 0L &&
      nzchar(trimws(as.character(cache_dir[[1L]])))) {
    cache_dir <- normalizePath(
      as.character(cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    cache_path <- file.path(cache_dir, cache_file)

    if (file.exists(cache_path) && !isTRUE(refresh_cache)) {
      cache <- tryCatch(
        readRDS(cache_path),
        error = function(e) tibble::tibble(.cache_key = character())
      )

      if (!".cache_key" %in% names(cache)) {
        cache$.cache_key <- .bf_tax_norm_key(
          paste(cache$input_name, cache$clean_name, sep = "||")
        )
      }
    }
  }

  out <- vector("list", length(input))
  need <- logical(length(input))

  for (i in seq_along(input)) {
    hit <- match(keys[[i]], cache$.cache_key)
    if (!is.na(hit)) {
      row <- cache[hit, setdiff(names(cache), ".cache_key"), drop = FALSE]
      row$input_name <- input[[i]]
      out[[i]] <- row
    } else {
      need[[i]] <- TRUE
    }
  }

  if (any(need)) {
    new_rows <- vector("list", sum(need))
    j <- 0L
    idx <- which(need)

    for (i in idx) {
      j <- j + 1L
      if (!isTRUE(quiet)) message("Resolving taxonomy: ", input[[i]])

      new_rows[[j]] <- bf_resolve_gbif_taxonomy(
        input[[i]],
        strict = strict,
        quiet = quiet,
        ...
      )
    }

    for (k in seq_along(idx)) out[[idx[[k]]]] <- new_rows[[k]]

    if (!is.null(cache_path)) {
      cache_add <- dplyr::bind_rows(new_rows)
      cache_add$.cache_key <- keys[idx]
      cache_new <- dplyr::bind_rows(cache, cache_add)
      cache_new <- cache_new[!duplicated(cache_new$.cache_key, fromLast = TRUE), , drop = FALSE]
      saveRDS(cache_new, cache_path)
    }
  }

  ans <- dplyr::bind_rows(out)

  # Preserve input order and ensure all expected columns are present.
  expected <- names(.bf_tax_empty_record("x", "x"))
  missing_cols <- setdiff(expected, names(ans))
  if (length(missing_cols)) for (nm in missing_cols) ans[[nm]] <- NA
  ans <- ans[, expected, drop = FALSE]

  if (isTRUE(strict_taxonomy)) {
    bad <- ans$taxon_match_status == "unresolved" & !ans$is_non_taxon & !is.na(ans$input_name) & nzchar(ans$input_name)
    if (any(bad, na.rm = TRUE)) {
      stop(
        "Taxonomy resolution failed for: ",
        paste(unique(ans$input_name[bad]), collapse = "; "),
        call. = FALSE
      )
    }
  }

  ans
}

#' Attach GBIF taxonomy to a data frame
#'
#' @param x Data frame.
#' @param name_col Column containing taxon names.
#' @param prefix Prefix for attached taxonomy columns. Use `""` for no prefix.
#' @param ... Passed to `bf_resolve_gbif_taxonomy_batch()`.
#'
#' @return A data frame with the same rows as `x` and GBIF taxonomy-resolution
#'   columns joined by `name_col`. The attached columns are produced by
#'   [bf_resolve_gbif_taxonomy_batch()] and, by default, are prefixed with
#'   `prefix`. They include resolved names, GBIF usage keys, rank,
#'   classification fields, match diagnostics and resolver provenance. The
#'   output is used to keep original records together with their taxonomy audit
#'   information.
#'
#' @export
bf_attach_gbif_taxonomy <- function(x,
                                    name_col = "species",
                                    prefix = "tax_",
                                    ...) {
  bf_require_packages(c("dplyr", "tibble"), context = "taxonomy helper")
  if (!is.data.frame(x)) stop("`x` must be a data frame.", call. = FALSE)
  if (!name_col %in% names(x)) stop("Name column not found: ", name_col, call. = FALSE)

  lookup <- bf_resolve_gbif_taxonomy_batch(unique(as.character(x[[name_col]])), ...)
  join_col <- paste0(".bf_join_", name_col)
  lookup[[join_col]] <- lookup$input_name

  keep <- setdiff(names(lookup), c("input_name"))
  lookup <- lookup[, c(join_col, keep), drop = FALSE]

  if (nzchar(prefix)) {
    names(lookup)[names(lookup) != join_col] <- paste0(prefix, names(lookup)[names(lookup) != join_col])
  }

  out <- x
  out[[join_col]] <- as.character(out[[name_col]])
  out <- dplyr::left_join(out, lookup, by = join_col)
  out[[join_col]] <- NULL
  out
}

#' Check taxonomic resolution quality
#'
#' @param taxonomy_tbl Output from `bf_resolve_gbif_taxonomy_batch()` or a data
#'   frame with equivalent columns.
#' @param min_confidence Optional minimum GBIF confidence threshold.
#'
#' @return A tibble containing the subset of rows from `taxonomy_tbl` with
#'   potential taxonomy-resolution issues. The output preserves the available
#'   taxonomy columns and adds an `issue` column describing the detected problem,
#'   such as `"non_taxon_skipped"`, `"unresolved"`,
#'   `"missing_gbif_taxon_key"`, `"low_confidence"`, `"non_exact_match"` or
#'   `"higher_rank_match"`. An empty tibble indicates that no rows were flagged
#'   by the selected checks.
#'
#' @export
bf_check_taxonomic_resolution <- function(taxonomy_tbl,
                                          min_confidence = 90) {
  bf_require_packages(c("dplyr", "tibble"), context = "taxonomy helper")
  if (!is.data.frame(taxonomy_tbl)) stop("`taxonomy_tbl` must be a data frame.", call. = FALSE)

  x <- taxonomy_tbl
  needed <- c("input_name", "taxon_match_status", "matchType", "confidence", "rank", "taxonKey_for_download", "is_non_taxon")
  missing <- setdiff(needed, names(x))
  if (length(missing)) for (nm in missing) x[[nm]] <- NA

  x |>
    dplyr::mutate(
      issue = dplyr::case_when(
        .data$is_non_taxon %in% TRUE ~ "non_taxon_skipped",
        .data$taxon_match_status == "unresolved" ~ "unresolved",
        is.na(.data$taxonKey_for_download) ~ "missing_gbif_taxon_key",
        !is.na(.data$confidence) & .data$confidence < min_confidence ~ "low_confidence",
        !is.na(.data$matchType) & !.data$matchType %in% c("EXACT", "LOOKUP") ~ "non_exact_match",
        !is.na(.data$rank) & .data$rank %in% c("GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM") ~ "higher_rank_match",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$issue))
}

#' Summarise taxonomy resolution outcomes
#'
#' @param taxonomy_tbl Output from `bf_resolve_gbif_taxonomy_batch()`.
#'
#' @return A one-row tibble summarising taxonomy-resolution outcomes. The output
#'   contains `n_input`, `n_matched`, `n_unresolved`, `n_non_taxon`,
#'   `n_higher_rank`, `pct_matched`, `resolved_by` and `match_types`. These
#'   fields report the number of input rows, how many received usable GBIF taxon
#'   keys, how many were unresolved or skipped, and which resolver and match-type
#'   categories were present.
#'
#' @export
bf_taxonomic_summary <- function(taxonomy_tbl) {
  bf_require_packages(c("dplyr", "tibble"), context = "taxonomy helper")
  if (!is.data.frame(taxonomy_tbl)) stop("`taxonomy_tbl` must be a data frame.", call. = FALSE)

  x <- taxonomy_tbl
  needed <- c("taxon_match_status", "resolved_by", "rank", "matchType", "is_non_taxon", "taxonKey_for_download")
  missing <- setdiff(needed, names(x))
  if (length(missing)) for (nm in missing) x[[nm]] <- NA

  n_input <- nrow(x)
  n_matched <- sum(!is.na(x$taxonKey_for_download), na.rm = TRUE)
  n_unresolved <- sum(x$taxon_match_status == "unresolved", na.rm = TRUE)
  n_non_taxon <- sum(x$is_non_taxon %in% TRUE, na.rm = TRUE)
  n_higher_rank <- sum(x$rank %in% c("GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM"), na.rm = TRUE)

  tibble::tibble(
    n_input = n_input,
    n_matched = n_matched,
    n_unresolved = n_unresolved,
    n_non_taxon = n_non_taxon,
    n_higher_rank = n_higher_rank,
    pct_matched = if (n_input > 0L) round(100 * n_matched / n_input, 1) else NA_real_,
    resolved_by = paste(sort(unique(stats::na.omit(x$resolved_by))), collapse = ";"),
    match_types = paste(sort(unique(stats::na.omit(x$matchType))), collapse = ";")
  )
}

# -----------------------------------------------------------------------------
# Backward-compatible wrappers
# -----------------------------------------------------------------------------

#' Check for NCBI Entrez API Key
#'
#' @return Invisibly returns a logical value of length one. The value is `TRUE`
#'   when the `ENTREZ_KEY` environment variable is present and non-empty, and
#'   `FALSE` otherwise. The function is called primarily for its side effect of
#'   reporting whether an NCBI Entrez API key is available for optional
#'   taxonomy-related workflows.
#'
#' @export
check_entrez_key <- function() {
  key <- Sys.getenv("ENTREZ_KEY")
  ok <- nzchar(key)

  if (!ok) {
    if (requireNamespace("cli", quietly = TRUE)) {
      cli::cli_alert_warning("No ENTREZ API key found. NCBI/taxize calls may be rate-limited.")
      cli::cli_alert_info("Set one with Sys.setenv(ENTREZ_KEY = 'YOUR_API_KEY') or taxize::use_entrez().")
    } else {
      message("No ENTREZ API key found. NCBI/taxize calls may be rate-limited.")
    }
    Sys.setenv(ENTREZ_KEY = "")
  } else {
    if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_success("Entrez API key detected.") else message("Entrez API key detected.")
  }

  invisible(ok)
}

#' Resolve and standardise species names
#'
#' Backward-compatible wrapper around `bf_resolve_gbif_taxonomy_batch()`.
#'
#' @param names Character vector of species names.
#' @param sources Retained for compatibility. Currently not used by the GBIF-first resolver.
#' @param use_taxize Retained for compatibility.
#' @param verbose Print progress messages.
#'
#' @return A tibble with one row per input name. The output includes
#'   backward-compatible columns `original_name`, `resolved_name`,
#'   `match_type`, `source` and `score`, followed by the full GBIF-first
#'   taxonomy-resolution columns returned by [bf_resolve_gbif_taxonomy_batch()].
#'   The table provides a compatibility interface for older workflows while
#'   retaining the current biofetchR taxonomy audit fields.
#'
#' @export
resolve_species_names <- function(names,
                                  sources = c("GBIF", "ITIS"),
                                  use_taxize = TRUE,
                                  verbose = TRUE) {
  tax <- bf_resolve_gbif_taxonomy_batch(
    names,
    use_worms_fallback = FALSE,
    quiet = !isTRUE(verbose)
  )

  dplyr::bind_cols(
    tibble::tibble(
      original_name = tax$input_name,
      resolved_name = tax$scientificName,
      match_type = tax$matchType,
      source = tax$resolved_by,
      score = tax$confidence
    ),
    tax
  )
}

#' Determine if a species is likely marine from taxonomy
#'
#' Uses GBIF taxonomy first. If needed, falls back to broad marine clades and a
#' small genus-name heuristic.
#'
#' @param species_name Scientific name.
#' @param verbose Print result.
#'
#' @return A logical value of length one. The value is `TRUE` when the supplied
#'   species name is identified as likely marine from resolved taxonomy clades or
#'   a small genus-name heuristic, and `FALSE` otherwise. When `verbose = TRUE`,
#'   the function also reports the classification result as a message.
#'
#' @export
is_marine_species <- function(species_name, verbose = FALSE) {
  marine_clades <- c(
    "Echinodermata", "Cnidaria", "Porifera", "Mollusca", "Arthropoda",
    "Crustacea", "Chondrichthyes", "Actinopterygii", "Elasmobranchii",
    "Cephalopoda", "Bivalvia", "Gastropoda", "Decapoda", "Bryozoa",
    "Ascidiacea", "Hydrozoa", "Anthozoa", "Foraminifera"
  )

  marine_genera <- c(
    "Mytilus", "Scomber", "Gadus", "Pomacentrus", "Sepia", "Engraulis",
    "Carcinus", "Mnemiopsis", "Undaria", "Didemnum", "Styela", "Asterias"
  )

  tax <- tryCatch(
    bf_resolve_gbif_taxonomy(species_name, quiet = TRUE),
    error = function(e) NULL
  )

  clades <- character(0)
  if (is.data.frame(tax) && nrow(tax)) {
    clades <- stats::na.omit(as.character(unlist(tax[, c("kingdom", "phylum", "class", "order", "family", "genus"), drop = FALSE])))
  }

  by_clade <- any(clades %in% marine_clades)
  gen <- bf_tax_extract_genus(species_name)[[1]]
  by_name <- !is.na(gen) && gen %in% marine_genera

  out <- isTRUE(by_clade || by_name)

  if (isTRUE(verbose)) {
    msg <- paste0(species_name, ": ", if (out) "likely marine" else "not identified as marine")
    if (requireNamespace("cli", quietly = TRUE)) {
      if (out) cli::cli_alert_success(msg) else cli::cli_alert_warning(msg)
    } else {
      message(msg)
    }
  }

  out
}
