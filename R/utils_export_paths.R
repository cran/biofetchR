################################################################################
# utils_export_paths.R
# -----------------------------------------------------------------------------
# biofetchR: universal export-path and manifest helpers
# -----------------------------------------------------------------------------
# PURPOSE
#   Provide a single deterministic filename contract for every occurrence-table
#   export produced by biofetchR. These helpers keep the physical CSV path, the
#   `gbif_summary$output_file` entry, the optional `export_manifest.csv` row and
#   regression-test expectations in sync.
#
# WHY THIS EXISTS
#   Earlier pipeline stages can legitimately write different kinds of outputs:
#   raw country-level occurrence tables, GADM-assigned tables, HydroBASINS-style
#   freshwater outputs, EEZ/marine-region exports, overlay-enriched exports,
#   cleaned-only exports, spatially thinned exports and native-evidence-filtered
#   outputs. If each pipeline constructs filenames independently, it is easy for
#   the package to write one file but record or test another expected filename.
#   This file centralises that logic.
#
# CORE RESPONSIBILITIES
#   1. Convert species names, region labels and workflow settings into safe,
#      deterministic filename components.
#   2. Infer the *actual* export stage from the object being written, not just
#      from requested pipeline settings.
#   3. Build complete export plans containing output paths and audit metadata.
#   4. Write occurrence CSVs from those plans.
#   5. Append stable rows to `export_manifest.csv`.
#   6. Build `gbif_summary` rows that point to the same resolved output file.
#
# DESIGN PRINCIPLES
#   - Do not duplicate filename logic inside individual pipelines.
#   - Keep default filenames short and readable.
#   - Add modifier tags only when settings materially change the output table.
#   - Avoid path separators, punctuation, whitespace and problematic encodings in
#     generated filenames.
#   - Preserve backwards compatibility through `.expected_output_file()`.
#
# DATA AND LICENSING NOTE
#   This script does not download, query, transform or redistribute external
#   biodiversity or spatial datasets. It only constructs paths and small audit
#   tables for outputs produced elsewhere in biofetchR. Licensing and citation
#   requirements therefore belong to the data-provider-specific pipeline and
#   overlay scripts that created the occurrence records being exported.
#
################################################################################

#' Return TRUE for blank scalar-like values
#'
#' Internal helper used by the export-path resolver. A value is treated as blank
#' when it is `NULL`, length zero, `NA`, or an empty string after trimming.
#'
#' @param x Object to test.
#'
#' @return A single logical value.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.bf_is_blank <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(TRUE)
  }

  x <- x[[1L]]

  if (is.na(x)) {
    return(TRUE)
  }

  !nzchar(trimws(as.character(x)))
}


#' Return a single non-blank value or a fallback
#'
#' @param x Object to coerce to a length-one character value.
#' @param fallback Value returned when `x` is blank.
#'
#' @return A length-one character value.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.bf_chr1 <- function(x, fallback = NA_character_) {
  if (.bf_is_blank(x)) {
    return(as.character(fallback)[[1L]])
  }

  as.character(x[[1L]])
}


#' Return a single non-blank value from several candidates
#'
#' @param ... Candidate values, ordered from highest to lowest priority.
#' @param fallback Value returned if all candidates are blank.
#'
#' @return A length-one character value.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.bf_first_nonblank <- function(..., fallback = NA_character_) {
  values <- list(...)

  for (value in values) {
    if (!.bf_is_blank(value)) {
      return(.bf_chr1(value))
    }
  }

  as.character(fallback)[[1L]]
}


#' Convert arbitrary text into a filesystem-safe slug
#'
#' This helper standardises species names, region identifiers and modifier tags
#' for use in output filenames. It avoids spaces, punctuation, path separators and
#' non-ASCII characters where possible.
#'
#' @param x Character vector to slugify.
#' @param max_chars Maximum number of characters retained per slug.
#' @param fallback Replacement used for blank values.
#'
#' @return A character vector of safe filename components.
#'
#' @examples
#' bf_safe_slug("Xenopus laevis")
#' bf_safe_slug("C\\u00f4te d'Ivoire")
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_safe_slug <- function(x, max_chars = 80L, fallback = "NA") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- fallback

  # Transliterate where possible so names such as C\\u00f4te d'Ivoire remain readable.
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_")
  x[is.na(x) | !nzchar(trimws(x))] <- fallback

  # Keep only letters and numbers, using underscores as stable separators.
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[is.na(x) | !nzchar(trimws(x))] <- fallback

  max_chars <- as.integer(max_chars[[1L]])
  if (is.na(max_chars) || max_chars < 1L) {
    max_chars <- 80L
  }

  too_long <- nchar(x) > max_chars
  x[too_long] <- substr(x[too_long], 1L, max_chars)
  x
}


#' Normalise a GADM level for filename use
#'
#' @param gadm_unit GADM level using biofetchR's 0-based convention: level 0 is
#'   country, level 1 is first administrative subdivision, and so on.
#'
#' @return A character scalar such as `"0"`, `"1"`, or `NA_character_`.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.bf_normalise_gadm_unit <- function(gadm_unit = NULL) {
  if (.bf_is_blank(gadm_unit)) {
    return(NA_character_)
  }

  x <- suppressWarnings(as.integer(gadm_unit[[1L]]))

  if (is.na(x)) {
    return(bf_safe_slug(gadm_unit, max_chars = 10L))
  }

  as.character(x)
}


#' Normalise a HydroBASINS level for filename use
#'
#' @param hydro_level HydroBASINS level, usually 04, 05, 06, etc.
#'
#' @return A two-character level string such as `"04"`, or `NA_character_`.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.bf_normalise_hydro_level <- function(hydro_level = NULL) {
  if (.bf_is_blank(hydro_level)) {
    return(NA_character_)
  }

  x <- suppressWarnings(as.integer(hydro_level[[1L]]))

  if (is.na(x)) {
    out <- bf_safe_slug(hydro_level, max_chars = 10L)
    out <- gsub("^[Ll]", "", out)
    return(out)
  }

  sprintf("%02d", x)
}


#' Infer the actual biofetchR export stage
#'
#' The export stage describes the structure of the file actually being written,
#' not merely the requested pipeline setting. This distinction matters because a
#' pipeline may request GADM joins but still write an intermediate raw-country
#' file first.
#'
#' Precedence is:
#'   1. explicit `export_stage` supplied by the caller;
#'   2. actual `workflow` value attached to the object being written;
#'   3. actual `region_type` value attached to the object being written;
#'   4. requested settings (`spatial_join_type`, `region_source`, `gadm_unit`);
#'   5. fallback to `"raw"`.
#'
#' @param export_stage Explicit stage, if already known.
#' @param workflow Actual workflow label in the data being written, e.g. `"raw"`,
#'   `"gadm"`, `"eez"`.
#' @param region_type Actual region type in the data being written, e.g.
#'   `"RAW_COUNTRY"`, `"GADM_L1"`.
#' @param region_source Requested region source, e.g. `"gadm"`, `"eez"`.
#' @param spatial_join_type Requested spatial-join mode, e.g. `"auto"`,
#'   `"gadm"`, `"eez"`.
#' @param gadm_unit GADM level using biofetchR's 0-based convention.
#' @param hydro_level HydroBASINS level.
#'
#' @return A lowercase filename-safe stage such as `"raw"`, `"gadm_l1"`,
#'   `"eez"`, or `"hbasin_l04"`.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_infer_export_stage <- function(
    export_stage = NULL,
    workflow = NULL,
    region_type = NULL,
    region_source = NULL,
    spatial_join_type = NULL,
    gadm_unit = NULL,
    hydro_level = NULL
) {
  explicit <- tolower(.bf_chr1(export_stage, fallback = ""))
  explicit <- gsub("[^a-z0-9]+", "_", explicit)
  explicit <- gsub("^_+|_+$", "", explicit)

  if (nzchar(explicit)) {
    return(explicit)
  }

  workflow_chr <- tolower(.bf_chr1(workflow, fallback = ""))
  workflow_chr <- gsub("[^a-z0-9]+", "_", workflow_chr)
  workflow_chr <- gsub("^_+|_+$", "", workflow_chr)

  if (workflow_chr %in% c("raw", "raw_country", "country")) {
    return("raw")
  }

  if (workflow_chr %in% c("gadm", "gadm_join", "gadm_grouped")) {
    lev <- .bf_normalise_gadm_unit(gadm_unit)
    return(if (!is.na(lev)) paste0("gadm_l", lev) else "gadm")
  }

  if (workflow_chr %in% c("eez", "marine", "marine_eez")) {
    return("eez")
  }

  if (workflow_chr %in% c("hbasin", "hbasins", "hydrobasin", "hydrobasins")) {
    lev <- .bf_normalise_hydro_level(hydro_level)
    return(if (!is.na(lev)) paste0("hbasin_l", lev) else "hbasin")
  }

  region_type_chr <- toupper(.bf_chr1(region_type, fallback = ""))
  region_type_chr <- gsub("[^A-Z0-9]+", "_", region_type_chr)
  region_type_chr <- gsub("^_+|_+$", "", region_type_chr)

  if (region_type_chr %in% c("RAW", "RAW_COUNTRY", "COUNTRY")) {
    return("raw")
  }

  if (grepl("^GADM(_LEVEL)?_?L?[0-9]+$", region_type_chr)) {
    lev <- sub(".*?([0-9]+)$", "\\1", region_type_chr)
    return(paste0("gadm_l", as.integer(lev)))
  }

  if (region_type_chr %in% c("EEZ", "MARINE_EEZ")) {
    return("eez")
  }

  if (grepl("^(HBASIN|HBASINS|HYDROBASIN|HYDROBASINS).*?([0-9]+)$", region_type_chr)) {
    lev <- sub(".*?([0-9]+)$", "\\1", region_type_chr)
    lev <- sprintf("%02d", as.integer(lev))
    return(paste0("hbasin_l", lev))
  }

  spatial_join_type_chr <- tolower(.bf_chr1(spatial_join_type, fallback = ""))
  region_source_chr <- tolower(.bf_chr1(region_source, fallback = ""))

  if (spatial_join_type_chr == "eez" || region_source_chr == "eez") {
    return("eez")
  }

  if (spatial_join_type_chr == "gadm" || region_source_chr == "gadm") {
    lev <- .bf_normalise_gadm_unit(gadm_unit)
    return(if (!is.na(lev)) paste0("gadm_l", lev) else "gadm")
  }

  if (region_source_chr %in% c("hbasin", "hbasins", "hydrobasin", "hydrobasins")) {
    lev <- .bf_normalise_hydro_level(hydro_level)
    return(if (!is.na(lev)) paste0("hbasin_l", lev) else "hbasin")
  }

  "raw"
}


#' Create a compact overlay filename tag
#'
#' Overlay names can be long and numerous. This helper creates a short readable
#' tag and records the full overlay list in `export_manifest.csv` instead.
#'
#' @param overlay_names Character vector of overlay/context-layer names.
#' @param max_items Maximum number of overlay names to include before adding a
#'   `plusN` suffix.
#' @param max_chars Maximum characters retained per individual overlay name.
#'
#' @return A character scalar such as `"ov-rivers-lakes-wdpa"`, or
#'   `NA_character_` when no overlays are supplied.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_overlay_tag <- function(overlay_names = NULL, max_items = 4L, max_chars = 20L) {
  if (is.null(overlay_names) || length(overlay_names) == 0L) {
    return(NA_character_)
  }

  overlay_names <- sort(unique(as.character(overlay_names)))
  overlay_names <- overlay_names[!is.na(overlay_names) & nzchar(trimws(overlay_names))]

  if (length(overlay_names) == 0L) {
    return(NA_character_)
  }

  clean <- bf_safe_slug(overlay_names, max_chars = max_chars, fallback = "overlay")

  max_items <- as.integer(max_items[[1L]])
  if (is.na(max_items) || max_items < 1L) {
    max_items <- 4L
  }

  if (length(clean) > max_items) {
    shown <- clean[seq_len(max_items)]
    extra <- length(clean) - max_items
    return(paste0("ov-", paste(shown, collapse = "-"), "-plus", extra))
  }

  paste0("ov-", paste(clean, collapse = "-"))
}


#' Build optional setting tags for output filenames
#'
#' These tags are deliberately conservative. Defaults are not tagged, which keeps
#' existing simple outputs readable and backwards-compatible. Tags are added when
#' settings materially change what the exported occurrence table contains.
#'
#' @param apply_cleaning Whether occurrence cleaning was applied. `FALSE` adds
#'   `"uncleaned"`; `TRUE` is treated as the default and adds no tag.
#' @param apply_thinning Whether spatial thinning was applied. `TRUE` adds a
#'   thinning tag such as `"thin5km"`.
#' @param thinning_dist_km Thinning distance in kilometres.
#' @param use_overlays Whether overlay/context columns were attached.
#' @param overlay_names Names of overlay/context layers.
#' @param prepare_taxonomy Whether taxonomy preparation was applied. `FALSE` adds
#'   `"notax"`; `TRUE` is treated as the default.
#' @param use_native_range Whether table-based native-range filtering was used.
#' @param use_native_web Whether web/SInAS native-range evidence was used.
#' @param extra_tags Additional user-defined tags to append.
#'
#' @return A character vector of filename-safe modifier tags.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_export_setting_tags <- function(
    apply_cleaning = NULL,
    apply_thinning = NULL,
    thinning_dist_km = NULL,
    use_overlays = NULL,
    overlay_names = NULL,
    prepare_taxonomy = NULL,
    use_native_range = NULL,
    use_native_web = NULL,
    extra_tags = NULL
) {
  tags <- character()

  if (identical(apply_cleaning, FALSE)) {
    tags <- c(tags, "uncleaned")
  }

  if (identical(apply_thinning, TRUE)) {
    dist <- suppressWarnings(as.numeric(.bf_chr1(thinning_dist_km, fallback = NA_character_)))
    if (!is.na(dist)) {
      dist_tag <- gsub("\\.", "p", as.character(dist))
      tags <- c(tags, paste0("thin", dist_tag, "km"))
    } else {
      tags <- c(tags, "thinned")
    }
  }

  if (identical(use_overlays, TRUE)) {
    ov_tag <- bf_overlay_tag(overlay_names)
    if (!is.na(ov_tag) && nzchar(ov_tag)) {
      tags <- c(tags, ov_tag)
    } else {
      tags <- c(tags, "overlays")
    }
  }

  if (identical(prepare_taxonomy, FALSE)) {
    tags <- c(tags, "notax")
  }

  if (identical(use_native_range, TRUE)) {
    tags <- c(tags, "native-table")
  }

  if (identical(use_native_web, TRUE)) {
    tags <- c(tags, "native-web")
  }

  if (!is.null(extra_tags) && length(extra_tags) > 0L) {
    extra_tags <- as.character(extra_tags)
    extra_tags <- extra_tags[!is.na(extra_tags) & nzchar(trimws(extra_tags))]
    tags <- c(tags, extra_tags)
  }

  tags <- unique(bf_safe_slug(tags, max_chars = 60L, fallback = "tag"))
  tags[!is.na(tags) & nzchar(tags)]
}


#' Resolve a complete biofetchR export plan
#'
#' This is the central filename contract for biofetchR. Call this immediately
#' before writing any occurrence CSV, then use the returned `output_file` both for
#' the write operation and for `gbif_summary$output_file`.
#'
#' The returned plan is intentionally verbose. In addition to the actual file
#' path, it stores the inferred export stage, safe filename components, modifier
#' tags and key workflow settings so that export manifests and summaries can be
#' reconstructed without guessing how the filename was built.
#'
#' @param output_dir Directory where the export should be written.
#' @param species Species name for the file being written.
#' @param input_region Input country/region code used for the GBIF query. For
#'   global marine workflows this can be `"GLOBAL"` or `NA`.
#' @param export_stage Explicit actual export stage, if already known.
#' @param workflow Actual workflow label in the data being written.
#' @param region_type Actual region type in the data being written.
#' @param region_source Requested spatial region source.
#' @param spatial_join_type Requested spatial join type.
#' @param gadm_unit GADM level using biofetchR's 0-based convention.
#' @param hydro_level HydroBASINS level.
#' @param region_id Optional output region identifier, used only when
#'   `split_by_region = TRUE`.
#' @param split_by_region If `TRUE`, include `region_id` in the filename. Leave
#'   `FALSE` when the file contains all regions for one species/input region.
#' @param apply_cleaning,apply_thinning,thinning_dist_km,use_overlays,overlay_names,prepare_taxonomy,use_native_range,use_native_web
#'   Settings used to construct compact modifier tags.
#' @param extra_tags Additional filename tags for specialised tests/workflows.
#' @param extension File extension, usually `"csv"`.
#'
#' @return A named list with `output_file`, `output_basename`, `export_stage`,
#'   and all metadata needed for summaries/manifests.
#'
#' @examples
#' bf_resolve_export_plan(
#'   output_dir = tempdir(),
#'   species = "Xenopus laevis",
#'   input_region = "FR",
#'   workflow = "raw",
#'   region_type = "RAW_COUNTRY"
#' )$output_basename
#'
#' bf_resolve_export_plan(
#'   output_dir = tempdir(),
#'   species = "Xenopus laevis",
#'   input_region = "FR",
#'   export_stage = "gadm_l1",
#'   use_overlays = TRUE,
#'   overlay_names = c("rivers", "lakes", "wdpa")
#' )$output_basename
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_resolve_export_plan <- function(
    output_dir,
    species,
    input_region = NULL,
    export_stage = NULL,
    workflow = NULL,
    region_type = NULL,
    region_source = NULL,
    spatial_join_type = NULL,
    gadm_unit = NULL,
    hydro_level = NULL,
    region_id = NULL,
    split_by_region = FALSE,
    apply_cleaning = NULL,
    apply_thinning = NULL,
    thinning_dist_km = NULL,
    use_overlays = NULL,
    overlay_names = NULL,
    prepare_taxonomy = NULL,
    use_native_range = NULL,
    use_native_web = NULL,
    extra_tags = NULL,
    extension = "csv"
) {
  if (.bf_is_blank(output_dir)) {
    stop("`output_dir` must be a non-empty directory path.", call. = FALSE)
  }

  if (.bf_is_blank(species)) {
    stop("`species` must be supplied when resolving an export path.", call. = FALSE)
  }

  extension <- bf_safe_slug(extension, max_chars = 10L, fallback = "csv")

  stage <- bf_infer_export_stage(
    export_stage = export_stage,
    workflow = workflow,
    region_type = region_type,
    region_source = region_source,
    spatial_join_type = spatial_join_type,
    gadm_unit = gadm_unit,
    hydro_level = hydro_level
  )

  species_tag <- bf_safe_slug(species, max_chars = 80L, fallback = "species")

  input_region_value <- .bf_chr1(input_region, fallback = "GLOBAL")
  input_region_tag <- bf_safe_slug(input_region_value, max_chars = 40L, fallback = "GLOBAL")

  parts <- c(stage, species_tag, input_region_tag)

  if (isTRUE(split_by_region) && !.bf_is_blank(region_id)) {
    parts <- c(parts, bf_safe_slug(region_id, max_chars = 80L, fallback = "region"))
  }

  setting_tags <- bf_export_setting_tags(
    apply_cleaning = apply_cleaning,
    apply_thinning = apply_thinning,
    thinning_dist_km = thinning_dist_km,
    use_overlays = use_overlays,
    overlay_names = overlay_names,
    prepare_taxonomy = prepare_taxonomy,
    use_native_range = use_native_range,
    use_native_web = use_native_web,
    extra_tags = extra_tags
  )

  if (length(setting_tags) > 0L) {
    parts <- c(parts, setting_tags)
  }

  output_basename <- paste0(paste(parts, collapse = "__"), ".", extension)
  output_file <- file.path(output_dir, output_basename)

  list(
    output_dir = output_dir,
    output_file = output_file,
    output_basename = output_basename,
    export_stage = stage,
    species = .bf_chr1(species),
    species_tag = species_tag,
    input_region = input_region_value,
    input_region_tag = input_region_tag,
    region_id = .bf_chr1(region_id, fallback = NA_character_),
    split_by_region = isTRUE(split_by_region),
    setting_tags = paste(setting_tags, collapse = ";"),
    overlay_tag = if (identical(use_overlays, TRUE)) bf_overlay_tag(overlay_names) else NA_character_,
    overlay_names = if (is.null(overlay_names)) NA_character_ else paste(sort(unique(as.character(overlay_names))), collapse = ";"),
    apply_cleaning = apply_cleaning,
    apply_thinning = apply_thinning,
    thinning_dist_km = thinning_dist_km,
    use_overlays = use_overlays,
    prepare_taxonomy = prepare_taxonomy,
    use_native_range = use_native_range,
    use_native_web = use_native_web,
    extension = extension
  )
}


#' Write an occurrence CSV using a resolved export plan
#'
#' Writes the supplied table to `export_plan$output_file` and appends simple
#' write diagnostics back onto the plan. This helper deliberately does not alter
#' the input table; callers should drop geometry, select columns or repair
#' encodings before calling it if those steps are required.
#'
#' @param x Data frame to write.
#' @param export_plan List returned by `bf_resolve_export_plan()`.
#'
#' @return The original export plan with write diagnostics appended.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_write_export_csv <- function(x, export_plan) {
  if (is.null(export_plan$output_file) || !nzchar(export_plan$output_file)) {
    stop("`export_plan$output_file` is missing.", call. = FALSE)
  }

  dir.create(dirname(export_plan$output_file), recursive = TRUE, showWarnings = FALSE)

  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(x, export_plan$output_file)
  } else {
    utils::write.csv(x, export_plan$output_file, row.names = FALSE)
  }

  export_plan$export_written <- file.exists(export_plan$output_file)
  export_plan$export_n_rows <- if (is.data.frame(x)) nrow(x) else NA_integer_
  export_plan$export_n_cols <- if (is.data.frame(x)) ncol(x) else NA_integer_
  export_plan$export_created_at <- as.character(Sys.time())

  export_plan
}


#' Convert an export plan to a one-row manifest table
#'
#' Converts the resolved export metadata into a rectangular audit row suitable
#' for `export_manifest.csv`. The manifest is intentionally broader than
#' `gbif_summary.csv`: it records filename components and workflow settings as
#' well as the final output path.
#'
#' @param export_plan List returned by `bf_resolve_export_plan()` and optionally
#'   updated by `bf_write_export_csv()`.
#' @param extra Named list of additional manifest fields.
#'
#' @return A one-row tibble.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_export_manifest_row <- function(export_plan, extra = list()) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required to build export manifests.", call. = FALSE)
  }

  row <- tibble::tibble(
    output_file = .bf_chr1(export_plan$output_file, fallback = NA_character_),
    output_basename = .bf_chr1(export_plan$output_basename, fallback = NA_character_),
    export_stage = .bf_chr1(export_plan$export_stage, fallback = NA_character_),
    species = .bf_chr1(export_plan$species, fallback = NA_character_),
    input_region = .bf_chr1(export_plan$input_region, fallback = NA_character_),
    region_id = .bf_chr1(export_plan$region_id, fallback = NA_character_),
    split_by_region = isTRUE(export_plan$split_by_region),
    setting_tags = .bf_chr1(export_plan$setting_tags, fallback = NA_character_),
    overlay_tag = .bf_chr1(export_plan$overlay_tag, fallback = NA_character_),
    overlay_names = .bf_chr1(export_plan$overlay_names, fallback = NA_character_),
    apply_cleaning = if (is.null(export_plan$apply_cleaning)) NA else isTRUE(export_plan$apply_cleaning),
    apply_thinning = if (is.null(export_plan$apply_thinning)) NA else isTRUE(export_plan$apply_thinning),
    thinning_dist_km = suppressWarnings(as.numeric(.bf_chr1(export_plan$thinning_dist_km, fallback = NA_character_))),
    use_overlays = if (is.null(export_plan$use_overlays)) NA else isTRUE(export_plan$use_overlays),
    prepare_taxonomy = if (is.null(export_plan$prepare_taxonomy)) NA else isTRUE(export_plan$prepare_taxonomy),
    use_native_range = if (is.null(export_plan$use_native_range)) NA else isTRUE(export_plan$use_native_range),
    use_native_web = if (is.null(export_plan$use_native_web)) NA else isTRUE(export_plan$use_native_web),
    export_written = if (is.null(export_plan$export_written)) file.exists(.bf_chr1(export_plan$output_file, fallback = "")) else isTRUE(export_plan$export_written),
    export_n_rows = suppressWarnings(as.integer(.bf_chr1(export_plan$export_n_rows, fallback = NA_character_))),
    export_n_cols = suppressWarnings(as.integer(.bf_chr1(export_plan$export_n_cols, fallback = NA_character_))),
    export_created_at = .bf_chr1(export_plan$export_created_at, fallback = as.character(Sys.time()))
  )

  if (length(extra) > 0L) {
    extra <- as.list(extra)
    for (nm in names(extra)) {
      row[[nm]] <- extra[[nm]][[1L]]
    }
  }

  row
}


#' Append one export-plan row to export_manifest.csv
#'
#' @param export_plan List returned by `bf_resolve_export_plan()` and optionally
#'   updated by `bf_write_export_csv()`.
#' @param manifest_file Full path to the manifest CSV. If `NULL`, the manifest is
#'   written to `file.path(export_plan$output_dir, "export_manifest.csv")`.
#' @param extra Named list of additional fields to append to the manifest row.
#'
#' @return Invisibly returns the manifest path.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_append_export_manifest <- function(export_plan, manifest_file = NULL, extra = list()) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required to append export manifests.", call. = FALSE)
  }

  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package `readr` is required to append export manifests.", call. = FALSE)
  }

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required to append export manifests.", call. = FALSE)
  }

  if (is.null(manifest_file) || !nzchar(manifest_file)) {
    manifest_file <- file.path(export_plan$output_dir, "export_manifest.csv")
  }

  dir.create(dirname(manifest_file), recursive = TRUE, showWarnings = FALSE)

  new_row <- bf_export_manifest_row(export_plan, extra = extra)

  # Manifest files are audit logs, not modelling inputs. Read and bind all
  # fields as character strings to avoid hidden append failures caused by
  # readr guessing different column types across rows, especially for columns
  # that are all-NA in one row but populated in another.
  new_row <- dplyr::mutate(
    tibble::as_tibble(new_row),
    dplyr::across(dplyr::everything(), as.character)
  )

  if (file.exists(manifest_file)) {
    old <- readr::read_csv(
      manifest_file,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE
    )

    old <- dplyr::mutate(
      tibble::as_tibble(old),
      dplyr::across(dplyr::everything(), as.character)
    )

    all_cols <- union(names(old), names(new_row))

    for (nm in setdiff(all_cols, names(old))) {
      old[[nm]] <- NA_character_
    }

    for (nm in setdiff(all_cols, names(new_row))) {
      new_row[[nm]] <- NA_character_
    }

    out <- dplyr::bind_rows(
      old[, all_cols, drop = FALSE],
      new_row[, all_cols, drop = FALSE]
    )
  } else {
    out <- new_row
  }

  readr::write_csv(out, manifest_file, na = "")

  invisible(manifest_file)
}


#' Build a standard gbif_summary row from an export plan
#'
#' This helper keeps the original `gbif_summary` columns but adds harmless
#' optional audit columns. Existing tests that require only the original columns
#' should continue to pass, while newer tests can check that summaries, manifests
#' and written files all refer to the same resolved export path.
#'
#' @param species Species name.
#' @param region_id Region identifier recorded in the summary row.
#' @param region_type Region type recorded in the summary row.
#' @param n_total Number of records before cleaning/filtering.
#' @param n_cleaned Number of records after cleaning.
#' @param n_thinned Number of records after thinning.
#' @param export_plan List returned by `bf_resolve_export_plan()`.
#' @param status Pipeline status, usually `"success"` or `"failed"`.
#' @param fail_stage Failure stage, if any.
#' @param fail_reason Failure reason, if any.
#' @param gbif_key GBIF download key, if available.
#' @param extra Named list of extra columns to append.
#'
#' @return A one-row tibble suitable for binding into `gbif_summary.csv`.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
bf_make_gbif_summary_row <- function(
    species,
    region_id,
    region_type,
    n_total,
    n_cleaned,
    n_thinned,
    export_plan = NULL,
    status = "success",
    fail_stage = NA_character_,
    fail_reason = NA_character_,
    gbif_key = NA_character_,
    extra = list()
) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required to build gbif_summary rows.", call. = FALSE)
  }

  output_file <- if (!is.null(export_plan)) .bf_chr1(export_plan$output_file, fallback = NA_character_) else NA_character_
  output_basename <- if (!is.null(export_plan)) .bf_chr1(export_plan$output_basename, fallback = NA_character_) else NA_character_
  export_stage <- if (!is.null(export_plan)) .bf_chr1(export_plan$export_stage, fallback = NA_character_) else NA_character_
  export_written <- if (!is.null(export_plan) && !.bf_is_blank(output_file)) file.exists(output_file) else FALSE

  row <- tibble::tibble(
    species = .bf_chr1(species, fallback = NA_character_),
    region_id = .bf_chr1(region_id, fallback = NA_character_),
    region_type = .bf_chr1(region_type, fallback = NA_character_),
    n_total = suppressWarnings(as.integer(.bf_chr1(n_total, fallback = NA_character_))),
    n_cleaned = suppressWarnings(as.integer(.bf_chr1(n_cleaned, fallback = NA_character_))),
    n_thinned = suppressWarnings(as.integer(.bf_chr1(n_thinned, fallback = NA_character_))),
    output_file = output_file,
    status = .bf_chr1(status, fallback = "success"),
    fail_stage = .bf_chr1(fail_stage, fallback = NA_character_),
    fail_reason = .bf_chr1(fail_reason, fallback = NA_character_),
    gbif_key = .bf_chr1(gbif_key, fallback = NA_character_),
    output_basename = output_basename,
    export_stage = export_stage,
    export_written = export_written
  )

  if (length(extra) > 0L) {
    extra <- as.list(extra)
    for (nm in names(extra)) {
      row[[nm]] <- extra[[nm]][[1L]]
    }
  }

  row
}


#' Backwards-compatible expected-output helper
#'
#' This wrapper keeps older pipeline code/tests working while routing the actual
#' filename resolution through the new universal export contract.
#'
#' @param ... Arguments passed to `bf_resolve_export_plan()`.
#'
#' @return Full output-file path.
#'
#'
#' @family internal export path helpers
#' @md
#' @keywords internal
#' @noRd
.expected_output_file <- function(...) {
  bf_resolve_export_plan(...)$output_file
}
