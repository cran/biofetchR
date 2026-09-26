################################################################################
# marine_pipeline.R
# -----------------------------------------------------------------------------
# biofetchR marine GBIF processing pipeline
# -----------------------------------------------------------------------------
#
# Purpose
#   This file contains the retained marine occurrence-processing workflow for
#   biofetchR. It is focused on data download, import, spatial overlay assignment,
#   coordinate cleaning, spatial thinning and tabular export.
#
# Main exported functions
#   - process_gbif_marine_pipeline()
#   - process_gbif_eez_pipeline()  # backwards-compatible wrapper
#
# Workflow overview
#   1. Validate and standardise the input species table.
#   2. Optionally apply the taxonomy gate before live GBIF submission.
#   3. Optionally attach species-level GRIIS evidence.
#   4. Optionally attach species-level native-origin evidence.
#   5. Group retained species into outer processing batches for organisation and
#      progress reporting.
#   6. Within each outer batch, process GBIF retrieval strictly species-by-species:
#      submit one species, wait for/import that download, reconcile its retrieval
#      outcome, and only then submit the next species.
#   7. Convert successful imports to WGS84 point sf objects.
#   8. Load a package-managed marine overlay, or use a supplied sf overlay.
#   9. Assign occurrence records to marine regions by spatial containment, with
#      an intersects fallback for boundary/coastal records.
#  10. Optionally clean coordinates and spatially thin records.
#  11. Export per species x marine-region occurrence tables plus processing and
#      GBIF retrieval audit files and/or return a combined in-memory table.
#
# Marine evidence interpretation
#   Marine recipient units are usually polygons such as EEZs, MEOW ecoregions,
#   LMEs, IHO sea areas, high seas, territorial seas, Longhurst provinces or
#   marine World Heritage polygons. GRIIS and native-origin joins are therefore
#   species-level by default, not recipient-polygon-specific, unless the user
#   supplies a separate recipient-specific evidence workflow.
#
# Dependency injection
#   The `deps` argument allows tests to inject mock summary, download, import,
#   export and thinning functions without changing the production pipeline.
#
# Encoding
#   External spatial resources can contain non-ASCII names. UTF-8 repair helpers
#   keep names such as Cura\\u00e7ao, S\\u00e3o Tom\\u00e9 and R\\u00e9union safe for console output and
#   CSV export across platforms.
#
# Data access and attribution
#   This pipeline can trigger live GBIF downloads and can use third-party marine
#   spatial overlays and native-evidence sources. biofetchR does not redistribute
#   those provider datasets from this script. Users should retain GBIF download
#   keys/DOIs and cite/attribute each provider dataset used in analysis or
#   publication according to the relevant source terms.
#
################################################################################

#' Normalise marine labels for matching
#'
#' Converts labels to lower-case underscore-separated identifiers. Used for
#' overlay aliases, workflow tags and other internal matching tasks.
#'
#' @param x Character value to normalise.
#'
#' @return Normalised character value.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_norm <- function(x) {
  x <- bf_repair_utf8_one(x)
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "_", x)
}

#' Coerce marine occurrence records to WGS84 sf points
#'
#' Accepts either an sf object or a data frame with `decimalLongitude` and
#' `decimalLatitude`, repairs text columns, assigns or transforms to EPSG:4326,
#' and drops empty geometries.
#'
#' @param x Occurrence table or sf object.
#'
#' @return An sf object in EPSG:4326, or `NULL` when `x` is `NULL`.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_ensure_sf <- function(x) {
  bf_require_packages("sf", context = "marine pipeline")

  if (is.null(x)) return(NULL)

  x <- bf_repair_utf8_df(x)

  if (!inherits(x, "sf")) {
    if (all(c("decimalLongitude", "decimalLatitude") %in% names(x))) {
      x <- sf::st_as_sf(
        x,
        coords = c("decimalLongitude", "decimalLatitude"),
        crs = 4326,
        remove = FALSE
      )
    } else {
      stop("Object is not sf and lacks decimalLongitude/decimalLatitude columns.", call. = FALSE)
    }
  }

  x <- bf_sf_wgs84(x)
  x <- x[!sf::st_is_empty(x), , drop = FALSE]
  x <- bf_repair_utf8_df(x)
  x
}

#' Make marine geometries valid where possible
#'
#' Uses `sf::st_make_valid()` when available, with a buffer-zero fallback for
#' older sf installations. Non-sf inputs are returned unchanged.
#'
#' @param x sf object or other object.
#'
#' @return Geometry-repaired object when possible.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_make_valid <- function(x) {
  x <- bf_sf_make_valid(x)
  x <- bf_repair_utf8_df(x)
  x
}

#' Restore longitude and latitude columns from sf geometry
#'
#' Ensures exported occurrence tables retain `decimalLongitude` and
#' `decimalLatitude` even after sf coercion or overlay joins.
#'
#' @param x sf point object.
#'
#' @return sf object with coordinate columns where possible.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_restore_coords <- function(x) {
  if (!inherits(x, "sf") || !nrow(x)) return(x)
  if (all(c("decimalLongitude", "decimalLatitude") %in% names(x))) {
    x <- bf_repair_utf8_df(x)
    return(x)
  }
  xy <- sf::st_coordinates(sf::st_geometry(x))
  x$decimalLongitude <- as.numeric(xy[, 1])
  x$decimalLatitude  <- as.numeric(xy[, 2])
  x <- bf_repair_utf8_df(x)
  x
}

#' Bind marine output groups safely
#'
#' Combines per-species/per-region outputs that may contain different optional
#' metadata columns. Incompatible non-geometry columns are coerced to character
#' before binding to avoid failures from mixed GBIF schemas.
#'
#' @param lst List of data frames or sf objects.
#'
#' @return Tibble containing all rows and the union of columns.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_bind_rows <- function(lst) {
  lst <- Filter(Negate(is.null), lst)
  if (!length(lst)) return(tibble::tibble(decimalLongitude = numeric(), decimalLatitude = numeric()))

  # Keep sf objects as data frames when different geometries/list columns would
  # otherwise make binding fragile in test mode.
  lst <- lapply(lst, function(x) {
    x <- bf_repair_utf8_df(x)
    if (inherits(x, "sf")) {
      x <- tibble::as_tibble(x)
    } else {
      x <- tibble::as_tibble(x)
    }
    x <- bf_repair_utf8_df(x)
    x
  })

  all_cols <- unique(unlist(lapply(lst, names)))
  lst <- lapply(lst, function(x) {
    missing <- setdiff(all_cols, names(x))
    if (length(missing)) for (nm in missing) x[[nm]] <- NA
    x[, all_cols, drop = FALSE]
  })

  # Convert incompatible columns to character except geometry/list columns.
  cls <- lapply(all_cols, function(nm) unique(vapply(lst, function(x) class(x[[nm]])[[1]], character(1))))
  names(cls) <- all_cols
  needs_chr <- names(cls)[vapply(cls, function(z) length(unique(z)) > 1L && !any(z %in% c("sfc_POINT", "sfc_POLYGON", "sfc_MULTIPOLYGON", "sfc_GEOMETRY", "sfc")), logical(1))]
  lst <- lapply(lst, function(x) {
    for (nm in needs_chr) x[[nm]] <- as.character(x[[nm]])
    x <- bf_repair_utf8_df(x)
    x
  })

  out <- dplyr::bind_rows(lst)
  out <- bf_repair_utf8_df(out)
  out
}

#' Create the default marine processing summary table
#'
#' Used when the package-level `initialize_summary()` helper is unavailable or
#' when tests inject only minimal dependencies.
#'
#' @return Empty tibble with the marine summary schema.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_default_summary <- function() {
  tibble::tibble(
    species = character(),
    region_id = character(),
    region_type = character(),
    n_total = integer(),
    n_cleaned = integer(),
    n_thinned = integer(),
    output_file = character(),
    status = character(),
    fail_stage = character(),
    fail_reason = character(),
    gbif_key = character()
  )
}

#' Append one row to the default marine summary table
#'
#' Fallback implementation used when `append_summary_row()` is unavailable.
#' The production pipeline then adds failure-stage and GBIF-key fields after the
#' append step.
#'
#' @param summary_tbl Existing summary table.
#' @param species Species name.
#' @param region_id Marine-region identifier or label.
#' @param region_type Overlay type.
#' @param n_total Number of imported records before cleaning.
#' @param n_cleaned Number of records retained after cleaning.
#' @param n_thinned Number of records retained after thinning.
#' @param output_file Output file path, or `NA` for in-memory outputs.
#' @param status Processing status.
#'
#' @return Updated summary table.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_default_append <- function(summary_tbl,
                                           species,
                                           region_id,
                                           region_type,
                                           n_total = 0,
                                           n_cleaned = NA_integer_,
                                           n_thinned = 0,
                                           output_file = NA_character_,
                                           status = "unknown") {
  row <- tibble::tibble(
    species = bf_repair_utf8_one(species),
    region_id = bf_repair_utf8_one(region_id),
    region_type = bf_repair_utf8_one(region_type),
    n_total = as.integer(n_total),
    n_cleaned = as.integer(n_cleaned),
    n_thinned = as.integer(n_thinned),
    output_file = bf_repair_utf8_one(output_file, fallback = NA_character_),
    status = bf_repair_utf8_one(status)
  )

  missing_in_row <- setdiff(names(summary_tbl), names(row))
  if (length(missing_in_row)) for (nm in missing_in_row) row[[nm]] <- NA

  missing_in_tbl <- setdiff(names(row), names(summary_tbl))
  if (length(missing_in_tbl)) for (nm in missing_in_tbl) summary_tbl[[nm]] <- NA

  out <- dplyr::bind_rows(summary_tbl, row[, names(summary_tbl), drop = FALSE])
  out <- bf_repair_utf8_df(out)
  out
}


#' Apply the marine taxonomy gate
#'
#' Cleans marine species names before GBIF submission, applies optional manual
#' fixes, rejects non-taxon placeholders and open/higher-rank names, writes audit
#' files when requested, and returns only species-level candidates. Accepted rows
#' are returned with the cleaned name in `taxonomy_name_col`, so downstream GBIF
#' downloads, summaries and filenames use the cleaned accepted name.
#'
#' @param df Input species table.
#' @param output_dir Directory for optional taxonomy audit outputs.
#' @param prepare_taxonomy Logical; run taxonomy preparation when `TRUE`.
#' @param manual_taxonomy_fixes Optional named character vector or data frame of
#'   manual corrections passed to `bf_apply_manual_taxonomy_fixes()`.
#' @param taxonomy_name_col Name of the column containing taxon names.
#' @param export_taxonomy_audit Logical; write `taxonomy_audit.csv`,
#'   `taxonomy_rejected.csv` and `taxonomy_summary.csv`.
#' @param quiet Logical; suppress progress messages.
#'
#' @return Tibble containing accepted rows with cleaned species names.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_apply_taxonomy_gate <- function(df,
                                                output_dir,
                                                prepare_taxonomy = FALSE,
                                                manual_taxonomy_fixes = NULL,
                                                taxonomy_name_col = "species",
                                                export_taxonomy_audit = TRUE,
                                                quiet = FALSE) {
  df <- bf_repair_utf8_df(df)

  if (!isTRUE(prepare_taxonomy)) return(tibble::as_tibble(df))

  for (fn in c(
    "bf_apply_manual_taxonomy_fixes",
    "bf_clean_taxon_names",
    "bf_tax_looks_non_taxon",
    "bf_tax_is_genus_level",
    "bf_tax_extract_genus"
  )) {
    if (!exists(fn, mode = "function", inherits = TRUE)) {
      stop(
        "`prepare_taxonomy = TRUE`, but required taxonomy helper is missing: ",
        fn,
        call. = FALSE
      )
    }
  }

  if (!taxonomy_name_col %in% names(df)) {
    stop("`taxonomy_name_col` was not found in `df`: ", taxonomy_name_col, call. = FALSE)
  }

  tax_raw <- tibble::as_tibble(df)
  tax_raw <- bf_repair_utf8_df(tax_raw)

  tax_fixed <- bf_apply_manual_taxonomy_fixes(
    tax_raw,
    manual_fixes = manual_taxonomy_fixes,
    name_col = taxonomy_name_col,
    overwrite_names = TRUE
  )
  tax_fixed <- bf_repair_utf8_df(tax_fixed)

  tax_audit <- tax_fixed |>
    dplyr::mutate(
      raw_name = tax_raw[[taxonomy_name_col]],
      fixed_name = .data[[taxonomy_name_col]],
      cleaned_name = bf_clean_taxon_names(.data[[taxonomy_name_col]]),
      was_manual_fixed = .data$raw_name != .data$fixed_name,
      is_non_taxon = bf_tax_looks_non_taxon(.data$fixed_name) | is.na(.data$cleaned_name),
      is_open_or_higher = dplyr::if_else(
        .data$is_non_taxon,
        FALSE,
        bf_tax_is_genus_level(.data$cleaned_name)
      ),
      genus = bf_tax_extract_genus(.data$cleaned_name),
      species_level_candidate = !.data$is_non_taxon &
        !.data$is_open_or_higher &
        !is.na(.data$cleaned_name) &
        grepl("^[A-Z][A-Za-z.-]+[[:space:]]+[a-z][A-Za-z.-]+", .data$cleaned_name),
      accepted = .data$species_level_candidate,
      rejection_reason = dplyr::case_when(
        .data$accepted ~ NA_character_,
        .data$is_non_taxon ~ "non_taxon_or_placeholder",
        .data$is_open_or_higher ~ "open_name_or_higher_rank",
        is.na(.data$cleaned_name) ~ "missing_cleaned_name",
        TRUE ~ "other_not_species_level"
      )
    )

  tax_audit <- bf_repair_utf8_df(tax_audit)

  taxonomy_rejected <- tax_audit |>
    dplyr::filter(!.data$accepted)

  taxonomy_rejected <- bf_repair_utf8_df(taxonomy_rejected)

  n_input <- nrow(tax_audit)
  n_accepted <- sum(tax_audit$accepted, na.rm = TRUE)

  taxonomy_summary <- tibble::tibble(
    n_input = n_input,
    n_accepted = n_accepted,
    n_rejected = sum(!tax_audit$accepted, na.rm = TRUE),
    n_manual_fixed = sum(tax_audit$was_manual_fixed, na.rm = TRUE),
    n_non_taxon = sum(tax_audit$is_non_taxon, na.rm = TRUE),
    n_open_or_higher = sum(tax_audit$is_open_or_higher, na.rm = TRUE),
    n_species_level_candidates = sum(tax_audit$species_level_candidate, na.rm = TRUE),
    pct_accepted = if (n_input > 0L) round(100 * n_accepted / n_input, 1) else NA_real_
  )

  taxonomy_summary <- bf_repair_utf8_df(taxonomy_summary)

  if (isTRUE(export_taxonomy_audit)) {
    readr::write_csv(tax_audit, file.path(output_dir, "taxonomy_audit.csv"))
    readr::write_csv(taxonomy_rejected, file.path(output_dir, "taxonomy_rejected.csv"))
    readr::write_csv(taxonomy_summary, file.path(output_dir, "taxonomy_summary.csv"))
  }

  accepted_idx <- which(tax_audit$accepted %in% TRUE)
  out <- tax_fixed[accepted_idx, , drop = FALSE]
  out[[taxonomy_name_col]] <- tax_audit$cleaned_name[accepted_idx]
  out <- bf_repair_utf8_df(out)

  if (!isTRUE(quiet)) {
    message(
      "biofetchR taxonomy gate: accepted ",
      nrow(out),
      " marine species row(s); rejected ",
      nrow(taxonomy_rejected),
      " row(s)."
    )
  }

  if (nrow(out) == 0L) {
    stop(
      "No rows remained after the marine taxonomy gate. Check taxonomy_rejected.csv in the output directory.",
      call. = FALSE
    )
  }

  tibble::as_tibble(out)
}

#' Apply the marine species-level GRIIS gate
#'
#' Attaches GRIIS evidence to a marine species table using species-only matching
#' (`require_country = FALSE`). This is appropriate for marine overlays because
#' recipient units are polygons rather than ISO country records. Optional filters
#' can retain all rows, only GRIIS-listed rows, or only rows flagged invasive.
#'
#' @param df Input species table.
#' @param output_dir Directory for optional GRIIS audit output.
#' @param use_griis_status Logical; attach GRIIS status when `TRUE`.
#' @param griis Optional standardised GRIIS table.
#' @param species_col Name of the species column in `df`.
#' @param griis_filter_mode One of `"audit_only"`, `"listed"` or `"invasive"`.
#' @param export_griis_audit Logical; write `griis_status_audit.csv`.
#' @param quiet Logical; suppress progress messages.
#'
#' @return Annotated or filtered tibble.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_apply_griis_gate <- function(df,
                                             output_dir,
                                             use_griis_status = FALSE,
                                             griis = NULL,
                                             species_col = "species",
                                             griis_filter_mode = c("audit_only", "listed", "invasive"),
                                             export_griis_audit = TRUE,
                                             quiet = FALSE) {
  griis_filter_mode <- match.arg(griis_filter_mode)

  df <- bf_repair_utf8_df(df)
  griis <- bf_repair_utf8_df(griis)

  if (!isTRUE(use_griis_status)) return(tibble::as_tibble(df))

  if (!exists("bf_attach_griis_status", mode = "function", inherits = TRUE)) {
    stop(
      "`use_griis_status = TRUE`, but `bf_attach_griis_status()` is not available. Add/load utils_griis.R first.",
      call. = FALSE
    )
  }

  out <- bf_attach_griis_status(
    df = df,
    species_col = species_col,
    griis = griis,
    require_country = FALSE,
    quiet = quiet
  )

  out <- bf_repair_utf8_df(out)

  if (isTRUE(export_griis_audit)) {
    readr::write_csv(out, file.path(output_dir, "griis_status_audit.csv"))
  }

  if (griis_filter_mode == "listed") {
    out <- out |> dplyr::filter(.data$griis_listed %in% TRUE)
  } else if (griis_filter_mode == "invasive") {
    out <- out |> dplyr::filter(.data$griis_invasive %in% TRUE)
  }

  out <- bf_repair_utf8_df(out)

  if (nrow(out) == 0L) {
    stop(
      "No marine rows remained after the GRIIS filter. Use `griis_filter_mode = 'audit_only'` first.",
      call. = FALSE
    )
  }

  tibble::as_tibble(out)
}

#' Apply the marine species-level native-origin gate
#'
#' Attaches native-origin evidence using species-only matching. When GRIIS
#' evidence is also present, the function can add reconciled origin-evidence
#' status columns. For marine workflows these columns represent species-level
#' origin evidence unless a separate recipient-specific workflow is supplied.
#'
#' @param df Input species table.
#' @param output_dir Directory for optional native-origin audit output.
#' @param use_native_status Logical; attach native-origin evidence when `TRUE`.
#' @param native_ranges Native-origin table or lookup.
#' @param species_col Species column in `df`.
#' @param native_species_col Species column in `native_ranges`.
#' @param native_filter_mode One of `"audit_only"` or `"origin_available_only"`.
#' @param reconcile_griis_native Logical; add reconciled evidence status when
#'   both GRIIS and native-origin columns are present.
#' @param export_native_audit Logical; write `native_status_audit.csv`.
#' @param quiet Logical; suppress progress messages.
#'
#' @return Annotated or filtered tibble.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_apply_native_gate <- function(df,
                                              output_dir,
                                              use_native_status = FALSE,
                                              native_ranges = NULL,
                                              species_col = "species",
                                              native_species_col = "species",
                                              native_filter_mode = c("audit_only", "origin_available_only"),
                                              reconcile_griis_native = TRUE,
                                              export_native_audit = TRUE,
                                              quiet = FALSE) {
  native_filter_mode <- match.arg(native_filter_mode)

  df <- bf_repair_utf8_df(df)
  native_ranges <- bf_repair_utf8_df(native_ranges)

  if (!isTRUE(use_native_status)) return(tibble::as_tibble(df))

  if (!exists("bf_attach_native_status", mode = "function", inherits = TRUE)) {
    stop(
      "`use_native_status = TRUE`, but `bf_attach_native_status()` is not available. Add/load utils_native_range.R first.",
      call. = FALSE
    )
  }

  out <- bf_attach_native_status(
    df = df,
    species_col = species_col,
    native_ranges = native_ranges,
    native_species_col = native_species_col,
    require_country = FALSE,
    strategy = "tiered",
    native_filter_mode = "audit_only",
    quiet = quiet
  )

  out <- bf_repair_utf8_df(out)

  has_griis <- all(c("griis_listed", "griis_invasive") %in% names(out))
  if (isTRUE(reconcile_griis_native) && has_griis) {
    if (!exists("bf_reconcile_griis_native_status", mode = "function", inherits = TRUE)) {
      stop(
        "`reconcile_griis_native = TRUE`, but `bf_reconcile_griis_native_status()` is not available.",
        call. = FALSE
      )
    }
    out <- bf_reconcile_griis_native_status(out)
    out <- bf_repair_utf8_df(out)
  }

  if (isTRUE(export_native_audit)) {
    readr::write_csv(out, file.path(output_dir, "native_status_audit.csv"))
  }

  if (native_filter_mode == "origin_available_only") {
    out <- out |> dplyr::filter(.data$native_has_origin %in% TRUE)
  }

  out <- bf_repair_utf8_df(out)

  if (nrow(out) == 0L) {
    stop(
      "No marine rows remained after the native-origin filter. Use `native_filter_mode = 'audit_only'` first.",
      call. = FALSE
    )
  }

  tibble::as_tibble(out)
}

#' Build species-level evidence metadata for marine outputs
#'
#' Extracts one row per species from taxonomy, GRIIS and native-origin evidence
#' columns so that those fields can be appended back to exported occurrence rows
#' after GBIF import and overlay assignment.
#'
#' @param df Annotated species table.
#' @param species_col Name of the species column.
#'
#' @return Tibble with one row per species and available evidence columns.
#'
#' @keywords internal
#' @noRd
.bf_marine_pipe_species_metadata <- function(df, species_col = "species") {
  df <- bf_repair_utf8_df(df)

  meta_cols <- names(df)[
    grepl("^(griis_|native_|biofetchr_origin_evidence_status$)", names(df)) |
      names(df) %in% c("raw_name", "fixed_name", "cleaned_name", "was_manual_fixed")
  ]

  if (!length(meta_cols)) {
    return(tibble::tibble(species = unique(as.character(df[[species_col]]))))
  }

  out <- df |>
    dplyr::transmute(
      species = as.character(.data[[species_col]]),
      dplyr::across(dplyr::all_of(meta_cols))
    ) |>
    dplyr::distinct(.data$species, .keep_all = TRUE)

  out <- bf_repair_utf8_df(out)
  out
}

#' Process marine GBIF occurrences using package-managed marine overlays
#'
#' `process_gbif_marine_pipeline()` is the main marine processing workflow in
#' biofetchR. It downloads/imports GBIF occurrence records for marine species,
#' optionally applies pre-download evidence gates, assigns records to a selected
#' marine spatial overlay, optionally cleans and thins the occurrence records,
#' and exports per species x marine-region occurrence tables.
#'
#' @details
#' The pipeline is deliberately ordered so that problematic or unsupported rows
#' are removed before live GBIF downloads are submitted:
#'
#' 1. The input species table is validated and deduplicated.
#' 2. If requested, the taxonomy gate applies manual fixes, cleans names and
#'    rejects non-taxa, open names and higher-rank placeholders.
#' 3. If requested, the GRIIS gate attaches species-level invasive-status
#'    evidence and can filter to listed or invasive species.
#' 4. If requested, the native-origin gate attaches species-level origin evidence
#'    and can filter to species with origin evidence.
#' 5. **Outer batching.** Retained species are grouped according to `batch_size`
#'    for progress and processing organisation.
#' 6. **Sequential GBIF retrieval within each batch.** One species is submitted
#'    to the configured GBIF backend, polled and imported completely before the
#'    next species in that same outer batch is submitted.
#' 7. **Retrieval reconciliation.** Each species receives an explicit GBIF
#'    retrieval outcome before downstream marine processing begins.
#' 8. **Point and overlay preparation.** Successful imports are converted to
#'    WGS84 point `sf` objects and a package-managed marine overlay is loaded once,
#'    or a supplied `overlay_sf` object is used.
#' 9. **Marine-region assignment.** Occurrences are assigned by
#'    `sf::st_within()`, with an `sf::st_intersects()` fallback for coastal or
#'    boundary records.
#' 10. **Cleaning and thinning.** Export groups are optionally cleaned and
#'     spatially thinned.
#' 11. **Export and auditing.** Results are written to CSV and/or returned as a
#'     combined in-memory table, with GBIF retrieval outcomes retained separately
#'     from downstream spatial-processing summaries.
#'
#'
#' @section GBIF download sequencing and retrieval auditing:
#' `batch_size` controls only the outer grouping of species. It does not define
#' the number of simultaneous asynchronous GBIF downloads.
#'
#' Within every outer batch, the pipeline deliberately uses:
#'
#' `submit species -> wait/poll -> import -> reconcile -> next species`.
#'
#' This ordering prevents a nominal batch containing several species from
#' submitting all of those GBIF jobs before any one has completed.
#'
#' Each retained species receives one row in `gbif_retrieval_audit.csv`. The
#' audit contains `species`, `request_id`, `request_type`, `gbif_key`,
#' `gbif_status`, `retrieval_status`, `n_records`, and `fail_reason`.
#' For the marine workflow, `request_type` is `"species"` and the request
#' identifier represents the species-level GBIF retrieval.
#'
#' Successfully completed downloads are recorded as `success` or `success_zero`.
#' Unresolved or failed states are preserved explicitly rather than being
#' interpreted as zero occurrences.
#'
#' Marine GRIIS and native-origin joins are species-level by default because
#' marine recipient units are usually polygons rather than ISO country records.
#' These columns should therefore be interpreted as species-level evidence unless
#' a separate recipient-specific origin table is supplied elsewhere in the
#' workflow.
#'
#' Supported overlay names are provided by `bf_marine_overlay_sources()` and
#' loaded by `bf_load_marine_regions_overlay()`. Typical retained overlays
#' include EEZs, Marine Ecoregions of the World, Large Marine Ecosystems, IHO sea
#' areas, high seas, territorial seas, internal waters, archipelagic waters,
#' Longhurst provinces and UNESCO marine World Heritage sites, depending on the
#' installed Marine Regions/mregions2 data products.
#'
#' The `deps` argument is for testing and advanced use. It can inject mock or
#' alternative implementations of summary, download, import, export or thinning
#' functions without changing the production pipeline.
#'
#' @section Output files:
#' When `store_in_memory = FALSE` or when export helpers are active, the pipeline
#' writes one CSV per retained species x marine-region group. The output
#' directory can also contain `gbif_summary.csv`, `gbif_retrieval_audit.csv`,
#' `gbif_unresolved_requests.csv`, and optional audit files for the taxonomy,
#' GRIIS, native-origin and native-web evidence gates. The retrieval audit records
#' one explicit GBIF outcome per retained species. Successful, empty, failed and
#' unresolved species therefore remain distinguishable independently of the
#' number of marine-region output groups produced downstream. Unresolved requests
#' retain their GBIF keys where available for targeted review or recovery.
#'
#' @section Data access, licensing and attribution:
#' This function can submit live GBIF downloads and can use third-party marine
#' overlay and native-evidence resources. biofetchR does not redistribute those
#' provider datasets from this script; it caches or exports user-side results
#' created during the workflow. Users should retain GBIF download keys/DOIs and
#' cite GBIF, Marine Regions or other overlay providers, GRIIS/native-evidence
#' sources and any supplied spatial layers according to the licence and citation
#' requirements of the exact data products used.
#'
#' @section Failure handling:
#' Download submission, polling, import, sf-conversion, overlay matching and export
#' failures are recorded explicitly where possible.
#'
#' GBIF retrieval outcomes are maintained independently in
#' `gbif_retrieval_audit.csv`. Unresolved or failed statuses copied to
#' `gbif_unresolved_requests.csv` include `taxon_key_failed`, `submit_failed`,
#' `submit_no_key`, `submit_timeout`, `invalid_key`, `pending_timeout`,
#' `download_failed`, `killed`, `cancelled`, `file_erased`, `import_failed`,
#' `split_failed`, and `internal_unreconciled`.
#'
#' `pending_timeout` means that the submitted GBIF download did not reach a
#' terminal state within the configured polling window; its key is retained and
#' the request is not interpreted as a zero-occurrence result.
#'
#' A final species-level completeness gate checks that every retained marine
#' species received a retrieval outcome. Any missing species is recorded as
#' `internal_unreconciled` rather than silently disappearing.
#'
#' Fatal configuration errors, such as missing required packages or unsupported
#' overlay names, still stop early before avoidable GBIF requests are submitted.
#'
#' @param df Data frame containing at least a `species` column. When
#'   `prepare_taxonomy = TRUE`, names are cleaned before GBIF submission and the
#'   cleaned accepted names are used downstream.
#' @param output_dir Directory for per-region CSV exports, `gbif_summary.csv`,
#'   `gbif_retrieval_audit.csv`, `gbif_unresolved_requests.csv`, and optional
#'   taxonomy/GRIIS/native-origin audit files.
#' @param user GBIF username.
#' @param pwd GBIF password.
#' @param email GBIF account email address.
#' @param batch_size Number of species grouped into each outer processing batch.
#'   GBIF submission and import are performed species-by-species within each
#'   batch so this value no longer determines the number of concurrent GBIF jobs.
#' @param apply_cleaning Logical; if `TRUE`, apply coordinate/quality cleaning
#'   through `thin_spatial_points(..., dist_km = 0)` when available.
#' @param apply_thinning Logical; if `TRUE`, spatially thin records after overlay
#'   assignment using `thin_spatial_points()`.
#' @param dist_km Numeric thinning distance in kilometres.
#' @param return_all_results Logical; if `TRUE` and `store_in_memory = TRUE`,
#'   return a combined table of processed occurrence rows.
#' @param export_summary Logical; if `TRUE`, write `gbif_summary.csv`,
#'   `gbif_retrieval_audit.csv`, and `gbif_unresolved_requests.csv`.
#' @param store_in_memory Logical; if `TRUE`, retain processed groups in memory;
#'   if `FALSE`, write groups to CSV and return invisibly.
#' @param use_planar Logical; if `TRUE`, disable spherical S2 predicates during
#'   this run for compatibility with planar overlay operations.
#' @param add_status Logical retained for backwards compatibility. When `TRUE`,
#'   an empty `occurrenceStatus` column is created before cleaning if absent.
#' @param overlay Marine overlay name or alias, resolved by
#'   `bf_marine_overlay_sources()` / `bf_load_marine_regions_overlay()`.
#' @param overlay_sf Optional preloaded `sf` polygon layer. Normally `NULL`, in
#'   which case the package-managed overlay loader is used.
#' @param overlay_cache_dir Cache directory for package-managed marine overlay
#'   downloads and processed overlay objects. If `NULL`, a cache folder is
#'   created under the explicitly supplied `output_dir`.
#' @param overlay_force_refresh Logical; rebuild or re-download overlay resources
#'   where supported.
#' @param strict_overlay_loading Logical; if `TRUE`, stop when the requested
#'   overlay cannot be loaded.
#' @param prepare_taxonomy Logical; if `TRUE`, run the taxonomy gate before GBIF
#'   download submission.
#' @param manual_taxonomy_fixes Optional named character vector or data frame of
#'   manual taxonomic corrections applied before name cleaning.
#' @param taxonomy_name_col Name of the input column containing taxon names.
#' @param export_taxonomy_audit Logical; if `TRUE`, write `taxonomy_audit.csv`,
#'   `taxonomy_rejected.csv` and `taxonomy_summary.csv` when taxonomy preparation
#'   is enabled.
#' @param use_griis_status Logical; if `TRUE`, attach species-level GRIIS status
#'   before GBIF submission.
#' @param griis Optional GRIIS table from `bf_read_griis()` or
#'   `bf_standardise_griis()`.
#' @param griis_filter_mode Character; one of `"audit_only"`, `"listed"` or
#'   `"invasive"`.
#' @param export_griis_audit Logical; if `TRUE`, write `griis_status_audit.csv`.
#' @param use_native_status Logical; if `TRUE`, attach species-level native-origin
#'   evidence before GBIF submission.
#' @param native_ranges Native-origin table or lookup accepted by
#'   `bf_attach_native_status()`. Leave `NULL` when `use_native_web = TRUE` so
#'   the pipeline can build species-level origin evidence from web/API sources
#'   after taxonomy and GRIIS filtering.
#' @param use_native_web Logical; if `TRUE`, fetch species-level native-origin
#'   evidence from package web/API helpers before applying the marine native
#'   evidence gate. Requires `use_native_status = TRUE` and `native_ranges = NULL`.
#' @param native_web_sources Character vector of native web sources passed to
#'   `bf_fetch_native_ranges_web()`, commonly `c("sinas", "gbif", "worms")`.
#' @param native_web_cache_dir Cache directory for native web/API responses. If
#'   `NULL`, a cache folder is created under the explicitly supplied `output_dir`.
#' @param native_web_force_refresh Logical; if `TRUE`, refresh native-web cache
#'   entries where supported.
#' @param native_web_sleep_sec Numeric delay in seconds between native-web
#'   requests.
#' @param native_web_sinas_main_path Optional local path to `SInAS_3.2.csv`
#'   when `"sinas"` is included in `native_web_sources`.
#' @param native_web_sinas_alllocations_path Optional local path to
#'   `AllLocations.xlsx`, `.csv` or `.tsv`.
#' @param native_web_sinas_fulltaxa_path Optional local path to
#'   `SInAS_3.2_FullTaxaList.csv`.
#' @param export_native_web_audit Logical; if `TRUE`, write native-web long,
#'   species, unmapped and summary audit files to `output_dir`.
#' @param native_species_col Species column in `native_ranges` when raw native
#'   range tables are supplied.
#' @param native_filter_mode Character; one of `"audit_only"` or
#'   `"origin_available_only"`.
#' @param reconcile_griis_native Logical; if `TRUE`, add reconciled evidence
#'   status when both GRIIS and native-origin columns are present.
#' @param export_native_audit Logical; if `TRUE`, write `native_status_audit.csv`.
#' @param quiet Logical; suppress console progress messages.
#' @param deps Optional named list of dependency overrides for tests. Supported
#'   entries include `initialize_summary`, `append_summary_row`, `download_fun`,
#'   `import_fun` and `thin_fun`.
#'
#' @return If `return_all_results = TRUE` and `store_in_memory = TRUE`, returns
#'   a tibble containing the combined processed marine occurrence records from
#'   all retained species x marine-region groups. The returned table includes
#'   occurrence columns imported from GBIF, restored `decimalLongitude` and
#'   `decimalLatitude` columns, marine-region assignment fields such as
#'   `region_id`, `region_name`, `region_type`, `workflow`,
#'   `marine_region_id`, `marine_region_name`, `marine_region_source`,
#'   `marine_regions_layer` and `geoname` where available, plus optional
#'   taxonomy, GRIIS and native-origin audit columns when those evidence gates
#'   are enabled. If `return_all_results = FALSE` or `store_in_memory = FALSE`,
#'   the function writes per-group CSV files and summary/audit outputs to
#'   `output_dir` and returns `invisible(NULL)`. In all modes, the primary side
#'   effects are GBIF download/import operations, marine overlay assignment,
#'   optional cleaning/thinning and writing workflow outputs.
#'
#' @examples
#' \donttest{
#' gbif_env <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
#'
#' if (interactive() && all(nzchar(Sys.getenv(gbif_env)))) {
#'   marine_input <- data.frame(
#'     species = c("Ficopomatus enigmaticus", "Carcinus maenas")
#'   )
#'
#'   process_gbif_marine_pipeline(
#'     df = marine_input,
#'     output_dir = file.path(tempdir(), "marine_gbif_outputs"),
#'     user = Sys.getenv("GBIF_USER"),
#'     pwd = Sys.getenv("GBIF_PWD"),
#'     email = Sys.getenv("GBIF_EMAIL"),
#'     overlay = "meow",
#'     prepare_taxonomy = TRUE,
#'     apply_cleaning = TRUE,
#'     apply_thinning = FALSE
#'   )
#' }
#' }
#'
#' @seealso
#' [download_gbif_batch()], [wait_and_import_gbif_safe()],
#' [thin_spatial_points()], [bf_load_marine_regions_overlay()]
#'
#' @family marine processing pipelines
#' @md
#' @export
process_gbif_marine_pipeline <- function(
    df,
    output_dir,
    user = Sys.getenv("GBIF_USER"),
    pwd = Sys.getenv("GBIF_PWD"),
    email = Sys.getenv("GBIF_EMAIL"),
    batch_size = 5,
    apply_cleaning = TRUE,
    apply_thinning = FALSE,
    dist_km = 5,
    return_all_results = TRUE,
    export_summary = TRUE,
    store_in_memory = TRUE,
    use_planar = TRUE,
    add_status = FALSE,
    overlay = "eez",
    overlay_sf = NULL,
    overlay_cache_dir = NULL,
    overlay_force_refresh = FALSE,
    strict_overlay_loading = TRUE,

    prepare_taxonomy = FALSE,
    manual_taxonomy_fixes = NULL,
    taxonomy_name_col = "species",
    export_taxonomy_audit = TRUE,

    use_griis_status = FALSE,
    griis = NULL,
    griis_filter_mode = c("audit_only", "listed", "invasive"),
    export_griis_audit = TRUE,

    use_native_status = FALSE,
    native_ranges = NULL,
    use_native_web = FALSE,
    native_web_sources = c("sinas", "gbif", "worms"),
    native_web_cache_dir = NULL,
    native_web_force_refresh = FALSE,
    native_web_sleep_sec = 0.25,
    native_web_sinas_main_path = NULL,
    native_web_sinas_alllocations_path = NULL,
    native_web_sinas_fulltaxa_path = NULL,
    export_native_web_audit = TRUE,
    native_species_col = "species",
    native_filter_mode = c("audit_only", "origin_available_only"),
    reconcile_griis_native = TRUE,
    export_native_audit = TRUE,

    quiet = FALSE,
    deps = NULL
) {
  bf_require_packages(c("sf", "dplyr", "tibble", "readr"), context = "marine pipeline")

  df <- bf_repair_utf8_df(df)
  griis <- bf_repair_utf8_df(griis)
  native_ranges <- bf_repair_utf8_df(native_ranges)
  overlay_sf <- bf_repair_utf8_df(overlay_sf)

  griis_filter_mode <- match.arg(griis_filter_mode)
  native_filter_mode <- match.arg(native_filter_mode)

  if (isTRUE(use_native_web) && !isTRUE(use_native_status)) {
    stop(
      "`use_native_web = TRUE` requires `use_native_status = TRUE`, because web-derived native-origin evidence is only used by the native evidence gate.",
      call. = FALSE
    )
  }

  if (isTRUE(use_native_web) && !is.null(native_ranges)) {
    stop(
      "`use_native_web = TRUE` currently requires `native_ranges = NULL`. Use either supplied native-range evidence or web-derived evidence, not both.",
      call. = FALSE
    )
  }

  if (!"species" %in% names(df)) {
    stop("`df` must contain a `species` column.", call. = FALSE)
  }

  if (missing(output_dir) || is.null(output_dir) || length(output_dir) == 0L ||
      !nzchar(trimws(as.character(output_dir[[1L]])))) {
    stop("`output_dir` must be supplied explicitly.", call. = FALSE)
  }

  output_dir <- normalizePath(
    as.character(output_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(overlay_cache_dir) || length(overlay_cache_dir) == 0L ||
      !nzchar(trimws(as.character(overlay_cache_dir[[1L]])))) {
    overlay_cache_dir <- file.path(
      output_dir,
      "_biofetchR_cache",
      "marine_overlays"
    )
  }

  overlay_cache_dir <- normalizePath(
    as.character(overlay_cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (is.null(native_web_cache_dir) || length(native_web_cache_dir) == 0L ||
      !nzchar(trimws(as.character(native_web_cache_dir[[1L]])))) {
    native_web_cache_dir <- file.path(
      output_dir,
      "_biofetchR_cache",
      "native_web"
    )
  }

  native_web_cache_dir <- normalizePath(
    as.character(native_web_cache_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (!exists("bf_marine_overlay_sources", mode = "function")) {
    stop("bf_marine_overlay_sources() is not available. Source/load utils_marine_context_overlays.R first.", call. = FALSE)
  }
  if (!exists("bf_load_marine_regions_overlay", mode = "function")) {
    stop("bf_load_marine_regions_overlay() is not available. Source/load utils_marine_context_overlays.R first.", call. = FALSE)
  }

  overlay <- .bf_marine_pipe_norm(overlay)
  supported <- bf_marine_overlay_sources()
  if (!overlay %in% supported) {
    stop(
      "Unsupported marine overlay: '", overlay, "'. Supported overlays are: ",
      paste(supported, collapse = ", "),
      call. = FALSE
    )
  }

  .dep <- function(key, fallback = NULL) {
    if (!is.null(deps) && is.list(deps) && !is.null(deps[[key]]) && is.function(deps[[key]])) return(deps[[key]])
    fallback
  }

  initialize_summary_fun <- .dep("initialize_summary", if (exists("initialize_summary", mode = "function")) initialize_summary else .bf_marine_pipe_default_summary)
  append_summary_fun     <- .dep("append_summary_row", if (exists("append_summary_row", mode = "function")) append_summary_row else .bf_marine_pipe_default_append)
  download_fun           <- .dep("download_fun", if (exists("download_gbif_batch", mode = "function")) download_gbif_batch else NULL)
  import_fun             <- .dep("import_fun", if (exists("wait_and_import_gbif_safe", mode = "function")) wait_and_import_gbif_safe else if (exists("wait_and_import_gbif", mode = "function")) wait_and_import_gbif else NULL)
  thin_fun               <- .dep("thin_fun", if (exists("thin_spatial_points", mode = "function")) thin_spatial_points else NULL)

  if (is.null(download_fun)) {
    stop("No GBIF download function is available. Provide deps$download_fun or define download_gbif_batch().", call. = FALSE)
  }
  if (is.null(import_fun)) {
    stop("No GBIF import function is available. Provide deps$import_fun or define wait_and_import_gbif[_safe]().", call. = FALSE)
  }

  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(!isTRUE(use_planar))

  df <- df |>
    dplyr::mutate(species = trimws(as.character(.data$species))) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species)) |>
    dplyr::distinct(.data$species, .keep_all = TRUE)

  df <- bf_repair_utf8_df(df)

  df <- .bf_marine_pipe_apply_taxonomy_gate(
    df = df,
    output_dir = output_dir,
    prepare_taxonomy = prepare_taxonomy,
    manual_taxonomy_fixes = manual_taxonomy_fixes,
    taxonomy_name_col = taxonomy_name_col,
    export_taxonomy_audit = export_taxonomy_audit,
    quiet = quiet
  ) |>
    dplyr::mutate(species = trimws(as.character(.data[[taxonomy_name_col]]))) |>
    dplyr::filter(!is.na(.data$species), nzchar(.data$species)) |>
    dplyr::distinct(.data$species, .keep_all = TRUE)

  df <- bf_repair_utf8_df(df)

  df <- .bf_marine_pipe_apply_griis_gate(
    df = df,
    output_dir = output_dir,
    use_griis_status = use_griis_status,
    griis = griis,
    species_col = "species",
    griis_filter_mode = griis_filter_mode,
    export_griis_audit = export_griis_audit,
    quiet = quiet
  ) |>
    dplyr::distinct(.data$species, .keep_all = TRUE)

  df <- bf_repair_utf8_df(df)

  if (isTRUE(use_native_web)) {
    if (!exists("bf_fetch_native_ranges_web", mode = "function", inherits = TRUE)) {
      stop(
        "`use_native_web = TRUE`, but `bf_fetch_native_ranges_web()` is not available. ",
        "Make sure `R/utils_native_range_web_sources.R` is present and loaded by the package.",
        call. = FALSE
      )
    }

    species_for_native_web <- unique(as.character(df$species))
    species_for_native_web <- species_for_native_web[!is.na(species_for_native_web) & nzchar(species_for_native_web)]

    native_web <- bf_fetch_native_ranges_web(
      species = species_for_native_web,
      sources = native_web_sources,
      cache_dir = native_web_cache_dir,
      force_refresh = native_web_force_refresh,
      sleep_sec = native_web_sleep_sec,
      quiet = quiet,
      sinas_main_path = native_web_sinas_main_path,
      sinas_alllocations_path = native_web_sinas_alllocations_path,
      sinas_fulltaxa_path = native_web_sinas_fulltaxa_path,
      return = "list"
    )

    if (isTRUE(export_native_web_audit)) {
      if (!exists("bf_write_native_web_outputs", mode = "function", inherits = TRUE)) {
        stop(
          "`export_native_web_audit = TRUE`, but `bf_write_native_web_outputs()` is not available.",
          call. = FALSE
        )
      }

      bf_write_native_web_outputs(
        native_web,
        output_dir = output_dir,
        prefix = "native_web"
      )
    }

    native_ranges <- bf_repair_utf8_df(native_web$species)
    native_species_col <- "species"

    if (!isTRUE(quiet)) {
      message(
        "biofetchR native web evidence: fetched species-level origin evidence for ",
        length(species_for_native_web),
        " marine species."
      )
    }
  }

  df <- .bf_marine_pipe_apply_native_gate(
    df = df,
    output_dir = output_dir,
    use_native_status = use_native_status,
    native_ranges = native_ranges,
    species_col = "species",
    native_species_col = native_species_col,
    native_filter_mode = native_filter_mode,
    reconcile_griis_native = reconcile_griis_native,
    export_native_audit = export_native_audit,
    quiet = quiet
  ) |>
    dplyr::distinct(.data$species, .keep_all = TRUE)

  # Repair encodings before species metadata, download submission and joins.
  # This protects the pipeline against invalid strings introduced through
  # taxonomy, GRIIS or native-range metadata.
  df <- bf_repair_utf8_df(df)

  species_metadata <- .bf_marine_pipe_species_metadata(df, species_col = "species")
  species_metadata <- bf_repair_utf8_df(species_metadata)

  if (!nrow(df)) {
    if (isTRUE(export_summary)) {
      empty_summary <- .bf_marine_pipe_default_summary()
      empty_summary <- bf_repair_utf8_df(empty_summary)
      readr::write_csv(
        empty_summary,
        file.path(output_dir, "gbif_summary.csv")
      )
    }

    if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
      return(tibble::tibble(decimalLongitude = numeric(), decimalLatitude = numeric()))
    }

    return(invisible(NULL))
  }

  species_list <- unique(df$species)
  batches <- split(species_list, ceiling(seq_along(species_list) / batch_size))

  pipeline_start_time <- bf_console_now()

  bf_console_pipeline_start(
    workflow = "marine",
    n_species = length(species_list),
    region_source = overlay,
    gadm_unit = NULL,
    overlay_mode = "single",
    cleaning = apply_cleaning,
    thinning = apply_thinning,
    quiet = quiet
  )

  bf_console_bullet(
    paste0("Loading marine overlay resource: ", toupper(bf_repair_utf8_one(overlay))),
    quiet = quiet
  )

  overlay_load_start <- bf_console_now()

  # ---------------------------------------------------------------------------
  # Load overlay once.
  # ---------------------------------------------------------------------------
  marine_overlay <- tryCatch(
    {
      if (!is.null(overlay_sf)) {
        .bf_marine_pipe_ensure_sf(overlay_sf)
      } else {
        bf_load_marine_regions_overlay(
          region_source = overlay,
          cache_dir = overlay_cache_dir,
          force_refresh = overlay_force_refresh,
          quiet = TRUE
        )
      }
    },
    error = function(e) {
      if (isTRUE(strict_overlay_loading)) {
        stop("Requested marine overlay failed to load: ", overlay, " - ", bf_compact_error(e), call. = FALSE)
      }
      NULL
    }
  )

  if (is.null(marine_overlay) || !inherits(marine_overlay, "sf") || !nrow(marine_overlay)) {
    stop("Marine overlay '", overlay, "' is unavailable or empty.", call. = FALSE)
  }

  marine_overlay <- .bf_marine_pipe_make_valid(.bf_marine_pipe_ensure_sf(marine_overlay))
  marine_overlay <- bf_repair_utf8_df(marine_overlay)

  if (!"geoname" %in% names(marine_overlay)) marine_overlay$geoname <- as.character(marine_overlay$marine_region_name)
  if (!"marine_region_id" %in% names(marine_overlay)) marine_overlay$marine_region_id <- as.character(seq_len(nrow(marine_overlay)))
  if (!"marine_region_name" %in% names(marine_overlay)) marine_overlay$marine_region_name <- as.character(marine_overlay$geoname)
  if (!"marine_region_source" %in% names(marine_overlay)) marine_overlay$marine_region_source <- overlay
  if (!"marine_regions_layer" %in% names(marine_overlay)) marine_overlay$marine_regions_layer <- overlay

  marine_overlay <- bf_repair_utf8_df(marine_overlay)

  bf_console_bullet(
    paste0(
      "Marine overlay loaded: ",
      toupper(bf_repair_utf8_one(overlay)),
      " | ",
      nrow(marine_overlay),
      " feature(s) | ",
      bf_console_elapsed(overlay_load_start),
      "."
    ),
    quiet = quiet
  )

  summary_tbl <- initialize_summary_fun()
  summary_tbl <- bf_repair_utf8_df(summary_tbl)

  for (nm in c("fail_stage", "fail_reason", "gbif_key")) {
    if (!nm %in% names(summary_tbl)) summary_tbl[[nm]] <- NA_character_
  }

  summary_tbl <- bf_repair_utf8_df(summary_tbl)

  retrieval_audit <- tibble::tibble(
    species = character(),
    request_id = character(),
    request_type = character(),
    gbif_key = character(),
    gbif_status = character(),
    retrieval_status = character(),
    n_records = integer(),
    fail_reason = character()
  )

  all_results <- list()
  result_i <- 0L

  append_row <- function(species, region_id, region_type,
                         n_total = 0, n_cleaned = NA_integer_, n_thinned = 0,
                         output_file = NA_character_, status = "unknown",
                         fail_stage = NA_character_, fail_reason = NA_character_, gbif_key = NA_character_) {
    species <- bf_repair_utf8_one(species)
    region_id <- bf_repair_utf8_one(region_id)
    region_type <- bf_repair_utf8_one(region_type)
    output_file <- bf_repair_utf8_one(output_file, fallback = NA_character_)
    status <- bf_repair_utf8_one(status)
    fail_stage <- bf_repair_utf8_one(fail_stage, fallback = NA_character_)
    fail_reason <- bf_repair_utf8_one(fail_reason, fallback = NA_character_)
    gbif_key <- bf_repair_utf8_one(gbif_key, fallback = NA_character_)

    summary_tbl <<- append_summary_fun(
      summary_tbl,
      species = species,
      region_id = region_id,
      region_type = region_type,
      n_total = n_total,
      n_cleaned = n_cleaned,
      n_thinned = n_thinned,
      output_file = output_file,
      status = status
    )

    summary_tbl <<- bf_repair_utf8_df(summary_tbl)

    i <- nrow(summary_tbl)
    summary_tbl$fail_stage[i] <<- fail_stage
    summary_tbl$fail_reason[i] <<- fail_reason
    summary_tbl$gbif_key[i] <<- gbif_key

    summary_tbl <<- bf_repair_utf8_df(summary_tbl)

    invisible(NULL)
  }

  append_retrieval_audit <- function(species,
                                   gbif_key = NA_character_,
                                   gbif_status = NA_character_,
                                   retrieval_status,
                                   n_records = NA_integer_,
                                   fail_reason = NA_character_) {
    species <- bf_repair_utf8_one(species)

    if (nrow(retrieval_audit)) {
      keep <- retrieval_audit$species != species
      retrieval_audit <<- retrieval_audit[keep, , drop = FALSE]
    }

    retrieval_audit <<- dplyr::bind_rows(
      retrieval_audit,
      tibble::tibble(
        species = species,
        request_id = species,
        request_type = "species",
        gbif_key = bf_repair_utf8_one(gbif_key, fallback = NA_character_),
        gbif_status = bf_repair_utf8_one(gbif_status, fallback = NA_character_),
        retrieval_status = bf_repair_utf8_one(retrieval_status),
        n_records = as.integer(n_records[[1L]]),
        fail_reason = bf_repair_utf8_one(fail_reason, fallback = NA_character_)
      )
    )

    retrieval_audit <<- bf_repair_utf8_df(retrieval_audit)
    invisible(NULL)
  }

  clean_and_thin <- function(x) {
    x <- bf_repair_utf8_df(x)

    n_total <- nrow(x)
    cleaned <- x
    n_cleaned <- n_total

    if (isTRUE(add_status) && !"occurrenceStatus" %in% names(cleaned)) {
      cleaned$occurrenceStatus <- NA_character_
    }

    cleaned <- bf_repair_utf8_df(cleaned)

    if (isTRUE(apply_cleaning) && is.function(thin_fun)) {
      tmp <- tryCatch(
        thin_fun(cleaned, dist_km = 0, filter_uncertain = TRUE, quiet = TRUE),
        error = function(e) cleaned
      )
      cleaned <- tmp
      cleaned <- bf_repair_utf8_df(cleaned)
      n_cleaned <- nrow(cleaned)
    }

    final <- cleaned
    final <- bf_repair_utf8_df(final)

    if (isTRUE(apply_thinning) && is.function(thin_fun) && is.finite(dist_km) && dist_km > 0) {
      tmp <- tryCatch(
        thin_fun(final, dist_km = dist_km, filter_uncertain = FALSE, quiet = TRUE),
        error = function(e) final
      )
      final <- tmp
      final <- bf_repair_utf8_df(final)
    }

    list(
      sf = final,
      n_total = n_total,
      n_cleaned = n_cleaned,
      n_thinned = nrow(final)
    )
  }

  export_group <- function(sf_group, species_name, region_id, region_name, gbif_key = NA_character_) {
    if (!nrow(sf_group)) return(invisible(NULL))

    species_name_safe <- bf_repair_utf8_one(species_name)
    overlay_safe <- bf_repair_utf8_one(overlay)
    region_id_safe <- bf_repair_utf8_one(region_id)
    region_name_safe <- bf_repair_utf8_one(region_name)
    gbif_key_safe <- bf_repair_utf8_one(gbif_key, fallback = NA_character_)

    sf_group <- bf_repair_utf8_df(sf_group)
    sf_group <- .bf_marine_pipe_restore_coords(.bf_marine_pipe_ensure_sf(sf_group))
    sf_group <- bf_repair_utf8_df(sf_group)

    proc <- clean_and_thin(sf_group)

    bf_console_bullet(
      paste0(
        "Preparing export group: ",
        species_name_safe,
        " | ",
        toupper(overlay_safe),
        " = ",
        region_id_safe
      ),
      quiet = quiet
    )

    bf_console_cleaning(
      n_before = proc$n_total,
      n_after = proc$n_cleaned,
      quiet = quiet
    )

    bf_console_thinning(
      n_before = proc$n_cleaned,
      n_after = proc$n_thinned,
      applied = apply_thinning,
      quiet = quiet
    )

    proc$sf <- .bf_marine_pipe_restore_coords(proc$sf)
    proc$sf <- bf_repair_utf8_df(proc$sf)

    proc$sf$species_region <- paste0(species_name_safe, "__", region_id_safe)
    proc$sf$region_id <- as.character(region_id_safe)
    proc$sf$region_name <- as.character(region_name_safe)
    proc$sf$region_type <- toupper(overlay_safe)
    proc$sf$workflow <- paste0("marine_", overlay_safe)

    proc$sf <- bf_repair_utf8_df(proc$sf)

    if (exists("species_metadata", inherits = TRUE) && is.data.frame(species_metadata)) {
      meta <- species_metadata[species_metadata$species %in% species_name_safe, , drop = FALSE]
      meta <- bf_repair_utf8_df(meta)
      if (nrow(meta)) {
        add_cols <- setdiff(names(meta), c("species", names(proc$sf)))
        if (length(add_cols)) {
          for (nm in add_cols) proc$sf[[nm]] <- meta[[nm]][[1]]
        }
      }
    }

    proc$sf <- bf_repair_utf8_df(proc$sf)

    # ---------------------------------------------------------------------------
    # Resolve and write grouped marine export through the universal export contract
    # ---------------------------------------------------------------------------
    # Marine outputs are global GBIF downloads assigned to one marine overlay. The
    # input region is therefore recorded as GLOBAL, while `region_id` records the
    # assigned marine polygon/unit.

    export_stage_value <- paste0("marine_", as.character(overlay_safe))

    export_plan <- bf_resolve_export_plan(
      output_dir = output_dir,
      species = species_name_safe,
      input_region = "GLOBAL",
      export_stage = export_stage_value,
      workflow = as.character(export_stage_value),
      region_type = toupper(as.character(overlay_safe)),
      region_source = as.character(overlay_safe),
      spatial_join_type = as.character(overlay_safe),
      region_id = as.character(region_id_safe),
      split_by_region = TRUE,
      apply_cleaning = apply_cleaning,
      apply_thinning = apply_thinning,
      thinning_dist_km = dist_km,
      use_overlays = FALSE,
      overlay_names = NULL,
      prepare_taxonomy = prepare_taxonomy,
      use_native_range = use_native_status,
      use_native_web = use_native_web
    )

    out_csv <- sf::st_drop_geometry(proc$sf)
    out_csv <- bf_repair_utf8_df(out_csv)

    export_plan <- tryCatch(
      {
        bf_write_export_csv(out_csv, export_plan)
      },
      error = function(e) {
        reason <- bf_compact_error(e)

        bf_console_bullet(
          paste0(
            "Universal export writer failed for ",
            species_name_safe,
            " / ",
            region_id_safe,
            "; falling back to readr::write_csv(). Reason: ",
            reason
          ),
          quiet = quiet
        )

        dir.create(dirname(export_plan$output_file), recursive = TRUE, showWarnings = FALSE)
        readr::write_csv(out_csv, export_plan$output_file)

        export_plan$export_written <- file.exists(export_plan$output_file)
        export_plan$export_n_rows <- nrow(out_csv)
        export_plan$export_n_cols <- ncol(out_csv)
        export_plan$export_created_at <- as.character(Sys.time())

        export_plan
      }
    )

    try(
      bf_append_export_manifest(
        export_plan,
        extra = list(
          overlay = as.character(overlay_safe),
          region_name = as.character(region_name_safe),
          gbif_key = as.character(gbif_key_safe)
        )
      ),
      silent = TRUE
    )

    output_file_path <- if (file.exists(export_plan$output_file)) {
      export_plan$output_file
    } else {
      NA_character_
    }

    out_data <- if (isTRUE(store_in_memory)) {
      proc$sf
    } else {
      output_file_path
    }

    append_row(
      species = species_name_safe,
      region_id = region_id_safe,
      region_type = toupper(overlay_safe),
      n_total = proc$n_total,
      n_cleaned = proc$n_cleaned,
      n_thinned = proc$n_thinned,
      output_file = output_file_path,
      status = "success",
      gbif_key = gbif_key_safe
    )

    if (exists("species_output_rows_console", inherits = TRUE)) {
      species_output_rows_console <<- species_output_rows_console + as.integer(proc$n_thinned)
    }

    if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
      result_i <<- result_i + 1L
      all_results[[result_i]] <<- out_data
    }

    invisible(NULL)
  }

  import_batch <- function(download_keys) {
    import_fun(download_keys)
  }

  submit_batch <- function(batch_df) {
    batch_df <- bf_repair_utf8_df(batch_df)

    # Test backend signature: download_fun(batch, user, pwd, email)
    # Some real backends use named arguments. Try both safely.
    out <- tryCatch(
      download_fun(batch = batch_df, user = user, pwd = pwd, email = email),
      error = function(e1) {
        tryCatch(
          download_fun(batch_df, user, pwd, email),
          error = function(e2) {
            tryCatch(
              download_fun(species = unique(batch_df$species), user = user, pwd = pwd, email = email),
              error = function(e3) stop(e1)
            )
          }
        )
      }
    )
    out
  }

  for (batch_species in batches) {
    batch_species <- bf_repair_utf8_chr(batch_species)

    batch_df <- df[df$species %in% batch_species, , drop = FALSE]
    batch_df <- bf_repair_utf8_df(batch_df)

    # `batch_size` still controls progress grouping, but GBIF submission/import is
    # intentionally species-by-species. This prevents a nominal batch of five
    # species from creating five simultaneous asynchronous downloads.
    imported <- list()
    gbif_key_map <- character(0)

    for (sp in batch_species) {
      sp <- bf_repair_utf8_one(sp)
      species_df <- batch_df[batch_df$species == sp, , drop = FALSE]
      species_df <- bf_repair_utf8_df(species_df)

      download_start_time <- bf_console_now()

      bf_console_bullet(
        paste0("Submitting marine GBIF request for ", sp, "."),
        quiet = quiet
      )

      one_keys <- tryCatch(
        submit_batch(species_df),
        error = function(e) e
      )

      if (inherits(one_keys, "condition")) {
        reason <- bf_compact_error(one_keys)
        append_row(
          species = sp,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "failed",
          fail_stage = "submit_download",
          fail_reason = reason
        )
        append_retrieval_audit(
          species = sp,
          retrieval_status = "submit_failed",
          fail_reason = reason
        )
        next
      }

      submission_audit <- attr(one_keys, "submission_audit", exact = TRUE)

      if (is.null(one_keys) || !length(one_keys)) {
        row_sp <- NULL
        if (is.data.frame(submission_audit) && "label" %in% names(submission_audit)) {
          hit <- which(as.character(submission_audit$label) == sp)
          if (length(hit)) row_sp <- submission_audit[hit[[1L]], , drop = FALSE]
        }

        sub_status <- if (!is.null(row_sp) && "submission_status" %in% names(row_sp)) {
          as.character(row_sp$submission_status[[1L]])
        } else {
          "submit_failed"
        }
        sub_reason <- if (!is.null(row_sp) && "message" %in% names(row_sp)) {
          as.character(row_sp$message[[1L]])
        } else {
          "GBIF request returned no usable download key."
        }
        sub_key <- if (!is.null(row_sp) && "gbif_key" %in% names(row_sp)) {
          as.character(row_sp$gbif_key[[1L]])
        } else {
          NA_character_
        }

        append_row(
          species = sp,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "failed",
          fail_stage = sub_status,
          fail_reason = sub_reason,
          gbif_key = sub_key
        )
        append_retrieval_audit(
          species = sp,
          gbif_key = sub_key,
          retrieval_status = sub_status,
          fail_reason = sub_reason
        )
        next
      }

      key_sp <- as.character(unlist(one_keys, use.names = FALSE))
      key_sp <- key_sp[!is.na(key_sp) & nzchar(key_sp)]
      key_sp <- if (length(key_sp)) key_sp[[1L]] else NA_character_

      bf_console_gbif_request(
        species = sp,
        elapsed = bf_console_elapsed(download_start_time),
        gbif_key = key_sp,
        quiet = quiet
      )

      import_start_time <- bf_console_now()
      one_imported <- tryCatch(
        import_batch(one_keys),
        error = function(e) e
      )

      if (inherits(one_imported, "condition")) {
        reason <- bf_compact_error(one_imported)
        append_row(
          species = sp,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "failed",
          fail_stage = "wait_import",
          fail_reason = reason,
          gbif_key = key_sp
        )
        append_retrieval_audit(
          species = sp,
          gbif_key = key_sp,
          retrieval_status = "import_failed",
          fail_reason = reason
        )
        next
      }

      import_audit <- attr(one_imported, "gbif_status", exact = TRUE)
      row_sp <- NULL
      if (is.data.frame(import_audit) && "label" %in% names(import_audit)) {
        hit <- which(as.character(import_audit$label) == sp)
        if (length(hit)) row_sp <- import_audit[hit[[1L]], , drop = FALSE]
      }

      if (!is.null(row_sp)) {
        retrieval_status_sp <- as.character(row_sp$retrieval_status[[1L]])
        gbif_status_sp <- as.character(row_sp$gbif_status[[1L]])
        n_records_sp <- suppressWarnings(as.integer(row_sp$n_records[[1L]]))
        reason_sp <- as.character(row_sp$message[[1L]])
        audit_key_sp <- as.character(row_sp$gbif_key[[1L]])

        append_retrieval_audit(
          species = sp,
          gbif_key = audit_key_sp,
          gbif_status = gbif_status_sp,
          retrieval_status = retrieval_status_sp,
          n_records = n_records_sp,
          fail_reason = reason_sp
        )

        if (!retrieval_status_sp %in% c("success", "success_zero")) {
          append_row(
            species = sp,
            region_id = overlay,
            region_type = toupper(overlay),
            status = "failed",
            fail_stage = retrieval_status_sp,
            fail_reason = reason_sp,
            gbif_key = audit_key_sp
          )
          next
        }
      }

      if (is.null(names(one_imported))) {
        names(one_imported) <- sp
      }

      if (!sp %in% names(one_imported)) {
        reason <- "Species request was neither returned by the importer nor represented as an explicit retrieval failure."
        append_retrieval_audit(
          species = sp,
          gbif_key = key_sp,
          retrieval_status = "internal_unreconciled",
          fail_reason = reason
        )
        append_row(
          species = sp,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "failed",
          fail_stage = "internal_unreconciled",
          fail_reason = reason,
          gbif_key = key_sp
        )
        next
      }

      imported[[sp]] <- bf_repair_utf8_df(one_imported[[sp]])
      gbif_key_map[[sp]] <- key_sp

      if (is.null(row_sp)) {
        n_records_sp <- as.integer(bf_console_nrow(imported[[sp]]))
        append_retrieval_audit(
          species = sp,
          gbif_key = key_sp,
          gbif_status = "SUCCEEDED",
          retrieval_status = if (n_records_sp > 0L) "success" else "success_zero",
          n_records = n_records_sp
        )
      }

      bf_console_gbif_import(
        n_records = as.integer(bf_console_nrow(imported[[sp]])),
        elapsed = bf_console_elapsed(import_start_time),
        quiet = quiet
      )
    }

    if (!length(imported)) {
      bf_console_bullet(
        paste0("GBIF import produced no successfully resolved species in batch: ", paste(batch_species, collapse = ", ")),
        quiet = quiet
      )
      next
    }

    names(imported) <- bf_repair_utf8_chr(names(imported))

    for (species_name in names(imported)) {
      species_name <- bf_repair_utf8_one(species_name)

      species_start_time <- bf_console_now()
      species_output_rows_console <- 0L

      species_index_console <- match(species_name, species_list)

      bf_console_species_start(
        species = species_name,
        index = species_index_console,
        total = length(species_list),
        region = toupper(bf_repair_utf8_one(overlay)),
        quiet = quiet
      )

      raw_sf <- imported[[species_name]]
      raw_sf <- bf_repair_utf8_df(raw_sf)

      gbif_key <- if (!is.null(gbif_key_map[[species_name]])) gbif_key_map[[species_name]] else NA_character_
      gbif_key <- bf_repair_utf8_one(gbif_key, fallback = NA_character_)

      if (is.null(raw_sf) || !nrow(raw_sf)) {
        append_row(
          species = species_name,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "empty",
          fail_stage = "import_empty",
          fail_reason = "GBIF import returned 0 rows.",
          gbif_key = gbif_key
        )

        bf_console_species_end(
          species = species_name,
          n_final = species_output_rows_console,
          start_time = species_start_time,
          quiet = quiet
        )

        next
      }

      pts <- tryCatch(.bf_marine_pipe_restore_coords(.bf_marine_pipe_ensure_sf(raw_sf)), error = function(e) NULL)
      pts <- bf_repair_utf8_df(pts)

      if (is.null(pts) || !nrow(pts)) {
        append_row(
          species = species_name,
          region_id = overlay,
          region_type = toupper(overlay),
          status = "failed",
          fail_stage = "ensure_sf",
          fail_reason = "Failed to coerce imported records to sf EPSG:4326.",
          gbif_key = gbif_key
        )

        bf_console_species_end(
          species = species_name,
          n_final = species_output_rows_console,
          start_time = species_start_time,
          quiet = quiet
        )

        next
      }
      pts$species <- species_name
      pts <- bf_repair_utf8_df(pts)

      bf_console_bullet(
        paste0(
          species_name,
          ": ",
          nrow(pts),
          " georeferenced marine GBIF records imported."
        ),
        quiet = quiet
      )

      overlay_join_start_time <- bf_console_now()

      bf_console_overlay_start(
        region_source = toupper(bf_repair_utf8_one(overlay)),
        gadm_unit = NULL,
        quiet = quiet
      )

      # Spatial containment join. For coastal/boundary points, fall back to intersects.
      joined <- tryCatch(
        suppressWarnings(sf::st_join(pts, marine_overlay, join = sf::st_within, left = FALSE)),
        error = function(e) NULL
      )

      joined <- bf_repair_utf8_df(joined)

      if (is.null(joined) || !nrow(joined)) {
        joined <- tryCatch(
          suppressWarnings(sf::st_join(pts, marine_overlay, join = sf::st_intersects, left = FALSE)),
          error = function(e) NULL
        )

        joined <- bf_repair_utf8_df(joined)
      }

      if (is.null(joined) || !nrow(joined)) {
        bf_console_overlay_result(
          n_records = nrow(pts),
          n_matched = 0L,
          id_col = "geoname",
          quiet = quiet
        )

        bf_console_bullet(
          paste0(
            "Marine overlay join completed in ",
            bf_console_elapsed(overlay_join_start_time),
            "."
          ),
          quiet = quiet
        )

        append_row(
          species = species_name,
          region_id = overlay,
          region_type = toupper(overlay),
          n_total = nrow(pts),
          n_cleaned = nrow(pts),
          n_thinned = 0,
          status = "empty",
          fail_stage = "no_overlay_match",
          fail_reason = "No imported points intersected the requested marine overlay.",
          gbif_key = gbif_key
        )

        bf_console_species_end(
          species = species_name,
          n_final = species_output_rows_console,
          start_time = species_start_time,
          quiet = quiet
        )

        next
      }

      if (!"geoname" %in% names(joined)) {
        joined$geoname <- if ("marine_region_name" %in% names(joined)) joined$marine_region_name else NA_character_
      }
      if (!"marine_region_id" %in% names(joined)) joined$marine_region_id <- joined$geoname
      if (!"marine_region_name" %in% names(joined)) joined$marine_region_name <- joined$geoname

      joined <- bf_repair_utf8_df(joined)

      joined <- joined[!is.na(joined$geoname) & nzchar(as.character(joined$geoname)), , drop = FALSE]
      joined <- bf_repair_utf8_df(joined)

      if (!nrow(joined)) {
        bf_console_overlay_result(
          n_records = nrow(pts),
          n_matched = 0L,
          id_col = "geoname",
          quiet = quiet
        )

        bf_console_bullet(
          paste0(
            "Marine overlay join completed in ",
            bf_console_elapsed(overlay_join_start_time),
            "."
          ),
          quiet = quiet
        )

        bf_console_species_end(
          species = species_name,
          n_final = species_output_rows_console,
          start_time = species_start_time,
          quiet = quiet
        )

        next
      }

      bf_console_overlay_result(
        n_records = nrow(pts),
        n_matched = nrow(joined),
        id_col = "geoname",
        quiet = quiet
      )

      bf_console_bullet(
        paste0(
          "Marine overlay join completed in ",
          bf_console_elapsed(overlay_join_start_time),
          "."
        ),
        quiet = quiet
      )

      split_ids <- as.character(joined$geoname)
      split_ids <- bf_repair_utf8_chr(split_ids)

      groups <- split(joined, split_ids)

      for (rid in names(groups)) {
        rid <- bf_repair_utf8_one(rid)

        region_name <- rid
        if ("marine_region_name" %in% names(groups[[rid]])) {
          region_name <- unique(as.character(groups[[rid]]$marine_region_name))[1]
        }

        region_name <- bf_repair_utf8_one(region_name)

        groups[[rid]] <- bf_repair_utf8_df(groups[[rid]])

        export_group(groups[[rid]], species_name, rid, region_name, gbif_key = gbif_key)
      }

      bf_console_species_end(
        species = species_name,
        n_final = species_output_rows_console,
        start_time = species_start_time,
        quiet = quiet
      )
    }
  }

  if (isTRUE(export_summary)) {
    expected_species <- tibble::tibble(
      species = unique(as.character(df$species))
    )
    observed_species <- retrieval_audit |>
      dplyr::distinct(.data$species)

    missing_species <- dplyr::anti_join(
      expected_species,
      observed_species,
      by = "species"
    )

    if (nrow(missing_species)) {
      for (i in seq_len(nrow(missing_species))) {
        append_retrieval_audit(
          species = missing_species$species[[i]],
          retrieval_status = "internal_unreconciled",
          fail_reason = "Final retrieval completeness gate found no recorded outcome for this species."
        )
      }
      warning(
        nrow(missing_species),
        " marine species lacked a retrieval outcome and were marked internal_unreconciled.",
        call. = FALSE
      )
    }

    unresolved_statuses <- c(
      "taxon_key_failed",
      "submit_failed",
      "submit_no_key",
      "submit_timeout",
      "invalid_key",
      "pending_timeout",
      "download_failed",
      "killed",
      "cancelled",
      "file_erased",
      "import_failed",
      "split_failed",
      "internal_unreconciled"
    )

    unresolved_requests <- retrieval_audit |>
      dplyr::filter(.data$retrieval_status %in% unresolved_statuses)

    summary_tbl <- bf_repair_utf8_df(summary_tbl)
    retrieval_audit <- bf_repair_utf8_df(retrieval_audit)
    unresolved_requests <- bf_repair_utf8_df(unresolved_requests)

    readr::write_csv(summary_tbl, file.path(output_dir, "gbif_summary.csv"))
    readr::write_csv(
      retrieval_audit,
      file.path(output_dir, "gbif_retrieval_audit.csv")
    )
    readr::write_csv(
      unresolved_requests,
      file.path(output_dir, "gbif_unresolved_requests.csv")
    )
  }

  bf_console_bullet(
    paste0(
      "Marine pipeline completed in ",
      bf_console_elapsed(pipeline_start_time),
      "."
    ),
    quiet = quiet
  )

  if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
    return(.bf_marine_pipe_bind_rows(all_results))
  }

  invisible(NULL)
}


#' Run the marine pipeline through the legacy EEZ wrapper
#'
#' Compatibility wrapper for older code that still calls
#' `process_gbif_eez_pipeline()`. All arguments are forwarded unchanged to
#' [process_gbif_marine_pipeline()]. The wrapper is retained so existing scripts
#' can continue to run while the package uses the more general marine-overlay
#' workflow internally.
#'
#' @param ... Arguments passed directly to [process_gbif_marine_pipeline()].
#'
#' @return Returns the same object as [process_gbif_marine_pipeline()]. When
#'   `return_all_results = TRUE` and `store_in_memory = TRUE`, this is a tibble
#'   containing combined processed marine occurrence records with marine-region
#'   assignment, coordinate, workflow and optional evidence-audit columns.
#'   Otherwise, the wrapper writes outputs through
#'   [process_gbif_marine_pipeline()] and returns `invisible(NULL)`. This
#'   function is a backwards-compatible wrapper and is called primarily for the
#'   side effects of running the marine GBIF pipeline through the legacy EEZ
#'   interface.
#'
#' @family marine processing pipelines
#' @md
#' @export
process_gbif_eez_pipeline <- function(...) {
  process_gbif_marine_pipeline(...)
}
