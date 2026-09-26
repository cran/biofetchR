################################################################################
# utils_pipeline_origin_gates.R
# ------------------------------------------------------------------------------
# biofetchR: internal pre-download origin-evidence gates
# ------------------------------------------------------------------------------
#
# PURPOSE
#   This file contains the shared origin-evidence gate used by the retained
#   biofetchR GBIF processing pipelines. It applies optional GRIIS and
#   native-range checks after taxonomy preparation and before any GBIF download
#   request is submitted.
#
# CANONICAL VERSION
#   This is the canonical version of the origin-gate helper. It includes the
#   newer native-web evidence route that can compile species-level native-origin
#   evidence from SInAS, GBIF Species API distributions and WoRMS before applying
#   the native-range gate.
#
# WHY THIS FILE EXISTS
#   GBIF downloads are slow, quota-limited and externally dependent. Rows that
#   fail an origin-evidence rule should therefore be filtered before they reach
#   the download backend. Centralising the GRIIS/native-range filtering logic here
#   keeps the terrestrial/freshwater and marine pipelines consistent and reduces
#   duplicated audit code.
#
# PIPELINE POSITION
#   The intended order is:
#
#     1. Input validation and taxonomic cleaning.
#     2. Optional GRIIS evidence attachment/filtering.
#     3. Optional native-web evidence compilation from SInAS/GBIF/WoRMS.
#     4. Optional native-range evidence attachment/filtering.
#     5. Optional reconciliation of GRIIS and native-origin evidence.
#     6. GBIF download submission for retained rows only.
#
# DESIGN NOTES
#   - This file is deliberately internal; it does not export public functions.
#   - All helper functions are documented with `@noRd` so roxygen can process
#     them without creating public help pages.
#   - Terrestrial/freshwater workflows normally use species + recipient-country
#     evidence (`require_country = TRUE`).
#   - Marine/global workflows can use species-only evidence by leaving recipient
#     columns unset and using `require_country = FALSE`.
#   - Native-web evidence is evaluated after any GRIIS filtering so web/API
#     requests are made only for rows still eligible for GBIF submission.
#   - Audit outputs are written when requested so users can inspect which rows
#     were retained or rejected before GBIF submission.
#
# MAIN INTERNAL ENTRY POINT
#   .bf_pipeline_apply_origin_gates()
#
# DATA-SOURCE AND ATTRIBUTION NOTE
#   This helper does not itself define the external data licences, but it can
#   trigger downstream helpers that use GRIIS, SInAS, GBIF Species API and WoRMS
#   evidence. Analyses using those outputs should cite the specific source
#   datasets/APIs and retained cache/audit files used in the workflow.
#
################################################################################


#' Write an origin-gate audit table
#'
#' Writes a data frame to CSV and creates the parent directory if needed. The
#' helper prefers `readr::write_csv()` when available and falls back to
#' [utils::write.csv()] otherwise. This keeps audit writing robust in minimal
#' installations while retaining tidy CSV output when `readr` is installed.
#'
#' @param x Data frame-like object to write.
#' @param path Output CSV path.
#'
#' @return Invisibly returns `path`.
#'
#' @md
#' @keywords internal
#' @family pipeline origin-gate helpers
#' @noRd
.bf_pipeline_origin_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(x, path)
  } else {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }

  invisible(path)
}

#' Test whether a supplied column name is usable
#'
#' Checks that a candidate column name is non-null, scalar, non-missing,
#' non-empty and present in a data frame. This is used before passing optional
#' ISO/country columns to GRIIS or native-range helpers.
#'
#' @param col Candidate column name.
#' @param dat Data frame in which the column should exist.
#'
#' @return Logical scalar.
#'
#' @md
#' @keywords internal
#' @family pipeline origin-gate helpers
#' @noRd
.bf_pipeline_origin_valid_col <- function(col, dat) {
  !is.null(col) &&
    length(col) == 1L &&
    !is.na(col) &&
    nzchar(col) &&
    col %in% names(dat)
}

#' Remove temporary origin-gate bookkeeping columns
#'
#' Drops columns whose names start with `.bf_origin_gate_`. These columns are
#' needed while filters and audit tables are being built but should not be
#' returned to downstream GBIF download code.
#'
#' @param dat Data frame produced by the origin-gate workflow.
#'
#' @return Data frame with temporary origin-gate columns removed.
#'
#' @md
#' @keywords internal
#' @family pipeline origin-gate helpers
#' @noRd
.bf_pipeline_origin_drop_internal_cols <- function(dat) {
  internal <- grepl("^\\.bf_origin_gate_", names(dat))
  dat[, !internal, drop = FALSE]
}

#' Apply optional GRIIS and native-range gates before GBIF download
#'
#' Applies the shared origin-evidence stage used by biofetchR pipelines. The
#' function can attach GRIIS evidence, filter by GRIIS status, attach
#' native-range evidence, filter by native/non-native status, reconcile GRIIS and
#' native-origin evidence where possible, and write audit tables describing all
#' retained and rejected rows.
#'
#' @details
#' This helper is intentionally placed after taxonomic cleaning and before GBIF
#' download submission. Rows rejected here are not sent to GBIF. This prevents
#' unnecessary external download requests for rows that are not relevant to the
#' selected invasion/origin-evidence criteria.
#'
#' Matching behaviour depends on `require_country`:
#'
#' - `require_country = TRUE` is appropriate for terrestrial/freshwater
#'   species-country workflows. GRIIS and native-range evidence are matched using
#'   species plus recipient country where the called helper supports that mode.
#' - `require_country = FALSE` is appropriate for species-only marine/global
#'   screening, where the recipient unit is not an ISO country.
#'
#' The function writes audit files when `export_origin_audit = TRUE`. These files
#' are useful for debugging and for documenting which records entered GBIF
#' download submission.
#'
#' @section Native-web evidence:
#' When `use_native_web = TRUE`, native-origin evidence is compiled after the
#' GRIIS gate and before the native-range gate. This ordering is intentional: it
#' avoids web/API requests for rows already rejected by GRIIS and then passes the
#' collapsed species-level evidence table into `bf_attach_native_status()`. The
#' native-web route requires `use_native_filter = TRUE` and `native_ranges = NULL`.
#'
#' @section Audit outputs:
#' Depending on the selected options, the helper can write GRIIS audit/rejection
#' files, native-web long/species/unmapped/summary files, native-range
#' audit/rejection/summary files, the combined gate summary and the final
#' accepted-for-GBIF table. These files are intended for transparent debugging of
#' pre-download filtering decisions.
#'
#' @section Data-source and attribution:
#' This helper can call downstream GRIIS and native-web evidence utilities. Those
#' utilities may use GRIIS, SInAS, GBIF Species API and WoRMS evidence. biofetchR
#' records the audit/cache outputs, but users remain responsible for citing the
#' specific data products, APIs and versions used in analysis.
#'
#' @param df Input data frame after taxonomic preparation.
#' @param output_dir Directory for origin-gate audit outputs and default
#'   pipeline-local caches. Must be supplied explicitly.
#' @param species_col Name of the species column in `df`.
#' @param iso2c_col Optional ISO2 recipient-country column.
#' @param iso3c_col Optional ISO3 recipient-country column.
#' @param country_col Optional recipient-country name column.
#' @param require_country Logical. If `TRUE`, origin helpers are called in
#'   species + country mode where possible. If `FALSE`, species-only matching is
#'   used.
#' @param use_griis_filter Logical. If `TRUE`, attach GRIIS status and optionally
#'   filter rows using `griis_filter_mode`.
#' @param griis_filter_mode Character. One of `"audit_only"`,
#'   `"listed_invasive"` or `"listed_any"`. `"audit_only"` attaches GRIIS
#'   evidence without filtering; `"listed_invasive"` keeps only rows flagged
#'   invasive by GRIIS; `"listed_any"` keeps any row listed in GRIIS.
#' @param griis Optional preloaded or standardised GRIIS table. If `NULL`, the
#'   downstream GRIIS helper may read/download GRIIS using `griis_cache_dir`.
#' @param griis_cache_dir Cache directory passed to `bf_attach_griis_status()`
#'   when GRIIS must be read internally. If `NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used when needed.
#' @param griis_force_refresh Logical. Force GRIIS refresh where supported.
#' @param use_native_web Logical. If `TRUE`, fetch species-level native-origin
#'   evidence from package web/API helpers after the GRIIS gate and before the
#'   native-range gate. This is disabled by default because it performs external
#'   requests unless cached responses already exist.
#' @param native_web_sources Character vector of native web sources passed to
#'   `bf_fetch_native_ranges_web()`, currently usually `c("sinas", "gbif", "worms")`.
#' @param native_web_cache_dir Cache directory used by native web-source helpers.
#'   If `NULL`, a subdirectory under `output_dir/_biofetchR_cache/` is used when
#'   `use_native_web = TRUE`.
#' @param native_web_force_refresh Logical. If `TRUE`, ignore existing native-web
#'   cache files where supported.
#' @param native_web_sleep_sec Numeric delay in seconds between native-web
#'   requests.
#' @param native_web_sinas_main_path Optional local path to `SInAS_3.2.csv`
#'   when `"sinas"` is included in `native_web_sources`.
#' @param native_web_sinas_alllocations_path Optional local path to
#'   `AllLocations.xlsx`, `.csv` or `.tsv`.
#' @param native_web_sinas_fulltaxa_path Optional local path to
#'   `SInAS_3.2_FullTaxaList.csv`.
#' @param export_native_web_audit Logical. If `TRUE`, write native-web long,
#'   species, unmapped and summary outputs using `bf_write_native_web_outputs()`.
#' @param use_native_filter Logical. If `TRUE`, attach native-range status and
#'   optionally filter rows using `native_filter_mode`.
#' @param native_filter_mode Character. One of `"audit_only"`,
#'   `"non_native_only"`, `"non_native_or_unknown"` or `"native_only"`.
#' @param native_ranges Native-range table or lookup passed to
#'   `bf_attach_native_status()`. Required when `use_native_filter = TRUE`.
#' @param native_species_col Species column in `native_ranges`.
#' @param native_strategy Native-range matching strategy passed to
#'   `bf_attach_native_status()`.
#' @param native_keep_unknown Logical. In `"non_native_or_unknown"` mode, retain
#'   rows with uncertain or missing origin evidence when `TRUE`.
#' @param native_max_origins Maximum number of native-origin entries passed to
#'   native-range helpers that support a `max_origins` argument.
#' @param reconcile_origin_evidence Logical. If `TRUE`, run
#'   `bf_reconcile_griis_native_status()` when available and relevant evidence
#'   columns are present.
#' @param export_origin_audit Logical. If `TRUE`, write gate-specific audit,
#'   rejected, summary and accepted-for-GBIF CSV files.
#' @param audit_prefix Filename prefix for audit outputs.
#' @param quiet Logical. If `TRUE`, suppress progress messages.
#'
#' @return A tibble containing rows retained for GBIF download submission, with
#'   temporary origin-gate bookkeeping columns removed.
#'
#' @md
#' @keywords internal
#' @family pipeline origin-gate helpers
#' @noRd
.bf_pipeline_apply_origin_gates <- function(
    df,
    output_dir,
    species_col = "species",

    # Recipient-country columns. For terrestrial/freshwater this should usually
    # be iso2c. For marine species-only workflows, leave these NULL and set
    # require_country = FALSE.
    iso2c_col = NULL,
    iso3c_col = NULL,
    country_col = NULL,
    require_country = TRUE,

    # GRIIS gate
    use_griis_filter = FALSE,
    griis_filter_mode = c("audit_only", "listed_invasive", "listed_any"),
    griis = NULL,
    griis_cache_dir = NULL,
    griis_force_refresh = FALSE,

    # Optional web/API native-origin evidence provider. This is evaluated after
    # any GRIIS filtering and before the native-range gate, so only retained
    # rows are queried.
    use_native_web = FALSE,
    native_web_sources = c("sinas", "gbif", "worms"),
    native_web_cache_dir = NULL,
    native_web_force_refresh = FALSE,
    native_web_sleep_sec = 0.25,
    native_web_sinas_main_path = NULL,
    native_web_sinas_alllocations_path = NULL,
    native_web_sinas_fulltaxa_path = NULL,
    export_native_web_audit = TRUE,

    # Native-range gate
    use_native_filter = FALSE,
    native_filter_mode = c(
      "audit_only",
      "non_native_only",
      "non_native_or_unknown",
      "native_only"
    ),
    native_ranges = NULL,
    native_species_col = "species",
    native_strategy = "tiered",
    native_keep_unknown = TRUE,
    native_max_origins = 40,

    # Reconcile evidence
    reconcile_origin_evidence = TRUE,

    # Output/audit behaviour
    export_origin_audit = TRUE,
    audit_prefix = "origin",

    quiet = FALSE
) {
  griis_filter_mode <- match.arg(griis_filter_mode)
  native_filter_mode <- match.arg(native_filter_mode)

  if (isTRUE(use_native_web) && !isTRUE(use_native_filter)) {
    stop(
      "`use_native_web = TRUE` requires `use_native_filter = TRUE`, because web-derived native-origin evidence is only used by the native-range gate.",
      call. = FALSE
    )
  }

  if (isTRUE(use_native_web) && !is.null(native_ranges)) {
    stop(
      "`use_native_web = TRUE` currently requires `native_ranges = NULL`. Use either supplied native-range evidence or web-derived evidence, not both.",
      call. = FALSE
    )
  }

  if (!is.data.frame(df)) {
    stop("`df` must be a data frame before applying origin gates.", call. = FALSE)
  }

  if (!species_col %in% names(df)) {
    stop(
      "`species_col` was not found in `df` before applying origin gates: ",
      species_col,
      call. = FALSE
    )
  }

  if (is.null(output_dir) || length(output_dir) == 0L ||
      !nzchar(trimws(as.character(output_dir[[1L]])))) {
    stop("`output_dir` must be supplied explicitly.", call. = FALSE)
  }

  output_dir <- normalizePath(
    as.character(output_dir[[1L]]),
    winslash = "/",
    mustWork = FALSE
  )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  .bf_pipeline_origin_cache_dir <- function(cache_dir, subdir) {
    if (is.null(cache_dir) || length(cache_dir) == 0L ||
        !nzchar(trimws(as.character(cache_dir[[1L]])))) {
      return(file.path(output_dir, "_biofetchR_cache", subdir))
    }

    normalizePath(
      as.character(cache_dir[[1L]]),
      winslash = "/",
      mustWork = FALSE
    )
  }

  griis_cache_dir <- .bf_pipeline_origin_cache_dir(
    griis_cache_dir,
    "griis"
  )

  native_web_cache_dir <- .bf_pipeline_origin_cache_dir(
    native_web_cache_dir,
    "native_web"
  )

  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required for pipeline origin gates.", call. = FALSE)
  }

  out <- tibble::as_tibble(df)
  out$.bf_origin_gate_rowid <- seq_len(nrow(out))

  summary_rows <- list()

  add_summary <- function(stage,
                          mode,
                          n_before,
                          n_after,
                          n_rejected,
                          detail = NA_character_) {
    summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
      stage = as.character(stage),
      mode = as.character(mode),
      n_before = as.integer(n_before),
      n_after = as.integer(n_after),
      n_rejected = as.integer(n_rejected),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  }

  add_summary(
    stage = "input_after_taxonomy_gate",
    mode = "none",
    n_before = nrow(out),
    n_after = nrow(out),
    n_rejected = 0L,
    detail = "Rows entering GRIIS/native origin gates"
  )

  # ---------------------------------------------------------------------------
  # Optional GRIIS status attach/filter
  # ---------------------------------------------------------------------------

  if (isTRUE(use_griis_filter)) {
    if (!exists("bf_attach_griis_status", mode = "function", inherits = TRUE)) {
      stop(
        "`use_griis_filter = TRUE`, but `bf_attach_griis_status()` is not visible.",
        call. = FALSE
      )
    }

    n_before_griis <- nrow(out)

    iso2_arg <- if (.bf_pipeline_origin_valid_col(iso2c_col, out)) iso2c_col else NULL
    iso3_arg <- if (.bf_pipeline_origin_valid_col(iso3c_col, out)) iso3c_col else NULL
    country_arg <- if (.bf_pipeline_origin_valid_col(country_col, out)) country_col else NULL

    out <- bf_attach_griis_status(
      df = out,
      species_col = species_col,
      iso2c_col = iso2_arg,
      iso3c_col = iso3_arg,
      country_col = country_arg,
      griis = griis,
      cache_dir = griis_cache_dir,
      force_refresh = griis_force_refresh,
      require_country = require_country,
      quiet = quiet
    )

    out$.bf_origin_gate_keep_griis <- switch(
      griis_filter_mode,
      audit_only = rep(TRUE, nrow(out)),
      listed_invasive = out$griis_invasive %in% TRUE,
      listed_any = out$griis_listed %in% TRUE
    )

    out$.bf_origin_gate_reject_griis <- ifelse(
      out$.bf_origin_gate_keep_griis,
      NA_character_,
      switch(
        griis_filter_mode,
        audit_only = NA_character_,
        listed_invasive = "not_griis_listed_invasive",
        listed_any = "not_griis_listed"
      )
    )

    if (isTRUE(export_origin_audit)) {
      .bf_pipeline_origin_write_csv(
        out,
        file.path(output_dir, paste0(audit_prefix, "_griis_audit.csv"))
      )

      griis_rejected <- out[!out$.bf_origin_gate_keep_griis, , drop = FALSE]

      .bf_pipeline_origin_write_csv(
        griis_rejected,
        file.path(output_dir, paste0(audit_prefix, "_griis_rejected.csv"))
      )
    }

    out <- out[out$.bf_origin_gate_keep_griis, , drop = FALSE]

    add_summary(
      stage = "griis_gate",
      mode = griis_filter_mode,
      n_before = n_before_griis,
      n_after = nrow(out),
      n_rejected = n_before_griis - nrow(out),
      detail = paste0(
        "GRIIS mode: ",
        griis_filter_mode,
        "; require_country=",
        require_country
      )
    )

    .bf_msg(
      "biofetchR GRIIS gate: retained ",
      nrow(out),
      " / ",
      n_before_griis,
      " row(s) using mode `",
      griis_filter_mode,
      "`.",
      quiet = quiet
    )
  } else {
    add_summary(
      stage = "griis_gate",
      mode = "disabled",
      n_before = nrow(out),
      n_after = nrow(out),
      n_rejected = 0L,
      detail = "GRIIS gate disabled"
    )
  }

  # ---------------------------------------------------------------------------
  # Optional web/API native-origin evidence fetch
  # ---------------------------------------------------------------------------
  # This step is intentionally placed after the GRIIS gate and before the
  # native-range gate. That order ensures web/API requests are made only for rows
  # still eligible for GBIF submission and that the returned species-level table
  # can be passed directly to bf_attach_native_status().

  if (isTRUE(use_native_web)) {
    if (!exists("bf_fetch_native_ranges_web", mode = "function", inherits = TRUE)) {
      stop(
        "`use_native_web = TRUE`, but `bf_fetch_native_ranges_web()` is not visible. ",
        "Make sure `R/utils_native_range_web_sources.R` is present and loaded by the package.",
        call. = FALSE
      )
    }

    species_for_native_web <- unique(as.character(out[[species_col]]))
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
          "`export_native_web_audit = TRUE`, but `bf_write_native_web_outputs()` is not visible.",
          call. = FALSE
        )
      }

      bf_write_native_web_outputs(
        native_web,
        output_dir = output_dir,
        prefix = paste0(audit_prefix, "_native_web")
      )
    }

    native_ranges <- native_web$species
    native_species_col <- "species"

    add_summary(
      stage = "native_web_evidence",
      mode = paste(native_web_sources, collapse = ";"),
      n_before = length(species_for_native_web),
      n_after = if (is.data.frame(native_ranges)) sum(native_ranges$native_has_origin %in% TRUE, na.rm = TRUE) else NA_integer_,
      n_rejected = 0L,
      detail = paste0(
        "Native web evidence fetched for ",
        length(species_for_native_web),
        " species; cache directory = ",
        native_web_cache_dir
      )
    )

    .bf_msg(
      "biofetchR native web evidence: fetched species-level origin evidence for ",
      length(species_for_native_web),
      " species.",
      quiet = quiet
    )
  } else {
    add_summary(
      stage = "native_web_evidence",
      mode = "disabled",
      n_before = nrow(out),
      n_after = nrow(out),
      n_rejected = 0L,
      detail = "Native web evidence disabled"
    )
  }

  # ---------------------------------------------------------------------------
  # Optional native-range status attach/filter
  # ---------------------------------------------------------------------------

  if (isTRUE(use_native_filter)) {
    if (!exists("bf_attach_native_status", mode = "function", inherits = TRUE)) {
      stop(
        "`use_native_filter = TRUE`, but `bf_attach_native_status()` is not visible.",
        call. = FALSE
      )
    }

    if (is.null(native_ranges)) {
      stop(
        "`use_native_filter = TRUE`, but `native_ranges` was not supplied.",
        call. = FALSE
      )
    }

    n_before_native <- nrow(out)

    iso2_arg <- if (.bf_pipeline_origin_valid_col(iso2c_col, out)) iso2c_col else NULL
    iso3_arg <- if (.bf_pipeline_origin_valid_col(iso3c_col, out)) iso3c_col else NULL
    country_arg <- if (.bf_pipeline_origin_valid_col(country_col, out)) country_col else NULL

    native_args <- list(
      df = out,
      species_col = species_col,
      iso2c_col = iso2_arg,
      iso3c_col = iso3_arg,
      country_col = country_arg,
      native_ranges = native_ranges,
      native_species_col = native_species_col,
      require_country = require_country,
      strategy = native_strategy,
      quiet = quiet
    )

    native_formals <- names(formals(bf_attach_native_status))

    if ("native_filter_mode" %in% native_formals) {
      native_args$native_filter_mode <- "audit_only"
    }

    if ("max_origins" %in% native_formals) {
      native_args$max_origins <- native_max_origins
    }

    out <- do.call(bf_attach_native_status, native_args)

    if (
      isTRUE(reconcile_origin_evidence) &&
      exists("bf_reconcile_griis_native_status", mode = "function", inherits = TRUE) &&
      any(c("griis_listed", "griis_invasive") %in% names(out))
    ) {
      out <- bf_reconcile_griis_native_status(out)
    }

    native_unknown <- (
      out$native_status %in% c(
        "origin_unknown",
        "recipient_country_missing",
        "cosmopolitan_or_uncertain",
        "origin_available_species_only"
      )
    ) | !(out$native_has_origin %in% TRUE)

    out$.bf_origin_gate_keep_native <- switch(
      native_filter_mode,
      audit_only = rep(TRUE, nrow(out)),
      non_native_only = out$native_is_non_native_recipient %in% TRUE,
      non_native_or_unknown = {
        if (isTRUE(native_keep_unknown)) {
          (out$native_is_non_native_recipient %in% TRUE) | native_unknown
        } else {
          out$native_is_non_native_recipient %in% TRUE
        }
      },
      native_only = out$native_is_native_recipient %in% TRUE
    )

    out$.bf_origin_gate_reject_native <- ifelse(
      out$.bf_origin_gate_keep_native,
      NA_character_,
      switch(
        native_filter_mode,
        audit_only = NA_character_,
        non_native_only = "not_confirmed_non_native",
        non_native_or_unknown = "not_non_native_or_unknown",
        native_only = "not_confirmed_native"
      )
    )

    if (isTRUE(export_origin_audit)) {
      .bf_pipeline_origin_write_csv(
        out,
        file.path(output_dir, paste0(audit_prefix, "_native_audit.csv"))
      )

      native_rejected <- out[!out$.bf_origin_gate_keep_native, , drop = FALSE]

      .bf_pipeline_origin_write_csv(
        native_rejected,
        file.path(output_dir, paste0(audit_prefix, "_native_rejected.csv"))
      )

      if (exists("bf_native_status_summary", mode = "function", inherits = TRUE)) {
        native_summary <- bf_native_status_summary(out)

        .bf_pipeline_origin_write_csv(
          native_summary,
          file.path(output_dir, paste0(audit_prefix, "_native_summary.csv"))
        )
      }
    }

    out <- out[out$.bf_origin_gate_keep_native, , drop = FALSE]

    add_summary(
      stage = "native_range_gate",
      mode = native_filter_mode,
      n_before = n_before_native,
      n_after = nrow(out),
      n_rejected = n_before_native - nrow(out),
      detail = paste0(
        "Native mode: ",
        native_filter_mode,
        "; require_country=",
        require_country,
        "; keep_unknown=",
        native_keep_unknown
      )
    )

    .bf_msg(
      "biofetchR native-range gate: retained ",
      nrow(out),
      " / ",
      n_before_native,
      " row(s) using mode `",
      native_filter_mode,
      "`.",
      quiet = quiet
    )
  } else {
    add_summary(
      stage = "native_range_gate",
      mode = "disabled",
      n_before = nrow(out),
      n_after = nrow(out),
      n_rejected = 0L,
      detail = "Native-range gate disabled"
    )
  }

  # ---------------------------------------------------------------------------
  # Final gate summary
  # ---------------------------------------------------------------------------

  add_summary(
    stage = "final_before_gbif_download",
    mode = "final",
    n_before = nrow(df),
    n_after = nrow(out),
    n_rejected = nrow(df) - nrow(out),
    detail = "Rows passed to GBIF download stage"
  )

  if (isTRUE(export_origin_audit)) {
    summary_tbl <- do.call(rbind, summary_rows)

    .bf_pipeline_origin_write_csv(
      summary_tbl,
      file.path(output_dir, paste0(audit_prefix, "_gate_summary.csv"))
    )

    .bf_pipeline_origin_write_csv(
      out,
      file.path(output_dir, paste0(audit_prefix, "_accepted_for_gbif.csv"))
    )
  }

  out <- .bf_pipeline_origin_drop_internal_cols(out)

  tibble::as_tibble(out)
}
