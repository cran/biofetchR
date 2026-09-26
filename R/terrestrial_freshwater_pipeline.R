################################################################################
# terrestrial_freshwater_pipeline.R
# ------------------------------------------------------------------------------
# biofetchR terrestrial and freshwater GBIF processing pipeline
# ------------------------------------------------------------------------------
#
# Package role
#   This file defines the main terrestrial/freshwater occurrence-processing
#   workflow used by biofetchR. Inputs remain species x recipient-country tables,
#   but accepted countries for the same species can be combined into one
#   asynchronous GBIF download, imported once, and split back to country-level
#   occurrence objects before downstream overlay assignment, cleaning, thinning
#   and export.
#
# Data access and attribution
#   This pipeline can submit GBIF downloads and can use multiple third-party
#   spatial or raster resources through helper functions in other biofetchR
#   files. biofetchR does not redistribute those provider datasets from this
#   script. Users should retain GBIF download keys/DOIs and follow the licence,
#   citation and attribution requirements of every data source requested in a
#   given run, including GBIF, GADM, HydroSHEDS/HydroLAKES/HydroRIVERS,
#   TEOW/FEOW, RESOLVE, WDPA, Ramsar, GMBA, Natural Earth, WorldCover,
#   SoilGrids, Human Footprint, Global Dam Watch and any user-supplied layers.
#
# Core responsibilities
#   1. Validate species x ISO2 input tables.
#   2. Apply optional taxonomy preparation before GBIF requests are made.
#   3. Apply optional GRIIS and native-range origin-evidence gates before GBIF
#      requests are made.
#   4. Pre-check retained species x country combinations for georeferenced GBIF
#      records.
#   5. For each species, combine all retained positive countries into one
#      multi-country asynchronous GBIF download request.
#   6. Import each successful GBIF archive once and reconstruct country-specific
#      occurrence objects using GBIF `countryCode`.
#   7. Reconcile every retained species x country retrieval request to an explicit
#      success, zero, failure, or unresolved status.
#   8. Optionally enrich occurrence points with raster-derived context fields.
#   9. Assign points to requested terrestrial/freshwater vector overlays.
#  10. Apply coordinate cleaning and, optionally, spatial thinning.
#  11. Export per-species/per-region CSVs plus processing and GBIF retrieval
#      audits.
#
# Design principles
#   - Conservative ordering: taxonomy and origin filters run before GBIF
#     submission, so rejected rows are never submitted to GBIF.
#   - Retrieval efficiency: species-country combinations remain the logical audit
#     units, while several countries for one species can share one GBIF download
#     key and one imported archive.
#   - Retrieval completeness: every retained species x country request receives an
#     explicit retrieval outcome; missing requests are marked
#     `internal_unreconciled` rather than disappearing silently.
#   - Timeout integrity: a polling timeout remains an unresolved download with its
#     GBIF key retained and is never converted to a zero-occurrence result.
#   - Fail-soft by default: overlay/context failures are recorded in the summary
#     rather than aborting an entire run, unless strict_* arguments are TRUE.
#   - Reproducible auditing: taxonomy, origin and GBIF retrieval decisions can be
#     written to CSV.
#   - Testable internals: the `deps` argument allows dependency injection for unit
#     tests without changing the public API.
#   - Data-processing scope: this file produces processed occurrence tables,
#     grouped CSV outputs and audit files only.
#
# Supported terrestrial/freshwater vector overlays
#   - gadm: administrative units from GADM.
#   - ne_admin1: Natural Earth first-level administrative units.
#   - ne_urban: Natural Earth urban areas.
#   - resolve2017: RESOLVE 2017 terrestrial ecoregions.
#   - teow: Terrestrial Ecoregions of the World.
#   - feow: Freshwater Ecoregions of the World.
#   - basins: HydroBASINS catchment polygons.
#   - lakes: HydroLAKES polygons.
#   - rivers: HydroRIVERS line features, assigned by nearest-feature snapping.
#   - gmba: Global Mountain Biodiversity Assessment mountain inventory.
#   - wdpa: World Database on Protected Areas polygons.
#   - ramsar: Ramsar wetlands.
#   - gdw_barriers: Global Dam Watch barrier features.
#   - biosphere_reserve: UNESCO biosphere-reserve features.
#   - gloric: GLORIC river-class features.
#   - hydrowaste: HydroWASTE wastewater-treatment features.
#   - gdw_reservoirs: Global Dam Watch reservoir polygons.
#   - global_mining: global mining-area polygons.
#
# Optional raster/context enrichment
#   Raster contexts add attributes to occurrence points and do not create new
#   grouped outputs by themselves. Currently supported values are:
#   - worldcover
#   - soilgrids
#   - human_footprint
#
################################################################################

#' Process terrestrial and freshwater GBIF occurrences
#'
#' Download, import, clean, assign and export GBIF occurrence records for
#' terrestrial and freshwater workflows. Input is supplied as a species x ISO2
#' country table. Species-country combinations remain the logical retrieval,
#' auditing and output units, but eligible countries for one species are combined
#' into a single asynchronous GBIF occurrence download whenever the built-in
#' multi-country backend is used.
#'
#' @description
#' `process_gbif_terrestrial_freshwater_pipeline()` is the main high-level
#' terrestrial/freshwater occurrence workflow in biofetchR. It applies optional
#' taxonomy and origin-evidence gates before GBIF submission, pre-checks retained
#' species-country requests for georeferenced records, submits one combined
#' multi-country download per species, imports that archive once, and reconstructs
#' country-specific occurrence objects before the existing context, overlay,
#' cleaning, thinning and export stages.
#'
#' Retrieval provenance is retained separately from downstream spatial-output
#' summaries so every original request can be classified explicitly as successful,
#' empty, failed, or unresolved.
#'
#' @details
#' The pipeline runs in the following order:
#'
#' 1. **Input validation and standardisation.** The input table must contain a
#'    species column and an ISO2 country column. Species names are coerced to
#'    character and ISO2 codes are upper-cased.
#' 2. **Optional taxonomy gate.** When `prepare_taxonomy = TRUE`, names are
#'    cleaned, manual fixes can be applied, invalid/non-species entries can be
#'    rejected, and accepted names are used downstream.
#' 3. **Optional origin-evidence gates.** GRIIS and native-range evidence can be
#'    attached and used to retain or reject species-country rows before GBIF is
#'    queried.
#' 4. **Country-level GBIF availability pre-check.** Each retained species-country
#'    combination is checked for georeferenced GBIF records. A confirmed zero is
#'    recorded explicitly as `precheck_zero` and is not included in the download.
#' 5. **Multi-country GBIF submission.** For each species, all countries retained
#'    after the pre-check are submitted together to the built-in backend as one
#'    asynchronous GBIF download using a country-IN predicate.
#' 6. **GBIF polling, import and country reconstruction.** The shared download is
#'    polled and imported once. Its records are split locally by GBIF
#'    `countryCode`, recreating one result object for every requested country.
#' 7. **Retrieval reconciliation.** Every retained species-country request is
#'    reconciled against the submission/import audit. Successful empty subsets,
#'    failures, timeouts and internally missing outcomes remain explicit.
#' 8. **Optional raster/context enrichment.** Requested raster/context variables
#'    are extracted to occurrence points as additional attributes.
#' 9. **Vector overlay assignment.** Points are assigned to requested
#'    terrestrial/freshwater spatial units using polygon containment or
#'    nearest-feature snapping, depending on overlay type.
#' 10. **Coordinate cleaning and thinning.** Coordinate-quality cleaning is
#'     applied when `apply_cleaning = TRUE`; distance-based thinning is applied
#'     when `apply_thinning = TRUE` and `dist_km > 0`.
#' 11. **Export and auditing.** Grouped occurrence CSVs, processing summaries and
#'     GBIF retrieval-audit files are written to `output_dir`. If
#'     `store_in_memory = TRUE`, grouped records can also be returned as a
#'     combined tibble.
#'
#' @section GBIF multi-country retrieval architecture:
#' The input and audit unit remains a species x ISO2-country combination. The
#' asynchronous GBIF download unit is different: when the built-in
#' `download_gbif_batch_gadm()` backend is used, all retained countries for one
#' species are combined in a single country-IN request.
#'
#' Consequently, several countries for the same species normally carry the SAME
#' GBIF download key. This is intentional. The completed archive is imported only
#' once and [wait_and_import_gbif()] splits it back to country-specific objects
#' using GBIF `countryCode`.
#'
#' This design reduces unnecessary asynchronous download jobs while preserving
#' country-level outputs and provenance. A successful country split containing no
#' records is retained as `success_zero`; it is not silently omitted.
#'
#' @section Input table:
#' `df` must contain at least:
#'
#' - `species`: submitted taxon name or common-name/manual-fix input, depending
#'   on the taxonomy settings.
#' - `iso2c`: recipient ISO2 country code used for terrestrial/freshwater GBIF
#'   country predicates.
#'
#' Additional columns are retained where possible and can be carried into
#' summaries and audit files as metadata.
#'
#' @section Supported vector overlays:
#' Use `region_source` to request one or more overlays. Supported values are
#' `"gadm"`, `"ne_admin1"`, `"ne_urban"`, `"resolve2017"`, `"teow"`,
#' `"feow"`, `"basins"`, `"lakes"`, `"rivers"`, `"gmba"`, `"wdpa"`,
#' `"ramsar"`, `"gdw_barriers"`, `"biosphere_reserve"`, `"gloric"`,
#' `"hydrowaste"`, `"gdw_reservoirs"`, and `"global_mining"`.
#'
#' Polygon overlays assign points by spatial intersection. Line or point-like
#' resources such as rivers, barriers, HydroWASTE, biosphere reserves and GLORIC
#' use nearest-feature assignment with overlay-specific maximum snapping
#' distances.
#'
#' @section Origin-evidence gates:
#' GRIIS and native-range evidence are optional and are evaluated before GBIF
#' download submission. For terrestrial/freshwater workflows, country-aware
#' matching is the default: GRIIS and native-range evidence are interpreted for
#' each species in the recipient country when the relevant data are available.
#' Audit files can be written by setting `export_origin_audit = TRUE`.
#'
#' @section Cleaning and thinning:
#' Coordinate cleaning and thinning are delegated to `thin_spatial_points()`.
#' Cleaning removes invalid, impossible, missing and zero-zero coordinates and can
#' filter records with high coordinate uncertainty. Thinning is distance-based and
#' prioritises higher-quality records before removing nearby clustered records.
#'
#' @section Output files:
#' The pipeline can write:
#'
#' - grouped occurrence CSVs for each species-region combination;
#' - `gbif_summary.csv`, recording processing status, row counts, download keys and
#'   failure stages;
#' - `gbif_retrieval_audit.csv`, recording one explicit GBIF retrieval outcome
#'   for each taxonomy-approved species x country request;
#' - `gbif_unresolved_requests.csv`, containing only unresolved or failed GBIF
#'   retrieval requests for targeted review/recovery;
#' - `taxonomy_audit.csv`, `taxonomy_rejected.csv` and `taxonomy_summary.csv`
#'   when taxonomy auditing is enabled;
#' - `origin_evidence_audit.csv`, `origin_evidence_rejected.csv` and
#'   `origin_evidence_summary.csv` when origin-evidence auditing is enabled.
#'
#'
#' @section Retrieval audit and completeness:
#' `gbif_retrieval_audit.csv` is deliberately separate from
#' `gbif_summary.csv`. The processing summary can contain multiple rows per
#' original request after spatial overlay assignment, whereas the retrieval audit
#' contains one row per retained GBIF request unit.
#'
#' The retrieval audit contains:
#'
#' \itemize{
#'   \item `species`;
#'   \item `request_id` (ISO2 country);
#'   \item `request_type` (`"country"`);
#'   \item `gbif_key`;
#'   \item `gbif_status`;
#'   \item `retrieval_status`;
#'   \item `n_records`;
#'   \item `fail_reason`.
#' }
#'
#' Resolved retrieval outcomes include `success`, `success_zero`, and
#' `precheck_zero`. Unresolved or failed states written to
#' `gbif_unresolved_requests.csv` include `taxon_key_failed`, `submit_failed`,
#' `submit_no_key`, `submit_timeout`, `invalid_key`, `pending_timeout`,
#' `download_failed`, `killed`, `cancelled`, `file_erased`, `import_failed`,
#' `split_failed`, and `internal_unreconciled`.
#'
#' The final completeness gate compares the retained expected species-country
#' requests with the retrieval ledger. Any request lacking an explicit outcome is
#' added as `internal_unreconciled` rather than being silently omitted.
#'
#' @section Data access, licences and attribution:
#' This function can submit live GBIF occurrence downloads and can call helper
#' functions that load third-party vector and raster resources. biofetchR does
#' not redistribute those provider datasets from this pipeline, but users are
#' responsible for complying with the licence, citation and attribution terms of
#' every dataset requested in a run. Retain GBIF download keys and cite the
#' corresponding GBIF download DOI where GBIF-mediated records are used.
#'
#' Licence-sensitive or attribution-sensitive resources may include, depending
#' on the arguments used, GBIF, GADM, Natural Earth, HydroSHEDS/HydroBASINS,
#' HydroLAKES, HydroRIVERS, TEOW, FEOW, RESOLVE ecoregions, WDPA, Ramsar, GMBA,
#' Global Dam Watch, HydroWASTE, global mining layers, ESA WorldCover,
#' SoilGrids, Human Footprint and user-supplied spatial or raster layers. Always
#' cite the exact product, version, year and provider used in the workflow.
#'
#' @section Failure handling:
#' Overlay and raster-context failures are handled according to
#' `strict_overlay_loading` and `strict_raster_context`. In fail-soft mode these
#' failures are recorded in `gbif_summary.csv` and processing continues where
#' possible.
#'
#' GBIF retrieval failures use a separate request-level ledger. Submission
#' failures, invalid keys, download terminal states, polling timeouts, import
#' failures, split failures and completeness failures are written to
#' `gbif_retrieval_audit.csv`. Unresolved states are also copied to
#' `gbif_unresolved_requests.csv`.
#'
#' A polling timeout is recorded as `pending_timeout` with the GBIF key retained.
#' It is not treated as evidence of zero occurrence records. Likewise, a request
#' that cannot be reconciled to either an imported result or an explicit failure
#' is recorded as `internal_unreconciled`.
#'
#' Strict overlay/context modes can still stop the pipeline immediately for
#' configuration or resource failures.
#'
#' @section Testing and dependency injection:
#' The `deps` argument allows tests to inject mock download, import, export,
#' summary, overlay and thinning helpers. Normal users should usually leave
#' `deps = NULL`; package tests can use it to exercise the pipeline without
#' submitting live GBIF downloads or requiring large remote overlay resources.
#'
#' @param df Data frame containing at least `species` and `iso2c`.
#' @param output_dir Directory where occurrence CSVs, `gbif_summary.csv`,
#'   `gbif_retrieval_audit.csv`, `gbif_unresolved_requests.csv`, and optional
#'   taxonomy/origin audit files are written.
#' @param user GBIF username.
#' @param pwd GBIF password.
#' @param email GBIF account email address.
#' @param region_source Character vector of terrestrial/freshwater vector overlay
#'   sources. See **Supported vector overlays** for accepted values.
#' @param raster_context Optional character vector of raster/context enrichments.
#'   Accepted values are `"worldcover"`, `"soilgrids"`, and
#'   `"human_footprint"`. These enrich occurrence points with attributes and do
#'   not create grouped outputs by themselves.
#' @param overlay_mode Overlay combination mode. `"single"` processes each
#'   requested overlay independently. `"dual_separate"` keeps requested overlays
#'   separate. `"dual_intersection"` exports groups based on combined overlay
#'   labels when multiple overlays are requested.
#' @param gadm_unit Integer GADM administrative level used when
#'   `region_source` includes `"gadm"`. Supported values are `0`, `1`, and `2`.
#' @param cache_dir_gadm Directory for cached GADM objects. If `NULL`, a
#'   subdirectory under `output_dir/_biofetchR_cache/` is used.
#'
#' @param teow_cache_dir,teow_method,teow_url TEOW loading and source settings.
#' @param teow_force_refresh,teow_clip_to_countries TEOW refresh and clipping
#'   settings. If `teow_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param feow_cache_dir,feow_method,feow_url FEOW loading and source settings.
#' @param feow_force_refresh,feow_clip_to_countries FEOW refresh and clipping
#'   settings. If `feow_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param lakes_cache_dir,lakes_force_refresh HydroLAKES loading and caching
#'   settings.
#' @param lakes_only,lakes_min_area_km2 HydroLAKES filtering settings. If
#'   `lakes_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param rivers_cache_dir,rivers_cache,rivers_cache_format HydroRIVERS caching
#'   settings.
#' @param rivers_force_refresh,rivers_regions HydroRIVERS refresh and regional
#'   loading settings.
#' @param rivers_min_strahler,rivers_min_discharge_cms HydroRIVERS filtering
#'   settings.
#' @param rivers_max_snap_m Maximum snapping distance in metres for assigning
#'   points to nearest HydroRIVERS reaches. If `rivers_cache_dir = NULL`, a
#'   subdirectory under `output_dir/_biofetchR_cache/` is used.
#'
#' @param basins_cache_dir,basins_level HydroBASINS loading settings.
#' @param basins_with_lakes,basins_clip_to_countries HydroBASINS filtering and
#'   clipping settings. If `basins_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param gmba_cache_dir,gmba_layer,gmba_extent GMBA loading settings.
#' @param gmba_force_refresh,gmba_clip_to_countries GMBA refresh and clipping
#'   settings. If `gmba_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param ne_cache_dir,ne_force_refresh Natural Earth loading and caching
#'   settings.
#' @param ne_admin1_scale Natural Earth administrative-unit scale. If
#'   `ne_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param resolve_cache_dir,resolve_force_refresh RESOLVE 2017 loading and
#'   caching settings.
#' @param resolve_clip_to_countries RESOLVE 2017 clipping setting. If
#'   `resolve_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param wdpa_cache_dir,wdpa_force_refresh WDPA loading and caching settings.
#' @param wdpa_exclude_marine WDPA filtering setting controlling whether marine
#'   protected areas are excluded.
#' @param wdpa_require_opt_in,wdpa_opt_in WDPA opt-in controls. WDPA loading is
#'   opt-in by default because this source can be large and may require user
#'   awareness of provider terms. If `wdpa_cache_dir = NULL`, a subdirectory
#'   under `output_dir/_biofetchR_cache/` is used.
#'
#' @param ramsar_cache_dir,ramsar_force_refresh Ramsar loading and caching
#'   settings.
#' @param ramsar_require_opt_in,ramsar_opt_in Ramsar opt-in controls. Ramsar
#'   loading is opt-in by default. If `ramsar_cache_dir = NULL`, a subdirectory
#'   under `output_dir/_biofetchR_cache/` is used.
#'
#' @param gdw_cache_dir,gdw_force_refresh Global Dam Watch barrier loading and
#'   caching settings.
#' @param gdw_barriers_source_path,gdw_barriers_source_url Optional local path or
#'   URL for Global Dam Watch barrier data.
#' @param gdw_barrier_max_snap_m Maximum snapping distance in metres for Global
#'   Dam Watch barrier assignment. If `gdw_cache_dir = NULL`, a subdirectory
#'   under `output_dir/_biofetchR_cache/` is used.
#'
#' @param biosphere_cache_dir,biosphere_force_refresh UNESCO biosphere-reserve
#'   loading and caching settings.
#' @param biosphere_source_path,biosphere_source_url Optional local path or URL
#'   for biosphere-reserve data.
#' @param biosphere_max_snap_m Maximum snapping distance in metres for
#'   biosphere-reserve assignment. If `biosphere_cache_dir = NULL`, a
#'   subdirectory under `output_dir/_biofetchR_cache/` is used.
#'
#' @param gloric_cache_dir,gloric_force_refresh GLORIC loading and caching
#'   settings.
#' @param gloric_source_path,gloric_source_url Optional local path or URL for
#'   GLORIC data.
#' @param gloric_max_snap_m Maximum snapping distance in metres for GLORIC
#'   assignment. If `gloric_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param hydrowaste_cache_dir,hydrowaste_force_refresh HydroWASTE loading and
#'   caching settings.
#' @param hydrowaste_source_path,hydrowaste_source_url Optional local path or URL
#'   for HydroWASTE data.
#' @param hydrowaste_max_snap_m Maximum snapping distance in metres for
#'   HydroWASTE assignment. If `hydrowaste_cache_dir = NULL`, a subdirectory
#'   under `output_dir/_biofetchR_cache/` is used.
#'
#' @param gdw_reservoirs_cache_dir,gdw_reservoirs_force_refresh Global Dam Watch
#'   reservoir loading and caching settings.
#' @param gdw_reservoirs_source_path,gdw_reservoirs_source_url Optional local
#'   path or URL for Global Dam Watch reservoir data. If
#'   `gdw_reservoirs_cache_dir = NULL`, a subdirectory under
#'   `output_dir/_biofetchR_cache/` is used.
#'
#' @param global_mining_cache_dir,global_mining_force_refresh Global mining-area
#'   loading and caching settings.
#' @param global_mining_source_path,global_mining_source_url Optional local path
#'   or URL for global mining-area data. If `global_mining_cache_dir = NULL`, a
#'   subdirectory under `output_dir/_biofetchR_cache/` is used.
#' @param strict_overlay_loading Logical. If `TRUE`, a requested overlay-loading
#'   failure stops the pipeline. If `FALSE`, the failure is recorded in
#'   `gbif_summary.csv` and processing continues where possible.
#' @param strict_raster_context Logical. If `TRUE`, a requested raster/context
#'   enrichment failure stops the pipeline. If `FALSE`, the failure is recorded
#'   and processing continues with unenriched points.
#' @param raster_cache_dir Directory for raster/context caches. If `NULL`, a
#'   subdirectory under `output_dir/_biofetchR_cache/` is used.
#' @param worldcover_source Path, URL or supported raster object for WorldCover
#'   extraction.
#' @param worldcover_year Year label used for WorldCover output metadata.
#' @param worldcover_buffer_m Buffer radius in metres for WorldCover extraction.
#'   `0` uses point extraction.
#' @param soilgrids_sources Named vector/list of paths, URLs or supported raster
#'   objects keyed by SoilGrids variable name.
#' @param soilgrids_vars Character vector of SoilGrids variable names to extract;
#'   must match names in `soilgrids_sources`.
#' @param soilgrids_buffer_m Buffer radius in metres for SoilGrids extraction.
#'   `0` uses point extraction.
#' @param human_footprint_source Path, URL or supported raster object for Human
#'   Footprint extraction.
#' @param human_footprint_buffer_m Buffer radius in metres for Human Footprint
#'   extraction. `0` uses point extraction.
#' @param batch_size Number of species grouped into each outer processing batch.
#'   Species are still processed sequentially within those batches. With the
#'   built-in terrestrial/freshwater backend, all eligible countries for one
#'   species are combined into one GBIF download; `batch_size` therefore does not
#'   represent the number of concurrent country-level GBIF jobs.
#' @param dist_km Minimum distance in kilometres used when
#'   `apply_thinning = TRUE`.
#' @param apply_thinning Logical. If `TRUE`, apply distance-based thinning after
#'   coordinate cleaning.
#' @param apply_cleaning Logical. If `TRUE`, apply coordinate-quality cleaning
#'   before export.
#' @param return_all_results Logical. If `TRUE` and `store_in_memory = TRUE`,
#'   return a combined tibble of all grouped records.
#' @param export_summary Logical. If `TRUE`, write `gbif_summary.csv`,
#'   `gbif_retrieval_audit.csv`, and `gbif_unresolved_requests.csv`.
#' @param store_in_memory Logical. If `TRUE`, retain grouped output rows in
#'   memory for return.
#' @param use_planar Logical. If `TRUE`, disable spherical geometry for selected
#'   spatial operations during this run.
#' @param use_overlays Logical. If `FALSE`, skip vector overlay assignment and
#'   export records grouped at the country/raw-occurrence level.
#' @param prepare_taxonomy Logical. If `TRUE`, run the taxonomy gate before GBIF
#'   download submission.
#' @param manual_taxonomy_fixes Optional named character vector or data frame of
#'   manual taxonomy fixes applied before automated cleaning.
#' @param taxonomy_name_col Name of the column containing submitted taxon names.
#' @param export_taxonomy_audit Logical. If `TRUE`, write taxonomy audit files to
#'   `output_dir`.
#' @param quiet Logical. If `TRUE`, reduce console messages.
#' @param deps Optional named list of dependency overrides for tests.
#' @param use_griis Logical. If `TRUE`, attach GRIIS evidence before GBIF
#'   submission.
#' @param filter_griis_invasive Logical. If `TRUE`, retain only rows flagged as
#'   invasive by GRIIS in the recipient country. Requires `use_griis = TRUE`.
#' @param griis Optional pre-read or standardised GRIIS table.
#' @param griis_cache_dir Cache directory used when GRIIS is loaded internally.
#'   If `NULL`, a subdirectory under `output_dir/_biofetchR_cache/` is used.
#' @param griis_force_refresh Logical. If `TRUE`, force GRIIS refresh when loaded
#'   internally.
#' @param griis_require_country Optional logical controlling whether GRIIS joins
#'   require recipient-country evidence. Defaults to `TRUE` in this pipeline.
#' @param use_native_range Logical. If `TRUE`, attach native-range evidence before
#'   GBIF submission.
#' @param native_ranges Native-range evidence table used by
#'   `bf_attach_native_status()`. Leave `NULL` when `use_native_web = TRUE` so
#'   the pipeline can build species-level native-origin evidence from web/API
#'   sources after taxonomy and GRIIS filtering.
#' @param use_native_web Logical. If `TRUE`, fetch native-origin evidence from
#'   package web/API helpers before applying the native-range gate. Requires
#'   `use_native_range = TRUE` and `native_ranges = NULL`.
#' @param native_web_sources Character vector of native web sources passed to
#'   `bf_fetch_native_ranges_web()`, commonly `c("sinas", "gbif", "worms")`.
#' @param native_web_cache_dir Cache directory for native web/API responses. If
#'   `NULL`, a subdirectory under `output_dir/_biofetchR_cache/` is used.
#' @param native_web_force_refresh Logical. If `TRUE`, refresh native-web cache
#'   entries where supported.
#' @param native_web_sleep_sec Numeric delay in seconds between native-web
#'   requests.
#' @param native_web_sinas_main_path Optional local path to `SInAS_3.2.csv`
#'   when `"sinas"` is included in `native_web_sources`.
#' @param native_web_sinas_alllocations_path Optional local path to
#'   `AllLocations.xlsx`, `.csv` or `.tsv`.
#' @param native_web_sinas_fulltaxa_path Optional local path to
#'   `SInAS_3.2_FullTaxaList.csv`.
#' @param export_native_web_audit Logical. If `TRUE`, write native-web long,
#'   species, unmapped and summary audit files to `output_dir`.
#' @param native_species_col Species-name column in `native_ranges`.
#' @param native_require_country Optional logical controlling whether native-range
#'   joins require recipient-country evidence. Defaults to `TRUE` in this
#'   pipeline.
#' @param native_filter_mode Native-range filter mode. Accepted values are
#'   `"audit_only"`, `"non_native_only"`, `"non_native_or_unknown"`, and
#'   `"native_only"`.
#' @param native_keep_unknown Logical. Used with `native_filter_mode =
#'   "non_native_only"`; if `TRUE`, rows with unknown origin are retained.
#' @param reconcile_origin_evidence Logical. If `TRUE`, reconcile attached GRIIS
#'   and native-range evidence into a combined origin-evidence status where the
#'   required helper is available.
#' @param export_origin_audit Logical. If `TRUE`, write origin-evidence audit
#'   files to `output_dir`.
#'
#' @return If `return_all_results = TRUE` and `store_in_memory = TRUE`, returns
#'   a tibble containing the combined processed terrestrial/freshwater occurrence
#'   records from all retained species x spatial-recipient groups. The returned
#'   table has geometry dropped and includes restored `decimalLongitude` and
#'   `decimalLatitude` columns, GBIF occurrence fields, grouping fields such as
#'   species, recipient region identifiers and region type, and any overlay,
#'   raster-context, taxonomy, GRIIS or native-origin audit columns created by
#'   the selected workflow options. If no rows remain after taxonomy or
#'   origin-evidence filtering, the returned tibble contains the stable empty
#'   coordinate columns `decimalLongitude` and `decimalLatitude`. If
#'   `return_all_results = FALSE` or `store_in_memory = FALSE`, the function
#'   writes grouped occurrence CSVs, `gbif_summary.csv` and any requested audit
#'   files to `output_dir` and returns `invisible(NULL)`. In all modes, the main
#'   side effects are GBIF download/import operations, optional raster/context
#'   enrichment, vector overlay assignment, coordinate cleaning, optional spatial
#'   thinning and writing workflow outputs.
#'
#' @family terrestrial and freshwater processing pipelines
#'
#' @seealso
#' [thin_spatial_points()], [download_gbif_batch_gadm()],
#' [bf_prepare_taxa_for_gbif()], [bf_attach_griis_status()],
#' [bf_attach_native_status()]
#'
#' @md
#' @export
process_gbif_terrestrial_freshwater_pipeline <- function(
    df, output_dir, user, pwd, email,
    # Keep the default lightweight and predictable. Users can request any of the
    # supported overlays explicitly via region_source, but the pipeline should
    # not attempt to download/load every global overlay by default.
    region_source            = "gadm",
    raster_context           = character(0),
    overlay_mode             = c("single","dual_separate","dual_intersection"),
    gadm_unit                = 1,
    cache_dir_gadm           = NULL,
    teow_cache_dir           = NULL,
    teow_method              = c("auto","mapme","direct"),
    teow_url                 = NULL,
    teow_force_refresh       = FALSE,
    teow_clip_to_countries   = TRUE,
    feow_cache_dir           = NULL,
    feow_method              = c("auto","mapme","direct"),
    feow_url                 = NULL,
    feow_force_refresh       = FALSE,
    feow_clip_to_countries   = TRUE,
    lakes_cache_dir          = NULL,
    lakes_force_refresh      = FALSE,
    lakes_only               = TRUE,
    lakes_min_area_km2       = 1,
    rivers_cache_dir         = NULL,
    rivers_cache             = c("disk","memory"),
    rivers_cache_format      = c("rds","gpkg","both"),
    rivers_force_refresh     = FALSE,
    rivers_regions           = NULL,
    rivers_min_strahler      = NULL,
    rivers_min_discharge_cms = NULL,
    rivers_max_snap_m        = 1000,
    basins_cache_dir         = NULL,
    basins_level             = 12,
    basins_with_lakes        = FALSE,
    basins_clip_to_countries = TRUE,
    gmba_cache_dir           = NULL,
    gmba_layer               = c("all","basic"),
    gmba_extent              = c("standard","broad"),
    gmba_force_refresh       = FALSE,
    gmba_clip_to_countries   = TRUE,
    ne_cache_dir             = NULL,
    ne_force_refresh         = FALSE,
    ne_admin1_scale          = "10m",
    resolve_cache_dir        = NULL,
    resolve_force_refresh    = FALSE,
    resolve_clip_to_countries= TRUE,
    wdpa_cache_dir           = NULL,
    wdpa_force_refresh       = FALSE,
    wdpa_exclude_marine      = TRUE,
    wdpa_require_opt_in      = TRUE,
    wdpa_opt_in              = FALSE,
    ramsar_cache_dir         = NULL,
    ramsar_force_refresh     = FALSE,
    ramsar_require_opt_in    = TRUE,
    ramsar_opt_in            = FALSE,
    gdw_cache_dir            = NULL,
    gdw_force_refresh        = FALSE,
    gdw_barriers_source_path = NULL,
    gdw_barriers_source_url  = NULL,
    gdw_barrier_max_snap_m   = 1000,
    biosphere_cache_dir      = NULL,
    biosphere_force_refresh  = FALSE,
    biosphere_source_path    = NULL,
    biosphere_source_url     = NULL,
    biosphere_max_snap_m     = 10000,
    gloric_cache_dir         = NULL,
    gloric_force_refresh     = FALSE,
    gloric_source_path       = NULL,
    gloric_source_url        = NULL,
    gloric_max_snap_m        = 1000,
    hydrowaste_cache_dir     = NULL,
    hydrowaste_force_refresh = FALSE,
    hydrowaste_source_path   = NULL,
    hydrowaste_source_url    = NULL,
    hydrowaste_max_snap_m    = 1000,
    gdw_reservoirs_cache_dir = NULL,
    gdw_reservoirs_force_refresh = FALSE,
    gdw_reservoirs_source_path   = NULL,
    gdw_reservoirs_source_url    = NULL,
    global_mining_cache_dir  = NULL,
    global_mining_force_refresh = FALSE,
    global_mining_source_path   = NULL,
    global_mining_source_url    = NULL,

    # Normal use remains fail-soft. Tests and development checks can set these
    # to TRUE so a requested overlay/raster-context failure is surfaced
    # immediately instead of being silently skipped.
    strict_overlay_loading   = FALSE,
    strict_raster_context    = FALSE,

    raster_cache_dir         = NULL,
    worldcover_source        = NULL,
    worldcover_year          = 2021,
    worldcover_buffer_m      = 0,
    soilgrids_sources        = NULL,
    soilgrids_vars           = NULL,
    soilgrids_buffer_m       = 0,
    human_footprint_source   = NULL,
    human_footprint_buffer_m = 0,
    batch_size               = 5,
    dist_km                  = 5,
    apply_thinning           = FALSE,
    apply_cleaning           = TRUE,
    return_all_results       = TRUE,
    export_summary           = TRUE,
    store_in_memory          = TRUE,
    use_planar               = TRUE,
    use_overlays             = TRUE,

    prepare_taxonomy         = FALSE,
    manual_taxonomy_fixes    = NULL,
    taxonomy_name_col        = "species",
    export_taxonomy_audit    = TRUE,

    quiet                    = FALSE,
    deps                     = NULL,

    # ---------------------------------------------------------------------------
    # Optional origin-evidence gates
    # ---------------------------------------------------------------------------
    use_griis = FALSE,
    filter_griis_invasive = FALSE,
    griis = NULL,
    griis_cache_dir = NULL,
    griis_force_refresh = FALSE,
    griis_require_country = NULL,

    use_native_range = FALSE,
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
    native_require_country = NULL,
    native_filter_mode = c(
      "audit_only",
      "non_native_only",
      "non_native_or_unknown",
      "native_only"
    ),
    native_keep_unknown = TRUE,

    reconcile_origin_evidence = TRUE,
    export_origin_audit = TRUE
) {

  # ---------------------------------------------------------------------------
  # Dependency injection helper
  # ---------------------------------------------------------------------------
  # `deps` is primarily for tests. Each named entry can replace a production
  # helper such as the GBIF downloader, importer, exporter, summary appender,
  # overlay loader or thinning function. This avoids live GBIF calls and large
  # remote downloads during unit tests while keeping the production API stable.
  .dep <- function(key, fallback = NULL) {
    if (!is.null(deps) && is.list(deps) && !is.null(deps[[key]]) && is.function(deps[[key]])) return(deps[[key]])
    fallback
  }

  init_summary_fun <- .dep("initialize_summary", if (exists("initialize_summary", mode="function")) initialize_summary else function() {
    tibble::tibble(
      species=character(), region_id=character(), region_type=character(),
      n_total=integer(), n_cleaned=integer(), n_thinned=integer(),
      output_file=character(), status=character()
    )
  })

  append_summary_fun <- .dep("append_summary_row", if (exists("append_summary_row", mode="function")) append_summary_row else function(summary_tbl, ...) {
    row <- list(...)
    tibble::add_row(summary_tbl, !!!row)
  })

  download_fun <- .dep("download_fun", if (exists("download_gbif_batch_gadm", mode="function")) download_gbif_batch_gadm else NULL)
  import_fun   <- .dep("import_fun", NULL)

  load_gadm_fun        <- .dep("load_gadm", if (exists("load_gadm", mode="function")) load_gadm else NULL)
  load_teow_fun        <- .dep("load_teow", if (exists("bf_load_teow", mode="function")) bf_load_teow else NULL)
  load_feow_fun        <- .dep("load_feow", if (exists("bf_load_feow", mode="function")) bf_load_feow else NULL)
  load_lakes_fun       <- .dep("load_lakes", if (exists("bf_load_lakes", mode="function")) bf_load_lakes else NULL)
  load_rivers_fun      <- .dep("load_rivers", if (exists("bf_load_rivers", mode="function")) bf_load_rivers else NULL)
  load_basins_fun      <- .dep("load_basins", if (exists("bf_load_basins", mode="function")) bf_load_basins else NULL)
  load_gmba_fun        <- .dep("load_gmba", NULL)
  load_ne_urban_fun    <- .dep("load_ne_urban", if (exists("bf_load_ne_urban", mode="function")) bf_load_ne_urban else NULL)
  load_ne_admin1_fun   <- .dep("load_ne_admin1", if (exists("bf_load_ne_admin1", mode="function")) bf_load_ne_admin1 else NULL)
  load_resolve2017_fun <- .dep("load_resolve2017", if (exists("bf_load_resolve2017", mode="function")) bf_load_resolve2017 else NULL)
  load_wdpa_fun        <- .dep("load_wdpa", if (exists("bf_load_wdpa", mode="function")) bf_load_wdpa else NULL)
  load_ramsar_fun      <- .dep("load_ramsar", if (exists("bf_load_ramsar", mode="function")) bf_load_ramsar else NULL)
  load_gdw_barriers_fun<- .dep("load_gdw_barriers", if (exists("bf_load_gdw_barriers", mode="function")) bf_load_gdw_barriers else NULL)
  load_biosphere_fun   <- .dep("load_biosphere_reserve", if (exists("bf_load_biosphere_reserve", mode="function")) bf_load_biosphere_reserve else NULL)
  load_gloric_fun      <- .dep("load_gloric", if (exists("bf_load_gloric", mode="function")) bf_load_gloric else NULL)
  load_hydrowaste_fun  <- .dep("load_hydrowaste", if (exists("bf_load_hydrowaste", mode="function")) bf_load_hydrowaste else NULL)
  load_gdw_reservoirs_fun <- .dep("load_gdw_reservoirs", if (exists("bf_load_gdw_reservoirs", mode="function")) bf_load_gdw_reservoirs else NULL)
  load_global_mining_fun <- .dep("load_global_mining", if (exists("bf_load_global_mining", mode="function")) bf_load_global_mining else NULL)
  enrich_raster_fun    <- .dep("enrich_raster_context", if (exists("bf_enrich_raster_context", mode="function")) bf_enrich_raster_context else NULL)

  thin_fun   <- .dep("thin_fun", if (exists("thin_spatial_points", mode="function")) thin_spatial_points else NULL)

  taxonomy_prepare_fun <- .dep(
    "taxonomy_prepare_fun",
    if (exists("bf_prepare_taxa_for_gbif", mode = "function")) bf_prepare_taxa_for_gbif else NULL
  )

  if (!all(c("species", "iso2c") %in% names(df))) {
    stop("`df` must contain columns: 'species' and 'iso2c'.", call. = FALSE)
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

  .bf_pipeline_cache_dir <- function(cache_dir, subdir) {
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

  cache_dir_gadm <- .bf_pipeline_cache_dir(cache_dir_gadm, "gadm")
  teow_cache_dir <- .bf_pipeline_cache_dir(teow_cache_dir, "teow")
  feow_cache_dir <- .bf_pipeline_cache_dir(feow_cache_dir, "feow")
  lakes_cache_dir <- .bf_pipeline_cache_dir(lakes_cache_dir, "lakes")
  rivers_cache_dir <- .bf_pipeline_cache_dir(rivers_cache_dir, "rivers")
  basins_cache_dir <- .bf_pipeline_cache_dir(basins_cache_dir, "basins")
  gmba_cache_dir <- .bf_pipeline_cache_dir(gmba_cache_dir, "gmba")
  ne_cache_dir <- .bf_pipeline_cache_dir(ne_cache_dir, "natural_earth")
  resolve_cache_dir <- .bf_pipeline_cache_dir(resolve_cache_dir, "resolve2017")
  wdpa_cache_dir <- .bf_pipeline_cache_dir(wdpa_cache_dir, "wdpa")
  ramsar_cache_dir <- .bf_pipeline_cache_dir(ramsar_cache_dir, "ramsar")
  gdw_cache_dir <- .bf_pipeline_cache_dir(gdw_cache_dir, "gdw_barriers")
  biosphere_cache_dir <- .bf_pipeline_cache_dir(biosphere_cache_dir, "biosphere_reserve")
  gloric_cache_dir <- .bf_pipeline_cache_dir(gloric_cache_dir, "gloric")
  hydrowaste_cache_dir <- .bf_pipeline_cache_dir(hydrowaste_cache_dir, "hydrowaste")
  gdw_reservoirs_cache_dir <- .bf_pipeline_cache_dir(gdw_reservoirs_cache_dir, "gdw_reservoirs")
  global_mining_cache_dir <- .bf_pipeline_cache_dir(global_mining_cache_dir, "global_mining")
  raster_cache_dir <- .bf_pipeline_cache_dir(raster_cache_dir, "raster_context")
  griis_cache_dir <- .bf_pipeline_cache_dir(griis_cache_dir, "griis")
  native_web_cache_dir <- .bf_pipeline_cache_dir(native_web_cache_dir, "native_web")

  sources <- unique(match.arg(region_source,
                              c("gadm","teow","feow","lakes","rivers","basins","gmba",
                                "ne_urban","resolve2017","ne_admin1","wdpa","ramsar",
                                "gdw_barriers","biosphere_reserve","gloric",
                                "hydrowaste","gdw_reservoirs",
                                "global_mining"),
                              several.ok = TRUE))
  raster_context <- unique(if (length(raster_context)) match.arg(raster_context, c("worldcover","soilgrids","human_footprint"), several.ok = TRUE) else character(0))
  overlay_mode <- match.arg(overlay_mode)
  native_filter_mode <- match.arg(native_filter_mode)
  teow_method  <- match.arg(teow_method)
  feow_method  <- match.arg(feow_method)
  rivers_cache <- match.arg(rivers_cache)
  rivers_cache_format <- match.arg(rivers_cache_format)
  gmba_layer <- match.arg(gmba_layer)
  gmba_extent <- match.arg(gmba_extent)

  if ("gadm" %in% sources && !gadm_unit %in% c(0,1,2)) stop("`gadm_unit` must be 0, 1, or 2.")

  if (!isTRUE(use_overlays)) {
    sources <- character(0)
    overlay_mode <- "single"
    if (!isTRUE(quiet)) message("biofetchR: use_overlays=FALSE -> no vector spatial joins; exporting raw occurrences.")
  }

  feow_method_resolved <- switch(feow_method, auto="auto", mapme="arcgis", direct="download", feow_method)

  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(!isTRUE(use_planar))

  df <- df |>
    dplyr::mutate(
      species = as.character(.data$species),
      iso2c   = toupper(trimws(as.character(.data$iso2c)))
    ) |>
    dplyr::filter(!(.data$species %in% c("", NA_character_)),
                  !(.data$iso2c %in% c("", NA_character_)))

  if (isTRUE(prepare_taxonomy)) {
    if (is.null(taxonomy_prepare_fun) || !is.function(taxonomy_prepare_fun)) {
      stop(
        "`prepare_taxonomy = TRUE`, but `bf_prepare_taxa_for_gbif()` is not available.",
        call. = FALSE
      )
    }

    tax_prep <- taxonomy_prepare_fun(
      df = df,
      name_col = taxonomy_name_col,
      iso2_col = "iso2c",
      manual_fixes = manual_taxonomy_fixes,
      require_species_level = TRUE,
      deduplicate = TRUE
    )

    if (isTRUE(export_taxonomy_audit)) {
      readr::write_csv(
        tax_prep$audit,
        file.path(output_dir, "taxonomy_audit.csv")
      )

      readr::write_csv(
        tax_prep$rejected,
        file.path(output_dir, "taxonomy_rejected.csv")
      )

      readr::write_csv(
        tax_prep$summary,
        file.path(output_dir, "taxonomy_summary.csv")
      )
    }

    if (!isTRUE(quiet)) {
      message(
        "biofetchR taxonomy gate: accepted ",
        nrow(tax_prep$accepted),
        " species-country row(s); rejected ",
        nrow(tax_prep$rejected),
        " row(s)."
      )
    }

    accepted_df <- tax_prep$accepted

    # User-facing downstream names must be the cleaned accepted names, not the raw
    # submitted names. Keep raw/fixed/cleaned names in taxonomy_audit.csv, but make
    # df$species clean before GBIF submission, summaries, filenames and console logs.
    if ("cleaned_name" %in% names(accepted_df)) {
      accepted_df$species <- accepted_df$cleaned_name
    } else if ("accepted_name" %in% names(accepted_df)) {
      accepted_df$species <- accepted_df$accepted_name
    } else if ("scientificName" %in% names(accepted_df)) {
      accepted_df$species <- accepted_df$scientificName
    }

    df <- accepted_df |>
      dplyr::mutate(
        species = as.character(.data$species),
        iso2c = toupper(trimws(as.character(.data$iso2c)))
      ) |>
      dplyr::filter(
        !is.na(.data$species),
        nzchar(.data$species),
        !is.na(.data$iso2c),
        nzchar(.data$iso2c)
      ) |>
      dplyr::distinct(.data$species, .data$iso2c, .keep_all = TRUE)

    if (!nrow(df)) {
      if (isTRUE(export_summary)) {
        readr::write_csv(
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
          ),
          file.path(output_dir, "gbif_summary.csv")
        )
      }

      if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
        return(tibble::tibble(decimalLongitude = numeric(), decimalLatitude = numeric()))
      }

      return(invisible(NULL))
    }
  }

  # ---------------------------------------------------------------------------
  # Optional GRIIS / native-range origin-evidence gates
  # ---------------------------------------------------------------------------
  # These gates must remain here: after taxonomy preparation and before species
  # queue construction / GBIF download submission. Rows rejected here should not
  # reach download_gbif_batch_gadm().

  if (isTRUE(filter_griis_invasive) && !isTRUE(use_griis)) {
    stop("`filter_griis_invasive = TRUE` requires `use_griis = TRUE`.", call. = FALSE)
  }

  if (isTRUE(use_native_web) && !isTRUE(use_native_range)) {
    stop(
      "`use_native_web = TRUE` requires `use_native_range = TRUE`, because web-derived native-origin evidence is only used by the native-range gate.",
      call. = FALSE
    )
  }

  if (isTRUE(use_native_web) && !is.null(native_ranges)) {
    stop(
      "`use_native_web = TRUE` currently requires `native_ranges = NULL`. Use either supplied native-range evidence or web-derived evidence, not both.",
      call. = FALSE
    )
  }

  if (is.null(griis_require_country)) {
    griis_require_country <- TRUE
  }

  if (is.null(native_require_country)) {
    native_require_country <- TRUE
  }

  .origin_iso2_col <- if ("iso2c" %in% names(df)) "iso2c" else NULL
  .origin_iso3_col <- if ("iso3c" %in% names(df)) "iso3c" else NULL
  .origin_country_col <- if ("country" %in% names(df)) "country" else NULL

  if (isTRUE(use_griis) || isTRUE(use_native_range)) {
    df$.bf_origin_rowid <- seq_len(nrow(df))

    if (isTRUE(use_griis)) {
      if (!exists("bf_attach_griis_status", mode = "function")) {
        stop(
          "`use_griis = TRUE`, but `bf_attach_griis_status()` is not available. ",
          "Make sure `R/utils_griis.R` is present and loaded by the package.",
          call. = FALSE
        )
      }

      df <- bf_attach_griis_status(
        df = df,
        species_col = "species",
        iso2c_col = .origin_iso2_col,
        iso3c_col = .origin_iso3_col,
        country_col = .origin_country_col,
        griis = griis,
        cache_dir = griis_cache_dir,
        force_refresh = griis_force_refresh,
        require_country = isTRUE(griis_require_country),
        quiet = quiet
      )
    }

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

      native_ranges <- native_web$species
      native_species_col <- "species"

      if (!isTRUE(quiet)) {
        message(
          "biofetchR native web evidence: fetched species-level origin evidence for ",
          length(species_for_native_web),
          " species."
        )
      }
    }

    if (isTRUE(use_native_range)) {
      if (!exists("bf_attach_native_status", mode = "function")) {
        stop(
          "`use_native_range = TRUE`, but `bf_attach_native_status()` is not available. ",
          "Make sure `R/utils_native_range.R` is present and loaded by the package.",
          call. = FALSE
        )
      }

      .native_filter_mode_for_helper <- switch(
        native_filter_mode,
        audit_only = "audit_only",
        non_native_only = if (isTRUE(native_keep_unknown)) "keep_non_native_or_unknown" else "keep_non_native",
        non_native_or_unknown = "keep_non_native_or_unknown",
        native_only = "keep_native"
      )

      df <- bf_attach_native_status(
        df = df,
        species_col = "species",
        iso2c_col = .origin_iso2_col,
        iso3c_col = .origin_iso3_col,
        country_col = .origin_country_col,
        native_ranges = native_ranges,
        native_species_col = native_species_col,
        require_country = isTRUE(native_require_country),
        native_filter_mode = .native_filter_mode_for_helper,
        quiet = quiet
      )
    }

    if (isTRUE(reconcile_origin_evidence) && exists("bf_reconcile_griis_native_status", mode = "function")) {
      has_griis_cols <- all(c("griis_listed", "griis_invasive") %in% names(df))
      has_native_cols <- all(c("native_status", "native_is_native_recipient", "native_is_non_native_recipient") %in% names(df))

      if (has_griis_cols || has_native_cols) {
        df <- bf_reconcile_griis_native_status(df)
      }
    }

    .origin_keep <- rep(TRUE, nrow(df))

    if (isTRUE(use_griis) && isTRUE(filter_griis_invasive)) {
      .origin_keep <- .origin_keep & (df$griis_invasive %in% TRUE)
    }

    if (isTRUE(use_native_range) && native_filter_mode != "audit_only") {
      if (!"native_filter_keep" %in% names(df)) {
        stop(
          "Native-range filtering was requested, but `native_filter_keep` was not created.",
          call. = FALSE
        )
      }

      .origin_keep <- .origin_keep & (df$native_filter_keep %in% TRUE)
    }

    .origin_rejection_reason <- rep(NA_character_, nrow(df))

    if (isTRUE(use_griis) && isTRUE(filter_griis_invasive)) {
      .origin_rejection_reason[!(df$griis_invasive %in% TRUE)] <- "not_griis_invasive"
    }

    if (isTRUE(use_native_range) && native_filter_mode != "audit_only") {
      idx <- is.na(.origin_rejection_reason) & !(df$native_filter_keep %in% TRUE)
      .origin_rejection_reason[idx] <- paste0("native_filter_", native_filter_mode)
    }

    origin_audit <- df
    origin_audit$origin_gate_keep <- .origin_keep
    origin_audit$origin_gate_rejection_reason <- .origin_rejection_reason

    origin_rejected <- origin_audit[!origin_audit$origin_gate_keep, , drop = FALSE]

    origin_summary <- tibble::tibble(
      n_input = nrow(origin_audit),
      n_retained = sum(origin_audit$origin_gate_keep, na.rm = TRUE),
      n_rejected = nrow(origin_rejected),
      use_griis = isTRUE(use_griis),
      filter_griis_invasive = isTRUE(filter_griis_invasive),
      use_native_range = isTRUE(use_native_range),
      use_native_web = isTRUE(use_native_web),
      native_web_sources = if (isTRUE(use_native_web)) paste(native_web_sources, collapse = ";") else NA_character_,
      native_filter_mode = native_filter_mode,
      native_keep_unknown = isTRUE(native_keep_unknown),
      n_griis_listed = if ("griis_listed" %in% names(origin_audit)) sum(origin_audit$griis_listed, na.rm = TRUE) else NA_integer_,
      n_griis_invasive = if ("griis_invasive" %in% names(origin_audit)) sum(origin_audit$griis_invasive, na.rm = TRUE) else NA_integer_,
      n_native_in_recipient = if ("native_is_native_recipient" %in% names(origin_audit)) sum(origin_audit$native_is_native_recipient, na.rm = TRUE) else NA_integer_,
      n_non_native_in_recipient = if ("native_is_non_native_recipient" %in% names(origin_audit)) sum(origin_audit$native_is_non_native_recipient, na.rm = TRUE) else NA_integer_,
      n_origin_unknown = if ("native_status" %in% names(origin_audit)) sum(origin_audit$native_status == "origin_unknown", na.rm = TRUE) else NA_integer_
    )

    if (isTRUE(export_origin_audit)) {
      readr::write_csv(origin_audit, file.path(output_dir, "origin_evidence_audit.csv"))
      readr::write_csv(origin_rejected, file.path(output_dir, "origin_evidence_rejected.csv"))
      readr::write_csv(origin_summary, file.path(output_dir, "origin_evidence_summary.csv"))
    }

    if (!isTRUE(quiet)) {
      message(
        "biofetchR origin gate: retained ",
        sum(.origin_keep, na.rm = TRUE),
        " row(s); rejected ",
        nrow(origin_rejected),
        " row(s)."
      )
    }

    df <- df[.origin_keep, , drop = FALSE]
    df$origin_gate_keep <- NULL
    df$origin_gate_rejection_reason <- NULL
    df$.bf_origin_rowid <- NULL

    if (!nrow(df)) {
      if (isTRUE(export_summary)) {
        readr::write_csv(
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
          ),
          file.path(output_dir, "gbif_summary.csv")
        )
      }

      if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
        return(tibble::tibble(decimalLongitude = numeric(), decimalLatitude = numeric()))
      }

      return(invisible(NULL))
    }
  }

  species_metadata <- df |>
    dplyr::distinct(.data$species, .data$iso2c, .keep_all = TRUE)

  species_metadata_cols <- setdiff(names(species_metadata), c("species", "iso2c"))

  species_list <- unique(df$species)
  batches      <- split(species_list, ceiling(seq_along(species_list) / batch_size))
  countries    <- unique(df$iso2c)

  summary_tbl <- init_summary_fun()
  if (!"fail_stage" %in% names(summary_tbl)) summary_tbl$fail_stage <- NA_character_
  if (!"fail_reason" %in% names(summary_tbl)) summary_tbl$fail_reason <- NA_character_
  if (!"gbif_key" %in% names(summary_tbl)) summary_tbl$gbif_key <- NA_character_
  all_results <- list()

  # One row per taxonomy-approved species x country retrieval request. This is
  # deliberately separate from gbif_summary.csv because downstream overlay
  # outputs can have multiple rows per original country request.
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

  # ---------------------------------------------------------------------------
  # Local status, summary and error helpers
  # ---------------------------------------------------------------------------
  # These helpers are local rather than package-level because they rely on the
  # current function arguments, summary table and overlay settings. They keep the
  # main loop readable and ensure that fail-soft branches still write a useful
  # summary row.
  .msg    <- function(...) { if (!isTRUE(quiet)) message(...) }
  .errmsg <- function(e) if (inherits(e, "condition")) conditionMessage(e) else as.character(e)

  .append_row <- function(species, region_id, region_type,
                          n_total=0, n_cleaned=0, n_thinned=0,
                          output_file=NA_character_, status="unknown",
                          fail_stage=NA_character_, fail_reason=NA_character_, gbif_key=NA_character_) {
    summary_tbl <<- append_summary_fun(
      summary_tbl,
      species=species, region_id=region_id, region_type=region_type,
      n_total=n_total, n_cleaned=n_cleaned, n_thinned=n_thinned,
      output_file=output_file, status=status
    )
    i <- nrow(summary_tbl)
    summary_tbl$fail_stage[i] <<- fail_stage
    summary_tbl$fail_reason[i] <<- fail_reason
    summary_tbl$gbif_key[i] <<- gbif_key
    invisible(NULL)
  }

  .append_retrieval_audit <- function(species,
                                      request_id,
                                      gbif_key = NA_character_,
                                      gbif_status = NA_character_,
                                      retrieval_status,
                                      n_records = NA_integer_,
                                      fail_reason = NA_character_) {
    species <- as.character(species[[1L]])
    request_id <- as.character(request_id[[1L]])

    # Upsert rather than append blindly: every species x country request should
    # have exactly one final retrieval outcome.
    if (nrow(retrieval_audit)) {
      keep <- !(
        retrieval_audit$species == species &
          retrieval_audit$request_id == request_id
      )
      retrieval_audit <<- retrieval_audit[keep, , drop = FALSE]
    }

    retrieval_audit <<- dplyr::bind_rows(
      retrieval_audit,
      tibble::tibble(
        species = species,
        request_id = request_id,
        request_type = "country",
        gbif_key = as.character(gbif_key[[1L]]),
        gbif_status = as.character(gbif_status[[1L]]),
        retrieval_status = as.character(retrieval_status[[1L]]),
        n_records = as.integer(n_records[[1L]]),
        fail_reason = as.character(fail_reason[[1L]])
      )
    )

    invisible(NULL)
  }


  .region_type_label <- function(src) {
    switch(as.character(src),
           gadm              = paste0("GADM_level_", gadm_unit),
           teow              = "TEOW",
           feow              = "FEOW",
           lakes             = "LAKES",
           rivers            = "RIVERS",
           basins            = paste0("BASINS_L", basins_level),
           gmba              = paste0("GMBA_", gmba_layer, "_", gmba_extent),
           ne_urban          = "NE_URBAN",
           resolve2017       = "RESOLVE2017",
           ne_admin1         = "NE_ADMIN1",
           wdpa              = "WDPA",
           ramsar            = "RAMSAR",
           gdw_barriers      = "GDW_BARRIERS",
           hydrowaste        = "HYDROWASTE",
           gdw_reservoirs    = "GDW_RESERVOIRS",
           global_mining     = "GLOBAL_MINING",
           biosphere_reserve = "BIOSPHERE_RESERVE",
           gloric            = "GLORIC",
           toupper(as.character(src)))
  }
  .region_type_combo <- function(src_vec, sep = "+") paste(vapply(src_vec, .region_type_label, character(1)), collapse = sep)
  .region_type_for_fail <- function() if (isTRUE(use_overlays) && length(sources)) .region_type_combo(sources, sep = "+") else "RAW_COUNTRY"

  .write_summary_now <- function() {
    if (isTRUE(export_summary)) {
      try(
        readr::write_csv(summary_tbl, file.path(output_dir, "gbif_summary.csv")),
        silent = TRUE
      )
    }
    invisible(NULL)
  }

  .record_overlay_load_failure <- function(src,
                                           reason,
                                           stage = "load_overlay") {
    .append_row(
      species = paste0("__overlay_loader__", src),
      region_id = paste(countries, collapse = ";"),
      region_type = .region_type_label(src),
      status = "failed",
      fail_stage = stage,
      fail_reason = reason
    )

    if (isTRUE(strict_overlay_loading)) {
      .write_summary_now()
      stop("Requested overlay failed to load: ", src, " - ", reason, call. = FALSE)
    }

    invisible(NULL)
  }


  .load_overlay_or_record <- function(src,
                                      fun,
                                      args,
                                      crop_after = FALSE,
                                      crop_tag = NULL) {
    if (is.null(fun) || !is.function(fun)) {
      .record_overlay_load_failure(
        src,
        "Loader function is not available.",
        stage = "missing_loader"
      )
      return(NULL)
    }

    out <- tryCatch(
      do.call(fun, args),
      error = function(e) {
        .record_overlay_load_failure(
          src,
          .errmsg(e),
          stage = "load_overlay"
        )
        NULL
      }
    )

    if (is.null(out)) {
      return(NULL)
    }

    if (!inherits(out, "sf") || !nrow(out)) {
      .record_overlay_load_failure(
        src,
        "Loader returned no sf rows.",
        stage = "empty_overlay"
      )
      return(NULL)
    }

    out <- tryCatch(
      .make_valid(.ensure_sf_wgs84(out)),
      error = function(e) {
        .record_overlay_load_failure(
          src,
          .errmsg(e),
          stage = "validate_overlay"
        )
        NULL
      }
    )

    if (is.null(out)) {
      return(NULL)
    }

    if (isTRUE(crop_after)) {
      out <- .crop_overlay_to_countries(
        out,
        countries,
        tag = if (is.null(crop_tag)) toupper(src) else crop_tag
      )
    }

    if (!inherits(out, "sf") || !nrow(out)) {
      .record_overlay_load_failure(
        src,
        "Overlay became empty after validation/cropping.",
        stage = "empty_overlay"
      )
      return(NULL)
    }

    out
  }

  .prepare_overlay_or_record <- function(src, out, fun) {
    if (is.null(out)) {
      return(NULL)
    }

    prepared <- tryCatch(
      fun(out),
      error = function(e) {
        .record_overlay_load_failure(
          src,
          .errmsg(e),
          stage = "prepare_overlay"
        )
        NULL
      }
    )

    if (is.null(prepared)) {
      return(NULL)
    }

    if (!inherits(prepared, "sf") || !nrow(prepared)) {
      .record_overlay_load_failure(
        src,
        "Overlay preparation returned no sf rows.",
        stage = "empty_overlay"
      )
      return(NULL)
    }

    prepared
  }
  if (requireNamespace("countrycode", quietly = TRUE)) {
    iso2_all <- countries
    iso2_all[iso2_all == "UK"] <- "GB"
    iso3 <- suppressWarnings(countrycode::countrycode(
      iso2_all, "iso2c", "iso3c",
      custom_match = c("XK" = "XKX", "NA" = "NAM"),
      warn = FALSE
    ))
    bad <- is.na(iso3) | !nzchar(iso3)
    if (any(bad)) {
      bad_iso2 <- unique(countries[bad])
      bad_rows <- df |> dplyr::filter(.data$iso2c %in% bad_iso2)
      for (i in seq_len(nrow(bad_rows))) {
        .append_row(
          species = as.character(bad_rows$species[i]),
          region_id = as.character(bad_rows$iso2c[i]),
          region_type = .region_type_for_fail(),
          status = "failed",
          fail_stage = "invalid_iso2c",
          fail_reason = "ISO2 code not recognised; skipping this country."
        )
      }
      df <- df |> dplyr::filter(!(.data$iso2c %in% bad_iso2))
      countries <- unique(df$iso2c)
      if (!length(countries) || nrow(df) == 0) {
        if (isTRUE(export_summary)) readr::write_csv(summary_tbl, file.path(output_dir, "gbif_summary.csv"))
        if (isTRUE(return_all_results) && isTRUE(store_in_memory)) return(tibble::tibble(decimalLongitude = numeric(), decimalLatitude = numeric()))
        return(invisible(NULL))
      }
    }
  }

  pipeline_start_time <- bf_console_now()

  species_vec_console <- unique(as.character(df$species))

  bf_console_pipeline_start(
    workflow = "terrestrial/freshwater",
    n_species = length(species_vec_console),
    region_source = sources,
    gadm_unit = gadm_unit,
    overlay_mode = overlay_mode,
    cleaning = apply_cleaning,
    thinning = apply_thinning,
    quiet = quiet
  )

  # ---------------------------------------------------------------------------
  # Local spatial-normalisation helpers
  # ---------------------------------------------------------------------------
  # These are intentionally kept inside the pipeline because they depend on the
  # selected geometry mode, output schema and export-repair rules used here.
  # They standardise CRS handling, geometry validity, longitude/latitude columns
  # and mixed-schema binding before records are joined to overlays or exported.
  .make_valid <- function(x) {
    if (!inherits(x,"sf")) return(x)
    if ("st_make_valid" %in% getNamespaceExports("sf")) return(suppressWarnings(sf::st_make_valid(x)))
    suppressWarnings(sf::st_buffer(x, 0))
  }

  .geom_name <- function(x) {
    nm <- tryCatch(attr(x, "sf_column"), error=function(e) NULL)
    if (!is.null(nm) && nm %in% names(x)) return(nm)
    cand <- intersect(c("geometry","geom","wkb_geometry"), names(x))
    if (length(cand)) { sf::st_geometry(x) <- cand[1]; return(cand[1]) }
    stop("No geometry column found.")
  }

  .ensure_sf_wgs84 <- function(x) {
    if (is.null(x) || !nrow(x)) return(x)
    if (!inherits(x, "sf") && all(c("decimalLongitude","decimalLatitude") %in% names(x))) {
      x <- sf::st_as_sf(x, coords=c("decimalLongitude","decimalLatitude"), crs=4326, remove=FALSE)
    }
    crs_now <- tryCatch(sf::st_crs(x), error=function(e) NA)
    if (is.na(crs_now) || is.null(crs_now)) x <- sf::st_set_crs(x, 4326)
    if (!is.na(sf::st_crs(x)) && sf::st_crs(x) != sf::st_crs(4326)) x <- sf::st_transform(x, 4326)
    x <- x[!sf::st_is_empty(x), , drop=FALSE]
    x
  }

  .ensure_lonlat_cols <- function(x) {
    if (!inherits(x,"sf")) return(x)
    x <- .ensure_sf_wgs84(x)
    if (all(c("decimalLongitude","decimalLatitude") %in% names(x))) {
      x$decimalLongitude <- suppressWarnings(as.numeric(x$decimalLongitude))
      x$decimalLatitude <- suppressWarnings(as.numeric(x$decimalLatitude))
      return(x)
    }
    coords <- tryCatch(sf::st_coordinates(x), error=function(e) NULL)
    if (is.null(coords) || !nrow(coords)) return(x)
    x$decimalLongitude <- as.numeric(coords[,1])
    x$decimalLatitude <- as.numeric(coords[,2])
    x
  }

  .strip_internal_cols <- function(df) {
    if (is.null(df) || !is.data.frame(df)) return(df)
    keep <- !startsWith(names(df), ".")
    df[, keep, drop = FALSE]
  }

  .restore_coords_df <- function(out_df, proc_sf) {
    need <- c("decimalLongitude","decimalLatitude")
    if (all(need %in% names(out_df))) return(out_df)
    if (".bf_rowid" %in% names(out_df) && ".bf_rowid" %in% names(proc_sf)) {
      key <- proc_sf |>
        sf::st_drop_geometry() |>
        dplyr::select(.bf_rowid, decimalLongitude, decimalLatitude)
      out_df <- dplyr::left_join(out_df, key, by = ".bf_rowid")
      return(out_df)
    }
    if (nrow(out_df) == nrow(proc_sf)) {
      out_df$decimalLongitude <- sf::st_drop_geometry(proc_sf)$decimalLongitude
      out_df$decimalLatitude <- sf::st_drop_geometry(proc_sf)$decimalLatitude
    }
    out_df
  }

  .restore_coords_file <- function(path, proc_sf) {
    if (!file.exists(path)) return(path)
    df_disk <- tryCatch(readr::read_csv(path, show_col_types = FALSE), error=function(e) NULL)
    if (is.null(df_disk)) return(path)
    df_fix <- .restore_coords_df(df_disk, proc_sf)
    df_fix <- .strip_internal_cols(df_fix)
    if (!all(c("decimalLongitude","decimalLatitude") %in% names(df_fix))) return(path)
    tmp <- paste0(path, ".tmp")
    tryCatch({
      readr::write_csv(df_fix, tmp)
      file.rename(tmp, path)
    }, error=function(e) try(unlink(tmp), silent=TRUE))
    path
  }

  .clean_and_thin <- function(sf_df) {
    n_total <- nrow(sf_df)
    cleaned <- sf_df
    n_cleaned <- NA_integer_
    if (isTRUE(apply_cleaning) && is.function(thin_fun)) {
      tmp <- tryCatch(thin_fun(sf_df, dist_km = 0, filter_uncertain = TRUE, quiet = TRUE), error=function(e) NULL)
      if (!is.null(tmp)) { cleaned <- tmp; n_cleaned <- nrow(tmp) }
    }
    final <- cleaned
    if (isTRUE(apply_thinning) && is.finite(dist_km) && dist_km > 0 && is.function(thin_fun)) {
      tmp2 <- tryCatch(thin_fun(final, dist_km = dist_km, filter_uncertain = FALSE, quiet = TRUE), error=function(e) NULL)
      if (!is.null(tmp2)) final <- tmp2
    }
    list(sf = final, n_total = n_total, n_cleaned = n_cleaned, n_thinned = nrow(final))
  }

  .bf_bind_rows_safe <- function(lst) {
    lst <- lapply(lst, function(x) tibble::as_tibble(x))
    all_cols <- unique(unlist(lapply(lst, names)))
    lst <- lapply(lst, function(df) {
      miss <- setdiff(all_cols, names(df))
      if (length(miss)) for (m in miss) df[[m]] <- NA
      df[, all_cols, drop = FALSE]
    })
    cls_map <- lapply(all_cols, function(col) unique(vapply(lst, function(df) class(df[[col]])[1], character(1))))
    names(cls_map) <- all_cols
    needs_chr <- vapply(cls_map, function(cl) {
      cl <- cl[!is.na(cl)]
      length(unique(cl)) > 1 || any(cl %in% c("integer64", "bit64"))
    }, logical(1))
    cols_to_chr <- names(needs_chr)[needs_chr]
    lst <- lapply(lst, function(df) {
      for (nm in names(df)) if (is.factor(df[[nm]])) df[[nm]] <- as.character(df[[nm]])
      if (length(cols_to_chr)) for (col in cols_to_chr) df[[col]] <- as.character(df[[col]])
      df
    })
    dplyr::bind_rows(lst)
  }

  .join_poly_label <- function(pts_sf, poly_sf, label_col) {
    if (is.null(poly_sf) || !inherits(poly_sf,"sf") || !nrow(poly_sf)) return(rep(NA_character_, nrow(pts_sf)))
    if (!label_col %in% names(poly_sf)) return(rep(NA_character_, nrow(pts_sf)))
    pts <- .ensure_sf_wgs84(pts_sf)
    poly <- .ensure_sf_wgs84(poly_sf)
    idx <- tryCatch(suppressWarnings(sf::st_intersects(pts, poly)), error=function(e) NULL)
    if (is.null(idx)) return(rep(NA_character_, nrow(pts_sf)))
    vals <- as.character(poly[[label_col]])
    out <- rep(NA_character_, length(idx))
    for (i in seq_along(idx)) {
      hits <- idx[[i]]
      if (length(hits)) out[i] <- vals[hits[1]]
    }
    out
  }

  .join_gadm <- function(pts_sf, gadm_sf) {
    out <- list(id = rep(NA_character_, nrow(pts_sf)), name = rep(NA_character_, nrow(pts_sf)), geoname = rep(NA_character_, nrow(pts_sf)))
    if (is.null(gadm_sf) || !inherits(gadm_sf, "sf") || !nrow(gadm_sf)) return(out)
    pts <- .ensure_sf_wgs84(pts_sf)
    poly <- .ensure_sf_wgs84(gadm_sf)
    if (!all(c("id_gadm","name_gadm","geoname_gadm") %in% names(poly))) stop("GADM overlay missing required columns: id_gadm / name_gadm / geoname_gadm")
    vals_id <- as.character(poly$id_gadm)
    vals_name <- as.character(poly$name_gadm)
    vals_geo <- as.character(poly$geoname_gadm)
    first_hit <- function(idx, vals) {
      outv <- rep(NA_character_, length(idx))
      for (i in seq_along(idx)) {
        h <- idx[[i]]
        if (length(h)) outv[i] <- vals[h[1]]
      }
      outv
    }
    idx1 <- tryCatch(suppressWarnings(sf::st_within(pts, poly)), error = function(e) NULL)
    if (!is.null(idx1)) {
      out$id <- first_hit(idx1, vals_id)
      out$name <- first_hit(idx1, vals_name)
      out$geoname <- first_hit(idx1, vals_geo)
    }
    need <- is.na(out$id)
    if (any(need)) {
      idx2 <- tryCatch(suppressWarnings(sf::st_covered_by(pts[need, , drop=FALSE], poly)), error = function(e) NULL)
      if (!is.null(idx2)) {
        out$id[need] <- first_hit(idx2, vals_id)
        out$name[need] <- first_hit(idx2, vals_name)
        out$geoname[need] <- first_hit(idx2, vals_geo)
        need <- is.na(out$id)
      }
    }
    if (any(need)) {
      idx3 <- tryCatch(suppressWarnings(sf::st_intersects(pts[need, , drop=FALSE], poly)), error = function(e) NULL)
      if (!is.null(idx3)) {
        out$id[need] <- first_hit(idx3, vals_id)
        out$name[need] <- first_hit(idx3, vals_name)
        out$geoname[need] <- first_hit(idx3, vals_geo)
      }
    }
    out
  }

  .join_nearest_id <- function(pts_sf, target_sf, id_col, snap_m = 1000) {
    if (is.null(target_sf) || !inherits(target_sf, "sf") || !nrow(target_sf)) return(rep(NA_character_, nrow(pts_sf)))
    if (!id_col %in% names(target_sf)) return(rep(NA_character_, nrow(pts_sf)))
    pts <- .ensure_sf_wgs84(pts_sf)
    tgt <- .ensure_sf_wgs84(target_sf)
    pts_m <- tryCatch(sf::st_transform(pts, 3857), error=function(e) pts)
    tgt_m <- tryCatch(sf::st_transform(tgt, 3857), error=function(e) tgt)
    idx <- tryCatch(sf::st_nearest_feature(pts_m, tgt_m), error=function(e) NULL)
    if (is.null(idx)) return(rep(NA_character_, nrow(pts_sf)))
    d <- tryCatch(suppressWarnings(sf::st_distance(sf::st_geometry(pts_m), sf::st_geometry(tgt_m[idx, ]), by_element = TRUE)), error=function(e) NULL)
    id <- as.character(tgt[[id_col]][idx])
    if (!is.null(d)) id[as.numeric(d) > snap_m] <- NA_character_
    id
  }

  .join_nearest_attr <- function(pts_sf, target_sf, attr_col, snap_m = 1000) {
    if (is.null(target_sf) || !inherits(target_sf, "sf") || !nrow(target_sf)) return(rep(NA_character_, nrow(pts_sf)))
    if (!attr_col %in% names(target_sf)) return(rep(NA_character_, nrow(pts_sf)))
    pts <- .ensure_sf_wgs84(pts_sf)
    tgt <- .ensure_sf_wgs84(target_sf)
    pts_m <- tryCatch(sf::st_transform(pts, 3857), error=function(e) pts)
    tgt_m <- tryCatch(sf::st_transform(tgt, 3857), error=function(e) tgt)
    idx <- tryCatch(sf::st_nearest_feature(pts_m, tgt_m), error=function(e) NULL)
    if (is.null(idx)) return(rep(NA_character_, nrow(pts_sf)))
    d <- tryCatch(suppressWarnings(sf::st_distance(sf::st_geometry(pts_m), sf::st_geometry(tgt_m[idx, ]), by_element = TRUE)), error=function(e) NULL)
    val <- as.character(target_sf[[attr_col]][idx])
    if (!is.null(d)) val[as.numeric(d) > snap_m] <- NA_character_
    val
  }

  .crop_overlay_to_countries <- function(overlay, iso2, tag="overlay") {
    if (is.null(overlay) || !inherits(overlay,"sf") || !nrow(overlay) || !length(iso2)) return(overlay)
    tryCatch({
      if (!requireNamespace("rnaturalearth", quietly=TRUE)) return(overlay)
      world <- rnaturalearth::ne_countries(scale="small", returnclass="sf")
      sel <- world[world$iso_a2 %in% toupper(iso2), , drop=FALSE]
      if (!nrow(sel)) return(overlay)
      bb <- sf::st_as_sfc(sf::st_bbox(sf::st_transform(sel, 4326)))
      pre <- nrow(overlay)
      cropped <- suppressWarnings(sf::st_crop(.ensure_sf_wgs84(overlay), sf::st_bbox(bb)))
      if (!isTRUE(quiet)) message(tag, " crop: ", nrow(cropped), " / ", pre, " features.")
      if (nrow(cropped) > 0) cropped else overlay
    }, error=function(e) overlay)
  }

  .load_gmba_inventory <- function(cache_dir, layer, extent, force_refresh, quiet_dl=TRUE) {
    bf_require_packages(c("sf"), context = "terrestrial/freshwater pipeline")

    if (is.function(load_gmba_fun)) {
      out <- load_gmba_fun(
        cache_dir = cache_dir,
        layer = layer,
        extent = extent,
        force_refresh = force_refresh,
        quiet = quiet_dl
      )
      if (!inherits(out, "sf") || !nrow(out)) {
        stop("External GMBA loader returned no sf rows.", call. = FALSE)
      }
      return(out)
    }

    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive=TRUE, showWarnings=FALSE)

    layer  <- match.arg(layer,  c("all", "basic"))
    extent <- match.arg(extent, c("standard", "broad"))

    base <- "https://data.earthenv.org/mountains"
    fname <- switch(
      paste(layer, extent, sep="_"),
      "all_standard"   = "GMBA_Inventory_v2.0_standard.zip",
      "all_broad"      = "GMBA_Inventory_v2.0_broad.zip",
      "basic_standard" = "GMBA_Inventory_v2.0_standard_basic.zip",
      "basic_broad"    = "GMBA_Inventory_v2.0_broad_basic.zip"
    )

    if (is.null(fname) || !nzchar(fname)) {
      stop("Unsupported GMBA layer/extent combination: ", layer, " / ", extent, call. = FALSE)
    }

    url <- paste(base, extent, fname, sep="/")
    tgt <- file.path(cache_dir, fname)
    exdir <- file.path(cache_dir, tools::file_path_sans_ext(fname))

    bf_download_cached(
      url = url,
      dest = tgt,
      force_refresh = force_refresh,
      quiet = quiet_dl,
      min_bytes = 1024,
      validate_not_html = TRUE
    )

    bf_unpack_archive(
      path = tgt,
      exdir = exdir,
      force_refresh = force_refresh,
      quiet = quiet_dl
    )

    shp <- list.files(exdir, pattern="\\.shp$", recursive=TRUE, full.names=TRUE)
    if (!length(shp)) stop("GMBA download completed but no .shp file was found in: ", exdir, call. = FALSE)

    # Prefer the main inventory shapefile if sidecar/admin layers are present.
    shp_base <- basename(shp)
    score <- as.integer(grepl("GMBA|Inventory|mountain", shp_base, ignore.case = TRUE))
    shp <- shp[order(score, decreasing = TRUE)]

    gmba <- sf::st_read(shp[1], quiet=TRUE)

    id_cand <- intersect(
      tolower(names(gmba)),
      tolower(c("GMBA_ID", "GMBA_V2_ID", "GMBA_V2", "ID", "RANGE_ID", "ID_", "gmba_id", "id"))
    )
    if (length(id_cand)) {
      orig <- names(gmba)[match(id_cand[1], tolower(names(gmba)))]
      gmba$gmba_id <- as.character(gmba[[orig]])
    } else {
      gmba$gmba_id <- as.character(seq_len(nrow(gmba)))
    }

    gcol <- .geom_name(gmba)
    gmba <- gmba[, c("gmba_id", gcol), drop=FALSE]
    .make_valid(.ensure_sf_wgs84(gmba))
  }

  .gbif_has_geo_records <- function(species, iso2) {
    if (identical(iso2, "NA")) return(TRUE)
    if (!requireNamespace("rgbif", quietly=TRUE)) return(TRUE)
    n <- tryCatch({
      res <- rgbif::occ_search(scientificName=species, country=iso2, hasCoordinate=TRUE, limit=0)
      as.integer(res$meta$count)
    }, error=function(e) NA_integer_)
    if (is.na(n)) return(TRUE)
    n > 0L
  }

  # ---------------------------------------------------------------------------
  # Overlay label configuration
  # ---------------------------------------------------------------------------
  # Each overlay exposes a different identifier column. The label map below
  # defines the column used to split joined occurrences into exported groups.
  # Point and line overlays are handled by nearest-feature helpers, while polygon
  # overlays are handled by containment/intersection helpers.
  label_cols <- list(
    gadm              = "id_gadm",
    teow              = "geoname_teow",
    feow              = "geoname_feow",
    lakes             = "geoname_lakes",
    rivers            = "hyriv_id",
    basins            = "hybas_id",
    gmba              = "gmba_id",
    ne_urban          = "geoname_ne_urban",
    resolve2017       = "geoname_resolve2017",
    ne_admin1         = "id_ne_admin1",
    wdpa              = "wdpa_id",
    ramsar            = "ramsar_id",
    gdw_barriers      = "gdw_barrier_id",
    biosphere_reserve = "biosphere_id",
    gloric            = "gloric_id",
    hydrowaste        = "hydrowaste_id",
    gdw_reservoirs    = "gdw_reservoir_id",
    global_mining     = "global_mining_id"
  )

  overlay_gadm <- overlay_teow <- overlay_feow <- overlay_lakes <- overlay_rivers <- overlay_basins <- overlay_gmba <- NULL
  overlay_ne_urban <- overlay_resolve2017 <- overlay_ne_admin1 <- overlay_wdpa <- overlay_ramsar <- NULL
  overlay_gdw_barriers <- overlay_biosphere <- overlay_gloric <- NULL
  overlay_hydrowaste <- NULL
  overlay_gdw_reservoirs <- overlay_global_mining <- NULL

  if (isTRUE(use_overlays) && length(sources)) {
    bf_console_bullet(
      paste0(
        "Loading terrestrial overlay resource(s): ",
        paste(sources, collapse = ", ")
      ),
      quiet = quiet
    )
  }

  if (isTRUE(use_overlays)) {
    if ("gadm" %in% sources) {
      if (is.null(load_gadm_fun)) .msg("biofetchR: load_gadm() not available; skipping GADM.") else {
        overlay_gadm <- tryCatch(load_gadm_fun(iso2c=countries, level=gadm_unit, cache_dir=cache_dir_gadm, quiet=TRUE), error=function(e) NULL)
        if (is.list(overlay_gadm) && length(overlay_gadm) && all(vapply(overlay_gadm, function(o) inherits(o, "sf"), logical(1)))) {
          overlay_gadm <- do.call(rbind, lapply(overlay_gadm, function(o) .ensure_sf_wgs84(o)))
        }
        overlay_gadm <- tryCatch(.ensure_sf_wgs84(overlay_gadm), error=function(e) NULL)
        if (inherits(overlay_gadm, "sf") && nrow(overlay_gadm)) {
          nms <- names(overlay_gadm); ln <- tolower(nms)
          name_col <- NULL
          for (cand in c(sprintf("name_%d", gadm_unit), sprintf("nl_name_%d", gadm_unit), "geoname_gadm", "gadm_name", "name", "GADM_name")) {
            hit <- which(ln == tolower(cand)); if (length(hit)) { name_col <- nms[hit[1]]; break }
          }
          id_col <- NULL
          for (cand in c(sprintf("gid_%d", gadm_unit), sprintf("id_%d", gadm_unit), "gid", "GID")) {
            hit <- which(ln == tolower(cand)); if (length(hit)) { id_col <- nms[hit[1]]; break }
          }
          if (!is.null(name_col)) {
            overlay_gadm$name_gadm <- as.character(overlay_gadm[[name_col]])
            overlay_gadm$id_gadm <- if (!is.null(id_col)) as.character(overlay_gadm[[id_col]]) else sprintf("GADM_L%d_%s", gadm_unit, seq_len(nrow(overlay_gadm)))
            overlay_gadm$geoname_gadm <- overlay_gadm$name_gadm
            gcol <- .geom_name(overlay_gadm)
            overlay_gadm <- overlay_gadm[, c("id_gadm", "name_gadm", "geoname_gadm", gcol), drop=FALSE]
            overlay_gadm <- .make_valid(overlay_gadm)
          } else overlay_gadm <- NULL
        } else overlay_gadm <- NULL
      }
    }


    if ("teow" %in% sources) {
      overlay_teow <- .load_overlay_or_record(
        src = "teow",
        fun = load_teow_fun,
        args = list(
          cache_dir = teow_cache_dir,
          method = teow_method,
          source_url = teow_url,
          force_refresh = teow_force_refresh,
          quiet = TRUE
        ),
        crop_after = isTRUE(teow_clip_to_countries),
        crop_tag = "TEOW"
      )

      overlay_teow <- .prepare_overlay_or_record(
        "teow",
        overlay_teow,
        function(x) {
          if (!"geoname" %in% names(x)) {
            nm <- intersect(
              c("ECO_NAME", "eco_name", "ECOREGION", "ecoregion"),
              names(x)
            )
            x$geoname <- if (length(nm)) {
              as.character(x[[nm[[1L]]]])
            } else {
              paste0("ECO_", seq_len(nrow(x)))
            }
          }

          names(x)[names(x) == "geoname"] <- "geoname_teow"
          gcol <- .geom_name(x)
          x[, c("geoname_teow", gcol), drop = FALSE]
        }
      )
    }

    if ("feow" %in% sources) {
      overlay_feow <- .load_overlay_or_record(
        src = "feow",
        fun = load_feow_fun,
        args = list(
          cache_dir = feow_cache_dir,
          method = feow_method_resolved,
          source_url = feow_url,
          force_refresh = feow_force_refresh,
          quiet = TRUE
        ),
        crop_after = isTRUE(feow_clip_to_countries),
        crop_tag = "FEOW"
      )

      overlay_feow <- .prepare_overlay_or_record(
        "feow",
        overlay_feow,
        function(x) {
          if (!"geoname" %in% names(x)) {
            nm <- intersect(
              c("ecoregion", "ECOREGION", "ECO_NAME", "FEOW_NAME", "NAME"),
              names(x)
            )
            x$geoname <- if (length(nm)) {
              as.character(x[[nm[[1L]]]])
            } else {
              paste0("FEOW_", seq_len(nrow(x)))
            }
          }

          names(x)[names(x) == "geoname"] <- "geoname_feow"
          gcol <- .geom_name(x)
          x[, c("geoname_feow", gcol), drop = FALSE]
        }
      )
    }

    if ("lakes" %in% sources) {
      overlay_lakes <- .load_overlay_or_record(
        src = "lakes",
        fun = load_lakes_fun,
        args = list(
          cache_dir = lakes_cache_dir,
          force_refresh = lakes_force_refresh,
          lakes_only = lakes_only,
          min_area_km2 = lakes_min_area_km2,
          quiet = TRUE
        ),
        crop_after = TRUE,
        crop_tag = "LAKES"
      )

      overlay_lakes <- .prepare_overlay_or_record(
        "lakes",
        overlay_lakes,
        function(x) {
          nm <- intersect(c("name", "Lake_name", "lake_name"), names(x))
          x$geoname_lakes <- if (length(nm)) {
            as.character(x[[nm[[1L]]]])
          } else {
            paste0("LAKE_", seq_len(nrow(x)))
          }

          gcol <- .geom_name(x)
          x[, c("geoname_lakes", gcol), drop = FALSE]
        }
      )
    }

    if ("rivers" %in% sources) {
      overlay_rivers <- .load_overlay_or_record(
        src = "rivers",
        fun = load_rivers_fun,
        args = list(
          cache_dir = rivers_cache_dir,
          cache = rivers_cache,
          cache_format = rivers_cache_format,
          force_refresh = rivers_force_refresh,
          regions = rivers_regions,
          min_strahler = rivers_min_strahler,
          min_discharge_cms = rivers_min_discharge_cms,
          quiet = TRUE
        ),
        crop_after = TRUE,
        crop_tag = "RIVERS"
      )

      overlay_rivers <- .prepare_overlay_or_record(
        "rivers",
        overlay_rivers,
        function(x) {
          if (!"hyriv_id" %in% names(x) && "HYRIV_ID" %in% names(x)) {
            x$hyriv_id <- x$HYRIV_ID
          }
          if (!"hyriv_id" %in% names(x)) {
            stop("HydroRIVERS overlay lacks `hyriv_id`/`HYRIV_ID`.", call. = FALSE)
          }

          gcol <- .geom_name(x)
          x[, c("hyriv_id", gcol), drop = FALSE]
        }
      )
    }

    if ("basins" %in% sources) {
      overlay_basins <- .load_overlay_or_record(
        src = "basins",
        fun = load_basins_fun,
        args = list(
          level = basins_level,
          with_lakes = basins_with_lakes,
          cache_dir = basins_cache_dir,
          quiet = TRUE
        ),
        crop_after = isTRUE(basins_clip_to_countries),
        crop_tag = "BASINS"
      )

      overlay_basins <- .prepare_overlay_or_record(
        "basins",
        overlay_basins,
        function(x) {
          if (!"hybas_id" %in% names(x) && "HYBAS_ID" %in% names(x)) {
            x$hybas_id <- x$HYBAS_ID
          }
          if (!"hybas_id" %in% names(x)) {
            stop("HydroBASINS overlay lacks `hybas_id`/`HYBAS_ID`.", call. = FALSE)
          }

          gcol <- .geom_name(x)
          x[, c("hybas_id", gcol), drop = FALSE]
        }
      )
    }

    if ("gmba" %in% sources) {
      overlay_gmba <- .load_overlay_or_record(
        "gmba",
        function(cache_dir, layer, extent, force_refresh, quiet = TRUE) {
          .load_gmba_inventory(
            cache_dir = cache_dir,
            layer = layer,
            extent = extent,
            force_refresh = force_refresh,
            quiet_dl = quiet
          )
        },
        list(
          cache_dir = gmba_cache_dir,
          layer = gmba_layer,
          extent = gmba_extent,
          force_refresh = gmba_force_refresh,
          quiet = TRUE
        ),
        crop_after = gmba_clip_to_countries,
        crop_tag = "GMBA"
      )
    }


    if ("ne_urban" %in% sources) {
      overlay_ne_urban <- .load_overlay_or_record(
        src = "ne_urban",
        fun = load_ne_urban_fun,
        args = list(
          cache_dir = ne_cache_dir,
          scale = "10m",
          force_refresh = ne_force_refresh,
          clip_to_countries = TRUE,
          iso2c = countries,
          quiet = TRUE
        ),
        crop_after = TRUE,
        crop_tag = "NE_URBAN"
      )
    }

    if ("ne_admin1" %in% sources) {
      overlay_ne_admin1 <- .load_overlay_or_record(
        src = "ne_admin1",
        fun = load_ne_admin1_fun,
        args = list(
          cache_dir = ne_cache_dir,
          scale = ne_admin1_scale,
          force_refresh = ne_force_refresh,
          clip_to_countries = TRUE,
          iso2c = countries,
          quiet = TRUE
        )
      )
    }

    if ("resolve2017" %in% sources) {
      overlay_resolve2017 <- .load_overlay_or_record(
        src = "resolve2017",
        fun = load_resolve2017_fun,
        args = list(
          cache_dir = resolve_cache_dir,
          iso2c = countries,
          force_refresh = resolve_force_refresh,
          clip_to_countries = resolve_clip_to_countries,
          quiet = TRUE
        ),
        crop_after = isTRUE(resolve_clip_to_countries),
        crop_tag = "RESOLVE2017"
      )
    }

    if ("wdpa" %in% sources) {
      overlay_wdpa <- .load_overlay_or_record(
        src = "wdpa",
        fun = load_wdpa_fun,
        args = list(
          iso2c = countries,
          cache_dir = wdpa_cache_dir,
          force_refresh = wdpa_force_refresh,
          quiet = TRUE,
          exclude_marine = wdpa_exclude_marine,
          require_opt_in = wdpa_require_opt_in,
          opt_in = wdpa_opt_in
        ),
        crop_after = TRUE,
        crop_tag = "WDPA"
      )
    }

    if ("ramsar" %in% sources) {
      overlay_ramsar <- .load_overlay_or_record(
        src = "ramsar",
        fun = load_ramsar_fun,
        args = list(
          iso2c = countries,
          cache_dir = ramsar_cache_dir,
          force_refresh = ramsar_force_refresh,
          quiet = TRUE,
          require_opt_in = ramsar_require_opt_in,
          opt_in = ramsar_opt_in
        ),
        crop_after = TRUE,
        crop_tag = "RAMSAR"
      )
    }

    if ("gdw_barriers" %in% sources) {
      overlay_gdw_barriers <- .load_overlay_or_record(
        "gdw_barriers",
        load_gdw_barriers_fun,
        list(
          cache_dir = gdw_cache_dir,
          force_refresh = gdw_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = gdw_barriers_source_path,
          source_url = gdw_barriers_source_url,
          quiet = TRUE
        )
      )
    }


    if ("biosphere_reserve" %in% sources) {
      overlay_biosphere <- .load_overlay_or_record(
        "biosphere_reserve",
        load_biosphere_fun,
        list(
          cache_dir = biosphere_cache_dir,
          force_refresh = biosphere_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = biosphere_source_path,
          source_url = biosphere_source_url,
          quiet = TRUE
        )
      )
    }

    if ("gloric" %in% sources) {
      overlay_gloric <- .load_overlay_or_record(
        "gloric",
        load_gloric_fun,
        list(
          cache_dir = gloric_cache_dir,
          force_refresh = gloric_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = gloric_source_path,
          source_url = gloric_source_url,
          quiet = TRUE
        )
      )
    }

    if ("hydrowaste" %in% sources) {
      overlay_hydrowaste <- .load_overlay_or_record(
        "hydrowaste",
        load_hydrowaste_fun,
        list(
          cache_dir = hydrowaste_cache_dir,
          force_refresh = hydrowaste_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = hydrowaste_source_path,
          source_url = hydrowaste_source_url,
          quiet = TRUE
        )
      )
    }

    if ("gdw_reservoirs" %in% sources) {
      overlay_gdw_reservoirs <- .load_overlay_or_record(
        "gdw_reservoirs",
        load_gdw_reservoirs_fun,
        list(
          cache_dir = gdw_reservoirs_cache_dir,
          force_refresh = gdw_reservoirs_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = gdw_reservoirs_source_path,
          source_url = gdw_reservoirs_source_url,
          quiet = TRUE
        )
      )
    }

    if ("global_mining" %in% sources) {
      overlay_global_mining <- .load_overlay_or_record(
        "global_mining",
        load_global_mining_fun,
        list(
          cache_dir = global_mining_cache_dir,
          force_refresh = global_mining_force_refresh,
          iso2c = countries,
          clip_to_countries = TRUE,
          source_path = global_mining_source_path,
          source_url = global_mining_source_url,
          quiet = TRUE
        )
      )
    }
  }

  # ---------------------------------------------------------------------------
  # Pre-download overlay availability gate
  # ---------------------------------------------------------------------------
  # If a user explicitly requests GADM overlay assignment, the pipeline cannot
  # produce valid GADM-grouped exports unless GADM polygons have loaded
  # successfully. This check happens before any GBIF occurrence downloads are
  # submitted, so first-time users are not left with wasted GBIF downloads and
  # an empty gbif_summary when the GADM/geodata provider is unavailable.
  #
  # In strict mode, stop immediately. In fail-soft mode, write one clear
  # pre-download failure row per retained species-country input and return an
  # empty result. This preserves the existing package convention that failures
  # are auditable in gbif_summary.csv, while still preventing unnecessary GBIF
  # submissions when the requested spatial output cannot be produced.
  if (
    isTRUE(use_overlays) &&
    "gadm" %in% sources &&
    (is.null(overlay_gadm) || !inherits(overlay_gadm, "sf") || !nrow(overlay_gadm))
  ) {
    reason <- paste0(
      "Requested region_source = 'gadm', but no GADM polygons were loaded. ",
      "This usually means load_gadm() failed, geodata::gadm() could not retrieve ",
      "the requested boundary data, or the GADM/geodata provider is temporarily ",
      "unavailable. GBIF downloads were not submitted because GADM-grouped ",
      "exports cannot be produced without the requested overlay."
    )

    fail_rows <- df |>
      dplyr::distinct(.data$species, .data$iso2c)

    for (i in seq_len(nrow(fail_rows))) {
      .append_row(
        species = as.character(fail_rows$species[[i]]),
        region_id = as.character(fail_rows$iso2c[[i]]),
        region_type = .region_type_label("gadm"),
        n_total = 0L,
        n_cleaned = 0L,
        n_thinned = 0L,
        output_file = NA_character_,
        status = "failed",
        fail_stage = "load_overlay",
        fail_reason = reason,
        gbif_key = NA_character_
      )
    }

    .write_summary_now()

    if (isTRUE(strict_overlay_loading)) {
      stop(reason, call. = FALSE)
    }

    warning(reason, call. = FALSE)

    bf_console_pipeline_end(
      workflow = "terrestrial",
      n_species = length(species_vec_console),
      n_records_out = 0L,
      start_time = pipeline_start_time,
      quiet = quiet
    )

    if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
      return(tibble::tibble(
        decimalLongitude = numeric(),
        decimalLatitude = numeric()
      ))
    }

    return(invisible(NULL))
  }

  .msg("overlay_gadm rows: ", if (!is.null(overlay_gadm) && "gadm" %in% sources) nrow(overlay_gadm) else 0)

  if (isTRUE(use_overlays) && length(sources)) {
    bf_console_bullet(
      "Overlay resources loaded; starting GBIF processing.",
      quiet = quiet
    )
  }

  .export_group <- function(sf_group, species_name, region_id, region_type, workflow_tag, region_name = NA_character_, gbif_key = NA_character_) {
    n_total <- nrow(sf_group)
    if (!n_total) return(invisible(NULL))
    sf_group <- .ensure_sf_wgs84(sf_group)
    sf_group <- .ensure_lonlat_cols(sf_group)
    if (!".bf_rowid" %in% names(sf_group)) sf_group$.bf_rowid <- seq_len(nrow(sf_group))
    proc <- .clean_and_thin(sf_group)

    bf_console_bullet(
      paste0(
        "Preparing export group: ",
        species_name,
        " | ",
        region_type,
        " = ",
        region_id
      ),
      quiet = quiet
    )

    bf_console_cleaning(
      n_before = proc$n_total,
      n_after = if (isTRUE(apply_cleaning)) proc$n_cleaned else proc$n_total,
      quiet = quiet
    )

    bf_console_thinning(
      n_before = if (isTRUE(apply_cleaning) && !is.na(proc$n_cleaned)) proc$n_cleaned else proc$n_total,
      n_after = proc$n_thinned,
      applied = apply_thinning,
      quiet = quiet
    )

    proc$sf <- .ensure_lonlat_cols(proc$sf)
    proc$sf$region_id <- as.character(region_id)
    proc$sf$region_name <- as.character(region_name)
    proc$sf$region_type <- as.character(region_type)
    proc$sf$workflow <- as.character(workflow_tag)
    if (identical(workflow_tag, "gadm")) proc$sf$GADM_level <- gadm_unit

    # ---------------------------------------------------------------------------
    # Resolve and write grouped export through the universal export contract
    # ---------------------------------------------------------------------------
    # This replaces the old local `.expected_output_file()` helper. The output CSV,
    # gbif_summary$output_file and export_manifest.csv now all use the same resolved
    # export plan.

    input_region_value <- if ("iso2c" %in% names(proc$sf)) {
      vals <- unique(as.character(proc$sf$iso2c))
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (length(vals)) paste(vals, collapse = "_") else "GLOBAL"
    } else {
      "GLOBAL"
    }

    export_stage_value <- if (identical(as.character(workflow_tag), "gadm")) {
      paste0("gadm_l", gadm_unit)
    } else {
      as.character(workflow_tag)
    }

    export_plan <- bf_resolve_export_plan(
      output_dir = output_dir,
      species = species_name,
      input_region = input_region_value,
      export_stage = export_stage_value,
      workflow = as.character(workflow_tag),
      region_type = as.character(region_type),
      region_source = as.character(workflow_tag),
      spatial_join_type = if (identical(as.character(workflow_tag), "gadm")) "gadm" else NA_character_,
      gadm_unit = if (identical(as.character(workflow_tag), "gadm")) gadm_unit else NULL,
      region_id = as.character(region_id),
      split_by_region = TRUE,
      apply_cleaning = apply_cleaning,
      apply_thinning = apply_thinning,
      thinning_dist_km = dist_km,
      use_overlays = isTRUE(use_overlays),
      overlay_names = as.character(workflow_tag),
      prepare_taxonomy = prepare_taxonomy,
      use_native_range = use_native_range,
      use_native_web = use_native_web
    )

    out_csv <- sf::st_drop_geometry(proc$sf)
    out_csv <- .restore_coords_df(out_csv, proc$sf)
    out_csv <- .strip_internal_cols(out_csv)

    export_plan <- tryCatch(
      {
        bf_write_export_csv(out_csv, export_plan)
      },
      error = function(e) {
        .msg(
          "Universal export writer failed for ",
          species_name,
          " / ",
          region_id,
          "; falling back to readr::write_csv(). Reason: ",
          .errmsg(e)
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
          region_name = as.character(region_name),
          gbif_key = as.character(gbif_key)
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
      out_csv
    } else {
      output_file_path
    }

    .append_row(
      species = species_name,
      region_id = as.character(region_id),
      region_type = as.character(region_type),
      n_total = n_total,
      n_cleaned = if (isTRUE(apply_cleaning)) proc$n_cleaned else NA_integer_,
      n_thinned = proc$n_thinned,
      output_file = output_file_path,
      status = "success",
      gbif_key = gbif_key
    )

    if (exists("species_output_rows_console", inherits = TRUE)) {
      species_output_rows_console <<- species_output_rows_console + as.integer(proc$n_thinned)
    }

    if (isTRUE(return_all_results) && isTRUE(store_in_memory) && is.data.frame(out_data)) {
      all_results[[paste(species_name, region_id, workflow_tag, sep = "_")]] <<- tibble::as_tibble(out_data)
    }
    invisible(NULL)
  }

  for (batch in batches) {
    for (species_name in batch) {

      species_start_time <- bf_console_now()
      species_output_rows_console <- 0L

      species_index_console <- match(species_name, species_vec_console)

      bf_console_species_start(
        species = species_name,
        index = species_index_console,
        total = length(species_vec_console),
        region = if (length(sources)) paste(sources, collapse = "+") else "raw",
        quiet = quiet
      )

      iso2_codes <- df |>
        dplyr::filter(.data$species == !!species_name) |>
        dplyr::pull(.data$iso2c) |>
        unique()

      ok_vec <- vapply(iso2_codes, function(cty) .gbif_has_geo_records(species_name, cty), logical(1))
      valid_countries <- iso2_codes[ok_vec]
      skipped_countries <- iso2_codes[!ok_vec]

      if (length(skipped_countries)) {
        for (cc in skipped_countries) {
          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "skipped",
            fail_stage = "precheck_no_geo",
            fail_reason = "GBIF precheck returned 0 georeferenced records (species x country)."
          )
          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            retrieval_status = "precheck_zero",
            n_records = 0L,
            fail_reason = "GBIF precheck returned 0 georeferenced records (species x country)."
          )
        }
      }
      if (!length(valid_countries)) next

      if (is.null(download_fun)) stop("No download function available (download_gbif_batch_gadm missing and deps$download_fun not provided).")

      download_start_time <- bf_console_now()

      bf_console_bullet(
        paste0(
          "Submitting GBIF request for ",
          species_name,
          " across ",
          length(valid_countries),
          " country/countries."
        ),
        quiet = quiet
      )

      download_keys <- tryCatch(
        download_fun(
          species = species_name,
          iso2_codes = valid_countries,
          user = user,
          pwd = pwd,
          email = email
        ),
        error = function(e) e
      )

      if (inherits(download_keys, "condition")) {
        for (cc in valid_countries) {
          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "failed",
            fail_stage = "submit_download",
            fail_reason = .errmsg(download_keys)
          )
          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            retrieval_status = "submit_failed",
            fail_reason = .errmsg(download_keys)
          )
        }

        bf_console_bullet(
          paste0("GBIF request failed for ", species_name, "."),
          quiet = quiet
        )
        next
      }

      submission_audit <- attr(download_keys, "submission_audit", exact = TRUE)

      if (is.null(download_keys) || !length(download_keys)) {
        for (cc in valid_countries) {
          row_cc <- NULL
          if (is.data.frame(submission_audit) && "label" %in% names(submission_audit)) {
            hit <- which(as.character(submission_audit$label) == cc)
            if (length(hit)) row_cc <- submission_audit[hit[[1L]], , drop = FALSE]
          }

          sub_status <- if (!is.null(row_cc) && "submission_status" %in% names(row_cc)) {
            as.character(row_cc$submission_status[[1L]])
          } else {
            "submit_failed"
          }
          sub_reason <- if (!is.null(row_cc) && "message" %in% names(row_cc)) {
            as.character(row_cc$message[[1L]])
          } else {
            "GBIF request returned no usable download key."
          }
          sub_key <- if (!is.null(row_cc) && "gbif_key" %in% names(row_cc)) {
            as.character(row_cc$gbif_key[[1L]])
          } else {
            NA_character_
          }

          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "failed",
            fail_stage = sub_status,
            fail_reason = sub_reason,
            gbif_key = sub_key
          )
          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            gbif_key = sub_key,
            retrieval_status = sub_status,
            fail_reason = sub_reason
          )
        }

        bf_console_bullet(
          paste0("GBIF request failed or returned no key for ", species_name, "."),
          quiet = quiet
        )
        next
      }

      bf_console_gbif_request(
        species = species_name,
        elapsed = bf_console_elapsed(download_start_time),
        gbif_key = unique(as.character(unlist(download_keys))),
        quiet = quiet
      )

      keys_named <- NULL
      if (is.character(download_keys)) keys_named <- download_keys
      if (is.list(download_keys)) {
        if (!is.null(names(download_keys))) {
          keys_named <- vapply(
            download_keys,
            function(z) if (length(z)) as.character(z[[1]]) else NA_character_,
            character(1)
          )
          names(keys_named) <- names(download_keys)
        } else {
          keys_named <- as.character(unlist(download_keys))
        }
      }
      if (is.null(keys_named)) keys_named <- as.character(unlist(download_keys))

      import_start_time <- bf_console_now()

      results_list <- tryCatch(
        {
          if (is.function(import_fun)) {
            import_fun(download_keys)
          } else if (exists("wait_and_import_gbif_safe", mode = "function")) {
            wait_and_import_gbif_safe(download_keys)
          } else if (exists("wait_and_import_gbif", mode = "function")) {
            wait_and_import_gbif(download_keys)
          } else {
            stop("No import function available (wait_and_import_gbif[_safe] missing and deps$import_fun not provided).")
          }
        },
        error = function(e) e
      )

      if (inherits(results_list, "condition")) {
        for (cc in valid_countries) {
          key_cc <- if (!is.null(names(keys_named)) && cc %in% names(keys_named)) {
            keys_named[[cc]]
          } else {
            NA_character_
          }

          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "failed",
            fail_stage = "wait_import",
            fail_reason = .errmsg(results_list),
            gbif_key = key_cc
          )
          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            gbif_key = key_cc,
            retrieval_status = "import_failed",
            fail_reason = .errmsg(results_list)
          )
        }
        next
      }

      import_audit <- attr(results_list, "gbif_status", exact = TRUE)
      result_names <- names(results_list)
      if (is.null(result_names)) result_names <- character(0)

      # Reconcile every requested country before downstream processing. A country
      # may be absent from `results_list` only when the import audit explicitly
      # explains why.
      for (cc in valid_countries) {
        row_cc <- NULL
        if (is.data.frame(import_audit) && "label" %in% names(import_audit)) {
          hit <- which(as.character(import_audit$label) == cc)
          if (length(hit)) row_cc <- import_audit[hit[[1L]], , drop = FALSE]
        }

        if (!is.null(row_cc)) {
          retrieval_status_cc <- as.character(row_cc$retrieval_status[[1L]])
          gbif_status_cc <- as.character(row_cc$gbif_status[[1L]])
          n_records_cc <- suppressWarnings(as.integer(row_cc$n_records[[1L]]))
          fail_reason_cc <- as.character(row_cc$message[[1L]])
          key_cc <- as.character(row_cc$gbif_key[[1L]])

          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            gbif_key = key_cc,
            gbif_status = gbif_status_cc,
            retrieval_status = retrieval_status_cc,
            n_records = n_records_cc,
            fail_reason = fail_reason_cc
          )

          if (!retrieval_status_cc %in% c("success", "success_zero")) {
            .append_row(
              species = species_name,
              region_id = cc,
              region_type = .region_type_for_fail(),
              status = "failed",
              fail_stage = retrieval_status_cc,
              fail_reason = fail_reason_cc,
              gbif_key = key_cc
            )
          }
        } else if (cc %in% result_names) {
          n_records_cc <- as.integer(bf_console_nrow(results_list[[cc]]))
          key_cc <- if (!is.null(names(keys_named)) && cc %in% names(keys_named)) {
            keys_named[[cc]]
          } else {
            NA_character_
          }

          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            gbif_key = key_cc,
            gbif_status = "SUCCEEDED",
            retrieval_status = if (n_records_cc > 0L) "success" else "success_zero",
            n_records = n_records_cc
          )
        } else {
          key_cc <- if (!is.null(names(keys_named)) && cc %in% names(keys_named)) {
            keys_named[[cc]]
          } else {
            NA_character_
          }
          reason_cc <- "Country request was neither returned by the importer nor represented in its retrieval audit."

          .append_retrieval_audit(
            species = species_name,
            request_id = cc,
            gbif_key = key_cc,
            retrieval_status = "internal_unreconciled",
            fail_reason = reason_cc
          )
          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "failed",
            fail_stage = "internal_unreconciled",
            fail_reason = reason_cc,
            gbif_key = key_cc
          )
        }
      }

      if (is.null(results_list) || !length(results_list)) {
        bf_console_bullet(
          paste0("GBIF import produced no successfully resolved country objects for ", species_name, "."),
          quiet = quiet
        )
        next
      }

      imported_record_counts_console <- vapply(
        results_list,
        function(x) as.integer(bf_console_nrow(x)),
        integer(1)
      )

      bf_console_gbif_import(
        n_records = sum(imported_record_counts_console, na.rm = TRUE),
        elapsed = bf_console_elapsed(import_start_time),
        quiet = quiet
      )

      if (is.null(names(results_list))) names(results_list) <- valid_countries[seq_len(min(length(valid_countries), length(results_list)))]

      for (cc in names(results_list)) {
        gbif_key_cc <- if (!is.null(names(keys_named)) && cc %in% names(keys_named)) keys_named[[cc]] else NA_character_
        raw_sf <- results_list[[cc]]

        if (is.null(raw_sf) || !nrow(raw_sf)) {
          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "empty",
            fail_stage = "import_empty",
            fail_reason = "GBIF download/import returned 0 rows.",
            gbif_key = gbif_key_cc
          )
          next
        }

        raw_sf <- tryCatch(.ensure_sf_wgs84(raw_sf), error=function(e) NULL)
        if (is.null(raw_sf) || !nrow(raw_sf)) {
          .append_row(
            species = species_name,
            region_id = cc,
            region_type = .region_type_for_fail(),
            status = "failed",
            fail_stage = "ensure_sf_wgs84",
            fail_reason = "Failed to coerce/transform to EPSG:4326.",
            gbif_key = gbif_key_cc
          )
          next
        }

        raw_sf <- .ensure_lonlat_cols(raw_sf)
        if (!"species" %in% names(raw_sf)) raw_sf$species <- species_name else raw_sf$species <- species_name

        # Attach species-country metadata produced by taxonomy / origin gates to
        # every occurrence before overlay splitting. This lets exported GBIF rows
        # retain GRIIS/native-range evidence without changing the download queue.
        if (exists("species_metadata", inherits = TRUE) && length(species_metadata_cols)) {
          meta_row <- species_metadata |>
            dplyr::filter(.data$species == !!species_name, .data$iso2c == !!cc) |>
            dplyr::slice_head(n = 1)

          if (nrow(meta_row)) {
            for (meta_col in species_metadata_cols) {
              if (!meta_col %in% names(raw_sf)) {
                raw_sf[[meta_col]] <- meta_row[[meta_col]][[1]]
              }
            }
          }
        }

        bf_console_bullet(
          paste0(
            species_name,
            " / ",
            cc,
            ": ",
            nrow(raw_sf),
            " georeferenced GBIF records imported."
          ),
          quiet = quiet
        )

        joined <- raw_sf

        if (length(raster_context) && is.function(enrich_raster_fun)) {
          joined <- tryCatch(
            enrich_raster_fun(
              sf_points = joined,
              raster_context = raster_context,
              cache_dir = raster_cache_dir,
              worldcover_source = worldcover_source,
              worldcover_year = worldcover_year,
              worldcover_buffer_m = worldcover_buffer_m,
              soilgrids_sources = soilgrids_sources,
              soilgrids_vars = soilgrids_vars,
              soilgrids_buffer_m = soilgrids_buffer_m,
              human_footprint_source = human_footprint_source,
              human_footprint_buffer_m = human_footprint_buffer_m,
              quiet = quiet
            ),
            error = function(e) {
              reason <- .errmsg(e)

              .append_row(
                species = species_name,
                region_id = cc,
                region_type = "RASTER_CONTEXT",
                status = "failed",
                fail_stage = "raster_context",
                fail_reason = reason,
                gbif_key = gbif_key_cc
              )

              if (isTRUE(strict_raster_context)) {
                .write_summary_now()
                stop("Raster-context enrichment failed for ", species_name, " / ", cc, ": ", reason, call. = FALSE)
              }

              .msg("Raster enrichment failed for ", species_name, " / ", cc, ": ", reason)
              joined
            }
          )
        }

        if (!isTRUE(use_overlays) || !length(sources)) {
          .export_group(joined, species_name, cc, "RAW_COUNTRY", "raw", region_name = cc, gbif_key = gbif_key_cc)
          next
        }

        overlay_join_start_time <- bf_console_now()

        bf_console_overlay_start(
          region_source = paste(sources, collapse = "+"),
          gadm_unit = if ("gadm" %in% sources) gadm_unit else NULL,
          quiet = quiet
        )

        if ("gadm" %in% sources) {
          joined$id_gadm <- joined$name_gadm <- joined$geoname_gadm <- NA_character_
          if (!is.null(overlay_gadm)) {
            g <- .join_gadm(joined, overlay_gadm)
            joined$id_gadm <- g$id
            joined$name_gadm <- g$name
            joined$geoname_gadm <- g$geoname
          }
          .msg("GADM labeled (", cc, "): ", sum(!is.na(joined$id_gadm)), " / ", nrow(joined))
        }

        if ("teow" %in% sources && !is.null(overlay_teow)) joined$geoname_teow <- .join_poly_label(joined, overlay_teow, "geoname_teow")
        if ("feow" %in% sources && !is.null(overlay_feow)) joined$geoname_feow <- .join_poly_label(joined, overlay_feow, "geoname_feow")
        if ("lakes" %in% sources && !is.null(overlay_lakes)) joined$geoname_lakes <- .join_poly_label(joined, overlay_lakes, "geoname_lakes")
        if ("basins" %in% sources && !is.null(overlay_basins)) joined$hybas_id <- .join_poly_label(joined, overlay_basins, "hybas_id")
        if ("rivers" %in% sources && !is.null(overlay_rivers)) joined$hyriv_id <- .join_nearest_id(joined, overlay_rivers, "hyriv_id", snap_m = rivers_max_snap_m)
        if ("gmba" %in% sources && !is.null(overlay_gmba)) joined$gmba_id <- .join_poly_label(joined, overlay_gmba, "gmba_id")
        if ("ne_urban" %in% sources && !is.null(overlay_ne_urban)) joined$geoname_ne_urban <- .join_poly_label(joined, overlay_ne_urban, "geoname_ne_urban")
        if ("resolve2017" %in% sources && !is.null(overlay_resolve2017)) joined$geoname_resolve2017 <- .join_poly_label(joined, overlay_resolve2017, "geoname_resolve2017")
        if ("ne_admin1" %in% sources && !is.null(overlay_ne_admin1)) {
          joined$id_ne_admin1 <- .join_poly_label(joined, overlay_ne_admin1, "id_ne_admin1")
          joined$geoname_ne_admin1 <- .join_poly_label(joined, overlay_ne_admin1, "geoname_ne_admin1")
        }
        if ("wdpa" %in% sources && !is.null(overlay_wdpa)) {
          joined$wdpa_id <- .join_poly_label(joined, overlay_wdpa, "wdpa_id")
          joined$wdpa_name <- .join_poly_label(joined, overlay_wdpa, "wdpa_name")
        }
        if ("ramsar" %in% sources && !is.null(overlay_ramsar)) {
          joined$ramsar_id <- .join_poly_label(joined, overlay_ramsar, "ramsar_id")
          joined$ramsar_name <- .join_poly_label(joined, overlay_ramsar, "ramsar_name")
        }
        if ("gdw_barriers" %in% sources && !is.null(overlay_gdw_barriers)) {
          joined$gdw_barrier_id <- .join_nearest_id(joined, overlay_gdw_barriers, "gdw_barrier_id", snap_m = gdw_barrier_max_snap_m)
          joined$gdw_barrier_name <- .join_nearest_attr(joined, overlay_gdw_barriers, "gdw_barrier_name", snap_m = gdw_barrier_max_snap_m)
        }
        if ("biosphere_reserve" %in% sources && !is.null(overlay_biosphere)) {
          geom_type_bio <- tryCatch(as.character(unique(sf::st_geometry_type(overlay_biosphere)))[1], error = function(e) "")
          if (grepl("POINT", geom_type_bio, ignore.case = TRUE)) {
            joined$biosphere_id <- .join_nearest_id(joined, overlay_biosphere, "biosphere_id", snap_m = biosphere_max_snap_m)
            joined$biosphere_name <- .join_nearest_attr(joined, overlay_biosphere, "biosphere_name", snap_m = biosphere_max_snap_m)
          } else {
            joined$biosphere_id <- .join_poly_label(joined, overlay_biosphere, "biosphere_id")
            joined$biosphere_name <- .join_poly_label(joined, overlay_biosphere, "biosphere_name")
          }
        }
        if ("gloric" %in% sources && !is.null(overlay_gloric)) {
          joined$gloric_id <- .join_nearest_id(joined, overlay_gloric, "gloric_id", snap_m = gloric_max_snap_m)
          joined$gloric_class <- .join_nearest_attr(joined, overlay_gloric, "gloric_class", snap_m = gloric_max_snap_m)
        }

        if ("hydrowaste" %in% sources && !is.null(overlay_hydrowaste)) {
          joined$hydrowaste_id <- .join_nearest_id(joined, overlay_hydrowaste, "hydrowaste_id", snap_m = hydrowaste_max_snap_m)
          joined$hydrowaste_name <- .join_nearest_attr(joined, overlay_hydrowaste, "hydrowaste_name", snap_m = hydrowaste_max_snap_m)
        }

        if ("gdw_reservoirs" %in% sources && !is.null(overlay_gdw_reservoirs)) {
          joined$gdw_reservoir_id <- .join_poly_label(joined, overlay_gdw_reservoirs, "gdw_reservoir_id")
          joined$gdw_reservoir_name <- .join_poly_label(joined, overlay_gdw_reservoirs, "gdw_reservoir_name")
        }

        if ("global_mining" %in% sources && !is.null(overlay_global_mining)) {
          joined$global_mining_id <- .join_poly_label(joined, overlay_global_mining, "global_mining_id")
          joined$global_mining_name <- .join_poly_label(joined, overlay_global_mining, "global_mining_name")
        }

        for (src_console in sources) {
          id_col_console <- label_cols[[src_console]]

          if (!is.null(id_col_console) && id_col_console %in% names(joined)) {
            bf_console_overlay_result(
              n_records = nrow(joined),
              n_matched = bf_console_count_matched(joined, id_col_console),
              id_col = id_col_console,
              quiet = quiet
            )
          }
        }

        bf_console_bullet(
          paste0(
            "Spatial overlay join completed in ",
            bf_console_elapsed(overlay_join_start_time),
            "."
          ),
          quiet = quiet
        )

        do_one_source <- function(src, dat) {
          col <- label_cols[[src]]
          if (is.null(col) || !col %in% names(dat)) return(invisible(NULL))
          x <- dat[!is.na(dat[[col]]), , drop=FALSE]
          if (!nrow(x)) return(invisible(NULL))
          rt <- .region_type_label(src)
          sp <- split(x, as.character(x[[col]]))
          for (rid in names(sp)) {
            rname <- rid
            if (src == "gadm" && "name_gadm" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$name_gadm)[1]
            if (src == "ne_admin1" && "geoname_ne_admin1" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$geoname_ne_admin1)[1]
            if (src == "wdpa" && "wdpa_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$wdpa_name)[1]
            if (src == "ramsar" && "ramsar_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$ramsar_name)[1]
            if (src == "gdw_barriers" && "gdw_barrier_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$gdw_barrier_name)[1]
            if (src == "biosphere_reserve" && "biosphere_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$biosphere_name)[1]
            if (src == "gloric" && "gloric_class" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$gloric_class)[1]
            if (src == "hydrowaste" && "hydrowaste_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$hydrowaste_name)[1]
            if (src == "gdw_reservoirs" && "gdw_reservoir_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$gdw_reservoir_name)[1]
            if (src == "global_mining" && "global_mining_name" %in% names(sp[[rid]])) rname <- unique(sp[[rid]]$global_mining_name)[1]
            .export_group(sp[[rid]], species_name, rid, rt, src, region_name = rname, gbif_key = gbif_key_cc)
          }
          invisible(NULL)
        }

        if (overlay_mode == "single" || length(sources) == 1) {
          do_one_source(sources[1], joined)
        } else if (overlay_mode == "dual_separate") {
          for (src in sources) do_one_source(src, joined)
        } else if (overlay_mode == "dual_intersection") {
          if (length(sources) != 2) {
            .msg("overlay_mode='dual_intersection' requires exactly 2 sources; falling back to dual_separate.")
            for (src in sources) do_one_source(src, joined)
          } else {
            s1 <- sources[1]; s2 <- sources[2]
            c1 <- label_cols[[s1]]; c2 <- label_cols[[s2]]
            if (!is.null(c1) && !is.null(c2) && c1 %in% names(joined) && c2 %in% names(joined)) {
              xi <- joined[!is.na(joined[[c1]]) & !is.na(joined[[c2]]), , drop=FALSE]
              if (nrow(xi)) {
                xi$combo <- paste0(xi[[c1]], " :: ", xi[[c2]])
                sp <- split(xi, xi$combo)
                for (rid in names(sp)) {
                  .export_group(sp[[rid]], species_name, rid, paste0(.region_type_label(s1), "x", .region_type_label(s2)), paste0(s1, "_", s2), region_name = rid, gbif_key = gbif_key_cc)
                }
              }
            }
          }
        }
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
    # Final completeness gate: every taxonomy-approved species x country request
    # must have one retrieval outcome, even when GBIF itself failed or timed out.
    expected_retrievals <- df |>
      dplyr::distinct(.data$species, .data$iso2c) |>
      dplyr::transmute(
        species = as.character(.data$species),
        request_id = as.character(.data$iso2c)
      )

    observed_retrievals <- retrieval_audit |>
      dplyr::distinct(.data$species, .data$request_id)

    missing_retrievals <- dplyr::anti_join(
      expected_retrievals,
      observed_retrievals,
      by = c("species", "request_id")
    )

    if (nrow(missing_retrievals)) {
      for (i in seq_len(nrow(missing_retrievals))) {
        .append_retrieval_audit(
          species = missing_retrievals$species[[i]],
          request_id = missing_retrievals$request_id[[i]],
          retrieval_status = "internal_unreconciled",
          fail_reason = "Final retrieval completeness gate found no recorded outcome for this request."
        )
      }
      warning(
        nrow(missing_retrievals),
        " species x country GBIF request(s) lacked a retrieval outcome and were marked internal_unreconciled.",
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

  if (isTRUE(return_all_results) && isTRUE(store_in_memory)) {
    if (!length(all_results)) {
      out <- tibble::tibble(
        decimalLongitude = numeric(),
        decimalLatitude = numeric()
      )

      bf_console_pipeline_end(
        workflow = "terrestrial",
        n_species = length(species_vec_console),
        n_records_out = 0L,
        start_time = pipeline_start_time,
        quiet = quiet
      )

      return(out)
    }

    out <- .bf_bind_rows_safe(all_results)
    out <- .strip_internal_cols(out)

    bf_console_pipeline_end(
      workflow = "terrestrial",
      n_species = length(species_vec_console),
      n_records_out = nrow(out),
      start_time = pipeline_start_time,
      quiet = quiet
    )

    return(out)
  }

  bf_console_pipeline_end(
    workflow = "terrestrial",
    n_species = length(species_vec_console),
    n_records_out = if ("n_thinned" %in% names(summary_tbl)) {
      sum(summary_tbl$n_thinned, na.rm = TRUE)
    } else {
      NA_integer_
    },
    start_time = pipeline_start_time,
    quiet = quiet
  )

  invisible(NULL)
}
