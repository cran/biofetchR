# tests/testthat/test-strict-overlay-loading.R

testthat::test_that(
  "strict overlay loading stops before GBIF for legacy overlays",
  {
    testthat::skip_if_not_installed("sf")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("readr")
    testthat::skip_if_not_installed("tibble")

    cases <- data.frame(
      source = c(
        "teow", "feow", "lakes", "rivers", "basins",
        "ne_urban", "ne_admin1", "resolve2017", "wdpa", "ramsar"
      ),
      dep_key = c(
        "load_teow", "load_feow", "load_lakes", "load_rivers", "load_basins",
        "load_ne_urban", "load_ne_admin1", "load_resolve2017",
        "load_wdpa", "load_ramsar"
      ),
      stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(cases))) {
      src <- cases$source[[i]]
      dep_key <- cases$dep_key[[i]]

      out_dir <- tempfile(paste0("bf_strict_", src, "_"))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      download_reached <- FALSE

      failing_loader <- local({
        source_name <- src
        function(...) {
          stop(
            paste0("synthetic ", source_name, " overlay failure"),
            call. = FALSE
          )
        }
      })

      mock_download <- function(...) {
        download_reached <<- TRUE
        stop("GBIF submission should not be reached in strict mode.", call. = FALSE)
      }

      deps <- list(download_fun = mock_download)
      deps[[dep_key]] <- failing_loader

      testthat::expect_error(
        process_gbif_terrestrial_freshwater_pipeline(
          df = tibble::tibble(
            species = "Xenopus laevis",
            iso2c = "FR"
          ),
          output_dir = out_dir,
          user = "mock_user",
          pwd = "mock_pwd",
          email = "mock@example.com",
          region_source = src,
          use_overlays = TRUE,
          strict_overlay_loading = TRUE,
          prepare_taxonomy = FALSE,
          apply_cleaning = FALSE,
          apply_thinning = FALSE,
          export_summary = TRUE,
          store_in_memory = FALSE,
          return_all_results = FALSE,
          quiet = TRUE,
          deps = deps
        ),
        regexp = paste0("Requested overlay failed to load: ", src)
      )

      testthat::expect_false(
        download_reached,
        info = paste("GBIF downloader was reached for", src)
      )

      summary_path <- file.path(out_dir, "gbif_summary.csv")
      testthat::expect_true(
        file.exists(summary_path),
        info = paste("gbif_summary.csv missing for", src)
      )

      if (file.exists(summary_path)) {
        summary_tbl <- readr::read_csv(
          summary_path,
          col_types = readr::cols(.default = readr::col_character()),
          show_col_types = FALSE,
          progress = FALSE
        )

        loader_rows <- summary_tbl[
          summary_tbl$species == paste0("__overlay_loader__", src),
          ,
          drop = FALSE
        ]

        testthat::expect_equal(
          nrow(loader_rows),
          1L,
          info = paste("expected one loader audit row for", src)
        )

        if (nrow(loader_rows) == 1L) {
          testthat::expect_equal(loader_rows$fail_stage[[1L]], "load_overlay")
          testthat::expect_match(
            loader_rows$fail_reason[[1L]],
            paste0("synthetic ", src, " overlay failure"),
            fixed = TRUE
          )
        }
      }
    }
  }
)

testthat::test_that(
  "fail-soft FEOW loading records the failure and can continue downstream",
  {
    testthat::skip_if_not_installed("sf")
    testthat::skip_if_not_installed("dplyr")
    testthat::skip_if_not_installed("readr")
    testthat::skip_if_not_installed("tibble")

    out_dir <- tempfile("bf_failsoft_feow_")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    download_reached <- FALSE

    failing_feow <- function(...) {
      stop("synthetic FEOW fail-soft failure", call. = FALSE)
    }

    # ISO2 "NA" bypasses the live GBIF pre-check in the terrestrial/freshwater
    # pipeline, so this regression test remains network-independent.
    mock_download <- function(...) {
      download_reached <<- TRUE
      stop("synthetic downstream submission stop", call. = FALSE)
    }

    testthat::expect_no_error(
      process_gbif_terrestrial_freshwater_pipeline(
        df = tibble::tibble(
          species = "Xenopus laevis",
          iso2c = "NA"
        ),
        output_dir = out_dir,
        user = "mock_user",
        pwd = "mock_pwd",
        email = "mock@example.com",
        region_source = "feow",
        use_overlays = TRUE,
        strict_overlay_loading = FALSE,
        prepare_taxonomy = FALSE,
        apply_cleaning = FALSE,
        apply_thinning = FALSE,
        export_summary = TRUE,
        store_in_memory = FALSE,
        return_all_results = FALSE,
        quiet = TRUE,
        deps = list(
          load_feow = failing_feow,
          download_fun = mock_download
        )
      )
    )

    testthat::expect_true(download_reached)

    summary_path <- file.path(out_dir, "gbif_summary.csv")
    testthat::expect_true(file.exists(summary_path))

    if (file.exists(summary_path)) {
      summary_tbl <- readr::read_csv(
        summary_path,
        col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      )

      loader_rows <- summary_tbl[
        summary_tbl$species == "__overlay_loader__feow",
        ,
        drop = FALSE
      ]

      testthat::expect_equal(nrow(loader_rows), 1L)

      if (nrow(loader_rows) == 1L) {
        testthat::expect_equal(loader_rows$fail_stage[[1L]], "load_overlay")
        testthat::expect_match(
          loader_rows$fail_reason[[1L]],
          "synthetic FEOW fail-soft failure",
          fixed = TRUE
        )
      }
    }
  }
)
