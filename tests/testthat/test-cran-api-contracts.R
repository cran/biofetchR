test_that("core public API functions are exported", {
  expected_exports <- c(
    "process_gbif_terrestrial_freshwater_pipeline",
    "process_gbif_marine_pipeline",
    "bf_attach_native_status",
    "bf_fetch_native_ranges_web",
    "load_gadm",
    "load_all_gadm",
    "gadm_join",
    "bf_load_ne_admin1",
    "bf_load_ne_urban",
    "bf_load_resolve2017",
    "bf_load_wdpa",
    "bf_load_ramsar",
    "bf_load_hydrowaste",
    "bf_load_gdw_barriers",
    "bf_load_global_mining",
    "bf_enrich_raster_context"
  )

  package_exports <- getNamespaceExports("biofetchR")
  missing_exports <- setdiff(expected_exports, package_exports)

  expect_true(
    length(missing_exports) == 0,
    label = paste(
      "Missing expected exports:",
      paste(missing_exports, collapse = ", ")
    )
  )
})

test_that("removed plotting API is not exported", {
  package_exports <- getNamespaceExports("biofetchR")

  expect_false(
    "list_plot_projections" %in% package_exports,
    label = "list_plot_projections should not be exported because plotting helpers were removed"
  )
})
