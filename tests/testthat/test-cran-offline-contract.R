test_that("CRAN offline contract is safe", {
  live_flag <- tolower(Sys.getenv("BIOFETCHR_RUN_LIVE_GBIF_DOWNLOADS", unset = "false"))
  expect_false(live_flag %in% c("1", "true", "yes", "y"))

  expect_true("bf_available_native_web_sources" %in% getNamespaceExports("biofetchR"))
  expect_true("bf_available_marine_overlays" %in% getNamespaceExports("biofetchR"))

  native_sources <- biofetchR::bf_available_native_web_sources()
  marine_overlays <- biofetchR::bf_available_marine_overlays()

  expect_true(is.data.frame(native_sources) || is.list(native_sources))
  expect_true(is.data.frame(marine_overlays) || is.list(marine_overlays))
})
