test_that("synthetic native-status workflow works without live downloads", {
  native_ranges <- data.frame(
    species = c("Rattus rattus", "Sturnus vulgaris"),
    native_origin_iso3 = c("IND;PAK", "DEU;FRA;GBR"),
    native_origin_flag = c("HAS_ORIGIN", "HAS_ORIGIN"),
    native_sources_used = c("synthetic", "synthetic"),
    stringsAsFactors = FALSE
  )

  recipient <- data.frame(
    species = c("Rattus rattus", "Sturnus vulgaris"),
    iso2c = c("GB", "DE"),
    stringsAsFactors = FALSE
  )

  attached <- biofetchR::bf_attach_native_status(
    recipient,
    species_col = "species",
    iso2c_col = "iso2c",
    native_ranges = native_ranges,
    native_filter_mode = "audit_only",
    quiet = TRUE
  )

  expect_equal(nrow(attached), nrow(recipient))
  expect_true(any(attached$species == "Rattus rattus" &
                  attached$native_status == "non_native_recipient"))
  expect_true(any(attached$species == "Sturnus vulgaris" &
                  attached$native_status == "native_recipient"))

  status_summary <- biofetchR::bf_native_status_summary(attached)
  expect_true(is.data.frame(status_summary))
  expect_gt(nrow(status_summary), 0)
})
