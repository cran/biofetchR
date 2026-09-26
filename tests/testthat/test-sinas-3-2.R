testthat::test_that("SInAS 3.2 is the default package-managed release", {
  urls <- bf_sinas_default_urls()

  testthat::expect_identical(
    unname(urls[["main_csv"]]),
    "https://zenodo.org/records/21933976/files/SInAS_3.2.csv?download=1"
  )
  testthat::expect_identical(
    unname(urls[["config_zip"]]),
    "https://zenodo.org/records/21933976/files/All_Config_Files_v3.2.zip?download=1"
  )
  testthat::expect_identical(
    unname(urls[["output_zip"]]),
    "https://zenodo.org/records/21933976/files/All_Output_Files_v3.2.zip?download=1"
  )
  testthat::expect_true(is.na(unname(urls[["fulltaxa_csv"]])))

  testthat::expect_identical(
    eval(formals(bf_sinas_default_urls)$record_id),
    "21933976"
  )
  testthat::expect_identical(
    eval(formals(bf_download_sinas_resources)$record_id),
    "21933976"
  )
  testthat::expect_identical(
    eval(formals(bf_fetch_native_ranges_sinas)$record_id),
    "21933976"
  )
  testthat::expect_identical(
    eval(formals(bf_fetch_native_ranges_web)$sinas_record_id),
    "21933976"
  )
})


testthat::test_that("SInAS 3.1.1 remains explicitly resolvable", {
  urls <- bf_sinas_default_urls("18220953")

  testthat::expect_match(
    urls[["main_csv"]],
    "SInAS_3.1.1.csv",
    fixed = TRUE
  )
  testthat::expect_match(
    urls[["config_zip"]],
    "All_Config_Files_SInAS_v3.1.1.zip",
    fixed = TRUE
  )
  testthat::expect_match(
    urls[["output_zip"]],
    "All_Output_Files_SInAS_v3.1.1.zip",
    fixed = TRUE
  )
})


testthat::test_that("SInAS 3.2 native evidence parses and feeds recipient classification", {
  testthat::skip_if_not_installed("readr")
  testthat::skip_if_not_installed("tibble")
  testthat::skip_if_not_installed("dplyr")

  td <- tempfile("biofetchR_sinas32_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)

  main_path <- file.path(td, "SInAS_3.2.csv")
  alllocations_path <- file.path(td, "AllLocations.csv")
  fulltaxa_path <- file.path(td, "SInAS_3.2_FullTaxaList.csv")

  readr::write_csv(
    tibble::tibble(
      location = c(
        "India",
        "France",
        "United Kingdom",
        "Atlantis",
        "Germany"
      ),
      locationID = c(
        "loc_ind",
        "loc_fra",
        "loc_gbr",
        "loc_atlantis",
        "loc_deu"
      ),
      taxon = c(
        "Rattus rattus",
        "Rattus rattus",
        "Carcinus maenas",
        "Carcinus maenas",
        "Sturnus vulgaris"
      ),
      taxonID = c("101", "101", "202", "202", "303"),
      establishmentMeans = c(
        "native",
        "introduced",
        "native",
        "native",
        "native"
      )
    ),
    main_path
  )

  readr::write_csv(
    tibble::tibble(
      locationID = c("loc_ind", "loc_fra", "loc_gbr", "loc_deu"),
      ISO3 = c("IND", "FRA", "GBR", "DEU"),
      country = c("India", "France", "United Kingdom", "Germany")
    ),
    alllocations_path
  )

  readr::write_csv(
    tibble::tibble(
      scientificName = c(
        "Rattus alexandrinus",
        "Rattus rattus",
        "Carcinus maenas",
        "Sturnus vulgaris"
      ),
      taxon = c(
        "Rattus rattus",
        "Rattus rattus",
        "Carcinus maenas",
        "Sturnus vulgaris"
      )
    ),
    fulltaxa_path
  )

  x <- bf_fetch_native_ranges_sinas(
    species = c(
      "Rattus rattus",
      "Rattus alexandrinus",
      "Carcinus maenas",
      "Sturnus vulgaris"
    ),
    main_path = main_path,
    alllocations_path = alllocations_path,
    fulltaxa_path = fulltaxa_path,
    record_id = "21933976",
    quiet = TRUE,
    return = "list"
  )

  testthat::expect_true(nrow(x$long) > 0L)
  testthat::expect_equal(nrow(x$species), 4L)

  rattus <- x$species[
    x$species$species == "Rattus rattus",
    ,
    drop = FALSE
  ]

  testthat::expect_equal(nrow(rattus), 1L)
  testthat::expect_match(
    rattus$native_origin_iso3[[1L]],
    "IND",
    fixed = TRUE
  )
  testthat::expect_false(
    grepl("FRA", rattus$native_origin_iso3[[1L]], fixed = TRUE)
  )

  alias <- x$species[
    x$species$species == "Rattus alexandrinus",
    ,
    drop = FALSE
  ]

  testthat::expect_equal(nrow(alias), 1L)
  testthat::expect_match(
    alias$native_origin_iso3[[1L]],
    "IND",
    fixed = TRUE
  )

  testthat::expect_true(
    any(
      x$unmapped$species == "Carcinus maenas" &
        x$unmapped$raw_native_area == "Atlantis"
    )
  )

  testthat::expect_identical(x$summary$sinas_version[[1L]], "3.2")
  testthat::expect_identical(
    x$summary$sinas_workflow_version[[1L]],
    "2.0"
  )
  testthat::expect_identical(
    x$summary$sinas_record_id[[1L]],
    "21933976"
  )
  testthat::expect_identical(
    x$summary$sinas_doi[[1L]],
    "10.5281/zenodo.21933976"
  )

  routed <- bf_fetch_native_ranges_web(
    species = c("Rattus rattus", "Carcinus maenas"),
    sources = "sinas",
    sinas_main_path = main_path,
    sinas_alllocations_path = alllocations_path,
    sinas_fulltaxa_path = fulltaxa_path,
    quiet = TRUE,
    return = "list"
  )

  testthat::expect_true(nrow(routed$long) > 0L)
  testthat::expect_identical(
    routed$summary$sinas_version[[1L]],
    "3.2"
  )
  testthat::expect_identical(
    routed$summary$sinas_record_id[[1L]],
    "21933976"
  )

  classified <- bf_attach_native_status(
    df = tibble::tibble(
      species = c("Rattus rattus", "Rattus rattus"),
      iso3c = c("IND", "FRA")
    ),
    species_col = "species",
    iso3c_col = "iso3c",
    native_ranges = x$species,
    native_species_col = "species",
    require_country = TRUE,
    native_filter_mode = "audit_only",
    quiet = TRUE
  )

  testthat::expect_true(
    classified$native_is_native_recipient[
      classified$iso3c == "IND"
    ][[1L]]
  )
  testthat::expect_true(
    classified$native_is_non_native_recipient[
      classified$iso3c == "FRA"
    ][[1L]]
  )
})
